# paper_rag/paper_search.py

from __future__ import annotations

import json
import os
import re
import threading
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:
    import faiss
except Exception:
    faiss = None

try:
    from rank_bm25 import BM25Okapi
except ImportError:
    BM25Okapi = None

from paper_rag.embedder_openai import OpenAIEmbedder


# ---- Optional blacklist for junk sections ----
DEFAULT_EXCLUDE_SECTIONS = {
    "Author Contributions",
    "Acknowledgements",
    "Acknowledgments",
    "Code availability.",
    "Code availability",
    "Data availability",
    "Competing interests",
    "Ethics statement",
}


@dataclass
class PaperChunk:
    doc_id: str
    chunk_id: str
    section_path: str
    text: str
    start_char: int
    end_char: int


class PaperSearchEngine:
    def __init__(
        self,
        rag_dir: Optional[str] = None,
        index_path: Optional[str] = None,
        meta_path: Optional[str] = None,
        embed_model: str = "text-embedding-3-small",
        min_score: float = 0.0,
        exclude_sections: Optional[set] = None,
    ):
        if faiss is None:
            raise RuntimeError("faiss not installed. pip install faiss-cpu")

        self.min_score = float(min_score)
        self.exclude_sections = exclude_sections or DEFAULT_EXCLUDE_SECTIONS

        _default_rag_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "data", "fetal_heart_v1",
        )
        self.rag_dir = rag_dir or os.environ.get("RAG_DIR", _default_rag_dir)
        self.index_path = index_path or os.environ.get(
            "RAG_INDEX", os.path.join(self.rag_dir, "index.faiss")
        )
        self.meta_path = meta_path or os.environ.get(
            "RAG_META", os.path.join(self.rag_dir, "metadata.jsonl")
        )

        self.embedder = OpenAIEmbedder(model=embed_model)

        self._index = None
        self._chunks: List[PaperChunk] = []
        self._bm25 = None

        self._load()

    def _load(self) -> None:
        if not os.path.exists(self.index_path):
            raise FileNotFoundError(self.index_path)
        if not os.path.exists(self.meta_path):
            raise FileNotFoundError(self.meta_path)

        self._index = faiss.read_index(self.index_path)

        with open(self.meta_path, "r", encoding="utf-8") as f:
            for line in f:
                obj = json.loads(line)
                self._chunks.append(
                    PaperChunk(
                        doc_id=obj.get("doc_id", ""),
                        chunk_id=obj["chunk_id"],
                        section_path=obj.get("section_path", "UNKNOWN"),
                        text=obj["text"],
                        start_char=int(obj.get("start_char", 0)),
                        end_char=int(obj.get("end_char", 0)),
                    )
                )

        # Build BM25 index in-memory from loaded chunks (fast, no disk storage)
        if BM25Okapi is not None and self._chunks:
            tokenized = [ch.text.lower().split() for ch in self._chunks]
            self._bm25 = BM25Okapi(tokenized)

    # ------------------------------------------------------------------
    # Dense (FAISS) search
    # ------------------------------------------------------------------

    def _raw_search(
        self,
        query: str,
        k: int,
        *,
        min_score: Optional[float] = None,
        exclude_sections: Optional[set] = None,
    ) -> List[Dict[str, Any]]:
        min_score_val = self.min_score if min_score is None else float(min_score)
        exclude_val = self.exclude_sections if exclude_sections is None else exclude_sections
        if not query.strip():
            return []

        vec = self.embedder.embed_texts([query]).astype(np.float32)
        faiss.normalize_L2(vec)

        scores, idxs = self._index.search(vec, k)
        scores = scores[0]
        idxs = idxs[0]

        hits = []
        for score, idx in zip(scores, idxs):
            if idx < 0:
                continue
            ch = self._chunks[idx]

            if float(score) < min_score_val:
                continue
            if ch.section_path in exclude_val:
                continue

            hits.append(
                {
                    "doc_id": ch.doc_id,
                    "chunk_id": ch.chunk_id,
                    "section_path": ch.section_path,
                    "score": float(score),
                    "text": ch.text,
                    "_row": int(idx),
                }
            )

        return hits

    # ------------------------------------------------------------------
    # BM25 keyword search
    # ------------------------------------------------------------------

    def _bm25_search(
        self,
        query: str,
        top_n: int = 50,
        exclude_sections: Optional[set] = None,
    ) -> List[Dict[str, Any]]:
        """Return top_n chunks ranked by BM25 keyword score."""
        if self._bm25 is None:
            return []
        if not query.strip():
            return []

        exclude_val = self.exclude_sections if exclude_sections is None else exclude_sections
        tokens = query.lower().split()
        scores = self._bm25.get_scores(tokens)  # numpy array, one per chunk

        ranked = np.argsort(scores)[::-1][:top_n]
        hits = []
        for row in ranked:
            if scores[row] <= 0:
                break
            ch = self._chunks[row]
            if ch.section_path in exclude_val:
                continue
            hits.append(
                {
                    "doc_id": ch.doc_id,
                    "chunk_id": ch.chunk_id,
                    "section_path": ch.section_path,
                    "score": float(scores[row]),
                    "text": ch.text,
                    "_row": int(row),
                }
            )
        return hits

    # ------------------------------------------------------------------
    # Reciprocal Rank Fusion
    # ------------------------------------------------------------------

    def _reciprocal_rank_fusion(
        self,
        dense_hits: List[Dict[str, Any]],
        bm25_hits: List[Dict[str, Any]],
        k: int = 60,
    ) -> List[Dict[str, Any]]:
        """Merge two ranked lists using RRF: score = sum of 1/(rank + k)."""
        rrf_scores: Dict[int, float] = {}
        all_hits: Dict[int, Dict] = {}

        for rank, h in enumerate(dense_hits):
            row = h["_row"]
            rrf_scores[row] = rrf_scores.get(row, 0.0) + 1.0 / (rank + k)
            all_hits[row] = h

        for rank, h in enumerate(bm25_hits):
            row = h["_row"]
            rrf_scores[row] = rrf_scores.get(row, 0.0) + 1.0 / (rank + k)
            if row not in all_hits:
                all_hits[row] = h

        merged = sorted(rrf_scores.items(), key=lambda x: x[1], reverse=True)
        return [
            {**all_hits[row], "score": score, "_row": row}
            for row, score in merged
        ]

    # ------------------------------------------------------------------
    # Post-processing helpers (unchanged)
    # ------------------------------------------------------------------

    def _apply_section_cap(self, hits, k, max_per_section):
        out = []
        per_sec = {}

        for h in hits:
            sec = h["section_path"]
            if per_sec.get(sec, 0) >= max_per_section:
                continue
            per_sec[sec] = per_sec.get(sec, 0) + 1
            out.append(h)
            if len(out) >= k:
                break

        return out

    def _add_context_window(self, selected, context_window, context_top_n):
        if context_window <= 0 or context_top_n <= 0:
            return selected

        seen = set()
        out = []

        def add(h):
            r = h["_row"]
            if r not in seen:
                seen.add(r)
                out.append(h)

        top_n = min(context_top_n, len(selected))

        for i, h in enumerate(selected):
            add(h)

            if i >= top_n:
                continue

            row = h["_row"]
            base_score = h["score"]

            for d in range(1, context_window + 1):
                for nbr in (row - d, row + d):
                    if nbr < 0 or nbr >= len(self._chunks):
                        continue
                    if nbr in seen:
                        continue

                    ch = self._chunks[nbr]
                    if ch.section_path != h["section_path"]:
                        continue

                    add(
                        {
                            "doc_id": ch.doc_id,
                            "chunk_id": ch.chunk_id,
                            "section_path": ch.section_path,
                            "score": base_score - 0.05 * d,
                            "text": ch.text,
                            "_row": nbr,
                            "_context_of": h["chunk_id"],
                        }
                    )

        return out

    # ------------------------------------------------------------------
    # Public search: hybrid dense + BM25 with RRF
    # ------------------------------------------------------------------

    def search(
        self,
        query: str,
        k: int = 8,
        max_per_section: int = 2,
        context_window: int = 0,
        context_top_n: int = 0,
        *,
        min_score: Optional[float] = None,
        exclude_sections: Optional[set] = None,
    ):
        # Dense retrieval — use min_score=0.0 so RRF sees the full ranked list
        dense_hits = self._raw_search(
            query,
            max(50, k * 10),
            min_score=0.0,
            exclude_sections=exclude_sections,
        )

        # BM25 keyword retrieval
        bm25_hits = self._bm25_search(query, top_n=50, exclude_sections=exclude_sections)

        # Merge via Reciprocal Rank Fusion
        raw = self._reciprocal_rank_fusion(dense_hits, bm25_hits)

        # Post-processing: section cap → context window → trim
        primary = self._apply_section_cap(raw, k, max_per_section)
        expanded = self._add_context_window(primary, context_window, context_top_n)

        out = []
        for h in expanded[:k]:
            h2 = dict(h)
            h2.pop("_row", None)
            out.append(h2)
        return out


_ENGINE = None
_LOCK = threading.Lock()


def paper_search(query, **kwargs):
    global _ENGINE
    if _ENGINE is None:
        init_paper_search()
    return _ENGINE.search(query, **kwargs)

def init_paper_search(
    rag_dir: Optional[str] = None,
    index_path: Optional[str] = None,
    meta_path: Optional[str] = None,
    embed_model: str = "text-embedding-3-small",
    *,
    min_score: float = 0.0,
    exclude_sections: Optional[set] = None,
) -> None:
    """
    Eagerly initialize the singleton search engine.
    Loads index.faiss + metadata.jsonl into memory and builds BM25 index.
    Does NOT rebuild embeddings or FAISS.
    """
    global _ENGINE
    with _LOCK:
        if _ENGINE is None:
            _ENGINE = PaperSearchEngine(
                rag_dir=rag_dir,
                index_path=index_path,
                meta_path=meta_path,
                embed_model=embed_model,
                min_score=min_score,
                exclude_sections=exclude_sections,
            )
