import re
from dataclasses import dataclass
from typing import List, Tuple

# Detect headings as standalone lines.
# Works with:
# - Markdown headings: "# Results"
# - Numbered headings: "1 Introduction", "2.3 Model training"
# - Common paper headings: "Abstract", "Methods", etc.
_HEADING_RE = re.compile(
    r"""
    ^\s*(
        \#\#HEADING\#\#\s+.+ |
        \#{1,6}\s+.+ |
        (?:\d+(?:\.\d+)*)\s+.+ |
        (?:Abstract|Introduction|Background|Methods|Materials and Methods|Results|Discussion|Conclusion|Limitations)\b.*
    )\s*$
    """,
    re.IGNORECASE | re.VERBOSE,
)


def _approx_tokens(s: str) -> int:
    # Cheap proxy. Good enough for chunk sizing.
    return max(1, len(s) // 4)

@dataclass
class Chunk:
    chunk_id: str
    section_path: str
    text: str
    start_char: int
    end_char: int

def chunk_with_headings(
    full_text: str,
    doc_id: str,
    target_tokens: int = 550,
    overlap_tokens: int = 80,
) -> List[Chunk]:
    """
    Chunk by headings, then split section text into paragraph-based chunks.
    Headings must appear on standalone lines in full_text.
    """
    lines = full_text.splitlines(keepends=True)

    sections: List[Tuple[str, str, int, int]] = []
    current_heading = "ROOT"
    buf: List[str] = []
    section_start = 0
    pos = 0

    for line in lines:
        if _HEADING_RE.match(line.strip()):
            # Flush previous section content
            if buf:
                sec_text = "".join(buf).strip()
                if sec_text:
                    sections.append((current_heading, sec_text, section_start, pos))
            # Start new section
            h = line.strip()
            if h.lower().startswith("##heading##"):
                h = h.split(" ", 1)[1].strip() if " " in h else "HEADING"
            current_heading = h
            buf = []
            section_start = pos
        else:
            buf.append(line)
        pos += len(line)

    # Flush last
    if buf:
        sec_text = "".join(buf).strip()
        if sec_text:
            sections.append((current_heading, sec_text, section_start, pos))

    chunks: List[Chunk] = []
    chunk_idx = 0

    for section_path, sec_text, sec_start, sec_end in sections:
        # Split into paragraphs
        paras = [p.strip() for p in re.split(r"\n\s*\n+|\n{1,}", sec_text) if p.strip()]
        if not paras:
            continue

        cur: List[str] = []
        cur_tokens = 0

        def emit(overlap_seed: str = ""):
            nonlocal chunk_idx, cur, cur_tokens
            if not cur:
                return
            text = "\n\n".join(cur).strip()
            if not text:
                cur = []
                cur_tokens = 0
                return

            chunks.append(
                Chunk(
                    chunk_id=f"{doc_id}::chunk_{chunk_idx:05d}",
                    section_path=section_path,
                    text=text,
                    start_char=sec_start,
                    end_char=sec_end,
                )
            )
            chunk_idx += 1

            # Overlap seed for next chunk
            if overlap_seed:
                cur = [overlap_seed]
                cur_tokens = _approx_tokens(overlap_seed)
            else:
                cur = []
                cur_tokens = 0

        for p in paras:
            p_tokens = _approx_tokens(p)

            # If a paragraph is huge, split roughly by sentences.
            if p_tokens > target_tokens:
                emit()
                sentences = re.split(r"(?<=[.!?])\s+", p)
                tmp: List[str] = []
                tmp_tokens = 0

                for s in sentences:
                    s = s.strip()
                    if not s:
                        continue
                    s_tokens = _approx_tokens(s)

                    if tmp_tokens + s_tokens > target_tokens and tmp:
                        text = " ".join(tmp).strip()
                        overlap = text[-overlap_tokens * 4 :]
                        chunks.append(
                            Chunk(
                                chunk_id=f"{doc_id}::chunk_{chunk_idx:05d}",
                                section_path=section_path,
                                text=text,
                                start_char=sec_start,
                                end_char=sec_end,
                            )
                        )
                        chunk_idx += 1
                        tmp = [overlap]
                        tmp_tokens = _approx_tokens(overlap)

                    tmp.append(s)
                    tmp_tokens += s_tokens

                if tmp:
                    chunks.append(
                        Chunk(
                            chunk_id=f"{doc_id}::chunk_{chunk_idx:05d}",
                            section_path=section_path,
                            text=" ".join(tmp).strip(),
                            start_char=sec_start,
                            end_char=sec_end,
                        )
                    )
                    chunk_idx += 1
                continue

            # Normal case: build chunk until budget, then emit with overlap
            if cur and (cur_tokens + p_tokens > target_tokens):
                prev_text = "\n\n".join(cur).strip()
                overlap = prev_text[-overlap_tokens * 4 :]
                emit(overlap_seed=overlap)

            cur.append(p)
            cur_tokens += p_tokens

        emit()

    return chunks
