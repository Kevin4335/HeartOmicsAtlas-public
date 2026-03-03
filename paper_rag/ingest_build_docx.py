import os
import hashlib
from paper_rag.docx_extract import extract_docx_with_headings, to_marked_text
from paper_rag.chunking import chunk_with_headings
from paper_rag.embedder_openai import OpenAIEmbedder
from paper_rag.faiss_store import save_jsonl, build_faiss_index, save_faiss

def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def main():
    docx_path = os.environ.get(
        "PAPER_DOCX",
        "/mnt/data/Copy of fetal_heart_manuscript_Jiajun_20260203.docx",
    )
    out_dir = os.environ.get("RAG_OUT_DIR", "data/rag/fetal_heart_paper_v1")

    if not os.path.exists(docx_path):
        raise FileNotFoundError(f"Missing DOCX at: {docx_path}")

    items = extract_docx_with_headings(docx_path)
    marked_text = to_marked_text(items)

    doc_hash = file_sha256(docx_path)[:16]
    doc_id = f"fetal_heart_paper::{doc_hash}"

    chunks = chunk_with_headings(
        full_text=marked_text,
        doc_id=doc_id,
        target_tokens=int(os.environ.get("CHUNK_TOKENS", "550")),
        overlap_tokens=int(os.environ.get("CHUNK_OVERLAP", "80")),
    )

    if len(chunks) < 30:
        raise RuntimeError(
            f"Only {len(chunks)} chunks. This usually means headings were not detected. "
            f"Fix the DOCX heading styles (Heading 1/2/3) or add explicit heading lines."
        )

    rows = []
    texts = []
    for c in chunks:
        rows.append(
            {
                "doc_id": doc_id,
                "chunk_id": c.chunk_id,
                "section_path": c.section_path,
                "start_char": c.start_char,
                "end_char": c.end_char,
                "text": c.text,
            }
        )
        texts.append(f"[{c.section_path}]\n{c.text}")

    embedder = OpenAIEmbedder()
    vecs = embedder.embed_texts(texts)

    index = build_faiss_index(vecs)

    os.makedirs(out_dir, exist_ok=True)
    meta_path = os.path.join(out_dir, "metadata.jsonl")
    faiss_path = os.path.join(out_dir, "index.faiss")

    save_jsonl(meta_path, rows)
    save_faiss(index, faiss_path)

    print(f"OK: doc_id={doc_id}")
    print(f"Wrote: {meta_path}")
    print(f"Wrote: {faiss_path}")
    print(f"Chunks: {len(rows)}")

if __name__ == "__main__":
    main()
