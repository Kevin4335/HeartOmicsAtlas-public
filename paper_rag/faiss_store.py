import json
import os
import numpy as np
import faiss


def save_jsonl(path: str, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


def build_faiss_index(vectors: np.ndarray):
    n, d = vectors.shape
    index = faiss.IndexFlatIP(d)
    faiss.normalize_L2(vectors)
    index.add(vectors)
    return index


def save_faiss(index, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    faiss.write_index(index, path)
