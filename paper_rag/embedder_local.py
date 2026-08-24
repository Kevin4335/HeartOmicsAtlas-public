import numpy as np
from typing import List

try:
    from sentence_transformers import SentenceTransformer
except ImportError:
    SentenceTransformer = None


class LocalEmbedder:
    """Local embedder using sentence-transformers. No API calls needed."""
    
    def __init__(self, model: str = "all-MiniLM-L6-v2"):
        if SentenceTransformer is None:
            raise RuntimeError(
                "sentence-transformers not installed. "
                "Run: pip install sentence-transformers"
            )
        self.model = SentenceTransformer(model)
        self.dimension = self.model.get_embedding_dimension()
    
    def embed_texts(self, texts: List[str]) -> np.ndarray:
        """Embed a list of texts and return as float32 numpy array."""
        vecs = self.model.encode(texts, convert_to_numpy=True)
        return vecs.astype(np.float32)
