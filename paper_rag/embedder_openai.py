import os
import numpy as np
from typing import List
from config import CHAT_KEY
from openai import OpenAI


class OpenAIEmbedder:
    def __init__(self, model: str = "text-embedding-3-small"):
        api_key = CHAT_KEY
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY not set in environment.")
        self.client = OpenAI(api_key=api_key)
        self.model = model

    def embed_texts(self, texts: List[str]) -> np.ndarray:
        resp = self.client.embeddings.create(
            model=self.model,
            input=texts
        )
        vecs = [x.embedding for x in resp.data]
        return np.array(vecs, dtype=np.float32)
