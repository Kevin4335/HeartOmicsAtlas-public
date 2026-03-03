from R_http import *
from utils import *
from mySecrets import hexToStr
from openai import OpenAI
import time
from queue import Queue
import json
from random import randint
from typing import Union, Literal, Tuple
from myBasics import binToBase64
from hashlib import sha256
from queue import Queue
from _thread import start_new_thread
from copy import deepcopy
import anthropic
import traceback
import secrets
import openai
import requests
from config import API_KEY
from config import CHAT_KEY
import os
from pathlib import Path
from typing import Tuple, Union, Literal, List, Dict, Any
import certifi
import re
from paper_rag.paper_search import paper_search, init_paper_search
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    import dotenv
    _root = Path(__file__).resolve().parent
    dotenv.load_dotenv(_root / ".env")
except ImportError:
    pass

__all__ = ['process_ai_chat']

PLOT_BACKEND_BASE = os.environ.get("PLOT_BACKEND_BASE", "http://128.84.41.80")
GLKB_LLM_AGENT_URL = os.environ.get("GLKB_LLM_AGENT_URL", "https://glkb.dcmb.med.umich.edu/api/frontend/llm_agent")

# ---------------------------------------------------------------------------
# Tool definitions for Claude native tool use
# ---------------------------------------------------------------------------

TOOLS = [
    {
        "name": "scRNA",
        "description": (
            "Show gene expression (UMAP, Violin, Dotplot, and split UMAP) in single cell RNA-seq. "
            "gene should be a gene ID (or genetic loci), case insensitive. "
            "type selects the organoid / cell population."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene ID or genetic locus, case insensitive.",
                },
                "type": {
                    "type": "string",
                    "enum": [
                        "Sinoid (SAN)",
                        "ACO (Atrial Cardioids)",
                        "VCO (Ventricular Cardioids)",
                        "SAN-PACO (SAN Paced Atrial Cardioids)",
                        "Mini-heart",
                    ],
                    "description": "The organoid or cell population to query.",
                },
            },
            "required": ["gene", "type"],
        },
    },
    {
        "name": "multiomics",
        "description": (
            "Show gene expression (UMAP, Violin plot) in snRNA-seq of multiomics, "
            "and also the chromosome accessibility (IGV coverage plot) in snATAC-seq of multiomics. "
            "gene should be a gene ID (or genetic loci), case insensitive."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene ID or genetic locus, case insensitive.",
                },
            },
            "required": ["gene"],
        },
    },
    {
        "name": "spatial_transcriptomics",
        "description": (
            "Show the gene's spatial expression in the human fetal heart, including "
            "spatial feature plot, UMAP, Violin plot, and Dot plot. "
            "gene should be a gene ID, case insensitive."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene ID, case insensitive.",
                },
            },
            "required": ["gene"],
        },
    },
    {
        "name": "static_images",
        "description": (
            "Show static default overview images to the user. "
            "All images are combined into single PNG files. name is case sensitive."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {
                    "type": "string",
                    "enum": [
                        "Sinoid (SAN)",
                        "ACO (Atrial Cardioids)",
                        "VCO (Ventricular Cardioids)",
                        "SAN-PACO (SAN Paced Atrial Cardioids)",
                        "Mini-heart",
                        "multiomics",
                        "spatial_transcriptomics",
                    ],
                    "description": (
                        "Which default image set to show. "
                        "Sinoid/ACO/VCO = scRNA defaults, SAN-PACO = SAN-PACO scRNA defaults, "
                        "Mini-heart = Mini-heart scRNA defaults, multiomics = multiomics defaults, "
                        "spatial_transcriptomics = spatial transcriptomics defaults."
                    ),
                },
            },
            "required": ["name"],
        },
    },
    {
        "name": "glkb_ai_assistant",
        "description": (
            "Ask a question to the GLKB (Genomic Literature Knowledge Base) AI assistant. "
            "The GLKB assistant can only see the question string — it has no access to "
            "chat history, so provide all necessary background and context in the question. "
            "The result will be displayed to the user directly."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "question": {
                    "type": "string",
                    "description": "The full question to send to GLKB, including all needed context.",
                },
            },
            "required": ["question"],
        },
    },
]

# ---------------------------------------------------------------------------
# Planner tool — used exclusively in Phase 1 to force a structured plan output
# ---------------------------------------------------------------------------

PLANNER_TOOL = [
    {
        "name": "create_plan",
        "description": (
            "Output a structured execution plan for the user's request. "
            "List every tool call that needs to be made. All listed tool calls will be "
            "executed in parallel, so do not list calls that depend on each other's output. "
            "If no tools are needed (pure text/question), output an empty steps list."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "intro_text": {
                    "type": "string",
                    "description": (
                        "A brief natural-language message to show the user while the tools run, "
                        "e.g. 'I'll generate the scRNA and spatial plots for INS now.'. "
                        "Leave empty string if no tools are needed."
                    ),
                },
                "steps": {
                    "type": "array",
                    "description": "List of tool calls to execute (all run in parallel).",
                    "items": {
                        "type": "object",
                        "properties": {
                            "tool": {
                                "type": "string",
                                "enum": [
                                    "scRNA",
                                    "multiomics",
                                    "spatial_transcriptomics",
                                    "static_images",
                                    "glkb_ai_assistant",
                                ],
                                "description": "Name of the tool to call.",
                            },
                            "input": {
                                "type": "object",
                                "description": "Input arguments for the tool, matching that tool's schema.",
                            },
                        },
                        "required": ["tool", "input"],
                    },
                },
            },
            "required": ["steps"],
        },
    }
]

PLANNER_SYSTEM_PROMPT = """You are the planner for HeartOmicsAtlas AI assistant.

Your ONLY job is to call the `create_plan` tool with a list of tool calls needed to answer the user's request.

## Available tools and their exact parameters

### scRNA
Visualize gene expression (UMAP, Violin, Dotplot, split UMAP) in single-cell RNA-seq.
Parameters:
  - gene: string — gene symbol, e.g. "INS", "HCN4", "TBX3"
  - type: one of EXACTLY these values (case-sensitive):
      "Sinoid (SAN)"
      "ACO (Atrial Cardioids)"
      "VCO (Ventricular Cardioids)"
      "SAN-PACO (SAN Paced Atrial Cardioids)"
      "Mini-heart"

### multiomics
Visualize snRNA-seq + snATAC-seq (chromatin accessibility) data.
Parameters:
  - gene: string — gene symbol, e.g. "INS", "HCN4"

### spatial_transcriptomics
Visualize spatial gene expression in the human fetal heart (spatial feature plot, UMAP, Violin, Dot plot).
Parameters:
  - gene: string — gene symbol, e.g. "INS", "HCN4"

### static_images
Show a default overview image for a dataset. Use when user asks to "show me the overview", "what does the data look like", or asks about a dataset without specifying a gene.
Parameters:
  - name: one of EXACTLY these values (case-sensitive):
      "Sinoid (SAN)"             — scRNA default for Sinoid/ACO/VCO
      "ACO (Atrial Cardioids)"   — scRNA default for Sinoid/ACO/VCO
      "VCO (Ventricular Cardioids)" — scRNA default for Sinoid/ACO/VCO
      "SAN-PACO (SAN Paced Atrial Cardioids)" — scRNA default for SAN-PACO
      "Mini-heart"               — scRNA default for Mini-heart
      "multiomics"               — multiomics default overview
      "spatial_transcriptomics"  — spatial transcriptomics default overview

### glkb_ai_assistant
Search the Genomic Literature Knowledge Base (33M+ PubMed abstracts) for external biology knowledge.
Parameters:
  - question: string — the full self-contained question with all needed context

## Rules

**General:**
- All steps run in parallel — only list independent tool calls.
- ALWAYS use exact enum values listed above for `type` and `name` parameters. Never invent values.

**Visualization tools (scRNA, multiomics, spatial_transcriptomics, static_images):**
- Call these when the user asks to visualize, show, plot, or display gene expression or data.
- If the user asks for a gene across multiple data types, list each as a separate step.

**glkb_ai_assistant:**
- Call this for questions about gene function, pathway biology, transcription factors, signaling mechanisms, or disease — anything that benefits from external literature context.
- Do NOT call it for purely website-specific questions like "what tabs does HeartOmicsAtlas have?" or "how do I use this website?".

**Combining tools — this is critical:**
- When a user asks a biology question about a gene or pathway, you should OFTEN plan BOTH glkb_ai_assistant AND visualization tools together. A good answer blends general biology knowledge with our dataset's specific findings.
- The synthesizer will combine GLKB's external knowledge with RAG paper evidence to produce a comprehensive answer.

**steps: [] (no tools):**
- Use ONLY for greetings, "thank you", clarifications, or purely website-usage questions.
- For biology questions where the RAG context alone is sufficient (e.g. "what are the main findings of this paper?"), steps: [] is OK — the synthesizer has the paper evidence.

## Examples

User: "What is the biological function of SHOX2?"
Plan: glkb_ai_assistant(question="What is the biological function of SHOX2 in cardiac development and the sinoatrial node?") + scRNA(gene="SHOX2", type="Sinoid (SAN)")
Reasoning: GLKB provides general biology of SHOX2; scRNA shows its expression pattern in our SAN data.

User: "Which transcription factors regulate SAN differentiation?"
Plan: glkb_ai_assistant(question="Which transcription factors regulate sinoatrial node SAN differentiation in the developing heart?")
Reasoning: broad biology question; RAG paper evidence will add our dataset-specific TF findings (TEAD1, YAP pathway, etc.) in the synthesizer.

User: "Show me the expression of HCN4 in SAN organoids"
Plan: scRNA(gene="HCN4", type="Sinoid (SAN)")
Reasoning: pure visualization request, no external knowledge needed.

User: "What does HCN4 do and show me its expression"
Plan: glkb_ai_assistant(question="What is the function of HCN4 in cardiac pacemaking?") + scRNA(gene="HCN4", type="Sinoid (SAN)")
Reasoning: user wants both explanation and visualization.

User: "Compare TEAD1 expression across scRNA, multiomics, and spatial data"
Plan: scRNA(gene="TEAD1", type="Sinoid (SAN)") + multiomics(gene="TEAD1") + spatial_transcriptomics(gene="TEAD1")
Reasoning: user wants the same gene across all modalities.

User: "Which signaling pathways govern SAN subdomain specification?"
Plan: glkb_ai_assistant(question="Which signaling pathways govern sinoatrial node SAN subdomain specification in the developing heart?")
Reasoning: broad biology + pathway question; RAG will add our Hippo pathway findings.

User: "What are the main findings of this paper?"
Plan: steps: []
Reasoning: purely about this paper; RAG paper evidence already in the synthesizer prompt answers this.

User: "Show me an overview of the multiomics data"
Plan: static_images(name="multiomics")
Reasoning: user asks for a default overview image.

User: "Hi, what can you do?"
Plan: steps: []
Reasoning: conversational, no tools needed.

Do NOT include any explanation in your response — only call create_plan.
"""

# ---------------------------------------------------------------------------
# System prompt (simplified — tool schemas handled by Claude native tool use)
# ---------------------------------------------------------------------------

PROMPT = """## 1. Introduction and Task

You are the AI assistant of HeartOmicsAtlas website. This website is for the display of some biology data, about human's fetal heart. Its main functions include showing the gene expression from the scRNA data and spatial distribution from spatial transcriptomics, as well as the chromosome accessibility from the snATAC data. You need to chat with users, answer their questions, and generate images (by calling the provided tools) based on their requirements.

HeartOmicsAtlas is an AI-powered, user-friendly, open-access platform designed to analyzing single cell RNA-seq (scRNA-seq), Multiomics, and Spatial Transcriptomics data. By integrating these datasets, HeartOmicsAtlas provides a comprehensive framework for systematically analyzing cellular responses and molecular changes in human fetal heart cells throughout development.

You can also ask questions to GLKB AI assistant. The Genomic Literature Knowledge Base (GLKB) is a comprehensive and powerful resource that integrates over 263 million biomedical terms and more than 14.6 million biomedical relationships. This collection is curated from 33 million PubMed abstracts and nine well-established biomedical repositories, offering an unparalleled wealth of knowledge for researchers and practitioners in the field. The AI assistant of GLKB can search the database and answer user's question based on the database. If user asks some questions in the biology field but not related to this website, you can use GLKB to answer user's question.

## 2. Tool Usage

When the user asks to visualize, show, plot, or display gene expression, multiomics, spatial data, or default images, you MUST use the provided tools. Do NOT write JSON or function call syntax in your text — use the actual tool calling mechanism provided to you. When you want to generate a visualization, invoke the appropriate tool directly.

- When the user asks for scRNA visualization of a specific gene and organoid type, call the `scRNA` tool.
- When the user asks for multiomics data for a gene, call the `multiomics` tool.
- When the user asks for spatial transcriptomics data for a gene, call the `spatial_transcriptomics` tool.
- When the user asks to see default/overview images, call the `static_images` tool.
- When the user asks a biology question not covered by the above, call the `glkb_ai_assistant` tool.

After each image tool call, you will receive the actual generated figure as an image. Examine it and provide a brief biological interpretation: describe what the expression pattern shows (e.g. which cell clusters express the gene, what the UMAP distribution looks like, notable features of the violin plot). Keep the interpretation concise and relevant to the user's question.

## 3. How to Answer Biology Questions

When answering biology questions (about genes, pathways, transcription factors, signaling, cell types, etc.):

1. **Start with general knowledge**: Provide a concise overview of the gene/pathway/concept from general biology. If GLKB results are available, use them.
2. **Then cite this paper's specific findings**: Transition with phrases like "In our dataset...", "In our fetal heart samples...", "In our HeartOmicsAtlas data...", "We identified...", "We found...". Highlight what this study specifically discovered or confirmed.
3. **Interpret any generated figures**: If visualization results are available, describe what the plots show — which cell clusters express the gene, expression levels, spatial patterns, UMAP distribution, etc.
4. **Highlight novel findings**: If the paper evidence contains unique discoveries (e.g., "TEAD1 and the YAP/Hippo pathway are critical for SAN specification"), emphasize these as key contributions of this study.

Example answer structure for "What is the biological function of SHOX2?":
- First: general biology of SHOX2 (homeobox TF, roles in embryonic development, cardiac conduction)
- Then: "In our dataset, we use SHOX2 as the main marker to identify SAN. Both in vivo fetal heart samples and in vitro hPSC-derived organoids show consistent high expression of SHOX2."
- Then: interpret the scRNA plot if one was generated

## 4. Behavioral Guidelines

- This is a chatting mode — remember previous messages and answer based on conversation history.
- Interact with the user directly (use "you" and "I").
- Always respond in natural language. NEVER output JSON arrays or objects.
- If you received image results, examine them and provide a brief biological interpretation.
- If user asks a vague question where the organoid type is unclear, ask them to choose (e.g., "Would you like to see this in Sinoid (SAN), ACO, VCO, SAN-PACO, or Mini-heart?").
- Answer all questions in the context of HeartOmicsAtlas, but you are still allowed to answer completely unrelated questions."""

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_PATH = os.path.join(BASE_DIR, "openai_logs.txt")

client = anthropic.Client(api_key=API_KEY)
init_paper_search()

rate_limit_records = []


def within_rate_limit():
    t = time.time()
    for i in range(len(rate_limit_records)-1, -1, -1):
        if (t - rate_limit_records[i] > 3600):
            rate_limit_records.pop(i)
    if(len(rate_limit_records) >= 100):
        return False
    rate_limit_records.append(t)
    return True


# ---------------------------------------------------------------------------
# History normalisation  (flat frontend format <-> Anthropic format)
# ---------------------------------------------------------------------------

def _flat_to_anthropic(flat_history: list) -> list:
    """Convert the frontend's [{role, content: str}] into Anthropic-compatible messages.

    Strips any non-user/assistant roles and ensures content is always a string.
    """
    msgs: list = []
    for m in flat_history:
        role = m.get("role")
        content = m.get("content", "")
        if role not in ("user", "assistant"):
            continue
        if not isinstance(content, str):
            content = str(content)
        if not content.strip():
            continue
        msgs.append({"role": role, "content": content})

    # Anthropic requires messages to alternate user/assistant.
    # Merge consecutive same-role messages.
    merged: list = []
    for m in msgs:
        if merged and merged[-1]["role"] == m["role"]:
            merged[-1]["content"] += "\n" + m["content"]
        else:
            merged.append(dict(m))

    # Must start with a user message
    while merged and merged[0]["role"] != "user":
        merged.pop(0)

    return merged


def _flatten_assistant_text(content_blocks) -> str:
    """Extract plain text from an Anthropic assistant response content list."""
    parts = []
    for block in content_blocks:
        if hasattr(block, "text"):
            parts.append(block.text)
        elif isinstance(block, dict) and block.get("type") == "text":
            parts.append(block.get("text", ""))
    return "\n\n".join(parts)


# ---------------------------------------------------------------------------
# RAG helpers  (unchanged)
# ---------------------------------------------------------------------------

def _latest_user_text(history: list) -> str:
    for msg in reversed(history):
        if isinstance(msg, dict) and msg.get("role") == "user":
            c = msg.get("content", "")
            if isinstance(c, str) and c.strip():
                return c.strip()
            if isinstance(c, list):
                for block in c:
                    if isinstance(block, dict) and block.get("type") == "text":
                        t = block.get("text", "").strip()
                        if t:
                            return t
    return ""

def _excerpt_around_term(text: str, term: str, window: int = 350) -> str:
    if not term:
        return text[:900] + ("..." if len(text) > 900 else "")

    m = re.search(re.escape(term), text, flags=re.IGNORECASE)
    if not m:
        return text[:900] + ("..." if len(text) > 900 else "")

    start = max(0, m.start() - window)
    end = min(len(text), m.end() + window)
    excerpt = text[start:end].strip()

    prefix = "... " if start > 0 else ""
    suffix = " ..." if end < len(text) else ""
    return prefix + excerpt + suffix


def _format_paper_evidence(hits: list, query: str) -> str:
    if not hits:
        return "No relevant paper excerpts found."

    term = ""
    q = (query or "").strip()
    if 1 <= len(q) <= 40 and " " not in q:
        term = q
    else:
        m = re.search(r"\bis\s+([A-Za-z0-9_-]{2,40})\s+mentioned\b", q, flags=re.IGNORECASE)
        if m:
            term = m.group(1)

    lines = []
    for i, h in enumerate(hits, 1):
        sec = h.get("section_path", "UNKNOWN")
        cid = h.get("chunk_id", "UNKNOWN")
        txt = (h.get("text", "") or "").strip()
        txt = _excerpt_around_term(txt, term, window=350)
        lines.append(f"[{i}] section={sec} chunk_id={cid}\n{txt}")
    return "\n\n".join(lines)


# ---------------------------------------------------------------------------
# Tool execution
# ---------------------------------------------------------------------------

# Timeout for fetching plot images from the R backend (seconds).
# R plot generation can be slow on first call; 120s gives it room to breathe.
_PLOT_FETCH_TIMEOUT = 120


def _fetch_png_bytes(url: str) -> Union[bytes, None]:
    """Fetch PNG bytes from the R plot server.

    Returns bytes on success, None on any error (timeout, HTTP error, etc.).
    """
    try:
        resp = requests.get(url, timeout=_PLOT_FETCH_TIMEOUT, verify=False)
        if resp.status_code == 200 and resp.content:
            return resp.content
        print(f"_fetch_png_bytes: HTTP {resp.status_code} from {url}")
        return None
    except Exception as e:
        print(f"_fetch_png_bytes: error fetching {url}: {e}")
        return None


def _image_tool_result(description: str, png_bytes: Union[bytes, None], url: str) -> Tuple[list, dict]:
    """Build a multimodal tool_result content list and a frontend display_msg.

    Claude receives the actual image as a base64 block (so it can interpret the figure).
    The frontend always gets the original URL as the image src — this keeps the JSON
    response small and avoids the React imgStatus key-collision bug that occurs when
    multi-megabyte base64 strings are used as state object keys.

    Returns:
        (tool_result_content, display_msg)
        - tool_result_content: list of content blocks for the Anthropic tool_result
        - display_msg: dict for the frontend {"type": "image", "content": url}
    """
    if png_bytes:
        b64 = binToBase64(png_bytes)
        tool_result_content = [
            {
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": "image/png",
                    "data": b64,
                },
            },
            {
                "type": "text",
                "text": description,
            },
        ]
    else:
        tool_result_content = [{"type": "text", "text": f"{description} (image could not be fetched)"}]

    # Always use the original URL for the frontend — keeps JSON small and imgStatus keys short
    display_msg = {"type": "image", "content": url}

    return tool_result_content, display_msg


def execute_tool(name: str, tool_input: dict) -> Tuple[Union[str, list], Union[dict, None]]:
    """Execute a single tool call.

    Returns:
        (tool_result_content, display_msg_or_None)
        - tool_result_content: str or list of content blocks fed back to Claude as tool_result.
          For image-producing tools this is a multimodal list [image_block, text_block] so
          Claude can see and reason about the generated figure.
        - display_msg: optional dict for the frontend ({"type": "image"/"text", "content": ...})
    """
    base = PLOT_BACKEND_BASE

    if name == "scRNA":
        gene = tool_input["gene"]
        scrna_type = tool_input["type"]
        if scrna_type in ("Sinoid (SAN)", "ACO (Atrial Cardioids)", "VCO (Ventricular Cardioids)"):
            subtype = (
                "sinoid" if scrna_type == "Sinoid (SAN)"
                else ("aco" if scrna_type == "ACO (Atrial Cardioids)" else "vco")
            )
            url = f"{base}:9027/genes/{gene}?subtype={subtype}"
        elif scrna_type == "SAN-PACO (SAN Paced Atrial Cardioids)":
            url = f"{base}:9028/genes/{gene}"
        elif scrna_type == "Mini-heart":
            url = f"{base}:9029/genes/{gene}"
        else:
            url = f"{base}:9027/genes/{gene}?subtype=sinoid"
        desc = f"scRNA plot for gene={gene}, type={scrna_type}"
        png_bytes = _fetch_png_bytes(url)
        tool_result_content, display_msg = _image_tool_result(desc, png_bytes, url)
        return tool_result_content, display_msg

    if name == "multiomics":
        gene = tool_input["gene"]
        url = f"{base}:9026/genes/{gene}"
        desc = f"Multiomics plot (snRNA + snATAC) for gene={gene}"
        png_bytes = _fetch_png_bytes(url)
        tool_result_content, display_msg = _image_tool_result(desc, png_bytes, url)
        return tool_result_content, display_msg

    if name == "spatial_transcriptomics":
        gene = tool_input["gene"]
        url = f"{base}:9025/genes/{gene}"
        desc = f"Spatial transcriptomics plot for gene={gene}"
        png_bytes = _fetch_png_bytes(url)
        tool_result_content, display_msg = _image_tool_result(desc, png_bytes, url)
        return tool_result_content, display_msg

    if name == "static_images":
        image_name = tool_input["name"]
        image_map = {
            "Sinoid (SAN)": "ACM_VCM_SAN_default_plots.png",
            "ACO (Atrial Cardioids)": "ACM_VCM_SAN_default_plots.png",
            "VCO (Ventricular Cardioids)": "ACM_VCM_SAN_default_plots.png",
            "SAN-PACO (SAN Paced Atrial Cardioids)": "SAN_PCO_default_plots.png",
            "Mini-heart": "mini_heart_default_plots.png",
            "multiomics": "Multiomics_default_plots.png",
            "spatial_transcriptomics": "Spatial_default_plots.png",
            "ACM_VCM_SAN": "ACM_VCM_SAN_default_plots.png",
            "SAN-PCO": "SAN_PCO_default_plots.png",
        }
        file_name = image_map.get(image_name, "ACM_VCM_SAN_default_plots.png")
        try:
            with open(os.path.join(BASE_DIR, "frontend", "src", "assets", "imgs", file_name), "rb") as f:
                png_bytes = f.read()
            b64 = binToBase64(png_bytes)
            tool_result_content = [
                {
                    "type": "image",
                    "source": {
                        "type": "base64",
                        "media_type": "image/png",
                        "data": b64,
                    },
                },
                {
                    "type": "text",
                    "text": f"Default overview image for '{image_name}' loaded successfully.",
                },
            ]
            # Use data URI for frontend here — static images have no external URL
            display_msg = {"type": "image", "content": "data:image/png;base64," + b64}
            return tool_result_content, display_msg
        except FileNotFoundError:
            return (
                f"Static image '{image_name}' not found on disk.",
                {"type": "text", "content": f"Static image {image_name} not found."},
            )

    if name == "glkb_ai_assistant":
        question = tool_input["question"]
        success, answer = glkb_chat(question)
        if not success:
            return (
                f"GLKB call failed. Error: {answer}",
                {"type": "text", "content": f"GLKB call failed.\n\nError:\n{answer}"},
            )
        return (
            f"GLKB answer: {answer}",
            {"type": "text", "content": f"Ask GLKB AI assistant: {question}\n\nAnswer:\n\n{answer}"},
        )

    return (f"Unknown tool '{name}'.", None)


# ---------------------------------------------------------------------------
# Parallel executor — runs all planned tool calls concurrently
# ---------------------------------------------------------------------------

def execute_plan_parallel(steps: List[dict]) -> List[Tuple[Union[str, list], Union[dict, None]]]:
    """Execute a list of planned tool calls concurrently.

    Each step is {"tool": name, "input": {...}}.
    Returns results in the same order as the input steps list.
    """
    if not steps:
        return []

    results: List[Union[Tuple, None]] = [None] * len(steps)

    with ThreadPoolExecutor(max_workers=min(8, len(steps))) as pool:
        futures = {
            pool.submit(execute_tool, s["tool"], s["input"]): i
            for i, s in enumerate(steps)
        }
        for future in as_completed(futures):
            idx = futures[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                tool_name = steps[idx]["tool"]
                print(f"execute_plan_parallel: tool '{tool_name}' raised: {e}")
                results[idx] = (f"Tool '{tool_name}' failed: {e}", None)

    return results  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Core Claude call: Planner → Parallel Executor → Synthesizer
# ---------------------------------------------------------------------------

def get_gpt_resp(history: list) -> Tuple[bool, str, list]:
    """Three-phase pipeline: Planner → Parallel Executor → Synthesizer.

    Phase 1 (Planner): A dedicated Claude call with tool_choice forced to
        create_plan. Claude outputs a structured list of tool calls.
    Phase 2 (Executor): All planned tool calls run concurrently via
        ThreadPoolExecutor. Images are fetched in parallel from R servers.
    Phase 3 (Synthesizer): Claude receives all tool results (including base64
        images) and generates the final natural-language response.

    Fallback: if the planner outputs steps=[], Phase 2 is skipped and Phase 3
        is a direct text-only call (no tools), e.g. for pure biology questions.

    Returns: (success, error_msg, frontend_messages)
    """
    try:
        user_q = _latest_user_text(history)
        MAX_QUERY_CHARS = 1200
        user_q = user_q[:MAX_QUERY_CHARS]

        # RAG — fetch paper evidence for both planner and synthesizer prompts
        hits = paper_search(
            user_q,
            k=10,
            max_per_section=2,
            context_window=1,
            context_top_n=2,
        )
        evidence_block = _format_paper_evidence(hits, user_q)
        paper_section = (
            "## 4. Paper Evidence (authoritative)\n"
            "You MUST use the excerpts below as the source of truth for questions about the paper.\n"
            "If the answer is not supported by the excerpts, say you cannot find it in the paper.\n"
            "Do not invent details.\n\n"
            f"{evidence_block}\n"
        )

        anthropic_msgs = _flat_to_anthropic(history)
        collected_messages: List[dict] = []

        # -------------------------------------------------------------------
        # Phase 1 — Planner
        # -------------------------------------------------------------------
        planner_system = PLANNER_SYSTEM_PROMPT + "\n\n" + paper_section

        planner_resp = client.messages.create(
            model="claude-sonnet-4-5-20250929",
            temperature=0,
            messages=anthropic_msgs,
            max_tokens=1024,
            system=planner_system,
            tools=PLANNER_TOOL,
            tool_choice={"type": "tool", "name": "create_plan"},
        )

        # Extract the create_plan call — it is guaranteed by tool_choice
        plan_input: dict = {}
        for block in planner_resp.content:
            if hasattr(block, "name") and block.name == "create_plan":
                plan_input = block.input
                break

        steps: List[dict] = plan_input.get("steps", [])
        intro_text: str = plan_input.get("intro_text", "").strip()

        print(f"======== Planner: {len(steps)} step(s): {[s['tool'] for s in steps]}")
        log_queue.put(json.dumps({"planner": plan_input}, ensure_ascii=False))

        if intro_text:
            collected_messages.append({"type": "text", "content": intro_text})

        # -------------------------------------------------------------------
        # Phase 2 — Parallel Executor (skipped when steps is empty)
        # -------------------------------------------------------------------
        # tool_results_for_claude: list of tool_result blocks to send to synthesizer
        # Fake tool_use_ids are generated since we bypassed the normal tool_use flow.
        tool_results_for_claude: List[dict] = []

        if steps:
            for s in steps:
                print(f"  -> planned: {s['tool']}({s['input']})")

            exec_results = execute_plan_parallel(steps)

            for i, (result_content, display_msg) in enumerate(exec_results):
                if display_msg:
                    collected_messages.append(display_msg)
                # We need a stable fake tool_use_id for each step so we can
                # build a valid tool_result message for the synthesizer.
                fake_id = f"plan_step_{i}"
                tool_results_for_claude.append({
                    "type": "tool_result",
                    "tool_use_id": fake_id,
                    "content": result_content,
                })

        # -------------------------------------------------------------------
        # Phase 3 — Synthesizer
        # -------------------------------------------------------------------
        synthesizer_system = PROMPT + "\n\n" + paper_section

        # Build the synthesizer message history.
        # We inject a fake assistant turn that "called" the planned tools,
        # followed by a user turn with all the tool results, so Claude
        # understands what was executed and can interpret the results.
        synth_msgs = list(anthropic_msgs)  # copy

        if steps and tool_results_for_claude:
            # Fake assistant turn: one tool_use block per step
            fake_tool_uses = [
                {
                    "type": "tool_use",
                    "id": f"plan_step_{i}",
                    "name": s["tool"],
                    "input": s["input"],
                }
                for i, s in enumerate(steps)
            ]
            synth_msgs.append({"role": "assistant", "content": fake_tool_uses})
            synth_msgs.append({"role": "user", "content": tool_results_for_claude})

        synth_resp = client.messages.create(
            model="claude-sonnet-4-5-20250929",
            temperature=0.2,
            messages=synth_msgs,
            max_tokens=3072,
            system=synthesizer_system,
            # No tools passed — synthesizer only generates text
        )

        print(f"======== Synthesizer stop_reason={synth_resp.stop_reason}, "
              f"blocks=[{', '.join(f'text({len(b.text)} chars)' if hasattr(b, 'text') else 'other' for b in synth_resp.content)}]")

        for block in synth_resp.content:
            if hasattr(block, "text") and block.text.strip():
                collected_messages.append({"type": "text", "content": block.text})

        # Build flat assistant text for the history returned to the frontend
        assistant_text_parts = [m["content"] for m in collected_messages if m["type"] == "text"]
        assistant_text = "\n\n".join(assistant_text_parts)
        if assistant_text.strip():
            history.append({"role": "assistant", "content": assistant_text})

        print(f"======== Claude success. {len(collected_messages)} display messages.")
        return (True, "", collected_messages)

    except Exception:
        err = traceback.format_exc()
        print(err)
        return (False, "Failed to get the response from the AI. Please copy your input, refresh the page and try again.", [])


# ---------------------------------------------------------------------------
# GLKB integration  (unchanged)
# ---------------------------------------------------------------------------

def glkb_chat(question: str) -> Tuple[bool, str]:
    """
    Safe-ish GLKB SSE caller.
    - No history (messages = [])
    - Handles SSE chunking
    - Adds basic retries + backoff
    - Adds caps to avoid runaway memory
    - Better error messages (HTTP body snippet, content-type, etc.)
    """
    URL = GLKB_LLM_AGENT_URL
    PREFIX = "[AGENT OUTPUT] FinalAnswerAgent | Output:"

    CONNECT_TIMEOUT_S = 10
    READ_TIMEOUT_S = 180
    MAX_ATTEMPTS = 3
    BACKOFF_S = 1.5

    MAX_CHUNKS = 1000
    MAX_TOTAL_CHARS = 2_000_000

    payload = {"question": question, "messages": []}

    last_err = None

    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with requests.Session() as session:
                with session.post(
                    URL,
                    json=payload,
                    stream=True,
                    timeout=(CONNECT_TIMEOUT_S, READ_TIMEOUT_S),
                    verify=False,
                    headers={
                        "Accept": "text/event-stream",
                        "Content-Type": "application/json",
                        "User-Agent": "glkb-python-client/1.0",
                    },
                ) as r:
                    if r.status_code < 200 or r.status_code >= 300:
                        body_snip = ""
                        try:
                            body_snip = (r.text or "")[:800]
                        except Exception:
                            body_snip = "<unable to read body>"
                        ct = r.headers.get("content-type", "")
                        return (
                            False,
                            f"HTTP {r.status_code}. Content-Type: {ct}. Body (first 800 chars): {body_snip}",
                        )

                    ct = (r.headers.get("content-type") or "").lower()
                    if "text/event-stream" not in ct:
                        body_snip = ""
                        try:
                            body_snip = (r.text or "")[:800]
                        except Exception:
                            body_snip = "<unable to read body>"
                        return (
                            False,
                            f"Expected SSE (text/event-stream) but got Content-Type: {ct}. Body (first 800 chars): {body_snip}",
                        )

                    chunks: List[str] = []
                    total_chars = 0

                    for raw_line in r.iter_lines(decode_unicode=True):
                        if raw_line is None:
                            continue

                        line = raw_line.strip()
                        if not line:
                            continue
                        if not line.startswith("data:"):
                            continue

                        data_str = line[len("data:"):].strip()
                        if not data_str:
                            continue
                        if data_str == "[DONE]":
                            break

                        try:
                            obj = json.loads(data_str)
                            step = obj.get("step")

                            if step == "Complete":
                                final_resp = obj.get("response")
                                final_text = ""

                                if isinstance(final_resp, str) and final_resp.strip():
                                    final_text = final_resp.strip()

                                refs = obj.get("references")
                                if isinstance(refs, list) and refs:
                                    ref_lines = []
                                    for i, r in enumerate(refs[:10], start=1):
                                        if not isinstance(r, list) or len(r) < 6:
                                            continue

                                        title = r[0] or "Untitled"
                                        url = r[1] or ""
                                        year = r[3] or ""
                                        journal = r[4] or ""
                                        authors = r[5] or []

                                        if isinstance(authors, list):
                                            authors = [a for a in authors if a.strip()]
                                            author_str = ", ".join(authors[:5])
                                            if len(authors) > 5:
                                                author_str += " et al."
                                        else:
                                            author_str = ""

                                        line = f"{i}. {title}\n   {author_str} ({year}) — {journal}\n   {url}"
                                        ref_lines.append(line)

                                    if ref_lines:
                                        final_text += "\n\nReferences:\n" + "\n".join(ref_lines)

                                return True, final_text.strip()

                        except json.JSONDecodeError:
                            continue

                        content = obj.get("content")
                        if not isinstance(content, str) or not content:
                            continue

                        idx = content.find(PREFIX)
                        if idx == -1:
                            continue

                        piece = content[idx + len(PREFIX):].lstrip()
                        if not piece:
                            continue

                        if len(chunks) >= MAX_CHUNKS:
                            break
                        if total_chars + len(piece) > MAX_TOTAL_CHARS:
                            remaining = MAX_TOTAL_CHARS - total_chars
                            if remaining > 0:
                                chunks.append(piece[:remaining])
                                total_chars += remaining
                            break

                        chunks.append(piece)
                        total_chars += len(piece)

                    if not chunks:
                        return (
                            False,
                            "No FinalAnswerAgent output found in SSE stream. "
                            "The agent may have failed, returned only traces, or the prefix format changed.",
                        )

                    result = "\n".join(chunks).strip()
                    return True, result

        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            last_err = f"{type(e).__name__}: {e}"
            if attempt < MAX_ATTEMPTS:
                time.sleep(BACKOFF_S * attempt)
                continue
            return False, f"Network/timeout after {MAX_ATTEMPTS} attempts: {last_err}"

        except requests.exceptions.RequestException as e:
            last_err = f"{type(e).__name__}: {e}"
            return False, f"HTTP/network error: {last_err}"

        except Exception:
            last_err = traceback.format_exc()
            return False, last_err

    return False, last_err or "Unknown error"


# ---------------------------------------------------------------------------
# HTTP handler  (API contract unchanged)
# ---------------------------------------------------------------------------

def process_ai_chat(request, path: str):
    print('AI chat')
    MAX_BODY_BYTES = 256_000

    cl = int(request.headers.get("Content-Length", "0") or "0")
    if cl <= 0 or cl > MAX_BODY_BYTES:
        request.send_response(413)
        request.send_header("Content-Length", 0)
        request.send_header("Access-Control-Allow-Origin", "*")
        request.end_headers()
        return

    user_input = request.rfile.read(cl).decode("utf-8", errors="replace")
    if within_rate_limit() == False:
        request.send_response(429)
        request.send_header('Connection', 'keep-alive')
        request.send_header('Content-Length', 0)
        request.send_header('Access-Control-Allow-Origin', '*')
        request.end_headers()
        request.wfile.write(b'')
        request.wfile.flush()
        return

    history: list = json.loads(user_input)
    MAX_TURNS = 30
    if isinstance(history, list) and len(history) > MAX_TURNS:
        history = history[-MAX_TURNS:]

    MAX_MSG_CHARS = 4000

    def _trim_msg_content(m):
        c = m.get("content", "")
        if isinstance(c, str) and len(c) > MAX_MSG_CHARS:
            m = dict(m)
            m["content"] = c[:MAX_MSG_CHARS] + "…[truncated]"
        return m

    history = [_trim_msg_content(m) for m in history if isinstance(m, dict)]

    log_queue.put(json.dumps({'history': history}, ensure_ascii=False))
    success, error_msg, messages = get_gpt_resp(history)

    if error_msg:
        request.send_response(500)
        error_msg = error_msg.encode('utf-8')
        request.send_header('Content-Length', len(error_msg))
        request.send_header('Connection', 'keep-alive')
        request.send_header('Access-Control-Allow-Origin', '*')
        request.end_headers()
        request.wfile.write(error_msg)
        request.wfile.flush()
        return

    request.send_response(200)
    resp_data = json.dumps({'history': history, 'messages': messages}, ensure_ascii=False)
    resp_data = resp_data.encode('utf-8')
    request.send_header('Content-Length', len(resp_data))
    request.send_header('Connection', 'keep-alive')
    request.send_header('Access-Control-Allow-Origin', '*')
    request.end_headers()
    request.wfile.write(resp_data)
    request.wfile.flush()
    return


# ---------------------------------------------------------------------------
# Async log writer
# ---------------------------------------------------------------------------

log_queue = Queue()


def write_logs():
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(f"{time.time()*1000:.3f}: Server started!\n")
        f.flush()
        while True:
            item = log_queue.get()
            f.write(f"{time.time()*1000:.3f}: {item}\n")
            f.flush()

start_new_thread(write_logs, ())
