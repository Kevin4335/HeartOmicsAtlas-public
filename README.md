<h1 align="center">HeartOmicsAtlas</h1>

<p align="center">
  <strong>An AI-powered, open-access web portal for exploring the human sinoatrial node.</strong>
</p>

<p align="center">
  Search genes · View expression across modalities · Ask the AI
</p>

---

## Info

HeartOmicsAtlas is the public codebase for a **web portal** that lets users explore single-cell and spatial omics data from the human fetal heart and related organoid models. It accompanies the manuscript *Unraveling Human Sinoatrial Node Development Using Fetal Heart and SAN-Paced Mini-Heart Models*.

| You can… | How |
|----------|-----|
| **Explore** | Search by gene and view expression (UMAPs, violin/dot plots, spatial maps, IGV) |
| **Compare** | Switch between scMultiomics, spatial transcriptomics, and scRNA-seq (ACM/VCM/SAN, SAN-PCO, Mini-heart) |
| **Ask** | Use the built-in AI assistant (Claude or OpenAI) to query the data and get plot suggestions |

---

## Data in the atlas

- **scMultiomics**
- **Spatial transcriptomics**
- **scRNA-seq** (ACM/VCM/SAN, SAN-PCO, Mini-heart)

---

## Tech stack

| Layer | Stack |
|-------|--------|
| **Frontend** | React, TypeScript, Vite, MUI, React Router, react-zoom-pan-pinch, react-markdown |
| **Backend** | Python (HTTP server), R (plot servers per modality) |
| **AI** | **Claude (Anthropic)** or **OpenAI** — configurable via `config.py` (`API_KEY` / `CHAT_KEY`). Optional RAG over the manuscript. |

---

## Quick start

### 1. R plot servers (ports 9025–9029)

From repo root, either use the helper script (set `BASE_DIR` in the script first):

```bash
./utils/restart_r_servers.sh
```

Or start each server manually:

```bash
cd resources-NEW/spatial_data     && Rscript Spatial_plot_function_new.R &
cd resources-NEW/multi_omics_data && Rscript Multiomics_plot_function_new.R &
cd resources-NEW/SAN_ACM_VCM      && Rscript ACM_VCM_SAN_plot_function_new.R &
cd resources-NEW/SAN-PCO          && Rscript SAN_PCO_plot_function_new.R &
cd resources-NEW/Mini-heart       && Rscript mini_heart_plot_function_new.R &
```

| Server      | Port |
|------------|------|
| Spatial    | 9025 |
| Multiomics | 9026 |
| ACM_VCM_SAN| 9027 |
| SAN-PCO    | 9028 |
| Mini-heart | 9029 |

### 2. Main web server

```bash
cd frontend && npm install && npm run build && cd ..
python server.py
```

Then open **http://localhost:8000**.  
For the AI chat, add a `config.py` with `API_KEY` and/or `CHAT_KEY` (Anthropic/OpenAI).

---

## Project layout

```
├── frontend/           React app (Explore Atlas, AI chat UI)
├── resources-NEW/      R plot servers (one per modality)
├── paper_rag/          Optional RAG for manuscript-backed AI answers
├── ai.py               AI chat backend (Claude/OpenAI, tools, GLKB)
├── server.py           Main HTTP server (frontend + /chat API)
└── utils/              restart_r_servers.sh, etc.
```

---

## Development

**Frontend only** (dev server, typically port 5173):

```bash
cd frontend && npm install && npm run dev
```

For full behavior (gene plots, AI), run the R servers and `python server.py` as above.

---

## Contact

**Chen Lab**, Weill Cornell Medicine — [shc2034@med.cornell.edu](mailto:shc2034@med.cornell.edu)
