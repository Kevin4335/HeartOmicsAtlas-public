from docx import Document
from typing import List, Tuple

def extract_docx_with_headings(docx_path: str) -> List[Tuple[str, str]]:
    """
    Returns [("heading"|"para", text), ...].
    Headings are detected primarily via Word paragraph styles.
    """
    doc = Document(docx_path)
    out: List[Tuple[str, str]] = []

    for p in doc.paragraphs:
        txt = (p.text or "").strip()
        if not txt:
            continue

        style_name = ""
        try:
            style_name = (p.style.name or "").lower()
        except Exception:
            style_name = ""

        # This catches Heading 1/2/3/... and any style containing "heading"
        is_heading = ("heading" in style_name) or (style_name == "title")

        out.append(("heading" if is_heading else "para", txt))

    return out

def to_marked_text(items: List[Tuple[str, str]]) -> str:
    """
    Convert into plain text where:
    - Every heading is on its own line and prefixed with ##HEADING##
    - Every element is separated by a blank line (critical for paragraph splitting)
    """
    lines: List[str] = []
    for kind, txt in items:
        txt = (txt or "").strip()
        if not txt:
            continue

        if kind == "heading":
            lines.append(f"##HEADING## {txt}")
        else:
            lines.append(txt)

        lines.append("")  # blank line boundary

    return "\n".join(lines).strip()
