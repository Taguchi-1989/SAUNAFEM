"""Convert governing_equations.md (Markdown + LaTeX math) to PDF via xelatex.

Usage:
    python scripts/md2pdf.py                          # default: docs/governing_equations.md
    python scripts/md2pdf.py docs/governing_equations.md -o output.pdf

Requirements:
    - MiKTeX (xelatex) installed
    - Python 3.11+
    - pip install markdown
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def md_to_latex(md_text: str) -> str:
    """Convert Markdown+LaTeX-math text to a full LaTeX document."""

    lines = md_text.split("\n")
    latex_lines: list[str] = []
    in_code_block = False
    in_table = False
    in_math_block = False
    skip_toc_section = False
    table_cols = 0
    table_header_pending = True

    def _close_table() -> None:
        nonlocal in_table, table_header_pending
        if in_table:
            latex_lines.append(r"\end{tabularx}")
            latex_lines.append(r"}")  # close \small
            latex_lines.append("")
            in_table = False
            table_header_pending = True

    for line in lines:
        if line.strip().startswith("```"):
            if in_code_block:
                latex_lines.append(r"\end{verbatim}")
                in_code_block = False
            else:
                _close_table()
                latex_lines.append(r"\begin{verbatim}")
                in_code_block = True
            continue

        if in_code_block:
            latex_lines.append(line)
            continue

        # Display math blocks ($$...$$) — pass through without escaping
        if line.strip() == "$$":
            if in_math_block:
                latex_lines.append("$")
                in_math_block = False
            else:
                latex_lines.append("$")
                in_math_block = True
            continue

        if in_math_block:
            latex_lines.append(line)  # raw LaTeX math — no escaping
            continue

        # Skip the manual TOC section (LaTeX generates its own)
        # Detect "## 目次" heading
        m_toc = re.match(r"^##\s+目次", line)
        if m_toc:
            skip_toc_section = True
            continue
        if skip_toc_section:
            # Skip until next ## heading or ---
            if re.match(r"^##\s+", line) or re.match(r"^-{3,}\s*$", line):
                skip_toc_section = False
                # Fall through to process this line
            else:
                continue

        # Headings
        m = re.match(r"^(#{1,6})\s+(.*)", line)
        if m:
            _close_table()
            level = len(m.group(1))
            title = _escape_tex(_inline_math_preserve(m.group(2)))
            cmds = [
                r"\section",
                r"\subsection",
                r"\subsubsection",
                r"\paragraph",
                r"\subparagraph",
                r"\subparagraph",
            ]
            latex_lines.append(f"{cmds[level - 1]}{{{title}}}")
            continue

        # Horizontal rule
        if re.match(r"^-{3,}\s*$", line):
            _close_table()
            latex_lines.append(r"\bigskip\hrule\bigskip")
            continue

        # Table rows — only detect if line starts/ends with |
        # (avoids false positive on inline math like $|x|$)
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|") and not stripped.startswith("```"):
            cells = [c.strip() for c in stripped.strip("|").split("|")]
            # Skip separator rows like |---|---|
            if all(re.match(r"^[-:]+$", c) for c in cells):
                continue
            if not in_table:
                table_cols = len(cells)
                # Use tabularx: first column fixed-width, rest auto-wrap
                if table_cols <= 3:
                    col_spec = "|".join(["L"] * table_cols)
                else:
                    # First column narrow, rest flexible
                    col_spec = "l|" + "|".join(["L"] * (table_cols - 1))
                latex_lines.append(r"{\small")
                latex_lines.append(r"\noindent")
                latex_lines.append(
                    r"\begin{tabularx}{\textwidth}{|" + col_spec + r"|}"
                )
                latex_lines.append(r"\hline")
                in_table = True
                table_header_pending = True
            escaped = [_escape_tex(_inline_math_preserve(c)) for c in cells]
            # Bold the header row (first row of each table)
            if table_header_pending:
                escaped = [r"\textbf{" + e + "}" for e in escaped]
                table_header_pending = False
            latex_lines.append(" & ".join(escaped) + r" \\")
            latex_lines.append(r"\hline")
            continue

        # Close table if we leave it
        if in_table and line.strip() == "":
            _close_table()
            latex_lines.append("")
            continue

        # Bullet lists
        m_bullet = re.match(r"^(\s*)-\s+(.*)", line)
        if m_bullet:
            content = _escape_tex(_inline_math_preserve(m_bullet.group(2)))
            latex_lines.append(rf"\begin{{itemize}} \item {content} \end{{itemize}}")
            continue

        # Numbered lists
        m_num = re.match(r"^(\s*)\d+\.\s+(.*)", line)
        if m_num:
            content = _escape_tex(_inline_math_preserve(m_num.group(2)))
            latex_lines.append(
                rf"\begin{{enumerate}} \item {content} \end{{enumerate}}"
            )
            continue

        # Empty line = paragraph break
        if line.strip() == "":
            latex_lines.append("")
            continue

        # Regular paragraph text
        latex_lines.append(_escape_tex(_inline_math_preserve(line)))

    _close_table()

    body = "\n".join(latex_lines)

    return _wrap_document(body)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_MATH_PLACEHOLDER = "\x00MATH"


def _inline_math_preserve(text: str) -> str:
    """Replace inline $...$ with placeholders so they survive escaping."""
    parts = []
    i = 0
    math_segments: list[str] = []
    while i < len(text):
        if text[i] == "$":
            j = text.find("$", i + 1)
            if j == -1:
                parts.append(text[i])
                i += 1
            else:
                placeholder = f"{_MATH_PLACEHOLDER}{len(math_segments)}\x00"
                math_segments.append(text[i : j + 1])
                parts.append(placeholder)
                i = j + 1
        elif text[i] == "`":
            j = text.find("`", i + 1)
            if j == -1:
                parts.append(text[i])
                i += 1
            else:
                code = text[i + 1 : j]
                placeholder = f"{_MATH_PLACEHOLDER}{len(math_segments)}\x00"
                math_segments.append(r"\texttt{" + _escape_tex_raw(code) + "}")
                parts.append(placeholder)
                i = j + 1
        else:
            parts.append(text[i])
            i += 1
    result = "".join(parts)
    # Restore math segments
    for idx, seg in enumerate(math_segments):
        result = result.replace(f"{_MATH_PLACEHOLDER}{idx}\x00", seg)
    return result


def _escape_tex_raw(text: str) -> str:
    """Escape TeX special characters (no math awareness)."""
    for ch, repl in [
        ("\\", r"\textbackslash{}"),
        ("{", r"\{"),
        ("}", r"\}"),
        ("&", r"\&"),
        ("%", r"\%"),
        ("#", r"\#"),
        ("_", r"\_"),
        ("~", r"\textasciitilde{}"),
        ("^", r"\textasciicircum{}"),
    ]:
        text = text.replace(ch, repl)
    return text


def _escape_tex(text: str) -> str:
    """Escape TeX special chars while preserving $...$ math and commands."""
    # Convert markdown links [text](url) → text only
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)

    # Split on $ math delimiters and \texttt{...} blocks — preserve both
    parts = re.split(r"(\$[^$]+\$|\\texttt\{[^}]*\})", text)
    result = []
    for part in parts:
        if part.startswith("$") and part.endswith("$"):
            result.append(part)
        elif part.startswith(r"\texttt{"):
            result.append(part)
        else:
            # Bold **text**
            part = re.sub(
                r"\*\*(.+?)\*\*", lambda m: r"\textbf{" + m.group(1) + "}", part
            )
            # Escape special chars (but not backslash commands we just inserted)
            out = []
            i = 0
            while i < len(part):
                if part[i] == "\\" and i + 1 < len(part) and part[i + 1].isalpha():
                    # LaTeX command — pass through until next non-alpha
                    j = i + 1
                    while j < len(part) and (part[j].isalpha() or part[j] in "{}"):
                        j += 1
                    out.append(part[i:j])
                    i = j
                elif part[i] in "&%#_":
                    out.append("\\" + part[i])
                    i += 1
                elif part[i] == "~":
                    out.append(r"\textasciitilde{}")
                    i += 1
                elif part[i] == "^":
                    out.append(r"\textasciicircum{}")
                    i += 1
                else:
                    out.append(part[i])
                    i += 1
            result.append("".join(out))
    return "".join(result)


def _wrap_document(body: str) -> str:
    """Wrap body in a full LaTeX document with CJK support."""
    return r"""\documentclass[a4paper,11pt]{article}

% --- CJK support (xelatex) ---
\usepackage{fontspec}
\usepackage{xeCJK}
\setCJKmainfont{Yu Mincho}
\setCJKsansfont{Yu Gothic}
\setCJKmonofont{Yu Gothic}

% --- Packages ---
\usepackage{amsmath,amssymb}
\usepackage{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage{fancyhdr}
\usepackage{xcolor}
\usepackage{enumitem}
\usepackage{tabularx}
\usepackage{array}

% Column type L: auto-wrapping left-aligned column for tabularx
\newcolumntype{L}{>{\raggedright\arraybackslash}X}

\geometry{margin=2.0cm}
\setlength{\headheight}{14pt}
\hypersetup{colorlinks=true, linkcolor=blue, urlcolor=blue}

\pagestyle{fancy}
\fancyhf{}
\rhead{SaunaFlow}
\lhead{支配方程式ドキュメント}
\cfoot{\thepage}

\title{\textbf{SaunaFlow 支配方程式ドキュメント}}
\author{SaunaFlow Project}
\date{\today}

\begin{document}
\maketitle
\tableofcontents
\newpage

""" + body + r"""

\end{document}
"""


def build_pdf(tex_path: Path, output_dir: Path) -> Path:
    """Run xelatex twice (for TOC) and return the PDF path."""
    for i in range(2):
        result = subprocess.run(
            [
                "xelatex",
                "-interaction=nonstopmode",
                "-output-directory",
                str(output_dir),
                str(tex_path),
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if result.returncode != 0 and i == 1:
            print("xelatex STDOUT:", result.stdout[-2000:], file=sys.stderr)
            print("xelatex STDERR:", result.stderr[-2000:], file=sys.stderr)
            sys.exit(1)

    pdf_name = tex_path.stem + ".pdf"
    return output_dir / pdf_name


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert Markdown+LaTeX to PDF")
    parser.add_argument(
        "input",
        nargs="?",
        default="docs/governing_equations.md",
        help="Input markdown file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output PDF path (default: same name as input, in docs/)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: {input_path} not found", file=sys.stderr)
        sys.exit(1)

    md_text = input_path.read_text(encoding="utf-8")
    latex_src = md_to_latex(md_text)

    # Determine output
    if args.output:
        output_pdf = Path(args.output)
    else:
        output_pdf = input_path.with_suffix(".pdf")

    # Build in temp dir, copy result
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        tex_file = tmp / "document.tex"
        tex_file.write_text(latex_src, encoding="utf-8")

        # Also save .tex next to the .md for reference
        tex_copy = input_path.with_suffix(".tex")
        tex_copy.write_text(latex_src, encoding="utf-8")
        print(f"LaTeX source saved: {tex_copy}")

        pdf_path = build_pdf(tex_file, tmp)

        shutil.copy2(pdf_path, output_pdf)
        print(f"PDF generated: {output_pdf}")


if __name__ == "__main__":
    main()
