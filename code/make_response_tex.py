#!/usr/bin/env python3
"""Generate response/response_to_reviewers.tex from reviewer_comments.txt.

Comments are transcribed verbatim from the decision letter so nothing is
paraphrased; each gets an empty \\response block for the authors to fill in.
Regenerating overwrites the file, so write responses in the .tex, not here.
"""
import re, os

SRC = "reviewer_comments.txt"
OUT = "../response/response_to_reviewers.tex"

UNI = {
    "\u2019": "'", "\u2018": "`", "\u201c": "``", "\u201d": "''",
    "\u2014": "---", "\u2013": "--", "\u2026": r"\dots{}",
    "\u00d7": r"$\times$", "\u03b1": r"$\alpha$", "\u00a0": " ",
}


def tex_escape(s):
    # Escape TeX specials on the RAW text first. Doing the unicode -> LaTeX
    # mapping first would get the inserted commands escaped in turn
    # (\dots{} -> \textbackslash{}dots{}), which is how "......" once ended up
    # rendering as literal "\{}dots{}".
    s = s.replace("\\", r"\textbackslash ")
    for ch in "&%$#_{}":
        s = s.replace(ch, "\\" + ch)
    s = s.replace("~", r"\textasciitilde ").replace("^", r"\textasciicircum ")
    # Only now substitute unicode, so these commands survive intact.
    for k, v in UNI.items():
        s = s.replace(k, v)
    return s


lines = open(SRC, encoding="utf8").read().split("\n")

# ---- split into reviewer blocks -------------------------------------------
rev_starts = [i for i, l in enumerate(lines) if re.match(r"^Reviewer: \d", l)]
end = next((i for i, l in enumerate(lines) if l.startswith("Editor's Comments")), len(lines))
blocks = []
for n, s in enumerate(rev_starts):
    e = rev_starts[n + 1] if n + 1 < len(rev_starts) else end
    blocks.append((lines[s].split(":")[1].strip(), lines[s + 1:e]))

HEAD = re.compile(r"^(Major|Minor)\s+(Comments|Concerns):?\s*$", re.I)
ITEM = re.compile(r"^(\d+)\.\s+(.*)$")


def parse(body):
    """-> (intro paragraphs, [(section, [items])])"""
    intro, sections, cur, item = [], [], None, None

    def flush():
        nonlocal item
        if item is not None:
            cur[1].append("\n\n".join(item).strip())
            item = None

    for ln in body:
        t = ln.strip()
        if not t or t == "Comments to the Author":
            if item is not None:
                item.append("")           # paragraph break inside an item
            continue
        h = HEAD.match(t)
        if h:
            flush()
            cur = (h.group(1).capitalize(), [])
            sections.append(cur)
            continue
        m = ITEM.match(t)
        if m:
            flush()
            if cur is None:               # reviewer with no Major/Minor header
                cur = ("", [])
                sections.append(cur)
            item = [m.group(2)]
            continue
        if item is not None:
            item.append(t)
        else:
            intro.append(t)
    flush()
    # collapse blank-line artefacts inside items
    sections = [(n, [re.sub(r"\n{3,}", "\n\n", i) for i in its]) for n, its in sections]
    return intro, sections


PREAMBLE = r"""% !TeX program = pdflatex
% ---------------------------------------------------------------
% Response to Reviewers -- BIB-26-1707, Briefings in Bioinformatics
% Reviewer comments transcribed verbatim from the decision letter
% of 01-Sep-2026. Write your replies in the \\response blocks.
% ---------------------------------------------------------------
\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[margin=1in]{geometry}
\usepackage{parskip}

% Overleaf reuses ONE output.aux for every document in the project. If the
% manuscript (which loads lineno) was compiled first, that aux holds hundreds
% of \@LN{..}{..} entries; this file would then die on them with no PDF.
% These no-ops make such a leftover aux harmless.
\makeatletter
\providecommand\@LN[2]{}
\providecommand\@LN@col[1]{}
\makeatother

\newenvironment{comment}{\par\medskip\begin{quote}\itshape}{\end{quote}}
\newenvironment{response}{\par\noindent\textbf{Response:}\quad}{\par\medskip}

\title{Response to Reviewers \\[0.3em]
  \large Manuscript BIB-26-1707, \textit{Briefings in Bioinformatics}}
\author{}
\date{}

\begin{document}
\maketitle

\noindent Dear Dr.~Ma, Dr.~Chen, and the reviewers,

\noindent We thank the reviewers for their careful reading of our manuscript.
% TODO: summarise the main changes here.
Below we respond to each comment in turn.

\noindent Sincerely,\\
The Authors

"""

parts = [PREAMBLE]
for num, body in blocks:
    intro, sections = parse(body)
    parts.append("\n\\section*{Reviewer %s}\n\\addcontentsline{toc}{section}{Reviewer %s}\n"
                 % (num, num))
    if intro:
        parts.append("\\begin{comment}\n%s\n\\end{comment}\n"
                     % tex_escape("\n\n".join(intro)))
        parts.append("\\begin{response}\n%% TODO\n\\end{response}\n")
    for sec_name, items in sections:
        if sec_name:
            parts.append("\n\\subsection*{%s Comments}\n" % sec_name)
        for i, it in enumerate(items, 1):
            # "Major"/"Minor" both start with M -- use an explicit prefix so
            # Major 1 and Minor 1 do not collide.
            prefix = {"Major": "", "Minor": "m"}.get(sec_name, "")
            tag = "R%s.%s%d" % (num, prefix, i)
            parts.append("\n\\noindent\\textbf{Comment %s}\n" % tag)
            parts.append("\\begin{comment}\n%s\n\\end{comment}\n" % tex_escape(it))
            parts.append("\\begin{response}\n%% TODO\n\\end{response}\n")

parts.append("\n\\end{document}\n")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

# Refuse to clobber written responses. Once anything has been typed into a
# \response block, this file is the source of truth, not the generator.
if os.path.exists(OUT):
    existing = open(OUT, encoding="utf8").read()
    blank = existing.count("\\begin{response}\n%% TODO") + existing.count("\\begin{response}\n%%%% TODO")
    written = existing.count("\\begin{response}") - blank
    if written > 0:
        raise SystemExit(
            "refusing to overwrite %s: %d response block(s) already contain text.\n"
            "Edit the .tex directly, or delete it first if you really want a blank one."
            % (OUT, written))

open(OUT, "w", encoding="utf8").write("".join(parts))

n_items = sum(len(i) for _, b in blocks for _, i in parse(b)[1])
print("wrote %s  (%d reviewers, %d comments)" % (OUT, len(blocks), n_items))
