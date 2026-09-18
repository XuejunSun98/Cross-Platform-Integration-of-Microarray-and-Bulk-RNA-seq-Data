#!/usr/bin/env python
"""Run Rank-In natively on Linux.

Rank-In is distributed only as a Windows executable (Rank-In.rar from
http://www.badd-cao.net/rank-in/), but that executable is a PyInstaller bundle
of a pure-Python 3.9 script (pandas/numpy/scipy only -- nothing Windows
specific). This runner executes the vendor's own code object straight out of
the bundle, so the code path is identical to the published tool: no
decompilation, no Wine, no reimplementation.

The program normally prompts for four paths interactively; we feed them on
stdin, which also replaces the processx stdin handshake used by the
integration_comparison_*Rank_in.R scripts.

Usage:
    python run_rankin.py <expr.txt> <class.txt> <out_result.txt> <out_deg.txt>

  expr.txt   tab-delimited, genes in rows, first column = gene id, header row
             = sample names  (e.g. simulation_2024/Rank_in/d_all_Rankin.txt)
  class.txt  tab-delimited, first column = sample_name, column "Class"
             (e.g. simulation_2024/Rank_in/sample.txt)
"""
import sys, os, marshal, io

PYC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "Rank-In.exe_extracted", "Rank-In.pyc")

if len(sys.argv) != 5:
    sys.exit(__doc__)

expr, cls, out_res, out_deg = sys.argv[1:5]

with open(PYC, "rb") as f:
    f.read(16)                      # CPython 3.9 .pyc header
    code = marshal.load(f)

sys.stdin = io.StringIO("\n".join([expr, cls, out_res, out_deg]) + "\n")

g = {"__name__": "__main__", "__builtins__": __builtins__, "__file__": "Rank-In.py"}
exec(code, g)
