# Source and computational certificates

This supplement contains the LaTeX source and bibliography for *Where the
quantum capacity sits in the arithmetical hierarchy*, plus exact finite
graph certificates and their replay programs. Paths below are relative
to the repository root. The paper uses `paper/abstract.txt` directly,
so its abstract and the submission text have one source.

## Build the paper

From `paper/`, run these commands in order with a TeX installation that
provides the article class, the packages in the preamble and plainurl.bst:

```
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Replay exact finite certificates

Use Python 3 with SymPy, networkx and psutil. Install python-sat only if
regenerating SAT formulas and proofs. Run one command at a time from
`certificates/`:

```
python -B bounded_run.py finite python -B replay.py finite
python -B bounded_run.py wrapper python -B check_wrapper.py
```

The finite replay reconstructs the 24 vertices and adjacency directly
from the parity equations. It checks the all-field rational dimension
dual, induced Mobius-ladder core, exact Pauli projections, clique bounds,
the 29-word witness and all 1,152 symmetry maps. It verifies the two
double cosets of empty-block patterns, of sizes 72 and 648, and the
question-fixing answer-transitive kernel used in the CNF reduction.
`dimension-dual.json` and `lower29.json` contain the exact witness data.
The wrapper check uses exact SymPy arithmetic on all 16 matrix units.

## Check the saved UNSAT proofs

The two DIMACS CNFs each have 11,496 variables and 50,449 clauses. They
encode the 28,512 strong-product conflicts, the representative empty
blocks, exactly five words in every first/second-question slice and one
symmetry-fixed occupied word. The saved Lingeling traces have 91,682
and 213,934 lines. Checking these files does not require a SAT solver.

`drat-trim.c` is the upstream checker source from
https://github.com/marijnheule/drat-trim, with its license intact.
Compile it from `certificates/` on a C toolchain with sys/time.h:

```
gcc -O2 -o drat-trim drat-trim.c
```

On a Windows GNU C toolchain lacking getc_unlocked, use:

```
gcc -O2 -Dgetc_unlocked=getc -o drat-trim.exe drat-trim.c
```

The substitution uses ordinary locked stdio and does not change the
checking algorithm. On a Unix shell, check the saved traces with:

```
python -B bounded_run.py proof-012345 ./drat-trim independent-upper30-012345.cnf independent-upper30-012345.drat
python -B bounded_run.py proof-013245 ./drat-trim independent-upper30-013245.cnf independent-upper30-013245.drat
```

On Windows, substitute `./drat-trim.exe` for `./drat-trim`. Require exit
status zero and `s VERIFIED` for both cases. SHA256SUMS.txt fixes the
released files; the two solver JSON files also record the CNF hashes.
The bound wrapper caps each process tree at 600 seconds and 1.4 billion
bytes of RSS, keeps numeric-library thread counts at one through the
replay programs, and writes its replay logs beside itself.

## Optional solver regeneration

Work in a copy of `certificates/`, because the following commands replace
the CNFs, traces and solver metadata with fresh outputs. With python-sat
and its Lingeling backend installed, run sequentially:

```
python -B bounded_run.py sat-012345 python -B replay.py 0
python -B bounded_run.py sat-013245 python -B replay.py 1
```

Then check both new traces as above. The replay explicitly saves and
flushes the trace before os._exit(0): this bypasses a Windows native
Lingeling cleanup failure observed after valid traces had been written.
The separate checker, rather than process exit alone, validates UNSAT.

## Verification scope

Both released traces have been checked with drat-trim. The graph
reduction, Boolean encoding and checker implementation remain in the
trust base. This is not a proof-assistant formalization. The checks
establish finite quantities only; they do not settle a Shannon-capacity
gap, a repetition conjecture or ordinary single-channel undecidability.
The computational proof files are preserved byte for byte. `replay.py`
changes only the original verifier's data paths and lower-witness file
layout so this archive can be replayed after extraction anywhere.
