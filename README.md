# Triples with equal sum and equal reciprocal sum: code and data

This repository reproduces every computation stated in the paper

> Palaash Gang, Akshaj Agadi, Arjun Veluri, Jerry Wang, Jeremy Barreto, Aarin Chouthaiwale, *Triples with equal sum and equal reciprocal sum*.

The paper studies unordered triples of positive integers with $1/a+1/b+1/c<1$ that share both their sum
and their reciprocal sum with another such triple, through the cubic pencil
$(x+y+z)(xy+yz+zx)=\Lambda xyz$. The repository contains:

- `theory/diophantine/`: the exact enumeration of all such classes with sum at most 4800
  (`enumerate_fast.py` with the C kernel `enumerate_core.c`), its data (`data/groups.csv`,
  `data/per_S.csv`), the checks of the curves (`variety_checks.py`), the families and dual pairs
  (`families.py`), the growth fits (`growth_fits.py`, `exponent_fits.py`) and the Mordell–Weil ranks
  (`ranks.py`, PARI/GP);
- `checks/classes4800/`: a second, independent enumeration to $S=4800$ (`enum_classes.c`) with the checks
  that rest on it (`check_enum.py`), the group-law checks (`check_groups.py`), the 2-isogeny descents
  (`check_descent.py`), the isosceles family (`check_isosceles.py`) and the conic-bundle family
  (`check_theorem5_family.py`);
- `checks/enumeration6000/`: a third enumeration, to $S=6000$ (`enum.c`, `check_enum.py`);
- `checks/note/`: the remaining checks of the paper (`check_note.py`, with the enumeration of
  $4801\le S\le6000$), and `count_sums.py`;
- `checks/point3P/`: the point $3P$ on $C_{155/12}$;
- `code/`: the reference enumerator for $S\le600$ that `enumerate_fast.py` is checked against;
- `figures/`: Figure 1 (`figures/src/F9.py`) with its style files;
- `tools/make_tables.py`: Tables 1 and 2 (`tables/`).

42 data files, 11.42 MB in all (`DATA.md`).

## Install and run

Requirements: Python 3.13 (3.9 or newer works for the checks; the byte-identical figure rebuilds need the
pinned versions), `make`, a C compiler (`cc`) and, for the PARI steps, `cypari2`.

```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/pip install --only-binary :all: -r requirements-pari.txt
make verify PYTHON=.venv/bin/python      # fast exact checks
make all PYTHON=.venv/bin/python         # everything
```

| target | what it runs | time |
|---|---|---|
| `make verify` | the fast exact checks | about 8 min |
| `make full` | the long computations | about 8 min |
| `make figures` | rebuild every figure from its data and compare with the committed file | about 1 min |
| `make tables` | regenerate every table and compare with the committed file | about 1 min |
| `make data-check` | compare every data file with `SHA256SUMS` | seconds |
| `make all` | all of the above | about 17 min |
| `make clean` | remove caches and compiled helpers | |

Each target prints PASS, FAIL or SKIP for every step and exits nonzero if any step fails. A step also
fails if it changes a committed data file: regenerating committed data is a reproduction test. Times were
measured on an Apple M-series laptop with 8 GB of memory. `python3 tools/steps.py --list` lists every step;
`TRACE.md` maps every checkable item of the paper to its script and target; `DATA.md` lists every data
file with its size, checksum, producing script and paper element.

## Optional components

- **PARI/GP** (`requirements-pari.txt`, cypari2 2.2.4 with PARI/GP 2.17.2): the ranks and the descents.
  Without it those two steps are skipped with a message. Set `PYTHON_PARI` to use another interpreter.
- **C compiler**: the enumerations. Without `cc` (or `CC`) those steps are skipped with a message.
- `make full` reruns the three enumerations from scratch and compares their output with the committed
  files; they use all cores and about 400 MB of memory.


## Licences

- **Code** (`*.py`, `*.c`, `*.sh`, `Makefile`, `tools/`): MIT licence, `LICENSE`.
- **Data, figures and text** (every file listed in `DATA.md`, and the documentation): Creative Commons
  Attribution 4.0 International, `LICENSE-DATA`.

Copyright (c) 2026 the authors: Palaash Gang, Akshaj Agadi, Arjun Veluri, Jerry Wang, Jeremy Barreto, Aarin Chouthaiwale.

## How to cite

Cite the paper and the archived version of this repository:

- Palaash Gang, Akshaj Agadi, Arjun Veluri, Jerry Wang, Jeremy Barreto, Aarin Chouthaiwale, *Triples with equal sum and equal reciprocal sum*, Preprint: to appear (arXiv).
- Palaash Gang, Akshaj Agadi, Arjun Veluri, Jerry Wang, Jeremy Barreto, Aarin Chouthaiwale, *Code and data for: Triples with equal sum and equal reciprocal sum*. The archived version of this repository, with its DOI, is on Zenodo (see the DOI badge and link on the GitHub page).

`CITATION.cff` and `.zenodo.json` hold the same metadata, with the authors' ORCID iDs and affiliations.
The repository is https://github.com/pa1aash/equal-sum-reciprocal-triples.
