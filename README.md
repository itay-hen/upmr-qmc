# UPMR-QMC

Reference implementation of the quantum Monte Carlo algorithm of

> I. Hen, *A Quantum Monte Carlo Algorithm for Arbitrary Physical Systems* (UPMR: universal permutation matrix representation).

The Hamiltonian is written as H = D_0 + Σ_t P_t D_t (translations P_t times diagonal amplitudes
D_t evaluated on the departure state), and the partition function is sampled as a sum over closed
walks whose weights are products of amplitudes times exponential divided differences of the
visited energies. No time discretization is involved.

## Contents

| file | what it is |
|---|---|
| `upmr_gen.cpp` | **model-agnostic sampler.** Reads a Hamiltonian as operator strings and derives the decomposition, inverse labels, and subcycle templates automatically. Supports fermions (Jordan–Wigner signs on the fly), spins of any S, and bosons with a cutoff, mixed freely; complex amplitudes; canonical (fixed-N) or grand-canonical sampling. |
| `gen_models.py` | model front end: writes spec files for `upmr_gen` and computes an independent ED reference from explicit sparse matrices. Contains the spin–boson, three-label, and fermion–spin–boson models of the paper. |
| `upmr_qmc.cpp` | earlier hand-coded sampler for the spin–boson benchmarks (Sec. "Benchmark application"). |
| `ed_ref.py`, `ed_lanczos.py` | ED references for the spin–boson tables (full ED; Lanczos in parity sectors for N_s = 4). |
| `make_tables.py`, `make_tables_mixed.py` | regenerate the LaTeX tables of the paper from the raw outputs. |
| `jobs.txt`, `jobsD.txt`, `gen/jobs.txt` | exact command lines behind every table. |
| `runs/`, `gen/runs/`, `ed_ref.out`, `gen/ed_mixed.out` | raw outputs used in the paper. |

## Build and quick check

    make            # builds upmr_gen and upmr_qmc (C++17, no dependencies)
    make check      # spin-boson pair at beta=5 vs ED, about a minute

The Python scripts need `numpy` and `scipy`.

## Defining a model

In `gen_models.py`, operators are built from elementary ones with `el(op, component)` and
combined with `+`, `*`, scalars and `.dag()`:

| species | component | elementary operators |
|---|---|---|
| fermion mode | `('F',)` | `cd`, `c`, `n` (Jordan–Wigner order = order of the F components) |
| spin S | `('S', 2S+1)` | `Sp`, `Sm`, `Sz` (local value v = 0..2S, S^z = v − S) |
| boson | `('B', nmax)` | `bd`, `b`, `nb` |

A model is a dict with `comps`, `terms` (list of `(tag, operator)`), `obs` (diagonal observables),
`root` (initial basis state), `N` (fermion number for canonical sampling, or `None`), and
`offdiag_tags`. Then

    python3 gen_models.py mixed --kw '{"pbc":false}' --betas 1,2,4 --spec gen/open.spec   # spec + ED
    ./upmr_gen --spec gen/open.spec --beta 2 --sub 1 --nonlocal 1 --sweeps 2e6 --therm 5e4 --seed 102

### Spec format (what `upmr_gen` reads)

    comp F | comp S <d> | comp B <nmax>             one line per component, in order
    canonical 0|1                                    fix the fermion number of the root
    root v0 v1 ...                                   initial basis state
    term <tag> <re> <im> <k> op c op c ...           operator string (ops left to right)
    obs <name> <re> <k> op c ...                     diagonal observable

### Sampler options

`--beta --sweeps --therm --bins --seed`; moves: `--sub 1` (subcycle replacement with automatically
generated templates, `--tmax` = max template length, default 4), `--nonlocal 1` (re-rooting,
inner-loop reversal, segment exchange), `--worm 1 --mu <mu>`. Basis-state shifts, swaps, and pair
insertions/deletions are always on.

Output: average phase, and for each quantity value and jackknife error — energy `E`, expansion
order `q` (⟨q⟩ = −β⟨H_off⟩), the off-diagonal energy of each term tag `H_<tag>` (from vertex
counts), and each diagonal observable (root estimator `<name>_root` and ω-weighted `<name>`).
Every run also checks the Hermitian pairing d_t̄(π_t z) = conj d_t(z) on random states and
recomputes the weight of every weight-preserving nonlocal move (`weight_check_fail` must be 0).

## Reproducing the paper's tables

    bash -c 'while read n a; do ./upmr_gen $a > gen/runs/$n.out; done < gen/jobs.txt'   # fermion-spin-boson
    bash -c 'while read n a; do ./upmr_qmc $a > runs/$n.out; done < jobs.txt'          # spin-boson
    python3 make_tables_mixed.py; python3 make_tables.py

Timings in the raw outputs were measured on a shared, heavily loaded workstation.
