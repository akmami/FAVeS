# Uniqueness analyses

The core is `match-length.cpp`: it computes the **exact** uniqueness curve `p(k)` for all `k` in one suffix-array pass, plus the maximal-match (repeat-length) distribution and base composition (`Q`, `k_c`).

The second important code is `sketch-seeds.cpp`; it sketches genomes and prints uniqueness and coverage information of sketches on various genomes with different relaxation radii.

## Quick start

```
cd analyses
make divsufsort                 # build libdivsufsort locally
make matchlen  GENOME=human     # exact p(k)=CDF(L) for one genome
make matchlenAll                # match-length for every genome in GENOMES
make plots                      # fit (alpha,beta) + render figures
```

`divsufsort` needs `cmake`, `make`, and `git`/`curl` (on clusters: `module load cmake git`). 
It builds a static libdivsufsort into `../third_party/libdivsufsort`; `match-length` auto-detects that prefix and falls back to a system install.

## Targets

| Target                | What it does                                                                            |
|-----------------------|-----------------------------------------------------------------------------------------|
| `make divsufsort`     | fetch + build libdivsufsort (static, 32+64-bit) locally                                 |
| `make matchlenAll`    | match-length for each genome                                                            |
| `make sketch-all`     | build + run technique for all specied over its grid -> `out/<m>/<g>.json`               |
| `make sweep`          | `kmer`/`minimizer`/`fracmh`/`syncmer`/`blend` over `KSWEEP` -> `out/sweep/<g>.<m>.json` |
| `make plots`          | fit `(α,β)`, render `uniqueness_curves.<s>.pdf` + `uniqueness_repeat_spectrum.<s>.pdf`  |

`<m>`: {`kmer fracmh minimizer syncmer blend strobemer strobemer2 lcp`}.
`<g>`: {`ecoliK12 mtub sacCer3 pf3D7 dm6 human`}.

## Analysis -> command

| # | Analysis                         | Command                                                                              |
|---|----------------------------------|--------------------------------------------------------------------------------------|
| 1 | Exact `p(k)=CDF(L)` identity     | `make matchlen` (compare `p_unique` to `make kmer`)                                  |
| 2 | Crossover `k_c = ln N / ln(1/Q)` | `k_c` field in the JSON vs where `p(k)` leaves the null                              |
| 3 | Plateau = repeat content         | `make plots` fits `1 − α·S_rep(k)`; compare `α` to RepeatMasker                      |
| 4 | Cross-species ordering           | `make matchlenAll GENOMES="ecoliK12 mtub sacCer3 pf3D7 dm6 human"` then `make plots` |
| 5 | Composition control              | `Q`, `base_counts`, `gc`, `sigma_eff` in the JSON                                    |
| 6 | Sketching invariance             | `make sweep` (kmer vs minimizer vs FracMinHash)                                      |
| 7 | Repeat-length resolution         | `match_length_hist` in the JSON; `make plots` -> `*_repeat_spectrum.pdf`             |

## Parameterization

```
make matchlen    GENOME=human KMAX=1024
make sweep       GENOME=dm6 KSWEEP="15 17 19 21 23" SWEEP_W=15
make blend       GENOME=sacCer3 blend_PARAMS="21 11 32 3 | 15 10 50 7"
make matchlenAll GENOMES="ecoliK12 sacCer3 human" JOBS=32
```

- **Genomes**: `<name>_FASTA` variables map a name to a path; pick one with `GENOME=<name>`, or a batch with `GENOMES="a b c"`.
- **Per-technique grids**: `<m>_PARAMS`, tuples separated by `|`; each tuple is the arguments passed to `sketch-seeds` after the FASTA.
- **Method registry**: add a technique by copying one 5-line block (`<m>_MACRO/_LIB/_SRC/_CC/_PARAMS`) and appending the name to `METHODS`.