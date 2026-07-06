# FAVeS

**FAVeS** - a lightweight SNP caller that anchors reads to a reference using **BLEND fuzzy seeds**, then confirms candidate variants with gapped alignment (WFA2).

FAVeS leans on two properties of BLEND sketches:

- **Uniqueness** — seeds that occur too often in the reference are dropped, so reads are anchored only by rare, position-specific seeds. This keeps mapping specific and avoids ambiguous placements.
- **Error tolerance (fuzziness)** — BLEND builds each seed by SimHash-combining several neighboring k-mers, so a seed keeps the *same* value even when one of those k-mers carries a mismatch. A single SNP or sequencing error therefore does not break the anchor, and the read still maps across the variant.

---

## Dependencies

- `gcc` with AVX2 / SSE4.1 support
- `zlib` (`-lz`), `pthread`, `libm`
- **WFA2-lib** — bundled in `WFA2-lib/` (built by the Makefile)
- **BLEND** sketch source in `sketch/` (built into `lib/libblend.a`)

---

## Build

You can run following commands for installation and build the project:

```bash

git clone https://github.com/akmami/FAVeS.git
make install
```

---

## Usage

```bash
./faves -f reference.fasta -q reads.fastq -o variants.bed [options]
```

### Required

| Flag | Long | Description |
|------|------|-------------|
| `-f` | `--fasta` | Reference FASTA file |
| `-q` | `--fastq` | Reads FASTQ file |
| `-o` | `--output` | Output BED file |

### Seeding (BLEND) options

| Flag | Long | Default | Description |
|------|------|---------|-------------|
| `-k` | `--kmer` | 21 | K-mer size |
| `-w` | `--window` | 11 | Window size |
| `-b` | `--blend-bits` | 50 | Number of hash bits |
| `-n` | `--n-neighbors` | 5 | Neighbors combined per fuzzy seed |
| `-r` | `--radius` | 4 | Seed span radius for alignment |

### Variant calling options

| Flag | Long | Default | Description |
|------|------|---------|-------------|
| `-c` | `--consensus` | 10 | Minimum supporting reads for a SNP |
| `-d` | `--consensus-frac` | 0.5 | Minimum support as a fraction of local depth |

### Runtime options

| Flag | Long | Default | Description |
|------|------|---------|-------------|
| `-t` | `--threads` | 4 | Worker threads (1–1024) |
| `-p` | `--progress` | off | Show progress |
| `-v` | `--verbose` | off | Verbose messages |
| `-h` | `--help` | — | Show help |

---

### Output (BED)

Each line reports one SNP:

```
chrom    start(0-based)    end(1-based)    ref    alt    support
```

`support` is the number of reads that agreed on the variant (the consensus count).

---

## Tuning notes

- **Specificity vs. sensitivity:** larger `-k` and stricter uniqueness filtering give more specific anchors; smaller values map more reads at the cost of ambiguity.
- **Fuzziness:** `-n` controls how many neighboring k-mers are blended per seed — more neighbors means greater error tolerance (seeds survive more mismatches) but coarser localization.
- **Confidence:** raise `-c` / `-d` for high-confidence calls on deep data; lower them for shallow coverage.

## License

`FAVeS` is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/akmami/FAVeS/blob/main/LICENSE).