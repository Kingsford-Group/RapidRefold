# RapidRefold

Fast incremental RNA secondary structure prediction under point mutations.

RapidRefold exploits the observation that a single-nucleotide mutation only
invalidates a subset of the dynamic programming (DP) subproblems in the Zuker
MFE algorithm. Instead of recomputing the full O(n^3) fold from scratch, it
selectively updates only the affected DP entries, achieving 2-80x speedups
depending on sequence length and optimization mode.

## Features

- **Two modes:** exact incremental (guaranteed identical to full recomputation)
  and top-K heuristic (faster, >92% structural accuracy)
- **Binary cache:** fold once, refold many mutations instantly
- **Parallel execution:** distribute mutations across threads
- **Saturation mutagenesis experiment:** systematic benchmark comparing exact
  vs heuristic across K=1,5,10 with structural accuracy metrics
- **Cross-platform:** Linux and macOS

## Quick Start

```bash
git clone --recursive <this-repo-url>
cd RapidRefold
make setup      # Build ViennaRNA dependency (one-time, ~5 min)
make            # Build RapidRefold

# Create a cache from a FASTA file
./build/rapidrefold cache -s examples/sample.fa -k 5 -o my_cache.bin

# Refold under a single mutation (position 5 → C)
./build/rapidrefold fold my_cache.bin -m 5 C

# Refold under multiple mutations from a file
./build/rapidrefold fold my_cache.bin -M mutations.txt -j 4 -o results.csv

# Run the saturation mutagenesis experiment (paper benchmark)
./build/fulltopk_experiment examples/sample.fa -o experiment_results.csv
```

## Installation

### Prerequisites

- C compiler (GCC or Clang)
- GNU Make
- Git (for submodule)
- Autotools (automake, autoconf, libtool) — needed to build ViennaRNA
- `xxd` — needed to generate hex template files for ViennaRNA

On macOS:
```bash
xcode-select --install
brew install automake autoconf libtool
```

On Ubuntu/Debian:
```bash
sudo apt-get install build-essential autoconf automake libtool xxd
```

On Fedora/RHEL:
```bash
sudo dnf install autoconf automake libtool gcc gcc-c++ vim-common
```

### Building

```bash
# Clone with submodules
git clone --recursive <repo-url>
cd RapidRefold

# Build ViennaRNA from submodule (one-time setup)
make setup

# Build RapidRefold
make
```

## Usage

### `rapidrefold cache` — Build a refolding cache

```
rapidrefold cache -s <input> [-k K | -inc] -o <output>
```

| Flag | Description |
|------|-------------|
| `-s FILE` | Input FASTA file or raw RNA sequence |
| `-k K` | Top-K heuristic mode with K candidates (default: 5) |
| `-i` / `-inc` | Exact incremental mode (no approximation) |
| `-o FILE` | Output binary cache file |

**Examples:**
```bash
# Top-K mode (fast, K=5 candidates cached per DP entry)
./build/rapidrefold cache -s sequence.fa -k 5 -o cache.bin

# Exact mode (slower but guaranteed correct)
./build/rapidrefold cache -s sequence.fa -inc -o cache_exact.bin

# Raw sequence input
./build/rapidrefold cache -s "GGGAAACCCUUUAAAGGGCCC" -k 10 -o cache.bin
```

### `rapidrefold fold` — Refold under mutations

```
rapidrefold fold <cache.bin> [-m POS BASE | -M FILE] [-j N] [-o FILE]
```

| Flag | Description |
|------|-------------|
| `<cache.bin>` | Binary cache file from the `cache` command |
| `-m POS BASE` | Single mutation (1-based position, nucleotide). Repeatable. |
| `-M FILE` | Mutation list file |
| `-j N` | Number of parallel worker threads (default: 1) |
| `-o FILE` | Output CSV file (default: stdout) |

**Examples:**
```bash
# Single mutation
./build/rapidrefold fold cache.bin -m 5 C

# Multiple mutations on the command line
./build/rapidrefold fold cache.bin -m 5 C -m 10 U -m 23 A

# Mutations from file, 4 threads, output to CSV
./build/rapidrefold fold cache.bin -M mutations.txt -j 4 -o results.csv
```

### Mutation File Format

Plain text, one mutation per line. Lines starting with `#` are comments.

```
# Position (1-based)  Base
5 C
10 U
23 A
47 G
```

### Output CSV Format

```csv
pos,orig,new,method,mfe,structure,time_ms
5,G,C,topk,-12.30,(((....)))...,0.142
10,A,U,topk,-11.80,.(((...)))...,0.098
```

| Column | Description |
|--------|-------------|
| `pos` | 1-based mutation position |
| `orig` | Original nucleotide |
| `new` | Mutated nucleotide |
| `method` | `incremental` (exact) or `topk` (heuristic) |
| `mfe` | Minimum free energy (kcal/mol) |
| `structure` | Predicted secondary structure (dot-bracket) |
| `time_ms` | Wall-clock time for this mutation (milliseconds) |

### `fulltopk_experiment` — Saturation mutagenesis benchmark

The main experiment tool for the paper. Performs saturation mutagenesis on
input sequences: for every position, tests all 3 possible single-nucleotide
substitutions. Compares exact incremental MFE (ground truth) against the
full top-K heuristic at K=1, K=5, K=10.

```
fulltopk_experiment input.fa [options]
```

| Flag | Description |
|------|-------------|
| `-o FILE` | Output CSV (default: stdout) |
| `-s PCT` | Stratified sample of positions (1-100%, default: 100) |
| `-h` | Show help |

**Output CSV columns:**
`seq_id, length, pos, orig, new, method, mfe, time_ms, mfe_diff, sensitivity, ppv, f1`

**Methods reported:**
- `initial_fold` — baseline MFE fold timing (median of 3 repeats)
- `topk_cache_KN` — time to build top-K cache for K=N
- `gt_mfe` — ground-truth exact incremental MFE
- `topk_full_KN` — full top-K heuristic result with K=N

**Examples:**
```bash
# Full saturation mutagenesis
./build/fulltopk_experiment sequences.fa -o results.csv

# 10% stratified sampling (faster for long sequences)
./build/fulltopk_experiment sequences.fa -s 10 -o results.csv
```

The output provides per-mutation structural accuracy (sensitivity, PPV, F1)
comparing each top-K prediction against the exact ground truth, enabling
systematic evaluation of the accuracy-speed tradeoff across K values.

## Performance

Typical speedups over full recomputation (`vrna_mfe`):

| Length | Exact Incremental | Top-K (K=5) |
|--------|-------------------|-------------|
| 100 nt | 2-5x | 5-10x |
| 300 nt | 5-15x | 15-40x |
| 500 nt | 10-30x | 30-60x |
| 1000 nt | 20-50x | 50-80x |

Top-K heuristic accuracy (K=5): structural F1 > 0.92 across benchmarks.

## Algorithm

### Exact Incremental

When position `p` mutates, only DP entries `c[i][j]` where `i <= p <= j` are
affected. The algorithm:

1. Updates sequence encodings and pair-type tables
2. Recomputes affected `c[i][j]` entries (pair decomposition)
3. Propagates changes through `fML[i][j]` (multibranch stems)
4. Recomputes the exterior loop `f5`
5. Backtracks to produce the new structure

Composable optimizations:
- **Lazy setup** (`VRNA_MFE_INC_LAZY`): O(n) ptype/hc update instead of O(n^2)
- **Sparse fill** (`VRNA_MFE_INC_SPARSE`): skip fML rows with no value changes

### Top-K Heuristic

Pre-caches the K best bifurcation split points and interior loop candidates
from the initial fold. During refolding, replaces O(n) scans with O(K) lookups:

- Bifurcation: K best split points for each fML entry
- MB loops: K best closing pair decompositions
- Interior loops: K best (k,l) inner pairs per (i,j)

This reduces the per-mutation complexity from O(n^2) to O(n^2 * K/n) ≈ O(n*K).

## Testing

```bash
make test
```

Runs:
- 6 basic correctness tests (exact vs full recomputation)
- Exhaustive 100bp test (300 mutations, exact match)
- Top-K accuracy evaluation (300 mutations, reports match %)
- Composable flag combination tests (7 combinations)

## ViennaRNA Submodule Modifications

During `make setup`, the build applies the following patch to the ViennaRNA
source tree (`extern/ViennaRNA/src/cephes/const.c`):

| Original | Patched |
|----------|---------|
| `double NAN = 1.0/0.0 - 1.0/0.0;` | `double NAN = __builtin_nan("");` |
| `double INFINITY = 1.0/0.0;` | `double INFINITY = __builtin_inf();` |

**Why:** The original cephes code computes NAN and INFINITY via floating-point
division-by-zero arithmetic. Modern compilers (GCC and Clang) on both Linux and
macOS treat this as undefined behavior under strict optimization, producing
warnings or incorrect constant folding. The patch replaces these with compiler
builtins that are well-defined on all platforms.

Additionally, the build only compiles the ViennaRNA library (`src/ViennaRNA/`)
rather than the full project (which includes CLI tools requiring `gengetopt`).
This reduces build dependencies and avoids failures on systems without
`gengetopt` installed.

## License

This project uses ViennaRNA (https://github.com/ViennaRNA/ViennaRNA) as a
dependency, which is distributed under its own license. See
`extern/ViennaRNA/COPYING` for details.
