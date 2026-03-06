
# Code for “Boosting Streaming Algorithms for Minimum Cost Submodular Cover with Theoretical Guarantees”

This repository contains the code used to evaluate streaming algorithms for the **Minimum Cost Submodular Cover (MCSC)** problem. It includes the implementations of the proposed methods, baseline algorithms, graph preprocessing tools, synthetic graph generators, and executables for several objective settings used in the experiments.

## Overview

The code is organized around a common experimental pipeline for graph-based MCSC instances. The repository supports multiple application settings, including coverage threshold, revenue threshold, maximum cut, and influence threshold. It also provides utilities to preprocess raw edge-list data into the binary graph format required by the executables.

## Repository Structure

.
├── data/                    # Input graphs in binary format
├── run/                     # Compiled executables
├── src/
│   ├── main.cpp             # Main entry point for experiments
│   ├── algs.cpp             # Proposed methods and baseline algorithms
│   ├── obj.cpp              # Objective functions
│   ├── mygraph.cpp          # Graph structures and I/O
│   ├── preprocess.cpp       # General graph preprocessing
│   ├── preprocess_ic.cpp    # Preprocessing for influence threshold instances
│   ├── gen_er.cpp           # Erdős–Rényi graph generator
│   ├── gen_ba.cpp           # Barabási–Albert graph generator
│   ├── logger.cpp           # Logging utilities
│   └── binheap.cpp          # Auxiliary heap structure
├── Makefile
└── LICENSE

## Implemented Methods

### Proposed methods

The three proposed algorithms in the paper are implemented in the codebase as follows:

* **SingStream** — code option: `ALG1`
* **ThreeStr** — code option: `ALG2`
* **MultiStr** — code option: `ALG3`

### Baselines

The repository also includes the following baseline methods:

* **SINGLE**
* **MULTI**
* **GREEDY**

## Supported Objective Modes

The `Makefile` builds different executables for different objective settings.

* `maxcov` — maximum coverage / coverage threshold setting
* `maxcut` — maximum cut setting
* `revmax` — revenue threshold, monotone mode
* `revmax_nm` — revenue threshold, non-monotone mode
* `ic` — influence threshold setting with OpenMP support

## Build Instructions

The project uses `g++` with C++14. The default optimization flags are `-std=c++14 -Wall -O2`. The `ic` target additionally uses OpenMP.

To build all main executables and tools:

```bash
make
```

This command creates the `run/` directory and builds:

* `run/maxcov`
* `run/maxcut`
* `run/revmax`
* `run/revmax_nm`
* `run/ic`
* `run/preproc`
* `run/preproc_ic`

To build individual targets:

```bash
make maxcov
make maxcut
make revmax
make revmax_nm
make ic
make preproc
make preproc_ic
make er
make ba
```

To remove compiled binaries:

```bash
make clean
```

## Data Preparation

### General preprocessing

The executable `preproc` converts an unweighted edge list into the binary graph format used by the main programs.

Build:

```bash
make preproc
```

Usage:

```bash
./run/preproc <input_edge_list> <output_binary_graph>
```

Input assumptions:

* the input file is an unweighted edge list,
* the graph is simplified,
* isolates are removed,
* vertices are renumbered before writing the binary file.

### Preprocessing for influence threshold

The executable `preproc_ic` prepares graph data for the influence threshold setting.

Build:

```bash
make preproc_ic
```

Usage:

```bash
./run/preproc_ic edges.txt output.bin [input_undirected(0/1)=1] [seed=42]
```

Expected input format:

* each line contains `u v [w]`,
* `u` and `v` are nonnegative integer node IDs,
* `w` is optional and defaults to `1.0`.

Processing rules:

* self-loops are removed,
* duplicate edges keep the first occurrence,
* incoming weights are normalized,
* node weights are generated randomly in `(0,1)`,
* if `input_undirected=1`, both directions are added.

## Datasets

All datasets used in the experiments are taken from **SNAP**:

[https://snap.stanford.edu/data/](https://snap.stanford.edu/data/)

The experimental datasets are listed below.

| Application         | Dataset  | Vertices |   Edges | Type       |
| ------------------- | -------- | -------: | ------: | ---------- |
| Revenue Threshold   | Facebook |     4039 |   88234 | Undirected |
| Revenue Threshold   | GrQc     |     5242 |   14496 | Undirected |
| Revenue Threshold   | Enron    |    36692 |  183831 | Undirected |
| Coverage Threshold  | AstroPh  |    18772 |  198110 | Undirected |
| Coverage Threshold  | Hept     |    34546 |  421578 | Undirected |
| Coverage Threshold  | Stanford |   281903 | 2312497 | Directed   |
| Influence Threshold | Email    |     1005 |   25571 | Undirected |

## Running Experiments

The main executables are built under `run/`. The examples below assume that the current working directory is `run/`. If you run from the project root, replace `./maxcov` by `./run/maxcov`, and similarly for the other binaries.

### Main arguments

The main program supports the following arguments.

* `-g <file>` : input graph in binary format
* `-a <algorithm>` : algorithm name
* `-t <T>` : threshold value
* `-e <epsilon>` : accuracy parameter, default `0.1`
* `-N <repetitions>` : number of repetitions
* `-o <file>` : output file name
* `-k <value>` : cardinality constraint
* `-b <value>` : budget
* `-q` : quiet mode
* `-l` : lazy mode
* `--nthreads <value>` : number of threads

### Algorithm names accepted by the code

For the proposed methods, use the following options:

* `-a ALG1` for **SingStream**
* `-a ALG2` for **ThreeStr**
* `-a ALG3` for **MultiStr**

For baselines, use:

* `-a SINGLE`
* `-a MULTI`
* `-a GREEDY`

### Example commands

Revenue threshold:

```bash
./revmax -g ../data/fb.bin -t 4272 -a SINGLE -e 0.2 >> ../output/revmax02.csv
```

Coverage threshold:

```bash
./maxcov -g ../data/google.bin -t 17514 -a SINGLE >> ../output/maxcov.csv
```

Influence threshold:

```bash
./ic -g ../data/email.bin -t 101 -a SINGLE >> ../output/email.csv
```

Examples using the proposed methods:

```bash
./revmax -g ../data/fb.bin -t 4272 -a ALG1 -e 0.2 >> ../output/singstream_rev.csv
./revmax -g ../data/fb.bin -t 4272 -a ALG2 -e 0.2 >> ../output/threestr_rev.csv
./revmax -g ../data/fb.bin -t 4272 -a ALG3 -e 0.2 >> ../output/multistr_rev.csv
```

```bash
./maxcov -g ../data/google.bin -t 17514 -a ALG1 >> ../output/singstream_cov.csv
./ic -g ../data/email.bin -t 101 -a ALG2 >> ../output/threestr_ic.csv
```

## Synthetic Graph Generation

The repository also includes two graph generators.

### Erdős–Rényi graphs

Build:

```bash
make er
```

This produces:

```bash
./run/er
```

### Barabási–Albert graphs

Build:

```bash
make ba
```

This produces:

```bash
./run/ba
```

These generators are useful for creating synthetic instances for controlled experiments.

## Batch Experiments

The `Makefile` defines a batch target:

```bash
make single-pass
```

This target enters `exp/single-pass` and runs a sequence of experiments for selected datasets and binaries. According to the current build script, it uses the dataset identifiers:

* `er`
* `ba`
* `fb`
* `slashdot`
* `pokec`

and the binaries:

* `maxcut`
* `revmax`

## Output Format

The main algorithms print experimental results in a comma-separated format. The standard header used by most runs is:

```text
alg,T,f,q,cost,pass,memory,memoryUnused,time
```

The fields have the following meaning.

* `alg` : algorithm name
* `T` : threshold value
* `f` : achieved objective value
* `q` : number of objective evaluations / queries
* `cost` : cost of the returned solution
* `pass` : number of passes
* `memory` : memory used by maintained candidates or structures
* `memoryUnused` : additional reserved but unused memory
* `time` : running time in seconds

When output is redirected using `>>`, results from multiple runs can be appended to the same CSV file.

## Reproducibility Notes

For consistent experiments, it is recommended to:

* preprocess all datasets using the provided tools,
* keep the same threshold and epsilon settings across methods,
* run all methods on the same binary graph files,
* use fixed seeds when generating synthetic data,
* record outputs in separate files for each objective setting.

The influence-threshold preprocessing and the revenue-related settings include randomized components, so fixing seeds is useful when exact replication is needed.

## License

See the `LICENSE` file for licensing information.

```
```
