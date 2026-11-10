# The angular-metric close enough traveling salesman problem

Code, instances and results of the paper

> E. G. Birgin, L. C. Coelho and P. C. R. Ertel. *The angular-metric close enough traveling salesman problem: an adaptive large neighborhood search enhanced with block coordinate descent.* Submitted, 2026.

In the angular-metric close enough traveling salesman problem (Angle-CETSP), a closed route must pass through every disk of a given set in the plane, and its cost is the sum of its turning angles. The paper proposes the ALNS-BCD, an adaptive large neighborhood search in which the visit points of each visiting sequence are computed by a block coordinate descent (BCD) method whose single-disk subproblem is solved in closed form.

## Contents

```
code/alns/              the ALNS-BCD (Fortran 90 and C++17)
code/block-rules/       comparison of the block selection rules of the BCD (Section 5.1)
instances/              the 400 instances derived from the Angle-TSP instances (Section 5.2)
results/                results of every run reported in the paper
```

## Building and running the ALNS-BCD

Requirements: `g++` with C++17, `gfortran` and Python 3 (tested with GCC 12.3 and 13.3 on Ubuntu).

```bash
cd code/alns
bash build.sh
python3 run_alns.py ../../instances/angle-tsp-radii/angletsp-n050-01.txt 123456
```

`run_alns.py` runs the configuration of the paper (Section 5.3): N = 200,000 iterations, a time limit of one hour, and the parameters of the paper. Its docstring describes the output: a log with the cost of each new best solution and the statistics of the operators, and a file with the best solution (route and visit points). The ALNS-BCD is deterministic for a given seed. The five seeds used for all instances are 123456, 42, 34567, 99 and 2026, and the 400 Angle-TSP instances were also solved with the seeds 7, 2718, 31415, 1618 and 8128.

`bash build.sh --sleef` builds a faster version, in which the `atan2` of the closed-form block update is computed by the SLEEF library (it requires `libsleef`). Its results agree with those of the default build up to rounding, but they are not bit-identical, so a run may follow a different trajectory. The runs of the paper used the default build.

Instance file format: the number m of disks in the first line, followed by one line `x y r` per disk. A depot is a disk of radius zero. The program translates and scales each instance so that its centers fit in the unit square, which does not change the cost of any route, and writes the solutions in the original coordinates.

## Instances

`instances/angle-tsp-radii/` contains the 400 instances built in Section 5.2 from the Angle-TSP instances of Staněk et al. (2019): each point c_i becomes the center of a disk of radius r_i = 0.45 min_{j≠i} ||c_i − c_j||. The instance `angletsp-nNNN-KK` is the point set KK (01 to 10) with NNN points of Staněk et al. (2019).

The other instances of the paper are available in the literature:

* the 400 Angle-TSP instances (disks of radius zero), their optimal values for up to 75 points and the results of the heuristics of Staněk et al. (2019) are in the ancillary files of arXiv:1803.03681;
* the 62 two-dimensional CETSP instances of Mennell (2009) and the 720 instances of Behdani and Smith (2014) are those used by Coutinho et al. (2016).

## Results

All costs are total turning angles in radians, and all times are CPU times in seconds.

| file | content |
|---|---|
| `results/alns-bcd-angle-tsp.csv` | one row per run of the ALNS-BCD on the 400 Angle-TSP instances, 10 seeds (4,000 rows) |
| `results/alns-bcd-disks.csv` | one row per run of the ALNS-BCD on the 1,182 instances with disks, 5 seeds (5,910 rows) |
| `results/angle-tsp-bks.csv` | best known solution (BKS) of each Angle-TSP instance |
| `results/gts-comparison-by-size.csv` | the comparison with the granular tabu search of Cavagnini et al. (2024), for each of the 40 sizes |
| `results/hexaly-disks.csv` | one row per Hexaly run on the 1,182 instances with disks |
| `results/block-rules/radii-x1.csv`, `-x4.csv`, `-x10.csv` | the comparison of the block selection rules of Section 5.1 |

### `alns-bcd-angle-tsp.csv` and `alns-bcd-disks.csv`

| column | meaning |
|---|---|
| `set` | `angle-tsp`, `angle-tsp-radii`, `behdani-r0.25`, `behdani-r0.50`, `behdani-r1.00` or `mennell` |
| `instance`, `m`, `seed` | instance, number of disks (depot included) and seed |
| `iterations` | N = 200,000 |
| `initial_cost`, `final_cost` | cost of the initial solution and of the best solution of the run |
| `time_limit_reached` | 1 if the run was stopped by the time limit of one hour |
| `construction_time_s` | time of reading, preprocessing and building the initial solution |
| `search_time_s` | time of the search, including the final BCD |
| `total_time_s` | `construction_time_s + search_time_s`; for runs stopped by the time limit, the wall time |
| `last_improvement_iteration`, `last_improvement_time_s` | iteration and search time of the last new best solution |
| `improvements` | number of new best solutions |
| `max_iterations_without_improvement` | longest stretch of iterations without a new best solution |
| `<operator>_calls`, `<operator>_new_best` | number of calls of each operator and number of new best solutions it produced |

The operators are the removal operators `worst_angles`, `random`, `cluster`, `random_string`, `bad_string`, `good_string`, `worst_sum_string`, `savings`, `sisr_random_split`, `sisr_spatial` and `lens`, and the insertion operators `greedy`, `regret2`, `regret3` and `sisr_greedy_blinks` (Section 4 of the paper). The 33 runs stopped by the time limit, all on instances of Mennell (2009) with 400 to 1,001 disks, have no final statistics: `final_cost` is the last new best solution and the operator columns are empty.

### `angle-tsp-bks.csv`

`bks` is the BKS in radians, `bks_type` is `optimal` (instances with up to 75 points) or `heuristic`, `bks_source` gives the method of Staněk et al. (2019) that found it, `lower_bound` is the lower bound reported there, and `published_bks` is the value as published (turning angles multiplied by 1,000). The gap of a cost f is 100 (f − bks) / bks.

### `gts-comparison-by-size.csv`

For each size m, the mean gaps to the BKS (in percent) of the best of ten runs and of the mean value of ten runs, and the mean time of a run, for the granular tabu search GTS-angular r-LENS, as reported by Cavagnini et al. (2024), and for the ALNS-BCD.

### `hexaly-disks.csv`

`total_turning_angle` is the cost of the route returned by Hexaly 14.5, recomputed with the cost of the paper, and is the value used in the comparisons. `reported_objective` is the value reported by Hexaly, `hexaly_status` its final status, `stopped_by_stagnation` is 1 when the run stopped after ten minutes without improvement, and `max_violation` is the largest distance of a visit point to its disk, in the coordinates of the instance file.

### `block-rules/radii-x*.csv`

For each random instance of Section 5.1 (`code/block-rules/instances/`), the cost of the solution (`psi_*`), the number of block updates (`updates_*`), the number of cycles (`cycles_*`) and the running time (`time_*_s`) of the cyclic rule and of the queue. To reproduce them:

```bash
cd code/block-rules
bash run.sh
```

## References

* E. Behdani and J. C. Smith. An integer-programming-based approach to the close-enough traveling salesman problem. *INFORMS Journal on Computing* 26(3):415–432, 2014.
* R. Cavagnini, M. Schneider and A. Theiß. A tabu search with geometry-based sparsification methods for angular traveling salesman problems. *Networks* 83(1):30–52, 2024.
* W. P. Coutinho, R. Q. do Nascimento, A. A. Pessoa and A. Subramanian. A branch-and-bound algorithm for the close-enough traveling salesman problem. *INFORMS Journal on Computing* 28(4):752–765, 2016.
* W. K. Mennell. *Heuristics for Solving Three Routing Problems: Close-Enough Traveling Salesman Problem, Close-Enough Vehicle Routing Problem, Sequence-Dependent Team Orienteering Problem.* PhD thesis, University of Maryland, 2009.
* R. Staněk, P. Greistorfer, K. Ladner and U. Pferschy. Geometric and LP-based heuristics for angular travelling salesman problems in the plane. *Computers & Operations Research* 108:97–111, 2019. Instances and results: arXiv:1803.03681.
