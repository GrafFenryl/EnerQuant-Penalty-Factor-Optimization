# Unit Commitment Benchmark Instances and Results

This repository contains the final unit commitment problem (UCP) instances,
solution records, reference MILP solutions, summary statistics, and standalone
Python utilities used for the accompanying publication.

The repository is self-contained: the instance, QUBO, and MILP utilities run
with the files and dependencies described below. All example paths are relative
to this repository's root directory, which contains this README.

## Repository contents

```text
.
|-- README.md
|-- final_statistics.xlsx
|-- instances/
|   `-- p_ucp_u{units}_t{timesteps}_v{version}
|-- results/
|   |-- fix_multiobj/
|   |-- it_biobj_repair/
|   |-- it_multiobj_repair/
|   |-- nelder_mead/
|   |-- nelder_mead_repair/
|   |-- nelder_mead_random_repair/
|   |-- nelder_mead_random_repair_old/  # Historical archive
|   `-- gurobi/
`-- code/
    |-- qubo_tools.py
    |-- qubo_usage.ipynb
    |-- ucp_mip.py
    `-- ucp_mip_usage.ipynb
```

## Instances

The `instances` directory contains 185 curated deterministic UCP instances.
Every retained combination of units and timesteps has five versions, numbered
`v1` through `v5`.

| Units | Retained timesteps | Instances |
|---|---|---:|
| 2 | 3, 6, 9, 12, 15, 18, 21, 24 | 40 |
| 4 | 3, 6, 9, 12, 15, 18, 21, 24 | 40 |
| 6 | 3, 6, 9, 12, 15, 18, 21, 24 | 40 |
| 8 | 3, 6, 9, 12, 15, 18 | 30 |
| 10 | 3, 6, 9, 12 | 20 |
| 12 | 3, 6, 9 | 15 |
| **Total** | | **185** |

Instance names follow this convention:

```text
p_ucp_u{units}_t{timesteps}_v{version}
```

For example, `p_ucp_u2_t3_v1` contains two thermal units, three timesteps, and
is version 1 of that size combination. The files intentionally have no file
extension.

Each instance is a JSON object with three fields:

- `units`: unit names, minimum and maximum generation, minimum up/down times,
  start costs, and variable costs;
- `supply`: exogenous supply for each timestep;
- `demand`: demand for each timestep.

The three tables use the JSON representation produced by the original data
pipeline. The standalone readers in `code/qubo_tools.py` and `code/ucp_mip.py`
decode this representation directly.

The publication filenames are definitive; no historical name mapping is needed
to use the repository.

## Result files

The six methods included in the workbook are `fix_multiobj`, `it_biobj_repair`,
`it_multiobj_repair`, `nelder_mead`, `nelder_mead_repair`, and
`nelder_mead_random_repair`. Each method directory contains one result for each
of the 185 instances. The directory `results/nelder_mead_random_repair_old/`
contains historical results and is excluded from the workbook.

Heuristic result files follow this convention:

```text
results/{method}/p_ucp_{method}_u{units}_t{timesteps}_v{version}.json
```

Each result contains metadata and six selected solution states:

- `best_original_solution_overall`
- `best_repaired_solution_overall`
- `best_original_solution_best_iteration`
- `best_repaired_solution_best_iteration`
- `best_original_solution_last_iteration`
- `best_repaired_solution_last_iteration`

The states were selected with the TDO criteria: technical feasibility, demand
feasibility, and objective quality. If the same state first appears more than
once, its earliest occurrence is used.

Every selected solution contains:

- `iteration_number`;
- `QUBO_value`;
- `split_energies`;
- `penalty_factors`;
- `decision_variables`.

The QUBO components use the following labels:

| Label | Meaning |
|---|---|
| `cost matrix` | Generation and start-up objective |
| `a` | Demand-balance penalty |
| `b` | Unit-status/start coupling penalty |
| `c` | Minimum up-time penalty |
| `d` | Minimum down-time penalty |
| `e` | Unit-status/generation-level coupling penalty |

Some consolidated outputs report `c+d`. This is the sum of the minimum up-time
and minimum down-time energies; both components use the same penalty factor.

### Reference MILP solutions

Reference results are stored under `results/gurobi` and follow this convention:

```text
results/gurobi/p_ucp_gurobi_u{units}_t{timesteps}_v{version}.json
```

Each file contains an objective value and the continuous generation (`gen`) and
binary commitment (`on`, `start`, and `stop`) variables. These variables are
from the MILP formulation and therefore differ from the logarithmically
discretized `level` variables in the QUBO results.

The reference objectives are instance-level quantities: a given
`u{units}_t{timesteps}_v{version}` instance has the same reference objective for
every heuristic method.

The workbook calls this objective `Gurobi_optimal_value`, but the reference is
not a proven optimum for every instance. Thirteen instances retain a nonzero
reported solver gap, with a maximum of approximately 0.00986%. The workbook's
`Gurobi_time` and `Gurobi_mip_gap` describe the historical reference solves;
they do not establish the runtime or optimality certificate of later objective
updates.

## Final statistics workbook

`final_statistics.xlsx` contains one sheet, `Statistics`, with 1,110 data
rows: one for each of the 185 retained instances and six included heuristic methods.
The workbook combines:

- instance dimensions and constraint counts;
- logical and physical qubit counts;
- iteration, evaluation, and timing statistics;
- QUBO values, split energies, penalty factors, and objective gaps for four
  selected states: original/repaired overall best and original/repaired last
  iteration;
- reference MILP objective, runtime, and solver gap.

The best-iteration states remain available in the result JSON files, giving six
selected states per file, but are omitted from the workbook.

The Gurobi reference columns are copied from one canonical instance-level record and are
identical across methods for the same instance.

### Timing statistics

All timing values are in seconds:

- `time_embedding`: historical embedding time for the instance, shared across
  the five hardware methods and zero for `nelder_mead_random_repair`;
- `time_CPU`: median classical-host processing time from the optimized timing
  replay, excluding external sampler waiting and campaign initialization;
- `time_QPU`: historical provider-reported preparation and QPU execution time,
  zero for the local random method;
- `time_compute_total`: `time_CPU + time_QPU`;
- `time_postprocessing`: the median postprocessing time, already included in
  `time_CPU`;
- `time_communication_overhead`: estimated communication overhead **per
  evaluation**, reported separately from `time_compute_total`.

For the hardware methods, the communication estimate averages positive
observations of `(historical sampler wall time - provider-reported QPU time) /
number of evaluations`, pooled across the five hardware methods for each
instance. The resulting estimate is shared by those methods. For instances
without a positive observation, the estimate is the same-size peers' mean plus
a seeded normal perturbation using their sample standard deviation, clipped
at zero. The estimate is zero for `nelder_mead_random_repair`. It is not a
directly measured total runtime for the row.

### Solution metrics

`sum technical matrix` is the unweighted sum of `start matrix`,
`min-up/down matrix`, and `discretization matrix`. These split-energy components
are signed, so cancellation can produce a zero sum even when individual
components are nonzero. A zero sum alone does not certify technical feasibility.

In each selected-state block, `MIP gap (%)` is the absolute objective difference
relative to the reference, calculated as
`100 * abs(objective matrix / Gurobi_optimal_value - 1)`. Here `objective matrix`
is the unweighted cost component; the complete weighted `QUBO value` also
includes penalties and must not be substituted for it. A zero objective gap
alone does not certify feasibility.

This objective comparison is distinct from `Gurobi_mip_gap`, the solver's
relative incumbent-to-bound gap stored as a fraction (multiply by 100 to obtain
a percentage).

## Standalone QUBO utilities

`code/qubo_tools.py` requires only the Python standard library and exposes three
main functions:

- `build_sub_QUBO`: construct the unweighted objective and constraint
  sub-QUBOs;
- `build_QUBO`: apply the penalty factors from a selected result state and
  construct the full QUBO;
- `evaluate_solution`: evaluate a selected result state on a constructed QUBO.

Run the complete example in `code/qubo_usage.ipynb`, or use the functions
directly from the repository root. Start Jupyter with `jupyter notebook` from
the repository root. The notebooks automatically use the repository root as
their working directory, including when the kernel starts in `code/`:

```python
from pathlib import Path
import sys

sys.path.insert(0, str(Path("code").resolve()))

from qubo_tools import build_sub_QUBO, build_QUBO, evaluate_solution

instance = Path("instances/p_ucp_u2_t3_v1")
result = Path(
    "results/fix_multiobj/"
    "p_ucp_fix_multiobj_u2_t3_v1.json"
)
solution = "best_repaired_solution_overall"

parts = build_sub_QUBO(instance)
qubo = build_QUBO(instance, result, solution)
value = evaluate_solution(qubo, result, solution)

print(qubo.summary())
print(value)
```

`build_sub_QUBO` reproduces the publication formulation, including generation
step size 1, initially-off unit status, logarithmic generation discretization,
and minimum up/down-time lookahead windows that omit the current timestep
(offset zero). This window convention applies throughout the horizon, including
at its first timestep.

## Standalone MILP utilities

`code/ucp_mip.py` builds a solver-neutral mixed-integer linear UCP. The model
uses continuous generation variables and binary `on`, `start`, and `stop`
variables. It includes demand balance, generation bounds, commitment
transitions, and minimum up/down-time constraints.

The module supports:

- solving with `scipy.optimize.milp`, which uses HiGHS;
- constructing an equivalent `gurobipy.Model` for Gurobi.

Install the dependencies needed for the SciPy workflow:

```text
python -m pip install numpy scipy jupyter
```

Then run `code/ucp_mip_usage.ipynb`, or use the module directly from the
repository root:

```python
from pathlib import Path
import sys

sys.path.insert(0, str(Path("code").resolve()))

from ucp_mip import build_mip, solve_with_scipy

problem = build_mip("instances/p_ucp_u2_t3_v1")
solution = solve_with_scipy(
    problem,
    time_limit=300.0,
    mip_rel_gap=0.0,
)

print(solution.objective_function_value)
print(solution.mip_gap)
print(solution.solver_status, solution.solver_message)
```

To construct the same model for Gurobi, install `gurobipy` and use a valid
Gurobi license:

```python
from ucp_mip import build_gurobi_model, solution_from_gurobi

model, variables, problem = build_gurobi_model(
    "instances/p_ucp_u2_t3_v1"
)
model.Params.MIPGap = 0.0
model.optimize()
if model.SolCount == 0:
    raise RuntimeError(f"Gurobi returned no feasible solution (status {model.Status}).")

solution = solution_from_gurobi(problem, variables, model)
print(solution.objective_function_value)
print(solution.solver_status, solution.mip_gap)
```

Setting `MIPGap` to zero requests zero relative optimality gap. Check the
solver's termination status and achieved gap before treating a result as an
exact reference; setting the requested tolerance alone does not certify it.

## Software requirements

- Python 3.10 or newer;
- no third-party dependency for QUBO construction and evaluation;
- NumPy and SciPy for solving the MILP with HiGHS;
- Jupyter for running the example notebooks;
- optional `gurobipy` and a Gurobi license for constructing and solving the
  Gurobi model;
- optional `openpyxl` for programmatic access to the Excel workbook.

## Reproducibility notes

- Run the inline examples from the repository root. The notebooks set their
  working directory to the same root. All example file paths refer to `code/`,
  `instances/`, or `results/` below that directory.
- The standalone utilities reconstruct the QUBO and MILP formulations and
  evaluate exported solutions. The full penalty-optimization, sampling, timing
  replay, and workbook-export pipelines are not included.
- QUBO energies include the constant offset.
- The decision-variable names are stored as tuple strings to reproduce the
  original result format exactly.
- The QUBO and MILP use different generation representations: discrete `level`
  variables in the QUBO and continuous `gen` variables in the MILP.
- Use the complete instance key (units, timesteps, and version) when comparing
  results across methods.

## Citation

If you use these instances, results, or utilities, please cite the accompanying
publication. Its bibliographic details will accompany the published release.
