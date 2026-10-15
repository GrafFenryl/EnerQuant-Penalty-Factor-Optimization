"""Standalone MILP construction and solving utilities for publication UCPs.

The module has no dependency on the EnerQuant repository.  SciPy is imported
only by :func:`solve_with_scipy`; gurobipy is imported only by
:func:`build_gurobi_model`.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence


INSTANCE_PATTERN = re.compile(
    r"^p_ucp_u(?P<units>\d+)_t(?P<timesteps>\d+)_v(?P<version>\d+)$"
)


@dataclass(frozen=True)
class Unit:
    """Thermal-unit parameters used by the MILP."""

    name: str
    min_power: float
    max_power: float
    min_uptime: int
    min_downtime: int
    start_cost: float
    var_cost: float


@dataclass(frozen=True)
class UCPInstance:
    """Parsed publication instance."""

    name: str
    units: tuple[Unit, ...]
    supply: tuple[float, ...]
    demand: tuple[float, ...]

    @property
    def num_units(self) -> int:
        return len(self.units)

    @property
    def num_timesteps(self) -> int:
        return len(self.demand)


@dataclass(frozen=True)
class ConstraintRow:
    """One sparse linear row with lower and upper bounds."""

    coefficients: Mapping[int, float]
    lower_bound: float
    upper_bound: float
    name: str


@dataclass(frozen=True)
class MIPProblem:
    """Solver-neutral mixed-integer linear problem."""

    instance: UCPInstance
    variable_keys: tuple[tuple[str, int, int], ...]
    objective: tuple[float, ...]
    integrality: tuple[int, ...]
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]
    constraints: tuple[ConstraintRow, ...]

    @property
    def num_variables(self) -> int:
        return len(self.variable_keys)

    @property
    def num_constraints(self) -> int:
        return len(self.constraints)


@dataclass(frozen=True)
class MIPSolution:
    """Normalized solution returned by a MILP solver."""

    objective_function_value: float
    decision_variables: Mapping[str, float | int]
    raw_values: tuple[float, ...]
    solver_status: int
    solver_message: str
    mip_gap: float | None


def _decode_table(value: Any, field: str) -> Mapping[str, Any]:
    decoded = json.loads(value) if isinstance(value, str) else value
    if not isinstance(decoded, Mapping):
        raise ValueError(f"Instance field {field!r} is not a JSON object.")
    return decoded


def _ordered_column(table: Mapping[str, Any], column: str) -> list[Any]:
    values = table.get(column)
    if not isinstance(values, Mapping):
        raise ValueError(f"Missing tabular column {column!r}.")
    numeric = {int(index): value for index, value in values.items()}
    expected = list(range(len(numeric)))
    if sorted(numeric) != expected:
        raise ValueError(
            f"Column {column!r} has non-contiguous indices: {sorted(numeric)}."
        )
    return [numeric[index] for index in expected]


def read_instance(path: str | Path) -> UCPInstance:
    """Read one publication instance without repository-specific imports."""

    instance_path = Path(path)
    match = INSTANCE_PATTERN.fullmatch(instance_path.name)
    if match is None:
        raise ValueError(f"Invalid publication instance name: {instance_path.name!r}.")
    raw = json.loads(instance_path.read_text(encoding="utf-8-sig"))
    if set(raw) != {"units", "supply", "demand"}:
        raise ValueError(
            f"Unexpected fields in {instance_path}: {sorted(raw)}."
        )

    unit_table = _decode_table(raw["units"], "units")
    columns = {
        name: _ordered_column(unit_table, name)
        for name in (
            "name",
            "min_power",
            "max_power",
            "min_uptime",
            "min_downtime",
            "start_cost",
            "var_cost",
        )
    }
    sizes = {len(values) for values in columns.values()}
    if len(sizes) != 1:
        raise ValueError(f"Unit columns have inconsistent lengths in {instance_path}.")
    units = tuple(
        Unit(
            name=str(columns["name"][index]),
            min_power=float(columns["min_power"][index]),
            max_power=float(columns["max_power"][index]),
            min_uptime=int(columns["min_uptime"][index]),
            min_downtime=int(columns["min_downtime"][index]),
            start_cost=float(columns["start_cost"][index]),
            var_cost=float(columns["var_cost"][index]),
        )
        for index in range(next(iter(sizes)))
    )
    for index, unit in enumerate(units):
        if not 0 <= unit.min_power <= unit.max_power:
            raise ValueError(f"Invalid power bounds for unit {index}: {unit}.")
        if unit.min_uptime < 1 or unit.min_downtime < 1:
            raise ValueError(f"Invalid minimum up/down time for unit {index}: {unit}.")

    supply_table = _decode_table(raw["supply"], "supply")
    demand_table = _decode_table(raw["demand"], "demand")
    supply = tuple(float(value) for value in _ordered_column(supply_table, "supply"))
    demand = tuple(float(value) for value in _ordered_column(demand_table, "demand"))
    if not demand or len(supply) != len(demand):
        raise ValueError(
            f"Supply and demand lengths differ in {instance_path}: "
            f"{len(supply)} != {len(demand)}."
        )
    if (len(units), len(demand)) != (
        int(match.group("units")),
        int(match.group("timesteps")),
    ):
        raise ValueError(
            f"Filename dimensions disagree with data in {instance_path}."
        )
    return UCPInstance(instance_path.name, units, supply, demand)


def build_mip(instance: UCPInstance | str | Path) -> MIPProblem:
    """Build the exact linear UCP formulation used for the reference MILP."""

    parsed = read_instance(instance) if isinstance(instance, (str, Path)) else instance
    num_units = parsed.num_units
    num_timesteps = parsed.num_timesteps

    variable_keys = tuple(
        (kind, unit, timestep)
        for kind in ("gen", "on", "start", "stop")
        for unit in range(num_units)
        for timestep in range(num_timesteps)
    )
    position = {key: index for index, key in enumerate(variable_keys)}
    objective: list[float] = []
    integrality: list[int] = []
    lower_bounds: list[float] = []
    upper_bounds: list[float] = []
    for kind, unit_index, _ in variable_keys:
        unit = parsed.units[unit_index]
        objective.append(
            unit.var_cost
            if kind == "gen"
            else unit.start_cost
            if kind == "start"
            else 0.0
        )
        integrality.append(0 if kind == "gen" else 1)
        lower_bounds.append(0.0)
        upper_bounds.append(unit.max_power if kind == "gen" else 1.0)

    rows: list[ConstraintRow] = []

    def add(
        coefficients: Mapping[tuple[str, int, int], float],
        lower: float,
        upper: float,
        name: str,
    ) -> None:
        row = {
            position[key]: float(value)
            for key, value in coefficients.items()
            if value != 0.0
        }
        rows.append(ConstraintRow(row, float(lower), float(upper), name))

    for timestep in range(num_timesteps):
        required_generation = parsed.demand[timestep] - parsed.supply[timestep]
        add(
            {("gen", unit, timestep): 1.0 for unit in range(num_units)},
            required_generation,
            required_generation,
            f"power_equilibrium_t{timestep}",
        )

    for unit_index, unit in enumerate(parsed.units):
        for timestep in range(num_timesteps):
            add(
                {
                    ("gen", unit_index, timestep): 1.0,
                    ("on", unit_index, timestep): -unit.min_power,
                },
                0.0,
                math.inf,
                f"min_power_i{unit_index}_t{timestep}",
            )
            add(
                {
                    ("gen", unit_index, timestep): 1.0,
                    ("on", unit_index, timestep): -unit.max_power,
                },
                -math.inf,
                0.0,
                f"max_power_i{unit_index}_t{timestep}",
            )
            add(
                {
                    ("start", unit_index, timestep): 1.0,
                    ("stop", unit_index, timestep): 1.0,
                },
                -math.inf,
                1.0,
                f"start_stop_uniqueness_i{unit_index}_t{timestep}",
            )

        for timestep in range(1, num_timesteps):
            add(
                {
                    ("on", unit_index, timestep): 1.0,
                    ("on", unit_index, timestep - 1): -1.0,
                    ("start", unit_index, timestep): -1.0,
                    ("stop", unit_index, timestep): 1.0,
                },
                0.0,
                0.0,
                f"on_start_stop_consistency_i{unit_index}_t{timestep}",
            )
        add(
            {
                ("on", unit_index, 0): 1.0,
                ("start", unit_index, 0): -1.0,
            },
            0.0,
            0.0,
            f"initial_on_start_consistency_i{unit_index}",
        )
        add(
            {("stop", unit_index, 0): 1.0},
            0.0,
            0.0,
            f"initial_stop_status_i{unit_index}",
        )

        for timestep in range(num_timesteps):
            uptime_start = max(timestep - unit.min_uptime + 1, 0)
            uptime = {
                ("start", unit_index, step): 1.0
                for step in range(uptime_start, timestep + 1)
            }
            uptime[("on", unit_index, timestep)] = -1.0
            add(
                uptime,
                -math.inf,
                0.0,
                f"min_uptime_i{unit_index}_t{timestep}",
            )

            downtime_start = max(timestep - unit.min_downtime + 1, 0)
            downtime = {
                ("stop", unit_index, step): 1.0
                for step in range(downtime_start, timestep + 1)
            }
            downtime[("on", unit_index, timestep)] = 1.0
            add(
                downtime,
                -math.inf,
                1.0,
                f"min_downtime_i{unit_index}_t{timestep}",
            )

    return MIPProblem(
        instance=parsed,
        variable_keys=variable_keys,
        objective=tuple(objective),
        integrality=tuple(integrality),
        lower_bounds=tuple(lower_bounds),
        upper_bounds=tuple(upper_bounds),
        constraints=tuple(rows),
    )


def decision_variable_key(kind: str, unit: int, timestep: int) -> str:
    """Return publication-compatible tuple-string variable names."""

    return str((kind, int(unit), int(timestep)))


def decode_solution(
    problem: MIPProblem,
    values: Sequence[float],
    *,
    binary_tolerance: float = 1e-6,
    zero_tolerance: float = 1e-9,
) -> dict[str, float | int]:
    """Convert a solver vector into named, JSON-ready decision variables."""

    if len(values) != problem.num_variables:
        raise ValueError(
            f"Expected {problem.num_variables} values, received {len(values)}."
        )
    decoded: dict[str, float | int] = {}
    for index, ((kind, unit, timestep), raw_value) in enumerate(
        zip(problem.variable_keys, values)
    ):
        value = float(raw_value)
        if problem.integrality[index]:
            rounded = round(value)
            if abs(value - rounded) > binary_tolerance or rounded not in (0, 1):
                raise ValueError(
                    f"Non-binary value {value} for {(kind, unit, timestep)}."
                )
            normalized: float | int = int(rounded)
        else:
            normalized = 0.0 if abs(value) <= zero_tolerance else value
        decoded[decision_variable_key(kind, unit, timestep)] = normalized
    return decoded


def solve_with_scipy(
    problem: MIPProblem | UCPInstance | str | Path,
    *,
    time_limit: float | None = 300.0,
    mip_rel_gap: float = 0.0,
    presolve: bool = True,
    display: bool = False,
) -> MIPSolution:
    """Solve a UCP MILP with ``scipy.optimize.milp`` and HiGHS."""

    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_array

    mip = problem if isinstance(problem, MIPProblem) else build_mip(problem)
    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    for row_index, row in enumerate(mip.constraints):
        for column_index, value in row.coefficients.items():
            row_indices.append(row_index)
            column_indices.append(column_index)
            coefficients.append(value)
    matrix = coo_array(
        (coefficients, (row_indices, column_indices)),
        shape=(mip.num_constraints, mip.num_variables),
    ).tocsr()
    constraints = LinearConstraint(
        matrix,
        np.asarray([row.lower_bound for row in mip.constraints], dtype=float),
        np.asarray([row.upper_bound for row in mip.constraints], dtype=float),
    )
    options: dict[str, Any] = {
        "disp": bool(display),
        "presolve": bool(presolve),
        "mip_rel_gap": float(mip_rel_gap),
    }
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    result = milp(
        c=np.asarray(mip.objective, dtype=float),
        integrality=np.asarray(mip.integrality, dtype=np.uint8),
        bounds=Bounds(mip.lower_bounds, mip.upper_bounds),
        constraints=constraints,
        options=options,
    )
    if not result.success or result.x is None or result.fun is None:
        raise RuntimeError(
            f"HiGHS did not solve {mip.instance.name}: "
            f"status={result.status}, message={result.message!s}"
        )
    gap = getattr(result, "mip_gap", None)
    return MIPSolution(
        objective_function_value=float(result.fun),
        decision_variables=decode_solution(mip, result.x),
        raw_values=tuple(float(value) for value in result.x),
        solver_status=int(result.status),
        solver_message=str(result.message),
        mip_gap=None if gap is None else float(gap),
    )


def build_gurobi_model(
    problem: MIPProblem | UCPInstance | str | Path,
    *,
    model_name: str = "publication_ucp",
):
    """Create an equivalent ``gurobipy.Model`` and its ordered variables.

    Returns ``(model, variables, mip_problem)``.  Gurobi is an optional
    dependency and is not required for the SciPy workflow.
    """

    import gurobipy as gp
    from gurobipy import GRB

    mip = problem if isinstance(problem, MIPProblem) else build_mip(problem)
    model = gp.Model(model_name)
    variables = []
    for index, (kind, unit, timestep) in enumerate(mip.variable_keys):
        variable = model.addVar(
            lb=mip.lower_bounds[index],
            ub=mip.upper_bounds[index],
            vtype=GRB.CONTINUOUS if mip.integrality[index] == 0 else GRB.BINARY,
            name=f"{kind}[{unit},{timestep}]",
        )
        variables.append(variable)
    model.setObjective(
        gp.quicksum(
            coefficient * variables[index]
            for index, coefficient in enumerate(mip.objective)
            if coefficient != 0.0
        ),
        GRB.MINIMIZE,
    )
    for row in mip.constraints:
        expression = gp.quicksum(
            coefficient * variables[index]
            for index, coefficient in row.coefficients.items()
        )
        if row.lower_bound == row.upper_bound:
            model.addConstr(expression == row.lower_bound, name=row.name)
        else:
            if math.isfinite(row.lower_bound):
                model.addConstr(
                    expression >= row.lower_bound,
                    name=f"{row.name}_lower",
                )
            if math.isfinite(row.upper_bound):
                model.addConstr(
                    expression <= row.upper_bound,
                    name=f"{row.name}_upper",
                )
    model.update()
    return model, tuple(variables), mip


def solution_from_gurobi(
    problem: MIPProblem,
    variables: Sequence[Any],
    model: Any,
) -> MIPSolution:
    """Normalize an already-optimized Gurobi model into :class:`MIPSolution`."""

    values = tuple(float(variable.X) for variable in variables)
    gap = getattr(model, "MIPGap", None)
    return MIPSolution(
        objective_function_value=float(model.ObjVal),
        decision_variables=decode_solution(problem, values),
        raw_values=values,
        solver_status=int(model.Status),
        solver_message="Gurobi optimization completed",
        mip_gap=None if gap is None else float(gap),
    )


__all__ = [
    "ConstraintRow",
    "INSTANCE_PATTERN",
    "MIPProblem",
    "MIPSolution",
    "UCPInstance",
    "Unit",
    "build_gurobi_model",
    "build_mip",
    "decision_variable_key",
    "decode_solution",
    "read_instance",
    "solution_from_gurobi",
    "solve_with_scipy",
]
