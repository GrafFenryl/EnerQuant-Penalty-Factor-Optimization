"""Standalone QUBO construction and evaluation for the published UCP data.

The module reproduces the formulation used for the publication runs:

* logarithmic generation discretization with step size 1,
* zero unit status before the modelled horizon, and
* minimum up/down-time lookahead windows omit the current timestep (offset zero).

Only the Python standard library is required.  The module deliberately has no
imports from the original benchmark repository, so these utilities can be used
in an independent repository.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence


PathLike = str | Path
Variable = tuple[str, int, int, tuple[int, ...]]


PART_DESCRIPTIONS = {
    "cost matrix": "generation and start-up objective",
    "a": "demand-balance penalty",
    "b": "on/start coupling penalty",
    "c": "minimum up-time penalty",
    "d": "minimum down-time penalty",
    "e": "on/generation-level coupling penalty",
}


@dataclass
class SparseQUBO:
    """A sparse upper-triangular QUBO plus its constant offset.

    ``terms[(u, v)]`` is the coefficient of ``u * v``.  Linear terms use
    ``u == v``.  Variable names are the exact strings used in the result JSON
    files, for example ``("on", 0, 0, ())`` rendered as a string.
    """

    terms: dict[tuple[str, str], float] = field(default_factory=dict)
    offset: float = 0.0
    name: str = "QUBO"

    def add_term(self, u: str, v: str, coefficient: float) -> None:
        """Add one linear or quadratic coefficient in canonical order."""

        key = (u, v) if u <= v else (v, u)
        self.terms[key] = self.terms.get(key, 0.0) + float(coefficient)

    def add_offset(self, value: float) -> None:
        """Add a constant contribution."""

        self.offset += float(value)

    def add_scaled(self, other: "SparseQUBO", factor: float = 1.0) -> None:
        """Add ``factor * other`` in place."""

        factor = float(factor)
        self.offset += factor * other.offset
        for (u, v), coefficient in other.terms.items():
            self.add_term(u, v, factor * coefficient)

    def scaled(self, factor: float, *, name: str | None = None) -> "SparseQUBO":
        """Return a scaled copy."""

        result = SparseQUBO(name=name or self.name)
        result.add_scaled(self, factor)
        return result

    @property
    def variables(self) -> tuple[str, ...]:
        """Return all variables in deterministic order."""

        return tuple(sorted({variable for pair in self.terms for variable in pair}))

    @property
    def linear(self) -> dict[str, float]:
        """Return diagonal coefficients as ``variable -> coefficient``."""

        return {u: value for (u, v), value in self.terms.items() if u == v}

    @property
    def quadratic(self) -> dict[tuple[str, str], float]:
        """Return off-diagonal coefficients."""

        return {(u, v): value for (u, v), value in self.terms.items() if u != v}

    def evaluate(self, state: Mapping[str, int | float], *, strict: bool = True) -> float:
        """Evaluate a state, including the constant offset.

        With ``strict=True`` (the default), the state must contain exactly the
        variables used by this QUBO.  This catches mismatched instance/result
        files early.
        """

        variables = set(self.variables)
        state_variables = set(state)
        if strict:
            missing = sorted(variables - state_variables)
            extra = sorted(state_variables - variables)
            if missing or extra:
                raise ValueError(
                    "State/QUBO variable mismatch: "
                    f"missing={missing[:5]}, extra={extra[:5]}."
                )
        else:
            missing = []

        if missing:
            raise ValueError(f"State is missing QUBO variables: {missing[:5]}.")

        energy = self.offset
        for (u, v), coefficient in self.terms.items():
            energy += coefficient * float(state[u]) * float(state[v])
        return float(energy)

    def summary(self) -> dict[str, int | float | str]:
        """Return a compact, notebook-friendly description."""

        return {
            "name": self.name,
            "variables": len(self.variables),
            "linear_terms": len(self.linear),
            "quadratic_terms": len(self.quadratic),
            "offset": self.offset,
        }


def _load_json(value: PathLike | Mapping[str, Any], description: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    path = Path(value)
    try:
        # ``utf-8-sig`` also reads ordinary UTF-8 and handles the BOM present
        # in the publication instance files.
        with path.open("r", encoding="utf-8-sig") as file:
            loaded = json.load(file)
    except OSError as exc:
        raise OSError(f"Could not read {description} {path}: {exc}") from exc
    if not isinstance(loaded, Mapping):
        raise TypeError(f"{description.capitalize()} must contain a JSON object: {path}")
    return loaded


def _load_instance(instance: PathLike | Mapping[str, Any]) -> dict[str, Any]:
    raw = _load_json(instance, "instance")
    required = {"units", "supply", "demand"}
    missing = sorted(required - set(raw))
    if missing:
        raise ValueError(f"Instance is missing fields: {missing}.")

    def nested_object(name: str) -> Mapping[str, Any]:
        value = raw[name]
        if isinstance(value, str):
            value = json.loads(value)
        if not isinstance(value, Mapping):
            raise TypeError(f"Instance field {name!r} must encode a JSON object.")
        return value

    units = nested_object("units")
    demand_object = nested_object("demand")
    supply_object = nested_object("supply")
    unit_fields = (
        "min_power",
        "max_power",
        "min_uptime",
        "min_downtime",
        "start_cost",
        "var_cost",
    )
    missing_unit_fields = sorted(set(unit_fields) - set(units))
    if missing_unit_fields:
        raise ValueError(f"Instance units are missing fields: {missing_unit_fields}.")

    unit_indices = sorted(int(index) for index in units[unit_fields[0]])
    if unit_indices != list(range(len(unit_indices))):
        raise ValueError("Unit indices must be consecutive and start at zero.")
    for field_name in unit_fields:
        indices = sorted(int(index) for index in units[field_name])
        if indices != unit_indices:
            raise ValueError(f"Unit field {field_name!r} uses inconsistent indices.")

    def ordered_values(values: Mapping[str, Any]) -> list[float]:
        return [float(values[str(index)]) for index in sorted(map(int, values))]

    demand_column = demand_object.get("demand")
    supply_column = supply_object.get("supply")
    if not isinstance(demand_column, Mapping) or not isinstance(supply_column, Mapping):
        raise ValueError("Demand and supply fields must contain indexed columns.")
    demand = ordered_values(demand_column)
    supply = ordered_values(supply_column)
    if len(demand) != len(supply):
        raise ValueError("Demand and supply must have the same number of timesteps.")

    return {
        "min_power": ordered_values(units["min_power"]),
        "max_power": ordered_values(units["max_power"]),
        "min_uptime": [int(value) for value in ordered_values(units["min_uptime"])],
        "min_downtime": [int(value) for value in ordered_values(units["min_downtime"])],
        "start_cost": ordered_values(units["start_cost"]),
        "var_cost": ordered_values(units["var_cost"]),
        # The original publication QUBO adapter uses the stored demand directly.
        "demand": demand,
        "supply": supply,
    }


def _variable(kind: str, timestep: int, unit: int, *subindices: int) -> Variable:
    return kind, timestep, unit, tuple(subindices)


def _variable_name(variable: Variable) -> str:
    return str(variable)


def _add_linear(qubo: SparseQUBO, variable: Variable, coefficient: float) -> None:
    name = _variable_name(variable)
    qubo.add_term(name, name, coefficient)


def _add_quadratic(
    qubo: SparseQUBO,
    first: Variable,
    second: Variable,
    coefficient: float,
) -> None:
    qubo.add_term(_variable_name(first), _variable_name(second), coefficient)


def _add_square(
    qubo: SparseQUBO,
    coefficients: Sequence[tuple[Variable, float]],
    constant: float,
) -> None:
    """Add ``(constant + sum(coefficient * variable)) ** 2``."""

    qubo.add_offset(constant * constant)
    for index, (variable, coefficient) in enumerate(coefficients):
        _add_linear(qubo, variable, coefficient * coefficient + 2.0 * constant * coefficient)
        for other, other_coefficient in coefficients[index + 1 :]:
            _add_quadratic(qubo, variable, other, 2.0 * coefficient * other_coefficient)


def build_sub_QUBO(
    instance: PathLike | Mapping[str, Any],
    *,
    step: float = 1.0,
    previous_status: Sequence[int] | None = None,
) -> dict[str, SparseQUBO]:
    """Build the six unweighted publication sub-QUBOs from an instance.

    The returned keys match the publication energy/factor labels:
    ``cost matrix`` (objective), ``a`` (demand), ``b`` (start coupling),
    ``c`` (minimum up time), ``d`` (minimum down time), and ``e``
    (generation-level coupling).

    ``previous_status`` defaults to all-zero status, as in the experiments.
    Minimum up/down-time lookahead windows omit the current timestep (offset
    zero), as in the original run configuration. Starts at the first horizon
    timestep still contribute to constraints on future timesteps.
    """

    data = _load_instance(instance)
    if step <= 0:
        raise ValueError("step must be positive.")
    units = len(data["min_power"])
    timesteps = len(data["demand"])
    if previous_status is None:
        previous = [0] * units
    else:
        previous = [int(value) for value in previous_status]
        if len(previous) != units or any(value not in (0, 1) for value in previous):
            raise ValueError("previous_status must contain one binary value per unit.")

    parts = {
        label: SparseQUBO(name=PART_DESCRIPTIONS[label])
        for label in ("cost matrix", "a", "b", "c", "d", "e")
    }

    # Logarithmic generation discretization, exactly as in the experiments.
    level_volumes: list[list[float]] = []
    for minimum, maximum in zip(data["min_power"], data["max_power"]):
        y = int((maximum - minimum) / step) + 1
        if y < 1:
            raise ValueError("Each unit must satisfy max_power >= min_power.")
        levels = 0 if y == 1 else math.floor(math.log2(y - 1)) + 1
        volumes = [step * (2**level) for level in range(max(0, levels - 1))]
        if y > 1:
            volumes.append(step * (y - 2 ** (levels - 1)))
        level_volumes.append(volumes)

    objective = parts["cost matrix"]
    demand_part = parts["a"]
    discretization = parts["e"]
    for timestep in range(timesteps):
        demand_coefficients: list[tuple[Variable, float]] = []
        for unit in range(units):
            on = _variable("on", timestep, unit)
            start = _variable("start", timestep, unit)
            _add_linear(objective, start, data["start_cost"][unit])
            _add_linear(objective, on, data["var_cost"][unit] * data["min_power"][unit])
            demand_coefficients.append((on, data["min_power"][unit]))

            for level, volume in enumerate(level_volumes[unit]):
                level_variable = _variable("level", timestep, unit, level)
                _add_linear(objective, level_variable, data["var_cost"][unit] * volume)
                demand_coefficients.append((level_variable, volume))

                # level * (1 - on)
                _add_linear(discretization, level_variable, 1.0)
                _add_quadratic(discretization, level_variable, on, -1.0)

        _add_square(demand_part, demand_coefficients, -data["demand"][timestep])

    coupling = parts["b"]
    minimum_up = parts["c"]
    minimum_down = parts["d"]
    for unit in range(units):
        for timestep in range(timesteps):
            on = _variable("on", timestep, unit)
            start = _variable("start", timestep, unit)

            # Start/on coupling.  Previous status is zero in the publication.
            if timestep == 0:
                _add_linear(coupling, on, 1 - previous[unit])
                _add_linear(coupling, start, 1 + 2 * previous[unit])
                _add_quadratic(coupling, start, on, -2.0)
            else:
                on_before = _variable("on", timestep - 1, unit)
                _add_linear(coupling, on, 1.0)
                _add_quadratic(coupling, on, on_before, -1.0)
                _add_linear(coupling, start, 1.0)
                _add_quadratic(coupling, start, on_before, 2.0)
                _add_quadratic(coupling, start, on, -2.0)

            # include_0th_step_in_minupdown=False -> offsets start at 1.
            up_cutoff = min(data["min_uptime"][unit], timesteps - timestep)
            up_offsets = range(1, up_cutoff)
            up_count = max(0, up_cutoff - 1)
            if up_count:
                _add_linear(minimum_up, start, up_count)
                for offset in up_offsets:
                    future_on = _variable("on", timestep + offset, unit)
                    _add_quadratic(minimum_up, start, future_on, -1.0)

            down_cutoff = min(data["min_downtime"][unit], timesteps - timestep)
            down_offsets = range(1, down_cutoff)
            if down_cutoff > 1:
                if timestep == 0:
                    # derived_stop = previous * (1 - on_0)
                    if previous[unit]:
                        for offset in down_offsets:
                            future_on = _variable("on", timestep + offset, unit)
                            _add_linear(minimum_down, future_on, 1.0)
                            _add_quadratic(minimum_down, on, future_on, -1.0)
                else:
                    on_before = _variable("on", timestep - 1, unit)
                    for offset in down_offsets:
                        future_on = _variable("on", timestep + offset, unit)
                        _add_quadratic(minimum_down, on, future_on, -1.0)
                        _add_quadratic(minimum_down, on_before, future_on, 1.0)
                        _add_quadratic(minimum_down, start, future_on, 1.0)

    return parts


def _solution_object(
    result: PathLike | Mapping[str, Any],
    solution: str | Mapping[str, Any],
) -> Mapping[str, Any]:
    if isinstance(solution, Mapping):
        return solution
    result_object = _load_json(result, "result")
    if solution not in result_object:
        available = sorted(key for key in result_object if key != "metadata")
        raise KeyError(f"Unknown solution {solution!r}; choose from {available}.")
    solution_object = result_object[solution]
    if not isinstance(solution_object, Mapping):
        raise TypeError(f"Result entry {solution!r} must be a JSON object.")
    return solution_object


def build_QUBO(
    instance: PathLike | Mapping[str, Any],
    result: PathLike | Mapping[str, Any],
    solution: str | Mapping[str, Any] = "best_repaired_solution_overall",
) -> SparseQUBO:
    """Build the weighted QUBO using one published solution's factors.

    ``instance`` and ``result`` may be paths or already-loaded mappings.
    ``solution`` may be a result key or the solution mapping itself.  Both the
    separate ``c``/``d`` and combined ``c+d`` publication conventions are
    supported.
    """

    solution_object = _solution_object(result, solution)
    factors = solution_object.get("penalty_factors")
    if not isinstance(factors, Mapping):
        raise ValueError("The selected solution has no penalty_factors mapping.")

    required = {"cost matrix", "a", "b", "e"}
    missing = sorted(required - set(factors))
    if missing:
        raise ValueError(f"Penalty factors are missing required labels: {missing}.")
    has_separate = "c" in factors and "d" in factors
    has_combined = "c+d" in factors
    if has_separate == has_combined:
        raise ValueError("Penalty factors must use either c and d, or the combined c+d.")
    allowed = required | ({"c", "d"} if has_separate else {"c+d"})
    extra = sorted(set(factors) - allowed)
    if extra:
        raise ValueError(f"Unsupported penalty-factor labels: {extra}.")

    parts = build_sub_QUBO(instance)
    full = SparseQUBO(name="weighted publication QUBO")
    full.add_scaled(parts["cost matrix"], float(factors["cost matrix"]))
    full.add_scaled(parts["a"], float(factors["a"]))
    full.add_scaled(parts["b"], float(factors["b"]))
    full.add_scaled(parts["e"], float(factors["e"]))
    if has_separate:
        full.add_scaled(parts["c"], float(factors["c"]))
        full.add_scaled(parts["d"], float(factors["d"]))
    else:
        full.add_scaled(parts["c"], float(factors["c+d"]))
        full.add_scaled(parts["d"], float(factors["c+d"]))
    return full


def evaluate_solution(
    qubo: SparseQUBO,
    result: PathLike | Mapping[str, Any],
    solution: str | Mapping[str, Any] = "best_repaired_solution_overall",
) -> float:
    """Read a published solution and evaluate it on ``qubo``."""

    if not isinstance(qubo, SparseQUBO):
        raise TypeError("qubo must be a SparseQUBO returned by this module.")
    solution_object = _solution_object(result, solution)
    decision_variables = solution_object.get("decision_variables")
    if not isinstance(decision_variables, Mapping):
        raise ValueError("The selected solution has no decision_variables mapping.")
    return qubo.evaluate(decision_variables)


__all__ = [
    "PART_DESCRIPTIONS",
    "SparseQUBO",
    "build_sub_QUBO",
    "build_QUBO",
    "evaluate_solution",
]
