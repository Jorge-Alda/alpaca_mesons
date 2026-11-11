#!/usr/bin/env python3
"""Pole tracking of D(s) = s K - M^2 along a parameter path (task 2.4).

Run:   .venv/bin/python scripts/track_alp_meson_poles.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_track_alp_meson_poles.py

Conventions as in scripts/solve_alp_meson_spectrum.py (Alda:2026nrd;
memory/alp/shared/finite-alp-meson-mixing.md): Phi = (P_i, a), ALP last,
pole vectors v^(n) with v^(n)T K v^(n) = 1. The LO path inputs are the
neutral 6x6 blocks (P_3, P_8, P_0, P_6, P_7, a) of
derive_alp_lo_neutral_blocks.py, scanned in m_a0^2 (GeV^2) or in
epsilon = f_pi/f_a (dimensionless, f_pi fixed).

Method.
  * At every point the spectrum is obtained with
    diagnose_alp_meson_spectrum.solve_physical_spectrum; an unstable point
    raises UnstableVacuumError (the parameter value is attached as
    `parameter`), so no unphysical spectrum enters a track.
  * Branches are labelled at the first point by increasing pole. Between two
    consecutive points the K-weighted overlaps
        O_bn = |v_b(previous)^T K_mid v^(n)(current)|,  K_mid = (K_prev + K_cur)/2,
    are maximized by a one-to-one assignment (scipy linear_sum_assignment).
    Signs are chosen so that the assigned overlaps are positive (continuous
    vectors). Within an exactly degenerate group of the current point
    (spread below `exact_relative_tolerance` times max|m^2|, e.g. the K0
    doublet for flavour-diagonal sources) the basis is first rotated to
    match the previous vectors (orthogonal Procrustes rotation), because the
    solver's basis there is arbitrary.
  * The step is bisected (down to `minimal_step`) when (a) an assigned overlap
    is below `overlap_threshold`; (b) two branches exchange their order and
    their states are coupled; or (c) the largest pole change in the step
    exceeds `gap_step_factor` times the smallest gap between neighbouring
    poles that both carry an ALP component (adaptive refinement near small
    gaps). An exchange of order is accepted without refinement only when the
    two states are decoupled: |v_b1^T (s dK - dM^2) v_b2| <=
    `decoupling_tolerance` |s dK - dM^2| with the previous vectors, the change
    dK, dM^2 over the step, and s the mean of the two poles. This is the case
    for a state with zero ALP component crossing the ALP-like state (e.g.
    the exact m_K0^2 pole, or mesons that do not mix with the ALP). An
    exchange that persists at the minimal step between coupled states is
    recorded as "unresolved".
  * The ALP-like branch at each point is the branch with the largest |v_a|.
  * Avoided crossings: interior local minima of the gap between two
    branches that both carry an ALP component and are neighbours in the
    order of the ALP-carrying poles. The minimum is refined with scipy's
    bounded scalar minimization of the gap between the two poles that best
    overlap the two branch vectors at the sampled minimum.

Units: parameter as given (GeV^2 for m_a0^2, dimensionless for epsilon);
poles and gaps in GeV^2. Inputs are illustrative, not fits.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.optimize
import sympy as sp

from derive_alp_lo_neutral_blocks import NEUTRAL_FIELD_LABELS, closed_form_blocks, illustrative_lo_inputs
from diagnose_alp_meson_spectrum import (
    UnstableVacuumError,
    isolated_two_level_expectation,
    isolated_two_level_inputs,
    solve_physical_spectrum,
)


@dataclass
class TrackedPoles:
    """Poles and K-normalized vectors labelled by branch along a parameter path.

    poles[k, b] is the pole of branch b at parameter_values[k];
    pole_vectors[k, :, b] its vector. `crossings` lists accepted order
    exchanges between branches with their classification.
    """

    parameter_values: np.ndarray
    poles: np.ndarray
    pole_vectors: np.ndarray
    kinetic: list
    mass_squared: list
    field_labels: tuple
    alp_index: int
    crossings: list = field(default_factory=list)
    refinement_count: int = 0

    @property
    def branch_count(self) -> int:
        return self.poles.shape[1]

    def alp_components(self):
        """|v_a| for every point (rows) and branch (columns)."""
        return np.abs(self.pole_vectors[:, self.alp_index, :])

    def alp_like_branch(self):
        """Branch with the largest |v_a| at each point."""
        return np.argmax(self.alp_components(), axis=1)

    def unresolved_crossings(self):
        return [crossing for crossing in self.crossings if crossing["kind"] == "unresolved"]

    def order_exchanges_between_alp_carrying_branches(self, alp_component_tolerance=1e-8):
        """(point index, branch, branch) where two branches carrying an ALP component change order.

        Empty when no labels cross between coupled states (the requirement
        for a scan through an avoided crossing).
        """
        carrying = self.alp_components() > alp_component_tolerance
        exchanges = []
        for index in range(len(self.parameter_values) - 1):
            for first in range(self.branch_count):
                for second in range(first + 1, self.branch_count):
                    if not (carrying[index, first] and carrying[index, second]
                            and carrying[index + 1, first] and carrying[index + 1, second]):
                        continue
                    before = self.poles[index, first] - self.poles[index, second]
                    after = self.poles[index + 1, first] - self.poles[index + 1, second]
                    if before * after < 0:
                        exchanges.append((index, first, second))
        return exchanges

    def alp_like_switches(self):
        """(parameter, from branch, to branch) where the ALP-like branch changes."""
        alp_like = self.alp_like_branch()
        return [(float(self.parameter_values[index]), int(alp_like[index - 1]), int(alp_like[index]))
                for index in range(1, len(alp_like)) if alp_like[index] != alp_like[index - 1]]

    def minimal_gaps(self):
        """For each pair of branches: smallest |m_b1^2 - m_b2^2| on the sampled points and where."""
        result = {}
        for first in range(self.branch_count):
            for second in range(first + 1, self.branch_count):
                gaps = np.abs(self.poles[:, first] - self.poles[:, second])
                position = int(np.argmin(gaps))
                result[(first, second)] = (float(gaps[position]), float(self.parameter_values[position]))
        return result


def _scale(poles):
    return max(1e-300, float(np.max(np.abs(poles))))


def _solve(build_blocks, parameter, field_labels):
    kinetic, mass_squared = build_blocks(parameter)
    try:
        spectrum, _ = solve_physical_spectrum(kinetic, mass_squared, field_labels=field_labels)
    except UnstableVacuumError as error:
        error.parameter = parameter
        error.args = (f"at parameter {parameter}: {error.args[0]}",)
        raise
    return (np.asarray(kinetic, dtype=float), np.asarray(mass_squared, dtype=float),
            np.array(spectrum.pole_squared_masses, dtype=float), np.array(spectrum.pole_vectors, dtype=float))


def _align_degenerate_groups(poles, vectors, previous_vectors, metric, exact_relative_tolerance):
    """Rotate the basis of each exactly degenerate group to best match the previous vectors."""
    vectors = vectors.copy()
    tolerance = exact_relative_tolerance * _scale(poles)
    groups, current = [], [0]
    for index in range(1, len(poles)):
        if poles[index] - poles[current[-1]] <= tolerance:
            current.append(index)
        else:
            groups.append(current)
            current = [index]
    groups.append(current)
    for group in groups:
        if len(group) < 2:
            continue
        basis = vectors[:, group]
        overlaps = basis.T @ metric @ previous_vectors
        selected = np.argsort(-np.linalg.norm(overlaps, axis=0))[:len(group)]
        left, _, right = np.linalg.svd(overlaps[:, selected])
        vectors[:, group] = basis @ (left @ right)
    return vectors


def _match(previous, current, alp_index, settings):
    """Assign current poles to branches; returns (accept, ordered poles, ordered vectors, crossings)."""
    previous_kinetic, previous_mass, previous_poles, previous_vectors = previous
    kinetic, mass_squared, poles, vectors = current
    metric = (previous_kinetic + kinetic) / 2
    vectors = _align_degenerate_groups(poles, vectors, previous_vectors, metric, settings["exact_relative_tolerance"])
    signed = previous_vectors.T @ metric @ vectors
    rows, columns = scipy.optimize.linear_sum_assignment(-np.abs(signed))
    order = np.empty(len(rows), dtype=int)
    order[rows] = columns
    ordered_poles = poles[order]
    ordered_vectors = vectors[:, order] * np.sign(signed[rows, columns])[np.argsort(rows)]
    smallest_overlap = float(np.min(np.abs(signed[rows, columns])))
    scale = _scale(np.concatenate([poles, previous_poles]))
    exact = settings["exact_relative_tolerance"] * scale
    change_kinetic, change_mass = kinetic - previous_kinetic, mass_squared - previous_mass
    crossings, coupled_crossing = [], False
    for first in range(len(order)):
        for second in range(len(order)):
            if not (previous_poles[second] - previous_poles[first] > exact
                    and ordered_poles[first] - ordered_poles[second] > exact):
                continue
            squared_momentum = (ordered_poles[first] + ordered_poles[second]) / 2
            change = squared_momentum * change_kinetic - change_mass
            coupling = abs(previous_vectors[:, first] @ change @ previous_vectors[:, second])
            decoupled = coupling <= settings["decoupling_tolerance"] * max(1e-300, np.linalg.norm(change))
            coupled_crossing |= not decoupled
            crossings.append({"branches": (first, second), "coupling": float(coupling),
                              "kind": "decoupled" if decoupled else "unresolved"})
    alp_carrying = np.abs(ordered_vectors[alp_index, :]) > settings["alp_component_tolerance"]
    previous_carrying = np.abs(previous_vectors[alp_index, :]) > settings["alp_component_tolerance"]
    gaps = []
    for poles_here, carrying in ((previous_poles, previous_carrying), (ordered_poles, alp_carrying)):
        values = np.sort(poles_here[carrying])
        if len(values) > 1:
            gaps.append(float(np.min(np.diff(values))))
    smallest_gap = min(gaps) if gaps else np.inf
    largest_change = float(np.max(np.abs(ordered_poles - previous_poles)))
    accept = (smallest_overlap >= settings["overlap_threshold"] and not coupled_crossing
              and largest_change <= settings["gap_step_factor"] * smallest_gap)
    return accept, ordered_poles, ordered_vectors, crossings


def track_poles(parameter_values, build_blocks, field_labels=None, alp_index=-1, overlap_threshold=0.8,
                gap_step_factor=0.5, minimal_step=None, exact_relative_tolerance=1e-9,
                decoupling_tolerance=1e-10, alp_component_tolerance=1e-8, maximal_points=100000):
    """Track all poles along `parameter_values` (monotonic) with adaptive bisection; see module docstring.

    build_blocks(parameter) -> (K, M^2) as arrays. Returns TrackedPoles.
    """
    parameter_values = [float(value) for value in parameter_values]
    if minimal_step is None:
        minimal_step = 1e-9 * abs(parameter_values[-1] - parameter_values[0])
    settings = {"overlap_threshold": overlap_threshold, "gap_step_factor": gap_step_factor,
                "exact_relative_tolerance": exact_relative_tolerance, "decoupling_tolerance": decoupling_tolerance,
                "alp_component_tolerance": alp_component_tolerance}
    first = _solve(build_blocks, parameter_values[0], field_labels)
    size = first[0].shape[0]
    alp = alp_index % size
    accepted_parameters, accepted_poles, accepted_vectors = [parameter_values[0]], [first[2]], [first[3]]
    kinetic_list, mass_list = [first[0]], [first[1]]
    previous, crossings, refinements = first, [], 0
    for target in parameter_values[1:]:
        pending = [target]
        while pending:
            if len(accepted_parameters) >= maximal_points:
                raise RuntimeError("maximal number of points reached while tracking poles")
            candidate = pending[-1]
            current = _solve(build_blocks, candidate, field_labels)
            accept, poles, vectors, step_crossings = _match(previous, current, alp, settings)
            step = abs(candidate - accepted_parameters[-1])
            if not accept and step > minimal_step:
                pending.append((accepted_parameters[-1] + candidate) / 2)
                refinements += 1
                continue
            pending.pop()
            for crossing in step_crossings:
                crossings.append({**crossing, "interval": (accepted_parameters[-1], candidate)})
            previous = (current[0], current[1], poles, vectors)
            accepted_parameters.append(candidate)
            accepted_poles.append(poles)
            accepted_vectors.append(vectors)
            kinetic_list.append(current[0])
            mass_list.append(current[1])
    return TrackedPoles(
        parameter_values=np.array(accepted_parameters), poles=np.array(accepted_poles),
        pole_vectors=np.array(accepted_vectors), kinetic=kinetic_list, mass_squared=mass_list,
        field_labels=tuple(field_labels) if field_labels is not None else tuple(
            [f"P_{index}" for index in range(size - 1)] + ["a"]),
        alp_index=alp, crossings=crossings, refinement_count=refinements)


def avoided_crossings(tracked: TrackedPoles, build_blocks, alp_component_tolerance=1e-8, refine=True):
    """Interior gap minima between two ALP-carrying branches, refined by bounded minimization.

    A pair of branches that both carry an ALP component is examined if the
    two are neighbours in the order of the ALP-carrying poles at the sampled
    minimum, or if the ALP-like branch passes from one to the other somewhere
    on the path (`alp_like_switch`, the ALP crossing a meson; other weakly
    coupled poles may lie in between).
    """
    alp_components = tracked.alp_components()
    switches = {frozenset((start, end)) for _, start, end in tracked.alp_like_switches()}
    found = []
    points = len(tracked.parameter_values)
    for index in range(1, points - 1):
        carrying = [branch for branch in range(tracked.branch_count) if alp_components[index, branch] > alp_component_tolerance]
        ordered = sorted(carrying, key=lambda branch: tracked.poles[index, branch])
        pairs = [(lower, upper) for position, lower in enumerate(ordered) for upper in ordered[position + 1:]
                 if upper == ordered[position + 1] or frozenset((lower, upper)) in switches]
        for lower, upper in pairs:
            gaps = np.abs(tracked.poles[index - 1:index + 2, upper] - tracked.poles[index - 1:index + 2, lower])
            if not (gaps[1] < gaps[0] and gaps[1] <= gaps[2]):
                continue
            entry = {"branches": (int(lower), int(upper)), "parameter": float(tracked.parameter_values[index]),
                     "minimal_gap": float(gaps[1]), "sampled_minimal_gap": float(gaps[1]),
                     "alp_components": (float(alp_components[index, lower]), float(alp_components[index, upper])),
                     "alp_like_switch": frozenset((lower, upper)) in switches}
            if refine:
                reference = tracked.pole_vectors[index][:, [lower, upper]]
                metric = tracked.kinetic[index]

                def gap(parameter):
                    _, _, poles, vectors = _solve(build_blocks, parameter, tracked.field_labels)
                    overlaps = np.abs(reference.T @ metric @ vectors)
                    first = int(np.argmax(overlaps[0]))
                    overlaps[1, first] = -1.0
                    second = int(np.argmax(overlaps[1]))
                    return abs(poles[second] - poles[first])

                bounds = (tracked.parameter_values[index - 1], tracked.parameter_values[index + 1])
                result = scipy.optimize.minimize_scalar(gap, bounds=bounds, method="bounded",
                                                        options={"xatol": 1e-12 * max(1.0, abs(bounds[1] - bounds[0]))})
                if result.fun <= gaps[1]:
                    entry["parameter"], entry["minimal_gap"] = float(result.x), float(result.fun)
            found.append(entry)
    return found


# ------------------------------------------------------------- LO scans


def lo_blocks_builder(inputs, parameter_name):
    """build_blocks(parameter) for the LO 6x6 blocks with `parameter_name` varied.

    parameter_name: "bare_alp_mass_squared" (m_a0^2 in GeV^2) or "ratio"
    (epsilon = f_pi/f_a with f_pi fixed).
    """

    def build(parameter):
        varied = dict(inputs)
        if parameter_name == "ratio":
            varied["alp_decay_constant"] = varied["pion_decay_constant"] / sp.nsimplify(parameter)
        else:
            varied[parameter_name] = sp.nsimplify(parameter)
        kinetic, mass_squared = closed_form_blocks(**varied)
        return np.array(kinetic.evalf(20), dtype=float), np.array(mass_squared.evalf(20), dtype=float)

    return build


def fast_lo_blocks_builder(inputs, parameter_name):
    """Same as lo_blocks_builder but evaluates the exact blocks once and rescales in floats.

    K_Pa, M^2_Pa scale as epsilon and K_aa - 1, M^2_aa - m_a0^2 as epsilon^2;
    m_a0^2 enters only M^2_aa.
    """
    base_ratio = float(inputs["pion_decay_constant"] / inputs["alp_decay_constant"])
    bare = float(inputs["bare_alp_mass_squared"])
    kinetic, mass_squared = closed_form_blocks(**inputs)
    kinetic = np.array(kinetic.evalf(30), dtype=float)
    mass_squared = np.array(mass_squared.evalf(30), dtype=float)

    def build(parameter):
        new_kinetic, new_mass = kinetic.copy(), mass_squared.copy()
        if parameter_name == "ratio":
            factor = parameter / base_ratio
            new_kinetic[:-1, -1] *= factor
            new_kinetic[-1, :-1] *= factor
            new_mass[:-1, -1] *= factor
            new_mass[-1, :-1] *= factor
            new_kinetic[-1, -1] = 1 + (kinetic[-1, -1] - 1) * factor**2
            new_mass[-1, -1] = bare + (mass_squared[-1, -1] - bare) * factor**2
        else:
            new_mass[-1, -1] = mass_squared[-1, -1] - bare + parameter
        return new_kinetic, new_mass

    return build


def scan_summary(tracked, build_blocks):
    """Crossings, ALP-like branch changes and avoided crossings of a track."""
    changes = tracked.alp_like_switches()
    return {
        "points": len(tracked.parameter_values),
        "refinements": tracked.refinement_count,
        "decoupled crossings": [crossing for crossing in tracked.crossings if crossing["kind"] == "decoupled"],
        "unresolved crossings": tracked.unresolved_crossings(),
        "ALP-like branch changes (parameter, from, to)": changes,
        "avoided crossings": avoided_crossings(tracked, build_blocks),
        "order exchanges between ALP-carrying branches": tracked.order_exchanges_between_alp_carrying_branches(),
    }


def _with_ratio(inputs, ratio):
    return dict(inputs, alp_decay_constant=inputs["pion_decay_constant"] / ratio)


def derive():
    """Scans of the LO system in m_a0^2 (isolated two-level, flavour-diagonal, neutral) and in epsilon."""
    results = {}
    isolated = isolated_two_level_inputs()
    builder = fast_lo_blocks_builder(isolated, "bare_alp_mass_squared")
    tracked = track_poles(np.linspace(0.002, 1.5, 31), builder, field_labels=NEUTRAL_FIELD_LABELS)
    summary = scan_summary(tracked, builder)
    expected = float(isolated_two_level_expectation(isolated)["splitting"])
    crossings = summary["avoided crossings"]
    summary["checks"] = {
        "no unresolved crossings": not summary["unresolved crossings"],
        "no order exchange between ALP-carrying branches": not summary["order exchanges between ALP-carrying branches"],
        "one avoided crossing (with P_3)": len(crossings) == 1 and crossings[0]["alp_like_switch"],
        "its minimal gap equals 2|M^2_P'a,can|": len(crossings) == 1 and abs(crossings[0]["minimal_gap"] - expected) < 1e-9,
        "it is at the bare mass where M^2_aa,can = m_pi^2": len(crossings) == 1 and abs(
            crossings[0]["parameter"] - float(isolated["bare_alp_mass_squared"])) < 1e-6,
        "ALP-like branch crosses the decoupled K0 doublet, eta and eta' (4 decoupled crossings)":
            len(summary["decoupled crossings"]) == 4,
    }
    results["m_a0^2 scan, isolated two-level inputs (epsilon = 1/10)"] = summary
    for source_structure in ("flavor_diagonal", "neutral"):
        for ratio in (sp.Rational(1, 20), None):
            inputs = illustrative_lo_inputs(source_structure)
            if ratio is not None:
                inputs = _with_ratio(inputs, ratio)
            label = f"epsilon = {ratio if ratio is not None else 'illustrative 23/50'}"
            builder = fast_lo_blocks_builder(inputs, "bare_alp_mass_squared")
            tracked = track_poles(np.linspace(0.002 if ratio is not None else 0.2, 1.5, 31), builder,
                                  field_labels=NEUTRAL_FIELD_LABELS)
            summary = scan_summary(tracked, builder)
            summary["checks"] = {
                "no unresolved crossings": not summary["unresolved crossings"],
                "no order exchange between ALP-carrying branches": not summary["order exchanges between ALP-carrying branches"],
                "avoided crossings have positive gaps": all(entry["minimal_gap"] > 0 for entry in summary["avoided crossings"]),
                "every ALP-like branch switch has an avoided crossing": all(
                    any(set(entry["branches"]) == {start, end} for entry in summary["avoided crossings"])
                    for _, start, end in summary["ALP-like branch changes (parameter, from, to)"]),
            }
            results[f"m_a0^2 scan, {source_structure} sources, {label}"] = summary
        inputs = illustrative_lo_inputs(source_structure)
        ratio_builder = fast_lo_blocks_builder(inputs, "ratio")
        tracked = track_poles(np.linspace(0.005, 0.5, 34), ratio_builder, field_labels=NEUTRAL_FIELD_LABELS)
        summary = scan_summary(tracked, ratio_builder)
        alp_like = tracked.alp_like_branch()
        summary["checks"] = {
            "no unresolved crossings": not summary["unresolved crossings"],
            "no order exchange between ALP-carrying branches": not summary["order exchanges between ALP-carrying branches"],
            "ALP-like pole within 1e-3 GeV^2 of m_a0^2 at epsilon = 0.005": abs(
                tracked.poles[0, alp_like[0]] - float(inputs["bare_alp_mass_squared"])) < 1e-3,
            "|v_a| of the ALP-like pole within 1e-2 of 1 at epsilon = 0.005": abs(
                tracked.alp_components()[0, alp_like[0]] - 1) < 1e-2,
        }
        results[f"epsilon scan, {source_structure} sources, m_a0^2 = 3/10"] = summary
    return results


def main():
    for name, summary in derive().items():
        print(name)
        for label, passed in summary["checks"].items():
            print(f"  [{'pass' if passed else 'FAIL'}] {label}")
        print(f"  points {summary['points']}, bisections {summary['refinements']}")
        print(f"  ALP-like branch changes: {summary['ALP-like branch changes (parameter, from, to)']}")
        for entry in summary["avoided crossings"]:
            print(f"  avoided crossing branches {entry['branches']} at {entry['parameter']:.9g}: "
                  f"minimal gap {entry['minimal_gap']:.6g} GeV^2, |v_a| "
                  f"({entry['alp_components'][0]:.4f}, {entry['alp_components'][1]:.4f}), "
                  f"ALP-like switch {entry['alp_like_switch']}")
        print(f"  decoupled crossings: {[(entry['branches'], entry['interval']) for entry in summary['decoupled crossings']]}")


if __name__ == "__main__":
    main()
