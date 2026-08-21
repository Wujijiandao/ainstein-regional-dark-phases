#!/usr/bin/env python3
"""History-dependent basin activation audit for the regional dark-phase prototype.

This module implements only the simulation-facing bookkeeping layer.  It does NOT
claim a microscopic transition law.  A basin is allowed to activate when the
regional candidate criterion Gamma >= 1 is met while it is expanding and
single-stream.  In the default *absorbing* prototype, activation is retained as
history information thereafter.  At every snapshot, maximal active nodes in the
laminar hierarchy are selected, producing a unique antichain and preventing
nested double counting.

The purpose of the synthetic trajectory is to demonstrate two structural facts:
(1) an instantaneous Gamma classifier can chatter around the threshold, whereas
    a history-dependent phase label need not;
(2) for a fixed laminar hierarchy with absorbing activation, the union of active
    domains is monotone, so its filling fraction is nondecreasing.  Combined with
    a fixed positive intrinsic V-phase energy density, this falls under the
    manuscript's phantom-or-Lambda residual sign theorem.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple
import csv

from analysis_basin_selector import Basin, ancestors, is_antichain

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class SnapshotState:
    time: float
    gamma: Dict[str, float]
    expanding: Dict[str, bool]
    single_stream: Dict[str, bool]


def hierarchy_template() -> List[Basin]:
    """Fixed laminar hierarchy for a deterministic history audit."""
    return [
        Basin("U", None, 1.00, gamma=0.0),
        Basin("A", "U", 0.55, gamma=0.0),
        Basin("A1", "A", 0.25, gamma=0.0),
        Basin("A2", "A", 0.30, gamma=0.0),
        Basin("B", "U", 0.45, gamma=0.0),
        Basin("B1", "B", 0.20, gamma=0.0),
        Basin("B2", "B", 0.25, gamma=0.0),
    ]


def synthetic_history() -> List[SnapshotState]:
    """A threshold-crossing history with deliberate near-threshold chatter."""
    names = [b.name for b in hierarchy_template()]
    values = [
        # U,    A,    A1,   A2,   B,    B1,   B2
        [0.2, 0.80, 0.90, 0.85, 0.70, 0.95, 0.75],
        [0.2, 0.90, 1.10, 0.92, 0.76, 1.05, 0.82],
        [0.2, 0.95, 0.97, 1.08, 0.82, 0.98, 0.88],
        [0.2, 1.05, 1.03, 1.11, 0.90, 1.02, 0.94],
        [0.2, 0.97, 0.96, 1.04, 1.03, 0.97, 1.01],
        [0.2, 1.10, 1.02, 1.06, 0.99, 1.03, 0.98],
    ]
    out: List[SnapshotState] = []
    for i, row in enumerate(values):
        g = dict(zip(names, row))
        expanding = {n: (n != "U") for n in names}
        single = {n: True for n in names}
        # Root is deliberately excluded from activation: it represents the box.
        out.append(SnapshotState(float(i), g, expanding, single))
    return out


def instantaneous_active(state: SnapshotState, basins: Sequence[Basin]) -> Set[str]:
    return {
        b.name for b in basins
        if b.parent is not None
        and state.expanding.get(b.name, False)
        and state.single_stream.get(b.name, False)
        and state.gamma.get(b.name, float("-inf")) >= 1.0
    }


def absorbing_active_history(states: Sequence[SnapshotState], basins: Sequence[Basin]) -> List[Set[str]]:
    """Activate on first eligibility crossing and retain activation thereafter.

    This is intentionally the strongest-hysteresis limiting prototype.  A real
    phase theory may include deactivation or node death/merger rules; those are
    exactly the extra ingredients required if observations reject monotone filling.
    """
    active: Set[str] = set()
    history: List[Set[str]] = []
    for s in states:
        active |= instantaneous_active(s, basins)
        history.append(set(active))
    return history


def maximal_from_active(active: Set[str], basins: Sequence[Basin]) -> List[Basin]:
    by_name = {b.name: b for b in basins}
    selected: List[Basin] = []
    for b in basins:
        if b.name not in active:
            continue
        if not any(a in active for a in ancestors(b.name, by_name)):
            selected.append(b)
    return sorted(selected, key=lambda b: b.name)


def union_volume_of_active(active: Set[str], basins: Sequence[Basin]) -> float:
    """For a laminar hierarchy, maximal active nodes exactly represent the union."""
    return sum(b.volume for b in maximal_from_active(active, basins))


def switch_count(series: Sequence[Set[str]]) -> int:
    """Count per-node Boolean state flips across snapshots."""
    names = sorted(set().union(*series))
    if not series:
        return 0
    count = 0
    all_names = sorted({b.name for b in hierarchy_template() if b.parent is not None})
    prev = {n: n in series[0] for n in all_names}
    for s in series[1:]:
        now = {n: n in s for n in all_names}
        count += sum(now[n] != prev[n] for n in all_names)
        prev = now
    return count


def trajectory_rows() -> List[dict]:
    basins = hierarchy_template()
    states = synthetic_history()
    inst = [instantaneous_active(s, basins) for s in states]
    hist = absorbing_active_history(states, basins)
    rows = []
    for s, ia, ha in zip(states, inst, hist):
        imax = maximal_from_active(ia, basins)
        hmax = maximal_from_active(ha, basins)
        rows.append({
            "time": s.time,
            "instant_active": ";".join(sorted(ia)),
            "instant_maximal": ";".join(b.name for b in imax),
            "phi_instant": union_volume_of_active(ia, basins),
            "history_active": ";".join(sorted(ha)),
            "history_maximal": ";".join(b.name for b in hmax),
            "phi_history": union_volume_of_active(ha, basins),
        })
    return rows


def write_outputs() -> None:
    outdir = ROOT / "results"
    outdir.mkdir(exist_ok=True)
    rows = trajectory_rows()
    csv_path = outdir / "basin_dynamics_trajectory.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    basins = hierarchy_template()
    states = synthetic_history()
    inst = [instantaneous_active(s, basins) for s in states]
    hist = absorbing_active_history(states, basins)
    ph = [union_volume_of_active(a, basins) for a in hist]
    summary = [
        "History-dependent basin-state audit",
        "prototype=absorbing activation after first Gamma>=1 eligibility crossing",
        f"instantaneous_switch_count={switch_count(inst)}",
        f"history_switch_count={switch_count(hist)}",
        f"phi_history={','.join(f'{x:.6f}' for x in ph)}",
        f"phi_history_monotone={all(b >= a - 1e-15 for a,b in zip(ph,ph[1:]))}",
        f"final_maximal_active={[b.name for b in maximal_from_active(hist[-1], basins)]}",
        f"final_antichain={is_antichain(maximal_from_active(hist[-1], basins), basins)}",
        "interpretation=algorithmic limiting case only; deactivation/interface physics is not derived",
    ]
    (outdir / "basin_dynamics_summary.txt").write_text("\n".join(summary)+"\n", encoding="utf-8")


def main() -> None:
    write_outputs()
    print((ROOT / "results" / "basin_dynamics_summary.txt").read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
