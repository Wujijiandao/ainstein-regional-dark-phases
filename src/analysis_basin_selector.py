#!/usr/bin/env python3
"""Topological basin-antichain selector for the regional dark-phase prototype.

This module does not infer real cosmic-web basins from observations. It implements
and regression-tests the purely combinatorial part of the manuscript's domain rule:
for a laminar basin hierarchy, keep only *maximal eligible* basins. This produces
an antichain and therefore prevents recursive double counting of nested domains.

The physical V-phase eligibility predicate is represented abstractly here by three
fields attached to a basin: expanding, single_stream, and gamma >= 1.  In a real
calculation these would be measured from a reconstructed/simulated regional state.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Basin:
    name: str
    parent: Optional[str]
    volume: float
    gamma: float
    expanding: bool = True
    single_stream: bool = True


def is_v_eligible(b: Basin) -> bool:
    """Idealized V-phase eligibility used in the manuscript prototype."""
    return b.expanding and b.single_stream and b.gamma >= 1.0


def ancestors(name: str, basins: Dict[str, Basin]) -> List[str]:
    out: List[str] = []
    p = basins[name].parent
    while p is not None:
        out.append(p)
        p = basins[p].parent
    return out


def maximal_eligible(basins: Sequence[Basin]) -> List[Basin]:
    """Return eligible nodes having no eligible strict ancestor.

    In a rooted laminar hierarchy, this is exactly the set of maximal eligible
    domains under set inclusion.  It is unique and forms an antichain.
    """
    by_name = {b.name: b for b in basins}
    eligible = {b.name for b in basins if is_v_eligible(b)}
    selected: List[Basin] = []
    for b in basins:
        if b.name not in eligible:
            continue
        if not any(a in eligible for a in ancestors(b.name, by_name)):
            selected.append(b)
    return sorted(selected, key=lambda x: x.name)


def is_antichain(nodes: Sequence[Basin], basins: Sequence[Basin]) -> bool:
    by_name = {b.name: b for b in basins}
    names = {b.name for b in nodes}
    for n in names:
        if any(a in names for a in ancestors(n, by_name)):
            return False
    return True


def naive_eligible_volume(basins: Sequence[Basin]) -> float:
    """Sum of all eligible node volumes; intentionally double-counts nesting."""
    return sum(b.volume for b in basins if is_v_eligible(b))


def antichain_volume(basins: Sequence[Basin]) -> float:
    """Volume of maximal eligible nodes for a laminar hierarchy."""
    return sum(b.volume for b in maximal_eligible(basins))


def synthetic_hierarchy() -> List[Basin]:
    """A deterministic nested example used only to audit non-recursion.

    Root U is not eligible. A is eligible and contains eligible A1/A2, so the
    antichain rule keeps A and retires its nested children. B is ineligible but
    contains eligible B1, so B1 remains active.
    """
    return [
        Basin("U", None, 1.00, gamma=0.55, expanding=True, single_stream=True),
        Basin("A", "U", 0.55, gamma=1.20),
        Basin("A1", "A", 0.25, gamma=1.40),
        Basin("A2", "A", 0.30, gamma=1.10),
        Basin("B", "U", 0.45, gamma=0.80),
        Basin("B1", "B", 0.20, gamma=1.30),
        Basin("B2", "B", 0.25, gamma=0.70),
    ]


def write_summary(path: Path) -> None:
    basins = synthetic_hierarchy()
    selected = maximal_eligible(basins)
    naive = naive_eligible_volume(basins)
    clean = antichain_volume(basins)
    text = [
        "Topological basin antichain audit",
        f"eligible_nodes={[b.name for b in basins if is_v_eligible(b)]}",
        f"maximal_eligible={[b.name for b in selected]}",
        f"antichain={is_antichain(selected, basins)}",
        f"naive_nested_volume_sum={naive:.6f}",
        f"maximal_antichain_volume={clean:.6f}",
        "note=synthetic hierarchy only; no cosmological abundance is inferred",
    ]
    path.write_text("\n".join(text) + "\n", encoding="utf-8")


def main() -> None:
    out = ROOT / "results" / "basin_selector_summary.txt"
    out.parent.mkdir(exist_ok=True)
    write_summary(out)
    print(out.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
