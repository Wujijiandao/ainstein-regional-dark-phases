#!/usr/bin/env python3
from __future__ import annotations
import subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CORE = ["analysis_phase.py", "analysis_regional_state.py", "analysis_excursion_phase.py", "analysis_power_clock.py", "analysis_curvature_selector.py", "analysis_basin_selector.py", "analysis_conversion_sign.py", "analysis_basin_dynamics.py", "analysis_interface_coarsening.py", "analysis_phase_interface.py", "analysis_phase_front.py", "analysis_bias_matching.py", "analysis_synthetic_3d_basins.py", "analysis_zeldovich_basin_history.py"]
for script in CORE:
    print(f"\n=== {script} ===")
    subprocess.run([sys.executable, str(ROOT / "src" / script)], check=True, cwd=ROOT)

rar = ROOT / "data" / "RAR.mrt"
if rar.exists():
    print("\n=== analysis_emg.py ===")
    subprocess.run([sys.executable, str(ROOT / "src" / "analysis_emg.py")], check=True, cwd=ROOT)
else:
    print("\n=== analysis_emg.py: skipped ===")
    print("SPARC data/RAR.mrt is not bundled in the GitHub-safe repository.")
    print("Run: python scripts/fetch_sparc.py")

print("\nAll available analyses completed.")
