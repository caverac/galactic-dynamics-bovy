#!/usr/bin/env python
# pylint: disable=unused-import
# flake8: noqa=F401
"""Generate documentation assets from solution plots.

This script discovers all registered plot functions and generates their
corresponding figure assets for the documentation.

Usage:
    uv run python scripts/generate_assets.py
"""

from pathlib import Path

# Import chapter modules to trigger decorator registration
import galactic_dynamics_bovy.chapter02.effective_radius_gamma
import galactic_dynamics_bovy.chapter02.ic_2574
import galactic_dynamics_bovy.chapter02.nbody_1d
import galactic_dynamics_bovy.chapter02.rvir_mvir
import galactic_dynamics_bovy.chapter02.spherical_exponential
import galactic_dynamics_bovy.chapter02.vcirc_profiles
from galactic_dynamics_bovy.utils.assets import get_registered_assets

ASSETS_DIR = Path(__file__).parent.parent / "docs" / "assets" / "generated"


def main() -> None:
    """Generate all documentation assets."""
    ASSETS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Assets directory: {ASSETS_DIR}")

    assets = get_registered_assets()
    print(f"Found {len(assets)} registered assets\n")

    written = 0
    skipped = 0

    for output_name, (plot_func, _) in sorted(assets.items()):
        output_path = ASSETS_DIR / output_name

        # Check if file exists and get mtime before generation
        existed = output_path.exists()
        mtime_before = output_path.stat().st_mtime if existed else None

        plot_func(path=output_path)

        # Check if file was actually written
        if not existed or output_path.stat().st_mtime != mtime_before:
            print(f"  [WRITE] {output_name}")
            written += 1
        else:
            print(f"  [SKIP]  {output_name} (unchanged)")
            skipped += 1

    print(f"\nDone: {written} written, {skipped} unchanged")


if __name__ == "__main__":
    main()
