#!/usr/bin/env python
"""Generate documentation assets from solution plots.

This script generates all figures for the documentation by calling
plot functions from solution modules.

Usage:
    uv run python scripts/generate_assets.py
"""

from pathlib import Path

from galactic_dynamics_bovy.chapter02.ic_2574 import plot_rotation_curve
from galactic_dynamics_bovy.chapter02.nbody_1d import plot_force_comparison
from galactic_dynamics_bovy.chapter02.rvir_mvir import plot_rvir_mvir_delta
from galactic_dynamics_bovy.chapter02.spherical_exponential import plot_exponential_vcirc
from galactic_dynamics_bovy.chapter02.vcirc_profiles import plot_vcirc_profiles

ASSETS_DIR = Path(__file__).parent.parent / "docs" / "assets" / "generated"


def main() -> None:
    """Generate all documentation assets."""
    ASSETS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Assets directory: {ASSETS_DIR}")

    # Chapter 2
    print("Generating IC 2574 rotation curve...")
    plot_rotation_curve(ASSETS_DIR / "ic2574_rotation_curve.png")

    print("Generating NFW virial mass/radius plot...")
    plot_rvir_mvir_delta(ASSETS_DIR / "nfw_rvir_mvir_delta.png")

    print("Generating circular velocity profiles...")
    plot_vcirc_profiles(ASSETS_DIR / "vcirc_profiles.png")

    print("Generating exponential disk velocity curves...")
    plot_exponential_vcirc(ASSETS_DIR / "exponential_vcirc.png")

    print("Generating 1D N-body force comparison...")
    plot_force_comparison(path=ASSETS_DIR / "nbody_1d_force.png")

    print("Done.")


if __name__ == "__main__":
    main()
