#!/usr/bin/env python
"""Generate documentation assets from solution plots.

This script generates all figures for the documentation by calling
plot functions from solution modules.

Usage:
    uv run python scripts/generate_assets.py
"""

from pathlib import Path

ASSETS_DIR = Path(__file__).parent.parent / "docs" / "assets" / "generated"


def main() -> None:
    """Generate all documentation assets."""
    from galactic_dynamics_bovy.chapter02.ic_2574 import plot_rotation_curve

    ASSETS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Assets directory: {ASSETS_DIR}")

    # Chapter 2
    print("Generating IC 2574 rotation curve...")
    plot_rotation_curve(ASSETS_DIR / "ic2574_rotation_curve.png")
    print("Done.")


if __name__ == "__main__":
    main()
