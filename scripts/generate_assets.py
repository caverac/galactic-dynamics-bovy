#!/usr/bin/env python
# pylint: disable=unused-import
# flake8: noqa=F401
"""Generate documentation assets from solution plots.

This script discovers all registered plot functions and generates their
corresponding figure assets for the documentation.

Usage:
    uv run python scripts/generate_assets.py              # generate all
    uv run python scripts/generate_assets.py p04_09 p02   # only matching names
"""

import argparse
from pathlib import Path
import re

from rich.console import Console

# Import chapter modules to trigger decorator registration
import galactic_dynamics_bovy.chapter02.effective_radius_gamma
import galactic_dynamics_bovy.chapter02.fermi_gas_vcirc
import galactic_dynamics_bovy.chapter02.ic_2574
import galactic_dynamics_bovy.chapter02.nbody_1d
import galactic_dynamics_bovy.chapter02.rvir_mvir
import galactic_dynamics_bovy.chapter02.spherical_exponential
import galactic_dynamics_bovy.chapter02.vcirc_profiles
import galactic_dynamics_bovy.chapter03.planet_comet_scattering
import galactic_dynamics_bovy.chapter04.adiabatic_isochrone
import galactic_dynamics_bovy.chapter04.gr_precession
import galactic_dynamics_bovy.chapter04.gr_precession_s2
from galactic_dynamics_bovy.utils.assets import get_registered_assets

ASSETS_DIR = Path(__file__).parent.parent / "docs" / "assets" / "generated"

_CHAPTER_PATTERN = re.compile(r"\.chapter(\d+)\.")


def _chapter_sort_key(item: tuple[str, tuple]) -> tuple[int, str]:
    """Sort key to order assets by chapter number, then by name."""
    output_name, (_, module_name) = item
    match = _CHAPTER_PATTERN.search(module_name)
    chapter_num = int(match.group(1)) if match else 999
    return (chapter_num, output_name)


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate documentation plot assets.")
    parser.add_argument(
        "filters",
        nargs="*",
        metavar="PATTERN",
        help="Only generate assets whose filename contains one of these substrings. "
        "If omitted, all assets are generated.",
    )
    return parser.parse_args()


def _matches_any(name: str, filters: list[str]) -> bool:
    """Return True if *name* contains any of the filter substrings."""
    return any(f in name for f in filters)


def main() -> None:
    """Generate documentation assets, optionally filtered by name."""
    args = _parse_args()
    console = Console()

    ASSETS_DIR.mkdir(parents=True, exist_ok=True)
    console.print(f"Assets directory: {ASSETS_DIR}")

    assets = get_registered_assets()

    if args.filters:
        assets = {k: v for k, v in assets.items() if _matches_any(k, args.filters)}
        console.print(f"Filter: {', '.join(args.filters)} -> {len(assets)} matching assets\n")
    else:
        console.print(f"Found {len(assets)} registered assets\n")

    written = 0
    skipped = 0
    current_chapter = None

    for output_name, (plot_func, module_name) in sorted(assets.items(), key=_chapter_sort_key):
        match = _CHAPTER_PATTERN.search(module_name)
        chapter = int(match.group(1)) if match else None
        if chapter != current_chapter:
            current_chapter = chapter
            label = f"Chapter {chapter:02d}" if chapter else "Other"
            console.print(f"[bold]{label}:[/bold]")

        output_path = ASSETS_DIR / output_name

        # Check if file exists and get mtime before generation
        existed = output_path.exists()
        mtime_before = output_path.stat().st_mtime if existed else None

        plot_func(path=output_path)

        # Check if file was actually written
        if not existed or output_path.stat().st_mtime != mtime_before:
            console.print(f"  [green]\\[WRITE][/green] {output_name}")
            written += 1
        else:
            console.print(f"  [dim]\\[SKIP]  {output_name} (unchanged)[/dim]")
            skipped += 1

    console.print(f"\nDone: [green]{written} written[/green], [dim]{skipped} unchanged[/dim]")


if __name__ == "__main__":  # pragma: no cover
    main()
