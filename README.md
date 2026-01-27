# Galactic Dynamics Bovy

Solutions to problems from **Dynamics and Astrophysics of Galaxies** by Jo Bovy (Princeton University Press, Princeton Series in Astrophysics).

## Book Reference

> Bovy, J. (2026). *Dynamics and Astrophysics of Galaxies*. Princeton University Press.
> - Online version: [galaxiesbook.org](https://galaxiesbook.org/)
> - Publisher: [Princeton University Press](https://press.princeton.edu/books/paperback/9780691212876/dynamics-and-astrophysics-of-galaxies)

This textbook provides an in-depth introduction to the dynamics, formation, and evolution of galaxies. The solutions in this repository leverage [galpy](https://www.galpy.org/), the galactic dynamics library developed by the author.

## Installation

This project uses [uv](https://docs.astral.sh/uv/) for dependency management.

```bash
# Clone the repository
git clone https://github.com/yourusername/galactic-dynamics-bovy.git
cd galactic-dynamics-bovy

# Install dependencies
uv sync

# Or install with specific groups
uv sync --group dev --group docs
```

## Usage

```python
# Import the package
import galactic_dynamics_bovy

# Import solutions from specific chapters (once available)
# from galactic_dynamics_bovy.chapter01 import solution_function
```

## Development

```bash
# Run tests
uv run pytest

# Run linting
uv run black .
uv run flake8

# Type checking
uv run mypy galactic_dynamics_bovy

# Serve documentation locally
uv run mkdocs serve
```

## Project Structure

```
galactic-dynamics-bovy/
├── galactic_dynamics_bovy/     # Main source package
│   ├── __init__.py
│   └── chapter*/               # Chapter-specific solutions
├── tests/                      # Test suite
│   └── unit/
├── docs/                       # Documentation
│   └── assets/generated/       # Generated figures
├── scripts/                    # Utility scripts
├── pyproject.toml              # Project configuration
└── README.md
```

## License

MIT License
