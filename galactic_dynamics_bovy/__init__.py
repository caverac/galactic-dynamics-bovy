"""Solutions to Dynamics and Astrophysics of Galaxies by Jo Bovy.

This package contains Python implementations and solutions to problems
from the textbook "Dynamics and Astrophysics of Galaxies" by Jo Bovy
(Princeton University Press, Princeton Series in Astrophysics).

The solutions leverage galpy, the galactic dynamics library developed
by the author of the textbook.

Reference:
    Bovy, J. (2026). Dynamics and Astrophysics of Galaxies.
    Princeton University Press. https://galaxiesbook.org/
"""

import matplotlib

__version__ = "0.1.0"

# Configure matplotlib for publication-quality figures
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
matplotlib.rcParams["mathtext.fontset"] = "stix"

__all__ = ["__version__"]
