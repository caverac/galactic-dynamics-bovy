"""MkDocs macros for dynamic content generation."""

from galactic_dynamics_bovy import __version__


def define_env(env):
    """Define macros available in documentation.

    Parameters
    ----------
    env : MacrosPlugin
        The MkDocs macros plugin environment.
    """

    @env.macro
    def project_version():
        """Return the current project version."""
        return __version__
