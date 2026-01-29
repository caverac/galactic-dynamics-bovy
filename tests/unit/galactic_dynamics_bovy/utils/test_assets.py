"""Unit tests for asset generation utilities."""

# pylint: disable=import-outside-toplevel,protected-access

from pathlib import Path

from matplotlib import pyplot as plt

from galactic_dynamics_bovy.utils import assets
from galactic_dynamics_bovy.utils.assets import (
    figures_are_equal,
    get_registered_assets,
    register_asset,
    save_figure_if_changed,
)


class TestRegisterAsset:
    """Tests for register_asset decorator."""

    def test_registers_function(self) -> None:
        """Decorator should add function to registry."""
        original_registry = assets._registry.copy()

        try:

            @register_asset("test_output.png")
            def test_plot(*, path: Path | None = None) -> None:  # pylint: disable=unused-argument
                pass

            registry = get_registered_assets()
            assert "test_output.png" in registry
            assert registry["test_output.png"][0] is test_plot
        finally:
            # Restore original registry
            assets._registry = original_registry

    def test_preserves_function(self) -> None:
        """Decorator should return the original function unchanged."""
        original_registry = assets._registry.copy()

        try:

            @register_asset("another_test.png")
            def my_plot(*, path: Path | None = None) -> None:  # pylint: disable=unused-argument
                """My docstring."""
                return None

            assert my_plot.__name__ == "my_plot"
            assert my_plot.__doc__ == "My docstring."
        finally:
            assets._registry = original_registry


class TestFiguresAreEqual:
    """Tests for figures_are_equal function."""

    def test_returns_false_for_nonexistent_file(self, tmp_path: Path) -> None:
        """Should return False if file doesn't exist."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3])

        result = figures_are_equal(tmp_path / "nonexistent.png", fig)

        plt.close(fig)
        assert result is False

    def test_returns_true_for_identical_figure(self, tmp_path: Path) -> None:
        """Should return True for identical figures."""
        # Create and save a figure
        fig1, ax1 = plt.subplots()
        ax1.plot([1, 2, 3], [1, 4, 9])
        path = tmp_path / "test.png"
        fig1.savefig(path, dpi=100)
        plt.close(fig1)

        # Create identical figure
        fig2, ax2 = plt.subplots()
        ax2.plot([1, 2, 3], [1, 4, 9])

        result = figures_are_equal(path, fig2, dpi=100)

        plt.close(fig2)
        assert result is True

    def test_returns_false_for_different_figure(self, tmp_path: Path) -> None:
        """Should return False for different figures."""
        # Create and save a figure
        fig1, ax1 = plt.subplots()
        ax1.plot([1, 2, 3], [1, 4, 9])
        path = tmp_path / "test.png"
        fig1.savefig(path, dpi=100)
        plt.close(fig1)

        # Create different figure
        fig2, ax2 = plt.subplots()
        ax2.plot([1, 2, 3], [9, 4, 1])  # Different data

        result = figures_are_equal(path, fig2, dpi=100)

        plt.close(fig2)
        assert result is False

    def test_returns_false_for_different_size(self, tmp_path: Path) -> None:
        """Should return False for figures with different dimensions."""
        # Create and save a small figure
        fig1, ax1 = plt.subplots(figsize=(4, 3))
        ax1.plot([1, 2, 3])
        path = tmp_path / "test.png"
        fig1.savefig(path, dpi=100)
        plt.close(fig1)

        # Create larger figure
        fig2, ax2 = plt.subplots(figsize=(8, 6))
        ax2.plot([1, 2, 3])

        result = figures_are_equal(path, fig2, dpi=100)

        plt.close(fig2)
        assert result is False


class TestSaveFigureIfChanged:
    """Tests for save_figure_if_changed function."""

    def test_saves_new_file(self, tmp_path: Path) -> None:
        """Should save file if it doesn't exist."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3])
        path = tmp_path / "new_plot.png"

        result = save_figure_if_changed(fig, path, dpi=100)

        plt.close(fig)
        assert result is True
        assert path.exists()

    def test_skips_unchanged_file(self, tmp_path: Path) -> None:
        """Should not save if figure is unchanged."""
        # Create and save initial figure
        fig1, ax1 = plt.subplots()
        ax1.plot([1, 2, 3], [1, 4, 9])
        path = tmp_path / "test.png"
        fig1.savefig(path, dpi=100)
        plt.close(fig1)

        # Get initial modification time
        initial_mtime = path.stat().st_mtime

        # Create identical figure and try to save
        fig2, ax2 = plt.subplots()
        ax2.plot([1, 2, 3], [1, 4, 9])

        result = save_figure_if_changed(fig2, path, dpi=100)

        plt.close(fig2)
        assert result is False
        # File should not have been modified
        assert path.stat().st_mtime == initial_mtime

    def test_saves_changed_file(self, tmp_path: Path) -> None:
        """Should save if figure has changed."""
        # Create and save initial figure
        fig1, ax1 = plt.subplots()
        ax1.plot([1, 2, 3], [1, 4, 9])
        path = tmp_path / "test.png"
        fig1.savefig(path, dpi=100)
        plt.close(fig1)

        # Create different figure and save
        fig2, ax2 = plt.subplots()
        ax2.plot([1, 2, 3], [9, 4, 1])  # Different data

        result = save_figure_if_changed(fig2, path, dpi=100)

        plt.close(fig2)
        assert result is True

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        """Should create parent directories if needed."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3])
        path = tmp_path / "subdir" / "nested" / "plot.png"

        result = save_figure_if_changed(fig, path, dpi=100)

        plt.close(fig)
        assert result is True
        assert path.exists()
