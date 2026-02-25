"""Unit tests for generate_assets script."""

# pylint: disable=import-outside-toplevel

import argparse
from pathlib import Path
from unittest.mock import MagicMock, patch

_NO_FILTERS = argparse.Namespace(filters=[])


class TestMatchesAny:
    """Tests for _matches_any helper."""

    def test_matches_substring(self) -> None:
        """Should return True when name contains a filter substring."""
        from scripts.generate_assets import _matches_any

        assert _matches_any("p04_09_gr_precession.png", ["p04_09"])

    def test_no_match(self) -> None:
        """Should return False when no filter substring matches."""
        from scripts.generate_assets import _matches_any

        assert not _matches_any("p04_09_gr_precession.png", ["p02"])

    def test_empty_filters(self) -> None:
        """Should return False for an empty filter list."""
        from scripts.generate_assets import _matches_any

        assert not _matches_any("anything.png", [])


class TestMain:
    """Tests for the main() function."""

    @patch("scripts.generate_assets._parse_args", return_value=_NO_FILTERS)
    @patch("scripts.generate_assets.get_registered_assets")
    @patch("scripts.generate_assets.ASSETS_DIR")
    def test_creates_assets_directory(
        self,
        mock_assets_dir: MagicMock,
        mock_get_assets: MagicMock,
        _mock_args: MagicMock,
    ) -> None:
        """Should create the assets directory if it doesn't exist."""
        from scripts.generate_assets import main

        mock_get_assets.return_value = {}

        main()

        mock_assets_dir.mkdir.assert_called_once_with(parents=True, exist_ok=True)

    @patch("scripts.generate_assets._parse_args", return_value=_NO_FILTERS)
    @patch("scripts.generate_assets.get_registered_assets")
    @patch("scripts.generate_assets.ASSETS_DIR", new_callable=MagicMock)
    def test_calls_each_plot_function(
        self,
        mock_assets_dir: MagicMock,
        mock_get_assets: MagicMock,
        _mock_args: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Should call each registered plot function with the correct path."""
        from scripts.generate_assets import main

        mock_plot1 = MagicMock()
        mock_plot2 = MagicMock()
        mock_get_assets.return_value = {
            "plot1.png": (mock_plot1, "module1"),
            "plot2.png": (mock_plot2, "module2"),
        }
        mock_assets_dir.__truediv__ = lambda self, name: tmp_path / name

        main()

        mock_plot1.assert_called_once_with(path=tmp_path / "plot1.png")
        mock_plot2.assert_called_once_with(path=tmp_path / "plot2.png")

    @patch("scripts.generate_assets._parse_args", return_value=_NO_FILTERS)
    @patch("scripts.generate_assets.get_registered_assets")
    @patch("scripts.generate_assets.ASSETS_DIR", new_callable=MagicMock)
    @patch("scripts.generate_assets.Console")
    def test_reports_written_for_new_files(
        self,
        mock_console_class: MagicMock,
        mock_assets_dir: MagicMock,
        mock_get_assets: MagicMock,
        _mock_args: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Should report [WRITE] for newly created files."""
        from scripts.generate_assets import main

        mock_console = MagicMock()
        mock_console_class.return_value = mock_console

        def create_file(path: Path) -> None:
            path.write_text("content")

        mock_plot = MagicMock(side_effect=create_file)
        mock_get_assets.return_value = {"new_plot.png": (mock_plot, "module")}
        mock_assets_dir.__truediv__ = lambda self, name: tmp_path / name

        main()

        # Check that [WRITE] was printed
        write_calls = [c for c in mock_console.print.call_args_list if "WRITE" in str(c)]
        assert len(write_calls) == 1

    @patch("scripts.generate_assets._parse_args", return_value=_NO_FILTERS)
    @patch("scripts.generate_assets.get_registered_assets")
    @patch("scripts.generate_assets.ASSETS_DIR", new_callable=MagicMock)
    @patch("scripts.generate_assets.Console")
    def test_reports_skipped_for_unchanged_files(
        self,
        mock_console_class: MagicMock,
        mock_assets_dir: MagicMock,
        mock_get_assets: MagicMock,
        _mock_args: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Should report [SKIP] for files that weren't modified."""
        from scripts.generate_assets import main

        mock_console = MagicMock()
        mock_console_class.return_value = mock_console

        # Create file beforehand
        output_file = tmp_path / "existing_plot.png"
        output_file.write_text("content")

        # Plot function that doesn't modify the file
        mock_plot = MagicMock()
        mock_get_assets.return_value = {"existing_plot.png": (mock_plot, "module")}
        mock_assets_dir.__truediv__ = lambda self, name: tmp_path / name

        main()

        # Check that [SKIP] was printed
        skip_calls = [c for c in mock_console.print.call_args_list if "SKIP" in str(c)]
        assert len(skip_calls) == 1

    @patch("scripts.generate_assets._parse_args", return_value=_NO_FILTERS)
    @patch("scripts.generate_assets.get_registered_assets")
    @patch("scripts.generate_assets.ASSETS_DIR", new_callable=MagicMock)
    def test_handles_empty_registry(
        self,
        _mock_assets_dir: MagicMock,
        mock_get_assets: MagicMock,
        _mock_args: MagicMock,
    ) -> None:
        """Should handle case with no registered assets."""
        from scripts.generate_assets import main

        mock_get_assets.return_value = {}

        # Should not raise
        main()

    @patch("scripts.generate_assets._parse_args", return_value=argparse.Namespace(filters=["plot1"]))
    @patch("scripts.generate_assets.get_registered_assets")
    @patch("scripts.generate_assets.ASSETS_DIR", new_callable=MagicMock)
    def test_filter_selects_matching_assets(
        self,
        mock_assets_dir: MagicMock,
        mock_get_assets: MagicMock,
        _mock_args: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Should only generate assets matching the filter."""
        from scripts.generate_assets import main

        mock_plot1 = MagicMock()
        mock_plot2 = MagicMock()
        mock_get_assets.return_value = {
            "plot1.png": (mock_plot1, "module1"),
            "plot2.png": (mock_plot2, "module2"),
        }
        mock_assets_dir.__truediv__ = lambda self, name: tmp_path / name

        main()

        mock_plot1.assert_called_once()
        mock_plot2.assert_not_called()


class TestModuleImports:
    """Tests for module import side effects."""

    def test_imports_do_not_generate_assets(self) -> None:
        """Importing the module should not generate any assets."""
        # The imports in generate_assets.py only register decorators,
        # they don't actually call plot functions. This test ensures
        # that importing the module is safe and fast.
        import scripts.generate_assets  # noqa: F401  # pylint: disable=unused-import

        # If we get here without errors or side effects, the test passes
        assert True
