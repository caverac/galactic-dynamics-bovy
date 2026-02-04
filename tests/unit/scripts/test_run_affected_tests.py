"""Unit tests for run_affected_tests script."""

# pylint: disable=import-outside-toplevel

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


class TestGetTestFileForSource:
    """Tests for get_test_file_for_source function."""

    def test_test_file_returns_itself(self, tmp_path: Path) -> None:
        """A test file should return itself."""
        from scripts.run_affected_tests import get_test_file_for_source

        test_file = tmp_path / "test_something.py"
        result = get_test_file_for_source(test_file)
        assert result == test_file

    def test_non_python_file_returns_none(self) -> None:
        """Non-Python files should return None."""
        from scripts.run_affected_tests import get_test_file_for_source

        result = get_test_file_for_source(Path("README.md"))
        assert result is None

    def test_file_outside_package_returns_none(self) -> None:
        """Files outside the package should return None."""
        from scripts.run_affected_tests import get_test_file_for_source

        result = get_test_file_for_source(Path("some_other_package/module.py"))
        assert result is None

    def test_source_file_maps_to_test_file(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Source files should map to their test files."""
        from scripts.run_affected_tests import get_test_file_for_source

        # Create the test file structure
        test_dir = tmp_path / "tests" / "unit" / "galactic_dynamics_bovy" / "chapter02"
        test_dir.mkdir(parents=True)
        test_file = test_dir / "test_ic_2574.py"
        test_file.touch()

        # Change to tmp_path so relative paths work
        monkeypatch.chdir(tmp_path)

        source_file = Path("galactic_dynamics_bovy/chapter02/ic_2574.py")
        result = get_test_file_for_source(source_file)

        assert result == Path("tests/unit/galactic_dynamics_bovy/chapter02/test_ic_2574.py")

    def test_source_file_without_test_returns_none(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Source files without corresponding test files should return None."""
        from scripts.run_affected_tests import get_test_file_for_source

        monkeypatch.chdir(tmp_path)

        source_file = Path("galactic_dynamics_bovy/chapter02/ic_2574.py")
        result = get_test_file_for_source(source_file)

        assert result is None


class TestMain:
    """Tests for main function."""

    @patch("scripts.run_affected_tests.sys")
    def test_no_files_returns_zero(self, mock_sys: MagicMock) -> None:
        """Should return 0 when no files are provided."""
        from scripts.run_affected_tests import main

        mock_sys.argv = ["run_affected_tests.py"]
        result = main()
        assert result == 0

    @patch("scripts.run_affected_tests.subprocess")
    @patch("scripts.run_affected_tests.get_test_file_for_source")
    @patch("scripts.run_affected_tests.sys")
    def test_no_test_files_found_returns_zero(
        self,
        mock_sys: MagicMock,
        mock_get_test: MagicMock,
        mock_subprocess: MagicMock,
    ) -> None:
        """Should return 0 when no test files are found."""
        from scripts.run_affected_tests import main

        mock_sys.argv = ["run_affected_tests.py", "some_file.txt"]
        mock_get_test.return_value = None

        result = main()

        assert result == 0
        mock_subprocess.run.assert_not_called()

    @patch("scripts.run_affected_tests.subprocess")
    @patch("scripts.run_affected_tests.get_test_file_for_source")
    @patch("scripts.run_affected_tests.sys")
    def test_runs_pytest_for_found_tests(
        self,
        mock_sys: MagicMock,
        mock_get_test: MagicMock,
        mock_subprocess: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Should run pytest for found test files."""
        from scripts.run_affected_tests import main

        test_file = tmp_path / "test_something.py"
        test_file.touch()

        mock_sys.argv = ["run_affected_tests.py", "source.py"]
        mock_sys.executable = "/usr/bin/python"
        mock_get_test.return_value = test_file
        mock_subprocess.run.return_value = MagicMock(returncode=0)

        main()

        mock_subprocess.run.assert_called_once()
        call_args = mock_subprocess.run.call_args
        assert "-m" in call_args[0][0]
        assert "pytest" in call_args[0][0]
        assert "--no-cov" in call_args[0][0]

    @patch("scripts.run_affected_tests.subprocess")
    @patch("scripts.run_affected_tests.get_test_file_for_source")
    @patch("scripts.run_affected_tests.sys")
    def test_returns_pytest_exit_code(
        self,
        mock_sys: MagicMock,
        mock_get_test: MagicMock,
        mock_subprocess: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Should return pytest's exit code."""
        from scripts.run_affected_tests import main

        test_file = tmp_path / "test_something.py"
        test_file.touch()

        mock_sys.argv = ["run_affected_tests.py", "source.py"]
        mock_sys.executable = "/usr/bin/python"
        mock_get_test.return_value = test_file
        mock_subprocess.run.return_value = MagicMock(returncode=1)

        result = main()

        assert result == 1

    @patch("scripts.run_affected_tests.subprocess")
    @patch("scripts.run_affected_tests.get_test_file_for_source")
    @patch("scripts.run_affected_tests.sys")
    def test_deduplicates_test_files(
        self,
        mock_sys: MagicMock,
        mock_get_test: MagicMock,
        mock_subprocess: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Should not run the same test file twice."""
        from scripts.run_affected_tests import main

        test_file = tmp_path / "test_something.py"
        test_file.touch()

        mock_sys.argv = ["run_affected_tests.py", "source1.py", "source2.py"]
        mock_sys.executable = "/usr/bin/python"
        # Both source files map to the same test file
        mock_get_test.return_value = test_file
        mock_subprocess.run.return_value = MagicMock(returncode=0)

        main()

        # pytest should only be called once with one test file
        call_args = mock_subprocess.run.call_args[0][0]
        test_file_args = [arg for arg in call_args if "test_" in arg]
        assert len(test_file_args) == 1
