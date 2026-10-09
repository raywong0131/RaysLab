"""Check that staging stays separate without losing historical regressions."""
from pathlib import Path
import tomllib

from scripts.check_offline import COMSOL_EXAMPLES, discover_tests


def test_staging_discovery_preserves_archive_and_excludes_computation():
    root = Path(__file__).resolve().parents[3]
    discovered = set(discover_tests())
    assert Path(__file__).resolve() in discovered
    assert root / "tests/archive/test_farfield_fft.py" in discovered
    assert all(p.name not in COMSOL_EXAMPLES for p in discovered)
    assert not any(".work" in p.relative_to(root / "tests").parts for p in discovered)
    settings = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["pytest"]["ini_options"]
    assert "archive" in settings["norecursedirs"] and settings["python_files"] == ["test_*.py"]
    assert (root / "tmp/archive").is_dir() and (root / "tests/archive").is_dir()
