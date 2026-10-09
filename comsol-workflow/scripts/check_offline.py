"""Run offline checks in separate solver-import and scientific processes."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
COMSOL_EXAMPLES = {
    "test_simulation_spatial.py", "test_simulation_utils_spatial_comsol.py",
    "test_simulation_plotting.py",
}
SCIENTIFIC_PREFIXES = ("test_s4_", "test_s5_", "test_dipolar_", "test_boundary_")


def discover_tests():
    """Read archived regressions explicitly and active checks without work trees."""
    base = ROOT / "tests"
    archived = list((base / "archive").glob("test_*.py"))
    active = [p for p in base.rglob("test_*.py")
              if not {"archive", ".work", "__pycache__"}.intersection(p.relative_to(base).parts)]
    return sorted(p for p in archived + active if p.name not in COMSOL_EXAMPLES)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=("core", "scientific"), help=argparse.SUPPRESS)
    parser.add_argument("--temp-root", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        os.environ["MPLCONFIGDIR"] = str(args.temp_root / "matplotlib")
        import pytest

        tests = [str(p) for p in discover_tests()
                 if p.name.startswith(SCIENTIFIC_PREFIXES) == (args.worker == "scientific")]
        if not tests:
            raise RuntimeError(f"No offline checks found for {args.worker}; refusing an empty pass")
        return pytest.main(["-q", "-p", "no:cacheprovider",
                            "--import-mode=importlib",
                            "--basetemp", str(args.temp_root / args.worker), *tests])

    task_parent = ROOT / "tests/.work"
    task_parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix="cw-", dir=task_parent)).resolve()
    try:
        codes = [subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()),
                                "--worker", group, "--temp-root", str(temporary)],
                               cwd=ROOT).returncode for group in ("core", "scientific")]
        return max(codes)
    finally:
        if temporary.parent != task_parent.resolve() or not temporary.name.startswith("cw-"):
            raise ValueError("Unexpected offline-check cleanup path")
        shutil.rmtree(temporary)


if __name__ == "__main__":
    raise SystemExit(main())
