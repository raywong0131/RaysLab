#!/usr/bin/env python3
"""Create self-contained HTML wrappers for finite-cavity field PNGs."""

from __future__ import annotations

import argparse
from pathlib import Path

from comsol_workflow.field_plotting import save_png_plot_html


FIELD_FILES = {
    "Hz_Re_2d.png": "Re(Hz)",
    "Hz_Im_2d.png": "Im(Hz)",
    "Wem_2d.png": "Wem",
}


def resolve_export_dir(path: Path) -> Path:
    path = Path(path).resolve()
    if path.name == "11_simulation_exports":
        export_dir = path
    else:
        export_dir = path / "11_simulation_exports"
    if not export_dir.is_dir():
        raise FileNotFoundError(f"11_simulation_exports not found: {export_dir}")
    return export_dir


def export_field_html(path: Path) -> list[Path]:
    export_dir = resolve_export_dir(path)
    outputs = []
    for png_name, title in FIELD_FILES.items():
        png_path = export_dir / png_name
        html_path = png_path.with_suffix(".html")
        save_png_plot_html(
            html_path,
            png_path,
            document_title=title,
        )
        outputs.append(html_path)
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "path",
        type=Path,
        help="mode directory or its 11_simulation_exports directory",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for output in export_field_html(args.path):
        print(output)


if __name__ == "__main__":
    main()
