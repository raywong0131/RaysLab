import base64
from pathlib import Path
import re

from comsol_workflow.field_plotting import save_png_plot_html
from scripts.analysis.export_finite_field_html import export_field_html


def test_html_plot_is_self_contained_and_preserves_png_bytes(tmp_path: Path):
    png_path = tmp_path / "field.png"
    png_bytes = b"\x89PNG\r\n\x1a\nfield-image"
    png_path.write_bytes(png_bytes)
    html_path = tmp_path / "field.html"

    save_png_plot_html(
        html_path,
        png_path,
        document_title="Re(Hz) & field",
    )

    html = html_path.read_text(encoding="utf-8")
    assert "<title>Re(Hz) &amp; field</title>" in html
    assert "background: transparent" in html
    assert "max-width: 100%" in html
    assert "width: auto" in html
    match = re.search(r"data:image/png;base64,([A-Za-z0-9+/=]+)", html)
    assert match is not None
    assert base64.b64decode(match.group(1)) == png_bytes
    assert "https://" not in html
    assert "http://" not in html


def test_export_field_html_writes_all_three_field_documents(tmp_path: Path):
    export_dir = tmp_path / "mode1" / "11_simulation_exports"
    export_dir.mkdir(parents=True)
    for name in ("Hz_Re_2d.png", "Hz_Im_2d.png", "Wem_2d.png"):
        (export_dir / name).write_bytes(b"png")

    outputs = export_field_html(export_dir.parent)

    assert [path.name for path in outputs] == [
        "Hz_Re_2d.html",
        "Hz_Im_2d.html",
        "Wem_2d.html",
    ]
    assert all(path.is_file() for path in outputs)
