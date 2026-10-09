"""Check the real figure saver in an isolated directory; no scientific run."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from scripts.analysis import analyze_dipolar_singularity as output


def test_save_separates_png_and_pdf_and_tracks_both():
    with TemporaryDirectory(prefix="dipolar-output-") as folder:
        root = Path(folder)
        (root / "80_logs").mkdir()
        with patch.object(output, "WORK", root):
            fig, ax = output.plt.subplots(figsize=(2, 2), layout="constrained")
            ax.plot([0, 1], [0, 1])
            output.save(fig, "Fig.1_check", dpi=72)
        png = root / "10_overview/Fig.1_check.png"
        pdf = root / "11_pdf/Fig.1_check.pdf"
        assert png.read_bytes().startswith(b"\x89PNG")
        assert pdf.read_bytes().startswith(b"%PDF")
        assert list((root / "10_overview").iterdir()) == [png]
        record = json.loads((root / "80_logs/figure_manifest.json").read_text(encoding="utf8"))
        assert record["Fig.1_check"]["files"] == [
            "10_overview/Fig.1_check.png", "11_pdf/Fig.1_check.pdf"]
