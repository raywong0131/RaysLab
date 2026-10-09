"""Output paths shared by retained regression and preview programs."""
import os
from pathlib import Path
import uuid


def output_root():
    value = os.environ.get("COMSOL_WORKFLOW_TEST_OUTPUT_ROOT")
    if not value:
        value = str(Path(__file__).resolve().parent / ".work" / ("manual_" + uuid.uuid4().hex[:8]))
        os.environ["COMSOL_WORKFLOW_TEST_OUTPUT_ROOT"] = value
    return Path(value)
