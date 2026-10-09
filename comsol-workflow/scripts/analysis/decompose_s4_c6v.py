"""Historical entrypoint; implementation: scripts.analysis.decompose_boundary_c6v."""
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from comsol_workflow._compat import redirect

redirect(__name__, 'scripts.analysis.decompose_boundary_c6v')
