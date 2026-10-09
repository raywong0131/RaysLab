"""Historical entrypoint; implementation: scripts.run_sweep.run_boundary_surface_signed_scan."""
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from comsol_workflow._compat import redirect

redirect(__name__, 'scripts.run_sweep.run_boundary_surface_signed_scan')
