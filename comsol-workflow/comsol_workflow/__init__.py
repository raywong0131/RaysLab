"""Reusable geometry, simulation, and post-processing tools for COMSOL workflows.

The package root intentionally performs no eager imports so that lightweight
geometry and analysis modules can be used without importing ``mph`` or starting
COMSOL.  Import concrete capabilities from their defining modules.
"""

__all__: list[str] = []
