from __future__ import annotations

from dataclasses import dataclass

from .finite_geometry import FiniteGeometryPlan
from .simulation_utils import BoundarySpec, LayerSpec


@dataclass(frozen=True)
class FiniteComsolGeometryInput:
    boundary: BoundarySpec
    layers: tuple[LayerSpec, ...]


def compile_finite_geometry_input(
    plan: FiniteGeometryPlan,
    *,
    layer_height: str,
    refractive_index: str,
    label: str = "patterned_slab",
    quarter: bool = False,
    holes: list | tuple | None = None,
) -> FiniteComsolGeometryInput:
    """Compile a pure finite plan to the existing named-selection COMSOL backend."""
    footprint = plan.footprint
    if quarter:
        shape = (
            "quarter_hexagon_boundary"
            if footprint.shape == "hex"
            else "quarter_rectangle"
        )
    else:
        shape = footprint.boundary_shape
    boundary = BoundarySpec(
        shape,
        a=footprint.width,
        b=footprint.height if "rectangle" in shape else None,
    )
    resolved_holes = plan.holes if holes is None else holes
    layer = LayerSpec(
        height=layer_height,
        holes=resolved_holes,
        refractive_index=refractive_index,
        label=label,
    )
    return FiniteComsolGeometryInput(boundary=boundary, layers=(layer,))
