from __future__ import annotations

from collections.abc import Sequence

from jpype.types import JInt


MIN_MESH_SIZE_DIVISOR = 33.3
PML_SWEEP_ELEMENTS = 8
MANUAL_MESH_SIZE_PARAMETERS: tuple[tuple[float, float, float], ...] = (
    (1.3, 0.2, 1.0),
    (1.35, 0.3, 0.85),
    (1.4, 0.4, 0.7),
    (1.45, 0.5, 0.6),
    (1.5, 0.6, 0.5),
    (1.6, 0.7, 0.4),
    (1.7, 0.8, 0.3),
    (1.85, 0.9, 0.2),
    (2.0, 1.0, 0.1),
)


def resolve_manual_mesh_size_parameters(mesh_auto_size: int) -> tuple[float, float, float]:
    if type(mesh_auto_size) is not int or not 1 <= mesh_auto_size <= len(
        MANUAL_MESH_SIZE_PARAMETERS
    ):
        raise ValueError("mesh_auto_size must be an integer from 1 to 9")
    return MANUAL_MESH_SIZE_PARAMETERS[mesh_auto_size - 1]


def _mesh_size_expressions(refractive_index: str | None, min_size_divisor: float) -> tuple[str, str]:
    nr = "1" if refractive_index is None else str(refractive_index).strip()
    if not nr or nr == "1":
        hmax = "lda0/5"
    else:
        hmax = f"lda0/({nr})/5"
    return hmax, f"({hmax})/{min_size_divisor:g}"


def _set_custom_size(
    size_feature,
    hmax: str,
    hmin: str,
    mesh_size_parameters: tuple[float, float, float],
    *,
    set_active: bool,
) -> None:
    hgrad, hcurve, hnarrow = mesh_size_parameters
    size_feature.set("custom", "on")
    size_feature.set("hmax", hmax)
    size_feature.set("hmin", hmin)
    size_feature.set("hgrad", hgrad)
    size_feature.set("hcurve", hcurve)
    size_feature.set("hnarrow", hnarrow)
    if set_active:
        size_feature.set("hmaxactive", True)
        size_feature.set("hminactive", True)
        size_feature.set("hgradactive", True)
        size_feature.set("hcurveactive", True)
        size_feature.set("hnarrowactive", True)


def construct_manual_mesh(
    jmodel,
    layers: Sequence[object],
    *,
    mesh_auto_size: int = 5,
    mesh_tag: str = "mesh1",
    min_size_divisor: float = MIN_MESH_SIZE_DIVISOR,
    pml_sweep_elements: int = PML_SWEEP_ELEMENTS,
) -> None:
    """Configure a named-selection based finite-cavity mesh sequence."""
    mesh_size_parameters = resolve_manual_mesh_size_parameters(mesh_auto_size)
    mesh = jmodel.component("comp1").mesh(mesh_tag)

    air_hmax, air_hmin = _mesh_size_expressions("1", min_size_divisor)
    _set_custom_size(
        mesh.feature("size"),
        air_hmax,
        air_hmin,
        mesh_size_parameters,
        set_active=False,
    )

    for layer_idx, layer in enumerate(layers):
        size_tag = f"size_layer{layer_idx}"
        mesh.create(size_tag, "Size")
        size_feature = mesh.feature(size_tag)
        size_feature.selection().named(f"sel_layer_{layer_idx}_mat_dom")

        hmax, hmin = _mesh_size_expressions(getattr(layer, "refractive_index", None), min_size_divisor)
        _set_custom_size(
            size_feature,
            hmax,
            hmin,
            mesh_size_parameters,
            set_active=True,
        )

    mesh.create("ftet1", "FreeTet")
    mesh.feature("ftet1").selection().named("sel_physical_dom")

    mesh.create("swe1", "Sweep")
    sweep = mesh.feature("swe1")
    sweep.selection().named("sel_pml_dom")
    sweep.create("dis1", "Distribution")
    sweep.feature("dis1").set("numelem", JInt(int(pml_sweep_elements)))
