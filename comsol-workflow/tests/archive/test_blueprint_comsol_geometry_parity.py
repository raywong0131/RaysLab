from __future__ import annotations

import numpy as np
import pytest

from blueprint_core.geometry import build_hex_layout
from comsol_workflow.basis_utils import (
    fourier_to_standard_basis,
    get_fourier_basis_min_nonzero_values,
)
from comsol_workflow.cladding_shift_profile import (
    CladdingShiftGeometry,
    CladdingShiftProfile,
)
from comsol_workflow.finite_geometry import (
    FiniteGeometrySpec,
    build_finite_geometry_plan,
)
from comsol_workflow.geometry_utils import create_hexagon_design


def _geometry_mapping(
    shift_kind: str,
    profile: dict[str, object],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "length_unit": "um",
        "lattice": {
            "kind": "triangular_hexagonal_v1",
            "period_um": 0.82,
            "origin_um": [0.0, 0.0],
            "orientation_deg": 0.0,
        },
        "unit_cell": {
            "kind": "six_triangular_holes_fourier_v1",
            "radial_reference_fraction": 1.0 / 3.0,
            "edge_reference_nm": 230.0,
            "cavity": {"b0_nm": 245.0, "eta": 0.96, "zeta": 1.156},
            "cladding": {"b0_nm": 242.0, "eta": 0.98, "zeta": 0.93},
        },
        "finite": {
            "kind": "hex_shells_v1",
            "cavity_layers": 2,
            "cladding_layers": 2,
            "empty_buffer_periods": 0.5,
        },
        "cladding_shift": {
            "geometry": {"kind": shift_kind},
            "profile": profile,
            "x_factor": 0.05,
            "y_factor": 0.08,
        },
    }


def _current_holes(
    mapping: dict[str, object],
    region: str,
) -> list[np.ndarray]:
    lattice = mapping["lattice"]
    unit_cell = mapping["unit_cell"]
    compact = unit_cell[region]
    period = float(lattice["period_um"])
    radial_reference = period * float(unit_cell["radial_reference_fraction"])
    edge_reference_um = float(unit_cell["edge_reference_nm"]) / 1000.0
    edge_dc_scale = (
        float(compact["b0_nm"]) / float(unit_cell["edge_reference_nm"])
    ) ** 2

    radial_fourier = np.zeros(6)
    theta_fourier = np.zeros(6)
    edge_square_fourier = np.zeros(6)
    phi_fourier = np.zeros(6)
    radial_fourier[0] = radial_reference * float(compact["eta"])
    edge_square_fourier[0] = edge_reference_um**2 * edge_dc_scale
    edge_square_fourier[3] = (
        edge_reference_um**2
        * edge_dc_scale
        * (float(compact["zeta"]) ** 2 - 1.0)
        / 2.0
    )
    minimum = get_fourier_basis_min_nonzero_values(6)
    rows = np.column_stack(
        [
            fourier_to_standard_basis(radial_fourier / minimum, 6),
            fourier_to_standard_basis(theta_fourier / minimum, 6),
            np.sqrt(
                np.maximum(
                    fourier_to_standard_basis(edge_square_fourier / minimum, 6),
                    0.0,
                )
            ),
            fourier_to_standard_basis(phi_fourier / minimum, 6),
        ]
    )
    return list(create_hexagon_design(period, rows)[1])


@pytest.mark.parametrize(
    ("shift_kind", "profile"),
    [
        ("groups", {"kind": "uniform"}),
        (
            "ellipse",
            {"kind": "tanh_power", "scale_layers": 2.5, "power": 2.0},
        ),
    ],
)
def test_standalone_hex_geometry_matches_current_finite_plan(
    shift_kind: str,
    profile: dict[str, object],
) -> None:
    mapping = _geometry_mapping(shift_kind, profile)
    blueprint = build_hex_layout(mapping)
    current = build_finite_geometry_plan(
        FiniteGeometrySpec(
            "hex",
            period=0.82,
            cavity_layers=2,
            cladding_layers=2,
            empty_buffer_periods=0.5,
        ),
        _current_holes(mapping, "cavity"),
        _current_holes(mapping, "cladding"),
        shift_geometry=CladdingShiftGeometry(shift_kind),
        shift_profile=CladdingShiftProfile.from_mapping(profile),
        x_shift_factor=0.05,
        y_shift_factor=0.08,
    )

    by_parent: dict[tuple[str, int, int], list] = {}
    for feature in blueprint.features:
        by_parent.setdefault(
            (feature.region, feature.parent_i, feature.parent_j), []
        ).append(feature)

    assert len(blueprint.features) == len(current.records)
    for placement in current.placements:
        region = "cavity" if placement.region == "bulk" else "cladding"
        features = sorted(
            by_parent[(region, placement.m, placement.n)],
            key=lambda item: item.local_index,
        )
        np.testing.assert_allclose(
            features[0].shift_um,
            placement.shift,
            rtol=0.0,
            atol=2e-15,
        )
        for feature, polygon in zip(features, placement.polygons, strict=True):
            np.testing.assert_allclose(
                np.asarray(feature.polygon.exterior.coords[:-1]),
                polygon,
                rtol=0.0,
                atol=2e-15,
            )

    np.testing.assert_allclose(
        np.asarray(blueprint.boundary.exterior.coords[:-1]),
        current.footprint.vertices,
        rtol=0.0,
        atol=2e-15,
    )
