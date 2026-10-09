import numpy as np
import pytest

from scripts.run_main import run_strip_1d


def _record(region, center, polygon, modulation_shift=(0.0, 0.0)):
    point = run_strip_1d.LatticePoint(
        i=2,
        j=-1,
        shell=2,
        x=float(center[0]),
        y=float(center[1]),
    )
    return {
        "region": region,
        "point": point,
        "polygon": np.asarray(polygon, dtype=float),
        "modulation_shift": np.asarray(modulation_shift, dtype=float),
        "marker": "preserved",
    }


def test_strip_cut_is_configured_for_the_slanted_side_family():
    assert run_strip_1d.STRIP_CUT_ANGLE == 60.0
    assert run_strip_1d.ALLOWED_STRIP_CUT_ANGLES == (0.0, 60.0)


def test_cut_basis_matches_global_tangent_and_normal():
    tangent, normal = run_strip_1d.strip_cut_basis(60.0)
    np.testing.assert_allclose(tangent, [0.5, np.sqrt(3.0) / 2.0], atol=1e-12)
    np.testing.assert_allclose(normal, [-np.sqrt(3.0) / 2.0, 0.5], atol=1e-12)
    np.testing.assert_allclose(np.dot(tangent, normal), 0.0, atol=1e-12)


def test_fourier_row_index_follows_cut_local_finite_axis():
    assert run_strip_1d.strip_fourier_row_index_coefficients(0.0) == (0, 1)
    assert run_strip_1d.strip_fourier_row_index_coefficients(60.0) == (-1, 0)


def test_sixty_degree_bulk_strip_has_all_expected_fourier_rows():
    from comsol_workflow.lattice_fourier_postprocess import strip_bulk_rows

    records = run_strip_1d.strip_bulk_cell_records(60.0)
    rows = strip_bulk_rows(
        records,
        run_strip_1d.BULK_RADIUS,
        row_index_coefficients=run_strip_1d.strip_fourier_row_index_coefficients(60.0),
    )

    assert rows == list(range(-run_strip_1d.BULK_RADIUS, run_strip_1d.BULK_RADIUS + 1))


def test_strip_row_fields_use_one_batched_interpolation(monkeypatch):
    from comsol_workflow import lattice_fourier_postprocess as postprocess

    calls = []

    def fake_interpolate(_field_points, _field, query_points):
        calls.append(query_points.copy())
        return query_points[:, 0] + 1j * query_points[:, 1]

    monkeypatch.setattr(postprocess, "interpolate_complex_field", fake_interpolate)
    row_centers = np.array([[0.0, 0.0], [0.41, 0.7]])
    rho_points = np.array([[-0.3, 0.0], [0.3, 0.2]])
    result = postprocess.extract_strip_bulk_row_fields(
        np.zeros((3, 2)),
        np.zeros(3, dtype=complex),
        rho_points,
        row_centers,
        kx_fraction=0.25,
        period=0.82,
    )

    unwrapped = row_centers[:, None, :] + rho_points[None, :, :]
    wrapped, shifts = postprocess.wrap_periodic_x(unwrapped.reshape(-1, 2), 0.82)
    expected = (
        np.exp(-2j * np.pi * 0.25 * shifts)
        * (wrapped[:, 0] + 1j * wrapped[:, 1])
    ).reshape(2, 2)
    np.testing.assert_allclose(result, expected)
    assert len(calls) == 1


def test_global_to_local_transform_is_rigid_and_invertible():
    points = np.array([[2.0, -1.0], [3.0, 0.0], [-0.25, 4.0]])
    local = run_strip_1d.points_to_strip_frame(points, 60.0)
    restored = run_strip_1d.points_from_strip_frame(local, 60.0)

    np.testing.assert_allclose(restored, points, atol=1e-12)
    np.testing.assert_allclose(
        np.linalg.norm(local[1] - local[0]),
        np.linalg.norm(points[1] - points[0]),
        atol=1e-12,
    )


def test_whole_record_transform_applies_to_bulk_and_cladding():
    center = np.array([2.0, -1.0])
    offsets = np.array([[0.2, 0.0], [0.0, 0.3], [-0.2, -0.1]])
    records = [
        _record("bulk", center, center + offsets),
        _record("cladding", center, center + offsets, modulation_shift=(0.1, -0.2)),
    ]

    transformed = run_strip_1d.records_to_strip_frame(records, 60.0)
    expected_polygon = run_strip_1d.points_to_strip_frame(center + offsets, 60.0)
    expected_center = run_strip_1d.points_to_strip_frame(center, 60.0)

    for record in transformed:
        np.testing.assert_allclose(record["polygon"], expected_polygon, atol=1e-12)
        np.testing.assert_allclose(
            [record["point"].x, record["point"].y],
            expected_center,
            atol=1e-12,
        )
        assert record["marker"] == "preserved"
    np.testing.assert_allclose(
        transformed[1]["modulation_shift"],
        run_strip_1d.vector_to_strip_frame(np.array([0.1, -0.2]), 60.0),
        atol=1e-12,
    )


def test_transform_preserves_hole_offsets_relative_to_cell_center():
    center = np.array([2.0, -1.0])
    polygon = center + np.array([[0.2, 0.0], [0.0, 0.3], [-0.2, -0.1]])
    record = _record("bulk", center, polygon)

    transformed = run_strip_1d.records_to_strip_frame([record], 60.0)[0]
    local_center = np.array([transformed["point"].x, transformed["point"].y])

    np.testing.assert_allclose(
        transformed["polygon"] - local_center,
        run_strip_1d.vector_to_strip_frame(polygon - center, 60.0),
        atol=1e-12,
    )


def test_zero_degree_cut_preserves_all_coordinates():
    center = np.array([1.0, 2.0])
    polygon = np.array([[1.0, 2.0], [2.0, 2.0], [1.0, 3.0]])
    record = _record("cladding", center, polygon, modulation_shift=(0.1, 0.2))

    transformed = run_strip_1d.records_to_strip_frame([record], 0.0)[0]

    np.testing.assert_allclose(transformed["polygon"], polygon)
    np.testing.assert_allclose(transformed["modulation_shift"], [0.1, 0.2])
    np.testing.assert_allclose([transformed["point"].x, transformed["point"].y], center)


@pytest.mark.parametrize("angle", [30.0, 120.0])
def test_unsupported_cut_angle_is_rejected(angle):
    with pytest.raises(ValueError, match="STRIP_CUT_ANGLE"):
        run_strip_1d.records_to_strip_frame([], angle)


def test_whole_geometry_transform_happens_before_local_strip_clipping():
    center = np.array([0.0, 0.0])
    record = _record(
        "bulk",
        center,
        np.array([[0.45, -0.05], [0.55, -0.05], [0.45, 0.05]]),
    )

    clipped = run_strip_1d.strip_records_for_cut([record], 60.0)
    transformed = run_strip_1d.records_to_strip_frame([record], 60.0)
    expected = run_strip_1d.clip_polygon_to_vertical_strip(transformed[0]["polygon"])

    assert expected is not None
    np.testing.assert_allclose(clipped[0]["polygon"], expected)
