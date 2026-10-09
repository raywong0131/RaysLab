import os
import json
from dataclasses import dataclass

import numpy as np

import mph
from jpype.types import JArray, JDouble, JInt, JString

from ..energy_recovery import integrate_field, interpolate_field

def _isotropic_tensor(value):
    return [value, "0", "0", "0", value, "0", "0", "0", value]


@dataclass(frozen=True)
class SimulationConfig:
    slab_height: str = "200 [nm]"
    wavelength: str = "1550 [nm]"
    air_height: str = "lda0"
    air_layer_top_factor: float = 1.0
    pml_layer_top_factor: float = 1.5
    air_refractive_index: str = "1"
    air_extinction_coefficient: str = "0"
    slab_refractive_index: str = "3.3"
    slab_extinction_coefficient: str = "0"
    eigenfrequency_shift: str = "c_const/1.55[um]"
    eigenmode_count: int = 8
    air_cutplane_z: str = "0.75*H_air"
    selection_tolerance: str = "10 [nm]"
    mode_type: str = "TE"

    @staticmethod
    def _layer_distance_expr(factor):
        if factor == 1.0:
            return "H_air"
        return f"H_air*{factor:g}"

    @property
    def air_layer_distance(self):
        return self._layer_distance_expr(self.air_layer_top_factor)

    @property
    def pml_layer_distance(self):
        return self._layer_distance_expr(self.pml_layer_top_factor)


class SimulationRun:
    def __init__(self, config=None):
        self.config = config or SimulationConfig()
        self.client = mph.start()
        self.model = self.client.create("Model")
        self._plotting_initialized = False
    
    def __del__(self):
        self.clear()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear()
        return False

    def clear(self):
        if self.model is not None:
            self.client.remove(self.model)
            self.model = None

    def build_and_run(self, a, holes, k=None):
        """Build and run the COMSOL eigenfrequency simulation.

        Parameters
        ----------
        a : float
            Lattice constant in um.
        holes : list of array-like, shape (N_holes, M_verts, 2)
            Polygon vertex coordinates (in um) for each hole in the unit cell.
            The hole count and vertex count can vary by design.
        k : dict, optional
            Bloch k-vector as fractions of G, e.g. ``{'kx': 0.1, 'ky': 0.0}``.
        """
        jmodel = self.model.java  # com.comsol.model.Model
        config = self.config

        jmodel.param().set("a", f"{a} [um]")
        jmodel.param().set("H", config.slab_height)
        jmodel.param().set("H_air", config.air_height)
        jmodel.param().set("lda0", config.wavelength)
        jmodel.param().set("selection_tol", config.selection_tolerance)
        jmodel.param().set("G", "2*pi/a")
        
        if k is not None:
            jmodel.param().set("kx", f"{k['kx']}*G")
            jmodel.param().set("ky", f"{k['ky']}*G")
        else:
            jmodel.param().set("kx", f"0*G")
            jmodel.param().set("ky", f"0*G")

        jmodel.component().create("comp1", True)
        jmodel.component("comp1").geom().create("geom1", 3)
        jmodel.component("comp1").mesh().create("mesh1")

        jmodel.component("comp1").geom("geom1").lengthUnit("um")
        jmodel.component("comp1").geom("geom1").create("wp1", "WorkPlane")
        jmodel.component("comp1").geom("geom1").feature("wp1").set("unite", True)

        square_table = [
            ["a/2", "a/2"],
            ["-a/2", "a/2"],
            ["-a/2", "-a/2"],
            ["a/2", "-a/2"],
        ]

        jmodel.component("comp1").geom("geom1").feature("wp1").geom().create("pol1", "Polygon")
        jmodel.component("comp1").geom("geom1").feature("wp1").geom().feature("pol1").label("Square")
        jmodel.component("comp1").geom("geom1").feature("wp1").geom().feature("pol1").set("source", "table")
        jmodel.component("comp1").geom("geom1").feature("wp1").geom().feature("pol1").set(
            "table",
            JArray(JString, 2)(square_table),
        )

        hole_tables = {
            f"pol{i + 2}": (
                f"Hole_{i + 1}",
                [[f"{x} [um]", f"{y} [um]"] for (x, y) in hole],
            )
            for i, hole in enumerate(holes)
        }

        for tag, (label, table) in hole_tables.items():
            jmodel.component("comp1").geom("geom1").feature("wp1").geom().create(tag, "Polygon")
            jmodel.component("comp1").geom("geom1").feature("wp1").geom().feature(tag).label(label)
            jmodel.component("comp1").geom("geom1").feature("wp1").geom().feature(tag).set("source", "table")
            jmodel.component("comp1").geom("geom1").feature("wp1").geom().feature(tag).set(
                "table",
                JArray(JString, 2)(table),
            )

        jmodel.component("comp1").geom("geom1").create("ext1", "Extrude")
        jmodel.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "H/2", 0)
        jmodel.component("comp1").geom("geom1").feature("ext1").selection("input").set("wp1")

        jmodel.component("comp1").geom("geom1").create("wp2", "WorkPlane")
        jmodel.component("comp1").geom("geom1").feature("wp2").set("quickz", "H/2")
        jmodel.component("comp1").geom("geom1").feature("wp2").set("unite", True)

        top_square_tag = f"pol{len(hole_tables) + 2}"
        jmodel.component("comp1").geom("geom1").feature("wp2").geom().create(top_square_tag, "Polygon")
        jmodel.component("comp1").geom("geom1").feature("wp2").geom().feature(top_square_tag).label("Square 1")
        jmodel.component("comp1").geom("geom1").feature("wp2").geom().feature(top_square_tag).set("source", "table")
        jmodel.component("comp1").geom("geom1").feature("wp2").geom().feature(top_square_tag).set(
            "table",
            JArray(JString, 2)(square_table),
        )

        jmodel.component("comp1").geom("geom1").create("ext2", "Extrude")
        jmodel.component("comp1").geom("geom1").feature("ext2").set(
            "distance",
            JArray(JString, 1)([config.air_layer_distance, config.pml_layer_distance]),
        )
        jmodel.component("comp1").geom("geom1").feature("ext2").set(
            "scale",
            JArray(JDouble, 2)([[1.0, 1.0], [1.0, 1.0]]),
        )
        jmodel.component("comp1").geom("geom1").feature("ext2").set(
            "displ",
            JArray(JDouble, 2)([[0.0, 0.0], [0.0, 0.0]]),
        )
        jmodel.component("comp1").geom("geom1").feature("ext2").set(
            "twist",
            JArray(JInt, 1)([0, 0]),
        )
        jmodel.component("comp1").geom("geom1").feature("ext2").selection("input").set("wp2")

        jmodel.component("comp1").geom("geom1").run()

        # -----------------------------------------------------------------------
        # Position-based selections.
        #
        # Boundary selections are seeded by small balls at known points on the
        # square unit-cell exterior.  Each seed is restricted to the exterior
        # boundary selection before continuous-tangent grouping is enabled, which
        # prevents coplanar internal faces from being selected accidentally.
        #
        # Geometry z layout (values in um, evaluated from model parameters):
        #   z = 0                                      bottom symmetry plane
        #   z = H/2                                    top of slab / air start
        #   z = H/2+air_layer_top_factor*H_air         bottom of PML layer
        #   z = H/2+pml_layer_top_factor*H_air         top of PML
        #
        # Square wall centers in xy:
        #   pc1: right (a/2, 0) + left   (-a/2, 0)
        #   pc2: top   (0, a/2) + bottom (0, -a/2)
        # -----------------------------------------------------------------------

        a_val     = float(jmodel.param().evaluate("a",     "um"))
        H_val     = float(jmodel.param().evaluate("H",     "um"))
        H_air_val = float(jmodel.param().evaluate("H_air", "um"))
        selection_tol = float(jmodel.param().evaluate("selection_tol", "um"))

        z_top    = H_val / 2 + H_air_val * config.pml_layer_top_factor
        z_pml_lo = H_val / 2 + H_air_val * config.air_layer_top_factor
        z_air_mid = H_val / 2 + H_air_val * config.air_layer_top_factor / 2
        large    = a_val * 10    # "infinity" along unconstrained axes

        def _box_sel(name, label, dim, xmn, xmx, ymn, ymx, zmn, zmx, cond="inside"):
            """Create a named Box selection."""
            jmodel.component("comp1").selection().create(name, "Box")
            s = jmodel.component("comp1").selection(name)
            s.set("entitydim", JInt(dim))
            s.set("xmin", xmn); s.set("xmax", xmx)
            s.set("ymin", ymn); s.set("ymax", ymx)
            s.set("zmin", zmn); s.set("zmax", zmx)
            s.set("condition", cond)
            if label:
                s.label(label)

        def _exterior_ball_sel(name, label, x, y, z):
            """Seed an exterior boundary face and expand by tangent continuity."""
            jmodel.component("comp1").selection().create(name, "Ball")
            s = jmodel.component("comp1").selection(name)
            s.set("entitydim", JInt(2))
            s.set("inputent", "selections")
            s.set("input", JArray(JString, 1)(["sel_exterior_boundaries"]))
            s.set("condition", "intersects")
            s.set("groupcontang", "on")
            s.set("posx", x)
            s.set("posy", y)
            s.set("posz", z)
            s.set("r", selection_tol)
            if label:
                s.label(label)

        # Boundaries exterior to the full 3D geometry. COMSOL 6.3 does not
        # support the exterioroutside subtype, so the small seed balls and
        # tangent grouping below localize this to the intended outside faces.
        jmodel.component("comp1").selection().create("sel_exterior_boundaries")
        _s = jmodel.component("comp1").selection("sel_exterior_boundaries")
        _s.geom("geom1", JInt(3), JInt(2), JArray(JString, 1)(["exterior"]))
        _s.all()
        _s.label("exterior_boundaries")

        # Bottom symmetry plane and top PML scattering boundary.
        _exterior_ball_sel("sel_bottom", "bottom", 0.0, 0.0, 0.0)
        _exterior_ball_sel("sel_top_face", "top", 0.0, 0.0, z_top)

        # Silicon slab domains: all domains in the slab z-range, minus the hole
        # air domains defined by the input hole polygons.
        _box_sel("sel_slab_layer", "slab_layer", 3,
                 -large, large, -large, large, -selection_tol, H_val / 2 + selection_tol)

        hole_domain_selections = []
        for i, hole in enumerate(holes):
            hole_center = np.asarray(hole, dtype=float).mean(axis=0)
            sname = f"sel_hole_dom_{i}"
            jmodel.component("comp1").selection().create(sname, "Ball")
            _s = jmodel.component("comp1").selection(sname)
            _s.set("entitydim", JInt(3))
            _s.set("condition", "intersects")
            _s.set("posx", float(hole_center[0]))
            _s.set("posy", float(hole_center[1]))
            _s.set("posz", 0.0)
            _s.set("r", selection_tol)
            hole_domain_selections.append(sname)

        if hole_domain_selections:
            jmodel.component("comp1").selection().create("sel_hole_doms", "Union")
            _s = jmodel.component("comp1").selection("sel_hole_doms")
            _s.set("entitydim", JInt(3))
            _s.set("input", JArray(JString, 1)(hole_domain_selections))
            _s.label("hole_doms")

            jmodel.component("comp1").selection().create("sel_slab_dom", "Difference")
            _s = jmodel.component("comp1").selection("sel_slab_dom")
            _s.set("entitydim", JInt(3))
            _s.set("add", JArray(JString, 1)(["sel_slab_layer"]))
            _s.set("subtract", JArray(JString, 1)(["sel_hole_doms"]))
            _s.label("slab_dom")
        else:
            jmodel.component("comp1").selection().create("sel_slab_dom", "Union")
            _s = jmodel.component("comp1").selection("sel_slab_dom")
            _s.set("entitydim", JInt(3))
            _s.set("input", JArray(JString, 1)(["sel_slab_layer"]))
            _s.label("slab_dom")

        # PML domain: the outermost z layer (z from z_pml_lo to z_top).
        _box_sel("sel_pml_dom", "pml_dom", 3,
                 -large, large, -large, large, z_pml_lo - selection_tol, z_top + selection_tol)

        # ---- Periodic wall face selections ----
        # One seed per disconnected side; each seed expands over the continuous
        # tangent exterior wall, then opposite sides are unioned into a PC pair.
        half_width = a_val / 2
        wall_seeds = {
            "pc1": [
                ( half_width, 0.0, z_air_mid),
                (-half_width, 0.0, z_air_mid),
            ],
            "pc2": [
                (0.0,  half_width, z_air_mid),
                (0.0, -half_width, z_air_mid),
            ],
        }

        for pc_key, seeds in wall_seeds.items():
            sub_names = []
            for i, (x, y, z) in enumerate(seeds):
                sname = f"sel_{pc_key}_w{i}"
                _exterior_ball_sel(sname, None, x, y, z)
                sub_names.append(sname)
            jmodel.component("comp1").selection().create(f"sel_{pc_key}", "Union")
            jmodel.component("comp1").selection(f"sel_{pc_key}").set("entitydim", JInt(2))
            jmodel.component("comp1").selection(f"sel_{pc_key}").label(f"{pc_key}_walls")
            jmodel.component("comp1").selection(f"sel_{pc_key}").set(
                "input", JArray(JString, 1)(sub_names))

        jmodel.component("comp1").material().create("mat1", "Common")
        jmodel.component("comp1").material().create("mat2", "Common")
        jmodel.component("comp1").material("mat1").propertyGroup().create("RefractiveIndex", "Refractive Index")
        jmodel.component("comp1").material("mat2").selection().named("sel_slab_dom")
        jmodel.component("comp1").material("mat2").propertyGroup().create("RefractiveIndex", "Refractive Index")

        jmodel.component("comp1").coordSystem().create("pml1", "PML")
        jmodel.component("comp1").coordSystem("pml1").selection().named("sel_pml_dom")

        jmodel.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
        jmodel.component("comp1").physics("ewfd").create("sctr1", "Scattering", 2)
        jmodel.component("comp1").physics("ewfd").feature("sctr1").selection().named("sel_top_face")

        mode_type = config.mode_type.upper()
        if mode_type == "TM":
            bottom_feature = "pec_bottom"
            bottom_condition = "PerfectElectricConductor"
        elif mode_type == "TE":
            bottom_feature = "pmc_bottom"
            bottom_condition = "PerfectMagneticConductor"
        else:
            raise ValueError("mode_type must be 'TE' or 'TM'")

        jmodel.component("comp1").physics("ewfd").create(bottom_feature, bottom_condition, 2)
        jmodel.component("comp1").physics("ewfd").feature(bottom_feature).selection().named("sel_bottom")

        jmodel.component("comp1").physics("ewfd").create("pc1", "PeriodicCondition", 2)
        jmodel.component("comp1").physics("ewfd").feature("pc1").selection().named("sel_pc1")

        jmodel.component("comp1").physics("ewfd").create("pc2", "PeriodicCondition", 2)
        jmodel.component("comp1").physics("ewfd").feature("pc2").selection().named("sel_pc2")

        if k is not None:
            jmodel.component("comp1").physics("ewfd").feature("pc1").set("PeriodicType", "Floquet")
            jmodel.component("comp1").physics("ewfd").feature("pc1").set("kFloquet", JArray(JString, 1)(["kx", "ky", "0"]))

            jmodel.component("comp1").physics("ewfd").feature("pc2").set("PeriodicType", "Floquet")
            jmodel.component("comp1").physics("ewfd").feature("pc2").set("kFloquet", JArray(JString, 1)(["kx", "ky", "0"]))

        jmodel.component("comp1").material("mat1").label("Air")
        jmodel.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("n", "")
        jmodel.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("ki", "")
        jmodel.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set(
            "n",
            JArray(JString, 1)(_isotropic_tensor(config.air_refractive_index)),
        )
        jmodel.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set(
            "ki",
            JArray(JString, 1)(_isotropic_tensor(config.air_extinction_coefficient)),
        )

        jmodel.component("comp1").material("mat2").label("Mat")
        jmodel.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("n", "")
        jmodel.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("ki", "")
        jmodel.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("n", "")
        jmodel.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("ki", "")
        jmodel.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set(
            "n",
            JArray(JString, 1)(_isotropic_tensor(config.slab_refractive_index)),
        )
        jmodel.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set(
            "ki",
            JArray(JString, 1)(_isotropic_tensor(config.slab_extinction_coefficient)),
        )

        jmodel.study().create("std1")
        jmodel.study("std1").create("eig", "Eigenfrequency")
        jmodel.study("std1").feature("eig").set("shift", config.eigenfrequency_shift)
        jmodel.study("std1").feature("eig").set("neigsactive", True)
        jmodel.study("std1").feature("eig").set("neigs", JInt(config.eigenmode_count))

        jmodel.component("comp1").mesh("mesh1").run()

        jmodel.study("std1").createAutoSequences("all")

        jmodel.sol("sol1").runAll()

        # output 1: eigenfrequencies
        jmodel.result().numerical().create("gev1", "EvalGlobal")
        jmodel.result().numerical("gev1").label("Eigenfrequencies (ewfd)")
        jmodel.result().numerical("gev1").set("data", "dset1")
        # angular freq = omega + i damp, Q = omega / (2*damp)
        # freq = angular freq / (2*pi), so divide by 2*pi to get freq and damp in THz
        jmodel.result().numerical("gev1").set("expr", ["ewfd.omega/2/pi", "ewfd.damp/2/pi", "ewfd.Qfactor"])
        jmodel.result().numerical("gev1").set("unit", ["THz", "THz", "1"])

        jmodel.result().table().create("tbl1", "Table")
        jmodel.result().numerical("gev1").set("table", "tbl1")
        jmodel.result().numerical("gev1").run()
        jmodel.result().numerical("gev1").setResult()

        # output 2: 2D field
        # center
        jmodel.result().dataset().create("cpl1", "CutPlane")
        jmodel.result().dataset("cpl1").set("quickplane", "xy")
        # air
        jmodel.result().dataset().create("cpl2", "CutPlane")
        jmodel.result().dataset("cpl2").set("quickplane", "xy")
        jmodel.result().dataset("cpl2").set("quickz", config.air_cutplane_z)
        # yz
        jmodel.result().dataset().create("cpl3", "CutPlane")
        jmodel.result().dataset("cpl3").set("quickplane", "yz")
        # xz
        jmodel.result().dataset().create("cpl4", "CutPlane")
        jmodel.result().dataset("cpl4").set("quickplane", "xz")

        self.plane_datasets = {
            "center": "cpl1",
            "air":    "cpl2",
            "yz":     "cpl3",
            "xz":     "cpl4",
        }
        jmodel.result().numerical().create("int1", "Interp")

        return

    def setup_plotting(self):
        """Lazily initialize COMSOL plot groups and export objects.
        Safe to call multiple times; only runs once per SimulationRun instance.
        """
        if self._plotting_initialized:
            return
        jmodel = self.model.java

        jmodel.result().create("pg1", "PlotGroup2D")
        jmodel.result("pg1").label("2D Field (ewfd)")
        jmodel.result("pg1").run()
        jmodel.result("pg1").create("surf1", "Surface")
        jmodel.result("pg1").feature("surf1").set("expr", "ewfd.normH")
        jmodel.result("pg1").run()
        jmodel.result().export().create("plot1", "pg1", "surf1", "Plot")
        jmodel.result().export().create("img1", "pg1", "Image")
        jmodel.result().export("img1").set("target", "file")

        jmodel.result().create("pg2", "PlotGroup3D")
        jmodel.result("pg2").run()
        jmodel.result("pg2").label("3D Plot")
        jmodel.result("pg2").create("mslc1", "Multislice")
        jmodel.result("pg2").feature("mslc1").set("multiplanexmethod", "coord")
        jmodel.result("pg2").feature("mslc1").set("xcoord", JInt(0))
        jmodel.result("pg2").feature("mslc1").set("multiplaneymethod", "coord")
        jmodel.result("pg2").feature("mslc1").set("ycoord", JInt(0))
        jmodel.result("pg2").feature("mslc1").set("multiplanezmethod", "coord")
        jmodel.result("pg2").feature("mslc1").set("zcoord", JInt(0))
        jmodel.result("pg2").feature("mslc1").set("expr", "ewfd.normE")
        jmodel.result("pg2").run()
        jmodel.result().export().create("img2", "pg2", "Image")
        jmodel.result().export("img2").set("target", "file")

        self._plotting_initialized = True

    def export_eigenfrequencies(self, save_path):
        os.makedirs(save_path, exist_ok=True)
        jmodel = self.model.java  # com.comsol.model.Model
        jmodel.result().table("tbl1").save(os.path.join(save_path, "eigenfrequencies.txt"))
    
    def get_eigenfrequencies(self):
        jmodel = self.model.java
        jmodel.result().numerical("gev1").run()
        data = np.asarray(jmodel.result().numerical("gev1").getReal()).T  # (expr, N) -> (N, expr)
        return data

    def export_2d_fields(self, eigenmode_idx, expr, expr_name, save_path, plane="center", export_data=True, export_image=True):
        self.setup_plotting()
        os.makedirs(save_path, exist_ok=True)
        jmodel = self.model.java  # com.comsol.model.Model
        jmodel.result("pg1").run()
        jmodel.result("pg1").set("data", self.plane_datasets[plane])
        jmodel.result("pg1").set("looplevel", JArray(JInt, 1)([eigenmode_idx + 1]))
        jmodel.result("pg1").feature("surf1").set("expr", expr)
        jmodel.result("pg1").run()
        if export_data:
            jmodel.result().export("plot1").set("filename", os.path.join(save_path, f"{eigenmode_idx:02d}_{expr_name}_{plane}_2d.txt"))
            jmodel.result().export("plot1").run()
        if export_image:
            from ..field_plotting import save_hz_simulation_plot

            filename = os.path.join(save_path, f"{eigenmode_idx:02d}_{expr_name}_{plane}_2d.png")
            if save_hz_simulation_plot(self, eigenmode_idx, expr, filename, plane):
                return
            jmodel.result().export("img1").set("pngfilename", os.path.join(save_path, f"{eigenmode_idx:02d}_{expr_name}_{plane}_2d.png"))
            jmodel.result().export("img1").run()
    
    # Interp Doc: https://doc.comsol.com/6.3/docserver/#!/com.comsol.help.comsol/comsol_api_results.52.075.html
    def get_2d_fields(self, eigenmode_idx, expr, plane="center"):
        jmodel = self.model.java
        jmodel.result().numerical("int1").set("data", self.plane_datasets[plane])
        jmodel.result().numerical("int1").set("expr", JArray(JString, 1)([expr]))
        jmodel.result().numerical("int1").set("solnum", JArray(JInt, 1)([eigenmode_idx + 1]))
        jmodel.result().numerical("int1").run()

        coords = np.asarray(jmodel.result().numerical("int1").getCoordinates()).T # (2, N) -> (N, 2)
        values = np.asarray(jmodel.result().numerical("int1").getData()).reshape(-1) # (expr, solnum, coordinates) -> (N,)
        return coords, values

    def compute_polarization(self, eigenmode_idx):
        """Compute Fourier-projected polarization (cx, cy) for a given eigenmode.

        Recovers kx and ky from the model parameters.

        For the eigenmode at k = (kx, ky):
            cx = integrate(Ex * exp(i*kx*x + i*ky*y)) / sqrt(Area_x)
            cy = integrate(Ey * exp(i*kx*x + i*ky*y)) / sqrt(Area_y)

        where Ex, Ey are complex fields on the air plane and Area_x/y is the
        total integration-domain area of the respective mesh.

        Returns
        -------
        cx, cy : complex
            Jones-vector components for x and y polarisation.
        """
        jmodel = self.model.java

        # param().evaluate() returns base unit in SI (1/m)
        # coordinates from get_2d_fields
        # are in the model's length unit (um), so convert kx/ky to 1/um.
        kx = float(jmodel.param().evaluate("kx", "um^-1"))
        ky = float(jmodel.param().evaluate("ky", "um^-1"))

        # Ex field on air plane – real and imaginary parts on potentially
        # different meshes; interpolate im onto re's mesh before combining.
        pts_ex_re, ex_re = self.get_2d_fields(eigenmode_idx, "ewfd.Ex",       "air")
        pts_ex_im, ex_im = self.get_2d_fields(eigenmode_idx, "ewfd.Ex*(-i)",  "air")
        ex_comp = ex_re + 1j * interpolate_field(pts_ex_im, ex_im, pts_ex_re)

        # Ey field on air plane
        pts_ey_re, ey_re = self.get_2d_fields(eigenmode_idx, "ewfd.Ey",       "air")
        pts_ey_im, ey_im = self.get_2d_fields(eigenmode_idx, "ewfd.Ey*(-i)",  "air")
        ey_comp = ey_re + 1j * interpolate_field(pts_ey_im, ey_im, pts_ey_re)

        # Phase factor: kx [1/um] * x [um] and ky [1/um] * y [um] -> dimensionless
        phase_ex = np.exp(1j * (kx * pts_ex_re[:, 0] + ky * pts_ex_re[:, 1]))
        phase_ey = np.exp(1j * (kx * pts_ey_re[:, 0] + ky * pts_ey_re[:, 1]))

        area_x = float(np.real(integrate_field(pts_ex_re, np.ones(len(pts_ex_re)))))
        area_y = float(np.real(integrate_field(pts_ey_re, np.ones(len(pts_ey_re)))))

        cx = integrate_field(pts_ex_re, ex_comp * phase_ex) / np.sqrt(area_x)
        cy = integrate_field(pts_ey_re, ey_comp * phase_ey) / np.sqrt(area_y)

        return complex(cx), complex(cy)

    def export_3d_fields(self, eigenmode_idx, expr, expr_name, save_path):
        self.setup_plotting()
        os.makedirs(save_path, exist_ok=True)
        jmodel = self.model.java  # com.comsol.model.Model
        jmodel.result("pg2").run()
        jmodel.result("pg2").set("data", "dset1")
        jmodel.result("pg2").set("looplevel", JArray(JInt, 1)([eigenmode_idx + 1]))
        jmodel.result("pg2").feature("mslc1").set("expr", expr)
        jmodel.result("pg2").run()
        jmodel.result().export("img2").set("pngfilename", os.path.join(save_path, f"{eigenmode_idx:02d}_{expr_name}_3d.png"))
        jmodel.result().export("img2").run()
