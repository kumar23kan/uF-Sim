#!/usr/bin/env python3
"""
Local simulation runner for microfluidic CFD.
Reads a job .json exported from uf_sim_ui.py and runs the full pipeline:
  STL → Gmsh mesh → FEniCSx solver → ParaView output (.xdmf / .vtu)

Usage:
    conda activate uf-sim
    python3 run_simulation.py <job.json>
    python3 run_simulation.py                  # opens file picker
"""

import argparse, json, os, sys, datetime, zipfile
import numpy as np

# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Microfluidic CFD local runner")
    p.add_argument("job", nargs="?", help="Path to job .json file")
    return p.parse_args()


def load_job(path=None):
    if path is None:
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk(); root.withdraw()
            path = filedialog.askopenfilename(
                title="Select job .json file",
                filetypes=[("JSON job file", "*.json")]
            )
            root.destroy()
        except Exception:
            sys.exit("No job file provided. Usage: python3 run_simulation.py <job.json>")
    if not path or not os.path.isfile(path):
        sys.exit(f"Job file not found: {path}")
    with open(path) as f:
        return json.load(f)


# ── Mesh ──────────────────────────────────────────────────────────────────────

MESH_SIZES    = {"Coarse": 0.5, "Medium": 0.2, "Fine": 0.08, "Very Fine": 0.04}
BL_FRACTIONS  = {"Coarse": 0.10, "Medium": 0.08, "Fine": 0.06, "Very Fine": 0.05}
BL_LAYERS     = 5
BL_RATIO      = 1.3

INLET_ID = 1; OUTLET_ID = 2; WALL_ID = 3; FLUID_ID = 4
BC_IDS = {"inlet": INLET_ID, "outlet": OUTLET_ID, "wall": WALL_ID, "fluid": FLUID_ID}


def build_mesh(stl_path, resolution, out_dir):
    import gmsh
    msh_path = os.path.join(out_dir, "channel.msh")
    gmsh.initialize()
    gmsh.model.add("channel")
    gmsh.option.setNumber("General.Verbosity", 2)

    gmsh.merge(stl_path)

    angle = 40 * (3.14159265 / 180)
    gmsh.model.mesh.classifySurfaces(angle, True, True, 3.14159265)
    gmsh.model.mesh.createGeometry()

    surfaces = gmsh.model.getEntities(2)
    s_tags   = [s[1] for s in surfaces]
    print(f"  Surfaces found: {len(s_tags)}")

    loop = gmsh.model.geo.addSurfaceLoop(s_tags)
    vol  = gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()

    # Classify inlet / outlet / wall by longest axis bounding box
    bb      = gmsh.model.getBoundingBox(-1, -1)
    x0,y0,z0,x1,y1,z1 = bb
    extents = {'x': (x0,x1), 'y': (y0,y1), 'z': (z0,z1)}
    flow_ax = max(extents, key=lambda k: extents[k][1] - extents[k][0])
    lo, hi  = extents[flow_ax]
    ax_idx  = {'x': 0, 'y': 1, 'z': 2}[flow_ax]
    tol     = (hi - lo) * 0.06

    cross_dims = [extents[ax][1] - extents[ax][0] for ax in extents if ax != flow_ax]
    min_cross  = min(cross_dims)

    inlet_tags, outlet_tags, wall_tags = [], [], []
    for dim, tag in surfaces:
        sbb = gmsh.model.getBoundingBox(dim, tag)
        slo, shi = sbb[ax_idx], sbb[ax_idx + 3]
        if shi < lo + tol:
            inlet_tags.append(tag)
        elif slo > hi - tol:
            outlet_tags.append(tag)
        else:
            wall_tags.append(tag)

    if not inlet_tags:  inlet_tags  = [s_tags[0]]
    if not outlet_tags: outlet_tags = [s_tags[-1]]
    wall_tags = wall_tags or [t for t in s_tags if t not in inlet_tags + outlet_tags]

    gmsh.model.addPhysicalGroup(2, inlet_tags,  INLET_ID); gmsh.model.setPhysicalName(2, INLET_ID,  "inlet")
    gmsh.model.addPhysicalGroup(2, outlet_tags, OUTLET_ID);gmsh.model.setPhysicalName(2, OUTLET_ID, "outlet")
    if wall_tags:
        gmsh.model.addPhysicalGroup(2, wall_tags, WALL_ID);gmsh.model.setPhysicalName(2, WALL_ID,   "walls")
    gmsh.model.addPhysicalGroup(3, [vol],        FLUID_ID);gmsh.model.setPhysicalName(3, FLUID_ID,  "fluid")

    size    = MESH_SIZES.get(resolution, 0.2)
    bl_frac = BL_FRACTIONS.get(resolution, 0.08)
    effective_size = min(size, min_cross * 0.25)
    bl_thickness   = min_cross * bl_frac

    gmsh.option.setNumber("Mesh.MeshSizeMax", effective_size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", effective_size * 0.05)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    if wall_tags:
        field = gmsh.model.mesh.field
        bl_id = field.add("BoundaryLayer")
        field.setNumbers(bl_id, "FacesList", wall_tags)
        field.setNumber( bl_id, "Size",      bl_thickness)
        field.setNumber( bl_id, "Ratio",     BL_RATIO)
        field.setNumber( bl_id, "NbLayers",  BL_LAYERS)
        field.setNumber( bl_id, "Quads",     0)
        field.setAsBoundaryLayer(bl_id)
        total_bl = bl_thickness * sum(BL_RATIO**i for i in range(BL_LAYERS))
        print(f"  Boundary layer: {BL_LAYERS} layers, "
              f"first={bl_thickness:.4f}, total={total_bl:.4f} "
              f"({total_bl/min_cross*100:.1f}% of channel depth)")

    print(f"  Bulk mesh size: {effective_size:.4f} (channel depth: {min_cross:.4f})")
    print(f"  Generating mesh ({resolution})...")
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.write(msh_path)
    gmsh.finalize()
    print(f"  Mesh written: {msh_path}")
    return msh_path


# ── Solvers ───────────────────────────────────────────────────────────────────

def solve_stokes(domain, facet_tags, bc_ids, params, out_dir):
    from dolfinx import fem, io
    from dolfinx.fem import functionspace, Constant
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    from petsc4py import PETSc
    import basix.ufl, ufl

    mu_val    = float(params.get("Fluid Viscosity (Pa·s)",  0.001))
    inlet_vel = float(params.get("Inlet Velocity (m/s)",     0.001))

    P2 = basix.ufl.element("Lagrange", domain.topology.cell_name(), 2,
                             shape=(domain.geometry.dim,))
    P1 = basix.ufl.element("Lagrange", domain.topology.cell_name(), 1)
    W  = functionspace(domain, basix.ufl.mixed_element([P2, P1]))

    (u, p) = ufl.TrialFunctions(W)
    (v, q) = ufl.TestFunctions(W)
    mu = Constant(domain, PETSc.ScalarType(mu_val))
    f  = Constant(domain, PETSc.ScalarType((0.0, 0.0, 0.0)))

    a = (mu * ufl.inner(ufl.grad(u), ufl.grad(v))
         - ufl.inner(p, ufl.div(v))
         + ufl.inner(ufl.div(u), q)) * ufl.dx
    L = ufl.inner(f, v) * ufl.dx

    V_sub, _ = W.sub(0).collapse()
    Q_sub, _ = W.sub(1).collapse()

    u_in   = Constant(domain, PETSc.ScalarType((0.0, 0.0, inlet_vel)))
    in_dof = fem.locate_dofs_topological((W.sub(0), V_sub), 2, facet_tags.find(bc_ids["inlet"]))
    bc_in  = fem.dirichletbc(u_in, in_dof, W.sub(0))

    u_wall = Constant(domain, PETSc.ScalarType((0.0, 0.0, 0.0)))
    wall_f = facet_tags.find(bc_ids["wall"])
    bcs    = [bc_in]
    if len(wall_f) > 0:
        bcs.append(fem.dirichletbc(u_wall,
            fem.locate_dofs_topological((W.sub(0), V_sub), 2, wall_f), W.sub(0)))

    p_out  = Constant(domain, PETSc.ScalarType(0.0))
    op_dof = fem.locate_dofs_topological((W.sub(1), Q_sub), 2, facet_tags.find(bc_ids["outlet"]))
    bcs.append(fem.dirichletbc(p_out, op_dof, W.sub(1)))

    print("  Solving Stokes...")
    wh  = LinearProblem(a, L, bcs=bcs,
          petsc_options={"ksp_type": "minres", "pc_type": "fieldsplit",
                         "pc_fieldsplit_type": "schur"}).solve()

    u_h = wh.sub(0).collapse(); u_h.name = "velocity"
    p_h = wh.sub(1).collapse(); p_h.name = "pressure"

    ip = fem.locate_dofs_topological(Q_sub, 2, facet_tags.find(bc_ids["inlet"]))
    op = fem.locate_dofs_topological(Q_sub, 2, facet_tags.find(bc_ids["outlet"]))
    dP = float(np.mean(p_h.x.array[ip]) - np.mean(p_h.x.array[op])) if len(ip) else 0.0
    print(f"  Pressure drop: {dP:.4f} Pa")

    for fname, fn in [("velocity.xdmf", u_h), ("pressure.xdmf", p_h)]:
        with io.XDMFFile(MPI.COMM_WORLD, os.path.join(out_dir, fname), "w") as xf:
            xf.write_mesh(domain); xf.write_function(fn)

    return u_h, p_h, dP


def solve_mixing(domain, facet_tags, bc_ids, params, out_dir, u_h):
    from dolfinx import fem, io
    from dolfinx.fem import functionspace, Constant
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    from petsc4py import PETSc
    import ufl

    D_val  = float(params.get("Diffusion Coefficient (m²/s)",   1e-9))
    c1_val = float(params.get("Inlet 1 Concentration (mol/m³)", 1.0))
    c2_val = float(params.get("Inlet 2 Concentration (mol/m³)", 0.0))

    V_c   = functionspace(domain, ("Lagrange", 1))
    c, ph = ufl.TrialFunction(V_c), ufl.TestFunction(V_c)
    D     = Constant(domain, PETSc.ScalarType(D_val))
    h     = ufl.CellDiameter(domain)
    u_m   = ufl.sqrt(ufl.inner(u_h, u_h) + 1e-16)
    tau   = h / (2.0 * u_m)

    a = (D * ufl.inner(ufl.grad(c), ufl.grad(ph))
         + ufl.dot(u_h, ufl.grad(c)) * ph
         + tau * ufl.dot(u_h, ufl.grad(c)) * ufl.dot(u_h, ufl.grad(ph))) * ufl.dx
    L = Constant(domain, PETSc.ScalarType(0.0)) * ph * ufl.dx

    in_f   = facet_tags.find(bc_ids["inlet"])
    in_d   = fem.locate_dofs_topological(V_c, 2, in_f)
    coords = domain.geometry.x[in_d]
    y_mid  = (coords[:, 1].max() + coords[:, 1].min()) / 2.0
    d1 = in_d[coords[:, 1] >= y_mid]
    d2 = in_d[coords[:, 1] <  y_mid]

    print("  Solving advection-diffusion...")
    c_h = LinearProblem(a, L,
        bcs=[fem.dirichletbc(PETSc.ScalarType(c1_val), d1, V_c),
             fem.dirichletbc(PETSc.ScalarType(c2_val), d2, V_c)],
        petsc_options={"ksp_type": "gmres", "pc_type": "ilu"}).solve()
    c_h.name = "concentration"

    cv    = c_h.x.array
    c_ref = (c1_val + c2_val) / 2.0
    var   = float(np.mean((cv - c_ref) ** 2))
    m_var = ((c1_val - c2_val) / 2.0) ** 2
    eff   = 1.0 - var / m_var if m_var > 0 else 1.0
    print(f"  Mixing efficiency: {eff*100:.1f}%")

    with io.XDMFFile(MPI.COMM_WORLD, os.path.join(out_dir, "concentration.xdmf"), "w") as xf:
        xf.write_mesh(domain); xf.write_function(c_h)

    return c_h, eff


def solve_particle_tracking(domain, facet_tags, bc_ids, params, out_dir, u_h):
    from dolfinx import fem
    from dolfinx.fem import functionspace
    import dolfinx.geometry as dgeom
    import meshio
    from scipy.integrate import solve_ivp

    n_part  = int(float(params.get("Number of Particles", 100)))
    t_end   = 1.0
    bb_tree = dgeom.bb_tree(domain, domain.topology.dim)

    def vel_at(pt):
        pts  = np.array([pt], dtype=np.float64)
        cand = dgeom.compute_collisions_points(bb_tree, pts)
        hits = dgeom.compute_colliding_cells(domain, cand, pts)
        lnks = hits.links(0)
        if len(lnks) == 0:
            return np.zeros(3)
        return u_h.eval(pts, lnks[:1])[0]

    V_c   = functionspace(domain, ("Lagrange", 1))
    in_d  = fem.locate_dofs_topological(V_c, 2, facet_tags.find(bc_ids["inlet"]))
    ic    = domain.geometry.x[in_d]
    side  = max(1, int(np.sqrt(n_part)))
    xs    = np.linspace(ic[:, 0].min(), ic[:, 0].max(), side)
    ys    = np.linspace(ic[:, 1].min(), ic[:, 1].max(), side)
    z0    = float(ic[:, 2].mean())
    seeds = [[x, y, z0] for x in xs for y in ys][:n_part]

    print(f"  Tracking {len(seeds)} particles...")
    all_tracks = []
    for i, s0 in enumerate(seeds):
        try:
            sol = solve_ivp(lambda t, y: vel_at(y), [0, t_end], s0,
                            method="RK45", max_step=1e-4)
            all_tracks.append(sol.y.T)
        except Exception:
            pass
        if (i + 1) % 10 == 0:
            print(f"    {i+1}/{len(seeds)}", end="\r")
    print(f"\n  Completed {len(all_tracks)} tracks")

    pts_arr = np.vstack(all_tracks)
    segs, idx = [], 0
    for tr in all_tracks:
        n = len(tr)
        segs.extend([[idx+j, idx+j+1] for j in range(n-1)])
        idx += n

    meshio.write(os.path.join(out_dir, "particle_tracks.vtu"),
                 meshio.Mesh(points=pts_arr, cells=[("line", np.array(segs))]))
    print("  Particle tracks written.")
    return all_tracks


def solve_heat_transfer(domain, facet_tags, bc_ids, params, out_dir, u_h):
    from dolfinx import fem, io
    from dolfinx.fem import functionspace, Constant
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    from petsc4py import PETSc
    import ufl

    k_val   = float(params.get("Thermal Conductivity (W/m·K)",  0.6))
    rho_val = float(params.get("Fluid Density (kg/m³)",         1000.0))
    cp_val  = float(params.get("Specific Heat Cp (J/kg·K)",     4182.0))
    T_in    = float(params.get("Inlet Temperature (°C)",         20.0))
    T_wall  = float(params.get("Wall Temperature (°C)",          37.0))

    V_T   = functionspace(domain, ("Lagrange", 1))
    T, ps = ufl.TrialFunction(V_T), ufl.TestFunction(V_T)
    k     = Constant(domain, PETSc.ScalarType(k_val))
    rcp   = Constant(domain, PETSc.ScalarType(rho_val * cp_val))

    a = (k * ufl.inner(ufl.grad(T), ufl.grad(ps))
         + rcp * ufl.dot(u_h, ufl.grad(T)) * ps) * ufl.dx
    L = Constant(domain, PETSc.ScalarType(0.0)) * ps * ufl.dx

    in_d  = fem.locate_dofs_topological(V_T, 2, facet_tags.find(bc_ids["inlet"]))
    bcs   = [fem.dirichletbc(PETSc.ScalarType(T_in), in_d, V_T)]
    wall_f = facet_tags.find(bc_ids["wall"])
    if len(wall_f) > 0:
        wd = fem.locate_dofs_topological(V_T, 2, wall_f)
        bcs.append(fem.dirichletbc(PETSc.ScalarType(T_wall), wd, V_T))

    print("  Solving energy equation...")
    T_h = LinearProblem(a, L, bcs=bcs,
        petsc_options={"ksp_type": "gmres", "pc_type": "ilu"}).solve()
    T_h.name = "temperature"

    tv = T_h.x.array
    print(f"  Temperature range: {tv.min():.2f} – {tv.max():.2f} °C")

    with io.XDMFFile(MPI.COMM_WORLD, os.path.join(out_dir, "temperature.xdmf"), "w") as xf:
        xf.write_mesh(domain); xf.write_function(T_h)

    return T_h


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    job  = load_job(args.job)

    params      = job["parameters"]
    sim_type    = job["simulation_type"]
    solver_back = job["solver_backend"]
    stl_file    = job["stl_file"]

    if not os.path.isfile(stl_file):
        # Try same directory as job file
        stl_file = os.path.join(os.path.dirname(args.job or ""),
                                job["stl_filename"])
    if not os.path.isfile(stl_file):
        sys.exit(f"STL file not found: {job['stl_filename']}\n"
                 f"Place it in the same folder as the job .json")

    OUTPUT_DIR = os.path.join(os.path.dirname(stl_file),
                              "output", job["job_id"])
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print()
    print("=" * 55)
    print(f"  Job:        {job['job_id']}")
    print(f"  Simulation: {sim_type}")
    print(f"  STL:        {job['stl_filename']}")
    print(f"  Mesh:       {job['mesh_resolution']}")
    print(f"  Output:     {OUTPUT_DIR}")
    print("=" * 55)

    # ── Mesh ──────────────────────────────────────────────────────────────────
    print("\n[1/3] Meshing...")
    msh_path = build_mesh(stl_file, job["mesh_resolution"], OUTPUT_DIR)

    from dolfinx.io import gmshio
    from mpi4py import MPI

    domain, cell_tags, facet_tags = gmshio.read_from_msh(
        msh_path, MPI.COMM_WORLD, gdim=3
    )
    domain.topology.create_connectivity(domain.topology.dim - 1, domain.topology.dim)

    n_cells = domain.topology.index_map(domain.topology.dim).size_global
    n_verts = domain.topology.index_map(0).size_global
    print(f"  Mesh: {n_cells} cells, {n_verts} vertices")

    # ── Solve ─────────────────────────────────────────────────────────────────
    print(f"\n[2/3] Solving: {sim_type}...")
    results = {"job_id": job["job_id"], "simulation_type": sim_type}

    if solver_back == "stokes":
        u_h, p_h, dP = solve_stokes(domain, facet_tags, BC_IDS, params, OUTPUT_DIR)
        results.update({"pressure_drop_Pa": dP,
                        "files": ["velocity.xdmf", "pressure.xdmf"]})

    elif solver_back == "mixing":
        u_h, p_h, dP = solve_stokes(domain, facet_tags, BC_IDS, params, OUTPUT_DIR)
        c_h, eff = solve_mixing(domain, facet_tags, BC_IDS, params, OUTPUT_DIR, u_h)
        results.update({"pressure_drop_Pa": dP,
                        "mixing_efficiency_pct": round(eff * 100, 2),
                        "files": ["velocity.xdmf", "pressure.xdmf", "concentration.xdmf"]})

    elif solver_back == "particle_tracking":
        u_h, p_h, dP = solve_stokes(domain, facet_tags, BC_IDS, params, OUTPUT_DIR)
        tracks = solve_particle_tracking(domain, facet_tags, BC_IDS, params, OUTPUT_DIR, u_h)
        results.update({"pressure_drop_Pa": dP,
                        "tracked_particles": len(tracks),
                        "files": ["velocity.xdmf", "particle_tracks.vtu"]})

    elif solver_back == "heat_transfer":
        u_h, p_h, dP = solve_stokes(domain, facet_tags, BC_IDS, params, OUTPUT_DIR)
        T_h = solve_heat_transfer(domain, facet_tags, BC_IDS, params, OUTPUT_DIR, u_h)
        results.update({"pressure_drop_Pa": dP,
                        "files": ["velocity.xdmf", "pressure.xdmf", "temperature.xdmf"]})

    elif solver_back == "dean_flow":
        print("  Dean flow: Stokes linear approximation.")
        u_h, p_h, dP = solve_stokes(domain, facet_tags, BC_IDS, params, OUTPUT_DIR)
        results.update({"pressure_drop_Pa": dP,
                        "files": ["velocity.xdmf", "pressure.xdmf"]})

    else:
        sys.exit(f"Unknown solver backend: {solver_back}")

    # ── Save summary ──────────────────────────────────────────────────────────
    results["completed_at"] = datetime.datetime.now().isoformat()
    summary_path = os.path.join(OUTPUT_DIR, "results_summary.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n[3/3] Results saved to: {OUTPUT_DIR}")
    print()
    print("=" * 55)
    for k, v in results.items():
        if k != "files":
            print(f"  {k}: {v}")
    print("  Output files:")
    for fn in results.get("files", []):
        print(f"    {OUTPUT_DIR}/{fn}")
    print()
    print("Open in ParaView:")
    print("  File → Open → select .xdmf files from the output folder above")
    print("=" * 55)


if __name__ == "__main__":
    main()
