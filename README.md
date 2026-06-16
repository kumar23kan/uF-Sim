# Microfluidic CFD Simulation

A two-part toolchain for running computational fluid dynamics (CFD) simulations on microfluidic channel geometries exported from Fusion 360.

---

## Overview

| Component | File | Purpose |
|---|---|---|
| Desktop UI | `uf_sim_ui.py` | Set parameters, choose simulation, export job config |
| Local runner | `run_simulation.py` | Run simulation locally using the `uf-sim` conda environment |
| Colab notebook | `microfluidic_cfd.ipynb` | Mesh, solve, and export results for ParaView (cloud option) |

The recommended workflow (local):
1. Export your microfluidic channel geometry from Fusion 360 as **STL**
2. Open the **desktop UI** to configure your simulation and export a job `.json`
3. Run `run_simulation.py` locally — it reads the `.json`, meshes, solves, and writes output
4. Open the results in **ParaView**

---

## Requirements

### Desktop UI only

- Python 3.9 or later with `tkinter` (included with most Python installs)

```bash
python3 uf_sim_ui.py
```

### Local simulation runner

Create the `uf-sim` conda environment once:

```bash
conda create -n uf-sim -c conda-forge \
    fenics-dolfinx mpich gmsh meshio h5py scipy python=3.11 -y
```

Then activate it before running simulations:

```bash
conda activate uf-sim
python3 run_simulation.py <job.json>
```

### Google Colab (cloud option)

Dependencies are installed automatically by Cell 2 of the notebook:

- [FEniCSx](https://fenicsproject.org/) — finite element solver
- [Gmsh](https://gmsh.info/) — 3D mesh generation from STL
- `meshio`, `h5py`, `scipy` — mesh I/O and particle integration

---

## Meshing Strategy

The notebook uses a two-layer mesh strategy inside `build_mesh()` (Cell 6):

**Bulk mesh**
- Size is clamped to 25% of the channel's shortest cross-section dimension, so thin channels always get enough cells regardless of the overall bounding box size.

**Boundary layer (auto-calculated)**
- Applied automatically to all wall surfaces.
- First-layer thickness = a fraction of the channel depth, scaled by resolution:

| Resolution | First-layer fraction | Layers | Growth ratio |
|---|---|---|---|
| Coarse | 10% of depth | 5 | 1.3× |
| Medium | 8% of depth | 5 | 1.3× |
| Fine | 6% of depth | 5 | 1.3× |
| Very Fine | 5% of depth | 5 | 1.3× |

This ensures the near-wall velocity gradient is resolved accurately for wall shear stress, pressure drop, and species concentration at walls — all critical in microfluidic channels.

---

## Step-by-step Instructions

### 1. Export STL from Fusion 360

1. Open your `.f3d` file in Fusion 360
2. Right-click the body/component in the browser → **Save As Mesh**
3. Format: **STL**, Units: **mm**, Refinement: High
4. Save the `.stl` file

---

### 2. Configure the simulation (Desktop UI)

```bash
python3 uf_sim_ui.py
```

| Section | What to do |
|---|---|
| **Geometry File** | Browse and select your `.stl` file |
| **Simulation Type** | Choose from the dropdown |
| **Parameters** | Edit the pre-filled values |
| **Mesh Resolution** | Coarse / Medium / Fine / Very Fine |
| **Output Format** | VTU or XDMF (both open in ParaView) |
| **Export Job File** | Click to save a `.json` config file |

The exported `.json` bundles all settings and is passed to Colab.

#### Simulation types

| Type | Physics solved |
|---|---|
| **Pressure Drop** | Stokes flow — velocity + pressure field, pressure drop value |
| **Flow Mixing** | Advection-diffusion — species concentration, mixing efficiency |
| **Particle Tracking** | Lagrangian — particle trajectories through the flow field |
| **Heat Transfer** | Energy equation — temperature distribution |
| **Dean Flow** | Stokes approximation for curved/spiral channels |

---

### 3a. Run the simulation locally (recommended)

```bash
conda activate uf-sim
python3 run_simulation.py uFSim_<timestamp>.json
```

- If no argument is given, a file-picker dialog opens.
- Results are written to `output/<job_id>/` next to the STL file.
- A `results_summary.json` is also saved with pressure drop, mixing efficiency, etc.

Typical runtime on a modern laptop:

| Resolution | Approximate time |
|---|---|
| Coarse | 30–90 seconds |
| Medium | 2–5 minutes |
| Fine | 10–20 minutes |
| Very Fine | 30–60 minutes |

---

### 3b. Run the simulation on Google Colab (cloud option)

1. Open [Google Colab](https://colab.research.google.com)
2. Upload `microfluidic_cfd.ipynb` via **File → Upload notebook**
3. Run the cells in order:

| Cell | Action | Notes |
|---|---|---|
| **Cell 2** | Install FEniCSx | Run once per session — takes ~5-10 min |
| **Cell 3** | Imports | Confirm versions print without error |
| **Cell 4** | Mount Drive + upload job | Upload your `.json` when prompted |
| **Cell 5** | Load config + upload STL | Upload your `.stl` when prompted |
| **Cell 6** | Generate mesh | Gmsh builds tetrahedral mesh + boundary layer from STL |
| **Cell 7-10** | Solver definitions | Functions are defined, not run yet |
| **Cell 11** | Run simulation | Automatically calls the right solver |
| **Cell 12** | Download results | Downloads a `.zip` of all output files |

> **Tip:** Enable GPU runtime in Colab for faster mesh generation:  
> Runtime → Change runtime type → T4 GPU

---

### 4. Visualise in ParaView

1. Download and install [ParaView](https://www.paraview.org/download/)
2. Extract the downloaded `.zip`
3. Open ParaView → **File → Open**
4. Select the output file for your simulation:

| Simulation | File to open |
|---|---|
| Pressure Drop | `pressure.xdmf`, `velocity.xdmf` |
| Flow Mixing | `concentration.xdmf` |
| Particle Tracking | `particle_tracks.vtu` |
| Heat Transfer | `temperature.xdmf` |

5. Click **Apply** in the Properties panel
6. Use the dropdown at the top toolbar to colour by field:
   - `pressure`, `velocity`, `concentration`, `temperature`
7. For particle tracks: add a **Tube** filter to give the lines width

---

## Boundary Conditions and Mesh

The mesh builder auto-detects boundaries from the STL bounding box:

| Boundary | Condition |
|---|---|
| **Inlet** | Face at the minimum extent of the longest axis — uniform inlet velocity |
| **Outlet** | Face at the maximum extent — zero pressure (outflow) |
| **Walls** | All remaining surfaces — no-slip |

If auto-detection fails for your geometry (e.g. channels not axis-aligned), edit the `build_mesh()` function in Cell 6 of the notebook and manually assign physical group tags.

---

## File Structure

```
Research/uF_Sim/
├── uf_sim_ui.py            # Desktop parameter UI
├── run_simulation.py       # Local simulation runner
├── microfluidic_cfd.ipynb  # Google Colab simulation notebook
├── generate_notebook.py    # Script that regenerates the notebook
└── README.md               # This file
```

---

## Troubleshooting

**FEniCSx installation fails in Colab**  
FEM-on-Colab may be temporarily unavailable. The notebook falls back to `condacolab` automatically — this will restart the runtime. Re-run from Cell 3 after the restart.

**Mesh generation fails with "no volume"**  
The STL may have open surfaces or self-intersections. In Fusion 360:
- Use **Inspect → Check** to find issues
- Re-export with **Watertight** body only

**Solver diverges or gives zero velocity**  
Check the inlet direction. The notebook assumes flow along the longest axis (+Z by default). If your channel is oriented differently, adjust the inlet velocity vector in `solve_stokes()` Cell 7:
```python
u_in = Constant(domain, PETSc.ScalarType((inlet_vel, 0.0, 0.0)))  # X-axis flow
```

**ParaView shows empty scene**  
Open the `.xdmf` file (not the `.h5`). Both files must be in the same folder.
