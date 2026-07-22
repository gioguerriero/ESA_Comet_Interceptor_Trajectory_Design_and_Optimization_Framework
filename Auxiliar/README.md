# `Auxiliar/` — Core Utilities & Shared Physics

## Table of Contents

- [Overview](#overview)
- [Dynamics — Equations of Motion](#dynamics--equations-of-motion)
- [Reference-Frame Conversions](#reference-frame-conversions)
- [Flyby & Manoeuvre Modelling](#flyby--manoeuvre-modelling)
- [Orbital-Mechanics Libraries](#orbital-mechanics-libraries)
- [Trajectory Assembly & Export](#trajectory-assembly--export)
- [Visualisation](#visualisation)
- [`Checks and tests/` — Validation & Analysis](#checks-and-tests--validation--analysis)

---

## Overview

`Auxiliar/` is the **core utility layer**: it implements the shared physics and low-level tools that every other module depends on. Whereas the optimization folders each own a phase of the pipeline, this folder holds the pieces they all reuse — the equations of motion, the reference-frame conversions, the lunar-flyby model, and standard orbital-mechanics libraries — so that the dynamical model is defined **once** and stays consistent across the whole framework.

Functions are grouped below by role. The `Checks and tests/` sub-folder contains standalone diagnostic scripts (not part of the automated pipeline) documented in its own section.

---

## Dynamics — Equations of Motion

These functions are the right-hand sides passed to MATLAB's `ode45`. They define the dynamical models on which every arc is propagated.

### `CR3BP`
- **Purpose:** Equations of motion of the Circular Restricted Three-Body Problem in the rotating (synodic), non-dimensional frame — the backbone dynamics of the whole search.
- **Inputs:** `t` (time, unused — autonomous system), `s` (state `[x y z vx vy vz]`), `mu` (system mass parameter).
- **Outputs:** `ds` (state derivative).
- **Usage:** Propagation of every CR3BP arc in `Moon2Comet/`, `Halo2Moon/`, and `Opt Manager/`.

### `CR3BP_STM`
- **Purpose:** CR3BP equations augmented with the **State Transition Matrix** (6 states + a 6×6 STM, 42 total), for sensitivity/differential-correction work.
- **Inputs:** `t`, `s` (augmented `[state; flattened STM]`), `mu`, `n` (mean motion, always 1; unused).
- **Outputs:** `ds` (augmented derivative).
- **Usage:** Available for STM-based analyses (e.g. halo generation / manifold work).

### `NBODY_J2000`
- **Purpose:** Spacecraft acceleration under **Sun + Earth** point-mass gravity using real SPICE ephemerides (Sun-centred ECLIPJ2000). The dynamics of the first refinement step.
- **Inputs:** `t` (time since `initial_epoch`), `S` (Sun-centred state), `initial_epoch` (ET), `c` (constants: `G`, `mEarth`, `mSun`).
- **Outputs:** `dS`.
- **Usage:** Multiple-shooting propagation in [`../CR3BP 2 Ephemeris Sun-Earth/`](../CR3BP%202%20Ephemeris%20Sun-Earth/).
- **Dependencies:** SPICE (`cspice_spkezr`).

### `NBODY_J2000_full_ephe`
- **Purpose:** Same as above but with **Sun + Earth + Moon** gravity — the full-ephemeris dynamics used to validate the flyby with explicit lunar attraction.
- **Inputs:** `t`, `S`, `initial_epoch`, `c` (adds `mMoon`).
- **Outputs:** `dS`.
- **Usage:** Propagation in [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/).
- **Dependencies:** SPICE.

### `NBODY_J2000_mod`
- **Purpose:** Homotopy dynamical model that blends CR3BP-like (Sun+Earth) and full (Sun+Earth+Moon) dynamics through a weight `eps` (`eps=0` → no Moon, `eps=1` → full).
- **Inputs:** `t`, `S`, `initial_epoch`, `eps` (homotopy weight).
- **Outputs:** `dS`.
- **Usage:** Continuation between fidelity levels.
- **Dependencies:** SPICE.

### `twoBody`
- **Purpose:** Keplerian two-body equations of motion in an inertial frame.
- **Inputs:** `t`, `S`, `muReal` (central-body GM).
- **Outputs:** `dsdt`.
- **Usage:** Local two-body propagation, e.g. the flyby-hyperbola reconstruction in `verify_flyby_geometry`.

---

## Reference-Frame Conversions

The search runs in the synodic frame; the ephemeris refinement runs in inertial J2000. These functions convert states between the two, using the *instantaneous* Sun–Earth geometry at each epoch (so length/velocity scales are epoch-dependent).

### `car2synodic` / `synodic2car`
- **Purpose:** Convert between an inertial Cartesian state and the CR3BP synodic frame (`car2synodic`: inertial → synodic; `synodic2car`: synodic → inertial). Primaries Sun/Earth.
- **Inputs:** the state, `time` (non-dimensional rotation angle), `mu`, `theta_0` (initial angle, default 0).
- **Outputs:** the converted state.
- **Usage:** Frame changes in `Moon2Comet/get_DSM_info` and related routines.

### `sun_J2000_to_synodic` / `synodic2sun_J2000`
- **Purpose:** Exact, epoch-dependent conversion between Sun-centred ECLIPJ2000 states and the non-dimensional Sun–Earth synodic frame, including the length-rate (`Ldot`) correction. They are exact inverses of one another.
- **Inputs:** the state(s) (`N×6`), `epochs` (ET), `mu`.
- **Outputs:** the converted state(s); `synodic2sun_J2000` also returns a `scales` struct (per-epoch `L`, `V`, `omega`, rotation matrices) reused for delta-v conversions.
- **Usage:** The bridge between the CR3BP search and the ephemeris refinement; used pervasively in `run_refinement`, the `variables_organizer` functions, and the diagnostic tools.
- **Dependencies:** SPICE (batch `cspice_spkezr`).

---

## Flyby & Manoeuvre Modelling

### `vinf_rotation`
- **Purpose:** The **lunar-flyby model**. Builds a `v_inf` vector by rotating the Moon's synodic velocity direction in-plane (`fpa`) and out-of-plane, then scaling by the magnitude — the instantaneous, zero-SOI representation of the gravity assist.
- **Inputs:** `synodic_moon` (Moon synodic state), `Vinf` (magnitude), `fpa` (in-plane angle), `out_of_plane` (angle).
- **Outputs:** `vinf_rotated` (3×1 vector).
- **Usage:** Central to both optimization phases and to the diagnostics — every flyby is realised through this rotation.

### `flyby_periapsis_state`
- **Purpose:** Reconstruct the flyby-hyperbola **periapsis state** (non-dimensional synodic) from the patched-conics `v_inf` vectors — a first approximation of the "true" flyby used when exporting/validating.
- **Inputs:** `vinf_in`, `vinf_out` (Moon-relative), `moon_state`, `epoch`, `muMoon`, `mu`, `c_const`.
- **Outputs:** `S_syn` (periapsis state).
- **Usage:** `write_python_inputs`, `verify_flyby_geometry`.
- **Dependencies:** `bplane_from_vinf` (in `Full_ephemeris_conversion/`), `sun_J2000_to_synodic`.

### `direction`
- **Purpose:** Build a unit delta-v direction in the local velocity (TNH) frame from an in-plane and an out-of-plane steering angle.
- **Inputs:** `state` (`[r; v]`), `alpha` (in-plane), `phi` (out-of-plane).
- **Outputs:** `d_hat` (unit vector).
- **Usage:** Steering-direction construction for manoeuvre parametrisation.

---

## Orbital-Mechanics Libraries

> These three are standard, well-established library functions (third-party origin) with their own extensive headers; they are used as black boxes.

### `lambert`
- **Purpose:** Robust Lambert solver (Izzo + Lancaster/Blanchard/Gooding) — solves the two-point boundary-value problem for a ballistic arc.
- **Inputs:** `r1`, `r2` (position vectors), `tf` (time of flight), `m` (revolutions), `GM_central`.
- **Outputs:** `V1`, `V2` (terminal velocities), extremal distances, exit flag.
- **Usage:** The heliocentric DSM2 leg in `Moon2Comet/get_DSM_info` and the GA objective.

### `astroConstants` / `getAstroConstants`
- **Purpose:** Return standard astrodynamic/planetary constants (`astroConstants` by numeric identifier; `getAstroConstants` a string-based wrapper).
- **Inputs:** an identifier vector / name strings.
- **Outputs:** the requested constant value(s).
- **Usage:** Constant lookups in `main.m` and elsewhere (e.g. Sun/Earth masses).

---

## Trajectory Assembly & Export

### `extract_mission_points`
- **Purpose:** Extract the key mission states (synodic, non-dimensional) and the segment times of flight from an optimised ephemeris solution — pre-injection, injection, DSM1, pre/post-flyby, DSM2, and the total TOF.
- **Inputs:** `x_opt` (optimised design vector), `S_halo_syn`, `target_position`, `epoch_comet_flyby`, `c`.
- **Outputs:** `out` (struct of key states + TOFs + arrival date).
- **Usage:** Post-analysis / verification of a refined solution (run manually).
- **Dependencies:** `CR3BP`, `NBODY_J2000`, `synodic2sun_J2000`.

### `write_python_inputs` / `write_python_inputs_serot`
- **Purpose:** Export a refined trajectory to a Python-readable text file — states, durations, and the encounter date. `write_python_inputs` writes the synodic non-dimensional frame; `write_python_inputs_serot` writes the dimensional Earth-centred SEROT frame.
- **Inputs:** `out_cr3bp` (a `variables_organizer` output), `c_const`, `filename` (optional).
- **Outputs:** none (writes a `.txt` file).
- **Usage:** Hand-off to external Python/Blender tooling.
- **Dependencies:** `flyby_periapsis_state`, `sun_J2000_to_synodic` / `synodic2sun_J2000`, SPICE time conversion.

### `propagate_S_pre_inj` *(script)*
- **Purpose:** Standalone script that propagates a manually-supplied pre-injection state with `NBODY_J2000`, converts back to the synodic frame, and plots it against the halo orbit. A verification/inspection aid.
- **Inputs:** `S_pre_inj_manual` in the workspace (edit parameters at the top).
- **Outputs:** a figure and printed states.

---

## Visualisation

### `plot_comets_synodic`
- **Purpose:** Plot the encounter positions of a set of comets in the non-dimensional synodic frame, optionally overlaying the unstable-manifold arcs.
- **Inputs:** `comets` (cell of comet structs), `mu`, plus optional Name-Value pairs (`PlotManifolds`, `S_halo`, `UnstableDir`, `C`, ...).
- **Outputs:** the comet positions (and a figure).
- **Usage:** Called in `main.m` to visualise target geometry; `Dependencies:` `CR3BP`, SPICE.

### `plot_trajectory_heliocentric`
- **Purpose:** Plot the full mission trajectory in the Sun-centred ECLIPJ2000 frame, colour-coded by phase, with event markers and body orbits.
- **Inputs:** trajectory states/times, event epochs and states, `selected_comet`.
- **Outputs:** a 3D figure.
- **Usage:** CR3BP-based trajectory visualisation inside `run_refinement`.
- **Dependencies:** SPICE.

### `plot_paper_trajectory`
- **Purpose:** Produce the publication-quality dual-panel mission figure (heliocentric overview + halo-departure zoom).
- **Inputs:** `global_results`, solution index `k`, `S_halo`, `c`, `selected_comet`, `unstable_dir`, `refine_params`.
- **Outputs:** a figure.
- **Usage:** Figure generation for reporting (run manually).
- **Dependencies:** `run_refinement`.

---

## `Checks and tests/` — Validation & Analysis

Standalone tools used to validate the pipeline and analyse a saved result set. They are **not** part of the automated run; they are executed manually, typically against a `Results/<comet>_runN/` dataset.

### `verify_flyby_geometry`
- **Purpose:** Reconstruct and visualise the lunar-flyby geometry (hyperbola, periapsis, B-plane, `v_inf` vectors) in non-dimensional Moon-relative coordinates, and numerically check that the asymptotic velocities reproduce `v_inf,in`/`v_inf,out`.
- **Inputs:** `out_cr3bp` (a `variables_organizer` output), `c_const`.
- **Outputs:** `info` (reconstructed geometry + verification residuals) and a figure.
- **Dependencies:** `bplane_from_vinf`, `twoBody`.

### `flyby_bending_check`
- **Purpose:** For each search solution, compare the **actual** flyby bending angle against the **maximum feasible** bending at the minimum allowed periapsis, flagging any physically unrealisable flyby.
- **Inputs:** `global_results`, `c`, optional Name-Value (`HminKm`, `Plot`, `Verbose`).
- **Outputs:** `res` (per-solution bending, differences, feasibility counts).
- **Dependencies:** `moon_state`, `vinf_rotation`.

### `earth_altitude_check`
- **Purpose:** For each search solution, reconstruct the full CR3BP trajectory and compute the **minimum Earth flyby altitude** (with a coarse-then-refined closest-approach search), to screen against passes too close to Earth.
- **Inputs:** `global_results`, `c`, `S_halo`, `unstable_dir`, `eps_vel_ms`, optional Name-Value.
- **Outputs:** `res` (per-solution min altitude, threshold counts, max heliocentric distance).
- **Dependencies:** `CR3BP`, `state_finder`, `vinf_rotation`, `moon_state`.

### `vinf_escape_study`
- **Purpose:** Energy-based validation of the swingby benefit. Computes the **ballistic-equivalent escape velocity** at a control radius by removing the propulsive energy of the manoeuvres, for comparison against reference literature values.
- **Inputs:** `global_results`, `c`, `S_halo`, `unstable_dir`, `eps_vel_ms`, `control_radius`, `selected_comet`, optional Name-Value.
- **Outputs:** `res` (signed/absolute `v_inf,ball`, statistics, anomaly report) and a PDF figure.
- **Dependencies:** `CR3BP`, `state_finder`, `vinf_rotation`, `moon_state`, `lambert`.

### `compare_vinf_sampling`
- **Purpose:** Visualise/compare the rectangular (meshgrid) vs. polar (uniform-areal-density) sampling of `v_inf` directions on the deflection cone — the sampling scheme used by `build_moon_manifold_db`.
- **Inputs / Outputs:** self-contained demonstration (figures).

### `rotation_test` *(script)*
- **Purpose:** Sanity-check script for `vinf_rotation` under the current sign conventions, using a dummy synodic Moon state.
- **Inputs / Outputs:** self-contained (figures).

---

## Related Folders

- [`../Opt Manager/`](../Opt%20Manager/), [`../Moon2Comet/`](../Moon2Comet/), [`../Halo2Moon/`](../Halo2Moon/) — consumers of the dynamics and flyby model.
- [`../CR3BP 2 Ephemeris Sun-Earth/`](../CR3BP%202%20Ephemeris%20Sun-Earth/), [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/) — consumers of the N-body dynamics and frame conversions.
- [`../kernels/`](../kernels/) — the SPICE kernels required by every `cspice_*` call used here.
