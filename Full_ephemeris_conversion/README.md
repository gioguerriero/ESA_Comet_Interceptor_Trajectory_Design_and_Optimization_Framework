# `Full_ephemeris_conversion/` — Full-Ephemeris Refinement (Moon Gravity + TCMs)

## Table of Contents

- [Overview](#overview)
- [Optimization Pipeline](#optimization-pipeline)
- [Design Variables](#design-variables)
- [Constraints](#constraints)
- [Cost Function](#cost-function)
- [Folder Structure](#folder-structure)

---

## Overview

This module is the **final fidelity upgrade** of the workflow. It takes the Sun+Earth ephemeris solution produced by [`../CR3BP 2 Ephemeris Sun-Earth/`](../CR3BP%202%20Ephemeris%20Sun-Earth/) and re-optimises it in a **fully realistic model in which the Moon's gravity is explicitly integrated** (`NBODY_J2000_full_ephe` in [`../Auxiliar/`](../Auxiliar/)). Adding the lunar attraction perturbs the previously ballistic flyby, so the trajectory would no longer close.

To absorb this perturbation, two **trajectory-correction manoeuvres (TCMs)** are introduced, a few days **before** the flyby (**TCM1**) and a few days **after** it (**TCM2**). The flyby itself is now targeted geometrically through its **B-plane** rather than by an instantaneous `v_inf` rotation. Because the perturbation splits naturally at the flyby, the refinement is decomposed into two independent sub-problems — **pre-flyby** and **post-flyby** — each solved by its own local optimiser.

> As noted in the project overview, this stage is not fully robust across all solutions; it demonstrates that refined, real-ephemeris trajectories with Moon gravity exist and preserve the CR3BP ΔV/ToF, but a general-purpose pipeline is left to future work.

---

## Optimization Pipeline

Driven by [`run_refinement`](../Opt%20Manager/run_refinement.m), the stage runs as follows.

### 1. Initialization — B-plane targeting of the flyby

Starting from the Sun+Earth solution (`out_cr3bp`), the pre-flyby state is back-propagated to the **TCM1** epoch (a few days before the flyby). `bplane_from_vinf` computes the desired **B-plane target point** (`B_T`, `B_R`) from the CR3BP `v_inf` vectors, and `bplane_tcm` solves for the **TCM1 impulse** that drives the perturbed flyby onto that B-plane target — a Newton iteration with a numerical Jacobian that repeatedly propagates the state to the flyby in the full model.

### 2. Post-flyby refinement

The segment from the flyby to the comet (TCM2 → DSM2 → comet) is refined. Two interchangeable formulations are provided:

- **Fast path** — `Refinement post flyby/refinement_post_flyby`: a direct `fmincon` on TCM2 + DSM2 + DSM2 epoch, targeting the comet position.
- **Global path** — `Refinement global post flyby/refinement_global_post_flyby`: a more general formulation (B-plane targets + control point + comet arrival velocity), designed for a Monotonic-Basin-Hopping-style search, enforcing the flyby altitude explicitly.

### 3. Pre-flyby refinement

`Refinement pre flyby/refinement_pre_flyby` refines the segment from injection to TCM1 (injection + DSM1 + DSM1 epoch), enforcing a **smooth patch** with the post-flyby leg.

### 4. Post-processing

`variables_organizer_refined` merges the pre- and post-flyby solutions (together with the TCM1 impulse) into the canonical **refined `out` struct**, containing every arc, epoch, and manoeuvre (`departure`, `injection`, `dsm1`, `tcm1`, `flyby`, `tcm2`, `dsm2`, `comet`) in synodic and J2000 form. `plot_full_trajectory` and `plot_traj_post_flyby` visualise the result.

### 5. Outputs

The stage returns the fully refined `out` trajectory struct plus the convergence flags of the two sub-optimisations (`success_pre`, `success_post`), used by `main.m` to decide whether to plot/accept the refined point.

---

## Design Variables

The refinement uses **three small design vectors**, one per sub-problem. TCM1 is *not* a design variable of any `fmincon`: it is computed as a **function of the B-plane targets** by `bplane_tcm`.

### Pre-flyby — `refinement_pre_flyby` (7 variables)

| Var | Symbol | Physical meaning | Role |
|---|---|---|---|
| `x(1:3)` | `ΔV_inj` | Injection impulse [km/s] | Departure from the halo |
| `x(4:6)` | `ΔV_DSM1` | First deep-space manoeuvre [km/s] | Deflection toward the Moon |
| `x(7)` | `Δepoch_DSM1` | DSM1 epoch variation [days] | Fine-tunes the DSM1 timing |

### Post-flyby, fast — `refinement_post_flyby` (7 variables)

| Var | Symbol | Physical meaning | Role |
|---|---|---|---|
| `x(1:3)` | `ΔV_TCM2` | Post-flyby cleanup impulse [km/s] | Absorbs residual Moon-gravity perturbation |
| `x(4:6)` | `ΔV_DSM2` | Second deep-space manoeuvre [km/s] | Final phasing to the comet |
| `x(7)` | `epoch_DSM2 / 1e8` | DSM2 epoch (scaled) | Times DSM2 |

### Post-flyby, global — `refinement_global_post_flyby` (12 variables)

| Var | Symbol | Physical meaning | Role |
|---|---|---|---|
| `x(1)`, `x(2)` | `B_T / 5000`, `B_R / 5000` | B-plane target components (scaled) | Target the flyby geometry (→ TCM1 via `bplane_tcm`) |
| `x(3:5)` | `ΔV_DSM2 / 0.1` | DSM2 impulse (scaled) | Final phasing |
| `x(6)` | `Δepoch_DSM2` | DSM2 epoch variation [days] | Times DSM2 |
| `x(7:9)` | `v_comet / 10` | Comet arrival velocity (scaled) | Control-point end velocity for the backward leg |
| `x(10:12)` | `ΔV_TCM2 / 0.1` | Post-flyby cleanup impulse (scaled) | Absorbs the perturbation |

All impulses are scaled to order 1 for good gradient conditioning.

---

## Constraints

Each sub-problem has its own nonlinear-constraint function returning `[c, ceq]`.

### `NC_refinement_pre_flyby`

| Constraint | Type | Interpretation | Implementation |
|---|---|---|---|
| **Maneuver spacing** | `c` | DSM1 kept clear of injection and TCM1 | timing inequalities on `epoch_dsm1` |
| **Smooth patch** | `ceq` | The injection-forward arc and the TCM1-backward arc must meet in position and velocity | `S_pre_patch − S_post_patch`, scaled (~100 km / 10 m/s) |

### `NC_refinement_post_flyby` (fast path)

Propagates TCM1 → TCM2 → DSM2 → comet in the full model and enforces:

| Constraint | Type | Interpretation | Implementation |
|---|---|---|---|
| **Maneuver spacing** | `c` | DSM2 kept clear of the other manoeuvres | timing inequalities on `epoch_dsm` |
| **Comet rendezvous** | `ceq` | Trajectory must reach the comet position at its epoch | `(S_comet(1:3) − comet_pos)/1e6` |

> A root-level `NC_refinement_post_flyby.m` provides a documented, alternative form of the same constraints; the copy inside `Refinement post flyby/` (which also carries the timing inequalities) is the one used by the fast-path driver.

### `NC_refinement_global_post_flyby`

| Constraint | Type | Interpretation | Implementation |
|---|---|---|---|
| **Flyby altitude** | `c` | The B-plane-targeted flyby must stay above the minimum lunar altitude | `(−h_flyby + h_min_moon)/0.5 ≤ 0` (via `bplane_tcm`) |
| **Maneuver spacing** | `c` | DSM2 kept clear of the flyby/comet | timing inequalities |
| **Smooth patch** | `ceq` | The forward TCM2 arc and the backward comet arc must meet in position/velocity at a patch point | `(S_patch_post − S_patch_pre)`, scaled (~1 km / 0.1 m/s) |

The **smooth-patch equalities are the closure conditions** that make the perturbed trajectory continuous once the Moon's gravity is present; the **flyby-altitude inequality** is what makes the geometric B-plane target physically realisable.

---

## Cost Function

Each sub-problem minimises the **total magnitude of the manoeuvres it controls**, keeping the time of flight fixed:

| Sub-problem | Objective | Meaning |
|---|---|---|
| Pre-flyby | `OF_refinement_pre_flyby` = `‖ΔV_inj‖ + ‖ΔV_DSM1‖` | Minimise the pre-flyby propellant. |
| Post-flyby (fast) | `OF_refinement_post_flyby` = `10⁻²·(‖ΔV_TCM2‖ + ‖ΔV_DSM2‖)` | Minimise the post-flyby propellant (scaled). |
| Post-flyby (global) | `OF_refinement_global_post_flyby` = `‖ΔV_TCM1‖ + ‖ΔV_TCM2‖ + ‖ΔV_DSM2‖` | Minimise all controlled impulses, including the B-plane-derived TCM1. |

The design intent is to **close the trajectory in the full-ephemeris model at minimum additional ΔV**: the TCMs exist only to compensate the lunar perturbation, so keeping them (and the DSMs) small is exactly the objective. Time of flight is not penalised because it is inherited from the Sun+Earth solution being refined.

---

## Folder Structure

| File / sub-folder | Role |
|---|---|
| `bplane_from_vinf.m` | Computes B-plane parameters (`B_T`, `B_R`, `r_p`, axes) from incoming/outgoing `v_inf` (Vallado Algorithm 79). |
| `bplane_tcm.m` | Solves for the TCM that drives the flyby onto a target B-plane point (Newton iteration with numerical Jacobian, full-ephemeris propagation). |
| `variables_organizer_refined.m` | Merges pre-/post-flyby solutions + TCM1 into the canonical refined `out` struct. |
| `plot_full_trajectory.m` | Plots the complete refined mission trajectory in the synodic frame. |
| `plot_traj_post_flyby.m` | Plots the guess vs. solution of the post-flyby refinement. |
| `NC_refinement_post_flyby.m` | Root-level (documented) form of the post-flyby constraints. |
| `Refinement pre flyby/` | Pre-flyby sub-problem: `refinement_pre_flyby`, `OF_refinement_pre_flyby`, `NC_refinement_pre_flyby` (+ a guess-plot helper). |
| `Refinement post flyby/` | Fast post-flyby sub-problem: `refinement_post_flyby`, `OF_refinement_post_flyby`, `NC_refinement_post_flyby`. |
| `Refinement global post flyby/` | Global post-flyby sub-problem: `refinement_global_post_flyby`, `OF_refinement_global_post_flyby`, `NC_refinement_global_post_flyby`. |

## Related folders

- [`../CR3BP 2 Ephemeris Sun-Earth/`](../CR3BP%202%20Ephemeris%20Sun-Earth/) — the previous fidelity step, whose `out_cr3bp` is refined here.
- [`../Opt Manager/`](../Opt%20Manager/) — `run_refinement`, which drives both refinement stages.
- [`../Auxiliar/`](../Auxiliar/) — `NBODY_J2000_full_ephe`, `synodic2sun_J2000`/`sun_J2000_to_synodic`, `flyby_periapsis_state`.
