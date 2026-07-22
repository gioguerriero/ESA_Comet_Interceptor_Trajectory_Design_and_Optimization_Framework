# `GODOT Refinement/` — Full-Ephemeris Refinement (Python / ESA GODOT)

## Purpose

This folder contains a **Python** implementation of a full-ephemeris trajectory refinement built on top of **ESA's [GODOT](https://godot.io.esa.int/) library** (the Agency's operational-grade orbit/trajectory propagation and optimisation toolkit). Its role is to take a trajectory produced by the MATLAB optimisation pipeline in this repository and **refine it in a high-fidelity, real-ephemeris dynamical model**, beyond the fidelity reached by the MATLAB-side ephemeris refinement (`Full_ephemeris_conversion/`).

Conceptually, this is the **final rung of the fidelity ladder**: the MATLAB framework performs the trajectory *search* and the CR3BP/intermediate-ephemeris optimisation; GODOT then re-propagates and refines the selected solution with full, operational-quality dynamics.

## Interface with the MATLAB pipeline

The two environments communicate through a single, human-readable **text file** exported on the MATLAB side by [`write_python_inputs`](../Auxiliar/write_python_inputs.m) (and its SEROT variant `write_python_inputs_serot`). That file carries everything needed to reconstruct the trajectory as an initial guess:

- the key mission states (post-injection, pre-injection, post-DSM1, flyby periapsis, post-DSM2);
- the durations of each segment;
- the comet encounter date.

The Python code in this folder reads that file and uses it to set up the GODOT refinement problem. See [§8 of the project README](../README.md#8-full-ephemeris-refinement-bridge-godot) for the description of the bridge from the MATLAB side.

```
MATLAB pipeline ──(write_python_inputs → .txt)──▶ GODOT Refinement (Python)
```

## Status — work in progress

> **Important:** this code is **not currently functional**. In its present state it **raises errors and does not complete a refinement**.

It is included in the repository deliberately, because:

- the **overall structure is correct** — the problem set-up, the import of the MATLAB interface file, and the intended GODOT refinement workflow are in place and follow a sound design;
- it was **left unfinished for lack of time**, not because the approach is wrong.

As such, it is meant to be **picked up and completed** by a future developer rather than treated as a finished tool. A reasonable starting point is to get the import/propagation set-up running against a known-good exported trajectory, then progressively enable the refinement.

## Related folders

- [`../Auxiliar/`](../Auxiliar/) — `write_python_inputs` / `write_python_inputs_serot`, which produce the interface file consumed here.
- [`../Full_ephemeris_conversion/`](../Full_ephemeris_conversion/) — the MATLAB-side full-ephemeris refinement (with Moon gravity and TCMs) that this Python code is intended to supersede in fidelity.
- Root [`../README.md`](../README.md) — project overview; the GODOT bridge is described in §8.
