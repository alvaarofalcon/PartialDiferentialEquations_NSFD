# SEIQR Reaction–Diffusion (MATLAB) — 1D and 2D Implementations

## Overview
This repository contains two MATLAB implementations of a **spatial SEIQR reaction–diffusion model** (Susceptible–Exposed–Infected–Quarantined–Recovered):

- **1D version** (one spatial dimension): evolution of `S, E, I, Q, R` along a 1D grid.
- **2D version** (two spatial dimensions): evolution of `S, E, I, Q, R` on a 2D grid (heatmap visualization).

Both versions share the same modeling idea:
- Each compartment is treated as a **spatial density field**.
- **Diffusion** spreads each compartment in space.
- **Reaction terms** couple the compartments (infection, progression, quarantine, recovery, removal, etc.).
- **Neumann (zero-flux) boundary conditions** are imposed so there is no flow across the domain boundary.
- The simulation stores snapshots every `save_interval` steps to reduce memory usage and to enable visualization at selected times.

---

## Repository Structure
This repo is organized into two main folders:

- `1D/`  
  MATLAB code for the **1D SEIQR reaction–diffusion model** (profiles vs. space).

- `2D/`  
  MATLAB code for the **2D SEIQR reaction–diffusion model** (2D maps / heatmaps).

Each folder contains:
- a main entry-point script/function (typically `SEIQR_main.m`)
- a solver/integration routine (typically `modelo_SEIQR.m`)
- optional helper utilities (depending on the folder)

> File names may be the same in both folders; run each version from inside its own folder (or add only that folder to the MATLAB path) to avoid name collisions.

---

## Model Components (Common to 1D and 2D)
### Compartments
- `S` — Susceptible
- `E` — Exposed
- `I` — Infected
- `Q` — Quarantined
- `R` — Recovered

### Parameters
Model parameters are defined in a `params` struct (names as used in the code), typically including:
- Diffusion coefficients: `lambda_S, lambda_E, lambda_I, lambda_Q, lambda_R`
- Reaction rates: `beta1, beta2, mu, delta, gamma, alfa, rho`
- Recruitment term: `Lambda`
- Discretization: `h` (spatial step), `k` (time step)

> Unless explicitly rescaled, time and space are treated in **model units** (dimensionless), as is common in numerical PDE experiments.

---

## Numerical Method (High-Level)
### Spatial diffusion
- **1D:** nearest-neighbor coupling along the line (left/right neighbors).
- **2D:** 4-neighbor stencil (Von Neumann neighborhood).

### Boundary conditions
Both versions impose **Neumann (zero-flux)** boundary conditions, meaning:
\[
\frac{\partial u}{\partial n} = 0 \quad \text{on the boundary}
\]
so there is **no net flow across the boundary**. In discrete form, boundaries are enforced by reflecting interior values into the boundary cells.

### Time stepping and storage
- The solver advances the system in time for `total_steps`.
- Snapshots are stored every `save_interval` steps in cell arrays such as:
  `S_hist, E_hist, I_hist, Q_hist, R_hist`.

---

## Visualization
- **1D version:** typically plots spatial profiles `u(x)` (line plots) at a selected time.
- **2D version:** displays 2D heatmaps using `imagesc` for `S, E, I, Q, R` at a selected time, usually with a shared colorbar/scale.

---

## How to Run
### Run the 1D version
1. Open MATLAB and change directory into `1D/`.
2. Run the main file:
```matlab
SEIQR_main
