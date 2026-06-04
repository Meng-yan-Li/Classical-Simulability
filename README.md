# Classical-Simulability

This repository contains MATLAB code for computing and analysing the "critical visibility" of families of quantum states with respect to classical simulation strategies. The critical visibility quantifies the largest white-noise admixture `v` in `[0,1]` for which a noisy ensemble `v * rho + (1-v)/d * I` admits a decomposition compatible with classical simulability constraints. The code implements both canonical and aggregated SDP relaxations and provides example state families and precomputed numeric results.

**Contents**

- **Files:**
  - [main_CS.m](main_CS.m) : Orchestrator script. Runs a set of experiments that compute critical visibilities for 2-dimensional state families and higher-dimensional state families, runs classical-simulability checks, and saves results into `Num_Data`.
  - [Critical_visibility_aggregate.m](Critical_visibility_aggregate.m) : Aggregated SDP hierarchy implementation (YALMIP). Builds an aggregated relaxation and returns the optimal critical visibility and solver status.
  - [Critical_visibility_canonical.m](Critical_visibility_canonical.m) : Canonical-lift SDP formulation (YALMIP). Returns `CV` plus structural variables (`Tau`, `Pi`, `c`) evaluated after solving the SDP.
  - [Critical_visibility_canonical_cvx.m](Critical_visibility_canonical_cvx.m) : CVX-based canonical formulation providing an alternative solver interface (CVX + MOSEK or other compatible solver).
  - [generate_2d_state_family.m](generate_2d_state_family.m) : Produces common qubit state families: BB84, Six-state, Trine, and Tetrahedral (SIC) families.
  - [generate_quantum_sets.m](generate_quantum_sets.m) : Produces higher-dimensional test sets (e.g., qutrit MUBs, SIC states, hybrid diagonal families, and random families).
- **Data folder:** `Num_Data/` contains saved `.mat` files produced by example runs. Example files included: `SimulationForBB84.mat`, `SimulationForSixstate.mat`, `WitnessForBB84.mat`, etc.

**Purpose & Approach**

- **Goal:** quantify when a family of quantum states can be simulated classically (by deterministic outcome assignments combined with shared randomness) under noise. The project computes the maximal visibility `v` such that noisy states remain classically simulable.
- **Approach:** express simulability conditions as semidefinite constraints on structural variables (`Tau`, `Pi`, `c`) and maximize `v` via SDP. Two closely related formulations are implemented: a canonical-lift approach (fine-grained per-`mu` variables) and an aggregated hierarchy that reduces some redundancy for numerical scaling.

**How to run**

- **Prerequisites:** MATLAB (R2016b or later recommended), YALMIP (for YALMIP-based scripts), MOSEK (recommended) or CVX (for the CVX variant). Add solver/toolbox paths to MATLAB's search path.

- **Run examples:** set current folder to the repository root, then run:

  main_CS.m

  The script runs through example sections: 2D families (BB84/Six-state/Trine/Tetrahedral), higher-dimensional families via `generate_quantum_sets`, optional classical-simulation checks, and witness computations. Results are saved under `Num_Data/`.

**Key variables and outputs**

- `CV` : the computed critical visibility (scalar in [0,1]).
- `Tau`, `Pi`, `c` : structural variables used in the SDP; canonical routine returns evaluated numeric values after solving.
- `W_opt` : dual witness operators (from CVX/YALMIP duals) used to verify or construct entanglement/measurement witnesses.
- Saved `.mat` files under `Num_Data/` contain the state families, computed CVs, runtime statistics, and optionally witness data.

**Notes for developers / reproducibility**

- The YALMIP routines use `sdpsettings` and pass options to MOSEK (e.g., `opts.threads`, `opts.gap`). Tune these options in the top of `main_CS.m` before running large experiments.
- For very large `d` or many states `n`, the SDP sizes grow quickly with the chosen hierarchy level `ell`. Use the aggregated formulation (`Critical_visibility_aggregate.m`) to reduce memory where possible.
- The CVX variant (`Critical_visibility_canonical_cvx.m`) provides an alternative modelling path; use it if you prefer CVX syntax or different solver options.

**Files of interest**

- Example entrypoint: [main_CS.m](main_CS.m)
- Core solver implementations: [Critical_visibility_canonical.m](Critical_visibility_canonical.m), [Critical_visibility_aggregate.m](Critical_visibility_aggregate.m), [Critical_visibility_canonical_cvx.m](Critical_visibility_canonical_cvx.m)
- State generators: [generate_2d_state_family.m](generate_2d_state_family.m), [generate_quantum_sets.m](generate_quantum_sets.m)
- Precomputed results: `Num_Data/` (see `Num_Data/SimulationForBB84.mat`, `Num_Data/WitnessForBB84.mat`, etc.)

