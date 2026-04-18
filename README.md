# VBI-modal-analysis

## Overview

This MATLAB program computes the dynamic response and coupled eigenvalues of single-degree-of-freedom sprung masses moving on a simply supported Bernoulli-Euler beam. The primary objective is to demonstrate how the coordinate definition parameter γ affects the system eigenvalues (coupled frequencies and damping ratios).

Two scenarios are provided:

- **Single-vehicle scenario**: One vehicle crosses the bridge. Used to validate the formulation and visualize the effect of γ on both vehicle and bridge coupled eigenvalues.
- **Periodic multi-vehicle scenario**: An idealized infinite sequence of equally spaced vehicles travels across the bridge at constant speed. This enables eigenvalue analysis under sustained vehicle loading and produces bridge acceleration signals of arbitrary duration.

## File Structure

Both scenarios share the same file naming convention (`f00`–`f06`) and are placed in separate folders.

| File | Description |
| --- | --- |
| `f00_mainProgram.m` | Main script — run this file to execute the full analysis |
| `f01_loadInputParams.m` | Defines bridge, vehicle, solver, and fleet parameters |
| `f02_buildElementMatrices.m` | Builds beam element matrices and computes uncoupled frequencies |
| `f03_assembleGlobalMatrices.m` | Assembles time-varying global matrices for the coupled VBI system |
| `f03a_hermiteShape.m` | Computes Hermite shape functions and their derivatives |
| `f03b_VBIElement.m` | Generates VBI element matrices based on the formulation with γ |
| `f04_solveEigenValue.m` | Solves the coupled eigenvalue problem at each time step |
| `f05_solveNewmarkBeta.m` | Solves the dynamic response using the Newmark-β method |
| `f06_postProcess.m` | Post-processing and visualization |

## How to Use

1. Open MATLAB and set the working directory to the desired scenario folder.
2. Open `f01_loadInputParams.m` and set the key parameters:
   - `sysParams.solver.gamma = 0;` for the conventional coordinate selection
   - `sysParams.solver.gamma = 1;` for the proposed coordinate selection
   - For the multi-vehicle scenario, additionally set `d_spacing` (vehicle spacing), `T_total` (simulation duration), and `x0_lead` (initial position of the leading vehicle).
3. Run `f00_mainProgram.m`.


## Requirements

- MATLAB (tested with R2022a)
- No additional toolboxes required

## Reference

To be added after publication.

## Author

Shi LuoKe (史罗克), 2026

## Acknowledgment

This code was developed by the author during a visiting PhD period at the Norwegian University of Science and Technology (NTNU). The program is based on a conventional VBI code originally provided by Xu Hao, and the theoretical interpretation of the formulation was completed under the guidance of the author's supervisor, Prof. Y.B. Yang. The author would like to express special gratitude to Daniel Cantero for his hospitality, support, and valuable suggestions during the visit at NTNU.

To improve code readability, the code was organized and formatted with the assistance of Claude (Anthropic). The correctness of the code has been verified by the author.
