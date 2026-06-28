# Requirements

## Software

- MATLAB R2022a or newer

## MATLAB Toolboxes

- Optimization Toolbox

## Third-Party Packages

### YALMIP

YALMIP is required for formulating the semidefinite optimization problem.

Repository:

https://github.com/yalmip/YALMIP

Installation:

1. Download YALMIP.
2. Extract the folder.
3. Add it to the MATLAB path.

```matlab
addpath(genpath('YALMIP'))
savepath
```

---

### MOSEK Solver

MOSEK is required for solving the Semidefinite Programming (SDP) optimization problem.

Website:

https://www.mosek.com/

Installation Steps

1. Download MOSEK.
2. Install the software.
3. Activate an Academic License (free for students).
4. Verify installation.

```matlab
mosekdiag
```

---


## Running

Open

```
src/robust_ris_simulation.m
```

Run the script.

The simulation automatically:

- Performs robust RIS optimization
- Solves the non-robust baseline
- Generates all plots
- Displays the numerical results table
