# Robust RIS Reflection Design for Vehicular Platoon Communications under Angle Uncertainty

A MATLAB-based simulation framework for designing and evaluating a **Robust Reconfigurable Intelligent Surface (RIS)** reflection strategy for mmWave vehicular platoon communications under angular uncertainty.

The project investigates how robust RIS beamforming improves communication reliability in high-mobility vehicular environments compared to conventional non-robust and random RIS configurations. The optimization problem is formulated as a Semidefinite Relaxation (SDR) problem and solved using convex optimization.

---

## Overview

Millimeter-wave (mmWave) vehicular communications offer high data rates but suffer from severe path loss, beam misalignment, and channel uncertainty caused by vehicle mobility.

This project proposes a **covariance-based robust RIS reflection design** that maximizes the minimum beamforming gain across all vehicles under bounded angular uncertainty.

The implementation evaluates system performance using Monte Carlo simulations and compares the proposed robust RIS design with conventional non-robust and random RIS beamforming methods.

---

## Features

- Robust RIS reflection optimization
- Semidefinite Relaxation (SDR) formulation
- Covariance-based angular uncertainty modeling
- Monte Carlo simulation framework
- Worst-case beamforming optimization
- Comparison with:
  - Robust RIS
  - Non-Robust RIS
  - Random RIS
- Publication-quality performance plots

---

## Methodology
<img width="1024" height="559" alt="image" src="https://github.com/user-attachments/assets/dfe979bf-406f-4ebb-8154-79cf1a6d0e76" />


---

## Performance Metrics

The simulation generates the following results:

- Broadcast Rate vs Angle Error
- Broadcast Rate vs Vehicle Speed
- Outage Probability vs Angle Error
- Worst-case SNR CDF
- Broadcast Rate vs RIS Size
- Performance Loss vs Angle Uncertainty
- Computation Time vs RIS Size
- RIS Gain Heatmap

---
---

## Requirements

- MATLAB R2022a or later
- YALMIP Toolbox
- MOSEK Solver
- Optimization Toolbox

See **requirements.md** for installation details.

---

## Running the Project

1. Clone the repository

```bash
git clone https://github.com/yourusername/Robust-RIS-Vehicular-Platoon-Communication.git
```

2. Open MATLAB.

3. Add the project folder to the MATLAB path.

4. Open

```
src/robust_ris_simulation.m
```

5. Run the script.

The program automatically:

- Solves the robust RIS optimization problem
- Computes the non-robust and random baselines
- Performs Monte Carlo simulations
- Generates all performance plots
- Displays the numerical results table

---

## Applications

- Intelligent Transportation Systems (ITS)
- Vehicular Platoon Communications
- 5G/6G Wireless Networks
- mmWave Communications
- Smart Road Infrastructure
- Reconfigurable Intelligent Surfaces (RIS)

---

## Future Improvements

- Deep learning-based RIS optimization
- Multi-RIS deployment
- Hardware implementation
- Real-time beam tracking
- UAV-assisted RIS communications

---

## Authors

Harshini S

---
# Contributors

Varsha V

Subiksha

----
## License

This project is licensed under the MIT License.
