# Adiabatic Controlled-Z (CZ) Gate in Capacitively Coupled Transmon Qutrits

Simulation and pulse-level optimization of a flux-tunable Controlled-Z (CZ) quantum gate. The gate is implemented using the second excited state ($|2\rangle$) of a capacitively coupled superconducting transmon pair, simulated dynamically using QuTiP.

📄 **[Read the Full Theoretical Derivation (PDF)](CZ_GATE_TRANSMON.pdf)**

## Physical Model & Overview

In superconducting circuit quantum electrodynamics (cQED), a high-fidelity two-qubit Controlled-Z gate can be realized by adiabatically tuning the magnetic flux $\Phi_{\text{ext}}(t)$ threading a DC-SQUID transmon. This brings the computational state $|11\rangle$ into resonance with the non-computational qutrit level $|02\rangle$.

1. **Circuit Quantization & RWA**: Starting from the classical capacitance matrix in the weak-coupling regime ($C_c \ll C_1, C_2$), the system is quantized in the transmon regime ($E_J / E_C \gg 50$) and truncated to the three lowest energy levels $\{|0\rangle, |1\rangle, |2\rangle\}$ per site (resulting in a 9-dimensional Hilbert space).
2. **Rotating Frame Decomposition**: Moving to a frame rotating at the first qubit's transition frequency $\omega_1$, the time-dependent Hamiltonian under the Rotating Wave Approximation (RWA) is expressed in standard control form:
   $$\hat{H}(t) = \hat{H}_0 + u(t)\hat{H}_u, \qquad u(t) = \sqrt{\cos\left(\pi \frac{\Phi_{\text{ext}}(t)}{\Phi_0}\right)}$$
3. **Adiabatic Flat-Top Pulse Shaping**: To suppress diabatic leakage into $|02\rangle$, the ideal square pulse is replaced by a smooth error-function ($\text{erf}$) envelope parameterized by total gate duration $T$ and ramp time $\sigma$:
   $$\Phi_{\text{ext}}(t) = \frac{\bar{\Phi}_{\text{ext}}}{2}\left[\text{erf}\left(\frac{t - 3\sigma}{\sigma}\right) - \text{erf}\left(\frac{t - T + 3\sigma}{\sigma}\right)\right]$$
4. **Virtual-Z Phase Compensation**: Single-qubit dynamical phases accumulated during the flux excursion are mathematically tracked and canceled via Virtual-Z rotations prior to evaluating the final gate fidelity $\mathcal{F}$ and leakage $\mathcal{L}$.

## Repository Structure

This project uses a flat directory structure for simplicity:

* `cz_simulator.py`: Core simulation engine (`TransmonCZSimulator`). It handles qutrit operators, Lab/Rotating frame Hamiltonians, propagator integration, and the gradient-free Nelder-Mead pulse optimization.
* `tutorial.ipynb`: Interactive Jupyter Notebook illustrating the energy spectrum, the $|11\rangle \leftrightarrow |02\rangle$ avoided crossing, and the optimization pipeline.
* `optimize_pulse.py`: Command-line interface (CLI) tool. It uses `argparse` to let you quickly find the optimal $(T, \sigma, \text{amp})$ pulse parameters for custom transmon hardware specifications without modifying the source code.
* `QCED_CZ_PRESENTATION.pdf`: Presentation slides containing the complete step-by-step analytical derivation of the model.

## Quickstart & Usage

### Prerequisites
Make sure you have the following libraries installed:
```bash
pip install numpy scipy qutip matplotlib
```

### 1. Interactive Walkthrough
To visualize the avoided crossing and understand the gate mechanics, open and run the Jupyter notebook:
```bash
jupyter notebook tutorial.ipynb
```

### 2. Automated Pulse Optimization (CLI)
Run the command-line optimizer with default transmon parameters ($E_{C1}=E_{C2}=0.3\text{ GHz}$, $E_{J1}=11.7\text{ GHz}$, $E_{J2}=13.5\text{ GHz}$, $g=0.015\text{ GHz}$):
```bash
python optimize_pulse.py
```
Or specify custom hardware parameters to tune your own simulated device:
```bash
python optimize_pulse.py --Ec1 0.29 --Ej1 12.7 --Ec2 0.30 --Ej2 14.2 --g 0.018 --sigma-min 5.0 --sigma-max 10.0 --ampl-tol 0.15
```