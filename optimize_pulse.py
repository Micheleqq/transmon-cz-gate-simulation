import argparse
from cz_simulator import TransmonCZSimulator

def main():
    parser = argparse.ArgumentParser(
        description="CLI tool to automatically optimize Flat-Top pulse parameters (T, sigma, amp) for a CZ gate."
    )
    parser.add_argument("--Ec1", type=float, default=0.3, help="Charging energy Qubit 1 (GHz)")
    parser.add_argument("--Ej1", type=float, default=11.7, help="Josephson energy Qubit 1 (GHz)")
    parser.add_argument("--Ec2", type=float, default=0.3, help="Charging energy Qubit 2 (GHz)")
    parser.add_argument("--Ej2", type=float, default=13.5, help="Josephson energy Qubit 2 (GHz)")
    parser.add_argument("--g", type=float, default=0.015, help="Capacitive coupling strength (GHz)")
    parser.add_argument("--sigma-min", type=float, default=4.0, help="Lower bound for ramp time sigma (ns)")
    parser.add_argument("--sigma-max", type=float, default=15.0, help="Upper bound for ramp time sigma (ns)")
    parser.add_argument("--amp-tol", type=float, default=0.15, help="Amplitude tuning tolerance fraction (default: 0.15 = 15%)")

    args = parser.parse_args()

    print("=== Initializing Transmon Qutrit CZ Simulator ===")
    sim = TransmonCZSimulator(Ec1=args.Ec1, Ej1=args.Ej1, Ec2=args.Ec2, Ej2=args.Ej2, g=args.g)
    print(f"Plasma Frequency Qubit 1 (omega_p1): {sim.omega1:.4f} GHz")
    print(f"Plasma Frequency Qubit 2 (omega_p2): {sim.omega2:.4f} GHz")
    print(f"Theoretical Target Flux (Phi/Phi0): {sim.phi_target:.4f}")
    print(f"Square Pulse Analytical Time (T)  : {sim.t_gate_square:.2f} ns\n")

    print(f"Running Nelder-Mead optimization over (T, sigma, amp) ...")
    result = sim.optimize_pulse(
        sigma_bounds=(args.sigma_min, args.sigma_max),
        amp_tol=args.amp_tol
    )

    print("\n=== Optimization Results ===")
    print(f"Optimal Gate Duration (T)    : {result['T_opt']:.2f} ns")
    print(f"Optimal Ramp Time (sigma)    : {result['sigma_opt']:.2f} ns")
    print(f"Optimal Flux Amplitude (amp) : {result['amp_opt']:.4f}")
    print(f"Max Gate Fidelity            : {result['fidelity']:.4f} %")
    print(f"Residual Leakage             : {result['leakage']:.4f} %")

if __name__ == "__main__":
    main()