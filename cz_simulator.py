import numpy as np
import qutip as qt
from scipy.special import erf
from scipy.optimize import minimize

class TransmonCZSimulator:
    """
    Simulation and pulse optimization engine for an adiabatic Controlled-Z (CZ) gate
    between a pair of capacitively coupled transmon qutrits, modulated via external 
    magnetic flux pulses.
    """
    def __init__(self, Ec1=0.3, Ej1=11.7, Ec2=0.3, Ej2=13.5, g=0.015, dim=3):
        self.Ec1 = Ec1
        self.Ej1 = Ej1
        self.Ec2 = Ec2
        self.Ej2 = Ej2
        self.g = g
        self.dim = dim

        # Characteristic plasma frequencies (GHz)
        self.omega1 = np.sqrt(8 * Ec1 * Ej1)
        self.omega2 = np.sqrt(8 * Ec2 * Ej2)

        # Theoretical operating point for |11> <-> |02> avoided crossing (Delta = Ec2)
        self.u_target = (self.omega1 - self.Ec1 + 2 * self.Ec2) / self.omega2
        if np.abs(self.u_target) > 1.0:
            raise ValueError(
                "Incompatible physical parameters: Unable to reach the |11> <-> |02> degeneracy point. "
                "Ensure that qubit 2 frequency is properly scaled relative to qubit 1."
            )
        self.phi_target = (1.0 / np.pi) * np.arccos(self.u_target**2)
        self.t_gate_square = np.pi / (np.sqrt(2) * self.g)

        self._build_operators()
        self._build_hamiltonians()

    def _build_operators(self):
        """Constructs uncoupled qutrit tensor operators in the 9-dimensional Hilbert space."""
        I = qt.qeye(self.dim)
        a = qt.destroy(self.dim)
        n = qt.num(self.dim)
        proj2 = qt.fock_dm(self.dim, 2)

        self.a1 = qt.tensor(a, I)
        self.a2 = qt.tensor(I, a)
        self.n1 = qt.tensor(n, I)
        self.n2 = qt.tensor(I, n)
        self.proj2_1 = qt.tensor(proj2, I)
        self.proj2_2 = qt.tensor(I, proj2)

    def _build_hamiltonians(self):
        """Builds both the Rotating Frame Hamiltonian (for time evolution) and Lab Frame (for spectra)."""
        self.H_int = self.g * (self.a1.dag() * self.a2 + self.a1 * self.a2.dag())

        # Rotating Frame Hamiltonian at frequency omega_1 = omega_p1 - Ec1
        H_loc1 = -self.Ec1 * self.proj2_1
        H_loc2 = (self.Ec1 - self.Ec2 - self.omega1) * self.n2 - self.Ec2 * self.proj2_2
        self.H_0 = H_loc1 + H_loc2 + self.H_int
        self.H_u = self.omega2 * self.n2

        # Lab Frame Hamiltonian
        H_loc1_lab = self.omega1 * self.n1 - self.Ec1 * self.proj2_1
        H_loc2_lab = -self.Ec2 * self.proj2_2
        self.H_0_lab = H_loc1_lab + H_loc2_lab + self.H_int

    @staticmethod
    def u_flattop(t, args):
        """
        Smooth flat-top flux pulse profile modulated via error-function (erf) ramps.
        Required signature: f(t, args) for QuTiP time-dependent solvers.
        """
        T = args['T']
        sigma = args['sigma']
        amp = args['amp']
        phi_t = amp * 0.5 * (erf((t - 3 * sigma) / sigma) - erf((t - T + 3 * sigma) / sigma))
        return np.sqrt(np.cos(np.pi * np.clip(phi_t, -0.5, 0.5)))

    def compute_spectrum(self, phi_max=0.15, num_points=300):
        """Calculates the adiabatic energy levels in the Lab Frame."""
        phi_array = np.linspace(0, phi_max, num_points)
        energies = np.zeros((num_points, self.dim**2))
        
        for i, phi in enumerate(phi_array):
            u = np.sqrt(np.cos(np.pi * phi))
            H_tot = self.H_0_lab + u * self.H_u
            evals = H_tot.eigenenergies()
            energies[i, :] = evals - evals[0]
            
        return phi_array, energies

    def evaluate_pulse(self, T, sigma, amp=None):
        """
        Simulates time propagation and returns (Fidelity %, Leakage %) 
        after single-qubit Virtual-Z phase compensation.
        """
        if amp is None:
            amp = self.phi_target

        args_pulse = {'T': float(T), 'sigma': float(sigma), 'amp': float(amp)}
        H_tot = [self.H_0, [self.H_u, self.u_flattop]]

        U_T = qt.propagator(H_tot, [0, args_pulse['T']], args=args_pulse)[-1]
        U_full = U_T.full()

        # Computational subspace {|00>, |01>, |10>, |11>} -> indices [0, 1, 3, 4]
        comp_indices = [0, 1, 3, 4]
        U_comp = U_full[np.ix_(comp_indices, comp_indices)]

        # Global phase removal relative to |00> state
        global_phase = np.angle(U_comp[0, 0])
        U_comp = U_comp * np.exp(-1j * global_phase)

        # Single-qubit phase correction (Virtual-Z gates)
        phi_01 = np.angle(U_comp[1, 1])
        phi_10 = np.angle(U_comp[2, 2])
        VZ = np.diag([
            1.0,
            np.exp(-1j * phi_01),
            np.exp(-1j * phi_10),
            np.exp(-1j * (phi_01 + phi_10))
        ])
        U_corr = VZ @ U_comp

        # Gate Fidelity and Leakage metrics
        U_cz_ideal = np.diag([1.0, 1.0, 1.0, -1.0])
        trace_val = np.trace(U_cz_ideal.conj().T @ U_corr)
        fidelity = (np.abs(trace_val)**2) / 16.0

        quad_trace = np.trace(U_comp.conj().T @ U_comp).real
        leakage = 1.0 - (quad_trace / 4.0)

        return fidelity * 100.0, max(0.0, leakage * 100.0)

    def optimize_pulse(
        self,
        sigma_bounds=(4.0, 15.0),
        T_bounds=(None, 250.0),
        amp_tol=0.15
    ):
        """
        Performs gradient-free Nelder-Mead optimization over (T, sigma, amp)
        to maximize Gate Fidelity while enforcing physical constraints (T >= 6*sigma)
        and a 15% amplitude tuning range around the theoretical target.
        """
        sigma_min, sigma_max = sigma_bounds
        
        # Natural lower bound constraint ensuring flat-top integrity: T >= 6*sigma
        natural_T_min = 6.0 * sigma_min
        T_min = natural_T_min if T_bounds[0] is None else max(T_bounds[0], natural_T_min)
        T_max = T_bounds[1]

        init_sigma = np.clip(6.0, sigma_min, sigma_max)
        init_T = np.clip(max(self.t_gate_square * 1.1, 6.0 * init_sigma + 10.0), T_min, T_max)

        amp_bounds = (
            self.phi_target * (1.0 - amp_tol),
            min(0.49, self.phi_target * (1.0 + amp_tol))
        )
        
        x0 = [init_T, init_sigma, self.phi_target]
        bounds = [(T_min, T_max), (sigma_min, sigma_max), amp_bounds]

        def cost_function(params):
            T_val, sigma_val, amp_val = params
            fid, _ = self.evaluate_pulse(T_val, sigma_val, amp=amp_val)
            return 100.0 - fid  # Minimizing Infidelity (100 - Fidelity %)

        res = minimize(
            cost_function,
            x0=x0,
            bounds=bounds,
            method='Nelder-Mead',
            options={'xatol': 1e-3, 'fatol': 1e-4, 'maxiter': 150}
        )

        opt_T, opt_sigma, opt_amp = res.x[0], res.x[1], res.x[2]
        opt_fid, opt_leak = self.evaluate_pulse(opt_T, opt_sigma, amp=opt_amp)

        return {
            'T_opt': opt_T,
            'sigma_opt': opt_sigma,
            'amp_opt': opt_amp,
            'fidelity': opt_fid,
            'leakage': opt_leak
        }