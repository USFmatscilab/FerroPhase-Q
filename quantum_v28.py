"""

-----------------------------------------------------------------------------------
Simulation of ferroelectric polarization dynamics

This code evolves the density matrix and analytical wavefunction of a
ferroelectric system under an applied electric field (AC or DC). The dynamics
include relaxation processes, which can be rescaled by 'gamma_rescale'.

Inputs:
  - Simulation parameters are read from PARAMS_dyn.json
  - Key options include temperature, field type/strength, time step,
    relaxation model, gamma rescaling, and plotting mode.

Outputs:
  1. simulation_output.txt   → time series of observables:
        * time [fs]
        * applied electric field [kV/cm]
        * polarization (analytical) [C/m²]
        * polarization (density matrix) [C/m²]
        * trace of density matrix (real)

  2. polarization_vs_time.png   → polarization expectation vs time
  3. polarization_vs_field.png  → polarization vs applied electric field (if AC field)

Plotting behavior:
  - plot_mode = "screen" → interactive plots shown + .png files saved
  - plot_mode = "none"   → no interactive plots, only .png files saved

All simulation parameters from PARAMS_dyn.json are echoed below for reproducibility.

-----------------------------------------------------------------------------------
"""



import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
#from scipy.sparse.linalg import eigsh
from scipy.sparse import issparse, csr_matrix, diags
import sys
import h5py
import os
import json
import time




# Physical constants
k_B = 8.617333262145e-5  # eV/K
hbar_eVfs = 0.6582  # ħ in eV·fs
start_time = time.time()
eps = 1e-8  # minimum populations below which analytical model does not work
use_analytical_model = True 


def read_interpolated_potential(filename="interpolated_potential.txt"):
    """Read external potential from file."""
    data = np.loadtxt(filename)
    x = data[:, 0]
    V = data[:, 1]
    return x, V


def read_hdf5_results(filename):
    """Read results from HDF5 file written by save_all_data."""
    with h5py.File(filename, "r") as f:
        results = {}

        # Root-level data
        results["field_values"] = f["field_values"][()]
        results["polarization_data"] = f["polarization_data"][()]
        results["z_star"] = f["z_star"][()]
        results["volume"] = f["volume"][()]
        results["num_solutions"] = f["num_solutions"][()]

        # Optional: read attributes
        results["mass"] = f.attrs["mass"]
        results["n"] = f.attrs["n"]
        results["num_sol_attr"] = f.attrs["num_sol"]
        results["Star"] = f.attrs["Star"]
        results["Vol"] = f.attrs["Vol"]

        # Read eigenvalues, eigenvectors, and potentials for each field
        eigen_group = f["eigen_data"]
        num_fields = len(eigen_group)
        eigenvalues = []
        eigenvectors = []
        potentials = []

        for i in range(num_fields):
            field_key = f"field_{i:03d}"
            group = eigen_group[field_key]
            eigenvalues.append(group["energies"][:])
            eigenvectors.append(group["eigenvectors"][:])
            potentials.append(group["potential"][:])

        results["eigenvalues"] = np.array(eigenvalues)
        results["eigenvectors"] = np.array(eigenvectors)
        results["potentials"] = np.array(potentials)

    return results




def comp_v_E(E):
    """Compute potential under electric field E (kV/cm)"""
    return v_int - z_star * E * x_int * 1.e-04

def build_hamiltonian(length, n, v, mass):
    """Construct the Hamiltonian matrix"""
    hbar = 1.054571596
    echarge = 1.602176462
    baremass = 9.10938188

    const = hbar ** 2 / (baremass * echarge) / mass
    deltax = length / n
    const /= deltax ** 2

    diagonal = v[1:-1] + const
    offdiagonal = const / 2.0

    diagonals = [diagonal, np.full(n - 2, -offdiagonal), np.full(n - 2, -offdiagonal)]
    return diags(diagonals, [0, -1, 1], format='csr')





def compute_psi_t(energy, phi, c_n0, t_fs, intrinsic_dynamics=True):
    """
    Compute time-evolved wavefunction ψ(t) = Σ c_n(0) e^(-iωₙt) ϕₙ
    with exponential decay
    
    Args:
        energies (array): Eigenvalues in eV (shape: (num_states,))
        phi (array): Eigenvectors (shape: (num_points, num_states))
        c_n0 (array): Initial coefficients (shape: (num_states,))
        t_fs (float): Time in femtoseconds
        
    Returns:
        psi_t (array): Wavefunction at time t (shape: (num_points,))
    """
    # Convert energies to angular frequencies (2π THz)
    hbar_eVfs = 0.6582119569  # ħ in eV·fs (conversion factor) 0.6582119569
    gamma_n=(energy-energy[0])/hbar_eVfs 
    if intrinsic_dynamics:
        omega_n = energy / hbar_eVfs
    else:
        omega_n = 0.0
    
    omega_n_prime=omega_n-1j*0.5*gamma_n
    
    
    # Compute phase factors for all states
    phase_factors = np.exp(-1j * omega_n_prime * t_fs)
    
    # Combine coefficients with phase factors and project to position basis
    psi_t = np.dot(phi, c_n0 * phase_factors)
    
    norm = np.linalg.norm(psi_t)
    psi_t = psi_t / norm
    
    return psi_t



def commutator(H, rho):
    """
    Compute the commutator (-i/ħ)[H, ρ] where H is sparse.
    
    Args:
        H (sparse matrix): Hamiltonian (shape [n, n])
        rho (array or sparse): Density matrix (shape [n, n])
        hbar (float): Reduced Planck constant (default 1.0)
    
    Returns:
        comm (array or sparse): Commutator (same type as rho)
    """
    # Ensure H is sparse (if not, convert to CSR)
    hbar = 0.6582119569 # ev fs 
    H = csr_matrix(H) if not issparse(H) else H
    
    # Compute [H, ρ] = Hρ - ρH
    if issparse(rho):
        H_rho = H.dot(rho)
        rho_H = rho.dot(H)
    else:
        H_rho = H.dot(rho)
        rho_H = rho.dot(H.toarray())  # Efficient for dense rho
    
    comm = (-1j / hbar) * (H_rho - rho_H)
    
    return comm





def rel_term(rho_t_loc, energy_loc, phi, Temp, relaxation_model,gamma_rescale):
    """
    Compute relaxation term for the density matrix in the eigenbasis of H(E).
    Includes population and coherence decay depending on the model.
    """

    hbar = 0.6582119569  # eV·fs
    k_B = 8.617333262145e-5  # eV/K
    
    beta = 1.0/(k_B*Temp)       
    gamma_n = gamma_rescale*(energy_loc - energy_loc[0]) / hbar  # Decay rates γ_n = (E_n - E_0)/ħ    
    num_states = len(energy_loc)
    
    # Transform to energy eigenbasis
    rho_eig = phi.conj().T @ rho_t_loc @ phi
    d_rho_eig = np.zeros_like(rho_eig, dtype=complex)

    # For rlx_equil_dens, use equilibrium density
    if relaxation_model == "DMB":
        exp_terms = np.exp(-beta * (energy_loc - energy_loc[0]))
        rho_eq_eig = np.diag(exp_terms / np.sum(exp_terms))
        for n in range(0,num_states):
            d_rho_eig[n, n] = -gamma_n[n] * (rho_eig[n, n] - rho_eq_eig[n, n])
        # Off-diagonal coherence decay
        for n in range(num_states):
            for m in range(num_states):
                if n != m:
                    d_rho_eig[n, m] -= 0.5 * (gamma_n[n] + gamma_n[m]) * rho_eig[n, m]    
                 

    # For rlx_lindblad_grn: ground state gets populations from excited states
    elif relaxation_model == "LLO":
        exp_terms = np.exp(-beta * (energy_loc - energy_loc[0]))
        for n in range(num_states):
            if n == 0:
                for m in range(1, num_states):
                    d_rho_eig[0, 0] +=  2.*(gamma_n[m]) * (rho_eig[m, m] - exp_terms[m] * rho_eig[0, 0])
            else:
                d_rho_eig[n, n] = -2.*(gamma_n[n]) * (rho_eig[n, n] - exp_terms[n] * rho_eig[0, 0])
        # Off-diagonal coherence decay 
        for n in range(num_states):
            for m in range(num_states):
                if n != m:
                    d_rho_eig[n, m] -=  (gamma_n[n] + gamma_n[m]) * rho_eig[n, m]  


    elif relaxation_model == "LO":
        num_states = len(energy_loc)
       
        # Compute transition rates γ_kn = |E_k - E_n| / ħ
        gamma_mat = np.zeros((num_states, num_states))

        for n in range(num_states):
            for k in range(n+1, num_states):  # only upper triangle
                deltaE = energy_loc[k] - energy_loc[n]
                val = np.abs(deltaE) / hbar
                gamma_mat[k, n] = val  # lower triangle
                gamma_mat[n, k] = val  # symmetric upper triangle
                gamma_mat[k, n] *= np.exp(-beta * deltaE)

        gamma_mat = gamma_rescale*gamma_mat 
        gamma_col_sum = np.sum(gamma_mat, axis=0)[:, None]  # column vector of row sums (shape: N×1)
               
    
        # --- Population terms ---
        for n in range(num_states):
            inflow = np.sum(gamma_mat[n, :] * np.diag(rho_eig))
            outflow = gamma_col_sum[n, 0]
            d_rho_eig[n, n] = 2.0 * (inflow - outflow * rho_eig[n, n])
    
        # --- Coherence terms ---
        for n in range(num_states):
            for m in range(n+1, num_states):  # only upper triangle
                val = - (gamma_col_sum[n, 0] + gamma_col_sum[m, 0]) * rho_eig[n, m]
                d_rho_eig[n, m] = val
                d_rho_eig[m, n] = np.conj(val)  # enforce Hermitian                      
    else:
        raise ValueError(f"Unknown relaxation model: {relaxation_model}")


    # Transform back to original basis
    d_rho_t = phi @ d_rho_eig @ phi.conj().T
    return d_rho_t



def predictor(rho_t,rv,rv1,rv2,delta_t):
    rho_0 = np.copy(rho_t)
    rho_t[1:-1,1:-1]+=(delta_t/12.)*(23.0*rv-16.0*rv1+5*rv2)
    rv2 = rv1.copy()
    rv1 = rv.copy()
    rho_t /= np.trace(rho_t)             # Force Tr(ρ) = 1
    return rho_t,rho_0,rv1,rv2

def corrector(rho_0,rv,rv1,rv2,delta_t):
    rho_t_loc = np.zeros_like(rho_0)
    rho_t_loc[1:-1,1:-1]=rho_0[1:-1,1:-1]+ (delta_t/12.)*(5.0*rv+8.0*rv1-rv2)
    rho_t_loc /= np.trace(rho_t_loc)             # Force Tr(ρ) = 1
    return rho_t_loc


def predictor_psi(psi,v_psi,v1_psi,v2_psi,delta_t):
    psi_0 = np.copy(psi)
    psi[1:-1]+=(delta_t/12.)*(23.0*v_psi-16.0*v1_psi+5*v2_psi)
    v2_psi=v1_psi
    v1_psi=v_psi
    norm = np.linalg.norm(psi)
    psi = psi / norm
    return psi,psi_0,v1_psi,v2_psi

def corrector_psi(psi,psi_0,v_psi,v1_psi,v2_psi,delta_t):
    psi[1:-1]=psi_0[1:-1]+ (delta_t/12.)*(5.0*v_psi+8.0*v1_psi-v2_psi)
    norm = np.linalg.norm(psi)
    psi = psi / norm
    return psi

def test_dt(H,delta_t):
    max_element_H = np.max(np.abs(H)) 
    delta_t_max=4.135/max_element_H
    if delta_t < delta_t_max:
        print(f"✅ Time step SAFE (δt = {delta_t:.3f} fs)")
        print(f"   Maximum stable time step = {delta_t_max:.3f} fs")
        print(f"   Ratio δt/δt_max = {delta_t/delta_t_max:.3f} (<1 is safe)")
        return True
    else:
        print(f"⚠️ WARNING: Time step UNSAFE (δt = {delta_t:.3f} fs)")
        print(f"   Maximum stable time step = {delta_t_max:.3f} fs")
        print(f"   Ratio δt/δt_max = {delta_t/delta_t_max:.3f} (≥1 may cause instability)")
        print("   Recommendation: Reduce time step or use implicit integration methods")
        return False
 

# Project density matrix onto the eigenbasis subspace
def project_density_matrix(rho_loc, phi_loc):
    """
    Projects density matrix rho onto subspace spanned by eigenvectors phi
    
    Args:
        rho: Density matrix in position basis (shape [n+1, n+1])
        phi: Eigenvectors (shape [n+1, num_sol])
    
    Returns:
        Projected and normalized density matrix
    """
    # Create projector P = Σ|φₙ›‹φₙ|
    projector = np.sum([np.outer(phi_loc[:,n], phi_loc[:,n].conj()) for n in range(phi_loc.shape[1])], axis=0)
    
    # Apply projection: ρ_proj = PρP†
    rho_proj = projector @ rho_loc @ projector.conj().T
    
    # Normalize to maintain Tr(ρ) = 1
    rho_proj /= np.trace(rho_proj)
    
    return rho_proj

def initialize_state(x_int, current_phi, conv_P2x, P_cen, delta_P):
    x_cen = conv_P2x * P_cen
    sigma = conv_P2x * delta_P

    # Gaussian wavepacket (real at this stage, but will be cast to complex later)
    psi_0 = np.exp(-0.5 * ((x_int - x_cen) / sigma)**2)
    psi_0 = psi_0.astype(np.complex128)   # ensure complex dtype
    psi_0 /= np.linalg.norm(psi_0)

    # Density matrix
    rho_0 = np.outer(psi_0.conj(), psi_0)
    rho_0 = project_density_matrix(rho_0, current_phi)

    # Expansion coefficients in eigenbasis
    c_n0 = current_phi.conj().T @ psi_0
    c_n0 /= np.linalg.norm(c_n0)

    # Reconstruct initial state from eigenbasis (guarantees consistency)
    psi_0 = current_phi @ c_n0

    return psi_0, c_n0, rho_0


def read_json_parameters(json_filename):
    """Read simulation parameters from JSON file."""
    with open(json_filename, "r") as f:
        params = json.load(f)
    
    # Provide defaults where appropriate
    params_out = {
        "Temp": params.get("Temp"),
        "field_type": params.get("field_type", "linear"),
        "E_max": params.get("E_max"),
        "E_min": params.get("E_min"),
        "num_periods": params.get("num_periods", 1.25),
        "nsteps": params.get("nsteps", 10000),
        "delta_t": params.get("delta_t", 0.085),
        "intrinsic_dynamics": params.get("intrinsic_dynamics", True),
        "P_cen": params.get("P_cen", 0.0),
        "delta_P": params.get("delta_P", 0.1),
        "n_out": params.get("n_out", 500),
        "relaxation_model": params.get("relaxation_model", "DMB"),
        "gamma_rescale": params.get("gamma_rescale", 1.0),
        "plot_mode": params.get("plot_mode", "screen"),
        "integrator": params.get("integrator", "simple") 
    }
    return params_out


# Read interpolated potential from file
x_int, v_int = read_interpolated_potential()
length = x_int[-1] - x_int[0]

# Read results from HDF5 file
data = read_hdf5_results("all_simulation_data.h5")
fields = data["field_values"]         
eigval_data = data["eigenvalues"]  
eigvec_data = data["eigenvectors"]   
# Access HDF5 parameters
z_star = data["z_star"]
Vol = data["volume"]
num_sol = data["num_solutions"]
n = data["n"] 
mass = data["mass"]


json_params = read_json_parameters("PARAMS_dyn.json") # read PARAMS file
Temp = json_params["Temp"]
field_type = json_params["field_type"]
E_max = json_params["E_max"]
E_min = json_params["E_min"]
num_periods = json_params["num_periods"]
nsteps = json_params["nsteps"]
delta_t = json_params["delta_t"]
intrinsic_dynamics = json_params["intrinsic_dynamics"]
P_cen = json_params["P_cen"]
delta_P = json_params["delta_P"]
n_out = json_params["n_out"]
relaxation_model = json_params["relaxation_model"]
gamma_rescale = json_params["gamma_rescale"]
integrator = json_params["integrator"] 
freq = num_periods*1.e+06/(nsteps*delta_t) #frequency of ac field in GHz

conv_P2x = Vol * 0.01 / (z_star * 1.6)




# === Generate field_sequence ===
if field_type == "linear":
    n_ramp = nsteps // 4
    field_sequence = np.concatenate([
        np.linspace(0, E_max, n_ramp),
        np.linspace(E_max, E_min, 2 * n_ramp),
        np.linspace(E_min, E_max, n_ramp)
    ])

    if len(field_sequence) < nsteps:
        field_sequence = np.pad(field_sequence, (0, nsteps - len(field_sequence)), mode='edge')

elif field_type == "ac":
    t_fs = np.arange(nsteps) * delta_t
    total_time_fs = nsteps * delta_t
    omega = 2 * np.pi * num_periods / total_time_fs
    field_sequence = E_max * np.sin(omega * t_fs)
elif field_type == "dc":
    field_sequence = np.full(nsteps, E_max)
else:
    raise ValueError(f"Unknown field_type: {field_type}")


E_current = 0.0 #field_sequence[0]  # Start at zero field


field_idx = np.argmin(np.abs(fields - E_current))
current_energy = eigval_data[field_idx]
current_phi = eigvec_data[field_idx]


#Initialize with Gaussian and project onto the current basis
psi_0, c_n0, rho_0 = initialize_state(x_int, current_phi, conv_P2x, P_cen, delta_P) 
rho_t=rho_0.astype(np.complex128)
psi_t = np.copy(psi_0) # time dependent state vector for analytical solution of Schrod eq.
# Compute gamma_n = (E_n - E_0)/ħ just for the output
gamma_n = (current_energy - current_energy[0]) / hbar_eVfs

# Avoid division by zero for the ground state
valid = gamma_n > 1e-12
tau_n = np.zeros_like(gamma_n)
tau_n[valid] = 1.0 / gamma_n[valid]
tau_n[~valid] = np.inf  # or set to 0, depending on preference

# Estimate min/max relaxation times (excluding ground state)
tau_min = np.min(tau_n[valid])
tau_max = np.max(tau_n[valid])

rel = rel_term(rho_t, current_energy, current_phi, Temp,relaxation_model,gamma_rescale)[1:-1, 1:-1]

if intrinsic_dynamics:
    v_E = comp_v_E(E_current)
    H = build_hamiltonian(length, n, v_E, mass)
    rv = commutator(H, rho_t[1:-1, 1:-1])+rel
    safeStep=test_dt(H,delta_t)
    if safeStep==False: 
       sys.exit("Error: safeStep is False. Execution stopped.")
else:
    rv=rel
rv1=np.zeros_like(rv); rv2=np.zeros_like(rv) 


output_file = open("simulation_output.txt", "w")

# Identify which script is producing the output
script_name = os.path.basename(sys.argv[0])

output_file.write("# ----------------------------------------------------------------------------------\n")
output_file.write(f"# Simulation output generated by: {script_name}\n")
output_file.write(f"# Generation date: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\n")
output_file.write("# ----------------------------------------------------------------------------------\n\n")


# Write parameters
output_file.write("# Simulation Parameters:\n")
for key, val in json_params.items():
    output_file.write(f"#   {key} = {val}\n")

mass_from_h5 = None
try:
    with h5py.File('all_simulation_data.h5', 'r') as f:
        if 'mass' in f.attrs:
            mass_from_h5 = f.attrs['mass']
        elif 'mass' in f:
            mass_from_h5 = f['mass'][()]
except Exception as e:
    print(f" Warning: Could not read mass from all_simulation_data.h5 ({e})")

# Append to parameter section
if mass_from_h5 is not None:
    output_file.write(f"#   mass = {mass_from_h5:.6f} a.u.\n")
else:
    output_file.write(f"#   mass = Not found in all_simulation_data.h5\n")    
    
if field_type == "ac":
    output_file.write(f"#\n# AC field frequency:\n")
    output_file.write(f"#   freq = {freq:.6f} GHz\n")    

output_file.write(f"#\n# Relaxation times (fs):\n")
output_file.write(f"#   tau_min = {tau_min:.6f} fs, tau_max = {tau_max:.6f} fs\n")

# Column header with units
output_file.write("#\n")
output_file.write("# Columns:\n")
output_file.write("#   time [fs], field_target [kV/cm], P_analytical [C/m²], P_density_matrix [C/mm²], Tr[rho] (real)\n")
output_file.write("# ------------------------------------------------------------------------------------------------------\n")





x_expect = np.zeros(nsteps)  # expectation value of position/polarization from analytical solution of SchEq
x_expect_num = np.zeros(nsteps) # expectation value of position/polarization from numer solution of SchEq
x_expect_rho = np.zeros(nsteps) # expectation value of position/polarization from  denisty matrix approach
time_points = np.zeros(nsteps)  # To store corresponding time values
E_vals = np.zeros(nsteps)  # To store electric field values
x_vs_E = np.zeros(nsteps)  # To store x_expect vs E
x_vs_E_rho = np.zeros(nsteps)  # To store x_expect vs E from desnity operator

prev_field_idx = -1  # initialize to detect field changes
t_fs_psi=0.0
for n_t in range(nsteps):

    
    # Your time-stepping code here
    t_fs = n_t * delta_t  # Current time in fs
    t_fs_psi+=delta_t
    time_points[n_t] = t_fs  # Store current time    
    E_current =  field_sequence[n_t]
    
    # Find closest field in loaded data
    field_idx = np.argmin(np.abs(fields - E_current))
    
    E_vals[n_t] = field_sequence[n_t] #fields[field_idx] #E_current  # Store current field value
    
    trace_val = np.real(np.trace(rho_t))
    print(f"Step: {n_t:6d} | Target Field: {E_current:10.4f} kV/cm | Actual Field: {fields[field_idx]:10.4f} kV/cm | Trace ρ: {trace_val:.6f}")
    if trace_val < 0 or trace_val > 1.001:
        raise ValueError(f"ERROR: Trace ρ is out of physical bounds (Trace ρ = {trace_val:.6f}). Terminating simulation.")
    if np.min(np.diag(rho_t)) < -5e-1:  # allow some negative from rounding < 10%
        print(f"WARNING: probability density became negative {np.min(np.diag(rho_t)):.4f}")
   
  

    
    # Rebuild Hamiltonian if field changed
    if field_idx != prev_field_idx:
 #       print("here")
        current_energy = eigval_data[field_idx]
        current_phi = eigvec_data[field_idx]
        
        if intrinsic_dynamics:
           v_E = comp_v_E(E_current)
           H = build_hamiltonian(length, n, v_E, mass)
        
        # Re-project density matrix onto new eigenbasis
        rho_t = project_density_matrix(rho_t, current_phi)       
        
        prev_field_idx = field_idx  # update for next step
        t_fs_psi =  delta_t #reset time for psi 0+delat_t

    # ANALYTICAL solution update
        if use_analytical_model:
            c_n0 = current_phi.conj().T @ psi_t
            c_n0=c_n0/np.linalg.norm(c_n0)
            if np.abs(c_n0[0]) < eps:
                print(
                    f"⚠️ Warning: Ground state population in wavefunction is very small: "
                    f"|c0|^2 = {np.abs(c_n0[0]):.3e}. Analytical model will be disabled.\n"
                    "Recommendation: decrease effective mass. Continuing without analytical modes"
                    )
                use_analytical_model = False
           
    
    if use_analytical_model:
        psi_t = compute_psi_t(current_energy, current_phi, c_n0, t_fs_psi,intrinsic_dynamics)
        psi_t2 = np.abs(psi_t)**2
        x_expect[n_t] = np.sum(x_int * psi_t2)
        x_vs_E[n_t] = x_expect[n_t]
    else:
    # Skip analytical update or set psi_t2 to zero/default
        psi_t2 = np.zeros_like(x_int)
        x_expect[n_t] = 0.0
        x_vs_E[n_t] = 0.0
   
    
    if integrator == "simple":
        # ---- Simple Forward Euler ----
        
        rel = rel_term(rho_t, current_energy, current_phi, Temp, relaxation_model, gamma_rescale)[1:-1, 1:-1]
        if intrinsic_dynamics:
            comm = commutator(H, rho_t[1:-1, 1:-1])
            rv = comm + rel
        else:
            rv = rel
        rho_t[1:-1, 1:-1] += rv * delta_t
        rho_t /= np.trace(rho_t)
    
    elif integrator == "predictor_corrector":
        # ---- Predictor–Corrector ----
        rho_t, rho_0, rv1, rv2 = predictor(rho_t, rv, rv1, rv2, delta_t)
    
        rel = rel_term(rho_t, current_energy, current_phi, Temp, relaxation_model, gamma_rescale)[1:-1, 1:-1]
        if intrinsic_dynamics:
            comm = commutator(H, rho_t[1:-1, 1:-1])
            rv = comm + rel
        else:
            rv = rel
    
        rho_t = corrector(rho_0, rv, rv1, rv2, delta_t)
    
    else:
        raise ValueError(f"Unknown integrator type: {integrator}")


    x_expect_rho[n_t] = np.sum(x_int * np.real(np.diagonal(rho_t)))
    x_vs_E_rho[n_t] = x_expect_rho[n_t]
 

    output_file.write(f"{n_t * delta_t:.6f} {E_current:.6f} {x_expect[n_t] / conv_P2x:.6f} {x_expect_rho[n_t] / conv_P2x  :.6f} {trace_val:.6f}\n")
# plot probability density for all two approaches for comparison    
    if json_params["plot_mode"] == "screen" and n_t % n_out == 0:
        # Convert x_int (position) to polarization:
        polarization = x_int / conv_P2x  # Now in C/m^2 units
        
        plt.plot(polarization, np.diagonal(rho_t), '-b', label='Density Matrix Diagonal')
        plt.plot(polarization, psi_t2, '-r', label='Probability Density (analytical solution)')
        # plt.plot(polarization, np.abs(psi)**2, '-g', label='Probability Density (numerical solution)')
    
        plt.title(f'Step {n_t}: Density Matrix vs Wavefunction (t = {n_t*delta_t:.2f} fs)')
        plt.xlabel('Polarization Expectation (C/m$^2$)')
        plt.ylabel('Probability Density')
        plt.legend()
        plt.show()


        
        
# Convert position expectation values to polarization
P_expect = x_expect / conv_P2x
P_expect_rho = x_expect_rho / conv_P2x
P_vs_E = x_vs_E / conv_P2x
P_vs_E_rho = x_vs_E_rho / conv_P2x
        

# Plot time-dependent polarization expectation value
plt.figure(figsize=(10, 6))
plt.plot(time_points, P_expect, 'r-', linewidth=2, label='Analytical Solution')
plt.plot(time_points, P_expect_rho, 'b-', linewidth=2, label='Density Matrix')
plt.xlabel('Time (fs)', fontsize=12)
plt.ylabel('<P>(t) (C/cm²)', fontsize=12)
plt.title(f'Polarization Expectation Value vs Time at T = {Temp:.1f} K', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("polarization_vs_time.png", dpi=300)   
if json_params["plot_mode"] == "screen":
    plt.show()
plt.close()



if field_type != "dc":
# Plot polarization vs electric field
    first_Emax_idx = np.argmax(E_vals >= E_max - 1e-6)  # Tolerance for floating point
    plt.figure(figsize=(12, 6))
    plt.plot(E_vals[first_Emax_idx:], P_vs_E[first_Emax_idx:], 'r-', linewidth=2, label='Analytical Solution')
    plt.plot(E_vals[first_Emax_idx:], P_vs_E_rho[first_Emax_idx:], 'b-', linewidth=2, label='Density Matrix')
    plt.xlabel('Electric Field (kV/cm)', fontsize=12)
    plt.ylabel('Polarization <P> (C/m²)', fontsize=12)
    plt.title(f'Polarization Expectation vs Applied Electric Field at T = {Temp:.1f} K', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig("polarization_vs_field.png", dpi=300)   # always save
    if json_params["plot_mode"] == "screen":
        plt.show()
    plt.close()
  

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Simulation completed in {elapsed_time:.2f} seconds")






