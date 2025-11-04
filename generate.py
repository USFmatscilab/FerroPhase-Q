


"""
===============================================================================
Hysteresis Loop Simulation for Ferroelectric Materials Using 1D Schrödinger Equation
===============================================================================

Description:
------------
This script simulates polarization hysteresis loops in ferroelectric materials
based on a 1D quantum mechanical model. It computes the average polarization ⟨P⟩
as a function of applied electric field E using a statistical mixture of
quantum states. The energy landscape U(P) is taken as input, and the system's
response is computed by solving the time-independent Schrödinger equation under
various field strengths.

Key Features:
-------------
- Uses statistical density matrix formalism to compute ⟨P⟩(E)
- Solves the 1D Schrödinger equation with field-dependent potential
- Computes thermal occupation of quantum states
- Generates and saves hysteresis loops
- Saves all data (fields, energies, wavefunctions, polarization) to HDF5
- Supports interactive analysis of eigenstates for any field value

Input:
------
- DataMagnitude.txt : Text file containing polarization vs potential energy
                      (U(P) in arbitrary units)

Output:
-------
- hysteresis_loop.png : Plot of ⟨P⟩ vs E (hysteresis loop)
- P_vs_E.txt          : Tabulated polarization vs field data
- all_simulation_data.h5 : HDF5 file with all simulation results
- eigenstates_field_XXX.png : Potential and first eigenstates for chosen field

Usage:
------
1. Make sure 'DataMagnitude.txt' is present in the same directory.
2. Run the script to compute and save the hysteresis loop data.
3. Use the interactive prompt to analyze eigenstates for selected field values.




===============================================================================


"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import h5py
import os
import json


# Constants and global definitions
k_B = 8.617333262145e-5  # eV/K
au2e = 1.66e4 / 9.109
hbar = 6.582119569

# read input file
with open('PARAMS.json', 'r') as f:
    params = json.load(f)

E_min = params['Emin']
E_max = params['Emax']
stepE = params['stepE']
n_states_default = params.get('n_states_to_plot', 3)  # Default to 3 if not specified
Temp = params.get('Temp', 300.0)
z_star = params['z_star']
Vol = params['Vol']
mass = params['mass']
num_sol = params['num_sol']
n = params['n']

conv_P2x = Vol * 0.01 / (z_star * 1.6) # converting polarization into nm


data = np.loadtxt("DataMagnitude.txt")
x = conv_P2x * data[:, 0]
length = x[-1] - x[0]
v = data[:, 1] - np.min(data[:, 1]) # shifts by the minimum



x_int = np.linspace(-length / 2, length / 2, n + 1)
tck = interpolate.splrep(x, v)
v_int = interpolate.splev(x_int, tck)
np.savetxt('interpolated_potential.txt',
           np.column_stack((x_int, v_int)),
           header='x(nm) V(x) (eV)',
           fmt='%.10f')

def comp_v_E(E):
    return v_int - z_star * E * 1.e-04 * x_int

def build_hamiltonian(length, n, v, mass):
    hbar = 1.054571596
    echarge = 1.602176462
    baremass = 9.10938188
    const = hbar ** 2 / (baremass * echarge) / mass
    deltax = length / n
    const /= deltax ** 2
    diagonal = v[1:-1] + const
    offdiagonal = const / 2.0
    diagonals = [diagonal, np.full(n - 2, -offdiagonal), np.full(n - 2, -offdiagonal)]
    Hmatrix = diags(diagonals, [0, -1, 1], format='csr')
    return Hmatrix

def solve(Hmatrix, num_sol):
    eigvals, eigvecs = eigsh(Hmatrix, k=num_sol, which='SA')
    sorted_indices = np.argsort(eigvals)
    e = eigvals[sorted_indices]
    phi = eigvecs[:, sorted_indices]
    phi = np.vstack([np.zeros(num_sol), phi, np.zeros(num_sol)])
    return e, phi

def comp_X_av(e_tmp, phi_tmp):
    Z = 0.
    X_av = 0.
    e_tmp=e_tmp-e_tmp[0]
    for i in range(num_sol):
        prob = np.exp(-e_tmp[i] / (k_B * Temp))
        Z += prob
        X_av += prob * np.sum(x_int * phi_tmp[:, i]**2)
    return X_av / Z

def save_all_data(field_values, all_energies, all_phi, X_E):
    with h5py.File('all_simulation_data.h5', 'w') as f:
        # Save simulation-wide parameters at the root level
        f.attrs['mass'] = mass
        f.attrs['n'] = n
        f.attrs['num_sol'] = num_sol
        f.attrs['Star'] = z_star
        f.attrs['Vol'] = Vol

        f.create_dataset('field_values', data=np.array(field_values, dtype='f8'))
        f.create_dataset('polarization_data', data=np.array(X_E, dtype='f8'))
        eigen_group = f.create_group('eigen_data')
        for i, E in enumerate(field_values):
            field_group = eigen_group.create_group(f'field_{i:03d}')
            field_group.create_dataset('field_value', data=np.float64(E))
            field_group.create_dataset('energies', data=np.array(all_energies[i], dtype='f8'))
            field_group.create_dataset('eigenvectors', data=np.array(all_phi[i], dtype='f8'))
            field_group.create_dataset('potential', data=np.array(comp_v_E(E), dtype='f8'))
        f.create_dataset('z_star', data=np.float64(z_star))
        f.create_dataset('volume', data=np.float64(Vol))
        f.create_dataset('temperature', data=np.float64(Temp))
        f.create_dataset('num_solutions', data=np.int32(num_sol))


def analyze_field_data_by_value(field_value_input, n_states_to_plot=3):
    """Analyze data for a specific field value (in kV/cm) and plot n_states_to_plot eigenfunctions"""
    try:
        with h5py.File('all_simulation_data.h5', 'r') as f:
            # Find closest matching field
            available_fields = list(f['eigen_data'].keys())
            all_field_values = [f[f'eigen_data/{field_name}/field_value'][()] for field_name in available_fields]
            closest_idx = np.argmin(np.abs(np.array(all_field_values) - field_value_input))
            field_name = available_fields[closest_idx]
            E = all_field_values[closest_idx]
            field_group = f[f'eigen_data/{field_name}']

            energies = field_group['energies'][:]
            phi = field_group['eigenvectors'][:]
            potential = field_group['potential'][:]

            print(f"\nAnalysis for field value: {E:.3f} kV/cm (closest to input {field_value_input:.3f} kV/cm)")
            print("First 10 energy levels (eV):")
            print(energies[:10])

            # Calculate thermal populations
            populations = np.exp(-(energies-energies[0]) / (k_B * Temp))
            populations /= np.sum(populations)

            print("\nTop 5 contributing states:")
            top5 = np.argsort(populations)[-5:][::-1]
            for i in top5:
                print(f"State {i}: Energy = {energies[i]:.4f} eV, Population = {populations[i]:.4f}")

            # Plot potential and wavefunctions using twin y-axes
            fig, ax1 = plt.subplots(figsize=(10, 6))

            ax1.set_xlabel('Position (nm)', fontsize=12)
            ax1.set_ylabel('Potential Energy (eV)', color='k', fontsize=12)
            ax1.plot(x_int, potential - np.min(potential), 'k-', label='Potential', linewidth=2)
            ax1.tick_params(axis='y', labelcolor='k')

            ax2 = ax1.twinx()
            ax2.set_ylabel('Wavefunction Intensity (arb. units)', color='b', fontsize=12)
            for i in range(min(n_states_to_plot, phi.shape[1])):
                ax2.plot(x_int, 5 * phi[:, i]**2 + energies[i], label=f'State {i}')
            ax2.tick_params(axis='y', labelcolor='b')

            fig.suptitle(f"Potential and Eigenstates at E = {E:.2f} kV/cm", fontsize=14)
            fig.tight_layout()
            fig.legend(loc='upper right')
            plt.grid(True)
            plt.savefig(f'eigenstates_field_{E:.2f}kVcm.png', dpi=150)
            plt.show()

            # Print all eigenvalues in eV, 2pi THz, and THz
            print("\nAll eigenvalues (index, eV, 2pi THz, THz):")
            ev_to_hz = 2.417989e14  # 1 eV = 2.417989e14 Hz
            for i, val in enumerate(energies):
                freq_hz = val * ev_to_hz
                freq_2pi_thz = freq_hz / 1e12         # Hz to THz (2pi included only in label)
                freq_thz = freq_hz / (2 * np.pi * 1e12)  # Hz to THz (no 2pi)
                print(f"{i:3d}  {val:10.6f}  {freq_2pi_thz:10.6f}  {freq_thz:10.6f}")

            return energies, phi, potential

    except Exception as e:
        print(f"Error analyzing field {field_value_input} kV/cm: {str(e)}")
        return None, None, None



    
field_values = np.linspace(E_min, E_max, stepE)

X_E = np.zeros((stepE, 2))
all_energies = np.zeros((stepE, num_sol))
all_phi = np.zeros((stepE, n + 1, num_sol))
highest_ke = 0.0  # track highest kinetic energy found (in eV)

for i, E in enumerate(field_values):
    v_E = comp_v_E(E)
    H = build_hamiltonian(length, n, v_E, mass)
    energy, phi = solve(H, num_sol)
    
    # --- Boltzmann weights at high T ---
    boltz_weights = np.exp(-(energy-energy[0]) / (k_B * 2000.0))  # ~infinite-T limit
    Z = np.sum(boltz_weights)
    probs = boltz_weights / Z

    
    # --- Identify significant states (cumulative prob ≥ 0.99) ---
    sorted_indices = np.argsort(probs)[::-1]
    cumulative_prob = 0.0
    selected_states = []
    for idx in sorted_indices:
        cumulative_prob += probs[idx]
        selected_states.append(idx)
        if cumulative_prob >= 0.99:
            break

    # --- Build kinetic Hamiltonian only (v_E = 0) ---
    H_kin = build_hamiltonian(length, n, np.zeros_like(v_E), mass)

    # --- Compute KE expectation for selected states ---
    for idx in selected_states:
        psi = phi[1:-1, idx]
        ke_val = np.real(psi.conj().T @ (H_kin @ psi))
        if ke_val > highest_ke:
            highest_ke = ke_val
#            print(highest_ke)

   
    
 # --- Polarization expectation ---   
    X_av = comp_X_av(energy, phi)

    X_E[i, :] = (E, X_av / conv_P2x)
    all_energies[i, :] = energy
    all_phi[i, :, :] = phi

    print(f'Progress: {100*(i+1)/stepE:.1f}% - E = {E:.1f} kV/cm, P = {100*X_av/conv_P2x:.2f} μC/cm²')

# same as building H constants, perhaps need better implementaion    
hbar = 1.054571596
echarge = 1.602176462
baremass = 9.10938188
deltax=length/n

lambda_dbr = 2*np.pi*hbar/np.sqrt(2.*baremass*mass*highest_ke*echarge)
# Number of grid points per wavelength
points_per_lambda = lambda_dbr / deltax

min_points_required = 10  # Suggested minimum number of points per wavelength
if points_per_lambda < min_points_required:
    print(f" WARNING: Spatial resolution may be insufficient for current mass!")
    print(f"   de Broglie wavelength = {lambda_dbr:.4f} nm")
    print(f"   Grid spacing dx = {deltax:.4f} nm")
    print(f"   Points per wavelength = {points_per_lambda:.2f} (minimum recommended = {min_points_required})")
    print(f"   Your current values are: mass = {mass:.4f}, n = {n}, dx = {deltax:.4f} nm, λ = {lambda_dbr:.4f} nm")
else:
    print(f"✅ Sampling is sufficient: {points_per_lambda:.2f} points per shortest wavelength")


    

save_all_data(field_values, all_energies, all_phi, X_E)

plt.figure(figsize=(10, 6))
plt.plot(X_E[:, 0], 100.0 * X_E[:, 1], 'b-', linewidth=2)
plt.xlabel('Electric Field (kV/cm)', fontsize=14)
plt.ylabel('Polarization (μC/cm²)', fontsize=14)
plt.title(f'Hysteresis Loop at {Temp} K', fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('hysteresis_loop.png', dpi=300)
plt.show()


np.savetxt('P_vs_E.txt', X_E,
header='Field(kV/cm) Polarization(microC/cm2)',
fmt='%.6f')

analyze_field_data_by_value(0.0)  # Default analysis for zero field
  


print("\nData analysis mode")
print("Enter electric field values in kV/cm to analyze (or 'q' to quit)")
while True:
  user_input = input("Field value (kV/cm) or 'q': ").strip()
  if user_input.lower() == 'q':
      break
  try:
      field_val = float(user_input)
      n_states_input = input("How many eigenfunctions to plot? (default 3): ").strip()
      if n_states_input == '':
          n_states = 3
      else:
          n_states = int(n_states_input)
      analyze_field_data_by_value(field_val, n_states)
  except ValueError:
      print("Please enter a valid numeric field value or 'q' to quit.")

