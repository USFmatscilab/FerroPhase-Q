# FerroPhase-Q
Quantum Ferroelectric Phase Transitions

1. generate.py (input PARAMS.json): generates input data
2. optional analyze.py (input PARAMS_anal.json) analyzes data and makes some plots
3. quantum_v20.py (input PARAMS_dyn.json) run calculations for several models, see explnation above


{
  "Temp": 900.0,
 "field_type": "ac",
 "E_max": 400.0,
 "E_min": -400.0,
 "num_periods": 1.25,
 "nsteps": 100000,
 "delta_t": 0.025,
 "intrinsic_dynamics": false,
 "P_cen": 0.0,
 "delta_P": 0.1,
 "n_out": 500,
 "relaxation_model": "rlx_lindblad_grn",
 "gamma_rescale": 1.0,
 "plot_mode": "screen",
  "integrator": "simple"
}





"Temp": 0.01, # Temp in K
"field_type": "ac", "linear"
"E_max" and "E_min": -400.0 range of fields to simulate
"num_periods": 1.25 # not properly tested
"nsteps": 100000, # number of steps
"delta_t": 0.085, # intergration time step. Use 0.00085 if intrinsic dynamics is turned on
"intrinsic_dynamics": false, # true or false. If false relaxation only will be used
"P_cen": 0.0, # initial condition, centering of the wave packet in polarization units C/m^2
"delta_P": 0.1, initial condition, width of the packet in polarization units C/m^2
"n_out": 500, # how often to plot probability densities
"relaxation_model": "rlx_lindblad_all" # different relaxation models
#"rlx_equil_dens"  relaxation to equilibrium desnity matrix
#"rlx_lindblad_grn" Lindblad relaxation to the ground state only
#"rlx_lindblad_all" Lindblad relaxation to all states
#"rlx_lindblad_grn_hybrid" Lindblad relaxation to the ground state only but no gain in ground state population, only happens through renormalization
"gamma_rescale": 1.0 rescalies damping parameter for rlx_lindblad_all and lx_lindblad_all, suggested use for rlx_lindblad_all is 1
"plot_mode": "screen"         // NEW: "screen" or "none"
 "integrator": "simple"   // options: "simple", "predictor_corrector" simple works better for large time steps
