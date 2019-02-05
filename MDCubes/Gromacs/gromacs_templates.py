gromacs_minimization = """
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = {nsteps:d}         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = {pbc}     ; Periodic Boundary Conditions in all 3 dimensions
"""


gromacs_nvt_npt = """
title		= NVT-NPT
define		= -DPOSRES	; position restrain the protein
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= {nsteps:d}		; number od md steps
dt		    = {timestep:f}		; md timestep
; Output control
nstxout		= {trajectory_steps:d}		; save coordinates
nstvout		= {trajectory_steps:d}		; save velocities 
nstenergy	= {trajectory_steps:d}		; save energies 
nstlog		= {report_steps:d}		; update log file 
; Bond parameters
continuation	        = no		; Restarting after NVT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= {temperature:f} 	  {temperature}	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = {pcoupl}	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = {pressure:f}		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc		= {pbc}		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= {gen_vel}		; Velocity generation is off 
"""