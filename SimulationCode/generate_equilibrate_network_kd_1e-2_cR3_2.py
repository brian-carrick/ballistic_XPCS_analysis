#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
#######################################
#                                     #
#-- Fracture Simulation of Networks --#
#------  Author: Akash Arora  --------#
#-- Postdoc, Olsen Group, ChemE, MIT--#
#                                     #
#######################################

 Overall Framework (Steps):
     1. Generate a Network following the algorithm published
        by AA Gusev, Macromolecules, 2019, 52, 9, 3244-3251
        if gen_net = 0, then it reads topology from user-supplied 
        network.txt file present in this folder
     
     2. Force relaxtion of network using Fast Inertial Relaxation Engine (FIRE) 
        to obtain the equilibrium positions of crosslinks (min-energy configuration)

     3. Compute Properties: Energy, Gamma (prestretch), and 
        Stress (all 6 componenets) 
     
     4. Deform the network (tensile) in desired direction by 
        supplying lambda_x, lambda_y, lambda_z

     5. Break bonds using kintetic theory of fracture (force-activated KMC). 
        Currently implemented algorithm is inspired by 
        Termonia et al., Macromolecules, 1985, 18, 2246

     6. Repeat steps 2-5 until the given extension (lam_total) is achieved Or    
        stress decreases below a certain (user-specified) value 
        indicating that material is completey fractured 
        (currently, only lam_total functionality is implemented)
"""


#-------------------------------------#
#       Simulation Parameters         #
#-------------------------------------#
# Parameters to vary
N    = 22 # Number of Kuhn segments
func = 5 #functionality of the network
tol = 0.1
tau_be = 0.1  # tau-Gillespie time step (s)
bond_exchange_steps = 100000  # Number of tau-Gillespie steps per iteration
k_A = 1  # Forward association rate constant


L = int(80/2) #length of simulation box
k_D = 1e-2  # Reverse dissociation rate constant

equilibration_network = 'equilibration_network_kd_1e-2_cR3_2.txt'


# Constants & calculations
rho  = 3 #??
l0   = 1
prob = 1
n_chains  = int(10000/8) #No. of chains in the simulation box (n/L3 defines the initial chain density)
n_links   = int(2*n_chains/func) #Stochiometric mixture, #links means the individual segments
print(n_links, n_chains)



Nb = N
K  = 1.0
r0 = 0.0
U0  = 1
tau = 1
erate = 1
lam_max = 11
max_itr = 10000
write_itr = 100     #output FIRE parameters ever x iterations
wrt_step = 1        #Output file every x iterations





import time
import math
import random
import netgen
import ioLAMMPS
import bondexchange
import matplotlib
import numpy as np
from relax import Optimizer
from numpy import linalg as LA
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
##import timeit

start = time.time()
print("hello")


random.seed(a=None, version=2)# Intialize random number generator
print('First random number of this seed: %d'%(random.randint(0, 10000))) 
# This is just to check whether different jobs have different seeds

netgen_flag = 1# if 0: read network from file, if 1: generate network using Gusev’s method
if(netgen_flag==0): #read from file
    
    [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
        atom_types, bond_types, mass, loop_atoms] = ioLAMMPS.readLAMMPS(equilibration_network, N, 0)


    
    # ===== CHECK ATOM VALENCE =====
    print(f'\n--- Checking atom valence ---')
    atom_bonds = [[] for _ in range(n_atoms)]
    
    for bond in bonds:
        link_1 = int(bond[2]) - 1
        link_2 = int(bond[3]) - 1
        atom_bonds[link_1].append(link_2)
        atom_bonds[link_2].append(link_1)
    
    for i in range(n_atoms):
        if len(atom_bonds[i]) > func:
            print(f'  Atom {i}: {len(atom_bonds[i])} bonds (OVERFLOW!)')
    
    print('--------------------------')
    
    print('xlo, xhi',xlo, xhi) 
    print('ylo, yhi',ylo, yhi) 
    print('zlo, zhi',zlo, zhi) 
    print('n_atoms', n_atoms) 
    print('n_bonds', n_bonds) 
    print('atom_types =', atom_types) 
    print('bond_types =', bond_types) 
    print('mass =', mass) 
    print('loop atoms =', len(loop_atoms))
    print('inter-atom bonds =', n_bonds - len(loop_atoms))
    print('--------------------------')

    # ===== DIAGNOSTIC =====
    print(f'\n--- RELOAD DIAGNOSTIC ---')
    print(f'bonds shape: {bonds.shape}')
    self_bonds = sum(1 for b in bonds if int(b[2]) == int(b[3]))
    inter_bonds = len(bonds) - self_bonds
    print(f'Self-bonds: {self_bonds}')
    print(f'Inter-atom bonds: {inter_bonds}')
    print(f'Total: {len(bonds)}')
    
    print('--------------------------')
   
elif(netgen_flag==1):#generate using Gussev's algorithm
   netgen.generate_network(prob, func, N, L, l0, n_chains, n_links)
   
   # ===== READ LOOP INFORMATION FROM FILES CREATED BY NETGEN =====
   loop_atoms_from_file = []
   try:
       with open('all_loops', 'r') as f:
           for line in f:
               parts = line.strip().split()
               if len(parts) >= 2:
                   atom_idx = int(parts[1])
                   if atom_idx != -1:
                       loop_atoms_from_file.append(atom_idx)
       print(f'Read {len(loop_atoms_from_file)} loop atoms from all_loops file')
   except FileNotFoundError:
       print('WARNING: all_loops file not found')
   
   # Now read the network
   [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
           atom_types, bond_types, mass, loop_atoms_dummy] = ioLAMMPS.readLAMMPS("network.txt", N, 0)
   
   # Use the loop atoms from the file, not from readLAMMPS
   loop_atoms = loop_atoms_from_file
   
   print(f'\n=== BEFORE writeLAMMPS ===')
   print(f'loop_atoms length: {len(loop_atoms)}')
   print(f'loop_atoms (first 10): {loop_atoms[:10] if len(loop_atoms) > 0 else "EMPTY"}')
   
   # Re-write to include self-bonds and mark loop atoms as type 2
   ioLAMMPS.writeLAMMPS("network.txt", xlo, xhi, ylo, yhi, zlo, zhi,
                        atoms, bonds, atom_types, bond_types, mass, loop_atoms)
   
   print(f'=== AFTER writeLAMMPS ===')
   
   # ===== CRITICAL: READ BACK THE UPDATED FILE =====
   [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
           atom_types, bond_types, mass, loop_atoms] = ioLAMMPS.readLAMMPS("network.txt", N, 0)
   
   print(f'\n=== AFTER SECOND readLAMMPS ===')
   print(f'loop_atoms length: {len(loop_atoms)}')
   print(f'loop_atoms (first 10): {loop_atoms[:10] if len(loop_atoms) > 0 else "EMPTY"}')

else:
   print('Invalid network generation flag')


fstr=open('stress','w')
fstr.write('#Lx, Ly, Lz, lambda, FE, deltaFE, st[0], st[1], st[2], st[3], st[4], st[5]\n') 

flen=open('strand_lengths','w')
flen.write('#lambda, ave(R), max(R)\n') 

fkmc=open('KMC_stats','w')
fkmc.write('#lambda, init bonds, final bonds\n')

fbe=open('bond_exchange_stats_kd_1e-2_cR3_2','w')  # NEW FILE FOR BOND EXCHANGE STATS
fbe.write('#step, alpha_inter, alpha_loop, alpha_dissoc, n_inter, n_loop, n_dissoc, n_bonds_total, n_bonds_bonded, F_total\n')

#-------------------------------------#
#   Open Bond Exchange Output File    #
#-------------------------------------#

fbe_detailed = open('bond_exchange_detailed_kd_1e-2_cR3_2.txt', 'w')
fbe_detailed.write('# Bond Exchange Detailed Statistics\n')
fbe_detailed.write('# Time[s]  n_inter  n_intra  n_dissoc  F_free/F_total  n_loops/n_total_chains  Energy  avg(r)/Nb  max(r)/Nb\n')
fbe_detailed.flush()

time_elapsed = 0.0  # Track cumulative time in seconds





print(n_chains, lam_max)
c=float(n_chains)/(L**3.0)
b=1
dim_conc=c*(b**3)*(N**1.5)
print('dim_conc',dim_conc)

#-------------------------------------#
#       First Force Relaxation        #
#-------------------------------------#

print('--------------------------')   
print('----First Force relaxation-------')   
print('--------------------------') 
mymin = Optimizer(atoms, bonds, xlo, xhi, ylo, yhi, zlo, zhi, K, r0, N, 'Mao')

print(f'\n--- AFTER Optimizer creation ---')
print(f'bonds variable length: {len(bonds)}')
print(f'mymin.bonds length: {len(mymin.bonds)}')
print(f'Are they the same object? {bonds is mymin.bonds}')

# Count self-bonds
self_bonds_in_bonds = sum(1 for b in bonds if int(b[2]) == int(b[3]))
self_bonds_in_mymin = sum(1 for b in mymin.bonds if int(b[2]) == int(b[3]))
print(f'Self-bonds in bonds: {self_bonds_in_bonds}')
print(f'Self-bonds in mymin.bonds: {self_bonds_in_mymin}')
print(f'loop_atoms length: {len(loop_atoms)}')

[e, Gamma] = mymin.fire_iterate(tol, max_itr, write_itr, 'log.txt')

print(f'\n--- AFTER FIRE ---')
print(f'bonds variable length: {len(bonds)}')
print(f'mymin.bonds length: {len(mymin.bonds)}')

self_bonds_in_bonds = sum(1 for b in bonds if int(b[2]) == int(b[3]))
self_bonds_in_mymin = sum(1 for b in mymin.bonds if int(b[2]) == int(b[3]))
print(f'Self-bonds in bonds: {self_bonds_in_bonds}')
print(f'Self-bonds in mymin.bonds: {self_bonds_in_mymin}')
print(f'loop_atoms length: {len(loop_atoms)}')

ioLAMMPS.writeLAMMPS(equilibration_network, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, mymin.zlo, 
                                  mymin.zhi, mymin.atoms, mymin.bonds, atom_types, bond_types, mass, loop_atoms)

dist = mymin.bondlengths()#calculate bond lengths
Lx0 = mymin.xhi-mymin.xlo# length in x direction
BE0 = e #what is this?
[pxx, pyy, pzz, pxy, pyz, pzx] = mymin.compute_pressure()
fstr.write('%7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f\n' \
                          %(mymin.xhi-mymin.xlo, mymin.yhi-mymin.ylo, mymin.zhi-mymin.zlo, 
                           (mymin.xhi-mymin.xlo)/Lx0, e, e-BE0, pxx, pyy, pzz, pxy, pyz, pzx)) # write to file 'stress'
fstr.flush() #clears the internal buffer of the file
 
flen.write('%7.4f  %7.4f  %7.4f\n'%((mymin.xhi-mymin.xlo)/Lx0, np.mean(dist[:,3])/N, np.max(dist[:,3])/N)) # write to file 'strand lengths()'
flen.flush()#clears the internal buffer of the file

fkmc.write('%7.4f  %5i  %5i\n'%((mymin.xhi-mymin.xlo)/Lx0, n_bonds, n_bonds)) #write to file kmc_stats
fkmc.flush()#clears the internal buffer of the file

ioLAMMPS.writeLAMMPS(equilibration_network, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, mymin.zlo, 
                                  mymin.zhi, mymin.atoms, mymin.bonds, atom_types, bond_types, mass, loop_atoms) #write to file in LAMMPS format

print('Gamma',Gamma/(N*b^2))
##print



#-------------------------------------#
#      Bond Exchange Simulation       #
#-------------------------------------#

print('--------------------------')   
print('----Bond Exchange Dynamics-------')   
print('--------------------------') 

Lx = mymin.xhi - mymin.xlo
Ly = mymin.yhi - mymin.ylo
Lz = mymin.zhi - mymin.zlo

# IMPORTANT: Pass atoms array so BondExchange can identify loop atoms by type 2
be_sim = bondexchange.BondExchange(mymin, bonds, N, b, k_A, k_D, func, 
                                 tau_be, Lx, Ly, Lz)


# ========== DIAGNOSTIC: CHECK INITIAL CONNECTION MATRIX ==========
print(f'\n--- Initial Network State Diagnostic ---')

f_i_init, F_total_init = be_sim.calculate_free_stickers()
F_max = be_sim.n_links * be_sim.func

print(f'Free stickers by atom (first 20): {f_i_init[:20]}')
print(f'Atoms with 0 free stickers: {np.sum(f_i_init == 0)}')
print(f'Atoms with all free stickers: {np.sum(f_i_init == be_sim.func)}')
print(f'Total free stickers: {F_total_init} / {F_max}')
print(f'Total bonded sites: {F_max - F_total_init}')
print(f'Expected number of bonds: {(F_max - F_total_init) // 2}')
print(f'Actual bonds in array: {len(bonds)}')

# Check each bond is in connection matrix
bonded_pairs = be_sim._find_bonded_pairs_from_conn()
print(f'Bonds found in connection matrix: {len(bonded_pairs)}')

# Verify self-bonds are in bonds array
self_bonds_in_array = 0
for bond in bonds:
    if int(bond[2]) == int(bond[3]):
        self_bonds_in_array += 1

print(f'Self-bonds in bonds array: {self_bonds_in_array}')
print(f'Loop atoms: {len(be_sim.loop_atoms_raw)}')
print(f'Expected self-bonds to equal loop atoms: {self_bonds_in_array} should equal {len(be_sim.loop_atoms_raw)}')
print('--- End Diagnostic ---\n')

for step in range(bond_exchange_steps):
    
    
    # ========== PERFORM TAU-GILLESPIE STEP ==========
    step_data = be_sim.tau_gillespie_step(step_num=step, verbose=False)
    
    
    # ========== FORCE EQUILIBRATION (FIRE) ==========
    try:
        [e, Gamma] = mymin.fire_iterate(tol, max_itr, write_itr, 'log.txt')
        [pxx, pyy, pzz, pxy, pyz, pzx] = mymin.compute_pressure()
    except Exception as fire_error:
        print(f'FIRE error at step {step}: {fire_error}')
        raise
    
    # ========== GATHER STATISTICS FOR OUTPUT ==========
    
    f_i_current, F_total_current = be_sim.calculate_free_stickers()
    F_max = be_sim.n_links * be_sim.func
    free_sticker_fraction = F_total_current / F_max
    
    n_loops_current = len(be_sim.loop_atoms_raw)
    n_total_chains = len(mymin.bonds) + n_loops_current
    loop_fraction = n_loops_current / n_chains
    
    dist = mymin.bondlengths()
    avg_r_over_Nb = np.mean(dist[:, 3]) / N if len(dist) > 0 else 0.0
    max_r_over_Nb = np.max(dist[:, 3]) / N if len(dist) > 0 else 0.0
    
    time_elapsed += tau_be
    
    fbe_detailed.write('%12.6e  %7i  %7i  %7i  %16.8f  %16.8f  %12.6e  %12.6f  %12.6f\n' \
                       %(time_elapsed, step_data['n_inter'], step_data['n_loop'], 
                         step_data['n_dissoc'], free_sticker_fraction, loop_fraction, 
                         e, avg_r_over_Nb, max_r_over_Nb))
    fbe_detailed.flush()
    
    # ========== WRITE OUTPUT ==========
    fstr.write('%7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %5i  %5i  %5i\n' \
                                     %(Lx, Ly, Lz, 
                                      Lx/Lx0, e, e-BE0, pxx, pyy, pzz, pxy, pyz, pzx, 
                                      step_data['n_inter'], step_data['n_loop'], step_data['n_dissoc']))
    fstr.flush()
    
    fbe.write('%7i  %.6e  %.6e  %.6e  %5i  %5i  %5i  %5i  %5i  %5i  %5i\n' \
                %(step, step_data['alpha_inter'], step_data['alpha_loop'], step_data['alpha_dissoc'],
                  step_data['n_inter'], step_data['n_loop'], step_data['n_dissoc'], 
                  step_data['n_bonds'], step_data['n_bonds_bonded'], step_data['F_total'], 
                  step_data['total_loop_chains']))
    fbe.flush()
    
    dist = mymin.bondlengths()
    flen.write('%7.4f  %7.4f  %7.4f\n'%(Lx/Lx0, np.mean(dist[:,3])/N, np.max(dist[:,3])/N))
    flen.flush()
    
    # Write restart file periodically
    if((step+1) % wrt_step == 0): 
        #filename = 'restart_network_%d.txt' %(step+1)
        filename = equilibration_network
        ioLAMMPS.writeLAMMPS(filename, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, 
                            mymin.zlo, mymin.zhi, mymin.atoms, mymin.bonds, 
                            atom_types, bond_types, mass, be_sim.loop_atoms_raw)
        print(f'\nStep {step+1}/{bond_exchange_steps}: n_bonds={step_data["n_bonds"]}, '
              f'n_loops={step_data["total_loop_chains"]}, free_stickers={F_total_current}/{F_max}, '
              f'alpha_total={step_data["alpha_total"]:.6e}\n')

#---------------------------------#
#   Close Bond Exchange Files     #
#---------------------------------#

fbe_detailed.close()


#---------------------------------#
#     Final Network Properties    #
#---------------------------------#

print('--------------------------')   
print('----Final Network Properties-------')   
print('--------------------------') 

[e, Gamma] = mymin.fire_iterate(tol, max_itr, write_itr, 'log.txt')
[pxx, pyy, pzz, pxy, pyz, pzx] = mymin.compute_pressure()

fstr.write('%7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %5i  %5i  %5i\n' \
                                 %(mymin.xhi-mymin.xlo, mymin.yhi-mymin.ylo, mymin.zhi-mymin.zlo, 
                              (mymin.xhi-mymin.xlo)/Lx0, e, e-BE0, pxx, pyy, pzz, pxy, pyz, pzx, 0, 0, 0))
fstr.flush()

dist = mymin.bondlengths()
flen.write('%7.4f  %7.4f  %7.4f\n'%((mymin.xhi-mymin.xlo)/Lx0, np.mean(dist[:,3])/N, np.max(dist[:,3])/N))
flen.flush()

# Print final statistics
final_stats = be_sim.get_statistics()
print('\n--- Final Bond Exchange Statistics ---')
print(f'Total intermolecular associations: {final_stats["total_associations_inter"]}')
print(f'Total intramolecular associations: {final_stats["total_associations_loop"]}')
print(f'Total dissociations: {final_stats["total_dissociations"]}')
print(f'Net bonds formed: {final_stats["net_bonds_change"]}')
print(f'Final free stickers: {F_total_current}/{F_max}')
print(f'Final loop chains: {len(be_sim.loop_atoms_raw)}')
print('--------------------------------------\n')

filename = equilibration_network
ioLAMMPS.writeLAMMPS(filename, mymin.xlo, mymin.xhi, mymin.ylo, mymin.yhi, mymin.zlo, mymin.zhi,
                                       mymin.atoms, mymin.bonds, atom_types, bond_types, mass, 
                                       be_sim.loop_atoms_raw)

fstr.close()
flen.close()
fbe.close()

end = time.time()
print(f'Total simulation time: {end - start:.2f} seconds')
