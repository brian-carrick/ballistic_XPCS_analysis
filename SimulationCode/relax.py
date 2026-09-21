#!/use/local/bin/env python
# -*- coding: utf-8 -*-
##
##-----------------------------------------------------------
## Fast Inertial Relaxation Engine (FIRE) Optimizer
## Ref: Bitzek et al, PRL, 97, 170201 (2006)
##
## Author: Akash Arora, Postdoc, Olsen Group, ChemE, MIT
## Implementation is inspired from LAMMPS and ASE Master Code
##-----------------------------------------------------------

import math
import time
import random
import numpy as np
import scipy.optimize as opt
from numpy import linalg as LA
from scipy.optimize import fsolve


class Optimizer(object): #class 

    def __init__(self, atoms, bonds, xlo, xhi, ylo, yhi, zlo, zhi, K, r0, N, ftype):

        self.atoms = atoms #crosslinks
        self.bonds = bonds #chains
##        print('bonds',bonds)
##        stop
        self.xlo = xlo #lowest, x
        self.xhi = xhi #highest, x
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi
        self.K = K # what is K??
        self.r0 = r0 #initial distance (length of chain)          
        self.N = N #number of Kuh  segments   
        self.ftype = ftype #force type

    def bondlengths(self):#calculate bond lengths- distance between two cross linkers
     
        atoms = self.atoms# crosslinks
        bonds = self.bonds # chains
        Lx = self.xhi - self.xlo #relative to xlo- ie. get the vector of the chain
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(self.atoms[:,0]) # this will always be N+1 right?- no - if atoms means crosslinks then not 
        n_bonds = len(self.bonds[:,0]) #this will always be N right?- no - if onds means chains then not 

        dist = np.zeros((n_bonds,4), dtype=float) #n_bondsX4 array

        for i in range (0, n_bonds):
              # ehre- the variable bonds is the same as the variable chains, and the structure of the variable bonds is as follows:
              #first column is the serial number, second column is all 1 (why is that??) #third column (column number 2) is link 1, fourth column is link 2
              lnk_1 = bonds[i,2]-1 #why -1?
              lnk_2 = bonds[i,3]-1 #why -1?
              delr = atoms[lnk_1,:] - atoms[lnk_2,:]
              
              delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx # why is this rounding off being done??
              delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
              delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                   
              dist[i,0:3] = delr # only 3 elemets- 0:3 means it will take element 0,1,2
              dist[i,3] = LA.norm(delr) #measure the distance(length) of the chain from the three bonds
    
        return dist


    def get_bondforce(self, r):
        """
        Calculate attractive bond force for a Gaussian chain (linear Hooke's law).
        Returns scalar force magnitude (negative = attractive).
        """
        K = self.K
        r0 = self.r0
        Nb = self.N
    
        x = (r - r0) / Nb
        fbkT = 3.0 * x
        fbond = -K * fbkT / r
        
        return fbond

    # def get_force(self):
        
    #     N = self.N
    #     b = 1.0  # Kuhn segment size
    #     atoms = self.atoms
    #     bonds = self.bonds
    #     ftype = self.ftype
    #     Lx = self.xhi - self.xlo
    #     Ly = self.yhi - self.ylo
    #     Lz = self.zhi - self.zlo
    #     n_atoms = len(atoms[:,0])
    #     n_bonds = len(bonds[:,0])
    
    #     e = 0.0 # Energy of chain
    #     Gamma = 0.0 #what is Gamma
    #     f = np.zeros((n_atoms,3), dtype = float) #empty array for force on each atom in all three directions
        
    #     # ==================== GAUSSIAN REPULSION PARAMETERS ====================
    #     epsilon = 45  # Repulsion strength
    #     Nb_sq = N * b**2  # R_e^2 = N*b^2
    #     r_c = 4*b * (2 * N / 3)**1.5
    #     gaussian_prefactor = (3.0 / (2.0 * np.pi * Nb_sq))**1.5  # Prefactor for Gaussian
        
    #     # ==================== BONDED FORCES ====================
    #     for i in range(0, n_bonds):

    #         lnk_1 = bonds[i,2]-1 #why is -1 done??
    #         lnk_2 = bonds[i,3]-1
    #         delr = atoms[lnk_1,:] - atoms[lnk_2,:] #the vector of the chain between link 1 and link 2
            
    #         delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx #why is this done?
    #         delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
    #         delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                
    #         r = LA.norm(delr)# the length of chain

    #         if (r > 0): #if there is any stretching, ie. when r not equal to 0
    #             fbond = self.get_bondforce(r) #get the bond forces
    #             e_stretch = (3/2) * (r-self.r0) * (r-self.r0)/N
    #             e = e + e_stretch 
    #         else: #when r=0
    #             fbond = 0.0
    #             e = e + 0.0
    
    #         Gamma = Gamma + r*r
    
    #         # apply force to each of 2 atoms        
    #         if (lnk_1 < n_atoms): #just to keep a check that the crosslinks which are outside our domain are not considered
    #             f[lnk_1,0] = f[lnk_1,0] + delr[0]*fbond
    #             f[lnk_1,1] = f[lnk_1,1] + delr[1]*fbond
    #             f[lnk_1,2] = f[lnk_1,2] + delr[2]*fbond
        
    #         if (lnk_2 < n_atoms):
    #             f[lnk_2,0] = f[lnk_2,0] - delr[0]*fbond
    #             f[lnk_2,1] = f[lnk_2,1] - delr[1]*fbond
    #             f[lnk_2,2] = f[lnk_2,2] - delr[2]*fbond
        
    #     # ==================== PAIRWISE GAUSSIAN REPULSION (ALL ATOMS) ====================
    #     for i in range(0, n_atoms):
    #         for j in range(i+1, n_atoms):
                
    #             delr = atoms[i,:] - atoms[j,:] #distance vector between atoms i and j
                
    #             delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx #apply periodic boundary conditions
    #             delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
    #             delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                
    #             r = LA.norm(delr) #distance between atoms
                
    #             if (r < r_c and r > 1e-10): #within cutoff and not identical
                    
    #                 # Isotropic Gaussian repulsion: U = epsilon * prefactor * exp(-3/2 * r^2 / (N*b^2))
    #                 arg = 1.5 * r**2 / Nb_sq
    #                 gaussian_envelope = gaussian_prefactor * np.exp(-arg)
                    
    #                 e_repel = epsilon * gaussian_envelope
    #                 e = e + e_repel
                    
    #                 # Repulsive force magnitude: f = epsilon * gaussian * 3 * r / (N*b^2)
    #                 f_repel_magnitude = epsilon * gaussian_envelope * 3.0 / Nb_sq
                    
    #                 # Apply repulsive force to both atoms
    #                 if (i < n_atoms):
    #                     f[i,0] = f[i,0] + f_repel_magnitude * delr[0] / r
    #                     f[i,1] = f[i,1] + f_repel_magnitude * delr[1] / r
    #                     f[i,2] = f[i,2] + f_repel_magnitude * delr[2] / r
                    
    #                 if (j < n_atoms):
    #                     f[j,0] = f[j,0] - f_repel_magnitude * delr[0] / r
    #                     f[j,1] = f[j,1] - f_repel_magnitude * delr[1] / r
    #                     f[j,2] = f[j,2] - f_repel_magnitude * delr[2] / r
        
    #     return f, e, Gamma         

    def get_force(self):
        
        N = self.N
        b = 1.0  # Kuhn segment size
        atoms = self.atoms
        bonds = self.bonds
        ftype = self.ftype
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(atoms[:,0])
        n_bonds = len(bonds[:,0])
    
        e = 0.0 # Energy of chain
        Gamma = 0.0 #what is Gamma
        f = np.zeros((n_atoms,3), dtype = float) #empty array for force on each atom in all three directions
        
        # ==================== GAUSSIAN REPULSION PARAMETERS ====================
        epsilon = 60  # Repulsion strength
        Nb_sq = N * b**2  # R_e^2 = N*b^2
        r_c = 4*b * (2 * N / 3)**1.5
        gaussian_prefactor = (3.0 / (2.0 * np.pi * Nb_sq))**1.5  # Prefactor for Gaussian
        
        # ==================== BONDED FORCES ====================
        for i in range(0, n_bonds):

            lnk_1 = bonds[i,2]-1 #why is -1 done??
            lnk_2 = bonds[i,3]-1
            delr = atoms[lnk_1,:] - atoms[lnk_2,:] #the vector of the chain between link 1 and link 2
            
            delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx #why is this done?
            delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
            delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                
            r = LA.norm(delr)# the length of chain

            if (r > 0): #if there is any stretching, ie. when r not equal to 0
                fbond = self.get_bondforce(r) #get the bond forces
                e_stretch = (3/2) * (r-self.r0) * (r-self.r0)/N
                e = e + e_stretch 
            else: #when r=0
                fbond = 0.0
                e = e + 0.0
    
            Gamma = Gamma + r*r
    
            # apply force to each of 2 atoms        
            if (lnk_1 < n_atoms): #just to keep a check that the crosslinks which are outside our domain are not considered
                f[lnk_1,0] = f[lnk_1,0] + delr[0]*fbond
                f[lnk_1,1] = f[lnk_1,1] + delr[1]*fbond
                f[lnk_1,2] = f[lnk_1,2] + delr[2]*fbond
        
            if (lnk_2 < n_atoms):
                f[lnk_2,0] = f[lnk_2,0] - delr[0]*fbond
                f[lnk_2,1] = f[lnk_2,1] - delr[1]*fbond
                f[lnk_2,2] = f[lnk_2,2] - delr[2]*fbond
        
        # ==================== BUILD CELL LIST ====================
        cell_size = r_c + 0.1  # Slightly larger than cutoff
        n_cells_x = max(1, int(Lx / cell_size))
        n_cells_y = max(1, int(Ly / cell_size))
        n_cells_z = max(1, int(Lz / cell_size))
        
        # Assign atoms to cells
        cells = {}
        for i in range(0, n_atoms):
            # Shift coordinates to origin, then find cell indices with PBC
            x_shifted = atoms[i,0] - self.xlo
            y_shifted = atoms[i,1] - self.ylo
            z_shifted = atoms[i,2] - self.zlo
            
            ix = int(x_shifted / cell_size) % n_cells_x
            iy = int(y_shifted / cell_size) % n_cells_y
            iz = int(z_shifted / cell_size) % n_cells_z
            
            cell_key = (ix, iy, iz)
            
            if cell_key not in cells:
                cells[cell_key] = []
            cells[cell_key].append(i)
        
        # ==================== PAIRWISE GAUSSIAN REPULSION (CELL LIST) ====================
        for cell_key, atom_list in cells.items():
            ix, iy, iz = cell_key
            
            # Check this cell and 26 neighboring cells
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    for dz in [-1, 0, 1]:
                        neighbor_key = ((ix+dx) % n_cells_x, (iy+dy) % n_cells_y, (iz+dz) % n_cells_z)
                        
                        if neighbor_key not in cells:
                            continue
                        
                        neighbor_list = cells[neighbor_key]
                        
                        for i in atom_list:
                            for j in neighbor_list:
                                if i >= j:  # Avoid duplicates and self-interaction
                                    continue
                                
                                delr = atoms[i,:] - atoms[j,:] #distance vector between atoms i and j
                                
                                delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx #apply periodic boundary conditions
                                delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
                                delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                                
                                r_sq = delr[0]**2 + delr[1]**2 + delr[2]**2  #distance squared (avoid sqrt if not needed)
                                
                                if (r_sq < r_c**2 and r_sq > 1e-20):  #within cutoff and not identical
                                    
                                    r = np.sqrt(r_sq)  #distance between atoms
                                    
                                    # Isotropic Gaussian repulsion: U = epsilon * prefactor * exp(-3/2 * r^2 / (N*b^2))
                                    arg = 1.5 * r_sq / Nb_sq
                                    gaussian_envelope = gaussian_prefactor * np.exp(-arg)
                                    
                                    e_repel = epsilon * gaussian_envelope
                                    e = e + e_repel
                                    
                                    # Repulsive force magnitude: f = epsilon * gaussian * 3 * r / (N*b^2)
                                    f_repel_magnitude = epsilon * gaussian_envelope * 3.0 / Nb_sq
                                    
                                    # Apply repulsive force to both atoms
                                    if (i < n_atoms):
                                        f[i,0] = f[i,0] + f_repel_magnitude * delr[0] / r
                                        f[i,1] = f[i,1] + f_repel_magnitude * delr[1] / r
                                        f[i,2] = f[i,2] + f_repel_magnitude * delr[2] / r
                                    
                                    if (j < n_atoms):
                                        f[j,0] = f[j,0] - f_repel_magnitude * delr[0] / r
                                        f[j,1] = f[j,1] - f_repel_magnitude * delr[1] / r
                                        f[j,2] = f[j,2] - f_repel_magnitude * delr[2] / r
        
        return f, e, Gamma
   
  
 
    def fire_iterate(self, ftol, maxiter, write_itr, logfilename):
      # based on FIRE algorithm, algorithm implemented exactly as done in paper
        tstart = time.time()

        ## Optimization parameters:
        eps_energy = 1.0e-8 #don't know what this is, but this is not used anywhere else, so doesn't matter!!
        delaystep = 5 #N_min
        dt_grow = 1.1 # f_inc
        dt_shrink = 0.5#f_dec
        alpha0 = 0.1 #alpha_start
        alpha_shrink = 0.99 #f_alpha
        tmax = 10.0#factor by which delta_t_max is more than delta_t_MD
        maxmove = 0.1 #maximum movement possible, right?? what does maxmove denote??
        last_negative = 0 #previous iteration step when P was negative

        dt = 0.005 #delta_t_MD
        dtmax = dt*tmax #delta_t_max
        alpha = alpha0
        last_negative = 0       
 
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        n_atoms = len(self.atoms[:,0])
        n_bonds = len(self.bonds[:,0])
        v = np.zeros((n_atoms,3), dtype = float) #array to keep track of velocity  of each crosslink(atom)

        n_bonds = len(self.bonds)
        dist = np.zeros((n_bonds,4), dtype=float)# first 3 columns are the distances in each direction, and last is the actual distance (morm)
# but this is not used, it is reassigned after next few lines

        [f,e,Gamma] = self.get_force()
        dist = self.bondlengths()

 
        fmaxitr = np.max(np.max(np.absolute(f))) #find maximum force in a particular direction
        # why is the force in a particular direction taken? what significance does it have?
        fnormitr = math.sqrt(np.vdot(f,f)) #vdot- the matrix is flattened # what is the significance of this??
##        logfile = open(logfilename,'w') 
##        logfile.write('FIRE: iter  Energy  fmax  fnorm  avg(r)/Nb  max(r)/Nb\n')
##        logfile.write('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f\n' %
##                              ('FIRE', 0, e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
##        logfile.flush()
        print('FIRE: iter  Energy  fmax  fnorm  avg(r)/Nb  max(r)/Nb')
        print('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f' %
                              ('FIRE', 0, e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))

        for itr in range (0, maxiter):
         
          vdotf = np.vdot(v,f)  #power P(t)=v.F
          if (vdotf > 0.0): #positive power
             vdotv = np.vdot(v,v)
             fdotf = np.vdot(f,f) 
             scale1 = 1.0 - alpha
             if (fdotf == 0.0): scale2 = 0.0 #if F is zero, then anyways there would not be any contribution from the force term
             else: scale2 = alpha * math.sqrt(vdotv/fdotf) #divide by sqrt(fdotf ) because f/sqrt(fdotf) is equal to unit vector in F direction
             v = scale1*v + scale2*f
              
             if (itr - last_negative > delaystep): # if number f steps since P was last negative is larger than N_min
                 dt = min(dt*dt_grow,dtmax)#increase time step
                 alpha = alpha*alpha_shrink #decrease alpha- means that rely less on present v and take more bias of F_(unitvecotr)*mod(v)
      
          else: #P is negative
             last_negative = itr
             dt = dt*dt_shrink #decrease delta t
             alpha = alpha0 #reset alpha_0
             v[:] = v[:]*0.0 #freeze system
      
          v = v + dt*f  #BIG DOUBT!! # why is this done??? This was not there in the algorithm?
          dr = dt*v #delta r- movement, vector- has 3 elements for each row
          normdr = np.sqrt(np.vdot(dr, dr))
          if (normdr > maxmove):
              dr = maxmove * dr / normdr # the maximum movement possible is max move, so cap it at maxmove, with the direction vector intact

          self.atoms = self.atoms + dr
          for i in range(0, n_atoms):
              self.atoms[i,0] = self.atoms[i,0] - math.floor((self.atoms[i,0]-self.xlo)/Lx)*Lx #like what was done earlier. but WHY??
              self.atoms[i,1] = self.atoms[i,1] - math.floor((self.atoms[i,1]-self.ylo)/Ly)*Ly
              self.atoms[i,2] = self.atoms[i,2] - math.floor((self.atoms[i,2]-self.zlo)/Lz)*Lz
          
  
          [f,e,Gamma] = self.get_force() #this is already done earlier. no need to put it here. taking up computation unnecessarily, right??
          fmaxitr = np.max(np.max(np.absolute(f)))
          fnormitr = math.sqrt(np.vdot(f,f))


          if((itr+1)%write_itr==0): # write_itr- number of iterations upto which we want to write to file
             dist = self.bondlengths()
##             logfile.write('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f\n' %
##                                  ('FIRE', itr+1, e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
##             logfile.flush()

             # Print on screen
             print('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f' %
                               ('FIRE', itr+1,  e, fmaxitr, fnormitr,  np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
        
   
          # Checking for convergence
          if (fnormitr < ftol): #convergence acheived, wite to file
             dist = self.bondlengths()
             tend = time.time() #end time
##             logfile.write('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f\n' %
##                                  ('FIRE', itr+1, e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
##             logfile.flush()
             print('%s: %5d  %9.6f  %9.6f  %9.6f  %9.4f  %9.4f' %
                               ('FIRE', itr+1,  e, fmaxitr, fnormitr, np.mean(dist[:,3])/self.N, np.max(dist[:,3])/self.N))
             print('Iterations converged, Time taken: %7.4f' %(tend-tstart))
             break
          elif (itr == maxiter-1):
             print('Maximum iterations reached')
     

##        logfile.close() 
        
        return e, Gamma       
                

    def compute_pressure(self):

        K = self.K #what is K??
        r0 = self.r0
        ftype = self.ftype
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        atoms = self.atoms
        bonds = self.bonds
        n_atoms = len(atoms[:,0])
        n_bonds = len(bonds[:,0])
       
        pxx = pyy = pzz = pxy = pyz = pzx = 0.0
        sigma = np.zeros((n_atoms,6), dtype=float)#stress tensor, with pij=pji
        inv_volume = 1.0/(Lx*Ly*Lz)#1/V
        for i in range(0, n_bonds):

            lnk_1 = bonds[i,2]-1#why is -1 done??
            lnk_2 = bonds[i,3]-1
            delr = atoms[lnk_1,:] - atoms[lnk_2,:]
            
            delr[0] = delr[0] - int(round(delr[0]/Lx))*Lx #periodic boundary condition
            delr[1] = delr[1] - int(round(delr[1]/Ly))*Ly
            delr[2] = delr[2] - int(round(delr[2]/Lz))*Lz
                 
            r = LA.norm(delr)
            if (r > 0.0): #if stretched
               fbond = self.get_bondforce(r) 
            else: fbond = 0.0 #if not stretched
            
            # apply pressure to each of the 2 atoms   
            # And for each of the 6 components     
            if (lnk_1 < n_atoms):
               sigma[lnk_1,0] = sigma[lnk_1,0] + 0.5*delr[0]*delr[0]*fbond # xx # why is stress~ force*r^2??? #why factor of 0.5??
               sigma[lnk_1,1] = sigma[lnk_1,1] + 0.5*delr[1]*delr[1]*fbond #yy
               sigma[lnk_1,2] = sigma[lnk_1,2] + 0.5*delr[2]*delr[2]*fbond #zz
               sigma[lnk_1,3] = sigma[lnk_1,3] + 0.5*delr[0]*delr[1]*fbond #xy
               sigma[lnk_1,4] = sigma[lnk_1,4] + 0.5*delr[1]*delr[2]*fbond #yz
               sigma[lnk_1,5] = sigma[lnk_1,5] + 0.5*delr[2]*delr[0]*fbond #zx
        
            if (lnk_2 < n_atoms):
               sigma[lnk_2,0] = sigma[lnk_2,0] + 0.5*delr[0]*delr[0]*fbond
               sigma[lnk_2,1] = sigma[lnk_2,1] + 0.5*delr[1]*delr[1]*fbond
               sigma[lnk_2,2] = sigma[lnk_2,2] + 0.5*delr[2]*delr[2]*fbond
               sigma[lnk_2,3] = sigma[lnk_2,3] + 0.5*delr[0]*delr[1]*fbond
               sigma[lnk_2,4] = sigma[lnk_2,4] + 0.5*delr[1]*delr[2]*fbond
               sigma[lnk_2,5] = sigma[lnk_2,5] + 0.5*delr[2]*delr[0]*fbond


        pxx = np.sum(sigma[:,0])*inv_volume #pressure=stress/volume # but the units don't match right?? pressure unit becomes force/length
        pyy = np.sum(sigma[:,1])*inv_volume
        pzz = np.sum(sigma[:,2])*inv_volume
        pxy = np.sum(sigma[:,3])*inv_volume
        pyz = np.sum(sigma[:,4])*inv_volume
        pzx = np.sum(sigma[:,5])*inv_volume

        return pxx, pyy, pzz, pxy, pyz, pzx


    def change_box(self, scale_x, scale_y, scale_z):

        xlo = self.xlo
        xhi = self.xhi
        ylo = self.ylo
        yhi = self.yhi
        zlo = self.zlo
        zhi = self.zhi
        atoms = self.atoms
        bonds = self.bonds
        n_atoms = len(atoms[:,0])
        n_bonds = len(bonds[:,0])

        xmid = (xlo+xhi)/2  
        ymid = (ylo+yhi)/2  
        zmid = (zlo+zhi)/2  

        new_xlo = xmid + scale_x*(xlo-xmid)# scaling is done this way, just by definition, new_xlo<xmid (since xlo<xmid)
        new_ylo = ymid + scale_y*(ylo-ymid)
        new_zlo = zmid + scale_z*(zlo-zmid)

        new_xhi = xmid + scale_x*(xhi-xmid)# new_xhi>xmid (since xhi>xmid)
        new_yhi = ymid + scale_y*(yhi-ymid)
        new_zhi = zmid + scale_z*(zhi-zmid)
        
        newLx = new_xhi - new_xlo
        newLy = new_yhi - new_ylo
        newLz = new_zhi - new_zlo
        for i in range(0, n_atoms):            
            atoms[i,0] = xmid + scale_x*(atoms[i,0]-xmid) #reassign the coordinates of the atoms (ie. cross-linkers)
            atoms[i,1] = ymid + scale_y*(atoms[i,1]-ymid)
            atoms[i,2] = zmid + scale_z*(atoms[i,2]-zmid)

        self.atoms = atoms #assign the modified values (after scaling the box accordingly
        self.xlo = new_xlo
        self.xhi = new_xhi
        self.ylo = new_ylo
        self.yhi = new_yhi
        self.zlo = new_zlo
        self.zhi = new_zhi
