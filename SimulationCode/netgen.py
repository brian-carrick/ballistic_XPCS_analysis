#!/use/local/bin/env python
# -*- coding: utf-8 -*-
##
##------------------------------------------------------
## Code for Generation of Network
## Follows algorithm published by 
## AA Gusev, Macromolecules, 2019, 52, 9, 3244-3251
##
## Author: Akash Arora, Postdoc, Olsen Group, ChemE, MIT
##------------------------------------------------------
import math
import numpy as np
from numpy import linalg as LA
import random
import ioLAMMPS

#random.seed(a=786564)


def sortSecond(val): 
    return val[1] #return the element at position 1 ie- second element

def bondlengths(n_chains, chains, links, Lx, Ly, Lz): #!!!!#seems like this is not used anywhere, the bondlength function from relax.py is used 
 
    dist = np.zeros((n_chains,4))#n_chains X 4 array
    dist[:,0:3] = chains 
    dist[:,3] = -1
        
    for i in range (0, n_chains):
        if(chains[i,2] !=-1):
      
          link_1 = chains[i,1]-1
          link_2 = chains[i,2]-1
          lk = links[link_1,:] - links[link_2,:]
          
          lk[0] = lk[0] - int(round(lk[0]/Lx))*Lx # why is this rounding off being done?
          lk[1] = lk[1] - int(round(lk[1]/Ly))*Ly
          lk[2] = lk[2] - int(round(lk[2]/Lz))*Lz
                
          dist[i,3] = LA.norm(lk)

    return dist


def find_neighbours(link_1, n_links, links, conn, cnt_length, mean_length, Lx, Ly, Lz):# does this find the neighbours which are empty or any neightbour?

    neigh = []
    lk    = np.zeros((n_links,3), dtype = float)
    
    lk[:,0] = links[:,0] - links[link_1,0] # x- coordinate of all cross linkers wrt the one under consideration   
    lk[:,1] = links[:,1] - links[link_1,1]    # y- coordinate of all cross linkers wrt the one under consideration
    lk[:,2] = links[:,2] - links[link_1,2]    # z- coordinate of all cross linkers wrt the one under consideration

    for i in range (0, n_links):
       
        lk[i,0] = lk[i,0] - int(round(lk[i,0]/Lx))*Lx #why is this being done???
        lk[i,1] = lk[i,1] - int(round(lk[i,1]/Ly))*Ly# this means that if the coordinates of the link exceeds the simulation box,
                                                        #then  consider the excess
        lk[i,2] = lk[i,2] - int(round(lk[i,2]/Lz))*Lz #is this like a periodic bc? but why is that done?
                
         
    dist = LA.norm(lk, axis=1)#LA is linalg, axis defines the axis for norm. axis=1 means 1 norm-
                                #ie. the distances of crosslinkers from the link under consideration
    dist_list = np.where(dist[:] <= cnt_length)[0] #only consider the crosslinks which are within the contour distance from the crosslink under consideration
                                                #this gives us the linkers to which a connection can be made, ie. the linkers which can count as neighbours
    # Checking what all linkers are empty 
    for i in range(0, np.size(dist_list)):
        lnk = dist_list[i]
        a = np.where(conn[lnk,:] == -1)[0]
        if(np.size(a) > 0): # if there is at least one crosslinker which is empty
           p = math.exp(-1.5*dist[lnk]**2/mean_length) # why is this being done?- this comes from gaussian distribution in 3D
           neigh.append([lnk, dist[lnk], p]) # why is this being appended?

    neigh.sort(key = sortSecond) #sort according to the 1 element(2nd element) of each array element present in neigh, ie.- sort according to dist[lnk]
                                                                                                                      #ie. -sort according to distances

    return neigh


def find_second_link(neigh):
    
    nbr  = np.zeros(len(neigh), dtype = int)
    prob = np.zeros(len(neigh), dtype = float)
    cum_prob = np.zeros(len(neigh)+1, dtype = float)

    for i in range(0, len(neigh)):
        nbr[i]  = neigh[i][0]#linel number, i dont think this is being used anywhere
        prob[i] = neigh[i][2]#probability calculated by the formula p = math.exp(-1.5*dist[lnk]**2/mean_length), as done in the function above

 
    prob = prob/np.sum(prob)#normalize
    cum_prob[0] = 0
    for i in range(0, len(prob)):
        cum_prob[i+1] = cum_prob[i] + prob[i] #calculate cumulative probability for each neighbour link
 
    rnd_num = random.uniform(0,1)
    a1 = np.where(cum_prob <= rnd_num)[0] #similar to kmc algorithm

    second_link = neigh[a1[len(a1)-1]][0] #keep track of the link number in the sorted neigh array

    return second_link



def generate_network(prob, func, N, L, l0, n_chains, n_links): #network is generated, no value returned.
                                                               #all important output is stored in network.txt file

    print('----------------------------')
    print('-----Generating Network-----')
    print('----------------------------')
        
    cnt_length = N*l0 #contour length , according to Mao model, contour length can be allowed to extend beyond Nb.
                     
    mean_length = N*l0*l0 # why is this value chosen? (this is l0 and not the number 10!!)
    Lx = L; Ly = L; Lz = L
    
    # Array sizes and initializations
    links   = np.zeros((n_links,3), dtype = float)  #keeps track of the coordinates of each cross linker
                                                    # these are the cross linkers assuminng that the system is fully reacted
                                                    # ie. there are no dangling chains at the ends (as opposed to the finite network case)
                                                    
    chains  = np.full((n_chains,3), -1, dtype = int)#return a n_chains X 3 array with all values as -1 filled
                                                    # the 3 indices are: 0- n_i (U_i can also be used in place of n_i), 1-linker1, 2-linker2 (index 1 nad 2 are the connectivities)
    conn    = np.full((n_links,func), -1, dtype = int)#each link is associated with func number of connections 
    
    
    for i in range (0, n_links): # loop over all the crosslinks
        x = random.uniform(0, Lx)# select x,y,z coordinates randomly
        y = random.uniform(0, Ly)
        z = random.uniform(0, Lz)
        links[i,:] = [x, y, z] #assign position to the crosslink randomly
    
    for i in range (0, n_chains):
                #***************IMPORTANT for mixed bond *********************    
        chains[i,0] = N   #position 0 is n_i (or U_i can also be used here alternatively.
                          #in case of mixed bonds- this will be different for the different kinds of bonds
        found = 0 #i think this variable is just for tracking, not anything particularly used in the computations
        while (found == 0): #link not found yet
              link_1 = random.randint(0, n_links-1) #select a link randomly
              a1 = np.where(conn[link_1,:] == -1)[0]#find the first position in the connections of link_1 have open positions/ unreacted dangling chains
              if(np.size(a1)>0): #if the above query returns some value(s) (which means there are some open positions), then set found to 1
                 found = 1 #this means that there is a link with unreacted positions left
    
        rnd_num = random.uniform(0,1) #uniform random number between 0 and 1
        if (prob > rnd_num): # prob is an input parameter, if this condition is satisfied, then assign the chains to the selected link
           chains[i,1] = link_1 #assign link_1 as the first link of the chain
           a1 = np.where(conn[link_1,:] == -1)[0]#find the first position in the connections of link_1 have open positions/ unreacted dangling chains
           conn[link_1, a1[0]] = n_links + 1  #denotes temporary filling, if not filled/replaced later, this will denote a dangling end
                                              #this works because there no actual link at n_links+1 position, so this would mean temp filling/dangling end
           # This ensures that link occ is filled
           # temporary filing by n_links + 1
           # later, it will be replaced by link 2 if it needs to be connected. 
           # Otherwise, it will be as it is, indicating that this site 
           # contains a dangling end.

           rnd_num = random.uniform(0,1)#uniform random number between 0 and 1
           neigh = find_neighbours(link_1, n_links, links, conn, cnt_length, mean_length, Lx, Ly, Lz)#find neighbours of cross linker 1- link_1
           if (prob > rnd_num and len(neigh) > 0): #means unfilled neighbours of link_1 is/are found,
                                                   # and condition is favourable for linking (according to probabiltiy)  
              link_2 = find_second_link(neigh)#chose the second link from among the nighbours
                                              #how to determine whether this second link has at least one open position??
    
              chains[i,2] = link_2 #assign the second link
              a2 = np.where(conn[link_2,:] == -1)[0] #check which is the first position of the second linker which is empty-
                                                      # why is this not being done before the step where the link is assigned to the chian?
                                                      #if it turns out that link 2 has no open positions, then there is no way to reverse that assignment? right???
              conn[link_1, a1[0]] = link_2 #assign link 2 to open position of link_1
              conn[link_2, a2[0]] = link_1 #assign link_1 to open position of link_2
    ##what does a1[0] and a2[0] mean? because a1 and a2 are anyways ()... [0], then again a1[0] won't make sense right??
    
    a1 = np.where(chains[:,1]==-1)[0] #find which chains have link 1 disconnected
    a2 = np.where(chains[:,2]==-1)[0] #find which chains have link2 disconnected
    n_free = np.size(a1) #a1 is the number of free chains. because for a chain- only if link 1 is connected can link 2 be connected
                        # so, if link 1 itself is disconnected, then link 2 is already discinnected. Hence, it is a free chain.
    n_dang = np.size(a2)-np.size(a1) #those chains which have link2 disconnected inlcude two types-
                                     #free(link 2 disconnected because link 1 is not connected in the first place)
                                     #and dangling ends (link 1 connected, but not link 2). 
    n_connected = n_chains - np.size(a2) # all chains with link 2 disconnected are either free chains or dangling ends. so, the rest are all connected
    
    n_loops = 0
    all_loop_atoms = []# includes all loops- including dumbells and primary loops
    for i in range (0, n_chains):
        if(chains[i,2] !=-1): #ie- link 2 not disconnected, means the chain is connected
           if(chains[i,1] == chains[i,2]): #the main concept of a primary loop is that link1 and link2 of a chain are both the same link(cross linker)
             n_loops = n_loops + 1
             all_loop_atoms.append(chains[i,1]) #add the link number for that primary loop into the array

    all_loop_atoms.sort() #this is counting only primary loops right?- No
    all_loop_atoms.append(-1) #why is this appending -1 being done?
 
    i = 0
    loop_atoms = [] 
    dumbbell_atoms = []
    while (i < len(all_loop_atoms)-1): #count within the primary loops
        if (all_loop_atoms[i] == all_loop_atoms[i+1]): #if the common linker of primary loop of chain 1 is the same as that of chain 2-
                                                        #then this would mean that the structure is dumbell
                                                        #BUT THIS IS NOT GENERAL. THIS WILL WORK ONLY FOR TETRAFUNCTIONAL NETWORKS(func=4)
                                                        #BUT THIS WILL NOT WORK FOR TRIFUCNTIONAL NETWORKS!!
           dumbbell_atoms.append(all_loop_atoms[i])
           i = i + 2 #because i+1 would be the same thing as i, hence move forward by i+1
        else:
           loop_atoms.append(all_loop_atoms[i]) #this means that it is only primary loop and not a dumbell
           i = i + 1
           
    #write to file      
    f1 = open('all_loops','w')
    for i in range(0, len(all_loop_atoms)):
        f1.write('%5i  %5i\n'%(i, all_loop_atoms[i])) 
    f1.close()

    f1 = open('primary_loops','w')
    for i in range(0, len(loop_atoms)):
        f1.write('%5i  %5i\n'%(i, loop_atoms[i])) 
    f1.close()

    f1 = open('dumbbells','w')
    for i in range(0, len(dumbbell_atoms)): 
        f1.write('%5i  %5i\n'%(i, dumbbell_atoms[i])) 
    f1.close()

    print('PRIMARY LOOP FRACTION', len(loop_atoms)/n_chains)
    
    # new_chain are the chains excluding all the loops, dangling, and sol (free) , kind of ideal chains which are there in the network
    new_chains = np.zeros((n_chains - n_loops - n_dang - n_free, 3), dtype = int)
    cnt = 0 
    for i in range (0, n_chains):
        if(chains[i,2]!=-1): #means link 2 is not disconnected- means link 2 is connected, means that it cannot be free or danlging end
           if(chains[i,1]!=chains[i,2]): #if it is not a loop (either primary or dumbell)
             new_chains[cnt,:] = chains[i,:] #add chain to the new chain array
             cnt = cnt + 1

    # Upto now, chains (or new chains) contain link with indices from 0 to 1199
    # But, for writing LAMMPS file as well as for calculating distances
    # chains are modified to have linkers with indices starting from 1 to 1200
    new_chains[:,1] = new_chains[:,1] + 1    
    new_chains[:,2] = new_chains[:,2] + 1    

    #write network characteristics to file
    print('----------------------------')
    print('--Network Characteristics---')
    print('----------------------------')
    file1=open("net_char","w")
    file1.write("Box Length = %6.4f\n" %(L))
    file1.write("Contour Length = %i\n" %(N))
    file1.write("No. of strands = %i\n" %(n_chains))
    file1.write("No. of linkers = %i\n" %(n_links))
    file1.write("Extent of reaction = %6.2f\n" %(prob))
    file1.write("Free chains = %i\n" %(n_free))
    file1.write("Dangling ends = %i\n" %(n_dang))
    file1.write("Connected chains = %i\n" %(n_connected))
    file1.write('No. of primary loops = %i\n' %(n_loops))
    file1.close()


##    dist = bondlengths(n_chains, chains, links, Lx, Ly, Lz)    
##    file1=open("chain_conn","w")
##    file2=open("links_conn","w")
##    for i in range(0, n_chains): 
##        file1.write("%4i  %4i  %4i  %6.4f\n" % (chains[i,0], chains[i,1], chains[i,2], dist[i,3]))
##    
##    for i in range(0, n_links):
##        file2.write("%4i  %4i  %4i\n" % (conn[i,0], conn[i,1], conn[i,2]))
##    
##    file1.close()
##    file2.close()


#write to file-0 the network obtained after generation
    xlo = 0; ylo = 0; zlo = 0
    xhi = L; yhi = L; zhi = L
    atom_types = 2; bond_types = 1; n_break = 0
    mass= np.full(atom_types, 1, dtype = float)
#    ioLAMMPS.writeLAMMPS(xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, n_loops, n_dang, n_free, n_break, links, chains, atom_types, bond_types, mass, loop_atoms)
#    print(len(new_chains[:,0]))
    ioLAMMPS.writeLAMMPSafternetgen('network.txt', xlo, xhi, ylo, yhi, zlo, zhi, links, new_chains, atom_types, bond_types, mass, loop_atoms)

##   [xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, 
##           atom_types, bond_types, mass] = netgen.generate_network(prob, func, N, L, l0, n_chains, n_links)
##
##    n_atoms = n_links
##    n_bonds = len(new_chains[:,0])  
    
#    return xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, links, new_chains, atoms_types, bond_types, mass
#    return links, new_chains, conn, n_loops, n_dang, n_free, n_break
    
#    return loop_atoms
