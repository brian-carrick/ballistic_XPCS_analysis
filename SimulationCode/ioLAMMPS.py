#!/use/local/bin/env python
# -*- coding: utf-8 -*-
##-----------------------------------------------------------
##
## Python script to write configuration in LAMMPS format
## Author: Akash Arora, Postdoc, Olsen Group, ChemE, MIT
##-----------------------------------------------------------

import math
import numpy as np

def readLAMMPS(filename, N, vflag):
    """
    Read LAMMPS file and extract network information
    Properly handles loop atoms identified by atom type 2
    """
    
    f1 = open(filename, "r")

    line1 = f1.readline()
    line2 = f1.readline()

    line3 = f1.readline()
    line3 = line3.strip()
    n_atoms = int(line3.split(" ")[0])  # Renamed from n_links to n_atoms for clarity
 
    line4 = f1.readline()
    line4 = line4.strip()
    atom_types = int(line4.split(" ")[0])

    line5 = f1.readline()
    line5 = line5.strip()
    n_bonds = int(line5.split(" ")[0])

    line6 = f1.readline()
    line6 = line6.strip()
    bond_types = int(line6.split(" ")[0])

    atoms_unsort = np.zeros((n_atoms, 4))  # id, type, x, y, z (will store id, type, x, y, z)
    atoms = np.zeros((n_atoms, 3), dtype=float)  # Just x, y, z coordinates
    bonds = np.full((n_bonds, 4), -1, dtype=int)  # bond_type, index, link_1, link_2
    mass = np.zeros(atom_types, dtype=float)

    line7 = f1.readline()
    line8 = f1.readline()
    line8 = line8.strip()
    xlo = float(line8.split(" ")[0])
    xhi = float(line8.split(" ")[1])

    line9 = f1.readline()
    line9 = line9.strip()
    ylo = float(line9.split(" ")[0])
    yhi = float(line9.split(" ")[1])

    line10 = f1.readline()
    line10 = line10.strip()
    zlo = float(line10.split(" ")[0])
    zhi = float(line10.split(" ")[1])

    for i in range(0, 3):
        f1.readline()
    
    for i in range(0, atom_types):
        line = f1.readline()
        line = line.strip()
        mass[i] = float(line.split(" ")[1])

    f1.close()

    # ==================== READ ATOMS ====================
    # Load atom data: id, type, x, y, z
    atoms_unsort = np.genfromtxt(filename, usecols=(0, 2, 3, 4, 5), skip_header=18, max_rows=n_atoms)
    
    # atoms_unsort columns: [id, type, x, y, z]
    atom_types_array = np.zeros(n_atoms, dtype=int)  # Store atom types for later
    
    for i in range(0, n_atoms):
        atom_id = int(atoms_unsort[i, 0])
        atom_type = int(atoms_unsort[i, 1])  # Now this is actually column 2 (atom_type)
        atoms[atom_id - 1, :] = atoms_unsort[i, 2:5]  # Store coordinates (x, y, z)
        atom_types_array[atom_id - 1] = atom_type  # Store atom type

    # ==================== READ BONDS ====================
    # Load bond data: bond_type, link_1, link_2
    bonds_unsort = np.genfromtxt(filename, usecols=(0, 1, 2, 3), skip_header=18 + n_atoms + 3, max_rows=n_bonds)
    
    # bonds_unsort columns: [bond_id, bond_type, link_1, link_2]
    for i in range(0, n_bonds):
        bond_id = int(bonds_unsort[i, 0])
        bonds[bond_id - 1, 0] = int(bonds_unsort[i, 1])  # bond_type
        bonds[bond_id - 1, 1] = bond_id  # bond_index
        bonds[bond_id - 1, 2] = int(bonds_unsort[i, 2])  # link_1 (1-indexed)
        bonds[bond_id - 1, 3] = int(bonds_unsort[i, 3])  # link_2 (1-indexed)

    # ==================== EXTRACT LOOP ATOMS FROM ATOM TYPES ====================
    # Loop atoms are identified by atom type 2
    loop_atoms = []
    for i in range(0, n_atoms):
        if atom_types_array[i] == 2:
            loop_atoms.append(i)  # 0-indexed
    
    print(f'\nReadLAMMPS summary:')
    print(f'  - {n_atoms} atoms')
    print(f'  - {n_bonds} bonds')
    print(f'  - {len(loop_atoms)} loop atoms (type 2)')
    print(f'  - Box: x=[{xlo:.4f}, {xhi:.4f}], y=[{ylo:.4f}, {yhi:.4f}], z=[{zlo:.4f}, {zhi:.4f}]')

    # ==================== RECONSTRUCT SELF-BONDS FROM LOOP ATOMS =====
    # If loop atoms exist but self-bonds aren't in bonds array, add them
    self_bonds_in_array = sum(1 for b in bonds if int(b[2]) == int(b[3]))
    
    if len(loop_atoms) > 0 and self_bonds_in_array == 0:
        print(f'Reconstructing {len(loop_atoms)} self-bonds from loop_atoms...')
        
        for loop_atom in loop_atoms:
            loop_atom_1indexed = int(loop_atom) + 1  # Convert to 1-indexed
            bond_id = len(bonds) + 1
            self_bond = np.array([[1, bond_id, loop_atom_1indexed, loop_atom_1indexed]], dtype=int)
            bonds = np.vstack([bonds, self_bond])
        
        n_bonds = len(bonds)
        print(f'Bonds after reconstruction: {n_bonds}')
    
    
    return xlo, xhi, ylo, yhi, zlo, zhi, n_atoms, n_bonds, atoms, bonds, atom_types, bond_types, mass, loop_atoms


def writeLAMMPS(filename, xlo, xhi, ylo, yhi, zlo, zhi, atoms, bonds, atom_types, bond_types, mass, loop_atoms):
    """
    Write network to LAMMPS format file
    The bonds array ALREADY contains self-bonds for loop atoms
    """
    
    n_atoms = len(atoms[:,0])
    
    # Convert loop_atoms to set for quick lookup (0-indexed)
    loop_atoms_set = set(int(atom) for atom in loop_atoms)
    
    n_bonds_total = len(bonds)  # Don't add loop_atoms - they're already in bonds!

    print(f'writeLAMMPS: {n_bonds_total} total bonds (including self-bonds)')
    
    func = 5
    max_bonds = (n_atoms * func) // 2
    print(f'writeLAMMPS: Maximum allowed bonds: {max_bonds}')
    
    if n_bonds_total > max_bonds:
        print(f'\n⚠ ERROR: Total bonds ({n_bonds_total}) exceeds maximum ({max_bonds})!')
        raise ValueError(f'Bond limit exceeded: {n_bonds_total} > {max_bonds}')
    
    # LAMMPS input file begins here
    file1=open(filename,"w")
    file1.write("# Gels-Network -- Configuration\n")
    file1.write("\n")
    file1.write("%i %s\n" % (n_atoms,"atoms"))
    file1.write("%i %s\n" % (atom_types,"atom types"))
    file1.write("%i %s\n" % (n_bonds_total,"bonds"))  # Don't add loop_atoms here
    file1.write("%i %s\n" % (bond_types,"bond types"))
    file1.write("\n")
    file1.write("%7.4f %7.4f %s %s\n" % (xlo,xhi,"xlo","xhi"))
    file1.write("%7.4f %7.4f %s %s\n" % (ylo,yhi,"ylo","yhi"))
    file1.write("%7.4f %7.4f %s %s\n" % (zlo,zhi,"zlo","zhi"))
    
    file1.write("\n")
    file1.write("Masses\n")
    file1.write("\n")
    for i in range(0, atom_types):
        file1.write("%i %6.4f\n" % (i+1,mass[i]))
    
    file1.write("\n")
    file1.write("Atoms\n")
    file1.write("\n")
    
    for i in range(0, n_atoms):
        # Check if this atom is a loop atom (type 2)
        if i in loop_atoms_set:
            atom_type = 2
        else:
            atom_type = 1
        
        file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(i+1, 1, atom_type, atoms[i,0], atoms[i,1], atoms[i,2]))
    
    # Bond writing starts here
    file1.write("\n")
    file1.write("Bonds\n")
    file1.write("\n")
    
    # Write ALL bonds (which already includes self-bonds)
    # Don't add them again!
    for i in range(0, len(bonds)):
        file1.write("{:4} {:4} {:4} {:4}\n".format(i+1, 1, bonds[i,2], bonds[i,3]))
    
    file1.close()



#def writeLAMMPSafternetgen(filename, xlo, xhi, ylo, yhi, zlo, zhi, n_links, n_chains, n_loops, n_dang, n_free, n_break, links, chains, atom_types, bond_types, mass):
def writeLAMMPSafternetgen(filename, xlo, xhi, ylo, yhi, zlo, zhi, links, chains, atom_types, bond_types, mass, loop_atoms):

#   n_atoms  = n_links
#   n_bonds  = n_chains - n_loops - n_dang - n_free -n_break
 
   n_atoms = len(links[:,0])  
   n_bonds = len(chains[:,0])  
   n_loops = len(loop_atoms)

   print(n_atoms, n_bonds, n_loops)
 
   # LAMMPS input file begins here
   file1=open(filename,"w")
   file1.write("# Gels-Network -- Initial Configuration\n")
   file1.write("\n")
   file1.write("%i %s\n" % (n_atoms,"atoms"))
   file1.write("%i %s\n" % (atom_types,"atom types"))
   file1.write("%i %s\n" % (n_bonds,"bonds"))
   file1.write("%i %s\n" % (bond_types,"bond types"))
   file1.write("\n")
   file1.write("%7.4f %7.4f %s %s\n" % (xlo,xhi,"xlo","xhi"))
   file1.write("%7.4f %7.4f %s %s\n" % (ylo,yhi,"ylo","yhi"))
   file1.write("%7.4f %7.4f %s %s\n" % (zlo,zhi,"zlo","zhi"))
   
   file1.write("\n")
   file1.write("Masses\n")
   file1.write("\n")
   for i in range(0, atom_types):
       file1.write("%i %6.4f\n" % (i+1,mass[i]))
   
   file1.write("\n")
   file1.write("Atoms\n")
   file1.write("\n")
   cnt = 0
   cnt_loops = 0
   for i in range(0, n_atoms):
       if(cnt_loops < n_loops and cnt == loop_atoms[cnt_loops]): #not yet reached the count for all loops and counter indicates loop, keep adding
          cnt = cnt + 1
          cnt_loops = cnt_loops + 1 #loop detected and written to file
          file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(cnt,1,2,links[i,0],links[i,1],links[i,2])) #why is there a column of 1s and 2s here?
                                                                                                         #maybe because 2 types of atoms -but why?
                                                                                                         #or is it that index=2 if loop detected, and 1 if not
       else:                                                                                             
          cnt = cnt + 1#here cnt_loops is not incremented because there is no loop detected
          file1.write("{:2} {:2} {:2} {:7} {:7} {:7}\n".format(cnt,1,1,links[i,0],links[i,1],links[i,2]))
   
  
   print(cnt, cnt_loops) 
   # Bond writing starts here
   file1.write("\n")
   file1.write("Bonds\n")
   file1.write("\n")
   cnt = 0
#   for i in range(0, n_chains):
#         if(chains[i,2]!=-1):
#            if(chains[i,1] != chains[i,2]):
#               cnt = cnt + 1
#               file1.write("{:4} {:4} {:4} {:4}\n".format(cnt,1,chains[i,1], chains[i,2]))
   
   for i in range(0, n_bonds):
       cnt = cnt + 1
       file1.write("{:4} {:4} {:4} {:4}\n".format(cnt,1,chains[i,1], chains[i,2])) #first column is the serial number, second column is all 1 (why is that??)
                                                                                   #maybe because only one type of bond??
                                                                                   #third column is link 1, fourth column is link 2
   
   file1.close()


   # ==================== WRITE LOOP FILES ====================
    
    # Analyze loop_atoms to separate primary loops and dumbbells
   loop_atom_counts = {}
   for atom in loop_atoms:
        atom_idx = int(atom)
        loop_atom_counts[atom_idx] = loop_atom_counts.get(atom_idx, 0) + 1
    
   # Separate into primary loops and dumbbells
   primary_loops = []
   dumbbell_atoms = []
    
   for atom_idx, count in loop_atom_counts.items():
       if count == 1:
           primary_loops.append(atom_idx)
       elif count == 2:
           dumbbell_atoms.append(atom_idx)
    
   # Create all_loops by combining primary loops and dumbbells (duplicates for dumbbells)
   all_loop_atoms = primary_loops.copy()
   for db_atom in dumbbell_atoms:
       all_loop_atoms.append(db_atom)
   all_loop_atoms.sort()
    
   # Write all_loops file
   f1 = open('all_loops', 'w')
   for i in range(0, len(all_loop_atoms)):
       f1.write('%5i  %5i\n' % (i, all_loop_atoms[i]))
   f1.close()
    
   # Write primary_loops file
   f1 = open('primary_loops', 'w')
   for i in range(0, len(primary_loops)):
       f1.write('%5i  %5i\n' % (i, primary_loops[i]))
   f1.close()
    
   # Write dumbbells file
   f1 = open('dumbbells', 'w')
   for i in range(0, len(dumbbell_atoms)):
       f1.write('%5i  %5i\n' % (i, dumbbell_atoms[i]))
   f1.close()
    
   print(f'Loop files updated:')
   print(f'  - all_loops: {len(all_loop_atoms)} total loop chains')
   print(f'  - primary_loops: {len(primary_loops)} atoms with single loop')
   print(f'  - dumbbells: {len(dumbbell_atoms)} atoms with dual loops')


def extract_loop_information(atoms_array, bonds_array):
    """
    Extract loop and connectivity information from LAMMPS arrays
    
    Parameters:
    -----------
    atoms_array : ndarray
        Atoms array from LAMMPS file (atom_id, type, ...)
    bonds_array : ndarray
        Bonds array from LAMMPS file (bond_id, type, link_1, link_2, ...)
    
    Returns:
    --------
    all_loop_atoms : list
        All atoms involved in loops (with duplicates for dumbbells)
    primary_loops : list
        Atoms with exactly one self-loop
    dumbbell_atoms : list
        Atoms with exactly two self-loops
    """
    
    # Find loop atoms - these are marked with type 2 in the Atoms section
    loop_atoms_set = set()
    for i, atom in enumerate(atoms_array):
        atom_type = int(atom[1])  # Second column is atom type
        if atom_type == 2:
            # Convert to 0-indexed (i is already 0-indexed in the array)
            loop_atoms_set.add(i)
    
    # Count how many self-loops each atom has
    self_loop_counts = {}
    for i in loop_atoms_set:
        self_loop_counts[i] = 0
    
    # Scan bonds to find self-loops (bonds where link_1 == link_2)
    for bond in bonds_array:
        link_1 = int(bond[2]) - 1  # Convert to 0-indexed
        link_2 = int(bond[3]) - 1  # Convert to 0-indexed
        
        # Check if this is a self-loop
        if link_1 == link_2 and link_1 in loop_atoms_set:
            self_loop_counts[link_1] += 1
    
    # Separate into primary loops and dumbbells
    primary_loops = []
    dumbbell_atoms = []
    
    for atom_idx in sorted(self_loop_counts.keys()):
        count = self_loop_counts[atom_idx]
        if count == 1:
            primary_loops.append(atom_idx)
        elif count == 2:
            dumbbell_atoms.append(atom_idx)
    
    # Create all_loop_atoms (includes duplicates for dumbbells)
    all_loop_atoms = primary_loops.copy()
    for db_atom in dumbbell_atoms:
        all_loop_atoms.append(db_atom)
    all_loop_atoms.sort()
    
    return all_loop_atoms, primary_loops, dumbbell_atoms