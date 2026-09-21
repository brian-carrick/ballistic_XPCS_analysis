"""
Bond Exchange Module for Reversible Bond Formation
Implements tau-Gillespie algorithm for stochastic chemical kinetics
Integrated with Gusev network generation format

Three independent reaction channels:
1. Intermolecular associations (between different atoms) - OPTIMIZED WITH CELL LIST
2. Intramolecular associations (loops on same atom)
3. Dissociations (bond breaking)
"""

import numpy as np
from numpy import linalg as LA
import random
import math


class BondExchange:
    """
    Class to manage reversible bond exchange using tau-Gillespie algorithm
    Maintains compatibility with Gusev network generation and Optimizer class
    Uses cell list for O(n_atoms) scaling of intermolecular associations
    
    Properly handles self-bonds (loops) that are already encoded in the bonds array
    """
    
    def __init__(self, optimizer, bonds_array, N, b, k_A, k_D, func, dt, Lx, Ly, Lz):
        # NO loop_atoms parameter - extract from bonds_array instead
        
        self.optimizer = optimizer
        self.N = N
        self.b = b
        self.k_A = k_A
        self.k_D = k_D
        self.func = func
        self.dt = dt
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        
        self.n_links = len(optimizer.atoms)
        self.bond_type = 1
        self.all_bonds = np.copy(bonds_array) if len(bonds_array) > 0 else np.zeros((0, 4), dtype=int)
        
        Nb_sq = N * b**2
        self.prob_const = (3.0 / (2.0 * np.pi * Nb_sq))**1.5
        
        self.conn = self._build_connection_matrix()
        
        # ===== EXTRACT loop_atoms FROM self-bonds in bonds_array =====
        self.loop_atoms_raw = []
        for bond in self.all_bonds:
            link_1 = int(bond[2]) - 1
            link_2 = int(bond[3]) - 1
            if link_1 == link_2:  # Self-bond
                self.loop_atoms_raw.append(link_1)
        
        self._analyze_loops()
        
        # Tracking statistics
        self.n_associations_inter = 0
        self.n_associations_loop = 0
        self.n_dissociations = 0
        self.alpha_inter_history = []
        self.alpha_loop_history = []
        self.alpha_dissoc_history = []
        self.bond_counter = len(self.all_bonds) + 1
        self.broken_bonds = []
        
        print(f'\nBondExchange initialized:')
        print(f'  - {self.n_links} atoms (junctions)')
        print(f'  - {len(self.all_bonds)} total bonds')
        print(f'  - {len(self.loop_atoms_raw)} loop chains (from self-bonds)')
        print(f'  - {self.n_loops} primary loops')
        print(f'  - {self.n_dumbbells} dumbbells')
        print(f'  - Total chains: {len(self.all_bonds)}')
        print()

    # def _build_connection_matrix(self):
    #     """
    #     Build connection matrix from ALL bonds (both inter-atom and self-bonds)
    #     """
    #     conn = np.full((self.n_links, self.func), -1, dtype=int)
        
    #     # Track which sticker slots have been used
    #     used_slots = {}
    #     for i in range(self.n_links):
    #         used_slots[i] = 0
        
    #     # Add all bonds
    #     for bond in self.all_bonds:
    #         link_1 = int(bond[2]) - 1  # Convert to 0-indexed
    #         link_2 = int(bond[3]) - 1  # Convert to 0-indexed
            
    #         if link_1 < 0 or link_2 < 0 or link_1 >= self.n_links or link_2 >= self.n_links:
    #             print(f'WARNING: Invalid bond link_1={link_1}, link_2={link_2}')
    #             continue
            
    #         # Check if we have available slots
    #         if used_slots[link_1] >= self.func:
    #             print(f'WARNING: Atom {link_1} already has {self.func} bonds!')
    #             continue
            
    #         if used_slots[link_2] >= self.func:
    #             print(f'WARNING: Atom {link_2} already has {self.func} bonds!')
    #             continue
            
    #         # Add bond to first available slot
    #         conn[link_1, used_slots[link_1]] = link_2
    #         conn[link_2, used_slots[link_2]] = link_1
            
    #         used_slots[link_1] += 1
    #         used_slots[link_2] += 1
        
    #     return conn
    
    def _build_connection_matrix(self):
        """
        Build connection matrix from ALL bonds (both inter-atom and self-bonds)
        """
        conn = np.full((self.n_links, self.func), -1, dtype=int)
        
        # Track which sticker slots have been used
        used_slots = {}
        for i in range(self.n_links):
            used_slots[i] = 0
        
        # Add all bonds
        for bond in self.all_bonds:
            link_1 = int(bond[2]) - 1  # Convert to 0-indexed
            link_2 = int(bond[3]) - 1  # Convert to 0-indexed
            
            if link_1 < 0 or link_2 < 0 or link_1 >= self.n_links or link_2 >= self.n_links:
                print(f'WARNING: Invalid bond link_1={link_1}, link_2={link_2}')
                continue
            
            # HANDLE SELF-LOOPS SPECIALLY
            if link_1 == link_2:
                # Self-loop uses 2 slots on the same atom
                if used_slots[link_1] + 1 >= self.func:
                    print(f'WARNING: Atom {link_1} cannot fit self-loop ({used_slots[link_1] + 1} slots needed, func={self.func})')
                    continue
                
                conn[link_1, used_slots[link_1]] = link_1
                conn[link_1, used_slots[link_1] + 1] = link_1
                used_slots[link_1] += 2
            
            else:
                # Regular inter-atom bond
                if used_slots[link_1] >= self.func:
                    print(f'WARNING: Atom {link_1} already has {self.func} bonds!')
                    continue
                
                if used_slots[link_2] >= self.func:
                    print(f'WARNING: Atom {link_2} already has {self.func} bonds!')
                    continue
                
                # Add bond to first available slot
                conn[link_1, used_slots[link_1]] = link_2
                conn[link_2, used_slots[link_2]] = link_1
                
                used_slots[link_1] += 1
                used_slots[link_2] += 1
        
        return conn

    
    def _find_bonded_pairs_from_conn(self):
        """
        Find all bonded pairs from connection matrix
        Correctly handles self-loops (counts each once)
        """
        bonded_pairs = []
        visited_self_loops = set()  # Track which atoms have already been added as loops
        
        for i in range(self.n_links):
            for bond_idx in range(self.func):
                j = self.conn[i, bond_idx]
                
                if j == -1 or j >= self.n_links:
                    continue
                
                # Handle self-loops (only count once per atom)
                if i == j:
                    if i not in visited_self_loops:
                        bonded_pairs.append({'i': i, 'j': j, 'idx_i': bond_idx})
                        visited_self_loops.add(i)
                else:
                    # For inter-atom bonds, avoid counting same bond twice
                    bond_key = tuple(sorted([i, j, bond_idx]))
                    # Check if this exact bond hasn't been added yet
                    bond_exists = False
                    for existing_pair in bonded_pairs:
                        if existing_pair['i'] == i and existing_pair['j'] == j:
                            bond_exists = True
                            break
                        if existing_pair['i'] == j and existing_pair['j'] == i:
                            bond_exists = True
                            break
                    
                    if not bond_exists and i < j:  # Only add once, when i < j
                        bonded_pairs.append({'i': i, 'j': j, 'idx_i': bond_idx})
        
        return bonded_pairs

    def _extract_loops_from_conn(self):
        """
        Extract loop chains from self-connections in conn matrix
        An atom with 2 self-connections (2 sticker slots to itself) = 1 loop chain
        An atom with 4 self-connections (dumbbell) = 2 loop chains
        """
        loop_atoms = []
        
        for i in range(self.n_links):
            # Count how many times this atom is connected to itself
            self_connections = np.sum(self.conn[i, :] == i)
            
            # Each pair of self-connections = 1 loop chain
            # (1 loop uses 2 sticker slots)
            n_loops_at_atom = self_connections // 2
            
            for _ in range(n_loops_at_atom):
                loop_atoms.append(i)
        
        return loop_atoms
    
    def _analyze_loops(self):
        """
        Analyze loop atom statistics
        loop_atoms_raw contains duplicates for atoms with multiple loops
        """
        # Count occurrences of each atom in loop_atoms list
        loop_atom_counts = {}
        for atom in self.loop_atoms_raw:
            atom_idx = int(atom)
            loop_atom_counts[atom_idx] = loop_atom_counts.get(atom_idx, 0) + 1
        
        # Separate into single loops and dumbbells
        self.loop_atoms_unique = []
        self.n_dumbbells = 0
        self.n_loops = 0
        
        for atom_idx, count in loop_atom_counts.items():
            self.loop_atoms_unique.append(atom_idx)
            if count == 1:
                self.n_loops += 1
            elif count == 2:
                self.n_dumbbells += 1
        
        self.total_loop_chains = self.n_loops + 2 * self.n_dumbbells
        
    
    def _sync_chains_and_bonds(self):
        """
        Update bonds array from connection matrix
        Returns both inter-atom bonds AND self-bonds
        """
        
        bonded_pairs = self._find_bonded_pairs_from_conn()
        
        bonds = []
        visited_bonds = set()
        bond_index = 1
        
        for pair in bonded_pairs:
            atom_i = pair['i']
            atom_j = pair['j']
            
            # Avoid counting same bond twice
            bond_key = tuple(sorted([atom_i, atom_j]))
            if bond_key in visited_bonds:
                continue
            visited_bonds.add(bond_key)
            
            # Add ALL bonds (both inter-atom and self-bonds)
            bond = np.array([self.bond_type, bond_index, atom_i + 1, atom_j + 1], dtype=int)
            bonds.append(bond)
            bond_index += 1
        
        bonds = np.array(bonds, dtype=int) if len(bonds) > 0 else np.zeros((0, 4), dtype=int)
        self.bond_counter = bond_index
        
        # Extract loop atoms from bonds
        new_loop_atoms = []
        for bond in bonds:
            link_1 = int(bond[2]) - 1
            link_2 = int(bond[3]) - 1
            if link_1 == link_2:
                new_loop_atoms.append(link_1)
        
        # ===== VALIDATION =====
        max_bonds = (self.n_links * self.func) // 2
        if len(bonds) > max_bonds:
            print(f'\n⚠ ERROR in _sync_chains_and_bonds:')
            print(f'Total bonds ({len(bonds)}) exceeds maximum ({max_bonds})')
            print(f'Inter-atom bonds: {len(bonds) - len(new_loop_atoms)}')
            print(f'Self-bonds/loops: {len(new_loop_atoms)}')
            
            # Count sticker usage
            f_i, F_total = self.calculate_free_stickers()
            bonded_stickers = self.n_links * self.func - F_total
            print(f'Sticker usage: {bonded_stickers} / {self.n_links * self.func}')
            
            raise ValueError(f'Bond limit exceeded during sync: {len(bonds)} > {max_bonds}')
        
        return bonds, new_loop_atoms
    

    
    
    def calculate_free_stickers(self):
        """
        Calculate number of free stickers on each atom
        
        Returns:
        --------
        f_i : ndarray
            Number of free stickers on each atom
        F_total : int
            Total number of free stickers in system
        """
        f_i = np.sum(self.conn == -1, axis=1)
        F_total = np.sum(f_i)
        return f_i, F_total
    
    
    def _periodic_distance(self, r_vec):
        """
        Apply periodic boundary conditions to distance vector
        """
        r_vec = np.copy(r_vec)
        r_vec[0] -= int(round(r_vec[0] / self.Lx)) * self.Lx
        r_vec[1] -= int(round(r_vec[1] / self.Ly)) * self.Ly
        r_vec[2] -= int(round(r_vec[2] / self.Lz)) * self.Lz
        return r_vec
    
    
    def _get_bond_distance(self, atom_i, atom_j):
        """
        Calculate distance between two atoms
        """
        delta_r = self.optimizer.atoms[atom_j, :3] - self.optimizer.atoms[atom_i, :3]
        delta_r = self._periodic_distance(delta_r)
        return LA.norm(delta_r)
    
    
    def _build_cell_list(self, cutoff_distance):
        """
        Build a cell list for O(n) neighbor finding
        """
        cell_size = cutoff_distance + 0.1
        n_cells_x = max(1, int(self.Lx / cell_size))
        n_cells_y = max(1, int(self.Ly / cell_size))
        n_cells_z = max(1, int(self.Lz / cell_size))
        
        cells = {}
        xlo = self.optimizer.xlo
        ylo = self.optimizer.ylo
        zlo = self.optimizer.zlo
        
        for i in range(self.n_links):
            x_shifted = self.optimizer.atoms[i, 0] - xlo
            y_shifted = self.optimizer.atoms[i, 1] - ylo
            z_shifted = self.optimizer.atoms[i, 2] - zlo
            
            ix = int(x_shifted / cell_size) % n_cells_x
            iy = int(y_shifted / cell_size) % n_cells_y
            iz = int(z_shifted / cell_size) % n_cells_z
            
            cell_key = (ix, iy, iz)
            
            if cell_key not in cells:
                cells[cell_key] = []
            cells[cell_key].append(i)
        
        return cells, cell_size, n_cells_x, n_cells_y, n_cells_z
    
    def _build_connection_matrix_from_bonds(self, bonds_array):
        """
        Rebuild connection matrix from bonds array
        Used after syncing to ensure consistency
        """
        conn = np.full((self.n_links, self.func), -1, dtype=int)
        
        # Add all bonds from the bonds array
        for bond in bonds_array:
            link_1 = int(bond[2]) - 1  # Convert to 0-indexed
            link_2 = int(bond[3]) - 1  # Convert to 0-indexed
            
            if link_1 < 0 or link_2 < 0 or link_1 >= self.n_links or link_2 >= self.n_links:
                continue
            
            # Find empty slots on both atoms
            empty_1 = np.where(conn[link_1, :] == -1)[0]
            empty_2 = np.where(conn[link_2, :] == -1)[0]
            
            if len(empty_1) > 0 and len(empty_2) > 0:
                conn[link_1, empty_1[0]] = link_2
                conn[link_2, empty_2[0]] = link_1
        
        return conn

    def tau_gillespie_step(self, step_num=0, verbose=True):
        """
        Perform one tau-leap time step of the tau-Gillespie algorithm
        with sequential draw-and-attempt for all reaction types.
        Bonds are converted directly from self.conn to avoid stale data.
        """

        if verbose:
            print(f'\n=== Step {step_num} START ===')
            print(f'Bonds before step: {len(self.optimizer.bonds)}')
        
        f_i, F_total = self.calculate_free_stickers()
        
        if verbose:
            print(f"Step {step_num} START: n_bonds={len(self.optimizer.bonds)}, F_total={F_total}, n_free_atoms={np.sum(f_i > 0)}")
        
        # ==================== PHASE 1: CALCULATE INTERMOLECULAR PROPENSITY & DRAW ====================
        
        alpha_inter = 0.0
        inter_pairs = []
        
        # cutoff_distance = 2.0 * np.sqrt(self.N * self.b**2)
        cutoff_distance = self.N * self.b
        cells, cell_size, n_cells_x, n_cells_y, n_cells_z = self._build_cell_list(cutoff_distance)
        
        Nb_sq = self.N * self.b**2
        
        for cell_key, atom_list in cells.items():
            ix, iy, iz = cell_key
            
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    for dz in [-1, 0, 1]:
                        neighbor_key = ((ix+dx) % n_cells_x, (iy+dy) % n_cells_y, (iz+dz) % n_cells_z)
                        
                        if neighbor_key not in cells:
                            continue
                        
                        neighbor_list = cells[neighbor_key]
                        
                        for i in atom_list:
                            if f_i[i] <= 0:
                                continue
                            
                            for j in neighbor_list:
                                if i >= j or f_i[j] <= 0:
                                    continue
                                
                                delta_r = self.optimizer.atoms[j, :3] - self.optimizer.atoms[i, :3]
                                delta_r = self._periodic_distance(delta_r)
                                
                                r_sq = np.sum(delta_r**2)
                                P_rij = self.prob_const * np.exp(-1.5 * r_sq / Nb_sq)
                                contrib = self.k_A * f_i[i] * f_i[j] * P_rij
                                alpha_inter += contrib
                                
                                inter_pairs.append({'i': i, 'j': j, 'P_rij': P_rij})
        
        # ==================== PHASE 2: CALCULATE INTRAMOLECULAR PROPENSITY & DRAW ====================
        
        alpha_loop = 0.0
        for i in range(self.n_links):
            if f_i[i] >= 2:
                alpha_loop += self.k_A * f_i[i] * (f_i[i] - 1) * self.prob_const
        
        # ==================== PHASE 3: CALCULATE DISSOCIATION PROPENSITY & DRAW ====================
        
        n_bonded_sites = self.n_links * self.func - F_total
        alpha_dissoc = self.k_D * n_bonded_sites / 2.0
        
        bonded_pairs_all = self._find_bonded_pairs_from_conn()
        
        # ==================== PHASE 4: TOTAL PROPENSITY & DRAW REACTION NUMBERS ====================
        
        alpha_total = alpha_inter + alpha_loop + alpha_dissoc
        
        self.alpha_inter_history.append(alpha_inter)
        self.alpha_loop_history.append(alpha_loop)
        self.alpha_dissoc_history.append(alpha_dissoc)
        
        if verbose:
            print(f"  Propensities: inter={alpha_inter:.6e}, loop={alpha_loop:.6e}, dissoc={alpha_dissoc:.6e}, total={alpha_total:.6e}")
        
        # Draw fixed numbers upfront
        n_inter_target = 0
        n_loop_target = 0
        n_dissoc_target = 0
        
        if alpha_total > 0:
            n_inter_target = np.random.poisson(alpha_inter * self.dt)
            n_loop_target = np.random.poisson(alpha_loop * self.dt)
            n_dissoc_target = np.random.poisson(alpha_dissoc * self.dt)
            
            if verbose:
                print(f"  Target reactions: inter={n_inter_target}, loop={n_loop_target}, dissoc={n_dissoc_target}")
        
        # ==================== PHASE 5: SEQUENTIAL INTERMOLECULAR ASSOCIATIONS ====================
        
        n_inter_executed = 0
        bonds_formed_inter = []
        candidate_inter_pairs = list(inter_pairs)
        
        while n_inter_executed < n_inter_target and len(candidate_inter_pairs) > 0:
            # RECALCULATE STATE
            f_i, F_total = self.calculate_free_stickers()
            
            # Build propensities from CURRENT f_i
            propensities = []
            valid_indices = []
            
            for idx, pair in enumerate(candidate_inter_pairs):
                i, j = pair['i'], pair['j']
                if f_i[i] > 0 and f_i[j] > 0:
                    prop = self.k_A * f_i[i] * f_i[j] * pair['P_rij']
                    propensities.append(prop)
                    valid_indices.append(idx)
            
            if not valid_indices:
                break
            
            propensities = np.array(propensities)
            propensities = propensities / np.sum(propensities)
            
            # Draw one pair
            choice = np.random.choice(len(valid_indices), p=propensities)
            pair_idx = valid_indices[choice]
            pair = candidate_inter_pairs[pair_idx]
            
            # ATTEMPT BOND
            success = self._associate_atoms(pair['i'], pair['j'])
            if success:
                n_inter_executed += 1
                self.n_associations_inter += 1
                dist = self._get_bond_distance(pair['i'], pair['j'])
                bonds_formed_inter.append({'atom_i': pair['i'] + 1, 'atom_j': pair['j'] + 1, 'distance': dist})
                if verbose and n_inter_executed <= 5:
                    print(f"    Inter: Atom {pair['i']+1:3d} -- Atom {pair['j']+1:3d} bonded (dist={dist:.4f})")
            else:
                # Remove this pair if it failed
                candidate_inter_pairs.pop(pair_idx)
        
        if verbose and n_inter_executed > 5:
            print(f"    ... and {n_inter_executed - 5} more intermolecular bonds formed")
        if verbose:
            print(f"  Executed intermolecular: {n_inter_executed} / {n_inter_target}")
        
        # ==================== PHASE 6: SEQUENTIAL INTRAMOLECULAR ASSOCIATIONS ====================
        
        n_loop_executed = 0
        bonds_formed_loop = []
        
        while n_loop_executed < n_loop_target:
            # RECALCULATE STATE
            f_i, F_total = self.calculate_free_stickers()
            
            # Find valid atoms (f_i >= 2)
            valid_loop_atoms = [i for i in range(self.n_links) if f_i[i] >= 2]
            
            if not valid_loop_atoms:
                break
            
            # Draw one atom
            atom_idx = np.random.choice(valid_loop_atoms)
            
            # ATTEMPT LOOP
            success = self._associate_atoms(atom_idx, atom_idx)
            if success:
                n_loop_executed += 1
                self.n_associations_loop += 1
                bonds_formed_loop.append({'atom': atom_idx + 1, 'f_i_after': f_i[atom_idx] - 2})
                if verbose and n_loop_executed <= 5:
                    print(f"    Loop: Atom {atom_idx+1:3d} self-loop formed")
            # If fails, just continue trying other atoms (don't remove)
        
        if verbose and n_loop_executed > 5:
            print(f"    ... and {n_loop_executed - 5} more intramolecular loops formed")
        if verbose:
            print(f"  Executed intramolecular: {n_loop_executed} / {n_loop_target}")
        
        # ==================== PHASE 7: SEQUENTIAL DISSOCIATIONS ====================
        
        n_dissoc_executed = 0
        bonds_broken_this_step = []
        
        while n_dissoc_executed < n_dissoc_target:
            # RECALCULATE STATE (get fresh bonded pairs)
            f_i, F_total = self.calculate_free_stickers()
            bonded_pairs_current = self._find_bonded_pairs_from_conn()
            
            if not bonded_pairs_current:
                break
            
            # Draw one pair uniformly
            pair_idx = np.random.randint(0, len(bonded_pairs_current))
            pair = bonded_pairs_current[pair_idx]
            
            atom_i = pair['i']
            atom_j = pair['j']
            idx_i = pair['idx_i']
            
            # VERIFY bond still exists (in case it was already broken)
            if self.conn[atom_i, idx_i] != atom_j:
                continue
            
            dist_before = self._get_bond_distance(atom_i, atom_j)
            
            # BREAK BOND
            bond_broken = False
            for idx_j in range(self.func):
                if self.conn[atom_j, idx_j] == atom_i:
                    self.conn[atom_i, idx_i] = -1
                    self.conn[atom_j, idx_j] = -1
                    n_dissoc_executed += 1
                    self.n_dissociations += 1
                    bond_broken = True
                    
                    bonds_broken_this_step.append({'atom_i': atom_i + 1, 'atom_j': atom_j + 1, 'distance': dist_before})
                    self.broken_bonds.append([step_num, atom_i + 1, atom_j + 1, dist_before])
                    
                    if verbose and n_dissoc_executed <= 5:
                        print(f"    Dissoc: Atom {atom_i+1:3d} -- Atom {atom_j+1:3d} bond broken (dist={dist_before:.4f})")
                    break
        
        if verbose and n_dissoc_executed > 5:
            print(f"    ... and {n_dissoc_executed - 5} more bonds dissociated")
        if verbose:
            print(f"  Executed dissociations: {n_dissoc_executed} / {n_dissoc_target}")

            # ==================== PHASE 7.5: DEBUG CONNECTIVITY MATRIX ====================

        if verbose:
            print("\n--- CONNECTIVITY MATRIX DIAGNOSTICS ---")
            
            # Check for invalid entries
            for i in range(self.n_links):
                bonded_count = np.sum(self.conn[i, :] >= 0)
                
                if bonded_count > self.func:
                    print(f"ERROR: Atom {i} has {bonded_count} bonds (exceeds func={self.func})")
                    print(f"  conn[{i}, :] = {self.conn[i, :]}")
                
                # Check for duplicates within a row
                bonded_atoms = self.conn[i, self.conn[i, :] >= 0]
                if len(bonded_atoms) > 0:
                    unique, counts = np.unique(bonded_atoms, return_counts=True)
                    for atom, count in zip(unique, counts):
                        if count > 1:
                            print(f"WARNING: Atom {i} bonded {count} times to atom {atom}")
                            print(f"  conn[{i}, :] = {self.conn[i, :]}")
            
            # Check for asymmetric bonds
            for i in range(self.n_links):
                for idx in range(self.func):
                    j = self.conn[i, idx]
                    if j >= 0 and j != i:  # Inter-atom bond
                        found_reverse = False
                        for idx_j in range(self.func):
                            if self.conn[j, idx_j] == i:
                                found_reverse = True
                                break
                        if not found_reverse:
                            print(f"ERROR: Asymmetric bond! conn[{i}, {idx}] = {j}, but {j} doesn't bond back to {i}")
            
            print("--- END DIAGNOSTICS ---\n")

        
        # ==================== PHASE 8: CONVERT FINAL CONN TO BONDS ====================

        if verbose:
            print(f'  Before conversion: {len(self.optimizer.bonds)} bonds')

        # Build bonds directly from the final self.conn state
        new_bonds = []
        bond_index = 1

        # All bonds: inter-atom (i < j) + self-loops (i == i)
        for i in range(self.n_links):
            for idx in range(self.func):
                j = self.conn[i, idx]
                if j >= 0 and i < j:  # Only count when i < j to avoid duplicates from being listed twice
                    bond = np.array([self.bond_type, bond_index, i + 1, j + 1], dtype=int)
                    new_bonds.append(bond)
                    bond_index += 1

        # Self-loops: count stickers bonded to self, divide by 2
        for i in range(self.n_links):
            self_loop_stickers = np.sum(self.conn[i, :] == i)
            if self_loop_stickers > 0:
                n_self_loops = self_loop_stickers // 2
                for _ in range(n_self_loops):
                    bond = np.array([self.bond_type, bond_index, i + 1, i + 1], dtype=int)
                    new_bonds.append(bond)
                    bond_index += 1

        # Convert to numpy array
        if len(new_bonds) > 0:
            self.optimizer.bonds = np.array(new_bonds, dtype=int)
        else:
            self.optimizer.bonds = np.zeros((0, 4), dtype=int)

        self.bond_counter = bond_index

        if verbose:
            print(f'  After conversion: {len(self.optimizer.bonds)} bonds')

        # ==================== PHASE 9: IDENTIFY LOOP ATOMS & VALIDATE VALENCY ====================

        loop_atoms_raw = []
        for i in range(self.n_links):
            # Check if atom has a self-loop
            has_self_loop = np.any(self.conn[i, :] == i)
            if has_self_loop:
                loop_atoms_raw.append(i)

        self.loop_atoms_raw = loop_atoms_raw

        if verbose:
            print(f'  Identified {len(loop_atoms_raw)} atoms with self-loops')

        # VALENCY CHECK: Only warn if an atom exceeds its functionality
        for i in range(self.n_links):
            bonded_count = np.sum(self.conn[i, :] >= 0)
            if bonded_count > self.func:
                print(f"ERROR: Atom {i+1} has {bonded_count} total bonds (exceeds func={self.func})")
                print(f"  conn[{i}, :] = {self.conn[i, :]}")


        # ==================== PHASE 10: FINAL STATE ====================
        
        f_i_final, F_total_final = self.calculate_free_stickers()
        
        if verbose:
            print(f'=== Step {step_num} END ===')
            print(f'Bonds after step: {len(self.optimizer.bonds)}')
            print(f"Step {step_num} END: n_bonds={len(self.optimizer.bonds)}, F_total={F_total_final}, n_loops={len(self.loop_atoms_raw)}")
            print()
        
        # ==================== DEBUG: VALIDATE FINAL CONN ====================
        
        for i in range(self.n_links):
            bonded_count = np.sum(self.conn[i, :] >= 0)
            if bonded_count > self.func:
                print(f"ERROR: Atom {i} has {bonded_count} bonds (exceeds func={self.func})")
                print(f"  conn row: {self.conn[i, :]}")
        
        # ==================== RETURN STATISTICS ====================
        
        step_data = {
            'n_inter': n_inter_executed,
            'n_loop': n_loop_executed,
            'n_dissoc': n_dissoc_executed,
            'n_inter_target': n_inter_target,
            'n_loop_target': n_loop_target,
            'n_dissoc_target': n_dissoc_target,
            'alpha_inter': alpha_inter,
            'alpha_loop': alpha_loop,
            'alpha_dissoc': alpha_dissoc,
            'alpha_total': alpha_total,
            'F_total': F_total_final,
            'n_bonds': len(self.optimizer.bonds),
            'n_bonds_bonded': self.n_links * self.func - F_total_final,
            'n_loops': self.n_loops,
            'n_dumbbells': self.n_dumbbells,
            'total_loop_chains': len(self.loop_atoms_raw),
            'bonds_broken': bonds_broken_this_step,
            'bonds_formed_inter': bonds_formed_inter,
            'bonds_formed_loop': bonds_formed_loop
        }
        
        return step_data
    
    
    def _associate_atoms(self, atom_i, atom_j):
        """
        Form a bond between atom i and atom j (or self-loop if i==j)
        """
        free_i = np.where(self.conn[atom_i, :] == -1)[0]
        free_j = np.where(self.conn[atom_j, :] == -1)[0]
        
        if len(free_i) == 0 or len(free_j) == 0:
            return False
        
        if atom_i == atom_j and len(free_i) < 2:
            return False
        
        idx_i = free_i[0]
        idx_j = free_j[0]
        
        self.conn[atom_i, idx_i] = atom_j
        self.conn[atom_j, idx_j] = atom_i
        
        return True
    
    
    def update_box_dimensions(self, Lx, Ly, Lz):
        """
        Update box dimensions if system is rescaled
        """
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
    
    
    def get_statistics(self):
        """
        Get cumulative statistics from all steps
        """
        stats = {
            'total_associations_inter': self.n_associations_inter,
            'total_associations_loop': self.n_associations_loop,
            'total_dissociations': self.n_dissociations,
            'total_bonds_formed': self.n_associations_inter + self.n_associations_loop,
            'net_bonds_change': self.n_associations_inter + self.n_associations_loop - self.n_dissociations,
            'broken_bonds': self.broken_bonds
        }
        return stats
