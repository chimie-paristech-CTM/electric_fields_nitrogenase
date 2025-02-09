from Bio.PDB import *
import numpy as np
import re
import os
import subprocess

import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def extract_fragment_around_active_site(pdb_file, enzyme_name, reactive_residue_name, cofactor_index, 
                                        radius=8, displacement=True, xyz_file="fragment.xyz"):
    """
    Extracts a fragment of the enzyme structure around the active site and saves the atomic coordinates to an XYZ file.

    This function identifies a region around a specified reactive residue and cofactor in the enzyme structure 
    within a given radius, and extracts the atomic coordinates for atoms within this region. The resulting fragment 
    is then saved in an XYZ format file.

    Args:
        pdb_file (str): Path to the input PDB file containing the enzyme structure.
        enzyme_name (str): Name of the enzyme, used for creating the directory and output file.
        reactive_residue_name (str): The name of the reactive residue (e.g., "CYS", "HIS") in the active site.
        cofactor_index (int): The index of the cofactor involved in the active site chemistry.
        radius (float, optional): The radius around the active site from which atoms will be extracted (default is 8 Å).
        displacement (bool, optional): Whether to displace the coordinates when saving them (default is True).
        xyz_file (str, optional): The name of the output XYZ file containing the atomic coordinates of the extracted fragment (default is "fragment.xyz").

    Returns:
        tuple: A tuple containing:
            - nearby_atom_list (list): A list of the atomic identifiers for atoms in the extracted fragment.
            - other_metal_atom_coord (list): A list of coordinates for metal atoms (if any) in the vicinity of the active site.
    """
    output_file = os.path.join(enzyme_name, xyz_file)
    fragment_atoms, other_metal_atom_coord = extract_fragment(pdb_file, cofactor_index, reactive_residue_name, radius=radius)
    nearby_atom_list = obtain_and_save_coord(fragment_atoms, output_file, displacement)

    return nearby_atom_list, other_metal_atom_coord


def get_coord_list(pdb_cofactor, path, reactive_residue_name, cofactor_index):
    """    
    Retrieves coordinates of Fe atoms from a PDB cofactor file for a specific reactive residue.

    Parameters:
    pdb_cofactor (pandas.DataFrame): DataFrame containing PDB data of the cofactor atoms.
    path (str): Path to directory where output files will be saved.
    reactive_residue_name (str): Name of the reactive residue for which Fe atom coordinates are retrieved.
    cofactor_index (int): The index of the cofactor to select.

    Returns:
    list: List of lists containing coordinates (x, y, z) of Fe atoms associated with the specified residue.

    Workflow:
    1. Extracts the unique chain ID from `pdb_cofactor`.
    2. Filters `pdb_cofactor` to retrieve Fe atoms belonging to the specified chain and residue.
    3. Saves the filtered Fe atom positions to 'Fe_positions.csv' in the specified `path`.
    4. Constructs and returns a list of lists containing coordinates of Fe atoms.
    """
    chain_id = pdb_cofactor['chain_id'].unique()[cofactor_index] # you can adjust this
    coord_list = []
    df_fe = pdb_cofactor[pdb_cofactor['Atom'] == 'Fe']
    df_fe = df_fe[df_fe['chain_id'] == f'{chain_id}']
    df_fe = df_fe[df_fe['AA'] == f'{reactive_residue_name}']
    df_fe.to_csv(os.path.join(path, 'Fe_positions.csv'))
    for _, row in df_fe.iterrows():
        coord = list([row['x'], row['y'], row['z']])
        coord_list.append(coord)

    return coord_list


def extract_fragment(pdb_file, cofactor_index, reactive_residue_name, radius=5):
    """
    Extracts a fragment of the enzyme structure around the active site and saves the atomic coordinates to an XYZ file.

    This function identifies a region around a specified reactive residue and cofactor in the enzyme structure 
    within a given radius, and extracts the atomic coordinates for atoms within this region. The resulting fragment 
    is then saved in an XYZ format file.

    Args:
        pdb_file (str): Path to the input PDB file containing the enzyme structure.
        enzyme_name (str): Name of the enzyme, used for creating the directory and output file.
        reactive_residue_name (str): The name of the reactive residue (e.g., "CYS", "HIS") in the active site.
        cofactor_index (int): The index of the cofactor involved in the active site chemistry.
        radius (float, optional): The radius around the active site from which atoms will be extracted (default is 8 Å).
        displacement (bool, optional): Whether to displace the coordinates when saving them (default is True).
        xyz_file (str, optional): The name of the output XYZ file containing the atomic coordinates of the extracted fragment (default is "fragment.xyz").

    Returns:
        tuple: A tuple containing:
            - nearby_atom_list (list): A list of the atomic identifiers for atoms in the extracted fragment.
            - other_metal_atom_coord (list): A list of coordinates for metal atoms (if any) in the vicinity of the active site.
    """
    fragment_atoms = []
    other_metal_atom_coord = None

    # get array of cofactor atoms
    cofactor_atoms = extract_cofactor_atoms(pdb_file, cofactor_index, reactive_residue_name.strip())
    cofactor_coord_array = np.array([atom.coord for atom in cofactor_atoms])

    for atom in cofactor_atoms:
        if atom.element == 'MO' or atom.element == 'V':
            other_metal_atom_coord = np.array([atom.coord]) # store the coordinates of non-iron metal atoms

    all_atoms = extract_all_atoms(pdb_file)

    # now iterate through all atoms and get the ones within the selected radius of the cofactor atoms
    for atom in all_atoms:
        squared_distances = np.sum((cofactor_coord_array - atom.coord) ** 2, axis=1)
        if np.any(np.sqrt(squared_distances) <= radius):
            fragment_atoms.append(atom)

    return fragment_atoms, other_metal_atom_coord


def obtain_and_save_coord(fragment_atoms, output_file, displacement=True):
    """    
    Obtains coordinates of atoms in a fragment and saves them to a specified output file.

    Parameters:
    fragment_atoms (list): A list of Atom objects representing atoms in the fragment.
    output_file (str): The path to the output file where coordinates will be saved.

    Returns:
    list: A list of numpy arrays representing coordinates of atoms in the fragment.
    """
    coord_list = []
    disp = fragment_atoms[0].coord
    with open(output_file, 'w') as f:
        f.write(f"{len(fragment_atoms)}\n\n")
        for atom in fragment_atoms:
            coord_list.append(np.array(atom.coord))
            if displacement:
                f.write(f"{atom.element} {atom.coord[0] - disp[0]} {atom.coord[1] - disp[1]} {atom.coord[2] - disp[2]}\n")
            else:
                f.write(f"{atom.element} {atom.coord[0]} {atom.coord[1]} {atom.coord[2]}\n")

    return coord_list


class Atom:
    def __init__(self, element, coord):
        self.element = element
        self.coord = coord


def extract_cofactor_atoms(pdb_file, cofactor_index, reactive_residue_name):
    """    
    Extracts atoms of a specific cofactor residue from a PDB file.

    Parameters:
    pdb_file (str): The path to the PDB file containing the enzyme structure.
    cofactor_index (int): The index of the specific cofactor residue to extract.
    reactive_residue_name (str): The name of the reactive residue to identify the cofactor.

    Returns:
    list: A list of Atom objects representing atoms of the specified cofactor residue.
    """
    parser = PDBParser()
    structure = parser.get_structure('enzyme', pdb_file)

    i = 0
    cofactor_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if reactive_residue_name in residue.id[0]:
                    if i != cofactor_index:
                        i += 1
                        break
                    else:
                        for atom in residue.get_atoms():
                            cofactor_atoms.append(Atom(atom.element, atom.coord))
                        i += 1

    return cofactor_atoms


def extract_all_atoms(pdb_file):
    """
    Extracts all atoms from a PDB file and returns them as a list of Atom objects.

    Parameters:
    pdb_file (str): The path to the PDB file to read atoms from.

    Returns:
    list: A list of Atom objects, each representing an atom extracted from the PDB file.
    """
    all_atoms = []

    with open(pdb_file, 'r') as file:
        for line in file:
            try:
                element = line.split()[-1]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except:
                print(line)
                continue
            coord = np.array([x, y, z])
            all_atoms.append(Atom(element, coord))

    return all_atoms
