from Bio.PDB import *
import numpy as np
import re
import os
import subprocess

import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def extract_fragment_around_active_site(enzyme_name, reactive_residue_name, cofactor_index, 
                                        output_xyz_file, radius=8, displacement=True, pdb_file=None):
    """
    Extracts the atomic fragment surrounding the active site of an enzyme within a specified radius.

    This function identifies atoms near the active site of a given enzyme and saves the coordinates
    of the fragment into an output file. It also optionally adjusts atom positions to account for
    structural displacement. Hydrogens are added to ensure a complete molecular structure before
    the extraction process.

    Parameters:
    ----------
    enzyme_name : str
        The name of the enzyme to process (e.g., a PDB file or identifier).
    reactive_residue_name : str
        The name of the reactive residue within the enzyme (e.g., "CYS").
    cofactor_index : int
        The index of the cofactor associated with the active site (typically in the enzyme structure file).
    output_xyz_file : str
        Path to the file where the extracted fragment will be saved in XYZ format.
    radius : float, optional
        The radius (in angstroms) around the active site within which atoms are included. Default is 8 Å.
    displacement : bool, optional
        If True, accounts for atomic displacement during fragment extraction. Default is True.

    Returns:
    -------
    nearby_atom_list : list
        A list of atoms within the specified radius of the active site.
    other_metal_atom_coord : list
        Coordinates of any additional metal atoms identified near the active site.
 
    """
    clean_up_and_add_hydrogens(enzyme_name, reactive_residue_name)
    nearby_atom_list, other_metal_atom_coord = extract_atoms_around_active_site(
        enzyme_name, reactive_residue_name, cofactor_index, output_xyz_file, radius, displacement, pdb_file)

    return nearby_atom_list, other_metal_atom_coord


def clean_up_and_add_hydrogens(enzyme_name, reactive_residue_name):
    """
    Cleans up the mutated PDB file, adds hydrogen atoms, and removes incorrectly added hydrogens 
    from the specified reactive residue.
    
    Parameters:
    enzyme_name (str): The name of the enzyme (also used as the directory name).
    reactive_residue_name (str): The name of the reactive residue to avoid adding hydrogen atoms to.
    """
    # Clean up the mutated pdb file (and add hydrogen atoms this time)
    input_file = f'{enzyme_name}_mutated.pdb'
    output_file = f'{enzyme_name}_mutated_clean.pdb' 

    command = f'obabel {os.path.join(enzyme_name, input_file)} -O {os.path.join(enzyme_name, output_file)} -h'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Print the output and error (if any)
    print(result.stdout)
    print(result.stderr)

    # Remove incorrectly added hydrogens added to the S atoms of the cofactor at the active site
    with open(os.path.join(enzyme_name, output_file), 'r') as f:
        lines = f.readlines()
    with open(os.path.join(enzyme_name, output_file), 'w') as f:
        for line in lines:
            try:
                if line.split()[2] == reactive_residue_name and line.split()[-1] == 'H':
                    continue
                else:
                    f.write(line)
            except:
                f.write(line)

    extract_atoms_pdb_file(enzyme_name, f'{enzyme_name}_mutated')


def extract_atoms_around_active_site(enzyme_name, reactive_residue_name, cofactor_index, 
                                     xyz_file="fragment.xyz", radius=8, displacement=True, pdb_file=None):
    """
    Extracts the atoms around the active site from a PDB file and saves it in an XYZ format.

    Parameters:
    enzyme_name (str): The name of the enzyme, used to construct file paths.
    reactive_residue_name (str): The name of the reactive residue to focus on within the enzyme.
    cofactor_index (int): The index of the cofactor to use in the extraction process.
    xyz_file (str): The name of the XYZ-file to which the fragment will be saved.
    radius (float): The radius to consider for the extraction of the fragment
    displacement (bool): Whether or not to center the fragment coordinates.

    Returns:
    tuple:
        nearby_atom_list (list): A list of atoms that are within the specified radius of the active site.
        other_metal_atom_coord (tuple): Coordinates of other metal atoms in the vicinity of the active site.

    Workflow:
    1. Defines the PDB file path based on the enzyme name.
    2. Sets the output file path for the fragment in XYZ format.
    3. Extracts fragment atoms and coordinates of nearby metal atoms from the PDB file.
    4. Saves the coordinates of the fragment atoms to the output file.
    5. Returns a list of nearby atoms and the coordinates of other metal atoms.
    """
    if pdb_file is not None:
        pass
    else:
        pdb_file = os.path.join(enzyme_name, f'{enzyme_name}_mutated_atoms.pdb')
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


def extract_atoms_pdb_file(path, enzyme_name):
    """
    Extracts ATOM and selected HETATM lines from a cleaned PDB file and writes them to separate output files.
    
    Parameters:
    path (str): The directory path where the input and output files are located.
    enzyme_name (str): The base name of the enzyme file (without extensions).
    
    The function generates two output files:
    - One containing all ATOM records and selected HETATM records (excluding water).
    - Another containing only the ATOM records for potential charge assignment.
    """
    input_file = os.path.join(path, f'{enzyme_name}_clean.pdb')
    output_file = os.path.join(path, f'{enzyme_name}_atoms.pdb')
    output_file_charge = os.path.join(path, f'{enzyme_name}_atoms_charge.pdb') 
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(output_file_charge, "w") as outfile_charge:
        for line in infile:
            if re.match(r'ATOM', line):
                outfile.write(line)
                outfile_charge.write(line)
            elif re.match(r'HETATM', line):
                residue_name = line[17:20].strip()
                if residue_name != 'HOH':
                    outfile.write(line)


def extract_fragment(pdb_file, cofactor_index, reactive_residue_name, radius=5):
    """    
    Extracts a fragment of atoms around a specific cofactor residue in a PDB file.

    Parameters:
    pdb_file (str): The path to the PDB file containing the enzyme structure.
    cofactor_index (int): The index of the specific cofactor residue to extract.
    reactive_residue_name (str): The name of the reactive residue to identify the cofactor.
    radius (float, optional): The radius in Angstroms to define the sphere around the cofactor for fragment extraction.
                              Default is 5 Angstroms.

    Returns:
    tuple: A tuple containing:
        - list: A list of Atom objects representing atoms within the specified radius around the cofactor.
        - numpy.ndarray or None: Coordinates of other metal atoms (like MO or V1) found in the cofactor.

    Workflow:
    1. Extracts atoms of the specific cofactor residue (`cofactor_atoms`) using `extract_cofactor_atoms`.
    2. Constructs a numpy array (`cofactor_coord_array`) containing coordinates of all `cofactor_atoms`.
    3. Searches through all atoms in the PDB file (`all_atoms`) using `extract_all_atoms`.
    4. Computes squared distances between each `all_atoms` and `cofactor_atoms`.
    5. Checks if any distance is within `radius`.
    6. If true, adds the `atom` to `fragment_atoms`.
    7. Stores coordinates of other metal atoms (like MO or V1) found in `cofactor_atoms` in `other_metal_atom_coord`.
    """
    fragment_atoms = []
    other_metal_atom_coord = None

    # get array of cofactor atoms
    cofactor_atoms = extract_cofactor_atoms(pdb_file, cofactor_index, reactive_residue_name)
    cofactor_coord_array = np.array([atom.coord for atom in cofactor_atoms])

    for atom in cofactor_atoms:
        if atom.element == 'MO' or atom.element == 'V1':
            other_metal_atom_coord = np.array([atom.coord]) # store the coordinates of non-iron metal atoms

    print(cofactor_coord_array)
    # now iterate through all atoms and get the ones within the selected radius of the cofactor atoms
    all_atoms = extract_all_atoms(pdb_file)
    for atom in all_atoms:
        print(atom.coord)
        squared_distances = np.sum((cofactor_coord_array - atom.coord) ** 2, axis=1)
        print(np.sqrt(squared_distances))
        print()
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
            element = line.split()[-1]  # Extract the element type
            x = float(line.split()[-6])  # Extract the x coordinate
            y = float(line.split()[-5])  # Extract the y coordinate
            z = float(line.split()[-4])  # Extract the z coordinate
            coord = np.array([x, y, z])
            all_atoms.append(Atom(element, coord))

    return all_atoms
