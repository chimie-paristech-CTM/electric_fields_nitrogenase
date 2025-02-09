import os
import pandas as pd
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import re


def prepare_and_mutate_enzyme(enzyme_name, element_active_site):
    """_summary_

    Args:
        enzyme_name (_type_): _description_
        reactive_residue_name (_type_): _description_
        element_active_site (_type_): _description_

    Returns:
        _type_: _description_
    """
    reactive_residue_name = preliminary_steps(enzyme_name)
    df_coord_fe, pdb_cofactor = mutate_enzyme(enzyme_name, reactive_residue_name, element_active_site)

    return reactive_residue_name, df_coord_fe, pdb_cofactor


def preliminary_steps(enzyme_name):
    """
    Performs preliminary clean-up steps for an enzyme PDB file and identifies the reactive residue.

    Parameters:
    enzyme_name (str): The name of the enzyme, used to construct file paths.

    Returns:
    str: The name of the reactive residue identified in the enzyme PDB file.

    Workflow:
    1. Sets the input and output file names based on the enzyme name.
    2. Cleans up the input PDB file and saves the cleaned structure to the output file.
    3. Identifies and returns the name of the reactive residue.
    """
    # set file names
    input_file = f'{enzyme_name}.pdb'
    output_file = f'{enzyme_name}_clean.pdb'

    # do a preliminary clean-up
    output_file = clean_up_pdb_file(enzyme_name, input_file, output_file)

    # get reactive residue name
    reactive_cofactor_names = set(find_reactive_cofactor(os.path.join(enzyme_name, input_file)))
    print(reactive_cofactor_names)
    reactive_residue_name = list(reactive_cofactor_names)[0] # For most of the selected enzymes, this appears OK, but always double check

    return reactive_residue_name


# for the enzymes analyzed, it looks like we are only dealing with ligating CYS residues, but for safety, we are keeping other options here as well
def mutate_enzyme(enzyme_name, reactive_residue_name, element_active_site):
    """
    Mutate specific residues in an enzyme to alanine and return the coordinates of Fe atoms associated with the mutated residues.

    This function performs the following steps:
    1. Extracts the protein matrix and cofactor from a PDB file of the enzyme.
    2. Identifies residues of specified atom types ("S" and "O") to be mutated to alanine.
    3. Applies the mutations to generate a new PDB file with the specified mutations.
    4. Returns the coordinates of Fe atoms associated with the mutated residues and the cofactor.

    Parameters:
    enzyme_name (str): The name of the enzyme, used to locate PDB files and output files.
    reactive_residue_name (str): The name of the reactive residue to be mutated.
    element_active_site (str): The element active site to focus the mutation process.

    Returns:
    df_coor_fe_final (pd.DataFrame): DataFrame containing the coordinates of Fe atoms associated with the mutated residues.
    pdb_cofactor (object): The cofactor extracted from the enzyme.
    """
    # extract protein matrix & cofactor
    extract_atoms_pdb_file(os.path.join(enzyme_name, f'{enzyme_name}_clean.pdb'), os.path.join(enzyme_name, f'{enzyme_name}_atoms.pdb'))
    pdb_protein_matrix = extract_protein_matrix(enzyme_name, enzyme_name)
    pdb_cofactor = extract_cofactor(enzyme_name, enzyme_name, reactive_residue_name=reactive_residue_name)

    residue_list_final = []
    df_coor_fe_list = []

    # find mutations to be made
    for atom_type in ["S", "O"]:
        residue_list, df_coor_fe = find_residues_to_be_mutated(pdb_protein_matrix, atom_type, pdb_cofactor, element_active_site)
        residue_list_final += residue_list
        df_coor_fe_list.append(df_coor_fe)

    df_coor_fe_final = pd.concat(df_coor_fe_list)

    print(f'Residues to be mutated: {residue_list_final}')

    input_file = f'{enzyme_name}_clean.pdb' #'test_with_hydrogens.pdb' #f'{enzyme_name}_clean.pdb'
    output_file = f'{enzyme_name}_mutated.pdb'

    # perform the mutations with the PDBFixer
    fixer = PDBFixer(os.path.join(enzyme_name, input_file))
    mutations = {}

    for residue_id, residue_type, chain_id in residue_list_final:
        if chain_id not in mutations.keys():
            mutations[chain_id] = [f'{residue_type}-{residue_id}-ALA']
        else:
            mutations[chain_id].append(f'{residue_type}-{residue_id}-ALA')

    for chain_id in mutations.keys():
        fixer.applyMutations(mutations[chain_id], chain_id)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(enzyme_name, output_file), 'w')) 

    extract_atoms_pdb_file(os.path.join(enzyme_name, f'{enzyme_name}_mutated.pdb'), os.path.join(enzyme_name, f'{enzyme_name}_mutated_atoms.pdb')) 

    # return the coordinates of the Fe atoms associated with mutated CYS residues
    return df_coor_fe_final, pdb_cofactor


def clean_up_pdb_file(path, input_file, output_file, chains_to_keep=None):
    """
    Cleans up a PDB file by performing several steps including removing nonstandard residues,
    adding missing atoms and hydrogens, and optionally retaining specified chains.

    Parameters:
    path (str): The directory path where the input and output files are located.
    input_file (str): The name of the input PDB file to be cleaned.
    output_file (str): The name of the output PDB file to save the cleaned structure.
    chains_to_keep (list of int, optional): List of chain indices to retain in the PDB file. If None, all chains are kept. Default is None.

    Returns:
    str: The name of the output file containing the cleaned PDB structure.

    Workflow:
    1. Initializes PDBFixer with the input PDB file.
    2. Identifies and removes missing residues at the beginning or end of chains.
    3. Replaces nonstandard residues with standard ones.
    4. Finds and adds missing atoms.
    5. Optionally removes chains not specified in `chains_to_keep`.
    6. Writes the cleaned structure to the output PDB file.
    """
    fixer = PDBFixer(os.path.join(path, input_file))
    fixer = find_missing_intermediate_residues(fixer)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    if chains_to_keep:
        nchains = len(list(fixer.topology.chains()))
        for index in range(nchains)[::-1]:
            if index not in chains_to_keep:
                try:
                    fixer.removeChains(range(index, index+1))
                except:
                    continue
    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(path, output_file), 'w'))

    return output_file


def find_missing_intermediate_residues(fixer):
    """
    Identifies and removes missing residues at the beginning or end of chains in a PDB structure, keeping only intermediate missing residues.

    Parameters:
    fixer (PDBFixer): An instance of PDBFixer containing the PDB structure to be processed.

    Returns:
    PDBFixer: The PDBFixer instance with updated missing residues, excluding those at the ends of chains.

    Workflow:
    1. Finds missing residues in the PDB structure using PDBFixer.
    2. Iterates through the identified missing residues.
    3. Removes entries corresponding to residues at the beginning or end of chains.
    4. Returns the PDBFixer instance with the filtered missing residues.
    """
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    keys_to_delete = []
    for key in keys:
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            keys_to_delete.append(key)
    for key in keys_to_delete:
        del fixer.missingResidues[key]

    return fixer


def find_reactive_cofactor(pdb_file):
    """
    Identifies reactive cofactors in a PDB file by searching for specific atom types.

    Parameters:
    pdb_file (str): The path to the PDB file to be searched.

    Returns:
    list: A list of residue names corresponding to the reactive cofactors identified in the PDB file.

    Workflow:
    1. Opens and reads the PDB file line by line.
    2. Identifies lines starting with 'HETATM'.
    3. Extracts and checks the atom type for specific reactive cofactors (e.g., atoms starting with 'MO' or 'V1').
    4. Appends the residue name to the list if a reactive cofactor is identified.
    5. Returns the list of identified reactive cofactor residue names.
    """
    selected_residues = []
    
    with open(pdb_file, 'r') as f:
        for line in f:     
            if line.startswith('HETATM'):
                atom_type = line.split()[1]
                if atom_type.startswith('MO') or atom_type == 'V1': #or atom_type == '\bFe\w*' or atom_type.startswith('FE'):
                    selected_residues.append(line.split()[2])     
    return selected_residues


def extract_atoms_pdb_file(input_file, output_file):
    """
    Extracts ATOM and selected HETATM lines from a cleaned PDB file and writes them to separate output files.

    Args:
        input_file (_type_): _description_
        output_file (_type_): _description_
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if re.match(r'ATOM', line):
                outfile.write(line)
            elif re.match(r'HETATM', line):
                residue_name = line[17:20].strip()
                if residue_name != 'HOH':
                    outfile.write(line) 


def extract_protein_matrix(path, enzyme_name):
    """
    Extracts the protein matrix from a PDB file and saves it as a CSV file.

    Parameters:
    path (str): The directory path where the input PDB file and output CSV file will be located.
    enzyme_name (str): The name of the enzyme, used to construct file paths.

    Returns:
    pandas.DataFrame: DataFrame containing the protein matrix extracted from the PDB file.

    Workflow:
    1. Reads the PDB file into a pandas DataFrame, skipping bad lines and specifying column names and separators.
    2. Filters the DataFrame to select rows where 'Structure' is 'ATOM', representing protein atoms.
    3. Converts coordinate columns ('x', 'y', 'z') to float type for numerical operations.
    4. Saves the filtered DataFrame as a CSV file named '{enzyme_name}_protein_matrix.csv' in the specified path.
    5. Returns the filtered DataFrame containing the protein matrix.
    """
    pdb = pd.read_csv(os.path.join(path, f'{enzyme_name}_atoms.pdb'), on_bad_lines='skip', 
                  encoding = 'utf-8', sep="\s+", low_memory=False, header=None)
    pdb.columns = ['Structure','site','type','AA','chain_id','position','x','y','z','n','m','Atom']
    pdb_protein_matrix = pdb.query('Structure == "ATOM"')
    pdb_protein_matrix['x'] = pdb_protein_matrix['x'].astype(float)
    pdb_protein_matrix['y'] = pdb_protein_matrix['y'].astype(float)
    pdb_protein_matrix['z'] = pdb_protein_matrix['z'].astype(float)
    pdb_protein_matrix.to_csv(os.path.join(path, f'{enzyme_name}_protein_matrix.csv'), index=False)

    return pdb_protein_matrix 


def extract_cofactor(path, enzyme_name, reactive_residue_name):
    """
    Extracts the coordinates and details of a specific cofactor residue from a PDB file and saves it as a CSV file.

    Parameters:
    path (str): The directory path where the input PDB file and output CSV file will be located.
    enzyme_name (str): The name of the enzyme, used to construct file paths.
    reactive_residue_name (str): The name of the reactive residue to extract from the PDB file.

    Returns:
    pandas.DataFrame: DataFrame containing the extracted cofactor residue information from the PDB file.

    Workflow:
    1. Reads the PDB file into a pandas DataFrame, skipping bad lines and specifying column names and separators.
    2. Filters the DataFrame to select rows where 'type' matches the provided 'reactive_residue_name'.
    3. Drops unnecessary columns ('n' and 'Atom') from the DataFrame.
    4. Converts coordinate columns ('x', 'y', 'z') to float type for numerical operations.
    5. Saves the filtered DataFrame as a CSV file named '{enzyme_name}_cofactor.csv' in the specified path.
    6. Returns the filtered DataFrame containing the cofactor residue information.
    """
    pdb = pd.read_csv(os.path.join(path, f'{enzyme_name}_atoms.pdb'), on_bad_lines='skip', 
                  encoding = 'utf-8', sep="\s+", low_memory=False, header=None)
    pdb.columns = ['Structure','site','type','AA','chain_id','position','x','y','z','n','m','Atom']
    pdb_cofactor = pdb.query(f'type == "{reactive_residue_name}"').drop(['n', 'Atom'], axis=1)

    pdb_cofactor.columns = ['Structure','type', 'AA', 'chain_id', 'id', 'x','y','z', 'n', 'Atom']

    pdb_cofactor['x'] = pdb_cofactor['x'].astype(float)
    pdb_cofactor['y'] = pdb_cofactor['y'].astype(float)
    pdb_cofactor['z'] = pdb_cofactor['z'].astype(float)
    pdb_cofactor.to_csv(os.path.join(path, f'{enzyme_name}_cofactor.csv'), index=False)

    return pdb_cofactor

def find_residues_to_be_mutated(pdb_protein_matrix, atom_type, pdb_cofactor, element_active_site):
    """
    Identify residues to be mutated based on proximity to a specific element in the active site and return their coordinates.

    This function performs the following steps:
    1. Extracts the coordinates of specific atom types from the protein matrix and the specified element from the cofactor.
    2. Calculates the distances between these atoms and the active site element.
    3. Identifies residues where the distance is less than or equal to 2.5 angstroms.
    4. Collects the coordinates of the active site element associated with these residues.

    Parameters:
    pdb_protein_matrix (pd.DataFrame): DataFrame containing the atomic coordinates and residue information of the protein.
    atom_type (str): The type of atom to consider for mutation (e.g., "S" for sulfur).
    pdb_cofactor (pd.DataFrame): DataFrame containing the atomic coordinates and residue information of the cofactor.
    element_active_site (str): The element in the active site to focus the mutation process (e.g., "Fe" for iron).

    Returns:
    res_to_mutate (list of lists): A list of residues to be mutated. Each entry is a list containing:
        - residue position (int)
        - residue type (str)
        - chain ID (str)
    df_coor_fe (pd.DataFrame): DataFrame containing the coordinates of the active site element associated with the residues to be mutated.
    """
    df_so = get_specific_element_from_pdb_df(pdb_protein_matrix, atom_type)
    df_fe = get_specific_element_from_pdb_df(pdb_cofactor, element_active_site)

    res_to_mutate = []
    coor_fe = {'x': [], 'y': [], 'z':[]}
    for i,row1 in df_so.iterrows():
        for j,row2 in df_fe.iterrows():
            dist = np.sqrt((row1['x'] - row2['x'])**2 + (row1['y'] - row2['y'])**2 +(row1['z'] - row2['z']) ** 2)
            if dist <= 2.5:
                res_to_mutate.append([int(row1['position']), row1['AA'], row1['chain_id']])
                coor_fe['x'].append (row2['x'])
                coor_fe['y'].append(row2['y'])
                coor_fe['z'].append(row2['z'])
    
    df_coor_fe = pd.DataFrame(coor_fe)

    return res_to_mutate, df_coor_fe


def get_specific_element_from_pdb_df(pdb_df, element):
    """
    Retrieves rows from a PDB DataFrame corresponding to a specific atom type.

    Parameters:
    pdb_df (pandas.DataFrame): DataFrame containing PDB data with columns including 'Atom' for atom types.
    element (str): The specific atom type to retrieve from the DataFrame.

    Returns:
    pandas.DataFrame: Subset of `pdb_df` containing rows where 'Atom' matches the specified `element`.
    """
    pdb_element = pdb_df.query(f'Atom == "{element}"')
    return pdb_element
