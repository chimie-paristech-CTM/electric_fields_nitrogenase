import os
import pandas as pd
import numpy as np
import re
import subprocess


def prepare_and_mutate_enzyme(enzyme_name, element_active_site):
    """_summary_

    Args:
        enzyme_name (_type_): _description_
        reactive_residue_name (_type_): _description_
        element_active_site (_type_): _description_

    Returns:
        _type_: _description_
    """
    # remove all the unnecessary lines and extract the name of the reactive cofactor
    reactive_residue_name = initial_clean_pdb(enzyme_name)
    
    # mutate the CYS residues directly connected to the reactive cofactor
    df_coord_fe, pdb_protein_matrix, pdb_cofactor, pdb_mutated_enzyme = mutate_enzyme(enzyme_name, reactive_residue_name, element_active_site) 
    # correct incorrect ASP residues to GLU (if CD is present) 
    correct_asp_and_ash(os.path.join(enzyme_name, f'{enzyme_name}_mutated.pdb'), os.path.join(enzyme_name, f'{enzyme_name}_mutated_corr.pdb')) 

    # TODO: remove alternative locations of side-chains


    # determine protonation states with propka and assign the right states in the original file (remove duplicate side chains)
    # + correct protonation for the HIS that need correcting
    protonated_pdb = protonate_enzyme(enzyme_name, os.path.join(enzyme_name, f'{enzyme_name}_mutated_corr.pdb'))
    correct_protonation_enzyme(pdb_protein_matrix, pdb_cofactor, protonated_pdb, 
                                os.path.join(enzyme_name, f'{enzyme_name}_mutated_protonated_corr2.pdb'))
    undo_rotations(os.path.join(enzyme_name, f'{enzyme_name}.pdb'), os.path.join(enzyme_name, f'{enzyme_name}_mutated_protonated_corr2.pdb'),
            os.path.join(enzyme_name, f'{enzyme_name}_final.pdb'))

    # TODO: Add hydrogen atoms to the final file with the help of TLEAP and save final geometry





    #protonated_pdb = protonate_enzyme(enzyme_name, pdb_mutated_enzyme)
    #assign_chains_on_reset(protonated_pdb, os.path.join(enzyme_name, f'{enzyme_name}_mutated_protonated_with_chains.pdb'))
    #correct_protonation_enzyme(pdb_protein_matrix, pdb_cofactor, os.path.join(enzyme_name, f'{enzyme_name}_mutated_protonated_with_chains.pdb'), 
    #                            os.path.join(enzyme_name, f'{enzyme_name}_mutated_protonated_with_chains_corr.pdb'))
    #extract_atoms_pdb_file(os.path.join(enzyme_name, f'{enzyme_name}_mutated_protonated_with_chains_corr.pdb'), 
    #                       os.path.join(enzyme_name, f'{enzyme_name}_final.pdb'))
    #extract_protein_matrix_pdb_file(protonated_pdb, os.path.join(enzyme_name, f'{enzyme_name}_final.pdb'))

    return reactive_residue_name, df_coord_fe, pdb_cofactor


def initial_clean_pdb(enzyme_name):
    """
    Extracts ATOM and HETATM lines from a cleaned PDB file and writes them to separate output files.
    """
    input_file = os.path.join(enzyme_name, f'{enzyme_name}.pdb')
    output_file =  os.path.join(enzyme_name, f'{enzyme_name}_clean.pdb')

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if re.match(r'ATOM', line):
                outfile.write(line)
            elif re.match(r'HETATM', line):
                outfile.write(line) 

    # get reactive residue name
    reactive_cofactor_names = set(find_reactive_cofactor(input_file))
    reactive_residue_name = list(reactive_cofactor_names)[0] # For most of the selected enzymes, this appears OK, but always double check

    return reactive_residue_name

# Histidine next to Mo should be HIE, the histidine close to S should be HIP
def correct_protonation_enzyme(protein_matrix_csv, cofactor_csv, input_file, output_file):
    """
    Args:
        enzyme_name (_type_): _description_
        reactive_residue_name (_type_): _description_

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    his_n = get_histidine_nitrogens_from_pdb_df(protein_matrix_csv)
    df_mo = get_specific_element_from_pdb_df(cofactor_csv, 'MO') 
    df_s = get_specific_element_from_pdb_df(cofactor_csv, 'S')

    change_to_hie = []
    for i,row1 in his_n.iterrows():
        for j,row2 in df_mo.iterrows():
            dist = np.sqrt((row1['x'] - row2['x'])**2 + (row1['y'] - row2['y'])**2 +(row1['z'] - row2['z']) ** 2)
            if dist <= 2.5:
                change_to_hie.append((row1['chain_id'], int(row1['position'])))
    
    change_to_hip = []
    for i,row1 in his_n.iterrows():
        for j,row2 in df_s.iterrows():
            dist = np.sqrt((row1['x'] - row2['x'])**2 + (row1['y'] - row2['y'])**2 +(row1['z'] - row2['z']) ** 2)
            if dist <= 4:
                if (row1['chain_id'], int(row1['position'])) not in change_to_hie:
                    change_to_hip.append((row1['chain_id'], int(row1['position'])))    

    print(change_to_hie, change_to_hip)
    change_protonation_states(change_to_hie, change_to_hip, input_file, output_file)


def get_histidine_nitrogens_from_pdb_df(pdb_file):
    """
    """
    histidine_labels = ["HIS", "HID", "HIE", "HIP"]
    his_n = pdb_file[(pdb_file["Atom"] == "N") & (pdb_file["AA"].isin(histidine_labels))]
    return his_n
 
def change_protonation_states(change_to_hie, change_to_hip, input_file, output_file):
    processed_lines = []

    # Read the PDB file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Extract relevant fields from the PDB line
            atom_name = line[12:16].strip()
            #residue_name = line[17:20].strip()
            chain_id = line[21].strip()
            residue_seq = int(line[22:26].strip())
            
            # Create a unique identifier for the residue
            residue_identifier = (chain_id, residue_seq)

            if residue_identifier in change_to_hie: 
                new_line = line[:17] + "HIE" + line[20:]
                processed_lines.append(new_line)
            elif residue_identifier in change_to_hip:
                new_line = line[:17] + "HIP" + line[20:]
                processed_lines.append(new_line)
            else:
                processed_lines.append(line)

    with open(output_file, 'w') as file:
        for line in processed_lines:
            file.write(line)


def undo_rotations(initial_file, processed_file, output_file):
    # first, you make a dict with all the residue identifiers as keys, and a list as value. This list should contain a dictionary in its own right,
    # with index within that residue, element type etc. as key, and the lines themselves as the value (only non-H atoms).
    # once you have this dictionary for both files, you determine which ones have changed, and then double check that order is the same -> replace coordinates if needed
    res_dict_orig = {}
    res_dict_new = {}
    rotated_residues = set()

    with open(initial_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        #try:
            if line.split()[-1] != 'H' and (line.startswith('HETATM') or line.startswith('ATOM')):
                residue_identifier = get_residue_identifier(line)
                if residue_identifier in res_dict_orig:
                    index = len(res_dict_orig[residue_identifier].keys())
                    res_dict_orig[residue_identifier][index] = (get_atom_info(line), line)
                else:
                    res_dict_orig[residue_identifier] = {0: (get_atom_info(line), line)}
        #except:
        #    continue

    with open(processed_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        try:
            if line.split()[-1] != 'H' and (line.startswith('HETATM') or line.startswith('ATOM')):
                residue_identifier = get_residue_identifier(line)
                if residue_identifier in res_dict_new:
                    index = len(res_dict_new[residue_identifier].keys())
                    res_dict_new[residue_identifier][index] = (get_atom_info(line), line)
                else:
                    res_dict_new[residue_identifier] = {0: (get_atom_info(line), line)}
        except:
            continue

    sorted_residues = sorted(res_dict_new.keys(), key=lambda x: (x[0], int(x[1])))

    with open(initial_file, 'r') as in_file, open(output_file, 'w') as out_file:
        for residue in sorted_residues:
            try:
                for i in range(len(res_dict_new[residue])):
                    if res_dict_new[residue][i][0] == res_dict_orig[residue][i][0]:
                        out_file.write(res_dict_new[residue][i][1])
                    else:
                        modified_line = modify_line(res_dict_new[residue][i][1], res_dict_orig[residue][i][1])
                        out_file.write(modified_line)
            except KeyError:
                out_file.write(res_dict_new[residue][i][1])


def modify_line(new_line, old_line):
    modified_line = new_line[:30] + old_line[30:54] + new_line[54:]
    return modified_line


def get_residue_identifier(line):
    chain_id = line[21].strip()
    residue_seq = line[22:26].strip() 
    # Create a unique identifier for the residue
    residue_identifier = (chain_id, residue_seq)

    return residue_identifier


def get_atom_info(line):
    coord = line[30:54] 
    element_type = line.split()[-1]
    atom_info = (element_type, coord)

    return atom_info  


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
    input_file = os.path.join(enzyme_name, f'{enzyme_name}_clean.pdb') #'test_with_hydrogens.pdb' #f'{enzyme_name}_clean.pdb'
    output_file = os.path.join(enzyme_name, f'{enzyme_name}_mutated.pdb')

    # extract protein matrix & cofactor
    pdb_protein_matrix = extract_protein_matrix(input_file, enzyme_name)
    pdb_cofactor = extract_cofactor(input_file, enzyme_name, reactive_residue_name=reactive_residue_name)

    residue_list_final = []
    df_coor_fe_list = []

    # find mutations to be made
    for atom_type in ["S", "O"]:
        residue_list, df_coor_fe = find_residues_to_be_mutated(pdb_protein_matrix, atom_type, pdb_cofactor, element_active_site)
        residue_list_final += residue_list
        df_coor_fe_list.append(df_coor_fe)

    df_coor_fe_final = pd.concat(df_coor_fe_list)

    print(f'Residues to be mutated: {residue_list_final}')

    mutate_sites(input_file, output_file, residue_list_final)

    # return the coordinates of the Fe atoms associated with mutated CYS residues
    return df_coor_fe_final, pdb_protein_matrix, pdb_cofactor, output_file
 

def mutate_sites(input_file, output_file, residue_list_final):
    """
    Processes a PDB file by retaining only the atoms in specified residues 
    that match those in alanine (ALA) and renaming the residue types to ALA.

    Args:
        input_file (str): Path to the input PDB file.
        output_file (str): Path to the output PDB file.
        residue_list_final (list): List of residue identifiers to process.

    Returns:
        list: Processed lines of the PDB file.
    """
    # Define the standard atom names for alanine (ALA)
    ala_atoms = {"N", "CA", "C", "O", "CB"}

    residue_identifiers_to_modify = [(sublist[2], sublist[0]) for sublist in residue_list_final]

    processed_lines = []

    # Read the PDB file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.startswith("ATOM"):
            # Extract relevant fields from the PDB line
            atom_name = line[12:16].strip()
            #residue_name = line[17:20].strip()
            chain_id = line[21].strip()
            residue_seq = int(line[22:26].strip())
            
            # Create a unique identifier for the residue
            residue_identifier = (chain_id, residue_seq)

            # Check if the residue is in the final list and if the atom matches ALA
            if residue_identifier in residue_identifiers_to_modify: 
                if atom_name in ala_atoms:
                    # Modify the residue name to ALA
                    new_line = line[:17] + "ALA" + line[20:]
                    processed_lines.append(new_line)
                else:
                    continue
            else:
                processed_lines.append(line)

    with open(output_file, 'w') as file:
        for line in processed_lines:
            file.write(line)


def correct_asp_and_ash(input_file, output_file):

    processed_lines = []

    with open(input_file, 'r') as file:
        lines = file.readlines()

    # first determine problematic residues
    asp_and_ash_res_to_change = set()

    for line in lines:
        if line.startswith("ATOM"):
            # Extract relevant fields from the PDB line
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21].strip()
            residue_seq = int(line[22:26].strip())
            
            if atom_name == 'CD' and (residue_name == 'ASP' or residue_name == 'ASH'):
                asp_and_ash_res_to_change.add((chain_id, residue_seq))

    for line in lines:
        if line.startswith("ATOM"):
            # Extract relevant fields from the PDB line
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21].strip()
            residue_seq = int(line[22:26].strip())
            
            if (chain_id, residue_seq) in asp_and_ash_res_to_change:
                new_line = line[:17] + "GLU" + line[20:]
                processed_lines.append(new_line)
            else:
                processed_lines.append(line)

    with open(output_file, 'w') as file:
        for line in processed_lines:
            file.write(line)      


def extract_protein_matrix(file_name, enzyme_name):
    """
    """
    pdb = pd.read_csv(file_name, on_bad_lines='skip', encoding = 'utf-8', sep="\s+", low_memory=False, header=None)
    pdb.columns = ['Structure','site','type','AA','chain_id','position','x','y','z','n','m','Atom']
    pdb_protein_matrix = pdb.query('Structure == "ATOM"')
    pdb_protein_matrix['x'] = pdb_protein_matrix['x'].astype(float)
    pdb_protein_matrix['y'] = pdb_protein_matrix['y'].astype(float)
    pdb_protein_matrix['z'] = pdb_protein_matrix['z'].astype(float)
    pdb_protein_matrix.to_csv(f'{enzyme_name}/{enzyme_name}_protein_matrix.csv', index=False)

    return pdb_protein_matrix 


def extract_cofactor(file_name, enzyme_name, reactive_residue_name):
    """
    """
    pdb = pd.read_csv(file_name, on_bad_lines='skip', encoding = 'utf-8', sep="\s+", low_memory=False, header=None)
    pdb.columns = ['Structure','site','type','AA','chain_id','position','x','y','z','n','m','Atom']
    pdb_cofactor = pdb.query(f'type == "{reactive_residue_name}"').drop(['m', 'Atom'], axis=1)
    pdb_cofactor.columns = ['Structure','site','type','AA','x','y','z','n','m','Atom']
    pdb_cofactor['x'] = pdb_cofactor['x'].astype(float)
    pdb_cofactor['y'] = pdb_cofactor['y'].astype(float)
    pdb_cofactor['z'] = pdb_cofactor['z'].astype(float)
    pdb_cofactor.to_csv(f'{enzyme_name}/{enzyme_name}_cofactor.csv', index=False)

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
    df_fe = get_specific_element_from_pdb_df(pdb_cofactor, element_active_site.upper())

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

def protonate_enzyme(enzyme_name, input_file, ph=7):
    """
    Args:
        input_file (_type_): _description_
    """
    output_pqr = os.path.join(enzyme_name, f'{enzyme_name}_mutated_corr_protonated.pqr')
    output_pdb = os.path.join(enzyme_name, f'{enzyme_name}_mutated_corr_protonated.pdb')

    try:
        # Step 1: Run pdb2pqr
        pdb2pqr_command = [
            "pdb2pqr",
            "--ff=AMBER",
            input_file,
            output_pqr,
            "--titration-state-method", "propka",
            f"--with-ph={ph}",
            "--ffout=AMBER",
            "--drop-water",
            "--keep-chain"
        ]
        print(f"Running pdb2pqr command: {' '.join(pdb2pqr_command)}")
        subprocess.run(pdb2pqr_command, check=True)
        
        # Step 2: Run obabel
        obabel_command = [
            "obabel",
            "-ipqr", output_pqr,
            "-opdb", "-O", output_pdb
        ]
        print(f"Running obabel command: {' '.join(obabel_command)}")
        subprocess.run(obabel_command, check=True)

        print("Commands executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing commands: {e}")

    
    return output_pdb


def assign_chains_on_reset(input_pdb, output_pdb):
    """
    Assign chain identifiers dynamically based on residue number resets in a PDB file.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to the output PDB file.
    """
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"  # Chain ID pool
    current_chain_index = 0
    last_residue_number = None

    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract residue number
                res_num = int(line[22:26].strip())

                # Detect residue number reset and increment chain ID
                if last_residue_number is not None and res_num < last_residue_number:
                    current_chain_index += 1
                    if current_chain_index >= len(chain_ids):
                        raise ValueError("Exceeded available chain identifiers (A-Z).")

                # Assign the current chain ID
                chain_id = chain_ids[current_chain_index]

                # Modify the chain ID (column 22) and write the updated line
                new_line = line[:21] + chain_id + line[22:]
                outfile.write(new_line)

                # Update the last residue number
                last_residue_number = res_num
            else:
                # Write non-ATOM/HETATM lines unchanged
                outfile.write(line)


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


def extract_protein_matrix_pdb_file(input_file, output_file):
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
