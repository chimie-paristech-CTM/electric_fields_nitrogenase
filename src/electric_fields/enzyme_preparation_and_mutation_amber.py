import os
import pandas as pd
import numpy as np
import re
import subprocess


def prepare_and_mutate_enzyme(enzyme_name, element_active_site, prot_state_his='HIP'):
    """
    """
    input_file = os.path.join(enzyme_name, f'{enzyme_name}.pdb')
    output_base_name =  os.path.join(enzyme_name, f'{enzyme_name}')

    # step 1: remove all the unnecessary lines and extract the name of the reactive cofactor
    reactive_residue_name = initial_clean_pdb(input_file, f'{output_base_name}_1.pdb')
 
    # step 2: mutate the CYS residues directly connected to the reactive cofactor
    df_coord_fe, pdb_protein_matrix, pdb_cofactor, pdb_mutated_enzyme = mutate_enzyme(f'{output_base_name}_1.pdb', 
                        f'{output_base_name}_2.pdb', enzyme_name, reactive_residue_name, element_active_site) 
    
    # step 3: correct incorrect ASP/ASN residues to GLU/GLN (if CD is present) 
    correct_asp_and_ash(f'{output_base_name}_2.pdb', f'{output_base_name}_3.pdb')

    # step 4: determine protonation states with propka and assign the right states in the original file (remove duplicate side chains)
    protonated_pdb = protonate_enzyme(f'{output_base_name}_3.pdb', f'{output_base_name}_4.pdb')

    # step 5: correct protonation for the HIS that need correcting
    correct_protonation_enzyme(f'{output_base_name}_4.pdb', f'{output_base_name}_5.pdb', pdb_protein_matrix, pdb_cofactor, prot_state_his)

    # step 6: undo rotations of the side chains                          
    undo_rotations(f'{output_base_name}_1.pdb', f'{output_base_name}_5.pdb', f'{output_base_name}_6.pdb')

    # step 7: add hydrogen atoms back to the main chain with tleap
    add_hydrogens_with_tleap(enzyme_name, f'{enzyme_name}_6.pdb', f'{enzyme_name}_7.pdb')

    # step 8: reassign chain indices
    assign_chains_on_reset(f'{output_base_name}_7.pdb', f'{output_base_name}_final.pdb')

    return reactive_residue_name, df_coord_fe, pdb_cofactor


def initial_clean_pdb(input_file, output_file):
    """
    Extracts ATOM and HETATM lines from a cleaned PDB file and writes them to separate output files.
    When alternative side-chain locations are provided, ignore the first one
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if re.match(r'ATOM', line):
                if len(line.split()[3]) == 4 and line.split()[3][0] == 'B':
                    continue
                # sometimes entry 3 and 4 are fused
                elif len(line.split()[2]) > 4 and 'B' in line.split()[2] and line.startswith('ATOM'):
                    continue 
                else:
                    outfile.write(line)
            elif re.match(r'HETATM', line):
                outfile.write(line) 

    # get reactive residue name
    reactive_cofactor_names = set([name for name in find_reactive_cofactor(input_file) if 'MO' not in name])
    reactive_residue_name = list(reactive_cofactor_names)[0].strip() # For most of the selected enzymes, this appears OK, but always double check

    return reactive_residue_name


# Histidine next to Mo should be HIE, the histidine close to S should be HIP (or HiD/HIE) -> there is no consensus!
def correct_protonation_enzyme(input_file, output_file, protein_matrix_csv, cofactor_csv, prot_state_his='HIP'):
    """
    Corrects the protonation states of histidine residues in an enzyme based on proximity to specific cofactors.

    This function identifies histidine residues that should be protonated as HIE (protonated on epsilon nitrogen) or HIP 
    (protonated on both nitrogens) based on their distances to molybdenum (MO) and sulfur (S) atoms in the cofactor. 
    It modifies the input file accordingly and writes the corrected structure to the output file.

    Args:
        input_file (str): Path to the input file containing the enzyme structure (e.g., PDB format).
        output_file (str): Path to the output file where the corrected enzyme structure will be saved.
        protein_matrix_csv (str): Path to a CSV file containing the protein matrix data (e.g., residue and atom coordinates).
        cofactor_csv (str): Path to a CSV file containing cofactor information (e.g., element types and coordinates).

    Returns:
        df_n_his: Dataframe containing the nitrogens with an uncertain protonation_state.
    """
    his_n = get_histidine_nitrogens_from_pdb_df(protein_matrix_csv)
    df_mo = get_specific_element_from_pdb_df(cofactor_csv, 'MO') #'V') 
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

    print(f'RESIDUE TO CHANGE TO HIE: {change_to_hie}, RESIDUE WITH UNCERTAIN PROTONATION STATE: {change_to_hip}')
    change_protonation_states(change_to_hie, change_to_hip, input_file, output_file, prot_state_his)


def get_histidine_nitrogens_from_pdb_df(pdb_file):
    """
    Extracts nitrogen atoms from histidine residues in a PDB DataFrame.

    This function filters the input DataFrame to identify nitrogen atoms (N) that belong to histidine residues. 
    It recognizes histidine variants including HIS, HID, HIE, and HIP.

    Args:
        pdb_file (pd.DataFrame): A pandas DataFrame representing a PDB file. The DataFrame must include columns 
                                 "Atom" (atom type) and "AA" (amino acid type).

    Returns:
        pd.DataFrame: A DataFrame containing rows corresponding to nitrogen atoms in histidine residues.
    """
    histidine_labels = ["HIS", "HID", "HIE", "HIP"]
    his_n = pdb_file[(pdb_file["Atom"] == "N") & (pdb_file["AA"].isin(histidine_labels))]
    return his_n
 

def change_protonation_states(change_to_hie, change_to_hip, input_file, output_file, prot_state_his='HIP'):
    """    
    """
    processed_lines = []

    # Read the PDB file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Extract relevant fields from the PDB line
            atom_name = line[12:16].strip()
            chain_id = line[21].strip()
            residue_seq = int(line[22:26].strip())
            
            # Create a unique identifier for the residue
            residue_identifier = (chain_id, residue_seq)

            if residue_identifier in change_to_hie: 
                new_line = line[:17] + "HIE" + line[20:]
                processed_lines.append(new_line)
            elif residue_identifier in change_to_hip:
                new_line = line[:17] + prot_state_his + line[20:]
                processed_lines.append(new_line)
            else:
                processed_lines.append(line)

    with open(output_file, 'w') as file:
        for line in processed_lines:
            file.write(line)


def undo_rotations(initial_file, processed_file, output_file):
    """    
    Restores the original atomic coordinates for non-hydrogen atoms in residues that have undergone rotation 
    during processing, while preserving the processed file's structure.

    This function compares two PDB files: the initial file and the processed file. It identifies residues that 
    have undergone rotation by comparing atom information (excluding hydrogen atoms) and restores the original 
    coordinates from the initial file for matching atoms.

    Args:
        initial_file (str): Path to the original PDB file before any processing.
        processed_file (str): Path to the processed PDB file where residues may have been rotated.
        output_file (str): Path to the output PDB file where restored atomic coordinates will be saved.

    Returns:
        None: Writes the corrected PDB structure to the output file.
    """
    # first, you make a dict with all the residue identifiers as keys, and a list as value. This list should contain a dictionary in its own right,
    # with index within that residue, element type etc. as key, and the lines themselves as the value (only non-H atoms).
    # once you have this dictionary for both files, you determine which ones have changed, and then double check that order is the same -> replace coordinates if needed
    res_dict_orig = {}
    res_dict_new = {}
    rotated_residues = set()

    with open(initial_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        try:
            if line.split()[-1] != 'H' and (line.startswith('HETATM') or line.startswith('ATOM')):
                residue_identifier = get_residue_identifier(line)
                if residue_identifier in res_dict_orig:
                    index = len(res_dict_orig[residue_identifier].keys())
                    res_dict_orig[residue_identifier][index] = (get_atom_info(line), line)
                else:
                    res_dict_orig[residue_identifier] = {0: (get_atom_info(line), line)}
        except:
            continue

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
    # Create a unique identifier for the atoms within the residue
    atom_info = (element_type, coord)

    return atom_info  


# for the enzymes analyzed, it looks like we are only dealing with ligating CYS residues, but for safety, we are keeping other options here as well
def mutate_enzyme(input_file, output_file, enzyme_name, reactive_residue_name, element_active_site):
    """
    Identifies and mutates specific residues in an enzyme structure, and extracts Fe atom coordinates 
    associated with the mutated residues.

    This function identifies residues in the enzyme structure that need to be mutated based on proximity 
    to a specified active site element. It performs the mutations, writes the modified structure to an output 
    file, and returns relevant structural data.

    Args:
        input_file (str): Path to the input PDB file containing the enzyme structure.
        output_file (str): Path to the output PDB file where the modified enzyme structure will be saved.
        enzyme_name (str): Name of the enzyme to be processed.
        reactive_residue_name (str): Name of the reactive residue (e.g., "CYS") in the active site.
        element_active_site (str): Element (e.g., "Fe") in the active site that dictates mutation criteria.

    Returns:
        tuple: A tuple containing:
            - pd.DataFrame: Coordinates of Fe atoms associated with the mutated residues.
            - pd.DataFrame: Protein matrix data extracted from the input file.
            - pd.DataFrame: Cofactor data extracted from the input file.
            - str: Path to the output file containing the mutated enzyme structure.
    """
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
    """    
    Corrects the residue names of problematic ASP and ASH residues in a PDB file.

    This function searches for aspartic acid residues (ASP or ASH) that have a 'CD' atom, 
    indicating a problematic protonation state or incorrect residue naming. It replaces these residues with 
    glutamic acid (GLU) in the output PDB file.

    Args:
        input_file (str): Path to the input PDB file containing the enzyme structure.
        output_file (str): Path to the output PDB file where the modified structure will be saved.

    Returns:
        None: Writes the corrected PDB structure to the output file.
    """
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
    Extracts the protein matrix from a PDB file and saves it to a CSV file.

    This function reads a PDB file, filters out only the ATOM entries, and extracts the coordinates 
    and other relevant information into a DataFrame. It then saves the extracted protein matrix 
    into a CSV file named after the enzyme.

    Args:
        file_name (str): Path to the input PDB file.
        enzyme_name (str): Name of the enzyme, used to save the extracted matrix into a CSV file.

    Returns:
        pd.DataFrame: A DataFrame containing the protein matrix, including atom coordinates and other information.
    """
    column_widths = [6, 5, 6, 4, 1, 8, 8, 8, 8, 6, 17, 2]
    column_names = ['Structure','site','type','AA','chain_id','position','x','y','z','n','m','Atom']
    pdb = pd.read_fwf(file_name, colspecs=[(sum(column_widths[:i]), sum(column_widths[:i+1])) for i in range(len(column_widths))], names=column_names, skiprows=0)
    pdb_protein_matrix = pdb.query('Structure == "ATOM"')
    pdb_protein_matrix['x'] = pdb_protein_matrix['x'].astype(float)
    pdb_protein_matrix['y'] = pdb_protein_matrix['y'].astype(float)
    pdb_protein_matrix['z'] = pdb_protein_matrix['z'].astype(float)
    pdb_protein_matrix.to_csv(f'{enzyme_name}/{enzyme_name}_protein_matrix.csv', index=False)

    return pdb_protein_matrix 


def extract_cofactor(file_name, enzyme_name, reactive_residue_name):
    """
    Extracts cofactor information from a PDB file for a specified reactive residue.

    This function reads a PDB file, filters the entries corresponding to the specified reactive residue 
    (e.g., a cofactor or metal ion), and extracts the relevant atomic information. The cofactor data is 
    then saved to a CSV file named after the enzyme.

    Args:
        file_name (str): Path to the input PDB file.
        enzyme_name (str): Name of the enzyme, used to save the extracted cofactor data into a CSV file.
        reactive_residue_name (str): The name of the reactive residue (e.g., "CYS", "MO") to extract from the PDB.
    """
    column_widths = [6, 5, 6, 4, 1, 8, 8, 8, 8, 6, 16, 3]
    column_names = ['Structure','site','type','AA','chain_id','position','x','y','z','n','m','Atom']
    pdb = pd.read_fwf(file_name, colspecs=[(sum(column_widths[:i]), sum(column_widths[:i+1])) for i in range(len(column_widths))], names=column_names, skiprows=0)
    pdb_cofactor = pdb.query(f'AA == "{reactive_residue_name.strip()}"')
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
                atom_type = line[11:17].strip()
                if atom_type.startswith('MO') or atom_type == 'V1': #or atom_type == '\bFe\w*' or atom_type.startswith('FE'):
                    selected_residues.append(line[17:21])    
    return selected_residues


def protonate_enzyme(input_file, output_pdb, ph=7):
    """
    Args:
        input_file (_type_): _description_
    """
    output_pqr = f"{output_pdb.strip('.pdb')}.pqr"

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

def add_hydrogens_with_tleap(enzyme_name: str, input_pdb: str, output_pdb: str, leaprc: str = "leaprc.protein.ff14SB"):
    """
    """
    home_dir = os.getcwd()
    os.chdir(enzyme_name)

    if not os.path.isfile(input_pdb):
        raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")

    tleap_input = "tleap.in"

    # Create tleap input script
    with open(tleap_input, "w") as f:
        f.write(f"source {leaprc}\n")
        f.write(f"enzyme = loadPdb \"{input_pdb}\"\n")
        f.write(f"savePdb enzyme \"{output_pdb}\"\n")
        f.write("quit\n")

    try:
        # Run tleap in a subprocess
        subprocess.run(["tleap", "-f", tleap_input], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running tleap: {e}")
    finally:
        # Clean up tleap input file
        if os.path.isfile(tleap_input):
            os.remove(tleap_input)

    os.chdir(home_dir)


def assign_chains_on_reset(input_pdb, output_pdb):
    """
    Assign chain identifiers dynamically based on residue number resets in a PDB file.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to the output PDB file.
    """
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"  # Chain ID pool
    current_chain_index = 0

    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if 'TER' in line:
                current_chain_index += 1
                if current_chain_index >= len(chain_ids):
                    raise ValueError("Exceeded available chain identifiers (A-Z).") 
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract residue number
                res_num = int(line[22:26].strip())

                # Assign the current chain ID
                chain_id = chain_ids[current_chain_index]

                # Modify the chain ID (column 22) and write the updated line
                new_line = line[:21] + chain_id + line[22:]
                outfile.write(new_line)

                # Update the last residue number
                last_residue_number = res_num
