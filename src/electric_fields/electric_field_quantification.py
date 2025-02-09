import os
import numpy as np
import pandas as pd
from .quantification_txt import QuantificationTxt
from .quantification_amoeba import QuantificationAmoeba
from .nearby_residue_identification import get_residues_near_cofactor
from Bio import PDB

from scipy.spatial.distance import cdist

field_conversion_factor = 51.4220675112


def get_discarded_atom_indices(pdb_filename, combined_pdb_file, reactive_residue_name, list_residue_indices_to_discard=[], radius_to_exclude=None):
    """
    Returns the indices of atoms in the discarded amino acid residues from the PDB file.
    
    Args:
    pdb_filename (str): Path to the PDB file.
    combined_pdb_file (str): Path to the PDB file with protein matrix + cofactor.
    reactive_residue_name (str): The name of the reactive residue.
    list_residue_indices_to_discard (list): List of residue indices (integers) to discard.
    radius_to_exclude (float, optional): If specified, residues within this radius of the cofactor will be discarded.
    
    Returns:
    list: List of atom indices (starting from 0) corresponding to the discarded residues.
    """
    if radius_to_exclude != None:
        list_residue_indices_to_discard = get_residues_near_cofactor(combined_pdb_file, reactive_residue_name, radius=radius_to_exclude) 
    
    # Parse the PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_filename)

    discarded_atom_indices = []

    # Iterate over all residues in the structure
    atom_index = 0  # We'll track the index for each atom
    
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.get_id()[1]  # Extract numerical residue index
                if residue_id in list_residue_indices_to_discard:
                    # If residue name is in the list of residues to discard, add atom indices
                    for atom in residue:
                        discarded_atom_indices.append(atom_index)
                        atom_index += 1
                else:
                    # Otherwise, just increment the atom index
                    for atom in residue:
                        atom_index += 1

    return discarded_atom_indices


def extract_charges_protein_matrix_amber(enzyme_name, file_name, reactive_residue_name, combined_pdb_file, 
    list_aa_to_discard=[], radius_to_exclude=None):
    """
    """
    atom_data = []
    atom_section = False

    with open(f'{enzyme_name}/{file_name}', 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                continue
            elif line.startswith("@<TRIPOS>") and atom_section:
                break  # Stop reading if we encounter a new section
            
            if atom_section:
                parts = line.split()
                if len(parts) >= 9:  # Ensure it is an atom line
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    element_type = parts[5]
                    residue_id = int(parts[6])
                    residue_name = parts[7]
                    charge = float(parts[8])
                    
                    atom_data.append([atom_id, atom_name, x, y, z, element_type, residue_id, residue_name, charge])

    # Create a DataFrame
    columns = ["atom ID", "atom Name", "x", "y", "z", "element Type", "residue ID", "residue Name", "charge"]
    n = pd.DataFrame(atom_data, columns=columns)

    if len(list_aa_to_discard) != 0:
        n = n[~n['residue ID'].isin(list_aa_to_discard)]
        n.to_csv(os.path.join(enzyme_name, 'retained_part_protein_matrix.csv'))
    elif radius_to_exclude is not None:
        list_aa_to_discard = get_residues_near_cofactor(combined_pdb_file, reactive_residue_name, radius=radius_to_exclude)
        n = n[~n['residue ID'].isin(list_aa_to_discard)]
        n.to_csv(os.path.join(enzyme_name, 'retained_part_protein_matrix.csv')) 
    n = n[['x','y','z', 'charge']]
    n.to_csv(os.path.join(enzyme_name, 'chrg.txt'), index=False, sep=' ', mode = 'w+', header = None)

    return n 


def quantify_electric_field(enzyme_name, reactive_residue_name, cofactor_index, pdb_cofactor, df_coord_fe, 
    nearby_atom_list, other_metal_atom_coord, fe_indices_to_retain=[2,6], mode='amber', atom_indices_to_discard=[]):
    """
    """
    coord_list, coord_list_filtered = get_coord_list(pdb_cofactor, enzyme_name, reactive_residue_name, cofactor_index, fe_indices_to_retain)
    cys_bonded_index = determine_which_mutated_aa_is_closest(coord_list, df_coord_fe)
    electric_field_mapping = compute_electric_fields_on_grid(coord_list_filtered, nearby_atom_list, enzyme_name, 
                                                             df_coord_fe, other_metal_atom_coord, cys_bonded_index, mode, atom_indices_to_discard)
    
    return electric_field_mapping, cys_bonded_index


def get_coord_list(pdb_cofactor, path, reactive_residue_name, cofactor_index, atom_indices_to_retain):
    """
    Extracts and filters atomic coordinates of Fe atoms from a given DataFrame, based on chain ID, 
    reactive residue type, and specified atom indices to retain.

    Parameters:
    ----------
    pdb_cofactor : pd.DataFrame
        A DataFrame containing information about the structure, including atom types, coordinates, 
        chain IDs, residue types, and sites.
    path : str
        Path to save the CSV file containing Fe atom positions.
    reactive_residue_name : str
        Name of the reactive residue (e.g., the residue associated with Fe atoms).
    cofactor_index : int
        Index of the chain ID corresponding to the reactive residue.
    atom_indices_to_retain : list
        A list of atom indices (extracted from the `site` column) to retain in the filtered coordinates.

    Returns:
    -------
    tuple
        - coord_list : list of lists
            A list of all Fe atom coordinates as [x, y, z].
        - coord_list_filtered : list of lists
            A list of Fe atom coordinates that match the specified indices in `atom_indices_to_retain`.    
    """
    chain_id = pdb_cofactor['chain_id'].unique()[cofactor_index] # you can adjust this
    coord_list, coord_list_filtered = [], []
    df_fe = pdb_cofactor[pdb_cofactor['Atom'] == 'FE']
    df_fe = df_fe[df_fe['chain_id'] == f'{chain_id}']
    df_fe = df_fe[df_fe['AA'] == f'{reactive_residue_name}']

    df_fe.to_csv(os.path.join(path, 'Fe_positions.csv'))
    for _, row in df_fe.iterrows():
        coord = list([row['x'], row['y'], row['z']])
        coord_list.append(coord)

    df_fe['index'] = df_fe['type'].apply(lambda x: int(x.strip('FE')))
    df_fe = df_fe[df_fe['index'].isin(atom_indices_to_retain)]

    for _, row in df_fe.iterrows():
        coord = list([row['x'], row['y'], row['z']])
        coord_list_filtered.append(coord)

    return coord_list, coord_list_filtered


def determine_which_mutated_aa_is_closest(fe_pos, aa_ligated_fe_pos):
    """
    Determines which mutated amino acid is closest to the given iron position.

    Parameters:
    fe_pos (numpy.ndarray): The coordinates of the iron position (shape: (3,)).
    aa_ligated_fe_pos (pandas.DataFrame or numpy.ndarray): The coordinates of the ligated iron positions (shape: (n, 3)).

    Returns:
    int: The index of the closest ligated iron position or None if no match is found.
    """
    aa_coor = aa_ligated_fe_pos.to_numpy()
    # Perform element-wise comparison to find matching points
    matching_indices = np.where((aa_coor[:, None] == fe_pos).all(axis=2))
    # If there are matching indices, return the index of the first match
    if matching_indices[0].size > 0:
        index = matching_indices[0][0]
        return index
    else:
        return None


def compute_electric_fields_on_grid(coord_list, nearby_atom_list, path, df_coord_fe, other_metal_atom_coord, 
    index=0, mode='amber', atom_indices_to_discard=[]):
    """    
    """
    if mode == 'amber':
        q = initialize_charge_distribution_amber(f'{path}/chrg')
    elif mode == 'amoeba':
        q = initialize_charge_distribution_amoeba(f'{path}/{path}_amber_opt.pdb', atom_indices_to_discard)

    electric_field_mapping = pd.DataFrame(columns=['idx', 'x', 'y', 'z', 'efx', 'efy', 'efz', 'ef_tot', 'oef_x', 'oef_y', 'oef_z', 'oef'])
    grid_points, fe_pos = [], []
    fe_pos_dict = {}
    
    for i, coord in enumerate(coord_list):
        spherical_grid = compute_spherical_grid(coord)
        grid_points.append(spherical_grid)
        fe_pos.append(np.array([[coord[0], coord[1], coord[2]] for _ in range(len(spherical_grid))]))
        fe_pos_dict[str(fe_pos[-1][-1])] = i

    grid_points = np.concatenate(grid_points)
    fe_pos = np.concatenate(fe_pos)
    nearby_atom_pos = np.array(nearby_atom_list)
    grid_points, fe_pos = filter_points(grid_points, fe_pos, nearby_atom_pos, df_coord_fe, other_metal_atom_coord, index)
    grid_points, fe_pos = add_interpolation_between_centers(grid_points, fe_pos, coord_list)

    for j, point in enumerate(grid_points):
        [efx, efy, efz, ef_tot, oef_x, oef_y, oef_z, oef] = calc_ef(q, point[0], point[1], point[2], fe_pos[j][0], fe_pos[j][1], fe_pos[j][2])
        value = pd.Series([f'{fe_pos_dict[str(fe_pos[j])]}_{j}', point[0], point[1], point[2], 
                efx, efy, efz, ef_tot, oef_x, oef_y, oef_z, oef], index=electric_field_mapping.columns)
        electric_field_mapping = electric_field_mapping.append(value, ignore_index=True)
    electric_field_mapping.to_csv(os.path.join(path, 'electric_field_mapping.csv'))

    return electric_field_mapping


def add_interpolation_between_centers(grid_points, fe_pos, coord_list, n_points=6):
    """
    Adds interpolated points between two centers in a 3D space and updates the grid and Fe positions.

    Parameters:
    ----------
    grid_points : array-like
        A collection of existing 3D grid points.
    fe_pos : array-like
        A collection of 3D positions associated with Fe atoms or other reference points.
    coord_list : list of lists or arrays
        A list containing two 3D coordinates (e.g., [[x1, y1, z1], [x2, y2, z2]]) 
        between which the interpolation is performed.
    n_points : int, optional
        The number of interpolated points to generate between the two centers. Default is 6.

    Returns:
    -------
    tuple of np.ndarray
        - Updated grid points as a NumPy array, including the interpolated points.
        - Updated Fe positions as a NumPy array, with new points associated with the second coordinate in `coord_list`.
    """
    x_range = np.linspace(coord_list[0][0], coord_list[1][0], n_points+2)
    y_range = np.linspace(coord_list[0][1], coord_list[1][1], n_points+2)
    z_range = np.linspace(coord_list[0][2], coord_list[1][2], n_points+2) 

    grid_points = list(grid_points)
    fe_pos = list(fe_pos)

    for x,y,z in zip(x_range[1:-1], y_range[1:-1], z_range[1:-1]):
        grid_points.append(np.array([x,y,z]))
        fe_pos.append(np.array([coord_list[1][0], coord_list[1][1], coord_list[1][2]]))

    return np.array(grid_points), np.array(fe_pos)


def initialize_charge_distribution_amber(name):
    q = QuantificationTxt(name, 0, 0, 0, 0, 0, 0, 1, 1, 1)
    return q

def initialize_charge_distribution_amoeba(pdb_file, multipoles_to_exclude=[]):
    q = QuantificationAmoeba(pdb_file, 0, 0, 0, 0, 0, 0, 1, 1, 1, multipoles_to_exclude= multipoles_to_exclude)
    return q


def calc_ef(q, x, y, z, c_x, c_y, c_z):
    """
    Calculates electric field components and total electric field at a point using a given quadrupole object.
    """
    q.point_x = x
    q.point_y = y
    q.point_z = z
    q.v1_x = x
    q.v1_y = y
    q.v1_z = z
    q.v2_x = c_x
    q.v2_y = c_y
    q.v2_z = c_z

    efx, efy, efz, ef_tot, oef = q.execute()
    oef_x, oef_y, oef_z = decompose_oef_vector(q.v1_x, q.v1_y, q.v1_z, q.v2_x, q.v2_y, q.v2_z, oef)

    return efx * field_conversion_factor, efy * field_conversion_factor, efz * field_conversion_factor, ef_tot * field_conversion_factor, \
         oef_x * field_conversion_factor, oef_y * field_conversion_factor, oef_z * field_conversion_factor, oef * field_conversion_factor


def decompose_oef_vector(x1, y1, z1, x2, y2, z2, oef):
    """
    Decomposes a vector into its x, y, and z components.

    Parameters:
        x1, y1, z1 (float): Coordinates of the first point.
        x2, y2, z2 (float): Coordinates of the second point.
        oef (float): Magnitude of the vector.

    Returns:
        tuple: The x, y, and z components of the vector.
    """
    # Calculate the direction vector
    direction_vector = np.array([x2 - x1, y2 - y1, z2 - z1])
    
    # Normalize the direction vector
    norm = np.linalg.norm(direction_vector)
    if norm == 0:
        raise ValueError("The two points are identical; cannot define a direction.")
    unit_vector = direction_vector / norm
    
    # Scale by the magnitude of the vector (oef)
    components = unit_vector * oef
    
    return tuple(components)


def compute_spherical_grid(center):
    """    
    Computes spherical grids around a specified center point with varying radii and number of points.

    Parameters:
    center (tuple or list): Coordinates (x, y, z) of the center point for the spherical grids.

    Returns:
    numpy.ndarray: Concatenated spherical grids containing coordinates (x, y, z) in a 3D space.

    Workflow:
    1. Iterates through predefined radii and corresponding number of points:
       - For each radius, computes a spherical grid using polar coordinates (phi, theta).
       - Generates mesh grids for phi and theta covering spherical surfaces.
       - Converts polar coordinates to Cartesian coordinates (x, y, z).
       - Flattens and stacks the coordinates to form a spherical grid.
    2. Concatenates all generated spherical grids into a single numpy array.
    """
    spherical_grids = []
    
    for radius, num_points in [(1.5, 400), (2.0, 400), (2.5, 400), (3.0, 400)]: 
        phi = np.linspace(0, np.pi, int(np.sqrt(num_points)))
        theta = np.linspace(0, 2 * np.pi, int(np.sqrt(num_points)))
        phi, theta = np.meshgrid(phi, theta)

        x = center[0] + radius * np.sin(phi) * np.cos(theta)
        y = center[1] + radius * np.sin(phi) * np.sin(theta)
        z = center[2] + radius * np.cos(phi)

        spherical_grid = np.column_stack((x.flatten(), y.flatten(), z.flatten()))
        spherical_grids.append(spherical_grid)

    return np.concatenate(spherical_grids, axis=0)


def mean_point_3d(coordinates):
    n = len(coordinates)
    sum_x = sum(x for x, _, _ in coordinates)
    sum_y = sum(y for _, y, _ in coordinates)
    sum_z = sum(z for _, _, z in coordinates)
    return np.array([sum_x / n, sum_y / n, sum_z / n])


def filter_points(points, fe_pos, nearby_atom_pos, df_coor_fe, other_metal_atom_coord, index=0, distance_threshold=0.15, min_atom_distance=1.4):
    """    
    Filters grid points based on proximity criteria to avoid overlap with other points and atoms.

    Parameters:
    points (numpy.ndarray): Array of grid points in Cartesian coordinates (x, y, z).
    fe_pos (numpy.ndarray): Coordinates of iron (Fe) atoms associated with each grid point.
    nearby_atom_pos (numpy.ndarray): Coordinates of nearby atoms to be considered in the filtering.
    df_coor_fe (pandas.DataFrame): DataFrame containing coordinates of Fe atoms bonded to CYS residues.
    other_metal_atom_coord (numpy.ndarray): Coordinates of other metal atoms for proximity comparison.
    index (int): Index of the Dataframe containing the Fe atoms bonded to CYS residues to retain. Defaults to 0.
    distance_threshold (float, optional): Threshold distance for determining proximity between grid points. Defaults to 0.3 Å.
    min_atom_distance (float, optional): Minimum distance allowed between grid points and nearby atoms. Defaults to 1.4 Å.

    Returns:
    tuple: A tuple containing:
        - numpy.ndarray: Filtered grid points after all proximity filters have been applied.
        - numpy.ndarray: Coordinates of Fe atoms associated with each filtered grid point.

    Workflow:
    1. Filters out grid points that are too close to each other using a distance threshold.
    2. Filters out grid points that are too close to any nearby atoms based on a minimum atom distance.
    3. Filters out grid points that are located within the interior of the cofactor.
    4. Filters out grid points that are too close to either the CYS/HIS/TYR-bonded Fe atoms or other metal atoms.
    """
    # first filter out the grid points that are close to each other
    distances = cdist(points, points)
    close_points = (distances <= distance_threshold)

    # Exclude self-comparisons
    np.fill_diagonal(close_points, False)

    for i in range(len(close_points)):
        close_points[i, i:] = False
    
    # Filter out points
    filtered_indices = np.logical_not(np.any(close_points, axis=1))
    filtered_points = points[filtered_indices]
    fe_pos = fe_pos[filtered_indices]

    print("Original number of points:", len(points))
    print("Filtered number of points:", len(filtered_points))
    
    # Second, filter out the grid points that are close to any of the nearby atoms
    distances = cdist(filtered_points, nearby_atom_pos)
    invalid_points = (distances <= min_atom_distance)

    # Exclude self-comparisons
    np.fill_diagonal(invalid_points, False)

    # Filter out points
    filtered_indices = np.logical_not(np.any(invalid_points, axis=1))
    filtered_points = filtered_points[filtered_indices]
    fe_pos = fe_pos[filtered_indices]

    # Third, filter out points at the interior of the cofactor
    cofactor_center = mean_point_3d(fe_pos)
    filtered_indices = []
    for i, point in enumerate(filtered_points):
        if np.linalg.norm(point - fe_pos[i]) + 0.5 < np.linalg.norm(point - cofactor_center):
            filtered_indices.append(i)

    filtered_points = filtered_points[filtered_indices]
    fe_pos = fe_pos[filtered_indices]

    # Finally, remove points too close to CYS/HIS/THYR-bonded metal or other metal
    filtered_indices = []
    fe_center_coord = np.array([df_coor_fe['x'][index], df_coor_fe['y'][index], df_coor_fe['z'][index]])

    for i, point in enumerate(filtered_points):
        if np.linalg.norm(point - fe_pos[i]) + 0.5 < np.linalg.norm(point - fe_center_coord) and \
        np.linalg.norm(point - fe_pos[i]) + 0.5 < np.linalg.norm(point - other_metal_atom_coord):
            filtered_indices.append(i)
    
    filtered_points = filtered_points[filtered_indices]
    fe_pos = fe_pos[filtered_indices] 

    print("Valid number of points:", len(filtered_points))

    return filtered_points, fe_pos
