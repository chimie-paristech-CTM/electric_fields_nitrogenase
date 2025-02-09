import pandas as pd
import math
import re
from Bio.PDB import *
import numpy as np
import os
import argparse

import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from electric_fields.fragment_extraction_amber import extract_fragment_around_active_site
from electric_fields.enzyme_preparation_and_mutation_amber import prepare_and_mutate_enzyme
from electric_fields.electric_field_quantification import quantify_electric_field, extract_charges_protein_matrix_amber, get_discarded_atom_indices
from electric_fields.final_plot_nitrogenase_quivers_html import nitrogenase_plot_with_quivers_html
from electric_fields.amber_optimization import optimize_amber, combine_residue_and_heteroatom_positions, get_charge_parameters_amber_opt


def get_args():
    """
    Parse command line arguments.

    Returns:
    - argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--enzyme-name', action='store', type=str, default='1m1n')
    parser.add_argument('--element-active-site', action='store', type=str, default='Fe')
    parser.add_argument('--cofactor-index', action='store', type=int, default=0)
    parser.add_argument('--prot-state-his', action='store', type=str, default='HIP')
    parser.add_argument('--fe-indices-to-retain', action='store', type=int, nargs='+', default=[2, 6])
    parser.add_argument('--residue-to-discard', action='store', type=int, default=None)
    parser.add_argument('--radius-to-discard', action='store', type=float, default=None)
    parser.add_argument('--charge-mode', action='store', type=str, default='amber')

    return parser.parse_args()


def count_missing_residues(pdb_file, max_residue):
    """
    Counts the number of missing residues in chain A up to a given residue number.
    
    Parameters:
        pdb_file (str): Path to the PDB file.
        max_residue (int): The residue number up to which missing residues are counted.
    
    Returns:
        int: Number of missing residues up to max_residue in chain A.
    """
    observed_residues = set()
    missing_residues = set(range(1, max_residue + 1))  # Assume a continuous sequence
    
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[21] == 'A':  # Ensure we are in chain A
                try:
                    residue_num = int(line[22:26].strip())
                    observed_residues.add(residue_num)
                except ValueError:
                    continue

    missing_residues -= observed_residues  

    return len(missing_residues)


def main(enzyme_name, element_active_site, cofactor_index, prot_state_his, fe_indices_to_retain, residue_to_discard, radius_to_exclude, charge_mode):
    """
    Main function to prepare, optimize, and analyze an enzyme, including calculating the electric field and visualizing the results.

    This function performs a series of steps to prepare the enzyme structure, mutate it as needed, optimize it using Amber,
    compute electric field data, and visualize the results. It also extracts relevant atomic coordinates and performs 
    post-processing on the enzyme structure.

    Args:
        enzyme_name (str): Name of the enzyme to be analyzed, used for file paths and outputs.
        element_active_site (str): The element (e.g., 'Fe', 'Cu') at the active site of the enzyme.
        cofactor_index (int): The index of the cofactor involved in the active site.
        prot_state_his (str): Protonation state of the critical HIS195 residue.
        fe_indices_to_retain (list): Indices of the Fe atoms around which to quantify the electric field.
        residue_to_discard (int): index of residue to discard (typically the adjacent HIS...)
        charge_mode (str): whether to use 'amber' charges, or 'amoeba' multipoles

    Returns:
        float: The minimum value of the electric field (oef) from the electric field mapping, extracted from the output CSV file.
    """
    # prepare enzyme
    reactive_residue_name, df_coord_fe, pdb_cofactor = prepare_and_mutate_enzyme(enzyme_name, element_active_site, prot_state_his)
    amber_opt_file = optimize_amber(enzyme_name)

    #amber_opt_file = os.path.join(enzyme_name, f'{enzyme_name}_amber_opt.pdb')

    # extract charges + compute electric field
    combined_pdb_file = combine_residue_and_heteroatom_positions(os.path.join(enzyme_name, f'{enzyme_name}.pdb'), amber_opt_file, enzyme_name)
    nearby_atom_list, other_metal_atom_coord = extract_fragment_around_active_site(combined_pdb_file, enzyme_name, reactive_residue_name, cofactor_index, 
                                        radius=8, displacement=False, xyz_file='fragment.xyz')

    if charge_mode == 'amber':
        mol2_file = get_charge_parameters_amber_opt(enzyme_name, f'{enzyme_name}_with_charges_amber_{prot_state_his}.mol2')
        if residue_to_discard != None:
            residue_to_discard_renumbered = residue_to_discard - count_missing_residues(os.path.join(enzyme_name, f'{enzyme_name}.pdb'), residue_to_discard)
            extract_charges_protein_matrix_amber(enzyme_name, mol2_file, reactive_residue_name, combined_pdb_file, 
            list_aa_to_discard=[residue_to_discard_renumbered])
        elif radius_to_exclude != None: 
            extract_charges_protein_matrix_amber(enzyme_name, mol2_file, reactive_residue_name, combined_pdb_file, 
            radius_to_exclude=radius_to_exclude)
        else:
            extract_charges_protein_matrix_amber(enzyme_name, mol2_file, reactive_residue_name, combined_pdb_file)

        _, _ = quantify_electric_field(enzyme_name, reactive_residue_name, cofactor_index, 
                            pdb_cofactor, df_coord_fe, nearby_atom_list, other_metal_atom_coord, fe_indices_to_retain, mode='amber')

    elif charge_mode == 'amoeba':
        if residue_to_discard != None:
            residue_to_discard_renumbered = residue_to_discard - count_missing_residues(os.path.join(enzyme_name, f'{enzyme_name}.pdb'), residue_to_discard)
            atom_indices_to_discard = get_discarded_atom_indices(amber_opt_file, combined_pdb_file, reactive_residue_name, list_residue_indices_to_discard=[residue_to_discard_renumbered])
        elif radius_to_exclude != None:
            atom_indices_to_discard = get_discarded_atom_indices(amber_opt_file, combined_pdb_file, reactive_residue_name, radius_to_exclude=radius_to_exclude)
        else:
            atom_indices_to_discard = []

        _, _ = quantify_electric_field(enzyme_name, reactive_residue_name, cofactor_index, 
                            pdb_cofactor, df_coord_fe, nearby_atom_list, other_metal_atom_coord, fe_indices_to_retain, mode='amoeba', atom_indices_to_discard=atom_indices_to_discard)

    # post-processing
    nearby_atom_list, other_metal_atom_coord = extract_fragment_around_active_site(combined_pdb_file, enzyme_name, reactive_residue_name, cofactor_index, 
                                        radius=2.4, displacement=False, xyz_file='fragment_small.xyz')
    
    base_name_plot = f'{enzyme_name}_quiver_plot_{prot_state_his}_{fe_indices_to_retain[0]}_{fe_indices_to_retain[1]}'

    if residue_to_discard != None:
        base_name_plot += f'_no_his{residue_to_discard}'
    elif radius_to_exclude != None:
        base_name_plot += f'_radius{radius_to_exclude}'

    if charge_mode == 'amber':
        base_name_plot += '_amber'
    elif charge_mode == 'amoeba':
        base_name_plot += '_amoeba'

    quiver_plot_html = f'{base_name_plot}.html'

    nitrogenase_plot_with_quivers_html(enzyme_name, os.path.join(enzyme_name, 'fragment_small.xyz'), os.path.join(enzyme_name, 'electric_field_mapping.csv'), 
                                os.path.join(enzyme_name, quiver_plot_html))
    df = pd.read_csv(os.path.join(enzyme_name, 'electric_field_mapping.csv'))

    print(df[['oef']].min())

    return df[['oef']].min()


if __name__ == '__main__':
    args = get_args()
    oef_max = main(args.enzyme_name, args.element_active_site, args.cofactor_index, args.prot_state_his, args.fe_indices_to_retain, 
        args.residue_to_discard, args.radius_to_discard, args.charge_mode)
