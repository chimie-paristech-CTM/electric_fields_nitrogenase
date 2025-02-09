import os
import subprocess
import parmed as pmd


def get_charge_parameters_amber_opt(enzyme_name, mol2_file='enzyme_with_charges.mol2'):
    """
    """
    home_dir = os.getcwd()
    os.chdir(enzyme_name)
    write_leap_input(f'{enzyme_name}_amber_opt')
    run_tleap()
    prmtop_file = f"{enzyme_name}_amber_opt.prmtop"
    inpcrd_file = f"{enzyme_name}_amber_opt.crd" 

    # Load the structure
    structure = pmd.load_file(prmtop_file, inpcrd_file)

    # Create a Mol2 file with charges included
    structure.save(mol2_file, overwrite=True)
    print(f"Mol2 file with charges saved as {mol2_file}")
    os.chdir(home_dir)

    return mol2_file


def combine_residue_and_heteroatom_positions(initial_pdb, amber_opt_pdb, enzyme_name):
    """
    Combines residue and heteroatom positions from two PDB files into a single output file.

    This function merges the positions of residues from the initial PDB file with the optimized heteroatom 
    positions from the AMBER-optimized PDB file. It appends the heteroatom information (excluding water) 
    to the optimized PDB file and saves the combined structure to a new PDB file.

    Args:
        initial_pdb (str): Path to the initial PDB file containing the starting residue and heteroatom information.
        amber_opt_pdb (str): Path to the AMBER-optimized PDB file containing the updated residue positions.
        enzyme_name (str): Name of the enzyme, used to save the combined structure into a PDB file.

    Returns:
        str: Path to the output PDB file containing the combined structure.
    """
    output_file = f'{enzyme_name}/{enzyme_name}_amber_opt_combined.pdb' 
    assign_chains_on_reset(amber_opt_pdb, output_file)
    with open(initial_pdb, 'r') as infile, open(output_file, 'a') as outfile:
        for line in infile:
            if line.startswith("HETATM") and line.split()[2] != 'HOH':
                outfile.write(line)
    
    return output_file


def optimize_amber(enzyme_name):
    """    
    Performs an AMBER-based optimization of the enzyme structure.

    This function automates the process of optimizing an enzyme structure using the AMBER force field. 
    It includes steps for extracting atomic information, preparing the structure, minimizing the energy, 
    and extracting the optimized PDB file.

    Args:
        enzyme_name (str): Name of the enzyme, used for locating the enzyme's directory and saving the output files.

    Returns:
        str: Path to the optimized PDB file produced by the AMBER optimization.
    """
    home_dir = os.getcwd()
    os.chdir(enzyme_name)
    extract_aa_atoms(enzyme_name, f'{enzyme_name}_final.pdb')
    prepare_input_structure(enzyme_name)
    minimize_structure(enzyme_name)
    extract_optimized_pdb(enzyme_name)
    os.chdir(home_dir)

    return f'{enzyme_name}/{enzyme_name}_amber_opt.pdb'


def extract_aa_atoms(enzyme_name, input_file):
    """
    Extracts atomic information for amino acid atoms from a PDB file, excluding hydrogen and OXT atoms.

    This function processes an input PDB file and writes out all lines corresponding to amino acid atoms 
    (i.e., atoms that are not hydrogen and do not belong to the OXT residue) to an output PDB file.

    Args:
        enzyme_name (str): Name of the enzyme, used to construct the output file name.
        input_file (str): Path to the input PDB file containing the enzyme structure.

    Returns:
        None: Writes the extracted atomic information to a new PDB file named "{enzyme_name}_amber.pdb".
    """

    # Define the input and output file paths
    output_file = f"{enzyme_name}_amber.pdb"  # Replace with your desired output .pdb file name

    # Open the input file for reading and the output file for writing
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            # Check if the line starts with "ATOM"
            if line.startswith("ATOM") and line.split()[-1] != 'H' and 'OXT' not in line:
            # Write the line to the output file
                outfile.write(line)


def prepare_input_structure(enzyme_name):
    write_leap_input(f'{enzyme_name}_amber')
    run_tleap()


def minimize_structure(enzyme_name):
    write_sander_input()
    write_minimization_bash_script(f'{enzyme_name}_amber')
    run_minimization('run_minimization.sh')


def extract_optimized_pdb(enzyme_name):
    create_extract_pdb_script(enzyme_name)
    make_and_run_extraction_script()


def write_leap_input(input_name):
    """
    Generate a LEaP input file for a specific PDB input.

    Parameters:
        input_name (str): Name of the input PDB file (without extension).
    """

    # LEaP input file content template
    leap_content = f"""source leaprc.protein.ff14SB

    protein = loadPdb "{input_name}.pdb"
    saveAmberParm protein {input_name}.prmtop {input_name}.crd

    # Save a PDB file with hydrogens added (optional)
    savePdb protein {input_name}_with_hydrogens.pdb

    # Quit tleap
    quit
    """
    # Define the output file name
    output_filename = "prepare_leap.in"

    # Write the content to the LEaP input file
    with open(output_filename, "w") as file:
        file.write(leap_content)

    return output_filename


def run_tleap():
    """
    Run the 'tleap' command.
    """
    try:
        # Run the tleap command as a subprocess
        result = subprocess.run(
            ["tleap", "-f", "prepare_leap.in"],
            capture_output=True,
            text=True,
            check=True
        )
        # Print the output from tleap
        print("TLEaP Output:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        # Print error message if tleap fails
        print("TLEaP Error:")
        print(e.stderr)


def write_sander_input():
    """
    Generate a minimization input file for sander.

    Parameters:
        input_name (str): Name of the input PDB file (without extension).
    """

    sander_content = f"""# Minimization input file for sander                                           
    &cntrl                                                                         
      imin = 1,         
      maxcyc = 1000,      
      ncyc = 250,        
      ntb = 0,  
      ntr = 1,         
      igb = 0,           
      cut = 12.0,                  
      restraint_wt = 10000.0,  
      restraintmask = '!@H=' 
    /
    """

    # Write the content to the LEaP input file
    with open("minimization.in", "w") as file:
        file.write(sander_content)


def write_minimization_bash_script(file_name):
    bash_script_content = f"""#!/bin/bash

    # Define input and output files
    INPUT=minimization.in
    PRMTOP={file_name}.prmtop       # Topology file from tleap
    INPCRD={file_name}.crd       # Coordinate file from tleap
    OUTCRD=minimized.rst     # Output coordinates after minimization
    OUTLOG=minimization.log  # Log file for minimization

    # Run sander minimization
    sander -O -i $INPUT -ref $INPCRD -p $PRMTOP -c $INPCRD -r $OUTCRD -o $OUTLOG"""

    # Write the content to the LEaP input file
    with open("run_minimization.sh", "w") as file:
        file.write(bash_script_content)


def run_minimization(script_name="run_minimization.sh"):
    """
    Make a shell script executable and then run it.

    Parameters:
        script_name (str): Name of the shell script to execute (default is "run_minimization.sh").
    """
    try:
        # Step 1: Make the script executable
        chmod_command = ["chmod", "+x", script_name]
        subprocess.run(chmod_command, check=True)

        # Step 2: Execute the script
        run_command = ["./" + script_name]
        result = subprocess.run(run_command, capture_output=True, text=True, check=True)

        # Print the output from the script
        print("Script Output:")
        print(result.stdout)

    except subprocess.CalledProcessError as e:
        # Print error if the script fails
        print(f"Error while executing '{script_name}':")
        print(e.stderr)

    except FileNotFoundError:
        print(f"Error: '{script_name}' not found or not in the current directory.")


def create_extract_pdb_script(enzyme_name, output_file="extract_pdb.sh"):
    """
    Create a shell script to extract a PDB file from a topology and restart file using cpptraj.

    Parameters:
        output_file (str): The name of the shell script to create (default is "extract_pdb.sh").
    """
    script_content = f"""#!/bin/bash

# Define input and output files
TOPOLOGY="{enzyme_name}_amber.prmtop"    # Topology file generated by tleap
RESTART="minimized.rst"   # Restart file from the minimization
OUTPUT_PDB="{enzyme_name}_amber_opt.pdb" # Final PDB output file

# Check if required files exist
if [[ ! -f $TOPOLOGY ]]; then
    echo "Error: Topology file $TOPOLOGY not found!"
    exit 1
fi

if [[ ! -f $RESTART ]]; then
    echo "Error: Restart file $RESTART not found!"
    exit 1
fi

# Create cpptraj input file
CPPTRAJ_INPUT="extract_geometry.in"
cat << EOF > $CPPTRAJ_INPUT
parm $TOPOLOGY
trajin $RESTART
trajout $OUTPUT_PDB pdb
EOF

# Run cpptraj
cpptraj -i $CPPTRAJ_INPUT > cpptraj.log

# Check if output was created
if [[ -f $OUTPUT_PDB ]]; then
    echo "PDB file created successfully: $OUTPUT_PDB"
else
    echo "Error: Failed to create PDB file. Check cpptraj.log for details."
fi

# Clean up temporary input file
rm $CPPTRAJ_INPUT
    """
    # Write the content to the output file
    with open(output_file, "w") as file:
        file.write(script_content)


def make_and_run_extraction_script(script_name="extract_pdb.sh"):
    """
    Make the specified script executable and execute it.

    Parameters:
        script_name (str): The name of the script to execute (default is "extract_pdb.sh").
    """
    try:
        # Step 1: Ensure the script is executable
        os.chmod(script_name, 0o755)  # Add execute permissions
        print(f"'{script_name}' has been made executable.")

        # Step 2: Execute the script
        result = subprocess.run(
            ["./" + script_name],
            capture_output=True,
            text=True,
            check=True
        )

        # Print the output from the script
        print("Script Output:")
        print(result.stdout)

    except FileNotFoundError:
        print(f"Error: Script '{script_name}' not found.")
    except subprocess.CalledProcessError as e:
        print(f"Error while executing '{script_name}':")
        print(e.stderr)


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
            elif line.startswith("ATOM"):
                # Extract residue number
                res_num = int(line[22:26].strip())

                # Assign the current chain ID
                chain_id = chain_ids[current_chain_index]

                # Modify the chain ID (column 22) and write the updated line
                new_line = line[:21] + chain_id + line[22:]
                outfile.write(new_line)

                # Update the last residue number
                last_residue_number = res_num
