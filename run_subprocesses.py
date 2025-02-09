import subprocess
import shlex


def run_calculations(commands):
    """
    Run a list of shell commands in sequence as subprocesses.
    
    Parameters:
        commands (list): A list of shell command strings to execute.
    """
    for i, command in enumerate(commands, start=1):
        print(f"Running command {i}/{len(commands)}: {command}")
        try:
            # Split the command string into a list for subprocess
            process = subprocess.run(shlex.split(command), check=True, text=True, capture_output=True)
            print(f"Output of command {i}:")
            print(process.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error while running command {i}: {e}")
            print(f"Error output: {e.stderr}")
        except Exception as e:
            print(f"Unexpected error while running command {i}: {e}")


if __name__ == "__main__":
    # Define the list of commands to run
    commands = []
    enzyme_residue_pair = [('3u7q', 195)]
    #enzyme_residue_pair = [('1m1n', 195), ('1g20', 195), ('1h1l',194), ('3u7q',195), ('8p8g', 195), ('4wza', 195), ('4wzb', 195), 
    #('3k1a', 195), ('2afh', 195), ('8e3u', 195), ('7mci', 195), ('6vxt', 195), ('6o7r', 195), ('6o7p', 195), ('6o7l', 195), ('6cdk', 195)]
    #enzyme_residue_pair = [('5vq4', 195), ('5vq3', 195), ('5koh', 195), ('4xpi', 195), ('4wes', 195), ('4tku', 195), 
    #        ('4nd8', 195), ('3min', 195), ('2min', 195), ('1qh8', 195), ('1qgu', 195)]

    #enzyme_residue_pair = [('1gpu', 194), ('5qv3', 186), ('5koh', 211), ('4wes', 186), ('1qh8', 194), ('5vq3', 186)]
    #enzyme_residue_pair = [('8p8g', 195)]

    # Vanadium
    #enzyme_residue_pair = [('5n6y', 180), ('6fea', 0), ('7ady', 0)]
    #enzyme_residue_pair = [('6fea', 1000)] 

    #____________
    #enzyme_residue_pair = [('3u7q',195), ('8e3u', 195), ('7mci', 195), ('2afh', 195)]
    #enzyme_residue_pair = [('8p8g',195)]

    for enzyme_name, residue_to_discard in enzyme_residue_pair:
        if residue_to_discard != 0:
            for prot_state in ['HID', 'HIE', 'HIP']:
                commands.append(f"python run_scripts/analyze_electric_fields_amber.py --enzyme-name {enzyme_name} --prot-state-his {prot_state} --fe-indices-to-retain 2 6")
            commands.append(f"python run_scripts/analyze_electric_fields_amber.py --enzyme-name {enzyme_name} --prot-state-his HID --fe-indices-to-retain 2 6 --residue-to-discard {residue_to_discard}")
            #commands.append(f"python run_scripts/analyze_electric_fields_amber.py --enzyme-name {enzyme_name} --prot-state-his HID --fe-indices-to-retain 4 5")
            #commands.append(f"python run_scripts/analyze_electric_fields_amber.py --enzyme-name {enzyme_name} --prot-state-his HID --fe-indices-to-retain 3 7")
        else: # protonation state is not relevant
            commands.append(f"python run_scripts/analyze_electric_fields_amber.py --enzyme-name {enzyme_name} --prot-state-his HID --fe-indices-to-retain 2 6") 
    print(commands)
    run_calculations(commands)
