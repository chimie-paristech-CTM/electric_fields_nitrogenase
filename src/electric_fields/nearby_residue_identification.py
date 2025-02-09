from Bio import PDB

def get_residues_near_cofactor(pdb_file, cofactor_name, radius=5.0):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    # Extract all atoms of the given cofactor
    cofactor_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                if res_name == cofactor_name:
                    cofactor_atoms.extend(residue.get_atoms())
                    found_cofactor = True 
            if found_cofactor:
                break
        if found_cofactor:
            break

    # Do the same for HCA
    hca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                if res_name == 'HCA':
                    hca_atoms.extend(residue.get_atoms())
                    found_hca = True 
            if found_hca:
                break
        if found_hca:
            break
    #******

    if not cofactor_atoms:
        raise ValueError(f"Cofactor {cofactor_name} not found in the PDB file.")
    
    # Find residues within the given radius
    ns = PDB.NeighborSearch(list(structure.get_atoms()))
    nearby_residues = set()
    
    for atom in cofactor_atoms:
        for neighbor in ns.search(atom.coord, radius):
            residue = neighbor.get_parent()
            if residue.id[0] == " ":  # Ignore heteroatoms except standard residues
                nearby_residues.add(residue)

    for atom in hca_atoms:
        for neighbor in ns.search(atom.coord, radius):
            residue = neighbor.get_parent()
            if residue.id[0] == " ":  # Ignore heteroatoms except standard residues
                nearby_residues.add(residue) 
    
    # Output residue indices
    residue_indices = [(res.parent.id, res.id[1]) for res in nearby_residues]
    print(len(residue_indices), [index[1] for index in residue_indices])
    return [index[1] for index in residue_indices]
