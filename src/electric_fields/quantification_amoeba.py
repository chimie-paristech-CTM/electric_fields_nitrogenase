import numpy as np
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *

Bohr_to_nm = 0.0529177
Bohr_to_ang = 0.529177
au_to_VA = 51.4220675112


class QuantificationAmoeba():
    """
    A class that represents a quantification calculation starting from an AMOEBA force field calculation

    Attributes:
    ----------
    name : string
        The name of the input and output files
    point_x : float
        The x-coordinate of the point at which the electric field will be quantified
    point_y : float
        The y-coordinate of the point at which the electric field will be quantified
    point_z : float
        The z-coordinate of the point at which the electric field will be quantified
    v1_x : float
        The x-coordinate of the first point making up the direction vector
    v1_y : float
        The y-coordinate of the first point making up the direction vector
    v1_z : float
        The z-coordinate of the first point making up the direction vector
    v2_x : float
        The x-coordinate of the second point making up the direction vector
    v2_y : float
        The y-coordinate of the second point making up the direction vector
    v2_z : float
        The z-coordinate of the second point making up the direction vector

    Methods:
    -------
    execute : executes the workflow associated to the quantification calculation
    """
    def __init__(self, pdb_file, point_x, point_y, point_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, multipoles_to_exclude):
        self.point_x = point_x / Bohr_to_ang 
        self.point_y = point_y / Bohr_to_ang
        self.point_z = point_z / Bohr_to_ang
        self.v1_x = v1_x / Bohr_to_ang
        self.v1_y = v1_y / Bohr_to_ang
        self.v1_z = v1_z / Bohr_to_ang
        self.v2_x = v2_x / Bohr_to_ang
        self.v2_y = v2_y / Bohr_to_ang
        self.v2_z = v2_z / Bohr_to_ang
        self.amoeba_parameters, self.positions = self.determine_amoeba_multipoles(pdb_file)
        self.multipoles_to_exclude = multipoles_to_exclude

    def get_direction_vector(self):
        return np.array([self.v2_x - self.v1_x, self.v2_y - self.v1_y, self.v2_z - self.v1_z])

    def determine_amoeba_multipoles(self, pdb_file):
        pdb, system, forcefield = load_system(pdb_file)
        multipole_force, induced_dipoles = compute_multipoles(system, pdb)
    
        # Create simulation context
        integrator = openmm.VerletIntegrator(0.001)
        platform = Platform.getPlatformByName('CPU')
        context = Context(system, integrator, platform)
        context.setPositions(pdb.positions)

        amoeba_parameters, positions = get_amoeba_parameters_in_reference_frame(multipole_force, induced_dipoles, context)

        return amoeba_parameters, positions

    def rescale_positions(self):
        self.point_x = self.point_x / Bohr_to_ang 
        self.point_y = self.point_y / Bohr_to_ang
        self.point_z = self.point_z / Bohr_to_ang
        self.v1_x = self.v1_x / Bohr_to_ang
        self.v1_y = self.v1_y / Bohr_to_ang
        self.v1_z = self.v1_z / Bohr_to_ang
        self.v2_x = self.v2_x / Bohr_to_ang
        self.v2_y = self.v2_y / Bohr_to_ang
        self.v2_z = self.v2_z / Bohr_to_ang 

    def execute(self):
        self.rescale_positions()
        direction_vector = self.get_direction_vector()
        efx, efy, efz = self.compute_electric_field()
        ef_tot = np.sqrt(efx ** 2 + efy ** 2 + efz ** 2)
        oef = self.compute_oef(direction_vector, np.array([efx, efy, efz]))

        return efx, efy, efz, ef_tot, oef
    
    def compute_oef(self, direction_vector, field_vector):
        direction_norm = np.linalg.norm(direction_vector)
        normalized_direction = direction_vector / direction_norm

        return np.dot(normalized_direction, field_vector)

    def compute_electric_field(self):
        """Computes the electric field at a given point due to AMOEBA charges and multipoles."""
    
        eval_point = np.array([self.point_x, self.point_y, self.point_z])

        field = np.zeros(3)

        for i in range(len(self.amoeba_parameters)):
            if i in self.multipoles_to_exclude:
                continue
            pos = self.positions[i]
            q = self.amoeba_parameters[i]['charge']
            dipole = self.amoeba_parameters[i]['dipole']
            induced_dipole = self.amoeba_parameters[i]['induced_dipole']
            Q = self.amoeba_parameters[i]['quadrupole']

            r_vec = eval_point - pos
            r = np.linalg.norm(r_vec)
            r_hat = r_vec / r

            # Electric field from charge (Coulomb's law)
            field += np.array(q * r_hat / r**2)

            # Dipole contribution
            p_dot_r = np.dot(dipole, r_hat)
            field += np.array((3 * p_dot_r * r_hat - dipole) / r**3)

            # Induced dipole contribution
            p_dot_r_ind = np.dot(induced_dipole, r_hat)
            field += np.array((3 * p_dot_r_ind * r_hat - induced_dipole) / r**3)
        
            # Quadrupole contribution (second-order term)
            # Compute scalar term: (r_hat^T Q r_hat)
            scalar_term = np.dot(r_hat, np.dot(Q, r_hat))

            # Compute vector term: Q * r_hat
            vector_term = np.array(np.dot(Q, r_hat))

            # Quadrupole electric field contribution
            field_quad = (5 * scalar_term * r_hat - 2 * vector_term) / r**4

            # Add quadrupole contribution to total field
            field += field_quad
        
        return field


def load_system(pdb_filename):
    """Loads a PDB file and sets up an AMOEBA force field system."""
    pdb = PDBFile(pdb_filename)
    forcefield = ForceField('amoeba2013.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None)
    return pdb, system, forcefield


def compute_multipoles(system, pdb):
    """Extracts atomic multipoles from the AMOEBA force field."""
    for force in system.getForces():
        if force.__class__.__name__ == 'AmoebaMultipoleForce':
            multipole_force = force
            break
    # Create an integrator and simulation
    integrator = LangevinIntegrator(300 * kelvin, 1/picosecond, 1 * femtoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    print('AMOEBA minimization... (may take some time)')
    # Minimize energy to converge induced dipoles
    simulation.minimizeEnergy(maxIterations=500)

    # Get induced dipoles (electronic polarization at each site)
    induced_dipoles = multipole_force.getInducedDipoles(simulation.context)

    return force, induced_dipoles


def get_amoeba_parameters_in_reference_frame(multipole_force, induced_dipoles, context):
    state = context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True)
    
    # unit conversion
    scaled_positions = []
    for position in positions:
        scaled_positions.append(position / (nanometer * Bohr_to_nm))

    amoeba_parameters = []

    for i in range(multipole_force.getNumMultipoles()):
        (charge, dipole, quadrupole, axis_type, z_atom, x_atom, y_atom, thole, damping, polarity) = multipole_force.getMultipoleParameters(i)
        induced_dipole = np.array(induced_dipoles[i])
        
        # unit conversions -- all distances in Bohr, and all charges in e
        charge = charge / elementary_charge
        dipole = dipole / (nanometer * elementary_charge * Bohr_to_nm)
        quadrupole = quadrupole / (nanometer**2 * elementary_charge * Bohr_to_nm ** 2)

        quadrupole = np.array(quadrupole).reshape(3,3)
        
        induced_dipole = np.array(induced_dipoles[i])
        rotation_matrix = get_multipole_rotation_matrix(scaled_positions, i, x_atom, y_atom, z_atom)

        amoeba_parameters.append({'charge': charge, 'dipole': rotate_dipole(dipole, rotation_matrix), 
            'induced_dipole': rotate_dipole(induced_dipole, rotation_matrix), 'quadrupole': transform_quadrupole(quadrupole, rotation_matrix)
            })

    return amoeba_parameters, np.array(scaled_positions)


def normalize(v):
    """ Normalize a vector to unit length. """
    return v / np.linalg.norm(v)


def get_multipole_rotation_matrix(positions, i, x_atom, y_atom, z_atom):
    central_pos = positions[i]

    if z_atom != -1:
        z_pos = positions[z_atom]
        e1 = normalize(z_pos - central_pos)  # Define local z-axis
    else:
        raise ValueError("Multipole frame requires at least a Z-axis reference.")

    if x_atom != -1:
        x_pos = positions[x_atom]
        e_temp = x_pos - central_pos
        e3 = normalize(np.cross(e1, e_temp))  # Local y-axis (ensuring perpendicularity)
        e2 = np.cross(e3, e1)  # Local x-axis (ensuring right-handed system)
    else:
        # If no X-axis reference, choose an arbitrary perpendicular direction
        e_temp = np.array([1.0, 0.0, 0.0]) if abs(e1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        e3 = normalize(np.cross(e1, e_temp))
        e2 = np.cross(e3, e1)

    # Construct the rotation matrix (columns are local frame vectors in global coordinates)
    R = np.column_stack((e1, e2, e3))

    return R


def rotate_dipole(dipole_local, R):
    """
    Transforms a dipole moment from the local frame to the global frame.

    Parameters:
    - dipole_local: (3,) NumPy array representing the dipole moment in the local frame
    - R: (3,3) Rotation matrix from local to global frame

    Returns:
    - dipole_global: (3,) NumPy array representing the dipole in the global frame
    """
    return np.dot(R, dipole_local)


def transform_quadrupole(quadrupole_local, R):
    """
    Transforms a quadrupole moment tensor from the local frame to the global frame.

    Parameters:
    - quadrupole_local: (3,3) NumPy array representing the quadrupole moment in the local frame
    - R: (3,3) Rotation matrix from local to global frame

    Returns:
    - quadrupole_global: (3,3) NumPy array representing the quadrupole in the global frame
    """
    return R @ quadrupole_local @ R.T # Use @ for matrix multiplication 
    
if __name__ == "__main__":
    pdb_filename = "3u7q_amber_opt.pdb"  # Replace with actual PDB filename
    q = QuantificationAmoeba(pdb_filename, 11.764049203807025,-6.902,55.27172586255095,11.764049203807025,-6.902,55.27172586255095,11.277,-6.902,53.853)
    efx, efy, efz, ef_tot, oef = q.execute()
    print(efx * au_to_VA, efy * au_to_VA, efz * au_to_VA, ef_tot * au_to_VA, oef * au_to_VA)
