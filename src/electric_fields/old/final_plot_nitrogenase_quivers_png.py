import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io import read
from ase.data import colors, covalent_radii
import os


def parse_xyz_with_ase(file_path):
    molecule = read(file_path)
    return molecule


def parse_electric_field(file_path):
    df = pd.read_csv(file_path)
    field_data = {
        'x': df['x'].values,
        'y': df['y'].values,
        'z': df['z'].values,
        'oef_x': df['oef_x'].values,
        'oef_y': df['oef_y'].values,
        'oef_z': df['oef_z'].values,
        'oef_magnitude': np.sqrt(df['oef_x']**2 + df['oef_y']**2 + df['oef_z']**2),
    }
    return field_data


def visualize_molecule_with_field_matplotlib(molecule, field_data, output_png):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot atoms (double the size, no transparency)
    atom_x = [atom.x for atom in molecule]
    atom_y = [atom.y for atom in molecule]
    atom_z = [atom.z for atom in molecule]
    atom_colors = [colors.jmol_colors[atom.number] for atom in molecule]

    # Plot bonds
    for i, atom1 in enumerate(molecule):
        for j, atom2 in enumerate(molecule):
            if i < j:
                distance = atom1.position - atom2.position
                if np.linalg.norm(distance) < (covalent_radii[atom1.number] + covalent_radii[atom2.number]) * 1.2:
                    ax.plot(
                        [atom1.x, atom2.x],
                        [atom1.y, atom2.y],
                        [atom1.z, atom2.z],
                        color='k', lw=1
                    )

    ax.scatter(atom_x, atom_y, atom_z, s=200, c=atom_colors, label="Atoms", alpha=1)  # Size doubled to 200

    # Plot electric field cones
    field_x = field_data['x']
    field_y = field_data['y']
    field_z = field_data['z']
    oef_x = field_data['oef_x']
    oef_y = field_data['oef_y']
    oef_z = field_data['oef_z']
    oef_magnitude = field_data['oef_magnitude']

    # Normalize vectors for cones
    #norm = np.sqrt(oef_x**2 + oef_y**2 + oef_z**2)
    u = oef_x * 100 #/ norm
    v = oef_y * 100 #/ norm
    w = oef_z * 100 #/ norm

    print(u,v,w)

    # Scale vectors for better visualization
    cone_scale = 0.5
    u = u * cone_scale
    v = v * cone_scale
    w = w * cone_scale

    print(oef_magnitude[2] / np.max(oef_magnitude))
    # Plot cones with color based on magnitude
    for i in range(len(field_x)):
        ax.quiver(
            field_x[i], field_y[i], field_z[i],
            u[i], v[i], w[i],
            color=plt.cm.Spectral_r(oef_magnitude[i] / 0.02),
            linewidth= oef_magnitude[i] / 0.03,
            length=cone_scale,
            normalize=True
        )

    # Customize plot appearance
    ax.set_axis_off()  # Turn off grid and axes

    # Add a color bar for cone magnitudes
    sm = plt.cm.ScalarMappable(cmap='Spectral_r', norm=plt.Normalize(vmin=0, vmax=oef_magnitude[i] / np.max(oef_magnitude)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.1, shrink=0.7)
    cbar.set_label('OEEF (a.u.)', fontsize=12)

    # Save plot
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close(fig)


def nitrogenase_plot_with_quivers_png(enzyme_name, xyz_file, field_file, output_file):
    molecule = parse_xyz_with_ase(xyz_file)
    field_data = parse_electric_field(field_file)
    visualize_molecule_with_field_matplotlib(molecule, field_data, output_file)


if __name__ == '__main__':
    nitrogenase_plot_with_cones_png(
        enzyme_name='1m1n',
        xyz_file='1m1n/fragment.xyz',
        field_file='1m1n/electric_field_mapping.csv',
        output_file='1m1n_molecule_with_field_cones.png'
    )
