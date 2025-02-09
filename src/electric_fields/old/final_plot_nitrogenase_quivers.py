import pandas as pd
import plotly.graph_objects as go
from ase.io import read
from ase.data import colors, covalent_radii
import numpy as np


def parse_xyz_with_ase(file_path):
    molecule = read(file_path)
    return molecule


def parse_electric_field(file_path):
    df = pd.read_csv(file_path)
    field_data = []
    for _, row in df.iterrows():
        field_point = {
            'x': row['x'],
            'y': row['y'],
            'z': row['z'],
            'efx': row['efx'],
            'efy': row['efy'],
            'efz': row['efz'],
            'oef': row['oef']
        }
        field_data.append(field_point)
    return field_data


def visualize_molecule_with_field_plotly(molecule, field_data, output_html):
    atom_radius = 0.25
    # Extract atom positions and elements
    atom_x = [atom.x for atom in molecule]
    atom_y = [atom.y for atom in molecule]
    atom_z = [atom.z for atom in molecule]
    atom_colors = [colors.jmol_colors[atom.number] for atom in molecule]
    
    # Create a trace for bonds
    bond_traces = []
    for i, atom1 in enumerate(molecule):
        for j, atom2 in enumerate(molecule):
            if i < j:
                distance = atom1.position - atom2.position
                if np.linalg.norm(distance) < (covalent_radii[atom1.number] + covalent_radii[atom2.number]) * 1.2:
                    prop_factor = np.sqrt(atom_radius**2 / ((atom1.x - atom2.x)**2 + (atom1.y - atom2.y)**2 + (atom1.z - atom2.z)**2))

                    if atom1.x < atom2.x:
                        x_start = atom1.x + prop_factor * abs(atom1.x - atom2.x)
                        x_end = atom2.x - prop_factor * abs(atom1.x - atom2.x)
                    else:
                        x_start = atom1.x - prop_factor * abs(atom1.x - atom2.x)
                        x_end = atom2.x + prop_factor * abs(atom1.x - atom2.x)
                    
                    if atom1.y < atom2.y:
                        y_start = atom1.y + prop_factor * abs(atom1.y - atom2.y)
                        y_end = atom2.y - prop_factor * abs(atom1.y - atom2.y)
                    else:
                        y_start = atom1.y - prop_factor * abs(atom1.y - atom2.y)
                        y_end = atom2.y + prop_factor * abs(atom1.y - atom2.y)

                    if atom1.z < atom2.z:
                        z_start = atom1.z + prop_factor * abs(atom1.z - atom2.z)
                        z_end = atom2.z - prop_factor * abs(atom1.z - atom2.z)
                    else:
                        z_start = atom1.z - prop_factor * abs(atom1.z - atom2.z)
                        z_end = atom2.z + prop_factor * abs(atom1.z - atom2.z)                    

                    bond_traces.append(go.Scatter3d(
                        x=[x_start, x_end],
                        y=[y_start, y_end],
                        z=[z_start, z_end],
                        mode='lines',
                        line=dict(color='rgb(50,50,50,0.5)', width=3)
                    ))

    # Create a scatter plot for atoms
    atom_trace = go.Scatter3d(
        x=atom_x,
        y=atom_y,
        z=atom_z,
        mode='markers',
        marker=dict(
            size=15,
            color=atom_colors,
            line=dict(width=3, color='rgb(50,50,50)'), opacity=1,
        ),
        text=[molecule[i].symbol for i in range(len(molecule))],
        hoverinfo='text'
    )
    
    # Extract field data
    field_x = [point['x'] for point in field_data]
    field_y = [point['y'] for point in field_data]
    field_z = [point['z'] for point in field_data]
    efx = [-point['efx'] for point in field_data]
    efy = [-point['efy'] for point in field_data]
    efz = [-point['efz'] for point in field_data]
    oef = [-point['oef'] for point in field_data]

    # Create a quiver plot for the electric field
    field_trace = go.Cone(
        x=field_x,
        y=field_y,
        z=field_z,
        u=efx,
        v=efy,
        w=efz,
        sizemode="scaled",
        sizeref=1,
        anchor="tail",
        colorscale='Spectral',  # Use spectral color scale
        cmin=0,          # Set the minimum for color scaling
        cmax=max(oef),          # Set the maximum for color scaling
        reversescale=True,
        colorbar=dict(title='EF (a.u.)', titleside='top', x=0.8, len=0.5),
    )
    
    # Combine all traces
    data = bond_traces + [field_trace] + [atom_trace]
    
    # Define layout
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Z'),
            bgcolor='white'
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        paper_bgcolor='white',  # Set the overall background color to white
        plot_bgcolor='white',  # Set the plot area background color to white
        showlegend=False
    )
 
    # Create figure
    fig = go.Figure(data=data, layout=layout)

    fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)
    
    # Save the figure as an HTML file
    fig.write_html(output_html, config=dict(displayModeBar=False))


def nitrogenase_plot_with_quivers(xyz_file, field_file):
    molecule = parse_xyz_with_ase(xyz_file)
    field_data = parse_electric_field(field_file)
    visualize_molecule_with_field_plotly(molecule, field_data, 'test.html')


if __name__ == '__main__':
    nitrogenase_plot_with_quivers('1m1n/fragment.xyz', '1m1n/electric_field_mapping.csv')
