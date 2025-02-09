import pandas as pd
import plotly.graph_objects as go
from ase.io import read
from ase.data import colors, covalent_radii
import numpy as np
import os


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
    field_x = np.array([point['x'] for point in field_data])
    field_y = np.array([point['y'] for point in field_data])
    field_z = np.array([point['z'] for point in field_data])
    oef = np.array([point['oef'] for point in field_data])

    field_x_pos = np.array([point['x'] for point in field_data if point['oef'] > 0])
    field_y_pos = np.array([point['y'] for point in field_data if point['oef'] > 0])
    field_z_pos = np.array([point['z'] for point in field_data if point['oef'] > 0])
    oef_pos = np.array([point['oef'] for point in field_data if point['oef'] > 0])

    print(oef, len(oef))
    print(oef_pos, len(oef_pos))

    # Create a scatter plot for the electric field
    field_trace = go.Scatter3d(
        x=field_x,  # Select points with negative values
        y=field_y,
        z=field_z,
        mode='markers',
        marker=dict(
            size=5,
            color=oef,
            colorscale='spectral',
            colorbar=dict(title='OEF magnitude (a.u.)', x=0.8, len=0.5),
            cmin=-0.02, #min(oef),
            cmax=0,
            opacity=1
        ),
        text=oef,
        hoverinfo='text'
    )

    positive_trace = go.Scatter3d(
        x=field_x_pos,  # Select points with positive values
        y=field_y_pos,
        z=field_z_pos,
        mode='markers',
        marker=dict(
            size=5,
            color=oef_pos,
            colorscale='greys', 
            cmin=0,
            cmax=max(oef),
            opacity=1  # Ensure full visibility
        ),
        hoverinfo='text'  # Disable hoverinfo for the overlay if unnecessary
    )

    # Combine all traces
    data = bond_traces + [field_trace] + [positive_trace] + [atom_trace]
    
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


def nitrogenase_plot(enzyme_name, xyz_file, field_file):
    molecule = parse_xyz_with_ase(xyz_file)
    field_data = parse_electric_field(field_file)
    visualize_molecule_with_field_plotly(molecule, field_data, os.path.join(enzyme_name, f'{enzyme_name}_plot.html'))


if __name__ == '__main__':
    nitrogenase_plot('3u7q/fragment.xyz', '3u7q/electric_field_mapping.csv')
