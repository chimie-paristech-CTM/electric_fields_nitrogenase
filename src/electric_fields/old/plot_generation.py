import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.io as pio
import pandas as pd
import os


def generate_plots(enzyme_name, electric_field_mapping, df_coor_fe, other_metal_atom_coord, index=0):
    """
    Generates vector, 3D, and quiver plots based on electric field mapping and Fe atom coordinates.

    Parameters:
    enzyme_name (str): Name of the enzyme or system for generating plots.
    electric_field_mapping (pd.DataFrame): DataFrame containing electric field mapping data,
                                           including 'x', 'y', 'z', 'oef', 'efx', 'efy', 'efz' columns.
    df_coor_fe (pd.DataFrame): DataFrame containing Fe atom coordinates with 'x', 'y', 'z' columns.
    other_metal_atom_coord (List): List of coordinates of the other metal atom which is part of the reactive cofactor.
    index (int): Index of the Dataframe containing the Fe atoms bonded to CYS residues to retain. Defaults to 0.

    Workflow:
    1. Generates a vector plot (3D plot with vector arrows representing electric field magnitude and direction).
    2. Generates a 3D plot with regular colored dots representing electric field magnitude.
    3. Generates a quiver plot (3D plot with quiver arrows representing electric field vectors).

    Outputs:
    - Creates '3d_plot_vect.html', '3d_plot.html', and '3d_plot_quivers.html' in the `enzyme_name` directory.    
    """
    # vector plot first
    output_file = f'{enzyme_name}/3d_plot_vect.html'
    coor_cys_fe = [df_coor_fe['x'][index], df_coor_fe['y'][index], df_coor_fe['z'][index]]
    #electric_field_mapping = pd.read_csv(os.path.join(enzyme_name, 'electric_field_mapping.csv'))
    fe_atoms = pd.read_csv(os.path.join(enzyme_name, 'Fe_positions.csv'))
    atoms = pd.read_csv(os.path.join(enzyme_name, f'{enzyme_name}_protein_matrix.csv'))

    generate_plot_vector(
        output_file, electric_field_mapping['x'], electric_field_mapping['y'], 
        electric_field_mapping['z'], electric_field_mapping['oef'], fe_atoms['x'], fe_atoms['y'], fe_atoms['z'],
        electric_field_mapping['efx'], electric_field_mapping['efy'], 
        electric_field_mapping['efz'], atoms['x'], atoms['y'], atoms['z'], 
        [coor_cys_fe[0]], [coor_cys_fe[1]], [coor_cys_fe[2]], [other_metal_atom_coord[0]], [other_metal_atom_coord[1]],
        [other_metal_atom_coord[2]]
    )

    # then 3D plot with regular colored dots
    output_file = f'{enzyme_name}/3d_plot.html'
    electric_field_mapping = pd.read_csv(os.path.join(enzyme_name, 'electric_field_mapping.csv'))
    fe_atoms = pd.read_csv(os.path.join(enzyme_name, 'Fe_positions.csv'))

    generate_plot(
        output_file, electric_field_mapping['x'], electric_field_mapping['y'], 
        electric_field_mapping['z'], electric_field_mapping['oef'], fe_atoms['x'], fe_atoms['y'], fe_atoms['z'], 
        [coor_cys_fe[0]], [coor_cys_fe[1]], [coor_cys_fe[2]], [other_metal_atom_coord[0]], [other_metal_atom_coord[1]],
        [other_metal_atom_coord[2]]
    )

    # finally quiver plot
    output_file = f'{enzyme_name}/3d_plot_quivers.html'

    generate_quiver_plot(output_file, electric_field_mapping, fe_atoms['x'], fe_atoms['y'], fe_atoms['z'],
                        [coor_cys_fe[0]], [coor_cys_fe[1]], [coor_cys_fe[2]], [other_metal_atom_coord[0]], 
                        [other_metal_atom_coord[1]], [other_metal_atom_coord[2]])  


# function to generate 3D plot with vectors and magnitude
def generate_plot_vector(output_file, x, y, z, values, x2, y2, z2, efx, efy, efz, xx, yy, zz, x_cys, y_cys, z_cys, x_other_metals, y_other_metals, z_other_metals):
    # Create a scatter plot
    trace1 = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=4,
            color=values,
            colorscale='spectral',
            cmin=min(values),
            cmax=0,
            opacity=0.8
        ),
        text=[f'oef: {oef:2f}' for oef in values],
        hoverinfo='text'
    )

    # Create a layout
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Z'),
        ),
        paper_bgcolor='rgba(0,0,0,0)',  # Transparent background
        plot_bgcolor='rgba(0,0,0,0)',
        title='Evenly Distributed Points on a Sphere'
    )

    # Create scatter plots for the second set of points
    trace2 = go.Scatter3d(
        x=x2,
        y=y2,
        z=z2,
        mode='markers',
        marker=dict(
            size=12,
            color='violet',
            opacity=0.8)     
    )
    # Plot vectors 
    trace3 = go.Cone(
        x=x,
        y=y,
        z=z,
        u=efx,
        v=efy,
        w=efz,
        sizeref=0.02,
        colorscale="Blues",         
        anchor="tail",
        name='Vectors',
        opacity=0.3 ,
        colorbar=dict(thickness=20, ticklen=4),
        sizemode="absolute",      
        #text=[f'oef: {oef:2f}' for oef in values],
        hoverinfo='text'
    )
    
    # Plot the Fe/Cu atom bonded to Cys/Ala  
    trace5 = go.Scatter3d(
        x=x_cys,
        y=y_cys,
        z=z_cys,
        mode='markers',
        marker=dict(
            size=12,
            color='blue',
            opacity=0.3)
    )

    # Plot the other metal atoms 
    trace6 = go.Scatter3d(
        x=x_other_metals,
        y=y_other_metals,
        z=z_other_metals,
        mode='markers',
        marker=dict(
            size=12,
            color='blue',
            opacity=0.3)
    )
    
    # Create a figure. Add trace4 to see the cofactor in the en molecule
    fig = go.Figure(data=[trace1, trace2, trace3, trace5, trace6], layout=layout)
    
    # Save the figure as an interactive HTML file
    pio.write_html(fig, f'{output_file}', auto_open=False)  


def generate_quiver_plot(output_file, electric_field_mapping, x_fe, y_fe, z_fe, x_cys, y_cys, z_cys, x_other_metals, y_other_metals, z_other_metals):
    # Add the 3D quiver plot
    trace1 = go.Cone(
        x=electric_field_mapping['x'],
        y=electric_field_mapping['y'],
        z=electric_field_mapping['z'],
        u=electric_field_mapping['efx'],
        v=electric_field_mapping['efy'],
        w=electric_field_mapping['efz'],
        sizemode="scaled",
        sizeref=2,
        anchor="tip"
    )

    # Create scatter plots for the metal centers
    trace2 = go.Scatter3d(
        x=x_fe,
        y=y_fe,
        z=z_fe,
        mode='markers',
        marker=dict(
            size=12,
            color='violet',
            opacity=0.8
        ),
        name='Second Set'
    )

    # Plot the Fe/Cu atom bonded to Cys/Ala  
    trace3 = go.Scatter3d(
        x=x_cys,
        y=y_cys,
        z=z_cys,
        mode='markers',
        marker=dict(
            size=12,
            color='blue',
            opacity=0.3)
    )

    # Plot the other metal atoms 
    trace4 = go.Scatter3d(
        x=x_other_metals,
        y=y_other_metals,
        z=z_other_metals,
        mode='markers',
        marker=dict(
            size=12,
            color='blue',
            opacity=0.3)
    )

    # Create the figure
    fig = go.Figure(data=[trace1, trace2, trace3, trace4])  

    # Update the layout
    fig.update_layout(
        scene=dict(
        xaxis=dict(title='X'),
        yaxis=dict(title='Y'),
        zaxis=dict(title='Z')
        ),
        paper_bgcolor='rgba(0,0,0,0)',  # Transparent background
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent plot background
        title='3D Electric Field Mapping'
    )  

    # Save the figure as an interactive HTML file
    pio.write_html(fig, f'{output_file}', auto_open=False) 


def generate_plot(output_file, x, y, z, values, x_fe, y_fe, z_fe, x_cys, y_cys, z_cys, x_other_metals, y_other_metals, z_other_metals):
    # Create a scatter plot
    trace1 = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=4,
            color=values,
            colorscale='spectral',
            cmin=min(values),
            cmax=0,
            opacity=0.8
        ),
        text=[f'oef: {oef:2f}' for oef in values],
        hoverinfo='text'

    )

    # Create a layout
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Z'),
        ),
        paper_bgcolor='rgba(0,0,0,0)',  # Transparent background
        plot_bgcolor='rgba(0,0,0,0)',
        title='Evenly Distributed Points on a Sphere'
    )

    # Create scatter plots for the metal centers
    trace2 = go.Scatter3d(
        x=x_fe,
        y=y_fe,
        z=z_fe,
        mode='markers',
        marker=dict(
            size=12,
            color='violet',
            opacity=0.8
        ),
        name='Second Set'
    )

    # Plot the Fe/Cu atom bonded to Cys/Ala  
    trace3 = go.Scatter3d(
        x=x_cys,
        y=y_cys,
        z=z_cys,
        mode='markers',
        marker=dict(
            size=12,
            color='blue',
            opacity=0.3)
    )

    # Plot the other metal atoms 
    trace4 = go.Scatter3d(
        x=x_other_metals,
        y=y_other_metals,
        z=z_other_metals,
        mode='markers',
        marker=dict(
            size=12,
            color='blue',
            opacity=0.3)
    )

    # Create a figure
    fig = go.Figure(data=[trace1, trace2, trace3, trace4], layout=layout)

    # Save the figure as an interactive HTML file
    pio.write_html(fig, output_file, auto_open=False) 
