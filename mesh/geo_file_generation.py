#%% imports
import geopandas as gpd
from pathlib import Path
import pandas as pd
import numpy as np
import copy

#%%
parent_directory = Path(r"C:\Users\03125327\github\fc_hydro_kalimantan_2022")
fn_pointers = parent_directory.joinpath(r'file_pointers.xlsx')
filenames_df = pd.read_excel(fn_pointers, header=2, dtype=str, engine='openpyxl')

# channel_network_nodes_fn = Path(filenames_df[filenames_df.Content ==
#                 'channel_network_nodes'].Path.values[0])
channel_network_lines_fn = Path(filenames_df[filenames_df.Content ==
                'channel_network_lines'].Path.values[0])

study_area_boundary_points_fn = Path(filenames_df[filenames_df.Content ==
                'study_area_boundary_points'].Path.values[0])
output_geo_fn = Path(filenames_df[filenames_df.Content ==
                'gmsh_geo'].Path.values[0])
output_geo_fn = parent_directory.joinpath(output_geo_fn)
#%% Read channel network points and study area points
def extract_coords_from_geodataframe(gdf:gpd.geodataframe)->list:
    coords = []
    for _, row in gdf.iterrows():
        coords.append((row.geometry.x, row.geometry.y))
    
    return coords

def get_points_from_linestring(linestring):
    xy = linestring.xy
    xcoords = xy[0].tolist()
    ycoords = xy[1].tolist()
    return [(xy[0][i], xy[1][i]) for i in range(len(xy[0]))]

def get_list_of_canal_lines(gdf: gpd.geodataframe)->list:
    canal_lines = []
    for _, row in channel_network_lines_gdf.iterrows():
        canal_lines.append(get_points_from_linestring(row.geometry)) 
    return canal_lines
    

# channel_network_nodes_gdf = gpd.read_file(channel_network_nodes_fn)
channel_network_lines_gdf = gpd.read_file(channel_network_lines_fn)
study_area_boundary_points_gdf = gpd.read_file(study_area_boundary_points_fn)

study_area_boundary_points_coords = extract_coords_from_geodataframe(study_area_boundary_points_gdf)
channel_network_lines = get_list_of_canal_lines(channel_network_lines_gdf)


# %% Generate a .geo file for gmsh

MESH_DIMENSIONS_GLOBAL = 500 # meters (?)
MESH_DIMENSIONS_AT_CANAL_POINTS = 50 # meters (?)
with open(output_geo_fn, "w") as f:
    # Define the mesh sizes
    f.write(f"lc = {MESH_DIMENSIONS_GLOBAL};\n") 
    f.write(f"lc2 = {MESH_DIMENSIONS_AT_CANAL_POINTS};\n") 
    f.write("IP = newp;\n") 
    f.write("IL = newl;\n") 
    f.write("IS = news;\n") 
    f.write("ILL = newll;\n") 

    # Write surface points
    point_number = 0 # keep track of point number
    for sa_coord in study_area_boundary_points_coords:
        f.write(f"Point(IP+{point_number}) = {{" + f"{sa_coord[0]}, {sa_coord[1]}, 0, lc}};\n")
        point_number += 1

    # Make the surface points into a line loop and a surface
    points_string = str()
    for i in range(0, point_number):
        points_string = points_string + f"IP+{i}, "
    points_string = points_string + "IP+0"
    f.write("Line(IL+0) = {" + points_string + "}; \n")
    f.write("Line Loop(ILL+0) = {IL+0}; \n")
    f.write("Plane Surface(IS) = {ILL:ILL+0}; \n")
    f.write('Physical Surface("Domain") = {IS}; \n')

    first_canal_point_number = copy.copy(point_number) # store this value for the next loop

    # Add points and lines in the canals with their mesh size
    line_count = 1
    latest_canal_point_number = first_canal_point_number
    for canal_line in canal_lines:
        for i, canal_coord in enumerate(canal_line):
            f.write(f"Point(IP+{latest_canal_point_number + i}) = {{" + f"{canal_coord[0]}, {canal_coord[1]}, 0, lc2}};\n")
            point_number += 1
    
        # Make a line with the canal points
        points_string = str()
        for i in range(latest_canal_point_number, point_number):
            points_string = points_string + f"IP+{i}, "
        f.write(f"Line(IL+{line_count}) =" + " {" + points_string[:-2] + "}; \n")

        line_count += 1
        latest_canal_point_number = point_number

    # Add canal points to surface
    for n in range(first_canal_point_number, first_canal_point_number + len(channel_network_nodes_coords) ):
        f.write("Point{IP+" + str(n) + "} In Surface{IS}; \n")
    
# %%
# %%
