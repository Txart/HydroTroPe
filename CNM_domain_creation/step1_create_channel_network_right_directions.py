
#%% imports
import numpy as np
from shapely.geometry import LineString, MultiPoint
from shapely.ops import split, linemerge
import shapely
import pandas as pd
import geopandas
import rasterio
from pathlib import Path
import pickle

import cwl_preprocess_data


# %% STEP 1. Produce canal network of right global direction
# Take a canal shapefile where junctions and ends are nodes,
# and each branch is an edge. Compute  direction of water flow
# accoding to the dem at each node. Save that canal network with the right
# water flow direction to a file.


#def explode_line_into_segments(line):
#    return list(map(LineString, zip(line.coords[:-1], line.coords[1:])))

def multilinestring_to_single(gdf, use_linemerge=True):
    if use_linemerge:
        return gdf['geometry'].apply(lambda x: linemerge(x))
    else:
        return gdf['geometry'].apply(lambda x: x[0])

def reverse_line(line):
    def _reverse(x, y, z=None):
        if z:
            return x[::-1], y[::-1], z[::-1]
        return x[::-1], y[::-1]

    return shapely.ops.transform(_reverse, line)

def reverse_if_wrong_direction(gdf_row):
    if gdf_row['beginnings_dem'] < gdf_row['ends_dem']:
        return reverse_line(gdf_row['geometry'])
    else:
        return gdf_row['geometry']

def sample_raster_from_geodataframe_of_points(raster_filename:str, point_coords:np.ndarray) ->np.ndarray:
    raster_src = rasterio.open(raster_filename)
    return np.array([raster_value[0] for raster_value in raster_src.sample(point_coords)])

def get_all_line_endpoints_from_geopackage_file(lines_gpkg):
    """Take a read geopackage file containing only lines 
    and spit out a list of the first and last points (endpoints) of that line.

    Args:
        lines_gpkg (geopandas geodataframe): holds the output of  geopandas.readfile("*.gpkg")

    Returns:
        (list): list of ends of lines in the geopackage file
    """
    all_line_endpoints = []
    for line in lines_gpkg.geometry:
        try:
            n1 = line.coords[0]
            n2 = line.coords[-1]
            all_line_endpoints.append(n1)
            all_line_endpoints.append(n2)
        except:
            raise Warning(" There was some problem with a line when trying to get its coords")
            continue

    return all_line_endpoints

def get_unique_nodes_in_line_gpkg(lines_gpkg):
    # nodes are the unique values of all_line_endpoints
    return list(set(get_all_line_endpoints_from_geopackage_file(lines_gpkg)))

def get_lists_of_node_coords_and_node_heights(raster_src, lines_gpkg):
    # get all unique line edpoints
    nodes_coords = get_unique_nodes_in_line_gpkg(lines_gpkg)
    
    nodes_heights = []
    nodes_coords_copy = nodes_coords.copy()
    for i, dem_height in enumerate(raster_src.sample(nodes_coords)):
        if np.isnan(dem_height[0]) or dem_height[0] < 0:  # Coordinates are outside the DEM
            # Remove from list of nodes
            nodes_coords_copy.remove(nodes_coords[i])
        else:
            nodes_heights.append(dem_height[0])
    nodes_coords = nodes_coords_copy.copy()
    del nodes_coords_copy

    return nodes_coords, nodes_heights

#%%
# Read files
print('Reading files...')
filepaths = pd.read_excel('file_pointers.xlsx')
fn_dem = Path(filepaths[filepaths.variable_name == 'fn_dem'].path.iloc[0])
fn_splitted_at_junctions = Path(filepaths[filepaths.variable_name == 'fn_splitted_at_junctions'].path.iloc[0])
# output file
fn_right_direction_lines = Path(filepaths[filepaths.variable_name == 'fn_right_direction_lines'].path.iloc[0])

print('Computing line directions')
raster_src = rasterio.open(fn_dem)
lines_gdf = geopandas.read_file(filename=fn_splitted_at_junctions)
lines_gdf = lines_gdf[~lines_gdf.is_empty] # Remove lines without any geometry
if type(lines_gdf.geometry[0]) == shapely.geometry.MultiLineString:
    lines_gdf['geometry'] = multilinestring_to_single(lines_gdf, use_linemerge=True)
nodes_coords, nodes_heights = get_lists_of_node_coords_and_node_heights(raster_src, lines_gdf)

lines_gdf['beginnings'] = lines_gdf['geometry'].apply(lambda x: x.coords[0])
lines_gdf['ends'] = lines_gdf['geometry'].apply(lambda x: x.coords[-1])
lines_gdf['beginnings_dem'] = sample_raster_from_geodataframe_of_points(raster_filename=fn_dem, point_coords=lines_gdf['beginnings'].to_numpy())
lines_gdf['ends_dem'] = sample_raster_from_geodataframe_of_points(raster_filename=fn_dem, point_coords=lines_gdf['ends'].to_numpy())

right_direction_lines = lines_gdf.apply(func=reverse_if_wrong_direction, axis=1)
right_direction_lines_gdf = geopandas.GeoDataFrame(right_direction_lines, columns=["geometry"])
right_direction_lines_gdf.to_file(fn_right_direction_lines, engine='fiona', crs=lines_gdf.crs)
print('STEP 1 completed! \n')
print(f'canal network with right directions saved to: {fn_right_direction_lines}')
# %%
