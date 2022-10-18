
#%% imports
import geopandas
import networkx as nx
import pickle
import pandas as pd
from pathlib import Path

import cwl_preprocess_data

#%% Read file
print('Reading file...')
filepaths = pd.read_excel('file_pointers.xlsx')
fn_dem = Path(filepaths[filepaths.variable_name == 'fn_dem'].path.iloc[0])
fn_regular_split = Path(filepaths[filepaths.variable_name == 'fn_regular_split'].path.iloc[0])
# output file
fn_graph_pickle = Path(filepaths[filepaths.variable_name == 'fn_graph_pickle'].path.iloc[0])

# %% STEP 3. Create the graph and save it as a pickle file
graph = cwl_preprocess_data.compute_channel_network_graph(gpkg_fn=fn_regular_split, fn_dtm=fn_dem)

#%% Add blocks to graph
# coords of blocks need to exactly coincide with channel network nodes' coords ( did that in QGIS)
fn_blocks = r"C:\Users\03125327\github\paper2\data\new_area\dams.gpkg"
gdf = geopandas.read_file(fn_blocks)
node_xs = nx.get_node_attributes(graph, 'x')
node_ys = nx.get_node_attributes(graph, 'y')
notfound = 0
node_coords = {}
for i in range(len(graph.nodes)):
    node_coords[i] = (node_xs[i], node_ys[i])

nx.set_node_attributes(graph, values=0, name='is_block') # initialize all nodes of channel network as having no block
for i, values in gdf.iterrows():
    coords_block = (values.geometry.x, values.geometry.y)
    try:
        network_node_with_block = list(node_coords.keys())[list(node_coords.values()).index(coords_block)]
        graph.nodes[network_node_with_block]['is_block'] = 1
    except:
        notfound += 1
    
print(f"A total of {notfound} blocks were not found in the search " )

# %% Dump graph into pickle
pickle.dump(graph, open(fn_graph_pickle, "wb"))
# %%
