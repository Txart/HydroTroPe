
#%% imports
import geopandas
import networkx as nx
import pickle
import pandas as pd
from pathlib import Path
import argparse

import cwl_preprocess_data

# %% Parse command-line arguments
parser = argparse.ArgumentParser(description='Perform Step 3 in generation of the graph from the canal shapefiles.')

parser.add_argument('--no-blocks', dest='block_opt', action='store_false', 
                    help='Do not consider blocks in the graph creation step. \n This option exists for when there exist no block shapefile. It is recommended to always create graphs with blocks because later you may choose to run the model without them.' )

args = parser.parse_args()
parser.ser_defaults(block_opt=True)
BLOCKS_OPT = args.block_opt

#%% Read file
print('Reading file...')
filepaths = pd.read_excel('file_pointers.xlsx')
fn_dem = Path(filepaths[filepaths.variable_name == 'fn_dem'].path.iloc[0])
fn_regular_split = Path(filepaths[filepaths.variable_name == 'fn_regular_split'].path.iloc[0])
fn_blocks = Path(filepaths[filepaths.variable_name == 'fn_blocks'].path.iloc[0])
# output file
fn_graph_pickle = Path(filepaths[filepaths.variable_name == 'fn_graph_pickle'].path.iloc[0])

# %% STEP 3. Create the graph and save it as a pickle file
graph = cwl_preprocess_data.compute_channel_network_graph(gpkg_fn=fn_regular_split, fn_dtm=fn_dem)

#%% Add blocks to graph
if BLOCKS_OPT == True:
    # coords of blocks need to exactly coincide with channel network nodes' coords ( did that in QGIS)
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
