import numpy as np
from tqdm import tqdm
import multiprocessing as mp
import networkx as nx

from classes.channel_network import ChannelNetwork

import math_preissmann, math_diff_wave
import cwl_utilities
import utilities

#%%
class CWLHydroParameters:
    def __init__(self, dt, dx, params_channel, 
                 downstream_diri_BC) -> None:
        self.g = 9.8  # m/s², Acceleration of gravity
        self.dt = dt  # s, timestep
        self.dx = dx  # m, delta x from finite dif.

        # Parameters for the function describing the n_manning channel network friction coefficient.
        self.porous_threshold_below_dem = float(params_channel['porous_threshold']) # position where n_manning = max, in metres from surface
        self.n_threshold = float(params_channel['n_threshold']) # Maximum value of n_manning, capped below the porous threshold.
        self.n1 = float(params_channel['n1']) # n_manning is prop to exp(-n1*y^n2), where y is elevation from ref. datum
        self.n2 = float(params_channel['n2']) 

        # weight parameter in Preissmann scheme
        self.a = 0.6  

        # Max number of iterations for...
        self.max_niter_newton = int(params_channel['max_niter_newton']) # the exact Newton method
        self.max_niter_inexact = int(params_channel['max_niter_newton_inexact']) # the inexact Newton method
        
        # Tolerance to test solution convergence
        self.rel_tol = 1e-7
        self.abs_tol = 1e-7

        # Weights to scale the solution of the Newton method to update the solution variables
        self.weight_A = 5e-2
        self.weight_Q = 5e-2

        # outer loop. If =1, it simulates dt seconds into the future. If =2, 2*dt, etc. So it is usually =1.
        self.ntimesteps = 1

        # Number of cpus to use by the cwl algorithm. If >1, multiprocess.
        self.ncpus = 1 

        # boolean for Downstream BC for the channel reaches. If True, ChannelNetwork.y_BC_below_DEM dirichlet values area applied at downstream nodes. If False, no flux Neumann are imposed.
        self.downstream_diri_BC = downstream_diri_BC

        pass


class AbstractChannelHydrology:
    """Abstract class. Always create one of its subclasses.
    """

    def __init__(self, cn: ChannelNetwork, cwl_params: CWLHydroParameters):
        self.cwl_params = cwl_params
        self.cn = cn
        
        # A list of channel_network classes of each component, to loop the hydro computation through them
        self.component_channel_networks = self._list_channel_networks_of_components()
        
        # List of nodes that are in channel network but don't appear in any of the components
        self.nodes_not_in_components = self._compute_nodes_not_in_components()
        
        # These are implemented in child classes
        self.cwl_model_algorithm = None  
        self.y_solution = None # Implemented after model run
        self.Q_solution = None

        pass

    def _list_channel_networks_of_components(self):
        # Return a list of channel network classes, one for each separated component
        # each one with the same params as the global one.
        return [ChannelNetwork(
                graph=g_com,
                params_channel=self.cn.params_channel,
                work_without_blocks=self.cn.work_without_blocks,
                is_components=True # crucial bit to compute the diffusive wave approx.
                ) for g_com in self.cn.component_graphs]
        
    def _compute_nodes_not_in_components(self):
        all_nodes_in_components = []
        for component_cn in self.component_channel_networks:
            all_nodes_in_components += list(component_cn.graph.nodes)
        nodes_not_in_components = list(set(range(0, self.cn.n_nodes)) - set(all_nodes_in_components))
        
        return nodes_not_in_components

    def _store_results(self, results):
        # Implemented by child class
        return None

    def _run_singleprocess(self):
        results = []
        for com_channel_network in tqdm(self.component_channel_networks):
            results.append(self.cwl_model_algorithm(
                self.cwl_params, com_channel_network))
        return results

    def _run_multiprocess(self):
        with mp.Pool(processes=self.cwl_params.ncpus) as pool:
            results = pool.starmap(self.cwl_model_algorithm, tqdm(
                [(self.cwl_params, com_channel_network) for com_channel_network in self.component_channel_networks]))

        return results

    def run(self):

        if self.cwl_params.ncpus == 1:
            results = self._run_singleprocess()
        elif self.cwl_params.ncpus > 1:
            results = self._run_multiprocess()

        # stores results in class attributes y_solution and, if preissmann, Q_solution
        self._store_results(results)

        return None

    def convert_y_nodedict_to_zeta_at_canals(self, y_nodedict):
        # Compute zeta = y - dem
        dem_dict = nx.get_node_attributes(self.cn.graph, 'DEM')
        zeta_at_canals = {node: y_value -
                          dem_dict[node] for node, y_value in y_nodedict.items()}
        return zeta_at_canals
    
    def convert_zeta_nodedict_to_y_at_canals(self, zeta_nodedict):
        # Compute zeta = y - dem
        dem_dict = nx.get_node_attributes(self.cn.graph, 'DEM')
        y_at_canals = {node: zeta_value +
                          dem_dict[node] for node, zeta_value in zeta_nodedict.items()}
        return y_at_canals
    
    def predict_q_for_next_timestep(self, theta_difference, seconds_per_timestep):
        # Returns volumetric lateral inflow per unit length that goes into canal in time dt (m^2/s)
        # theta difference (and not h difference!)
        # gives change in cell water volume per unit area between t and t + dt (m/timestep)
        # B (m) is a geometric quantity derived from geometric reasoning of canal and triangle shapes.
        # seconds_per_timestep changes units from timesteps to seconds.
        return self.cn.B * theta_difference / seconds_per_timestep
        
    def _find_zeta_in_ponding_canal_nodes(self, y_canal_nodedict):
        zeta_at_canals = self.convert_y_nodedict_to_zeta_at_canals(y_canal_nodedict)
        ponding_canal_nodes = {node:zeta for node, zeta in zeta_at_canals.items() if zeta > 0}
        return ponding_canal_nodes



class PreissmanModel(AbstractChannelHydrology):
    def __init__(self, cn: ChannelNetwork, cwl_params: CWLHydroParameters):
        super().__init__(cn, cwl_params)

        # This is the cwl algorithm
        self.cwl_model_algorithm = math_preissmann.simulate_one_component

    def _convert_solution_to_dictionary_of_nodes(self, results):
        y_results = {}
        Q_results = {}
        for channel_network, result in zip(self.component_channel_networks, results):
            y_sol = channel_network.from_nparray_to_nodedict(result[0])
            y_results = cwl_utilities.merge_two_dictionaries(y_results, y_sol)
            Q_sol = channel_network.from_nparray_to_nodedict(result[1])
            Q_results = cwl_utilities.merge_two_dictionaries(Q_results, Q_sol)

        return y_results, Q_results
    
    def _store_results(self, results):
        self.y_solution, self.Q_solution = self._convert_solution_to_dictionary_of_nodes(
            results)
        
        return None
    
class DiffWaveModel(AbstractChannelHydrology):
    def __init__(self, cn: ChannelNetwork, cwl_params: CWLHydroParameters, inexact_newton_raphson: bool):
        super().__init__(cn, cwl_params)

        # This is the cwl algorithm
        if inexact_newton_raphson:
            self.cwl_model_algorithm = math_diff_wave.simulate_one_component_inexact_newton_raphson
        if not inexact_newton_raphson:
            self.cwl_model_algorithm = math_diff_wave.simulate_one_component

    def _convert_solution_to_dictionary_of_nodes(self, results):
        y_results = {}
        for channel_network, result in zip(self.component_channel_networks, results):
            y_sol = channel_network.from_nparray_to_nodedict(result)
            y_results = cwl_utilities.merge_two_dictionaries(y_results, y_sol)

        return y_results
    
    def _store_results(self, results):
        self.y_solution = self._convert_solution_to_dictionary_of_nodes(
            results)
        
        return None
    
    
#%%

def set_up_channel_hydrology(model_type, cn: ChannelNetwork, cwl_params: CWLHydroParameters):
    """Returns the right subclass of the abstract CWLHydroParameters class

    Args:
        cwl_params (CWLHydroParameters): Parameters for the CWL computation
        cn (ChannelNetwork): Channel network stuff
        model_type (str, optional): Possible values: 'preissmann', 'diff-wave-implicit', 'diff-wave-explicit'

    Returns:
        [class]: subclass of CWLHydroParameters class
    """
    # cwl_model (str, optional):
    if model_type == 'preissmann':
        raise Warning('You have selected the Preissmann solver for the canal water level. It has fewer features than the others, and the code might not run.')
        return PreissmanModel(cn=cn, cwl_params=cwl_params)
    elif model_type == 'diff-wave-implicit':
        raise Warning('You have selected the exact Newton method for the solution of the diffusive wave approximation of the open channel flow equations. Expect this to be slower than the inexact counterpart.')
        return DiffWaveModel(cn=cn, cwl_params=cwl_params, inexact_newton_raphson=False)
    elif model_type == 'diff-wave-implicit-inexact':
        return DiffWaveModel(cn=cn, cwl_params=cwl_params, inexact_newton_raphson=True)
    elif model_type == 'diff-wave-explicit':
        raise NotImplementedError
    else:
        raise ValueError("I don't recognize that cwl hydro model type")
