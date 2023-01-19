class PeatlandHydroParameters:
    def __init__(self,
                 params_peat,
                 use_several_weather_stations) -> None:
        self.dx = float(params_peat['dx']) # only used if rectangular grid
        self.max_sweeps = int(params_peat['max_sweeps']) # max number of iterations for convergence of FiPy's numerical solution
        self.fipy_desired_residual = float(params_peat['fipy_desired_residual'])
        self.zeta_diri_bc = float(params_peat['zeta_diri_BC'])

        # Compute P-ET based on weighted average distance of several weather stations
        self.use_several_weather_stations = use_several_weather_stations

        pass
