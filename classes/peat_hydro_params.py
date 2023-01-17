class PeatlandHydroParameters:
    def __init__(self,
                 params_peat,
                 dt, dx,
                 use_several_weather_stations) -> None:
        self.dt = dt
        self.dx = dx
        self.max_sweeps = int(params_peat['max_sweeps']) # max number of iterations for convergence of FiPy's numerical solution
        self.fipy_desired_residual = float(params_peat['fipy_desired_residual'])

        # parameters for S and T functions. These must be set from code.
        self.s1 = None
        self.s2 = None
        self.t1 = None
        self.t2 = None

        # Compute P-ET based on weighted average distance of several weather stations
        self.use_several_weather_stations = use_several_weather_stations

        pass
