

# Import coupling module
from gastli.atm_models_interp import atm_models_interp

# Other Python modules
import numpy as np
import matplotlib.pyplot as plt

path_input = "/Users/acuna/Desktop/gastli_input_data/"

# Test initialization
myatmmodel = atm_models_interp(path_to_file=path_input)


# Test function
Tint_test = 452.
g_surf_test = 1e3
Teq_test = 300.
CO_test = 0.55
log_FeH_test = 2.4

myatmmodel.calc_interior_mass_fraction(Tint_test, g_surf_test, Teq_test, CO_test,\
                                                            log_FeH_test)

myatmmodel.calc_PTprofile(Tint_test,g_surf_test,Teq_test)

print("Tsurf [K] = ",myatmmodel.Tsurf)

