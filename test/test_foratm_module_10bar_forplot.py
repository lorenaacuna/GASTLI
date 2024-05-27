

# Import coupling module
from gastli.atm_models_interp import atm_models_interp

# Other Python modules
import numpy as np
import matplotlib.pyplot as plt

path_input = "/Users/acuna/Desktop/gastli_input_data/"

# Test function
Tints = np.array([50.,150.,250.,350.,449.,550., 650.,750.,850.,950.])
g_surf_test = 1e3
Teq_test = 300.
CO_test = 0.55
log_FeH_test = 2.4


for i_Tint, Tint in enumerate(Tints):
    file_name = "10bar_PTprofile_" + "{:2.0f}".format(Tint) + "K.dat"

    # Test initialization
    myatmmodel = atm_models_interp(path_to_file=path_input, name_grid="gastli_atm_grid_10bar.hdf5")

    myatmmodel.calc_interior_mass_fraction(Tint, g_surf_test, Teq_test, CO_test,\
                                                            log_FeH_test,P_surf=9.5)

    myatmmodel.calc_PTprofile(Tint,g_surf_test,Teq_test,P_surf=9.5)

    #print("Tsurf [K] = ",myatmmodel.Tsurf)

    data = np.zeros((len(myatmmodel.Pprofile),2))
    data[:,0] = myatmmodel.Pprofile
    data[:,1] = myatmmodel.Tprofile

    fmt = '%1.4e'
    np.savetxt(file_name, data, fmt=fmt)
