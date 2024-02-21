

# Import coupling module
import gastli.Thermal_evolution as therm

# Other Python modules
import numpy as np
import matplotlib.pyplot as plt

# Path to input files
# Dont forget the "/" at the end of the string
path_input = "/Users/acuna/Desktop/gastli_input_data/"

# Create thermal evolution class
my_therm_obj = therm.thermal_evolution(path_to_file=path_input)


# Input for interior
# 1 Mjup in Mearth units
Mjup = 318.
M_P = Mjup

# Equilibrium temperatures
Teqpl = 105.

# Core mass fraction
CMF = 0.01
Zenv = 0.04149

Tint_array = np.asarray([50.,90.,100.,110.,120.,130.,140.,150.,160.,200.,240.,340.,440.,540.,640.,740.,\
                         840.])
tolerance_array = np.asarray([1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,\
                              0.007,0.018,0.018,0.018,0.018])

# Run sequence of interior models at different internal temperatures
my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, FeH_flag=False, Zenv=Zenv, tolerance=tolerance_array)


# Recommended: save sequence of interior models in case thermal evol eq. solver stops
data = np.zeros((len(my_therm_obj.f_S),6))
data[:,0] = my_therm_obj.f_S
data[:,1] = my_therm_obj.s_mean_TE
data[:,2] = my_therm_obj.Tint_array
data[:,3] = my_therm_obj.Rtot_TE
data[:,4] = my_therm_obj.Rbulk_TE
data[:,5] = my_therm_obj.Tsurf_TE

fmt = '%1.4e','%1.4e','%1.4e','%1.4e','%1.4e','%1.4e'
np.savetxt('thermal_sequence_Jupiter.dat', data,header='f_S s_mean_TE Tint Rtot Rbulk Tsurf',comments='',fmt=fmt)



