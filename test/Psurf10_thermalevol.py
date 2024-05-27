

# Import coupling module
import gastli.Thermal_evolution as therm

# Other Python modules
import numpy as np
import matplotlib.pyplot as plt

# Path to input files
# Dont forget the "/" at the end of the string
path_input = "/Users/acuna/Desktop/gastli_input_data/"

# Create thermal evolution class
my_therm_obj = therm.thermal_evolution(path_to_file=path_input,pow_law_formass=0.30)


# Input for interior
# 1 Mjup in Mearth units
#Mjup = 318.
M_P = 100.

# Equilibrium temperatures
Teqpl = 300.

# Core mass fraction
CMF = 0.1

Tint_array = np.asarray([100.,110.,120.,130.,140.,150.,160.,200.,240.,340.,440.,540.,640.,740.,\
                         840.])


'''
Tint_array = np.asarray([50.,90.,100.,110.,120.,130.,140.,150.,160.,200.,240.,340.,440.,540.,640.,740.,\
                         840.])
'''

# Run sequence of interior models at different internal temperatures
my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, CO=0.55, log_FeH=2.4,Rguess=11.2*0.4,P_surf=9.5,\
                  name_grid="gastli_atm_grid_10bar.hdf5",j_max=99)


# Recommended: save sequence of interior models in case thermal evol eq. solver stops
data = np.zeros((len(my_therm_obj.f_S),7))
data[:,0] = my_therm_obj.f_S
data[:,1] = my_therm_obj.s_mean_TE
data[:,2] = my_therm_obj.s_top_TE
data[:,3] = my_therm_obj.Tint_array
data[:,4] = my_therm_obj.Rtot_TE
data[:,5] = my_therm_obj.Rbulk_TE
data[:,6] = my_therm_obj.Tsurf_TE

fmt = '%1.10e','%1.10e','%1.10e','%1.10e','%1.10e','%1.10e','%1.10e'
np.savetxt('thermal_sequence_10bar.dat', data,header='f_S s_mean_TE s_top_TE Tint Rtot Rbulk Tsurf',\
           comments='',fmt=fmt)



