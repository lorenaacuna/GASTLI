

# Import coupling module
import gastli.Thermal_evolution as therm

# Other Python modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd

# Path to input files
# Dont forget the "/" at the end of the string
path_input = "/Users/acuna/Desktop/gastli_input_data/"

# Create thermal evolution class
my_therm_obj = therm.thermal_evolution(path_to_file=path_input)

# Read in data saved in step 1
data = pd.read_csv('thermal_sequence_Jupiter.dat', sep='\s+',header=None,skiprows=1)
my_therm_obj.f_S = data[0]
my_therm_obj.s_mean_TE = data[1]
my_therm_obj.s_top_TE = data[2]
my_therm_obj.Tint_array = data[3]
my_therm_obj.Rtot_TE = data[4]
my_therm_obj.Rbulk_TE = data[5]
my_therm_obj.Tsurf_TE = data[6]

my_therm_obj.solve_thermal_evol_eq()

'''
Plots
'''

# Plot thermal evolution
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(2, 2, 1)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.S_solution/1e6, '-', color='black', linewidth=4.)
plt.plot(my_therm_obj.age_points, my_therm_obj.s_top_TE/1e6, 'o', color='grey', markeredgecolor='k')
plt.xlim(0.,15.)


plt.ylabel(r'Entropy [MJ/kg/K]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)


ax = fig.add_subplot(2, 2, 2)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.Tint_solution, '-', color='orange', linewidth=4.)
plt.plot(my_therm_obj.age_points, my_therm_obj.Tint_array, 'o', color='gold', markeredgecolor='k')


plt.ylabel(r'Internal temperature [K]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

plt.xlim(0.,15.)
plt.ylim(0.,1000.)


ax = fig.add_subplot(2, 2, 4)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.Rtot_solution, '-', color='blue', linewidth=4.)
plt.plot(my_therm_obj.age_points, my_therm_obj.Rtot_TE, 'o', color='dodgerblue', markeredgecolor='k')


plt.ylabel(r'Radius [$R_{Jup}$]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

plt.xlim(0.,15.)



ax = fig.add_subplot(2, 2, 3)

plt.plot(my_therm_obj.Tint_array, my_therm_obj.Rtot_TE, 'o-', color='limegreen', markeredgecolor='k')


plt.ylabel(r'Radius [$R_{Jup}$]', fontsize=14)
plt.xlabel(r'Internal temperature [K]', fontsize=14)

#plt.ylim(0.,10.)
plt.xlim(0.,1000.)



fig.savefig('thermal_evolution.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)


# Calculate Jupiter's internal temperature at age = 4.5 Gyrs

thermal_evol_function = interpolate.interp1d(my_therm_obj.t_Gyr, my_therm_obj.Tint_solution)

Tint_Jupiter = thermal_evol_function(4.5)
print("Jupiter's Tint [K] @ 4.5 Gyr = ", Tint_Jupiter)
