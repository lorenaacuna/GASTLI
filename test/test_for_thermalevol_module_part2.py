

# Import coupling module
import gastli.Thermal_evolution as therm

# Other Python modules
import numpy as np
import time
import matplotlib.pyplot as plt


# Path to input files
# Dont forget the "/" at the end of the string
path_input = "/Users/acuna/Desktop/gastli_input_data/"

# Create thermal evolution class
# Time initialization
start = time.time()
my_therm_obj = therm.thermal_evolution(path_to_file=path_input)

end = time.time()
print('Initialization time [s] = ', end-start)
print('')

# Input for interior
# 1 Mjup in Mearth units
Mjup = 318.
# M = 1 Mjup
M_P = Mjup
# Equilibrium temperatures
Teqpl = 105.
# Core mass fraction and envelope metal mass fraction
CMF = 0.05
Zenv = 0.04
# Internal temperature array
# Tint_array = np.linspace(110.,200.,5)

# If Tint = 600 K
#Error in Coupling.py: The number of interior-atmosphere iterations is greater than 20.
#The current relative difference between radii is 0.017744705208929703




Tint_array = np.asarray([100,500.,600.,700.,800.,849.])
tolerance_array = np.asarray([1e-3,1e-3,0.018,0.018,0.018,0.018])
#tolerance_array = np.asarray([1e-3,1e-3,1e-3,1e-3,1e-3,2e-2,2e-2,2e-2])

#def main(self, M_P, x_core, Teq, Tint_array, CO=0.55, log_FeH=0., Zenv=0.03, FeH_flag=True, Tguess=2000.,
#         tolerance=1e-3, \
#         t_Gyr=np.linspace(2.1e-6, 10., 100), S0=12.):


my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, Zenv=Zenv, FeH_flag=False,tolerance=tolerance_array)

# Recommended: save sequence of interior models in case thermal evol eq. solver fails

my_therm_obj.solve_thermal_evol_eq()

"""
my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, Zenv=Zenv, FeH_flag=False, tolerance=tolerance_array)
"""
'''
Plots
'''

# Plot thermal evolution
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(2, 2, 1)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.S_solution/1e6, '-', color='black', linewidth=4.)
plt.plot(my_therm_obj.age_points, my_therm_obj.s_mean_TE/1e6, 'o', color='grey', markeredgecolor='k')

plt.xlim(0.,10.)


plt.ylabel(r'Entropy [MJ/kg/K]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)


ax = fig.add_subplot(2, 2, 2)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.Tint_solution, '-', color='orange', linewidth=4.)
plt.plot(my_therm_obj.age_points, my_therm_obj.Tint_array, 'o', color='gold', markeredgecolor='k')


plt.ylabel(r'Internal temperature [K]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

plt.xlim(0.,10.)
plt.ylim(0.,1000.)


ax = fig.add_subplot(2, 2, 4)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.Rtot_solution, '-', color='blue', linewidth=4.)
plt.plot(my_therm_obj.age_points, my_therm_obj.Rtot_TE, 'o', color='dodgerblue', markeredgecolor='k')


plt.ylabel(r'Radius [$R_{Jup}$]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

plt.xlim(0.,10.)



ax = fig.add_subplot(2, 2, 3)

plt.plot(my_therm_obj.Tint_array, my_therm_obj.Rtot_TE, 'o-', color='limegreen', markeredgecolor='k')


plt.ylabel(r'Radius [$R_{Jup}$]', fontsize=14)
plt.xlabel(r'Internal temperature [K]', fontsize=14)

#plt.ylim(0.,10.)
plt.xlim(0.,1000.)



fig.savefig('thermal_evolution.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)


print("t_Gyr = ", my_therm_obj.t_Gyr)
print("Solution: Tint [K] = ",my_therm_obj.Tint_solution)
print("Solution: Total radius [Rjup] = ",my_therm_obj.Rtot_solution)
print("Rtot [Rjup] = ",my_therm_obj.Rtot_TE)
print("Tint [K] = ",Tint_array)
print("Rbulk [Rjup] = ",my_therm_obj.Rbulk_TE)
print("Tsurf [K] = ",my_therm_obj.Tsurf_TE)
print("Mean entropy [MJ/kg/K] = ",my_therm_obj.s_mean_TE/1e6)
print("Top entropy [MJ/kg/K] = ",my_therm_obj.s_top_TE/1e6)




