
# Import GASTLI thermal module
import gastli.Thermal_evolution as therm
import gastli.constants as cte
# Other Python modules
import numpy as np
import matplotlib.pyplot as plt
# Create thermal evolution class object
my_therm_obj = therm.thermal_evolution()
# Input for interior
M_P = 100.     # Earth units
# Equilibrium temperatures
Teqpl = 700.
# Core mass fraction
CMF = 0.2
log_FeH = 1.
Tint_array = np.asarray([50., 100., 200., 300.,400., 500., 600., 700., 800.])
P_surf_array = np.asarray([1e3, 1e3, 1e3, 1e3, 9.5, 9.5, 9.5, 9.5, 9.5])

my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, log_FeH=log_FeH,P_surf=P_surf_array)

my_therm_obj.solve_thermal_evol_eq(t_Gyr=np.linspace(2.1e-6, 15., 10000))

'''
Plots
'''

# Plot thermal evolution
fig = plt.figure(figsize=(19, 6))
ax = fig.add_subplot(1, 3, 1)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.S_solution/1e6, linestyle='solid', color='black',linewidth=4)
plt.plot(my_therm_obj.age_points, my_therm_obj.s_top_TE/1e6, 'o', color='grey', markeredgecolor='k')



plt.ylabel(r'Entropy [MJ/kg/K]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

ax.set_xscale('log')
#plt.xlim(1e-2,12.)
plt.ylim(0.03,0.12)



ax = fig.add_subplot(1, 3, 2)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.Tint_solution, linestyle='solid', color='black',linewidth=4)
plt.plot(my_therm_obj.age_points, my_therm_obj.Tint_array, 'o', color='grey', markeredgecolor='k')



plt.ylabel(r'T$_{int}$ [K]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

ax.set_xscale('log')
#plt.xlim(1e-2,12.)
plt.ylim(0.,1000)




ax = fig.add_subplot(1, 3, 3)

sigma = 5.67e-8
Lsun = 3.846e26
Lint = (4 * np.pi * sigma * (my_therm_obj.Rtot_TE*11.2*cte.constants.r_e)**2 * my_therm_obj.Tint_array**4)/Lsun
Lsolution = (4 * np.pi * sigma * (my_therm_obj.Rtot_solution*11.2*cte.constants.r_e)**2 *\
             my_therm_obj.Tint_solution**4)/Lsun



plt.plot(my_therm_obj.t_Gyr, Lsolution, linestyle='solid', color='black',linewidth=4)
plt.plot(my_therm_obj.age_points, Lint, 'o', color='grey', markeredgecolor='k')



plt.ylabel(r'Luminosity [L$_{sun}$]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

ax.set_yscale('log')
ax.set_xscale('log')
#plt.xlim(1e-3,12.)


fig.savefig('thermal_evolution_concatenate.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)
