
# Import GASTLI thermal module
import gastli.Thermal_evolution as therm
import gastli.constants as cte
# Other Python modules
import numpy as np
import matplotlib.pyplot as plt
# Path to input files
# Dont forget the "/" at the end of the string
path_input = "/Users/acuna/Desktop/gastli_input_data/"
# Create thermal evolution class object
my_therm_obj = therm.thermal_evolution(path_to_file=path_input)
# Input for interior
M_P = 100.     # Earth units
# Equilibrium temperatures
Teqpl = 700.
# Core mass fraction
CMF = 0.2
log_FeH = 1.
Tint_array = np.asarray([50., 100., 200., 300.])
# Run sequence of interior models at different internal temperatures (up to 300 K)
my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, log_FeH=log_FeH)

f_S_cold = my_therm_obj.f_S
s_mean_TE_cold = my_therm_obj.s_mean_TE
s_top_TE_cold = my_therm_obj.s_top_TE
Tint_array_cold = my_therm_obj.Tint_array
Rtot_TE_cold = my_therm_obj.Rtot_TE
Rbulk_TE_cold = my_therm_obj.Rbulk_TE
Tsurf_TE_cold = my_therm_obj.Tsurf_TE


# Run sequence of interior models at different internal temperatures (from 400 K)
Tint_array = np.asarray([400., 500., 600., 700., 800.])
my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, log_FeH=log_FeH, P_surf=9.5)

f_S_hot = my_therm_obj.f_S
s_mean_TE_hot = my_therm_obj.s_mean_TE
s_top_TE_hot = my_therm_obj.s_top_TE
Tint_array_hot = my_therm_obj.Tint_array
Rtot_TE_hot = my_therm_obj.Rtot_TE
Rbulk_TE_hot = my_therm_obj.Rbulk_TE
Tsurf_TE_hot = my_therm_obj.Tsurf_TE


# Concatenate
my_therm_obj.f_S = np.concatenate((f_S_cold,f_S_hot))
my_therm_obj.s_mean_TE = np.concatenate((s_mean_TE_cold,s_mean_TE_hot))
my_therm_obj.s_top_TE = np.concatenate((s_top_TE_cold,s_top_TE_hot))
my_therm_obj.Tint_array = np.concatenate((Tint_array_cold,Tint_array_hot))
my_therm_obj.Rtot_TE = np.concatenate((Rtot_TE_cold,Rtot_TE_hot))
my_therm_obj.Rbulk_TE = np.concatenate((Rbulk_TE_cold,Rbulk_TE_hot))
my_therm_obj.Tsurf_TE = np.concatenate((Tsurf_TE_cold,Tsurf_TE_hot))

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
