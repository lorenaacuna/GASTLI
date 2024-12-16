
# Import coupling module
import gastli.Thermal_evolution as therm
import gastli.constants as cte
# Other Python modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd

# Create thermal evolution class
my_therm_obj = therm.thermal_evolution()

# Read in data saved in step 1
data = pd.read_csv('thermal_sequence_HATP26b_CMF50_20xsolar.dat', sep='\s+',header=None,skiprows=1)
my_therm_obj.f_S = data[0]
my_therm_obj.s_mean_TE = data[1]
my_therm_obj.s_top_TE = data[2]
my_therm_obj.Tint_array = data[3]
my_therm_obj.Rtot_TE = data[4]
my_therm_obj.Rbulk_TE = data[5]
my_therm_obj.Tsurf_TE = data[6]

my_therm_obj.solve_thermal_evol_eq(t_Gyr=np.linspace(2.1e-6, 15., 10000))

# Plot thermal evolution
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)

plt.title(r"Envelope composition: 20 $\times$ solar")

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.Rtot_solution, '-', color='springgreen', linewidth=4.) #,label="CMF = 0.9")
plt.plot(my_therm_obj.age_points, my_therm_obj.Rtot_TE, 'o', color='springgreen', markeredgecolor='k')




plt.text(12.3,6.3,"CMF = 0.5")

yerr = np.zeros((2,1))
yerr[0,0] = 0.36
yerr[1,0] = 0.81


xerr = np.zeros((2,1))
xerr[0,0] = 4.9
xerr[1,0] = 3.


plt.errorbar([9.],[6.33], yerr, xerr,'X',color='black',label="HAT-P-26 b")

plt.legend()

plt.ylabel(r'Radius [$R_{\oplus}$]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

plt.xlim(0.,15.)
plt.ylim(3.,10.)


fig.savefig('thermal_evolution_HATP26b_20xsolar.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)


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
plt.xlim(1e-2,12.)
plt.ylim(0.03,0.12)



ax = fig.add_subplot(1, 3, 2)

plt.plot(my_therm_obj.t_Gyr, my_therm_obj.Tint_solution, linestyle='solid', color='black',linewidth=4)
plt.plot(my_therm_obj.age_points, my_therm_obj.Tint_array, 'o', color='grey', markeredgecolor='k')



plt.ylabel(r'T$_{int}$ [K]', fontsize=14)
plt.xlabel(r'Age [Gyrs]', fontsize=14)

ax.set_xscale('log')
plt.xlim(1e-2,12.)
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
plt.xlim(1e-2,12.)


fig.savefig('thermal_evolution_all.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)
