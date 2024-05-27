

# Import coupling module
import gastli.Coupling as cpl
# Import constants module
import gastli.constants as cte

# Other Python modules
import numpy as np
import time
import matplotlib.pyplot as plt


# Path to input files
# Dont forget the "/" at the end of the string
path_input = "/Users/acuna/Desktop/gastli_input_data/"

#my_coupling = cpl.coupling(path_to_file="/Users/acuna/Desktop/gastli_input_data/")

my_coupling = cpl.coupling(path_to_file="/Users/acuna/Desktop/gastli_input_data/",\
                           pow_law_formass=0.31)



# Input for interior
# 1 Mjup in Mearth units
Mjup = 318.
# M = 1 Mjup
M_P = Mjup
# Equilibrium temperatures
Teqpl = 105.
# Core mass fraction and envelope metal mass fraction
CMF = 0.01

Tintpl = 140.


# Case 2, Zenv is known
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=0.55, log_FeH=0.)


# Output
print("log(Fe/H) atm = ",my_coupling.myatmmodel.log_FeH)
print("C/O atm = ",my_coupling.myatmmodel.CO)
print("Zenv (for interior) = ",my_coupling.myatmmodel.Zenv_pl)
print("Total planet mass M [M_earth] = ",my_coupling.Mtot)
print("Temperature at 1000 bar [K] = ",my_coupling.T_surf)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup)
print("log10_g: Planet surface gravity (1000 bar) [cm/s2] = ",np.log10(my_coupling.g_surf_planet))
print("Total planet radius [R_jup] = ",my_coupling.Rtot)

#print(my_coupling.myplanet.intrf)

'''
Plots
'''


# Jupiter radius in Earth radii
Rjup_Rearth = 11.2
xmax = Rjup_Rearth*my_coupling.Rtot*1.1

# Plot interior profiles
fig = plt.figure(figsize=(6, 30))
ax = fig.add_subplot(5, 1, 1)

plt.plot(my_coupling.myplanet.r / cte.constants.r_e, my_coupling.myplanet.g, '-', color='lime')

plt.xlabel(r'Radius [$R_{\oplus}$]', fontsize=16)
plt.ylabel(r'Gravity acceleration [$m/s^{2}$]', fontsize=16)

plt.xlim(0, xmax)
plt.ylim(0, 1.1 * np.nanmax(my_coupling.myplanet.g))

ax = fig.add_subplot(5, 1, 2)

plt.plot(my_coupling.myplanet.r / cte.constants.r_e, my_coupling.myplanet.P / 1e9, '-', color='blue')

plt.xlabel(r'Radius [$R_{\oplus}$]', fontsize=16)
plt.ylabel('Pressure [GPa]', fontsize=16)

plt.xlim(0, xmax)
plt.ylim(0, 1.1 * np.amax(my_coupling.myplanet.P / 1e9))

ax = fig.add_subplot(5, 1, 3)

plt.plot(my_coupling.myplanet.r / cte.constants.r_e, my_coupling.myplanet.T, '-', color='magenta')

plt.xlabel(r'Radius [$R_{\oplus}$]', fontsize=16)
plt.ylabel('Temperature [K]', fontsize=16)

plt.xlim(0, xmax)
plt.ylim(0, 1.1 * np.amax(my_coupling.myplanet.T))

ax = fig.add_subplot(5, 1, 4)

plt.plot(my_coupling.myplanet.r / cte.constants.r_e, my_coupling.myplanet.rho, '-', color='red')

plt.xlabel(r'Radius [$R_{\oplus}$]', fontsize=16)
plt.ylabel(r'Density [$kg/m^{3}$]', fontsize=16)

plt.xlim(0, xmax)
plt.ylim(0, 1.1 * np.nanmax(my_coupling.myplanet.rho))


ax = fig.add_subplot(5, 1, 5)

plt.plot(my_coupling.myplanet.r / cte.constants.r_e, my_coupling.myplanet.entropy/1e6, '-', color='black')

plt.xlabel(r'Radius [$R_{\oplus}$]', fontsize=16)
plt.ylabel(r'Entropy [MJ/kg/K]', fontsize=16)

plt.xlim(0, xmax)
plt.ylim(0, 1.1 * np.nanmax(my_coupling.myplanet.entropy/1e6))


fig.savefig('interior_structure_profiles.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)


# Plot planet core and envelope
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)

r_core = my_coupling.myplanet.r[my_coupling.myplanet.intrf[1] - 1]\
         / my_coupling.myplanet.r[my_coupling.myplanet.intrf[2] - 1]
r_lm = my_coupling.myplanet.r[my_coupling.myplanet.intrf[2] - 1]\
       / my_coupling.myplanet.r[my_coupling.myplanet.intrf[2] - 1]

circle4 = plt.Circle((0.5, 0.5), r_core, color='teal')
circle3 = plt.Circle((0.5, 0.5), r_lm, color='mediumspringgreen')

ax.add_patch(circle3)
ax.add_patch(circle4)

plt.tick_params(axis='both', which='both', bottom=False, top=False, \
                labelbottom=False, right=False, left=False, labelleft=False)

plt.axis('equal')

fig.savefig('core_and_envelope.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)

# Plot atm. profiles
fig = plt.figure(figsize=(24, 6))
ax = fig.add_subplot(1, 4, 1)

plt.plot(my_coupling.T_atm_profile,my_coupling.P_atm_profile/1e5, '-', color='black')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel(r'Temperature [K]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)


ax = fig.add_subplot(1, 4, 2)

plt.plot(my_coupling.rho_atm_profile,my_coupling.P_atm_profile/1e5, '-', color='blue')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel(r'Density [kg/m$^{3}$]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)

ax = fig.add_subplot(1, 4, 3)

plt.plot(my_coupling.g_atm_profile,my_coupling.P_atm_profile/1e5, '-', color='orange')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel('Gravity acceleration [m/s$^{2}$]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)


ax = fig.add_subplot(1, 4, 4)

# Rjup = 7.149e7    # Jupiter radius in m

plt.plot(my_coupling.r_atm_profile/7.149e7,my_coupling.P_atm_profile/1e5, '-', color='red')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel('Radius [$R_{Jup}$]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)


fig.savefig('atmospheric_profiles.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)

