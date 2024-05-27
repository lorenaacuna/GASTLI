

# Import coupling module
import gastli.Coupling as cpl

# Other Python modules
import numpy as np
import time
import matplotlib.pyplot as plt


# Path to input files
# Dont forget the "/" at the end of the string
path_input = "/Users/acuna/Desktop/gastli_input_data/"

my_coupling = cpl.coupling(path_to_file="/Users/acuna/Desktop/gastli_input_data/")


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

Tintpl = 50.

# Case 2, Zenv is known
my_coupling.main(M_P, CMF, Teqpl, Tintpl, FeH_flag=False, Zenv=Zenv, limit_level=7.)

# Output
print("log(Fe/H) atm = ",my_coupling.myatmmodel.log_FeH)
print("C/O atm = ",my_coupling.myatmmodel.CO)
print("Zenv (for interior) = ",my_coupling.myatmmodel.Zenv_pl)
print("Total planet mass M [M_earth] = ",my_coupling.Mtot)
print("Temperature at 1000 bar [K] = ",my_coupling.T_surf)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup)
print("log10_g: Planet surface gravity (1000 bar) [cm/s2] = ",np.log10(my_coupling.g_surf_planet))
print("Total planet radius [R_jup] = ",my_coupling.Rtot)

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

#print(my_coupling.g_atm_profile)

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

