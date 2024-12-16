

# Import coupling module
import gastli.Coupling as cpl
# Other Python modules
import numpy as np
# Import constants module
import gastli.constants as cte
import matplotlib.pyplot as plt


# Initialise coupling class
my_coupling = cpl.coupling()

# Input for interior
# 1 Mjup in Mearth units
M_P = 318.
# Internal temperature
Tintpl = 99.
# Equilibrium temperature
Teqpl = 122.
# Core mass fraction
CMF = 0.01

# Envelope log-metallicity is solar
log_FeHpl = 0
# C/O ratio is solar
CO_planet = 0.55
# Run coupled interior-atmosphere model
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=CO_planet, log_FeH=log_FeHpl)

print("Case 1, log(Fe/H) is known")
# Composition input
print("log(Fe/H) atm [x solar] (input) = ",my_coupling.myatmmodel.log_FeH)
print("C/O atm (input) = ",my_coupling.myatmmodel.CO)
# Output
print("Zenv (output) = ",my_coupling.myatmmodel.Zenv_pl)
print("Total planet mass M [M_earth] = ",my_coupling.Mtot)
print("Temperature at 1000 bar [K] = ",my_coupling.T_surf)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup)
print("log10_g: Planet surface gravity (1000 bar) [cm/s2] = ",np.log10(my_coupling.g_surf_planet))
print("Total planet radius [R_jup] = ",my_coupling.Rtot)
tmm = my_coupling.Mtot*CMF + my_coupling.Mtot*(1-CMF)*my_coupling.myatmmodel.Zenv_pl
print("Total metal mass [M_earth] = ",tmm)

'''
Plots
'''


# Jupiter radius in Earth radii
Rjup_Rearth = 11.2
xmax = Rjup_Rearth*my_coupling.Rbulk_Rjup

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



# Plot atm. profiles
fig = plt.figure(figsize=(24, 6))
ax = fig.add_subplot(1, 4, 1)

plt.plot(my_coupling.myatmmodel.T_ode,my_coupling.myatmmodel.P_ode/1e5, '-', color='black')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel(r'Temperature [K]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)


ax = fig.add_subplot(1, 4, 2)

plt.plot(my_coupling.myatmmodel.rho_ode,my_coupling.myatmmodel.P_ode/1e5, '-', color='blue')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel(r'Density [kg/m$^{3}$]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)

ax = fig.add_subplot(1, 4, 3)

plt.plot(my_coupling.myatmmodel.g_ode,my_coupling.myatmmodel.P_ode/1e5, '-', color='orange')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel('Gravity acceleration [m/s$^{2}$]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)


ax = fig.add_subplot(1, 4, 4)

# Rjup = 7.149e7    # Jupiter radius in m

plt.plot(my_coupling.myatmmodel.r/7.149e7,my_coupling.myatmmodel.P_ode/1e5, '-', color='red')

plt.ylabel(r'Pressure [bar]', fontsize=16)
plt.xlabel('Radius [$R_{Jup}$]', fontsize=16)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,2e-2)


fig.savefig('atmospheric_profiles.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)

