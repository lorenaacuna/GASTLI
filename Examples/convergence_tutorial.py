

# Import coupling module
import gastli.Coupling as cpl
import gastli.constants as cte

# Other Python modules
# Other Python modules
import numpy as np
import matplotlib.pyplot as plt


# Create coupling class
# Time initialization time
my_coupling = cpl.coupling(j_max=99, pow_law_formass=0.31)



# Input for interior
# 1 Mjup in Mearth units
M_P = 50.

# Internal and equilibrium temperatures
Tintpl = 150.
Teqpl = 1000.

# Core mass fraction
CMF = 0.1


# Call to coupled interior-atm. model (and time it)


# Case 1, log(Fe/H) is known
# You must have FeH_flag=True, which is the default value
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=0.55, log_FeH=0.,Rguess=6.)



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

base = my_coupling.myplanet.intrf_hist[0,:]
core = my_coupling.myplanet.intrf_hist[1,:]
envelope = my_coupling.myplanet.intrf_hist[2,:]
surface = my_coupling.myplanet.intrf_hist[3,:]
x = my_coupling.myplanet.iter_num
r = my_coupling.myplanet.r/cte.constants.r_e

mask = core != 0

# Plot interior profiles
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)


plt.plot(x[mask], r[base[mask]-1], linestyle='solid', color='black')
plt.plot(x[mask], r[core[mask]-1], linestyle='solid', color='brown')
plt.plot(x[mask], r[envelope[mask]-1], linestyle='solid', color='deepskyblue')
plt.plot(x[mask], r[surface[mask]-2], linestyle='solid', color='grey')

ax.fill_between(x[mask], r[base[mask]-1], r[core[mask]-1], facecolor='brown',alpha=0.5)
ax.fill_between(x[mask], r[core[mask]-1], r[envelope[mask]-1], facecolor='deepskyblue',alpha=0.5)


plt.xlabel(r'Iteration #', fontsize=16)
plt.ylabel(r'Radius [$R_{\oplus}$]', fontsize=16)

plt.xlim(0, 100)


fig.savefig('convergence_tutorial.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)

