

# Import coupling module
import gastli.Coupling as cpl
# Other Python modules
import numpy as np
import os


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

print("")
print("Case 2, Zenv is known")
# Case 2, Zenv is known
# You must have FeH_flag=False
# Envelope metal mass fraction
Zenvpl = 0.013
my_coupling.main(M_P, CMF, Teqpl, Tintpl, FeH_flag=False, CO=CO_planet, Zenv=Zenvpl)

# Composition input
print("Zenv (input) = ",my_coupling.myatmmodel.Zenv_pl)
print("C/O atm (input) = ",my_coupling.myatmmodel.CO)


# Output
print("log(Fe/H) atm [x solar] (output) = ",my_coupling.myatmmodel.log_FeH)
print("Total planet mass M [M_earth] = ",my_coupling.Mtot)
print("Temperature at 1000 bar [K] = ",my_coupling.T_surf)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup)
print("log10_g: Planet surface gravity (1000 bar) [cm/s2] = ",np.log10(my_coupling.g_surf_planet))
print("Total planet radius [R_jup] = ",my_coupling.Rtot)

tmm = my_coupling.Mtot*CMF + my_coupling.Mtot*(1-CMF)*my_coupling.myatmmodel.Zenv_pl
print("Total metal mass [M_earth] = ",tmm)

