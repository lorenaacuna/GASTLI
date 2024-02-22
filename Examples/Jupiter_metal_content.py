

# Import coupling module
import gastli.Coupling as cpl

# Other Python modules
import numpy as np
import time


# Create coupling class
# Time initialization time
start = time.time()
my_coupling = cpl.coupling(path_to_file="/Users/acuna/Desktop/gastli_input_data/")

end = time.time()
print('Initialization time [s] = ', end-start)
print('')


# Input for interior
# 1 Mjup in Mearth units
Mjup = 318.
M_P = Mjup

# Internal and equilibrium temperatures
Tintpl = 150.
Teqpl = 105.

# Core mass fraction
CMF = 0.01


# Call to coupled interior-atm. model (and time it)
start = time.time()

# Case 1, log(Fe/H) is known
# You must have FeH_flag=True, which is the default value
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=0.55, log_FeH=0.)

end = time.time()
print('Model computation time [s] = ', end-start)
print('')

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
my_coupling.main(M_P, CMF, Teqpl, Tintpl, FeH_flag=False, CO=0.55, Zenv=0.04149)

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

