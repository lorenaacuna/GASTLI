

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
# M = 1 Mjup
M_P = Mjup
# Internal and equilibrium temperatures

Tintpl = 150.
Teqpl = 105.

# Core mass fraction and envelope metal mass fraction
CMF = 0.05
Zenv = 0.03


# Call to coupled interior-atm. model (and time it)
start = time.time()


my_coupling.main(M_P, CMF, 0.04, Teqpl, Tintpl)

end = time.time()
print('Model computation time [s] = ', end-start)
print('')

# Output
print("Total planet mass M [M_earth] = ",my_coupling.Mtot)
print("Temperature at 1000 bar [K] = ",my_coupling.T_surf)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup)
print("log10_g: Planet surface gravity (1000 bar) [cm/s2] = ",np.log10(my_coupling.g_surf_planet))

print("Total planet radius [R_jup] = ",my_coupling.Rtot)

