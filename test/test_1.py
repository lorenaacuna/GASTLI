

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
CMF = 0.01
Zenv = 0.04149


# Call to coupled interior-atm. model (and time it)
start = time.time()

#main(self,M_P,x_core,Teq,Tint,CO=0.55,log_FeH=0.,Zenv=0.03,FeH_flag=True,Tguess=2000.,Rguess=11.2,\
#             tolerance=1e-3):

# Case 1, log(Fe/H) is known
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=0.55,log_FeH=0.)

end = time.time()
print('Model computation time [s] = ', end-start)
print('')

# Output
print("log(Fe/H) atm = ",my_coupling.myatmmodel.log_FeH)
print("C/O atm = ",my_coupling.myatmmodel.CO)
print("Zenv (for interior) = ",my_coupling.myatmmodel.Zenv_pl)
print("Total planet mass M [M_earth] = ",my_coupling.Mtot)
print("Temperature at 1000 bar [K] = ",my_coupling.T_surf)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup)
print("log10_g: Planet surface gravity (1000 bar) [cm/s2] = ",np.log10(my_coupling.g_surf_planet))

print("Total planet radius [R_jup] = ",my_coupling.Rtot)

"""
# Case 2, Zenv is known
my_coupling.main(M_P, CMF, Teqpl, Tintpl,FeH_flag=False,Zenv=Zenv)

end = time.time()
print('Model computation time [s] = ', end-start)
print('')

# Output
print("log(Fe/H) atm = ",my_coupling.myatmmodel.log_FeH)
print("C/O atm = ",my_coupling.myatmmodel.CO)
print("Zenv (for interior) = ",my_coupling.myatmmodel.Zenv_pl)
print("Total planet mass M [M_earth] = ",my_coupling.Mtot)
print("Temperature at 1000 bar [K] = ",my_coupling.T_surf)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup)
print("log10_g: Planet surface gravity (1000 bar) [cm/s2] = ",np.log10(my_coupling.g_surf_planet))

print("Total planet radius [R_jup] = ",my_coupling.Rtot)

tmm = my_coupling.Mtot*CMF + my_coupling.Mtot*(1-CMF)*my_coupling.myatmmodel.Zenv_pl
print("Total metal mass [M_earth] = ",tmm)
"""
"""
# OLD TESTS
# Tests
my_coupling = coupling()
#(M_P,x_core,Zenv,Teq,Tint,Tguess=2000.,tolerance=1e-3)
my_coupling.main(318., 0.05, 0.04, 900., 50.)
print("Total planet radius [R_jup] = ",my_coupling.Rtot)
"""