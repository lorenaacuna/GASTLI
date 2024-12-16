

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
Tintpl = 700.
Teqpl = 1000.

# Core mass fraction
CMF = 0.5


# Case 1, log(Fe/H) is known
# You must have FeH_flag=True, which is the default value
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=0.55, log_FeH=2.4,Rguess=6.,P_surf=9.5)



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

