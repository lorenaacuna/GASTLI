

# Import coupling module
import gastli.Coupling as cpl
# Import constants module
import gastli.constants as cte


# Other Python modules
import numpy as np


# Create coupling class
# Time initialization time
my_coupling = cpl.coupling(path_to_file="/Users/acuna/Desktop/gastli_input_data/",\
                           name_grid="gastli_atm_grid_10bar.hdf5", j_max=99, pow_law_formass=0.30)


# Input for interior
# 1 Mjup in Mearth units
# Mjup = 318.
M_P = 30.       # 30 Mearth

# Internal and equilibrium temperatures
Tintpl = 450.
Teqpl = 300.

# Core mass fraction
CMF = 0.3


# Case 2, Zenv is known
# You must have FeH_flag=False
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=0.55, log_FeH=2.4,P_surf=9.5,Rguess=11.2*0.3)

'''
Model info
'''
file_name = 'Tint_450K_10bar'
file_out = open(file_name+'.dat','w')


# Composition input
print("Zenv (input) = ",my_coupling.myatmmodel.Zenv_pl,file=file_out)
print("C/O atm (input) = ",my_coupling.myatmmodel.CO,file=file_out)


# Output
print("log(Fe/H) atm [x solar] (output) = ",my_coupling.myatmmodel.log_FeH,file=file_out)
print("Total planet mass M [M_earth] = ",my_coupling.Mtot,file=file_out)
print("Temperature at Psurf [K] = ",my_coupling.T_surf,file=file_out)
print("Planet bulk radius [R_jup] = ",my_coupling.Rbulk_Rjup,file=file_out)
print("log10_g: Planet surface gravity (Psurf) [cm/s2] = ",np.log10(my_coupling.g_surf_planet),file=file_out)
print("Total planet radius [R_jup] = ",my_coupling.Rtot,file=file_out)

tmm = my_coupling.Mtot*CMF + my_coupling.Mtot*(1-CMF)*my_coupling.myatmmodel.Zenv_pl
print("Total metal mass [M_earth] = ",tmm,file=file_out)

'''
Interior profiles
'''

data = np.zeros((len(my_coupling.myplanet.r), 5))
data[:, 0] = my_coupling.myplanet.r/cte.constants.r_e
data[:, 1] = my_coupling.myplanet.g
data[:, 2] = my_coupling.myplanet.P
data[:, 3] = my_coupling.myplanet.T
data[:, 4] = my_coupling.myplanet.rho

press = my_coupling.myplanet.P
S = my_coupling.myplanet.entropy
index_find = np.array(np.where(press <= 1e3*1e5))
i_top = index_find[0]
a = S[i_top]

fmt = '%1.4e'
np.savetxt(file_name+'_interior_profiles.dat', data, fmt=fmt, header='r[R_E]  g[m/s2]  P[Pa]  T[K]  rho[kg/m3]')

'''
Atm profiles
'''

data = np.zeros((len(my_coupling.myatmmodel.T_ode), 5))
data[:, 0] = my_coupling.myatmmodel.P_ode
data[:, 1] = my_coupling.myatmmodel.T_ode
data[:, 2] = my_coupling.myatmmodel.rho_ode
data[:, 3] = my_coupling.myatmmodel.g_ode
data[:, 4] = my_coupling.myatmmodel.r/cte.constants.r_e

fmt = '%1.4e'
np.savetxt(file_name+'_atm_profiles.dat', data, fmt=fmt, header='P[Pa]  T[K]  rho[kg/m3]  g[m/s2]  r[R_E]')

'''
Convergence
'''

base = my_coupling.myplanet.intrf_hist[0,:]
core = my_coupling.myplanet.intrf_hist[1,:]
envelope = my_coupling.myplanet.intrf_hist[2,:]
surface = my_coupling.myplanet.intrf_hist[3,:]
x = my_coupling.myplanet.iter_num

data = np.zeros((len(base), 5))
data[:, 0] = base
data[:, 1] = core
data[:, 2] = envelope
data[:, 3] = surface
data[:, 4] = x

np.savetxt(file_name+'_convergence.dat', data, fmt='%i', header='base  core  envelope  surface  x')



