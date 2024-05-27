

# Import coupling module
from gastli.atm_models_interp import atm_models_interp

# Other Python modules
import numpy as np
import matplotlib.pyplot as plt

path_input = "/Users/acuna/Desktop/gastli_input_data/"

# Test initialization
myatmmodel = atm_models_interp(path_to_file=path_input)

"""
# Test function 1
Tint_test = 150.
g_surf_test = 1020.
Teq_test = 110.
CO_test = 0.55
log_FeH_test = 0.
myatmmodel.calc_interior_mass_fraction(Tint_test,g_surf_test,Teq_test,CO_test,log_FeH_test)

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)


plt.plot(myatmmodel.MMF_profile, myatmmodel.Pprofile, '-', color='black', linewidth=4.)


plt.ylabel(r'Pressure [bar]', fontsize=14)
plt.xlabel(r'Metal mass fraction, $Z_{env}$', fontsize=14)

ax.invert_yaxis()
ax.set_yscale('log')

plt.xlim(0.,1.)
plt.ylim(1e3,1e-5)

fig.savefig('metal_mass_fraction_profile.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)


# Test function 2, case 1: log(Fe/H) is known

myatmmodel.calc_PTprofile(Tint_test,g_surf_test,Teq_test)
print("Zenv = ",myatmmodel.Zenv_pl)
print("Tsurf [K] = ",myatmmodel.Tsurf)


fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)

plt.plot(myatmmodel.Tprofile, myatmmodel.Pprofile, '-', color='black', linewidth=4.)


plt.ylabel(r'Pressure [bar]', fontsize=14)
plt.xlabel(r'Temperature [K]', fontsize=14)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,1e-5)
plt.xlim(0.,2000.)

fig.savefig('PT_profile.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)

# Test function 3, case 1
myatmmodel.calc_thickness(1., 1e-3)
print("Total radius [R_Jup] = ",myatmmodel.total_radius)
print("Atmospheric thickness [R_Jup] = ",myatmmodel.z_ode)


"""

# Test function 2, case 2: Zenv is known, and log(Fe/H) is calculated internally to interpolate
# the PT profile
# Default C/O in this case is 0.55, but you can specify any other value you want

Tint_test = 150.
g_surf_test = 1020.
Teq_test = 110.


myatmmodel.calc_PTprofile(Tint_test,g_surf_test,Teq_test,Zenv=0.10,FeH_flag=False)
print("Zenv = ",myatmmodel.Zenv_pl)
print('metal. [x solar] =', myatmmodel.log_FeH)
print('C/O =', myatmmodel.CO)
print("Tsurf [K] = ",myatmmodel.Tsurf)


fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)

plt.plot(myatmmodel.Tprofile, myatmmodel.Pprofile, '-', color='black', linewidth=4.)


plt.ylabel(r'Pressure [bar]', fontsize=14)
plt.xlabel(r'Temperature [K]', fontsize=14)

ax.invert_yaxis()
ax.set_yscale('log')

plt.ylim(1e3,1e-5)
plt.xlim(0.,2000.)

fig.savefig('PT_profile_Zenv.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)

# Test function 3, case 2
myatmmodel.calc_thickness(1., 1e-3)
print("Total radius [R_Jup] = ",myatmmodel.total_radius)
print("Atmospheric thickness [R_Jup] = ",myatmmodel.z_ode)


"""
# OLD TESTS
# Test
myatmmodel = atm_models_interp()
#(T_int,g_surf,Zenv,Teq):
myatmmodel.calc_PTprofile(150.,1020.,0.001,110.)
print(myatmmodel.Tsurf)
print(myatmmodel.Psurf)
#(Rbulk,Matm_earthunits):
myatmmodel.calc_thickness(1.,1e-3)
print(myatmmodel.total_radius)
print(myatmmodel.z_ode)

myatmmodel = atm_models_interp()
myatmmodel.calc_PTprofile(670.,10**2.8,0.03,110.)

myatmmodel = atm_models_interp()
myatmmodel.calc_PTprofile(150.,1000.,0.03,105.)
myatmmodel.calc_thickness(0.967983603432466,0.03814827028578978)

"""
