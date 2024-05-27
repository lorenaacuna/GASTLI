

# Import interior module
from gastli.GASTLI import int_planet
import gastli.constants as cte

# Other Python modules
import numpy as np
import matplotlib.pyplot as plt

path_to_input = "/Users/acuna/Desktop/gastli_input_data/"
myplanet_test = int_planet(path_to_file=path_to_input,\
                           j_max=99, pow_law=0.30)

myplanet_test.setup_parameters()

M_P = 30.                 # mass in Mearth
CMF = 0.5                 # core mass fraction
Tsurf = 1250.             # surface temperature in K
Psurf = 1e1*1e5
#Psurf = 1e2*1e5
#Psurf = 1e3*1e5           # surface pressure in K
Zenv = 0.013              # metal mass fraction in envelope

myplanet_test.calc_radius(M_P,CMF,1-CMF,Tsurf,Psurf,Zenv)

print('Output radius [R_earth] = ', myplanet_test.R_P)


pressure = myplanet_test.P
temp = myplanet_test.T

data = np.zeros((len(pressure),2))
data[:,0] = pressure
data[:,1] = temp

fmt = '%1.4e'
np.savetxt('1000bar_test.dat', data, fmt=fmt)



base = myplanet_test.intrf_hist[0,:]
core = myplanet_test.intrf_hist[1,:]
envelope = myplanet_test.intrf_hist[2,:]
surface = myplanet_test.intrf_hist[3,:]
x = myplanet_test.iter_num
r = myplanet_test.r/cte.constants.r_e

'''
PLOT
'''

mask = core != 0

# Plot interior profiles
fig = plt.figure(figsize=(13, 6))
ax = fig.add_subplot(1, 2, 1)

plt.plot(x[mask], base[mask], linestyle='solid', color='black')
plt.plot(x[mask], core[mask], linestyle='solid', color='brown')
plt.plot(x[mask], envelope[mask], linestyle='solid', color='deepskyblue')
plt.plot(x[mask], surface[mask], linestyle='solid', color='grey')

ax.fill_between(x[mask], base[mask], core[mask], facecolor='brown',alpha=0.5)
ax.fill_between(x[mask], core[mask], envelope[mask], facecolor='deepskyblue',alpha=0.5)


plt.xlabel(r'Iteration #', fontsize=16)
plt.ylabel(r'Spatial index', fontsize=16)

plt.xlim(0, 100)
#plt.ylim(0, 1.1 * np.nanmax(my_coupling.myplanet.g))


ax = fig.add_subplot(1, 2, 2)

plt.plot(x[mask], r[base[mask]-1], linestyle='solid', color='black')
plt.plot(x[mask], r[core[mask]-1], linestyle='solid', color='brown')
plt.plot(x[mask], r[envelope[mask]-1], linestyle='solid', color='deepskyblue')
plt.plot(x[mask], r[surface[mask]-2], linestyle='solid', color='grey')

ax.fill_between(x[mask], r[base[mask]-1], r[core[mask]-1], facecolor='brown',alpha=0.5)
ax.fill_between(x[mask], r[core[mask]-1], r[envelope[mask]-1], facecolor='deepskyblue',alpha=0.5)


plt.xlabel(r'Iteration #', fontsize=16)
plt.ylabel(r'Radius [$R_{\oplus}$]', fontsize=16)

plt.xlim(0, 100)
#plt.ylim(0, 1.1 * np.nanmax(my_coupling.myplanet.g))


fig.savefig('layer_convergence_10bar_lowT.pdf', bbox_inches='tight', format='pdf', dpi=1000)
plt.close(fig)
