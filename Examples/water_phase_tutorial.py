



from matplotlib import pyplot as plt
import numpy as np
import gastli.water_curves as water_curv
import gastli.Coupling as cpl

my_coupling = cpl.coupling()


# Input for interior
# Mearth units
M_P = 50.
# Internal temperature
Tintpl = 50.
# Equilibrium temperature
Teqpl = 300.
# Core mass fraction
CMF = 0.5

# Envelope log-metallicity is solar
log_FeHpl = 2.4
# C/O ratio is solar
CO_planet = 0.55
# Run coupled interior-atmosphere model
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=CO_planet, log_FeH=log_FeHpl,Rguess=6.)

my_coupling_hot = cpl.coupling()
my_coupling_hot.main(M_P, CMF, 1000., Tintpl, CO=CO_planet, log_FeH=log_FeHpl,Rguess=6.)


water_phase_lines = water_curv.water_curves()


fig,ax = plt.subplots(nrows=1,ncols=1)

plt.title(r'M = 50 $M_{\oplus}$, CMF = 0.5, $T_{int}$ = 50 K, [Fe/H] = 250 x solar')
water_phase_lines.plot_water_curves(ax)

plt.plot(my_coupling.myplanet.T, my_coupling.myplanet.P, '-', color='blue',label=r'$T_{eq}$ = 300 K')
plt.plot(my_coupling.myatmmodel.T_ode, my_coupling.myatmmodel.P_ode, '-', color='blue')

plt.plot(my_coupling_hot.myplanet.T, my_coupling_hot.myplanet.P, '-', color='red',label='$T_{eq}$ = 1000 K')
plt.plot(my_coupling_hot.myatmmodel.T_ode, my_coupling_hot.myatmmodel.P_ode, '-', color='red')


plt.yscale('log')
plt.xscale('log')

plt.ylabel(r'Pressure [Pa]', fontsize=14)
plt.xlabel(r'Temperature [K]',fontsize=14)

xmin = 100
xmax = 2e4

plt.xlim((xmin,xmax))
plt.ylim((1,1e15))

plt.text(1000, 1e9, 'Supercritical')
plt.text(400, 5e10, 'Ice VII')
plt.text(300, 5e7, 'Liquid')
plt.text(1000, 100, 'Vapor')

plt.legend()

fig.savefig('phase_diagram.pdf',bbox_inches='tight',format='pdf', dpi=1000)
plt.close(fig)
