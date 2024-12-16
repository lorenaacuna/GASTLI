
# Import coupling module
import gastli.Coupling as cpl
# Other Python modules
import numpy as np

# Input for interior
## 1 Mjup in Mearth units
Mjup = 318.
mass_array = Mjup * np.arange(0.05 ,1.6 ,0.05)
n_mrel = len(mass_array)
## Internal temperature
Tintpl = 107.     # K
## Equilibrium temperature
Tstar = 5777.     # K
Rstar = 0.00465   # AU
ad = 5.2          # AU
Teq_4 = Tstar**4./4. * (Rstar /ad )**2.
Teqpl = Teq_4**0.25
# Core mass fraction
CMF = 0.
# Mass-radius curve output file
file_out = open('Jupiter_MRrel_CMF0_logFeH_0.dat' ,'w')
file_out.write('  M_int[M_E]  M_tot[M_E]  x_core  ')
file_out.write('T_surf[K]  R_bulk[R_J]  R_tot[R_J]  T_int[K]  Zenv  z_atm[R_J] ')
file_out.write("\n")
# For loop that changes the mass in each call of the coupling class
for k in range(0, n_mrel):
    M_P_model = mass_array[k]
    print('---------------')
    print('Mass [Mearth] = ', M_P_model)
    print('Model = ', k+ 1, ' out of ', n_mrel)
    print('---------------')
    # Create coupling class (this also resets parameters)
    my_coupling = cpl.coupling(pow_law_formass=0.31)
    # Case 1, log(Fe/H) is known
    # You must have FeH_flag=True, which is the default value
    my_coupling.main(M_P_model, CMF, Teqpl, Tintpl, CO=0.55, log_FeH=0.)
    # Save data
    file_out.write('%s %s' % ("  ", str(M_P_model)))
    file_out.write('%s %s' % ("  ", str(my_coupling.Mtot)))
    file_out.write('%s %s' % ("  ", str(CMF)))
    file_out.write('%s %s' % ("  ", str(my_coupling.T_surf)))
    file_out.write('%s %s' % ("  ", str(my_coupling.Rbulk_Rjup)))
    file_out.write('%s %s' % ("  ", str(my_coupling.Rtot)))
    file_out.write('%s %s' % ("  ", str(Tintpl)))
    file_out.write('%s %s' % ("  ", str(my_coupling.myatmmodel.Zenv_pl)))
    zatm_RJ = my_coupling.Rtot - my_coupling.Rbulk_Rjup
    file_out.write('%s %s' % ("  ", str(zatm_RJ)))
    file_out.write("\n")
    file_out.flush()
# End of for loops
file_out.close()
