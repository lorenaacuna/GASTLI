
# Import GASTLI thermal module
import gastli.Thermal_evolution as therm
# Other Python modules
import numpy as np
# Create thermal evolution class object
my_therm_obj = therm.thermal_evolution()
# Input for interior
M_P = 18.76     # Earth units
# Equilibrium temperatures
Teqpl = 1000.
# Core mass fraction
CMF = 0.5
log_FeH = np.log10(20.) # 20 x solar
Tint_array = np.asarray([50. ,60. ,70. ,80., 100., 110., 120., 130., 140., 160., 150., 160., 200., 240., 300.])
# Run sequence of interior models at different internal temperatures
my_therm_obj.main(M_P, CMF, Teqpl, Tint_array, log_FeH=log_FeH)
# Recommended: save sequence of interior models in case thermal evol eq. solver stops
Rjup = 11.2  # Jupiter radius in Earth units
data = np.zeros((len(my_therm_obj.f_S) ,7))
data[: ,0] = my_therm_obj.f_S
data[: ,1] = my_therm_obj.s_mean_TE
data[: ,2] = my_therm_obj.s_top_TE
data[: ,3] = my_therm_obj.Tint_array
data[: ,4] = my_therm_obj.Rtot_TE *Rjup
data[: ,5] = my_therm_obj.Rbulk_TE *Rjup
data[: ,6] = my_therm_obj.Tsurf_TE *Rjup
fmt = '%1.4e' ,'%1.4e' ,'%1.4e' ,'%1.4e' ,'%1.4e' ,'%1.4e' ,'%1.4e'
np.savetxt('thermal_sequence_HATP26b_CMF50_20xsolar.dat', data
           ,header='f_S s_mean_TE s_top_TE Tint Rtot Rbulk    Tsurf' ,comments='' ,fmt=fmt)
