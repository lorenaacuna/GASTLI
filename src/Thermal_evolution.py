

# Import coupling module
import Coupling as cpl
# Other Python modules
import numpy as np
import math
from scipy import integrate
from scipy import interpolate
from scipy.integrate import odeint




class thermal_evolution:
    """ Class defining objects to calculate the thermal evolution for a constant mass, composition and equilibrium
    temperature

    Args:
    :param pow_law_formass (Optional): power exponent for planet radius estimation in the interior model.
           Default is 0.32. Increase if planet is very massive (greater than 5 Jupiter masses aprox.)

    """

    def __init__(self,pow_law_formass=None):

        self.pow_law_formass = pow_law_formass

        # Create coupling class
        if self.pow_law_formass == None:
            self.my_coupling = cpl.coupling()
        else:
            self.my_coupling = cpl.coupling(pow_law_formass = self.pow_law_formass)



        # Constants for thermal evolution
        self.sigma = 5.67051e-8              # W/m2/K4
        self.R_jup = 11.2
        self.R_earth = 6.371e6               # m
        self.s_conv = 365 * 24 * 3600 * 1e9
        self.m_h = 1.660e-27                 # In kg
        self.kb = 1.380649e-23               # J
        self.M_earth = 5.972e24              # kg





    def main(self,M_P,x_core,Zenv,Teq,Tint_array,Tguess=2000.,tolerance=1e-3,\
             t_Gyr=np.linspace(2.1e-6, 10., 100),S0=12.):
        """
        Function that runs compute a series of interior structure models and solves the entropy differential equation
        to calculate the thermal evolution
        Args:
        :param M_P: Planet mass in Earth masses
        :param x_core: Core mass fraction
        :param Zenv: Metal mass fraction in envelope
        :param Teq: Equilibrium temperature in K (at zero Bond albedo)
        :param Tint_array: Array containing the internal temperatures at which the entropy (and age) are to be computed
        :param Tguess (optional): Initial guess for the surface temperature in K. Default is 2000 K.
        :param tolerance (optional): relative difference in radius between interior-atm. steps. Default is 0.001
        :param t_Gyr (optional): time array in Gyr used to solve the entropy differential equation.
               Default has 100 points between 2100 yrs and 10 Gyr.
        :param S0 (optional): Initial condition for the entropy, S(t=0) in k_b*m_h units. Default is 12.

        :return:
        self.Tint_solution: array containing the internal temperature in K corresponding to the ages in t_Gyr
        self.S_solution: array containing the entropy in J/kg/K corresponding to the ages in t_Gyr
        self.s_top_TE: array containing the entropy at 1000 bar of the interior models calculated as part of
                       the thermal evolution sequence. Units: J/kg/K.
        self.s_mean_TE: array containing the mean envelope entropy of the interior models calculated as part of
                       the thermal evolution sequence. Units: J/kg/K.
        self.Mtot_TE: array containing the total mass of the interior models calculated as part of
                       the thermal evolution sequence. Units: Earth masses.
        self.Rtot_TE: array containing the total radius of the interior models calculated as part of
                       the thermal evolution sequence. Units: Jupiter radii.
        self.Rbulk_TE: array containing the bulk radius of the interior models calculated as part of
                       the thermal evolution sequence. Units: Jupiter radii.
        self.Tsurf_TE: array containing the temperature at 1000 bar of the interior models calculated as part of
                       the thermal evolution sequence. Units: K.
        self.L_TE: array containing the luminosity of the interior models calculated as part of
                       the thermal evolution sequence. Units: J/s.
        self.f_S: array containing the entropy derivative dS/dt of the interior models calculated as part of
                       the thermal evolution sequence. In SI units (J/kg/K per s).
        """

        self.M_P = M_P
        self.x_core = x_core
        self.Zenv = Zenv
        self.Teq = Teq
        self.Tint_array = Tint_array
        self.t_Gyr = t_Gyr
        self.S0 = S0


        n_therm = len(self.Tint_array)

        self.s_top_TE = np.zeros(n_therm)
        self.s_mean_TE = np.zeros(n_therm)
        self.Mtot_TE = np.zeros(n_therm)
        self.Rtot_TE = np.zeros(n_therm)
        self.Rbulk_TE = np.zeros(n_therm)
        self.Tsurf_TE = np.zeros(n_therm)
        self.L_TE = np.zeros(n_therm)
        self.f_S = np.zeros(n_therm)


        for i in range(0,n_therm):
            # Call to interior model
            self.my_coupling.main(self.M_P, self.x_core, self.Zenv, self.Teq, self.Tint_array[i], Tguess, tolerance)

            # Entropy
            S = self.my_coupling.myplanet.entropy
            i_top = self.my_coupling.myplanet.intrf[2] - 1
            i_core_env_bound = self.my_coupling.myplanet.intrf[1] - 1

            self.s_top_TE[i] = S[i_top]
            self.s_mean_TE[i] = np.mean(S[i_core_env_bound:i_top])
            self.Mtot_TE[i] = self.my_coupling.Mtot
            self.Rtot_TE[i] = self.my_coupling.Rtot
            self.Rbulk_TE[i] = self.my_coupling.Rbulk_Rjup
            self.Tsurf_TE[i] = self.my_coupling.T_surf
            """
            print(self.Mtot_TE)
            print(self.Rtot_TE)
            print(self.Rbulk_TE)
            print(self.Tsurf_TE)
            print(self.s_mean_TE)
            """


            # Luminosity
            self.L_TE[i] = 4. * math.pi * self.sigma * (self.Rbulk_TE[i] * self.R_jup * self.R_earth) ** 2\
                   * self.Tint_array[i] ** 4


            # f = dt/dS function
            r = self.my_coupling.myplanet.r[0:i_top]
            rho = self.my_coupling.myplanet.rho[0:i_top]
            temp = self.my_coupling.myplanet.T[0:i_top]

            m = np.zeros_like(r)

            for j in range(0, i_top):
                r_prim = r[0:j]
                rho_prim = rho[0:j]
                func_m = 4 * math.pi * r_prim ** 2 * rho_prim
                integ_m = integrate.trapezoid(func_m, x=r_prim)
                m[j] = integ_m / (self.M_P * self.M_earth)

            self.f_S[i] = (-1) * (self.M_P*self.M_earth)/self.L_TE[i] * integrate.trapezoid(temp, x=m)


        # dS/dt
        inv_f_S = 1 / self.f_S
        dSdt_func_interp = interpolate.interp1d(self.s_mean_TE, inv_f_S, bounds_error=False, fill_value="extrapolate")

        def dSdt_func(y, t):
            '''
            Calculates the derivative dS/dt
            :param y: entropy
            :param t: time (function not dependent on it)
            :return:
            '''
            dSdt = dSdt_func_interp(y)
            return dSdt

        # Conversion from Gyr to seconds
        t = self.t_Gyr * self.s_conv

        # Initial condition
        y0 = self.S0 * self.kb / self.m_h

        # ODE solver
        self.S_solution = odeint(dSdt_func, y0, t)

        Tint_func = interpolate.interp1d(self.s_mean_TE, self.Tint_array, bounds_error=False, fill_value="extrapolate")
        self.Tint_solution = Tint_func(self.S_solution)

















