
import GASTLI
import constants as cte
import atm_models_interp

import numpy as np
import time
import math
import sys

#
class coupling:
    """ Class defining objects to run a normal interior-atmosphere coupling model
        (no thermal evolution)

    Args:
        pow_law_formass: power exponent for planet radius estimation in the interior model.
         Default is 0.32. Increase if planet is very massive (greater than 5 Jupiter masses aprox.)
    """

    def __init__(self,pow_law_formass=None):
        '''
        - Initialises interior and atm. models
        - Other constants
        '''

        self.pow_law_formass = pow_law_formass

        # Initialise interior model
        if self.pow_law_formass == None:
            self.myplanet = GASTLI.int_planet()
        else:
            self.myplanet = GASTLI.int_planet(pow_law = self.pow_law_formass)


        self.myplanet.setup_parameters()

        # Initialise atm models
        self.myatmmodel = atm_models_interp.atm_models_interp()

        # Important constants
        #self.M_P = 318.             # M in M_earth
        self.Mearth = 5.972e24       # In m



    def main(self,M_P,x_core,Zenv,Teq,Tint,Tguess=2000.,tolerance=1e-3):
        """
        Function that runs coupled interior-atmosphere model
        Args:
        :param M_P: Planet mass in Earth masses
        :param x_core: Core mass fraction
        :param Zenv: Metal mass fraction in envelope
        :param Teq: Equilibrium temperature in K (at zero Bond albedo)
        :param Tint: Internal temperature in K
        :param Tguess (optional): Initial guess for the surface temperature in K. Default is 2000 K.
        :param tolerance (optional): relative difference in radius between interior-atm. steps. Default is 0.001
        :return:
        ## From interior model (only) ##
        self.myplanet.R_P: Planet bulk (from center to Psurf) radius in Earth radii
        self.myplanet.rho_p: Planet bulk density in g/cm3
        self.myplanet.OtoH: O:H ratio in envelope
        self.myplanet.metal: Envelope metallicity
        ## From int.-atm. coupling ##
        self.T_surf: Surface temperature (1000 bar) in K
        self.Rbulk_Rjup: Planet bulk radius in Jupiter radii
        self.Matm_earthunits: Atmospheric mass in Earth masses
        self.Mtot: Total mass (bulk + atm.) in Earth masses
        self.Rtot: Converged total planet radius in Jupiter radii
        ## From atmospheric model (only) ##
        self.r_atm_profile: Atm. radius profile in m
        self.g_atm_profile: Atm. gravity profile in m/s2
        self.P_atm_profile: Atm. pressure profile in Pa
        self.T_atm_profile: Atm. temperature profile in K
        self.rho_atm_profile: Atm. density profile in kg/m3
        """

        self.M_P = M_P
        self.x_core = x_core
        self.Zenv = Zenv

        x_H2O = 1.-self.x_core           # Envelope mass fraction = 1.0 - CMF
        P_surf = 1e3 * 1e5               # Surface P in Pa

        # Initial guess
        self.T_surf = Tguess             # Surface T in K

        # Input for heat model
        self.Teq = Teq
        self.Tint = Tint


        # Difference between Rguess and Rinterior
        diff = 1.0
        Rguess = 0.


        Tsurf_arr = np.zeros(20)
        R_arr = np.zeros(20)

        counter = 0

        # Convergence loop starts here
        while diff > tolerance:

            if counter == len(Tsurf_arr):
                print('Error in Coupling.py: The number of interior-atmosphere iterations is greater than 20.')
                print('The current relative difference between radii is',diff)
                print('Increase the tolerance.')
                sys.exit(1)

            Tsurf_arr[counter] = self.T_surf
            R_arr[counter] = Rguess
            """
            print("Input M = ",self.M_P)
            print("Input x_core = ",self.x_core)
            print("Input x_h2o = ",x_H2O)
            print("Input T_surf = ",self.T_surf)
            print("Input P_surf = ",P_surf)
            print("Input Zenv = ",self.Zenv)
            """


            self.myplanet.calc_radius(self.M_P,self.x_core,x_H2O,self.T_surf,P_surf,self.Zenv)

            """
            print('')
            print('Output:')
            print('R_P [R_E] = ', self.myplanet.R_P)
            print('rho [g.cm-3]= ', self.myplanet.rho_p)
            print('Envelope O:H = ', self.myplanet.OtoH)
            print('Envelope [Fe/H] = ', self.myplanet.metal)
            print('')
            """

            '''
            # Arrays for thermal evolution
            r_heat = self.myplanet.r[0 : self.myplanet.intrf[2]-1]
            cv_heat = self.myplanet.cv[0 : self.myplanet.intrf[2]-1]
            rho_heat = self.myplanet.rho[0 : self.myplanet.intrf[2]-1]
            '''


            # Atm model
            self.g_surf_planet = 100 * 9.8 * self.M_P/self.myplanet.R_P**2   # In cm/s2

            self.myatmmodel.calc_PTprofile(self.Tint,self.g_surf_planet,self.Zenv,self.Teq)
            """
            print('Atm. models from prt')
            print('Input:')
            print('T_int [K] = ', self.Tint_evol)
            print('g_surf [cgs] = ', g_surf_planet)

            print('')
            print('Output:')
            print('T_surf [K] = ', self.myatmmodel.Tsurf)
            print('P_surf [bar] = ', self.myatmmodel.Psurf)
            print('')
            """
            Rinterior = self.myplanet.R_P
            diff = abs(Rinterior - Rguess) / Rinterior
            Rguess = Rinterior
            counter = counter + 1
            self.T_surf = self.myatmmodel.Tsurf


        print("")
        print("Convergence reached in surface temperature and bulk radius")
        print("")

        Tsurf_arr[counter] = self.T_surf
        R_arr[counter] = Rguess
        """
        print("Rbulk [R_E] = ", R_arr)
        print("Tsurf [K] = ", Tsurf_arr)
        print("# iterations = ", counter)
        """

        # Final radius
        Rjup = 7.149e7    # Jupiter radius in m
        Rearth = 6.371e6  # Earth radius in m

        self.Rbulk_Rjup = R_arr[counter]*Rearth/Rjup


        Pbase = P_surf                                  # In Pa
        Rbulk = self.Rbulk_Rjup*Rjup                    # In m
        g_surf = 9.8 * self.M_P/R_arr[counter]**2       # In m/s2
        Matm = Pbase*4*math.pi*Rbulk**2/g_surf          # In kg
        self.Matm_earthunits = Matm/self.Mearth
        self.Mtot = self.M_P + self.Matm_earthunits

        self.myatmmodel.calc_thickness(self.Rbulk_Rjup,self.Matm_earthunits)


        self.Rtot = self.myatmmodel.total_radius

        self.r_atm_profile = self.myatmmodel.r
        self.g_atm_profile = self.myatmmodel.g_ode
        self.P_atm_profile = self.myatmmodel.P_ode
        self.T_atm_profile = self.myatmmodel.T_ode
        self.rho_atm_profile = self.myatmmodel.rho_ode



        """
        print('')
        print('Atm z calculation:')
        print('--------------------')

        print('')
        print('Input:')
        print('Rbulk [Rjup] = ', self.Rbulk_jup)
        print('')
        print('z_atm [R_jup] = ', self.myatmmodel.z_ode)
        print('Matm [Mearth] = ', Matm/self.Mearth)
        print('Matm [Mjup] = ', Matm/(self.Mearth*M_P))


        print('Output:')
        print('Rtotal [Rjup] = ', self.Rtot)
        """


