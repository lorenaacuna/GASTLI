
from gastli.GASTLI import int_planet
from gastli.atm_models_interp import atm_models_interp

import numpy as np
import time
import math
import sys

class coupling:
    def __init__(self,path_to_file,name_grid=None,pow_law_formass=0.32,j_max=30):
        r"""Class defining objects to run one interior-atmosphere coupled model

        Args:
            path_to_file:
                path to input data directory
            pow_law_formass:
                power exponent to estimate the initial guess of the planet radius in the interior model.
                Default is 0.32. Increase if planet is very massive (greater than 5 Jupiter masses aprox). Decrease if
                core mass fraction is very low (< 0.03 approx.) and/or planet is low mass (15-20 Earth masses approx.)
        """

        '''
        - Initialises interior and atm. models
        - Other constants
        '''

        self.pow_law_formass = pow_law_formass
        self.path_to_file = path_to_file
        self.j_max = j_max

        # Initialise interior model
        self.myplanet = int_planet(path_to_file=self.path_to_file, pow_law=self.pow_law_formass, j_max=self.j_max)

        self.myplanet.setup_parameters()

        # Initialise atm models
        if name_grid == None:
            self.myatmmodel = atm_models_interp(path_to_file=self.path_to_file)
        else:
            self.myatmmodel = atm_models_interp(path_to_file=self.path_to_file, name_grid=name_grid)

        # Important constants
        #self.M_P = 318.             # M in M_earth
        self.Mearth = 5.972e24       # In m


    def main(self,M_P,x_core,Teq,Tint,CO=0.55,log_FeH=0.,Zenv=0.03,FeH_flag=True,Tguess=2000.,Rguess=11.2,\
             tolerance=1e-3,guillot=False,P_surf=1e3):
        r"""Function that runs coupled interior-atmosphere model
            ## Output parameters of atmosphere class ##
            To access the output parameters of the atmosphere class, such as the atmospheric thickness or profiles,
            add "myatmmodel." to the name of the parameter.
            Example: to get the pressure and temperature atm. profiles, get the attributes from the coupling class
            named "myatmmodel.P_ode" and "myatmmodel.T_ode".
            To see the list of parameters from the atmosphere class, see the docstring of function "calc_thickness"
            ## Output parameters of interior class ##
            To access the output parameters of the interior class, such as the interior profiles,
            add "myplanet." to the name of the parameter.
            Example: to get the gravity and radius interior profiles, get the attributes from the coupling class
            named "myplanet.g" and "myplanet.r".
            To see the list of parameters from the interior class, see the docstring of function "calc_radius"


        Args:
            M_P:
                Planet mass in Earth masses
            x_core:
                Core mass fraction
            Teq:
                Equilibrium temperature in K (at zero Bond albedo)
            Tint:
                Internal temperature in K
            CO (optional):
                C-to-O ratio. Default value is 0.55 (solar)
            FeH_flag (optional):
                equals True if the amount of metals is specified as log10(metallicity) in x solar
                units. If specified as metal mass fraction, set FeH_flag = False. In the former case, the metal mass
                fraction profile is extracted from easychem data (calculated with calc_interior_mass_fraction).
                In the latter, it is set constant to Zenv, and converted to log10(metallicity) with Fortney+13 relation
                (appendix A1) to interpolate the PT profiles from atmospheric grid
            Zenv (optional):
                atmospheric metal mass fraction (for FeH_flag = False). Default value is 0.03.
            log_FeH:
                log10(metallicity) in x solar units (for FeH_flag = True).
                Default value is zero (solar composition)
            Tguess (optional):
                Initial guess for the surface temperature in K. Default is 2000 K.
            Rguess (optional):
                Initial guess for planet radius in Earth radii. Default is 11.2 Earth radii (Jupiter's
                radius). Changing it may speed up calculations when FeH_flag = True.
            tolerance (optional):
                maximum relative difference in radius between interior-atm. steps. Default is 0.001
            guillot (optional):
                False if you do not want to use Guillot 2010 atm. profile
            P_surf:
                Boundary pressure between interior and atmosphere. Default is 1000 bars. For models with high Tint
                you may need to decrease it to 9.5 bars (if using the default atm. grid)

        Return:
            Rtot:
                Converged total planet radius in Jupiter radii
            Mtot:
                Total mass (bulk interior + atm.) in Earth masses
            T_surf:
                Final surface temperature (1000 bar) in K
            CMF_conv:
                Re-computed core mass fraction. It takes into account the mass of the atmosphere.
            Rbulk_Rjup:
                Planet bulk radius in Jupiter radii
            Matm_earthunits:
                Atmospheric mass in Earth masses

        """

        self.M_P = M_P
        self.x_core = x_core
        self.Zenv = Zenv
        self.log_FeH = log_FeH
        self.CO_pl = CO
        self.FeH_flag = FeH_flag
        self.Zenv = Zenv


        x_H2O = 1.-self.x_core           # Envelope mass fraction = 1.0 - CMF

        # Initial guess
        self.T_surf = Tguess             # Surface T in K

        # Temperatures
        self.Teq = Teq
        self.Tint = Tint


        # Difference between Rguess and Rinterior
        n_iterations = 20
        diff = 1.0
        Tsurf_arr = np.zeros(n_iterations+1)
        R_arr = np.zeros(n_iterations+1)
        counter = 0

        # Convergence loop starts here
        # Run for at least two iterations to avoid "false accurate" runs
        while diff > tolerance or counter<2:
            #print(diff)

            if counter == n_iterations:
                print('Warning in Coupling.py: The number of interior-atmosphere iterations is greater than ',\
                      n_iterations)
                print('The current relative difference between radii is',diff)
                #print('Increase the tolerance.')
                # For debugging
                """
                print("Rbulk array [R_E] = ", R_arr)
                print("Tsurf array [K] = ", Tsurf_arr)
                print("# iterations = ", counter)
                """
                # Re-initialise everything again
                ## Interior class
                self.pow_law_formass = self.pow_law_formass - 0.005
                print('Readjusting mass power law to ', self.pow_law_formass)
                #print(self.path_to_file)
                #print(self.pow_law_formass)
                self.myplanet = int_planet(path_to_file=self.path_to_file, pow_law=self.pow_law_formass)
                self.myplanet.setup_parameters()

                ## Arrays
                diff = 1.0
                Tsurf_arr = np.zeros(n_iterations+1)
                R_arr = np.zeros(n_iterations+1)
                counter = 0

                continue
                #sys.exit(1)

            #print("Tsurf_arr = ", Tsurf_arr)
            #print("counter = ", counter)

            Tsurf_arr[counter] = self.T_surf
            R_arr[counter] = Rguess

            #print("Tsurf = ",Tsurf_arr)
            #print("Rarr  = ", R_arr)
            #print('diff =', diff)
            #print("tolerance = ",tolerance)

            # Calculate Zenv for interior
            # if FeH_flag = True, then Zenv for interior module is calculated from easychem tables with log(Fe/H)
            # if FeH_flag is False, then Zenv for interior is the one provided by user
            if FeH_flag == True:
                self.g_surf_planet = 100 * 9.8 * self.M_P / Rguess ** 2  # In cm/s2
                self.myatmmodel.calc_interior_mass_fraction(self.Tint, self.g_surf_planet, self.Teq, self.CO_pl,\
                                                            self.log_FeH,P_surf=P_surf)
                self.Zenv = self.myatmmodel.MMF_surf



            self.myplanet.calc_radius(self.M_P,self.x_core,x_H2O,self.T_surf,P_surf*1e5,self.Zenv)

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


            # Update surface gravity
            self.g_surf_planet = 100 * 9.8 * self.M_P/self.myplanet.R_P**2   # In cm/s2


            # Atm model
            if FeH_flag == True:
                self.myatmmodel.calc_PTprofile(self.Tint,self.g_surf_planet,self.Teq,guillot=guillot,P_surf=P_surf)
            else:
                self.myatmmodel.calc_PTprofile(self.Tint, self.g_surf_planet, self.Teq, self.Zenv, FeH_flag=False,\
                                               CO_def=self.CO_pl,guillot=guillot,P_surf=P_surf)

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

        #print("Tsurf_arr = ", Tsurf_arr)
        #print("counter = ", counter)

        Tsurf_arr[counter] = self.T_surf
        R_arr[counter] = Rguess
        """
        print("Rbulk array [R_E] = ", R_arr)
        print("Tsurf array [K] = ", Tsurf_arr)
        print("# iterations = ", counter)
        """
        # Final radius
        Rjup = 7.149e7    # Jupiter radius in m
        Rearth = 6.371e6  # Earth radius in m

        self.Rbulk_Rjup = R_arr[counter]*Rearth/Rjup

        Pbase = P_surf*1e5                              # In Pa
        Rbulk = self.Rbulk_Rjup*Rjup                    # In m
        g_surf = 9.8 * self.M_P/R_arr[counter]**2       # In m/s2
        Matm = Pbase*4*math.pi*Rbulk**2/g_surf          # In kg
        self.Matm_earthunits = Matm/self.Mearth
        self.Mtot = self.M_P + self.Matm_earthunits
        self.CMF_conv = (M_P*x_core)/self.Mtot


        self.myatmmodel.calc_thickness(self.Rbulk_Rjup,self.Matm_earthunits)


        self.Rtot = self.myatmmodel.total_radius
        """
        self.r_atm_profile = self.myatmmodel.r
        self.g_atm_profile = self.myatmmodel.g_ode
        self.P_atm_profile = self.myatmmodel.P_ode
        self.T_atm_profile = self.myatmmodel.T_ode
        self.rho_atm_profile = self.myatmmodel.rho_ode
        """


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

