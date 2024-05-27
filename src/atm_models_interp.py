

import numpy as np
import pandas as pd
import math
from scipy import interpolate
from scipy.integrate import odeint
import h5py
from scipy.interpolate import RegularGridInterpolator
import sys
import os
import gastli.Guillot10 as Guillot10
import gastli.radius_hse as rh
import gastli.constants as cte

#import matplotlib.pyplot as plt

# Auxiliary functions
# Functions for interpolation of tables

def maxloc(array_in,value_in):
    r'''This function locates the index of the element x <= value_in in an array
    i.e
    index_final = maxloc(array_in,value_in)
    is equivalent to Fortran's
    index_final = MAXLOC(array_in, DIM = 1, MASK=(array_in <= value_in))

    Args:
        array_in:
            array with values
        value_in:
            value whose index we want to locate
    '''
    array_minus = abs(array_in-value_in)
    array_index_prov = np.argmin(array_minus)
    if array_in[array_index_prov] <= value_in:
        index_final = array_index_prov
    else:
        index_final = array_index_prov-1
    return index_final


def interp2d_opt(x_in,y_in,x_array,y_array,z_array):
    r'''2d interpolation

    Args:
        x_in:
            Value to evaluate in x-axis
        y_in:
            Value to evaluate in y-axis
        x_array:
            Array of x-axis grid. This one is the one that always changes i.e [a,b,c,...,a,b,c...]
            Dimension is Nx.
        y_array:
            Array of y-axis grid. This one is the constant one i.e [d,d,d,...,e,e,e...]
            Dimension is Ny
        z_array:
            Array with the parameter to interpolate. Dimension is Nx*Ny

    Returns:
        z(x_in,y_in):
            final interpolated value
    '''

    Nx = len(x_array)

    i_y1 = maxloc(y_array,y_in)
    i_y2 = i_y1 + 1
    i_x1 = maxloc(x_array,x_in)
    i_x2 = i_x1 + 1

    y1 = y_array[i_y1]
    y2 = y_array[i_y2]
    x1 = x_array[i_x1]
    x2 = x_array[i_x2]

    f11 = z_array[i_y1*Nx+i_x1]
    f12 = z_array[i_y2*Nx+i_x1]
    f21 = z_array[i_y1*Nx+i_x2]
    f22 = z_array[i_y2*Nx+i_x2]

    conste = 1/((x2-x1)*(y2-y1))

    interp = conste*( f11*(x2-x_in)*(y2-y_in) + f21*(x_in-x1)*(y2-y_in) +\
                 f12*(x2-x_in)*(y_in-y1) + f22*(x_in-x1)*(y_in-y1) )

    return interp


def vapor_pressure(Tin):
    r'''This function calculates the saturation curve for water (Wagner & Pruss 2002)

    Args:
        Tin:
            Temperature in K

    Returns:
        Pout:
            Saturation pressure at Tin in MPa
    '''
    Tc = 647.096   # K
    pc = 22.064    # MPa
    rhoc = 322     # kg/m3

    a1 = -7.85951783
    a2 = 1.84408259
    a3 = -11.7866497
    a4 = 22.6807411
    a5 = -15.9618719
    a6 = 1.80122502

    theta = Tin/Tc
    tau = 1 - theta

    term1 = a1*tau
    term2 = a2*tau**1.5
    term3 = a3*tau**3
    term4 = a4*tau**3.5
    term5 = a5*tau**4
    term6 = a6*tau**7.5


    x = (Tc/Tin) * (term1+term2+term3+term4+term5+term6)
    pout = pc * np.exp(x)

    return pout

def vapor_density(Tin):
    r'''This function calculates the density of water vapour in the saturation limit (Wagner & Pruss 2002)

    Args:
        Tin:
            Temperature in K

    Returns:
        rhoout:
            Density in kg/m3
    '''

    Tc = 647.096   # K
    pc = 22.064    # MPa
    rhoc = 322     # kg/m3

    c1 = -2.03150240
    c2 = -2.68302940
    c3 = -5.38626492
    c4 = -17.2991605
    c5 = -44.7586581
    c6 = -63.9201063

    theta = Tin/Tc
    tau = 1 - theta

    term1 = c1*tau**(2/6)
    term2 = c2*tau**(4/6)
    term3 = c3*tau**(8/6)
    term4 = c4*tau**(18/6)
    term5 = c5*tau**(37/6)
    term6 = c6*tau**(71/6)

    x = term1+term2+term3+term4+term5+term6
    rhoout = rhoc * np.exp(x)

    return rhoout

# Functions for metallicity to Fe/H
# Constants

O_to_H_sun = 4.898e-4         # From Lodders+ 03
mu_H = 1.
mu_H2O = 18.

def O_to_H_molecular(Z):
    r'''This function converts metal mass fraction into metallicity (or O:H ratio). See Fortney et al. 2013

    Args:
        Z:
            metal mass fraction in atmosphere

    Returns:
        OtoH:
            O:H ratio (metallicity [Fe/H])

    '''
    X = (1-Z)/1.379
    mu_X = 2*mu_H
    mu_Z = mu_H2O
    OtoH = (Z/X)*(mu_X/mu_Z)/( 2 + 2*( (Z/X)*(mu_X/mu_Z) ) )

    return OtoH


class atm_models_interp:
    def __init__(self,path_to_file,name_grid="gastli_default_atm_grid.hdf5"):
        r"""Class defining object to carry out interpolation of atmospheric models

        Args:
            path_to_file:
                path to input data directory
            name_grid (optional):
                name of grid with atmospheric models. The name of the default grid is
                "gastli_default_atm_grid.hdf5"
                The atmospheric grid must have the following dimensions:
                "C/O": C/O ratio
                "FeH": log10(Fe/H) in x solar units
                "logg": log10 of surface gravity in cm/s2
                "Teq": Equilibrium temperature in K
                "Tint": Internal temperature in K
                "pressure": Atmospheric pressure in bar
                There should be two 6-dimensional datasets in the hdf5 file: one named "PT_profiles" and the
                other "metal mass fractions".
                The data must be defined on a rectilinear grid, in other words, a regular grid structure
                If some models are missing, fill these with np.nan. The atm_models_interp class will raise a warning
                and a recommended value for a np.nan model
        """

        '''
        Description:
        - Initialises parameters
        - Loads in data from atmospheric models and water and H/He EOS
        '''
        # Planet parameters
        self.T_int_pl = 0.0         # K
        self.g_surf_pl = 0.0        # cgs
        self.Zenv_pl = 0.0          # Mass fraction
        self.Rbulk_pl = 0.0         # Rjupiter units

        # Dimensions
        self.npoints = 130

        # Load atmosphere data
        self.path_to_file = path_to_file

        file_atm = h5py.File(self.path_to_file+"Input/Atmospheric data/"+name_grid, 'r')

        self.data_set_atm = file_atm['PT_profiles'][()]
        self.data_set_metal_mass_fractions = file_atm['metal_mass_fractions'][()]

        self.CO_atm = file_atm['CO'][()]
        self.FeH_atm = file_atm['FeH'][()]
        self.Teq_atm = file_atm['Teq'][()]
        self.Tint_atm = file_atm['Tint'][()]
        self.logg_atm = file_atm['logg'][()]
        press_atm_file = file_atm['pressure'][()]

        # Create atm. data function: PT profile
        self.temperature_func = RegularGridInterpolator((self.CO_atm,self.FeH_atm,self.logg_atm,self.Teq_atm,self.Tint_atm,\
                                                        press_atm_file), self.data_set_atm, bounds_error=False,\
                                                        fill_value=None)


        # Create atm. data function: abundances
        self.metal_mass_fraction_func = RegularGridInterpolator((self.CO_atm,self.FeH_atm,self.logg_atm,self.Teq_atm,\
                                                        self.Tint_atm, press_atm_file), \
                                                        self.data_set_metal_mass_fractions, bounds_error=False,\
                                                        fill_value=None)


        # Load EOS data

        #Read Chabrier data
        Nt = 121
        Np = 441

        self.logTinput = np.zeros(Nt)
        self.logPinput = np.zeros(Np)

        for i in range(0, Nt):
            self.logTinput[i] = 2.0 + i * 0.05  # Decimal logarithm of temperature [K]

        for i in range(0, Np):
            self.logPinput[i] = -9.0 + i * 0.05  # Decimal logarithm of pressure [GPa]

        ntot = Nt * Np

        self.logrho = np.zeros(ntot)

        c = 0

        file = open(self.path_to_file+'Input/Chabrier/TABLEEOS_2021_TP_Y0275_v1', 'r')
        Lines = file.readlines()

        for j in range(1, Nt + 1):
            for i in range(1, Np + 1):
                index_arr = (j - 1) * Np + j + 1
                string = Lines[index_arr + i - 1]
                string_arr = map(float, string.split())
                float_arr = [float(i_n) for i_n in string_arr]
                self.logrho[c] = float_arr[2]  # Decimal logarithm of density [g/cm3]
                c = c + 1


        data = pd.read_csv(self.path_to_file+"Tables/aqua_eos_pt_v1_0.dat",\
                           sep='\s+', header=None, skiprows=19)
        self.rho_aqua = np.asarray(data[2])

        data = pd.read_csv(self.path_to_file+"Tables/P_AQUA_Pa.dat", sep='\s+', header=None)
        self.press_aqua = np.asarray(data[0])  # Pressure in Pa

        data = pd.read_csv(self.path_to_file+"Tables/T_AQUA_K.dat", sep='\s+', header=None)
        self.temp_aqua = np.asarray(data[0])  # Temperature in K

        data = pd.read_csv(self.path_to_file+"Input/Chabrier/HG23_Vmix_Smix.csv", sep=',', header=0)
        self.V_mix = np.asarray(data['Vmix'])

        data = pd.read_csv(self.path_to_file+"Input/Chabrier/logP_HG23_corr.dat", sep='\s+', header=None)
        self.logP_Vmix = np.asarray(data[0])

        data = pd.read_csv(self.path_to_file+"Input/Chabrier/logT_HG23_corr.dat", sep='\s+', header=None)
        self.logT_Vmix = np.asarray(data[0])


        #############################
        #############################
        # END Reading in files
        #############################
        #############################

    def calc_interior_mass_fraction(self,Tint,g_surf,Teq,CO,log_FeH,P_surf=1e3):
        r'''Function that calculates the metal mass fraction for a given log10(metallicity) [x solar].
            It uses data from easychem.

        Args:
            Tint:
                Internal temperature in K
            g_surf:
                Surface gravity in cm/s2
            Teq:
                Equilibrium temperature in K (at Bond albedo zero)
            CO:
                atm. C-to-O ratio
            log_FeH:
                log10(metallicity), where the atm. metallicity is in x solar units
            P_surf:
                surface pressure in bar

        Returns:
            Pprofile:
                pressure profile in bar
            MMF_profile:
                metal mass fraction profile for atm. thickness calculation
            MMF_surf:
                metal mass fraction in the surface (1000 bars) for interior model input
        '''

        self.T_int_pl = Tint
        self.g_surf_pl = g_surf
        self.Teq_pl = Teq
        self.CO = CO
        self.log_FeH = log_FeH
        self.P_surf = P_surf

        logg_pl = np.log10(self.g_surf_pl)

        # Interpolation
        # Order: FeH_atm, logg_atm, Teq_atm, Tint_atm, press_atm
        self.press_atm = np.logspace(-5, np.log10(P_surf), self.npoints)

        pts = np.zeros((self.npoints, 6))
        pts[:, 0] = self.CO
        pts[:, 1] = self.log_FeH
        pts[:, 2] = logg_pl
        pts[:, 3] = self.Teq_pl
        pts[:, 4] = self.T_int_pl
        pts[:, 5] = self.press_atm

        # Extrapolation warnings
        if self.CO < min(self.CO_atm) or self.CO > max(self.CO_atm):
            print("C/O is out of atmospheric grid limits. Extrapolating")

        if self.log_FeH < min(self.FeH_atm) or self.log_FeH > max(self.FeH_atm):
            print("Fe/H is out of atmospheric grid limits. Extrapolating")

        if self.Teq_pl < min(self.Teq_atm) or self.Teq_pl > max(self.Teq_atm):
            print("Equilibrium temperature is out of atmospheric grid limits. Extrapolating")

        if logg_pl < min(self.logg_atm) or logg_pl > max(self.logg_atm):
            print("Surface gravity is out of atmospheric grid limits. Extrapolating")


        self.Pprofile = self.press_atm
        self.MMF_profile = self.metal_mass_fraction_func(pts)
        self.MMF_surf = self.MMF_profile[-1]


        if np.isnan(self.MMF_surf):
            Teq_node = maxloc(self.Teq_atm, self.Teq_pl)
            logg_node = maxloc(self.logg_atm, logg_pl)
            FeH_node = maxloc(self.FeH_atm, self.log_FeH)
            CO_node = maxloc(self.CO_atm, self.CO)

            index_Teq = np.array([Teq_node,Teq_node])
            index_logg = np.array([logg_node,logg_node])
            index_FeH = np.array([FeH_node,FeH_node])
            index_CO = np.array([0, 1])

            maxTint_list = []


            for i_Teq in index_Teq:
                for i_logg in index_logg:
                    for i_FeH in index_FeH:
                        for i_CO in index_CO:
                            Tsurf_fornans = self.data_set_atm[i_CO, i_FeH, i_logg, i_Teq, :, -1]
                            mask = ~np.isnan(Tsurf_fornans)
                            Tint_masked = self.Tint_atm[mask]
                            maxTint_list.append(max(Tint_masked))

            maxTint = np.array(maxTint_list)
            Tint_limit = min(maxTint)

            print("No atmospheric models available for this case (np.nan in grid).")
            print("Decrease the interior temperature or increase the surface pressure")
            #print("Decrease the interior temperature below ", Tint_limit, "K")
            sys.exit(1)



    def calc_PTprofile(self,Tint,g_surf,Teq,Zenv=0.03,FeH_flag=True,CO_def=0.55,guillot=False,P_surf=1e3):
        r'''Calculation of pressure-temperature atmospheric profile by interpolating the grid of data of
            atmospheric models

        Args:
            Tint:
                Internal temperature in K
            g_surf:
                Surface gravity in cm/s2
            Teq:
                Equilibrium temperature in K (at Bond albedo zero)
            CO_def (optional):
                C-to-O ratio. Default value is 0.55 (solar)
            FeH_flag (optional):
                equals True if the amount of metals is specified as log10(metallicity) in x solar
                units. If specified as metal mass fraction, set FeH_flag = False. In the former case, the metal mass
                fraction profile is extracted from easychem data (calculated with calc_interior_mass_fraction).
                In the latter, it is set constant to Zenv, and converted to log10(metallicity) with Fortney+13 relation
                (appendix A1) to interpolate the PT profiles from atmospheric grid
            Zenv (optional):
                atmospheric metal mass fraction (for FeH_flag = False). Default value is 0.03.

        Returns:
            Pprofile:
                atm. pressure profile in bar
            Tprofile:
                atm. temperature profile in K
            Tsurf:
                Temperature at bottom of atmosphere in K
            Psurf:
                Pressure at bottom of atmosphere in bar. Normally set to 1000 bars
        '''

        self.T_int_pl = Tint
        self.g_surf_pl = g_surf
        self.Teq_pl = Teq

        logg_pl = np.log10(self.g_surf_pl)

        if FeH_flag == False:
            """
            if FeH_flag == False, you are providing the metal mass fraction 
            and log(Fe/H) needs to be calculated
            """
            # Calculate from Fortney+ 2013 approximation
            self.Zenv_pl = Zenv
            self.log_FeH = np.log10(O_to_H_molecular(self.Zenv_pl) / O_to_H_sun)
            # Set CO
            self.CO = CO_def
            # set metal mass fraction profile
            self.MMF_profile = np.ones(self.npoints)*self.Zenv_pl
            self.P_surf = P_surf
        else:
            """
            if FeH_flag == True, you are providing the log(Fe/H)
            and Zenv was calculated by calc_interior_mass_fraction
            """
            self.Zenv_pl = self.MMF_surf


        # Interpolation
        # Order: CO_atm, FeH_atm, logg_atm, Teq_atm, Tint_atm, press_atm
        self.press_atm = np.logspace(-5, np.log10(P_surf), self.npoints)

        pts = np.zeros((self.npoints, 6))
        pts[:, 0] = self.CO
        pts[:, 1] = self.log_FeH
        pts[:, 2] = logg_pl
        pts[:, 3] = self.Teq_pl
        pts[:, 4] = self.T_int_pl
        pts[:, 5] = self.press_atm
        """
        print("log_FeH = ",log_FeH_pl)
        print("log_g = ", logg_pl)
        print("Teq = ", self.Teq_pl)
        print("Tint = ", self.T_int_pl)
        """


        # Extrapolation warnings
        if self.CO < min(self.CO_atm) or self.CO > max(self.CO_atm):
            print("C/O is out of atmospheric grid limits. Extrapolating")

        if self.log_FeH < min(self.FeH_atm) or self.log_FeH > max(self.FeH_atm):
            print("Fe/H is out of atmospheric grid limits. Extrapolating")

        if self.Teq_pl < min(self.Teq_atm) or self.Teq_pl > max(self.Teq_atm):
            print("Equilibrium temperature is out of atmospheric grid limits. Extrapolating")

        if logg_pl < min(self.logg_atm) or logg_pl > max(self.logg_atm):
            print("Surface gravity is out of atmospheric grid limits. Extrapolating")


        '''
        Class' final output from atmospheric models
        '''

        if guillot==True:
            kappa_IR = 0.01
            gamma = 0.4

            gravity = 1e1 ** logg_pl

            Guillot10_model = Guillot10.guillot_global(self.press_atm, kappa_IR, gamma, gravity,\
                                                        self.T_int_pl, self.Teq_pl)

            self.Pprofile = self.press_atm
            self.Tprofile = Guillot10_model
            self.Tsurf = self.Tprofile[-1]
            self.Psurf = self.Pprofile[-1]

        else:
            self.Pprofile = self.press_atm
            self.Tprofile = self.temperature_func(pts)
            self.Tsurf = self.Tprofile[-1]
            self.Psurf = self.Pprofile[-1]

            if np.isnan(self.Tsurf):
                Teq_node = maxloc(self.Teq_atm, self.Teq_pl)
                logg_node = maxloc(self.logg_atm, logg_pl)
                FeH_node = maxloc(self.FeH_atm, self.log_FeH)
                CO_node = maxloc(self.CO_atm, self.CO)

                index_Teq = np.array([Teq_node,Teq_node])
                index_logg = np.array([logg_node,logg_node])
                index_FeH = np.array([FeH_node,FeH_node])
                index_CO = np.array([0, 1])

                maxTint_list = []

                #self.temperature_func = RegularGridInterpolator(
                #    (self.CO_atm, self.FeH_atm, self.logg_atm, self.Teq_atm, self.Tint_atm, \
                #     press_atm)


                for i_Teq in index_Teq:
                    for i_logg in index_logg:
                        for i_FeH in index_FeH:
                            for i_CO in index_CO:
                                """
                                print("i_CO =",i_CO)
                                print("i_FeH =", i_FeH)
                                print("i_logg =", i_logg)
                                print("i_Teq =", i_Teq)
                                """
                                Tsurf_fornans = self.data_set_atm[i_CO, i_FeH, i_logg, i_Teq, :, -1]
                                mask = ~np.isnan(Tsurf_fornans)
                                Tint_masked = self.Tint_atm[mask]
                                maxTint_list.append(max(Tint_masked))

                maxTint = np.array(maxTint_list)
                Tint_limit = min(maxTint)

                print("No atmospheric models available for this case (np.nan in grid).")
                print("Decrease the interior temperature or decrease the surface pressure")
                #print("Decrease the interior temperature below ", Tint_limit, "K")
                sys.exit(1)



    def calc_thickness(self,Rbulk,Matm_earthunits):
        r'''Calculates thickness of atmosphere

        Args:
            Rbulk:
                Planet bulk radius (centre to Psurf) in Jupiter radii
            Matm_earthunits:
                Atmospheric mass in Earth units

        Returns:
            r:
                Radius atm. profile in m
            P_ode:
                Pressure atm. profile in Pa
            g_ode:
                Gravity atm. profile in m/s2
            rho_ode:
                Density atm. profile in kg/m3
            T_ode:
                Temperature atm. profile in K
            total_radius:
                Total planet radius in Jupiter radii
            z_ode:
                Atmospheric thickness in Jupiter radii
        '''

        Rjup = 7.149e7                         # Jupiter radius in m
        Mjup = 1.899e27                        # Jupiter mass in kg
        Rearth = 6.371e6                       # Earth radius in m
        Mearth = 5.972e24                      # Earth mass in kg
        G = 6.674e-11                          # Grav. cte in m3/kg/s2 (SI)

        self.Rbulk_pl = Rbulk
        g0 = self.g_surf_pl                    # cm/s2
        Rbulk_inm = self.Rbulk_pl*Rjup         # In m

        Mbulk =  ( ((g0/100)/9.8) * (Rbulk_inm/Rearth)**2 )  # In Mearth units

        Z = self.MMF_surf
        mu = Z * 18. + (1 - Z) * 2.33



        MMWs_amu = mu * np.ones_like(self.Tprofile)
        pressures_cgs = self.Pprofile * 1e5 * 10.
        reference_pressure_cgs = 1000. * 1e5 * 10.
        reference_radius_cm = self.Rbulk_pl * 11.2 * cte.constants.r_e * 100.
        """
        print(self.Tprofile)
        print(MMWs_amu)
        print(g0)
        print(reference_pressure_cgs)
        print(reference_radius_cm)
        print(pressures_cgs)
        """
        radius_planet = rh.calc_radius_hydrostatic_equilibrium(self.Tprofile,
                                                               MMWs_amu,
                                                               g0,
                                                               reference_pressure_cgs,
                                                               reference_radius_cm,
                                                               pressures_cgs,
                                                               variable_gravity=True)

        self.r = radius_planet/100
        self.P_ode = self.Pprofile * 1e5
        self.g_ode = G * (Mbulk * Mearth) / self.r ** 2
        self.T_ode = self.Tprofile

        func_radius = interpolate.interp1d(self.Pprofile, self.r)
        self.total_radius = func_radius(20 * 1e-3)/Rjup
        # modified for Solar System planets
        #self.total_radius = func_radius(1.) / Rjup
        self.z_ode = self.total_radius - Rbulk_inm/Rjup


        #
        # Calculation of density profile
        #

        Y_til = 0.245
        Y_astk = 0.275

        density_profile = np.zeros(self.npoints)
        density_water = np.zeros(self.npoints)
        density_hhe = np.zeros(self.npoints)

        for i in range(0, self.npoints):
            P_forprof = self.Pprofile[i] * 1e5  # Pa
            T_forprof = self.Tprofile[i]        # K

            # rho from aqua
            # here: if P>Psat -> dens = dens_sat(T)

            if T_forprof < 647.096:
                if P_forprof >= 1e6 * vapor_pressure(T_forprof):
                    rho_aqua_forprof = vapor_density(T_forprof)
                else:
                    rho_aqua_forprof = interp2d_opt(T_forprof, P_forprof, self.temp_aqua, self.press_aqua,
                                                    self.rho_aqua)
            else:
                rho_aqua_forprof = interp2d_opt(T_forprof, P_forprof, self.temp_aqua, self.press_aqua, \
                                                self.rho_aqua)  # kg/m3

            density_water[i] = rho_aqua_forprof

            # rho from H/He EOS
            logP_forprof = math.log(P_forprof / 1e9, 10)
            logT_forprof = math.log(T_forprof, 10)
            logrho_cgs = interp2d_opt(logP_forprof, logT_forprof, self.logPinput, self.logTinput, self.logrho)
            rho_cd21 = (10 ** logrho_cgs)

            # HG23 Correction
            Vmix_value = interp2d_opt(logT_forprof, logP_forprof, self.logT_Vmix, self.logP_Vmix + 10.0, self.V_mix)
            inv_rho_corr = (1 / rho_cd21) + Vmix_value * (Y_astk * (1 - Y_astk) - Y_til * (1 - Y_til))  # In g/cm3
            rho_corr = 1 / inv_rho_corr
            rho_si = rho_corr * 1e3  # In kg/m3
            density_hhe[i] = rho_si

            # Note that H/He is enclosing the condensates from easychem
            # but we consider the contribution of condensates negligible
            # The thickness may increase slightly due to a lower density
            # compared to the condensate case
            inv_rho = self.MMF_profile[i] / rho_aqua_forprof + (1 - self.MMF_profile[i]) / rho_si

            density_profile[i] = 1 / inv_rho

        self.rho_ode = density_profile

        """
        Y_til = 0.245
        Y_astk = 0.275
        #
        # Calculation of density profile
        #
        density_profile = np.zeros(self.npoints)
        density_water = np.zeros(self.npoints)
        density_hhe = np.zeros(self.npoints)

        for i in range(0, self.npoints):
            P_forprof = self.Pprofile[i] * 1e5    # Pa
            T_forprof = self.Tprofile[i]          # K

            # rho from aqua
            # here: if P>Psat -> dens = dens_sat(T)
            
            if T_forprof<647.096:
                if P_forprof >= 1e6*vapor_pressure(T_forprof):
                    rho_aqua_forprof = vapor_density(T_forprof)
                else:
                    rho_aqua_forprof = interp2d_opt(T_forprof, P_forprof, self.temp_aqua, self.press_aqua,
                                                    self.rho_aqua)
            else:
                rho_aqua_forprof = interp2d_opt(T_forprof, P_forprof, self.temp_aqua, self.press_aqua,\
                                                self.rho_aqua)  # kg/m3

            density_water[i] = rho_aqua_forprof

            # rho from H/He EOS
            logP_forprof = math.log(P_forprof / 1e9, 10)
            logT_forprof = math.log(T_forprof, 10)
            logrho_cgs = interp2d_opt(logP_forprof, logT_forprof, self.logPinput, self.logTinput, self.logrho)
            rho_cd21 = (10 ** logrho_cgs)


            # HG23 Correction
            Vmix_value = interp2d_opt(logT_forprof, logP_forprof, self.logT_Vmix, self.logP_Vmix+10.0, self.V_mix)
            inv_rho_corr = (1 / rho_cd21) + Vmix_value * (Y_astk * (1 - Y_astk) - Y_til * (1 - Y_til))    # In g/cm3
            rho_corr = 1/inv_rho_corr
            rho_si = rho_corr * 1e3             # In kg/m3
            density_hhe[i] = rho_si

            # Note that H/He is enclosing the condensates from easychem
            # but we consider the contribution of condensates negligible
            # The thickness may increase slightly due to a lower density
            # compared to the condensate case
            inv_rho = self.MMF_profile[i] / rho_aqua_forprof + (1 - self.MMF_profile[i]) / rho_si

            density_profile[i] = 1 / inv_rho

        '''
        # Output profiles
        nplot = len(density_profile)
        points = np.linspace(1, nplot, num=nplot)

        data = np.zeros((nplot, 6))
        data[:, 0] = points
        data[:, 1] = self.Pprofile * 1e5
        data[:, 2] = self.Tprofile
        data[:, 3] = density_profile
        data[:, 4] = density_hhe
        data[:, 5] = density_water

      
        # Ouput file
        fmt = '%d', '%1.4e', '%1.4e', '%1.4e', '%1.4e', '%1.4e '
        np.savetxt('Output/profile_density_atm.dat', data, delimiter='      ',
                    header='  i  P[Pa]  T[K]  rho_mix[kg.m3]  rho_HHe[kg.m3]  rho_H2O[kg.m-3]  ', comments='', fmt=fmt)
        '''

        #
        # Prepare interpolation of P-T and P-rho
        #

        Rjup = 7.149e7                         # Jupiter radius in m
        Mjup = 1.899e27                        # Jupiter mass in kg
        Rearth = 6.371e6                       # Earth radius in m
        Mearth = 5.972e24                      # Earth mass in kg
        G = 6.674e-11                          # Grav. cte in m3/kg/s2 (SI)

        g0 = self.g_surf_pl/100                # m/s2
        Rbulk_inm = self.Rbulk_pl*Rjup         # In m

        Mbulk =  ( (g0/9.8) * (Rbulk_inm/Rearth)**2 )  # In Mearth units

        # Pressure in bars!
        func_rho = interpolate.interp1d(self.Pprofile, density_profile, bounds_error=False,fill_value="extrapolate")
        func_temp = interpolate.interp1d(self.Pprofile, self.Tprofile, bounds_error=False,fill_value="extrapolate")

        #limit_level = 10.     # 10 bars

        def thickness_func(y, t):
            '''
            Function that calculates pressure and radius derivatives with respect to the enclosed mass
            Args:
            :param y: [p, r] values to evaluate
            :param t: enclosed mass
            press units: in Pa
            mass units: in kg
            radius units: in m
            :return: dydt = [dP/dr,dr/dm]
            '''
            press_ft, radius_ft = y

            # dr/dm
            press_inbar = press_ft/1e5                   # From Pa to bar
            #limit_level = 10.                            # 2 bar
            limit_dens = func_rho(limit_level)
            #limit_dens = 1e-2
            #print(limit_dens)
            rho_ft = max(func_rho(press_inbar),limit_dens)
            radius_term = 1 / (4 * math.pi * rho_ft * radius_ft ** 2.)

            # dP/dm
            press_term = - G * (Mbulk*Mearth+t)/(4*math.pi*radius_ft**4.)


            dydt = [press_term, radius_term]
            return dydt


        Psurf = 1e3 * 1e5                     # In Pa
        y0 = [Psurf, Rbulk_inm]
        fmin = 1e-3                           # Matm in kg
        fmax = 1.3*Matm_earthunits*Mearth

        t = np.linspace(fmin, fmax, 100000)    # mass points in kg

        sol, infodict = odeint(thickness_func, y0, t, full_output=True)


        pressure_sol = sol[:, 0]              # Pa
        radius_sol = sol[:, 1]                # m



        #print('pressure_sol = ',pressure_sol/1e5)
        #print('radius_sol = ', radius_sol/Rjup)

        # Calculate density
        rho_ode_nofilter = func_rho(pressure_sol/1e5)

        # Mask for density>limit
        #mask_dens = (rho_ode_nofilter>1e-2)
        mask_dens = (pressure_sol > limit_level*1e5)

        # Mask pressure and radius accordingly
        pressure_sol_masked = pressure_sol[mask_dens]
        radius_sol_masked = radius_sol[mask_dens]
        t_masked = t[mask_dens]


        # Calculate extrapolation functions
        func_radius = interpolate.interp1d(pressure_sol_masked, radius_sol_masked, \
                                             bounds_error=False, fill_value="extrapolate")

        func_t = interpolate.interp1d(pressure_sol_masked, t_masked, \
                                             bounds_error=False, fill_value="extrapolate")


        self.P_ode = np.logspace(-5, 3, 1000)*1e5

        # Calculate gravity
        self.r = func_radius(self.P_ode)
        grav = G * (func_t(self.P_ode) + Mbulk * Mearth) / func_radius(self.P_ode) ** 2

        # Arrays
        self.g_ode = grav
        self.T_ode = func_temp(self.P_ode/1e5)
        self.rho_ode = func_rho(self.P_ode/1e5)

        # check
        dp_dr = np.diff(self.P_ode)/np.diff(self.r)
        if any(dp_dr>-1e-1):
            print('Error in atmospheric thickness calculation: your pressure level (limit_level) is too low')
            print('Your current limit_level =',limit_level," bar")
            print('Increase your value for limit_level')
            sys.exit(1)

        ##############################################################

        self.total_radius = func_radius(20*1e-3*1e5)/Rjup
        self.z_ode = self.total_radius - Rbulk_inm/Rjup
        """


        ####
        """
        mask_press = (pressure_sol > 0)

        pressure_sol_masked = pressure_sol[mask_press]
        radius_sol_masked = radius_sol[mask_press]
        t_masked = t[mask_press]

        grav = G * (t_masked + Mbulk * Mearth) / radius_sol_masked ** 2


        # Arrays
        self.r = radius_sol_masked
        self.P_ode = pressure_sol_masked
        self.g_ode = grav
        self.rho_ode = func_rho(self.P_ode/1e5)
        self.T_ode = func_temp(self.P_ode/1e5)


        ##############################################################

        func_pressure = interpolate.interp1d(pressure_sol/1e5, radius_sol/Rjup, \
                                             bounds_error=False,fill_value="extrapolate")

        self.total_radius = func_pressure(20*1e-3)
        self.z_ode = self.total_radius - Rbulk_inm/Rjup

        """

