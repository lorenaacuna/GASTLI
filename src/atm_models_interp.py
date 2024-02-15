

import numpy as np
import pandas as pd
import math
from scipy import interpolate
from scipy.integrate import odeint
import h5py
from scipy.interpolate import RegularGridInterpolator
import sys
import os


# Auxiliary functions
# Functions for interpolation of tables

def maxloc(array_in,value_in):
    '''
    This function locates the index of the element x <= value_in in an array
    i.e
    index_final = maxloc(array_in,value_in)
    is equivalent to Fortran's
    index_final = MAXLOC(array_in, DIM = 1, MASK=(array_in <= value_in))
    Args:
        array_in: array with values
        value_in: value whose index we want to locate
    '''
    array_minus = abs(array_in-value_in)
    array_index_prov = np.argmin(array_minus)
    if array_in[array_index_prov] <= value_in:
        index_final = array_index_prov
    else:
        index_final = array_index_prov-1
    return index_final


def interp2d_opt(x_in,y_in,x_array,y_array,z_array):
    '''
    2d interpolation
    :param x_in: Value to evaluate in x-axis
    :param y_in: Value to evaluate in y-axis
    :param x_array: Array of x-axis grid. This one is the one that always changes i.e [a,b,c,...,a,b,c...]
                    Dimension is Nx.
    :param y_array: Array of y-axis grid. This one is the constant one i.e [d,d,d,...,e,e,e...]
                    Dimension is Ny
    :param z_array: Array with the parameter to interpolate. Dimension is Nx*Ny
    :return: final interpolated value, z(x_in,y_in)
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
    '''
    This function calculates the saturation curve for water (Wagner & Pruss 2002)
    :param Tin: Temperature in K
    :return: Saturation pressure at Tin in MPa
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
    '''
    This function calculates the density of water vapour in the saturation limit (Wagner & Pruss 2002)
    :param Tin: Temperature in K
    :return: Density in kg/m3
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
    '''
    This function converts metal mass fraction into metallicity (or O:H ratio). See Fortney et al. 2013
    :param Z: metal mass fraction in atmosphere
    :return: O:H ratio (metallicity [Fe/H])
    '''
    X = (1-Z)/1.379
    mu_X = 2*mu_H
    mu_Z = mu_H2O
    OtoH = (Z/X)*(mu_X/mu_Z)/( 2 + 2*( (Z/X)*(mu_X/mu_Z) ) )

    return OtoH


class atm_models_interp:
    """ Class defining objects for carrying out interpolation of atmospheric models

    Args:
        path_grid (optional): path to atmospheric model grid. It must be in hdf5 format,
        and must have the following dimensions:
        "FeH": log10(Fe/H) in x solar units
        "logg": log10 of surface gravity in cm/s2
        "Teq": Equilibrium temperature in K
        "Tint": Internal temperature in K
        "pressure": Atmospheric pressure in bar
        The 5-dimensional dataset must be named "PT_profiles" in the hdf5 file
        The data must be defined on a rectilinear grid, in other words, a regular grid structure
        If some models are missing, fill these with np.nan. The atm_models_interp class will raise a warning
        and a recommended value for a np.nan model
    """

    def __init__(self,path_to_file):
        '''
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

        file_atm = h5py.File(self.path_to_file+"Input/Atmospheric data/CO_055_grid.hdf5", 'r')

        self.data_set_atm = file_atm['PT_profiles'][()]
        self.FeH_atm = file_atm['FeH'][()]
        self.Teq_atm = file_atm['Teq'][()]
        self.Tint_atm = file_atm['Tint'][()]
        self.logg_atm = file_atm['logg'][()]
        press_atm = file_atm['pressure'][()]

        # Create atm. data function
        self.temperature_func = RegularGridInterpolator((self.FeH_atm,self.logg_atm,self.Teq_atm,self.Tint_atm,\
                                                        press_atm), self.data_set_atm, bounds_error=False,\
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




    def calc_PTprofile(self,Tint,g_surf,Zatm,Teq):
        '''
        Calculation of pressure-temperature atmospheric profile by interpolating the grid of data of atmospheric models
        Args:
        :param Tint: Internal temperature in K
        :param g_surf: Surface gravity in cm/s2
        :param Zatm: Atmosphere metal mass fraction
        :param Teq: Equilibrium temperature in K (at Bond albedo zero)
        Return:
        self.Pprofile: atm. pressure profile in bar
        self.Tprofile: atm. temperature profile in K
        self.Tsurf: Temperature at bottom of atmosphere in K
        self.Psurf: Pressure at bottom of atmosphere in bar. Normally set to 1000 bars
        '''

        self.T_int_pl = Tint
        self.g_surf_pl = g_surf
        self.Zenv_pl = Zatm
        self.Teq_pl = Teq

        log_FeH_pl = np.log10(O_to_H_molecular(self.Zenv_pl)/O_to_H_sun)
        logg_pl = np.log10(self.g_surf_pl)

        # Interpolation
        # Order: FeH_atm, logg_atm, Teq_atm, Tint_atm, press_atm
        self.press_atm = np.logspace(-5, 3, self.npoints)

        pts = np.zeros((self.npoints, 5))
        pts[:, 0] = log_FeH_pl
        pts[:, 1] = logg_pl
        pts[:, 2] = self.Teq_pl
        pts[:, 3] = self.T_int_pl
        pts[:, 4] = self.press_atm
        """
        print("log_FeH = ",log_FeH_pl)
        print("log_g = ", logg_pl)
        print("Teq = ", self.Teq_pl)
        print("Tint = ", self.T_int_pl)
        """


        # Extrapolation warnings
        if log_FeH_pl < min(self.FeH_atm) or log_FeH_pl > max(self.FeH_atm):
            print("Fe/H is out of atmospheric grid limits. Extrapolating")

        if self.Teq_pl < min(self.Teq_atm) or self.Teq_pl > max(self.Teq_atm):
            print("Equilibrium temperature is out of atmospheric grid limits. Extrapolating")

        if logg_pl < min(self.logg_atm) or logg_pl > max(self.logg_atm):
            print("Surface gravity is out of atmospheric grid limits. Extrapolating")


        '''
        Class' final output from atmospheric models
        '''
        self.Pprofile = self.press_atm
        self.Tprofile = self.temperature_func(pts)
        self.Tsurf = self.Tprofile[-1]
        self.Psurf = self.Pprofile[-1]

        if np.isnan(self.Tsurf):
            Teq_node = maxloc(self.Teq_atm, self.Teq_pl)
            logg_node = maxloc(self.logg_atm, logg_pl)
            FeH_node = maxloc(self.FeH_atm, log_FeH_pl)

            index_Teq = np.array([Teq_node,Teq_node+1])
            index_logg = np.array([logg_node,logg_node+1])
            index_FeH = np.array([FeH_node,FeH_node+1])

            maxTint_list = []

            for i_Teq in index_Teq:
                for i_logg in index_logg:
                    for i_FeH in index_FeH:
                        Tsurf_fornans = self.data_set_atm[i_FeH, i_logg, i_Teq, :, -1]
                        mask = ~np.isnan(Tsurf_fornans)
                        Tint_masked = self.Tint_atm[mask]
                        maxTint_list.append(max(Tint_masked))

            maxTint = np.array(maxTint_list)
            Tint_limit = min(maxTint)

            print("No atmospheric models available for this case (np.nan in grid).")
            print("Decrease the interior temperature below ", Tint_limit, "K")
            sys.exit(1)



    def calc_thickness(self,Rbulk,Matm_earthunits):
        '''
        Calculates thickness of atmosphere
        Args:
        :param Rbulk: Planet bulk radius (centre to Psurf) in Jupiter radii
        :param Matm_earthunits: Atmospheric mass in Earth units
        :return:
        self.r: Radius atm. profile in m
        self.P_ode: Pressure atm. profile in Pa
        self.g_ode: Gravity atm. profile in m/s2
        self.rho_ode: Density atm. profile in kg/m3
        self.T_ode: Temperature atm. profile in K
        self.total_radius: Total planet radius in Jupiter radii
        self.z_ode: Atmospheric thickness in Jupiter radii
        '''


        self.Rbulk_pl = Rbulk

        Y_astk = 0.275
        Y_til = 0.245

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
            #"""
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

            inv_rho = self.Zenv_pl / rho_aqua_forprof + (1 - self.Zenv_pl) / rho_si

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
            press_inbar = press_ft/1e5                  # From Pa to bar
            limit_level = 2.                            # 2 bar
            limit_dens = func_rho(limit_level)
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

