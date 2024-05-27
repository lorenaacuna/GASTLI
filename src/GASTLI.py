
import gastli.constants as cte
import gastli.dimensions as dim
from gastli.fortinput import fort_input as fi
from gastli.gastli_interior import gastli_interior as interior

import numpy as np
#import time


class int_planet:
    def __init__(self, path_to_file, j_max=30, cnt_conv_max=3, conv_prec=1e-5, \
                 pow_law=0.32, chk_EOS=0, EOS_lim_P=[5e11, 5e11, 5e11, 5e11, \
                                                     5e11, 5e11, 5e11, 5e11, 5e11, 5e11], corEOS=1):
        r"""Class defining objects for carrying out interior structure calculations
        for a given set of mass, composition and surface conditions

        Args:
            path_to_file:
                path to input data directory
            j_max (Optional):
                Maximum number of iterations - must be < 100.
            cnt_conv_max (Optional):
                Maximum number of detected oscillations.
            conv_prec (Optional):
                Precision of the convergence condition.
            pow_law (Optional):
                Power exponent for planet radius estimation. Default is 0.32. Increase if planet is very massive
            chk_EOS (Optional):
                Index to check the validity range of the used EOS. Do not change
            EOS_lim_P (Optional, array):
                Upper limit in pressure. Do not change
            corEOS (Optional):
                Type of correction on the thermodynamical parameters of Vinet EOS.
                (=0: no correction / =1: range [1:1.5] / =2: range [1:5] / =3:
                range [1:10]). Do not change
        """


        # Arguments of __init__
        self.path_to_file = path_to_file
        self.j_max = j_max
        self.cnt_conv_max = cnt_conv_max
        self.conv_prec = conv_prec
        self.pow_law = pow_law
        self.chk_EOS = chk_EOS
        self.EOS_lim_P = EOS_lim_P
        self.corEOS = corEOS

        # Dimensions
        self.n_lay = dim.dimensions.n_lay
        self.n_EOS = dim.dimensions.n_eos
        self.n_pts = dim.dimensions.n_pts
        self.n_int = dim.dimensions.n_int
        self.n_mat_max = dim.dimensions.n_mat_max

        # Load constant file
        self.Ttp,self.Ptp = fi.read_constants(self.path_to_file)

        # Load layer files
        self.ilayer,self.use_lay,self.iEOS,self.EOS_th,self.del_T,self.n_mat,\
        self.x,self.cf_Mmol,self.cf_Z_eff,self.cf_rho_0,self.cf_T_0,\
        self.cf_K_0,self.cf_Kp_0,self.cf_gam_0,self.cf_q,self.cf_a_T,\
        self.cf_b_T,self.cf_c_T,self.cf_a_P,\
        self.cf_T_d0 = fi.read_layer_parameters(self.path_to_file,self.n_lay,self.n_mat_max)

        # Load EOS data

        # Rock EOS from SESAME data base
        self.logT_sesame,self.logP_sesame,self.logrho_sesame,self.logS_sesame,\
        self.dlrho_dlT_p_sesame,self.dlS_dlT_p_sesame,\
        self.dlrho_dlP_t_sesame = fi.read_eos_sesame(self.path_to_file)

        # Water EOS from Mazevet+19 (only the density, for optimization)
        self.P_maz_water,self.T_maz_water,\
        self.rho_maz_water = fi.read_mazevet_water(self.path_to_file)

        # H/He EOS from Chabrier+23
        self.logT_input,self.logP_input,self.logrho_input,\
        self.logrho_ch,self.logU_ch,self.logP_ch,self.logS_ch,\
        self.dlrho_dlT_P,self.dlrho_dlP_T,self.dlS_dlT_P,self.grad_ad,\
        self.grad_ad_PT = fi.read_eos_ch(self.path_to_file)


        # H/He EOS correction from Howard & Guillot 2023
        self.logP_HG,self.logT_HG,\
        self.Vmix_HG, self.Smix_HG = fi.read_hg23_corr(self.path_to_file)


        # Initialise parameters
        self.f_alloy = 0.0
        self.MgD = 0.0
        self.MgSi = 0.0

        self.M_P = 0.0
        self.x_core = 0.0
        self.x_H2O = 0.0
        self.T_surf = 0.0
        self.P_surf = 0.0

        self.x_H2Oc = 0.0
        self.x_corec = 0.0
        self.M_Pc = 0.0
        self.R_P = 0.0
        self.FeSi = 0.0
        self.rho_p = 0.0
        self.intrf = np.array(np.zeros(self.n_lay + 1))
        self.intrf_hist = np.zeros((self.n_lay + 1,100))
        self.iter_num = np.zeros(100)
        self.g = np.array(np.zeros(self.n_pts))
        self.T = np.array(np.zeros(self.n_pts))
        self.P = np.array(np.zeros(self.n_pts))
        self.rho = np.array(np.zeros(self.n_pts))
        self.cv = np.array(np.zeros(self.n_pts))
        self.r = np.array(np.zeros(self.n_pts))


        #############################
        #############################
        # END Reading in files
        #############################
        #############################

    def setup_parameters(self, f_alloy=cte.constants.f_alloy_e, \
                         MgD=cte.constants.mgd_e, \
                         MgSi=cte.constants.mgsi_e):
        r"""Set up parameters for chemistry of mantle and core
            Do not change any of the parameters below

        Args:
            f_alloy (Optional):
                if no Fe in core, not applicable
            MgD (Optional):
                if not Vinet Eq. for mantle, not applicable
            MgSi (Optional):
                if not Vinet Eq. for mantle, not applicable

        Returns:
            x:
                matrix with material properties
        """

        self.f_alloy = f_alloy
        self.MgD = MgD
        self.MgSi = MgSi

        self.Mmol, self.Z_eff, self.rho_0, self.T_0, self.K_0, self.Kp_0, \
        self.gam_0, self.q, self.a_T, self.b_T, self.c_T, self.a_P, self.T_d0, xout = \
        fi.calc_parameters(self.n_mat, self.x, self.cf_Mmol, self.cf_Z_eff, \
        self.cf_rho_0, self.cf_T_0, self.cf_K_0, self.cf_Kp_0, self.cf_gam_0, \
        self.cf_q, self.cf_a_T, self.cf_b_T, self.cf_c_T, self.cf_a_P, \
        self.cf_T_d0, self.f_alloy, self.MgD, self.MgSi, [self.n_lay, self.n_mat_max])

        self.x = xout

    def calc_radius(self, M_P, x_core, x_H2O, T_surf, P_surf, Zenv):
        r'''Function that calculates planet radius and interior structure

        Args:
            M_P:
            Planet mass in Earth masses
            x_core:
            Core mass fraction
            x_H2O:
            Envelope mass fraction. Defined as 1-x_core in this case
            T_surf:
            Outer surface (boundary) temperature in K
            P_surf:
            Outer surface (boundary) pressure in Pa

        Returns:
            x_H2Oc:
                Re-computed envelope mass fraction (for checking)
            x_corec:
                Re-computed core mass fraction (for checking)
            M_Pc:
                Re-computed interior mass in Earth masses (for checking)
            R_P:
                Planet radius in Earth radii
            FeSi:
                Fe-to-Si ratio (not applicable if no Fe core)
            rho_p:
                Planet density in g/cm3
            OtoH:
                Envelope O:H ratio
            metal:
                Envelope metallicity ([Fe/H] in x solar units)
            intrf:
                Array with indexes of interfaces between layers
            g:
                Gravity acceleration profile array, g(r), in cm/s2
            T:
                Temperature profile array, T(r), in K
            P:
                Pressure profile array, P(r), in Pa
            rho:
                Density profile array, rho(r), in kg/m3
            cv:
                Specific heat capacity at constant volume profile array, Cv(r), in J/kg/K
            entropy:
                Entropy profile array, S(r), in J/kg/K
            r:
                Radius profile array, r, in m
        '''
        self.M_P = M_P
        self.x_core = x_core
        self.x_H2O = x_H2O
        self.T_surf = T_surf
        self.P_surf = P_surf
        self.Zenv = Zenv
        """
        print("Input M = ", self.M_P)
        print("Input x_core = ", self.x_core)
        print("Input x_h2o = ", self.x_H2O)
        print("Input T_surf = ", self.T_surf)
        print("Input P_surf = ", self.P_surf)
        print("Input Zenv = ", self.Zenv)
        """
        print("")
        print("Running interior structure model")

        # Output parameters
        m_pout, x_coreout, x_h2oout, r_pout, fesiout, rhoout, rout, interfaceout, intrf_hist, iter_num,\
        gout, tout, pout, rhoarrout, cvout, Sout,\
        OtoHout, metalout = interior.gastli_interior_subroutine(self.j_max,\
        # Input: run variables & material data
        self.cnt_conv_max,self.conv_prec,self.pow_law,self.chk_EOS, \
        self.EOS_lim_P, self.corEOS, self.Ttp, self.Ptp, self.f_alloy,self.MgD, \
        self.MgSi, self.ilayer, self.use_lay, self.iEOS, self.EOS_th,self.n_mat, \
        self.del_T, self.x, self.Mmol, self.Z_eff, self.rho_0, self.T_0,self.K_0, \
        self.Kp_0, self.gam_0, self.q, self.a_T, self.b_T, self.c_T,self.a_P,self.T_d0, \
        # Input: EOS data
        # SESAME for rock
        self.logT_sesame,self.logP_sesame,self.logrho_sesame,self.logS_sesame,\
        self.dlrho_dlT_p_sesame,self.dlS_dlT_p_sesame,self.dlrho_dlP_t_sesame,\
        # Mazevet+19 for water
        self.P_maz_water,self.T_maz_water,self.rho_maz_water,\
        # Chabrier+ for H/He
        self.logT_input,self.logP_input,self.logrho_input,self.logrho_ch,\
        self.logU_ch,self.logP_ch,self.logS_ch,self.dlrho_dlT_P,self.dlrho_dlP_T,\
        self.dlS_dlT_P,self.grad_ad,self.grad_ad_PT,\
        # HG23 correction for H/He
        self.logP_HG,self.logT_HG,self.Vmix_HG,self.Smix_HG,\
        # Input: planet parameters
        self.M_P, self.x_core, self.x_H2O, self.T_surf,self.P_surf,\
        # NEW PLANET PARAMETERS INPUT
        self.Zenv)

        '''
        Class' final output from interior model
        '''
        self.x_H2Oc = x_h2oout
        self.x_corec = x_coreout
        self.M_Pc = m_pout
        self.R_P = r_pout
        self.FeSi = fesiout
        self.rho_p = rhoout
        self.OtoH = OtoHout
        self.metal = metalout
        self.intrf = interfaceout
        self.intrf_hist = intrf_hist
        self.iter_num = iter_num
        self.g = gout
        self.T = tout
        self.P = pout
        self.rho = rhoarrout
        self.cv = cvout
        self.entropy = Sout
        self.r = rout


