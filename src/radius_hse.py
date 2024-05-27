
import numpy as np
from gastli.fortspec import fort_spec as fs

kB = 1.3806488e-16
amu = 1.66053892e-24

def calc_radius_hydrostatic_equilibrium(temperatures_K,
                                        MMWs_amu,
                                        reference_gravity_cgs,
                                        reference_pressure_cgs,
                                        reference_radius_cm,
                                        pressures_cgs,
                                        variable_gravity=True):

    rho = pressures_cgs * MMWs_amu * amu / kB / temperatures_K
    radius = fs.calc_radius(pressures_cgs,
                            reference_gravity_cgs,
                            rho,
                            reference_pressure_cgs,
                            reference_radius_cm,
                            variable_gravity)
    return radius
