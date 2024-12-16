

# Import coupling module
import gastli.Coupling as cpl

# Other Python modules
import numpy as np

# Create coupling class
my_coupling = cpl.coupling(j_max=99,pow_law_formass=0.315)


# Input for interior
# 1 Mjup in Mearth units
Mjup = 318.
# M = 1 Mjup
M_P = Mjup
# Internal and equilibrium temperatures
Tintpl = 107.
Teqpl = 105.

# Core mass fraction and envelope metal mass fraction
CMF = 0.01



# Call to coupled interior-atm. model (and time it)
# Case 1, log(Fe/H) is known
my_coupling.main(M_P, CMF, Teqpl, Tintpl, CO=0.55,log_FeH=0.)


Zenv_expected = 0.0135335763819442
Mtot_expected = 318.0390773794476
T_surf_expected = 1328.9730303240333
Rbulk_Rjup_expected = 0.9738243869527924
Rtot_expected = 0.9799956272275188

def test_answer():
    assert my_coupling.myatmmodel.Zenv_pl == Zenv_expected
    assert my_coupling.Mtot == Mtot_expected
    assert my_coupling.T_surf == T_surf_expected
    assert my_coupling.Rbulk_Rjup == Rbulk_Rjup_expected
    assert my_coupling.Rtot == Rtot_expected


