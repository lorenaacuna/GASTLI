

# Import interior module
from gastli.GASTLI import int_planet

# Other Python modules
import numpy as np
import time

myplanet_test = int_planet()

n_lay_expected = int(3)
Ttp_expected = np.array([273.16,251.165,256.164,273.31,355.,647.096])
ilayer_expected = np.array([1,2,3])

myplanet_test.setup_parameters()

M_P = 318.
x_core = 0.10
T_surf = 1300.
P_surf = 1e3*1e5
Zenv = 0.04

start = time.time()
myplanet_test.calc_radius(M_P,x_core,1-x_core,T_surf,P_surf,Zenv)
end = time.time()
print('Calculation time [s] = ', end-start)

expected_radius = 10.188813052994131

def test_answer():
    assert myplanet_test.n_lay == n_lay_expected
    assert len(myplanet_test.Ttp) == len(Ttp_expected)
    assert all([a == b for a, b in zip(myplanet_test.Ttp, Ttp_expected)])
    assert len(myplanet_test.ilayer) == len(ilayer_expected)
    assert all([a == b for a, b in zip(myplanet_test.ilayer, ilayer_expected)])
    assert myplanet_test.R_P == expected_radius
