

# Import interior module
from gastli.GASTLI import int_planet

# Other Python modules
import numpy as np
import time

path_to_input = "/Users/acuna/Desktop/gastli_input_data/"
myplanet_test = int_planet(path_to_file=path_to_input)

print("Number of layers: ",myplanet_test.n_lay)
print("Temperature array with water data: ",myplanet_test.Ttp)
print("ilayer value: ",myplanet_test.ilayer)

myplanet_test.setup_parameters()

start = time.time()
myplanet_test.calc_radius(318.,0.10,1-0.10,1300.,1e3*1e5,0.04)
end = time.time()
print('Output radius [R_earth] = ', myplanet_test.R_P)
print('Calculation time [s] = ', end-start)
