# GASTLI
GAS gianT modeL for Interiors

Download zip file for the input data [here](https://www.dropbox.com/scl/fi/euhdgfhyotxtw0jqe4ydi/gastli_input_data.zip?rlkey=vw7epxyly844xfkf9yws4fai5&dl=0))

Create an environment for GASTLI: 
```
conda create -n GASTLI_packaging_test python=3.10
```

Activate environment:
```
conda activate GASTLI_packaging_test
```

Install GASTLI:
```
pip install .
```

If installed correctly, this should work in python:
```
import gastli.dimensions as dim
print(dim.dimensions.n_lay)
import gastli.constants as cte
print(cte.constants.f_alloy_e)
```
