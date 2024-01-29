# GASTLI
GAS gianT modeL for Interiors

Download zip file and create an environment for GASTLI: 
```
conda create -n GASTLI_packaging_test python=3.10
```

Activate environment:
```
conda activate GASTLI_packaging_test
```

Install pip:
```
conda install pip
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
