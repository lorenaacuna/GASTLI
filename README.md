# GASTLI
GAS gianT modeL for Interiors

Download zip file for the input data [here](https://www.dropbox.com/scl/fi/fagjvx4umni04f4cxqdhy/gastli_input_data.zip?rlkey=8y47oldcyokbxqo5z2e5izw9c&dl=0)

Create an environment for GASTLI: 
```
conda create -n for_GASTLI_env python=3.10
```

Activate environment:
```
conda activate for_GASTLI_env
```

Install GASTLI:
```
pip install .
```

You can check quickly if it installed correctly by running this in Python:
```
import gastli.dimensions as dim
print(dim.dimensions.n_lay)
import gastli.constants as cte
print(cte.constants.f_alloy_e)
```

Or by running the scripts in the test/ or Examples/ folders. The test/test_for_interior_module.py is particularly useful to know how fast your computer is to run the interior model.

We don't have the documentation yet (stay tuned), but the Examples/ scripts and the docstrings in the src/ code are your guide for the moment.

