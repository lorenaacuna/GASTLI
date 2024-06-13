============
Installation
============

Prerequisites for basic installation
====================================
To install GASTLI, you need to install:

- Python 3.9+,
- a fortran compiler, for example ``gfortran``,
- download the zip file for the input data `here <https://www.dropbox.com/scl/fi/pvctvgps31oxs0agtyb5b/gastli_input_data.zip?rlkey=2svm3y9usqgqgp3ab6lv5vxsa&st=yyltghlr&dl=0>`_, and place it in a directory in your computer that is easy to find

We recommend installing GASTLI via the pip install command

Installation of GASTLI via pip install
=============================================
To install GASTLI via pip install, open a terminal and create a new conda environment (for example, "GASTLI_env") by running:

.. code-block:: bash

    conda create -n GASTLI_env python=3.10

Activate the environment:

.. code-block:: bash

    conda activate GASTLI_env

Download the GASTLI github repository, and from the root directory, run: 

.. code-block:: bash

    pip install .

You can check quickly if GASTLI installed correctly by running this in Python:

.. code-block:: python

   import gastli.dimensions as dim
   print(dim.dimensions.n_lay)
   import gastli.constants as cte
   print(cte.constants.f_alloy_e)

The output should be:

.. code-block:: python

   3
   0.13
