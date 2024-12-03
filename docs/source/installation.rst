============
Installation
============

Prerequisites for basic installation
====================================

To install GASTLI, you need to:

- have a Python 3.9+ installation,
- have a fortran compiler, for example ``gfortran``,
- **download the zip file for the input data** `here <https://www.dropbox.com/scl/fi/p2kawqp8gtzh5psn21tjc/gastli_input_data.zip?rlkey=fc0mfxvpck5mukkqhk1f8hkad&st=ggsa4zmk&dl=0>`_, and place it in a directory in your computer that is easy to find

If you do not have Fortran installed, follow the Fortran installation instructions appropriate for your OS. We offer an installation guide for different OS below. Once Fortran is installed, you can then proceed with the GASTLI package installation. 

Fortran installation 
=====================

.. warning::

   Mac OS with Apple Silicon processors (M1, M2 and M3) have ARM64 as the native architecture. These systems have shown problems when installing a Python package that compiles a fortran backend because Python's default fortran compiler uses Intel architecture. One way to recognise this is to try installing GASTLI with pip (or from source) and testing the installation as indicated below. Users that get a message such as the following need to follow the instructions for Apple Silicon in the :doc:`macos_arm64_install` section:
   
   .. code-block:: python

    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    ImportError: dlopen(/Users/acuna/anaconda3/envs/GASTLI_norosetta/lib/python3.10/site-packages/gastli/dimensions.cpython-310-darwin.so, 0x0002): symbol not found in flat namespace (_f2pyinitdimensions_)


.. toctree::
   linux_install
   macos_intel_install
   macos_arm64_install


Installation of GASTLI's Python package
=======================================


Our recommendation is to use a Python virtual environment such as `venv <https://docs.python.org/3/library/venv.html>`_ or `conda <https://docs.anaconda.com/anaconda/install/index.html>`_, to prevent potential conflicts. 

If you decide to create a new conda environment (for example, "GASTLI_env"), run in the terminal:

.. code-block:: bash

    conda create -n GASTLI_env python=3.10

Then activate the environment:

.. code-block:: bash

    conda activate GASTLI_env

After activating the environment, we will be ready to install GASTLI with pip, ``pip install gastli==0.9.1``

.. warning::

   The latest stable version is 0.9.1, which corresponds to the main branch. Version 0.9.2 is the development branch. We recommend users to make sure that the install version 0.9.1 if they use pip ``pip install gastli==0.9.1``



Installation from source
------------------------

You can download the .zip file from the `GASTLI github repository <https://github.com/lorenaacuna/GASTLI>`_, or use git to clone the repository in your local computer as explained `here <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`_. Make sure it is the main branch (and not the development one). Then from the root directory, run in your terminal: 

.. code-block:: bash

    pip install .


Testing the installation
------------------------

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


