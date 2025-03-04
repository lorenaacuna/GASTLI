============
Installation
============

Prerequisites for basic installation
====================================

To install GASTLI, you need to:

- have a Python 3.9+ installation,
- have a fortran compiler, for example ``gfortran``,
- **download the zip file for the input data** `here <https://www.dropbox.com/scl/fi/p2kawqp8gtzh5psn21tjc/gastli_input_data.zip?rlkey=fc0mfxvpck5mukkqhk1f8hkad&st=ggsa4zmk&dl=0>`_, and place it in a directory in your computer that is easy to find. See the subsection "Setting the input data path" below for instructions on how to specify the path of this directory.

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

If instead you want to create a virtual environment, use the following command. This will create a new directory that contains all the files necessary to run the environment:

.. code-block:: bash

    python3 -m venv GASTLI_env

Then activate it with:

.. code-block:: bash

    source GASTLI_env/bin/activate


Once the environment is activated, we will be ready to install GASTLI's Python package. This can be done with pip:

.. code-block:: bash

    pip install gastli


Or for a particular version of the package:

.. code-block:: bash

    pip install gastli==0.9.3


The package can also be installed from source. This can be done by downloading the .zip file from the `GASTLI github repository <https://github.com/lorenaacuna/GASTLI>`_ while making sure you download the intended branch - main, or development - or using git to clone the repository in your local computer as explained `here <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`_:

.. code-block:: bash

    git clone https://github.com/lorenaacuna/GASTLI


To move to the root directory and install the package, run in your terminal: 

.. code-block:: bash

    cd GASTLI
    pip install .


Setting the input data path
---------------------------

After downloading the input data file (see "Prerequisites for basic installation"), place the ``gastli_input_data`` folder in any directory on your computer - just make sure you remember where, such as your Desktop. 

Next, set the environment variable ``GASTLI_input_data_path`` to your ``.bash_profile``, ``.bashrc``, or ``.zshrc`` file, depending on your operating system and shell type. This variable should contain the path to the ``gastli_input_data`` folder on your computer. Here are examples of how to do this:

- **On macOS**: Run the following command in your terminal: 

.. code-block:: bash

    echo 'export GASTLI_input_data_path="/PATH/TO/YOUR/FOLDER/gastli_input_data/"' >>~/.bash_profile

- **On Linux**: Run the following command in your terminal: 

.. code-block:: bash

    echo 'export GASTLI_input_data_path="/PATH/TO/YOUR/FOLDER/gastli_input_data/"' >>~/.bashrc

- **On a SLURM system**: Add the following line to your ``.sh`` script: 

.. code-block:: bash

    export GASTLI_input_data_path="/PATH/TO/YOUR/FOLDER/gastli_input_data/"

You only need to add this line to your bash profile once.

Alternatively, you can define the environment variable directly within your Python scripts. If you choose this method, add the following line to every script where you use GASTLI's functions or classes: 

.. code-block:: python

    os.environ['GASTLI_input_data_path']="/PATH/TO/YOUR/FOLDER/gastli_input_data/"


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


