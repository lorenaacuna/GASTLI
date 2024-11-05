Fortran installation in Mac OS (Intel)
-----------------------------------

For Mac OS Intel, we recommend to use `homebrew <https://brew.sh>`_. Trying to install gfortran (or gcc) from source can be complicated, and setup-related errors can be easily introduced. 

For a safe installation, start executing in your terminal:

.. code-block:: bash

    brew update
    brew upgrade
    brew doctor

The command ``brew doctor`` may show a list of fixes. We recommend to go through this list before moving on to the next step.

Finally, the fortran compiler can be installed by running in the terminal: 

.. code-block:: bash

    brew install gcc

