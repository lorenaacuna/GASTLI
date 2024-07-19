==========
Test suite
==========

To run the test suite, ``pytest`` needs to be installed. Pytest version 8.2.2 is already installed with GASTLI, as it is part of the requirements. In the main repository, the test files can be found in the directory ``test``. From this directory, running the command ``pytest`` will run all test files whose names start with ``test_``. There is one test per class: interior only (``test_GASTLI``), atmosphere only (``test_atm_models_interp``), interior-atmosphere coupling (``test_Coupling``) and thermal evolution (``test_Thermal_evolution``).

