=========================================
GASTLI: The GAS gianT modeL for Interiors
=========================================

Welcome to the **GASTLI** documentation! GASTLI is a user-friendly Python package developed to support both observers and expert modelers in modeling of exoplanet interiors with volatile-rich envelopes. 

Key features:
=============

- **Mass-radius and thermal evolution curves**: GASTLI enables you to compute mass-radius relations, thermal evolution curves and interior structure profiles using the ``Coupling`` and ``Thermal_evolution`` classes. 
- **Up-to-date equations of state (EOS) tables**: GASTLI integrates the latest non-ideal EOS for H/He, water and rock to calculate the density and Grüneisen parameter (adiabatic gradient), which is essential for accurate modelling of high-pressure planetary envelopes. Refer to `this paper <https://www.aanda.org/articles/aa/full_html/2024/08/aa50559-24/aa50559-24.html>`_ for available EOS. If you wish to use other EOS in your models, please check our :doc:`contributing` section.
- **Flexible mass, composition and irradiation input**: GASTLI supports a wide range of core mass fractions, envelope metal content (i.e. atmospheric metallicity), masses, and irradiation conditions - provided the atmospheric grid covers relevant surface gravities, atmospheric metallicity and equilibrium temperatures. The default atmospheric grid's applicability is covered in our :doc:`atm_grid` section.
- **Customizable atmospheric grid**: The default atmospheric grid corresponds to clear, 1D atmospheres in chemical equilibrium. If you wish to change these assumptions for your models (e.g., clouds, 3D effects, chemical disequilibrium) or our default atmospheric grid does not meet your required surface gravities, compositions and equilibrium temperatures, GASTLI allows you to input a custom atmospheric grid. For details on grid file format and requirements, check our :doc:`atm_grid` section.
- **Retrievals**: GASTLI is optimized to compute numerous models within a grid, allowing interpolation for retrievals on mass, radius, atmospheric metallicity, age, and/or internal temperature. See our :doc:`retrieval` section for a comprehensive guide on performing retrievals with GASTLI.


License and how to cite
=======================

GASTLI is open source under a BSD 3-Clause License. 

The following papers document GASTLI:

- The base package for gas giants is described in `Acuña et al. (2024) <https://arxiv.org/abs/2406.10032>`_.
- Interior-atmosphere coupling algorithm: `Acuña et al. (2021) <https://www.aanda.org/articles/aa/full_html/2021/03/aa39885-20/aa39885-20.html>`_.

Please cite these two papers if you make use of GASTLI.

.. _contact: acuna@mpia.de


Development
===========

GASTLI is based on the Marseille Super-Earth Interior Model (MSEI), which was adapted and optimized by Lorena Acuña, supported by Laura Kreidberg, Paul Mollière, Meng Zhai, Iva Momcheva and Morgan Fouenesneau (Max Planck Institute for Astronomy). MSEI's former developers include Olivier Mousis, Magali Deleuil, Artem Aguichine, and Bastien Brugger (Laboratoire d'Astrophysique de Marseille).


Contents
========

.. toctree::
   :maxdepth: 1
   :caption: Installation

   installation

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   tutorial
   atm_grid
   retrieval

.. toctree::
   :maxdepth: 1
   :caption: Code Documentation

   test_suite
   contributing
   
   


Release Notes
=============
.. include::  ../../CHANGELOG.rst


