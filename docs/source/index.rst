=========================================
GASTLI: The GAS gianT modeL for Interiors
=========================================

Welcome to the **GASTLI** documentation. GASTLI is a Python package for calculating 
interior structure models for gas giant exoplanets. GASTLI allows you to compute mass-radius 
relations, thermal evolution curves and interior structure profiles for a wide variety of irradiations,
masses and internal compositions (check out our :doc:`atm_grid` section).


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
   :caption: User Guide

   installation
   atm_grid
   test_suite

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorial
   retrieval

Release Notes
=============
.. include::  ../../CHANGELOG.rst


