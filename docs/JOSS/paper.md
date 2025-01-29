---
title: 'GASTLI: A Python package for coupled interior–atmosphere modelling of volatile-rich planets'
tags:
  - Python
  - astronomy
  - exoplanets
  - interiors
  - atmospheres
authors:
  - name: Lorena Acuña 
    corresponding: true # (This is how to denote the corresponding author)
    orcid: 0000-0002-9147-7925
    affiliation: "1" 
  - name: Laura Kreidberg
    affiliation: "1"
  - name: Meng Zhai
    orcid: 0000-0003-1207-3787
    affiliation: "1,2"
  - name: Paul Mollière
    affiliation: "1"
  - name: Morgan Fouesneau
    orcid: 0000-0001-9256-5516
    affiliation: "1"
affiliations:
 - name: Max Planck Institute for Astronomy, Königstuhl 17, 69117 Heidelberg, Germany
   index: 1
 - name: Chinese Academy of Sciences South America Center for Astronomy (CASSACA), National Astronomical Observatories, Chinese Academy of  Sciences, Beijing 100101, PR China
   index: 2
date: 8 August 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Gaseous exoplanets are far more diverse in mass, density and temperature than their Solar System counterparts. These range from temperate sub-Neptunes whose masses and sizes are below that of the Solar System ice giants, to extremely irradiated, Jupiter-mass exoplanets known as hot Jupiters. The bulk composition of sub-Neptunes and gas giants is key to understanding their formation and evolution pathways, including their dominant formation mechanism [@Thorngren16; @Teske19; @Ginzburg20; @Schneider21; @Mah_Bitsch23] and interactions between the irradiation of the host star and the planetary atmosphere [@OW13; @Acuna22; @Rogers24]. Their building blocks consist of rocks and water (metals), and hydrogen and helium (H/He), which are present in a wide variety of mass fractions in the exoplanet population [@ariel_int_rev]. For exoplanets that contain H/He, the planetary radius (or density) not only depends on mass, but also on the temperature due to the distance from its host star [@Burrows97], its age [@Baraffe03], the mass of the core [@Baraffe08], and the fraction of H/He in its envelope [@Rogers23]. The inference of bulk composition from observational data -- mass, radius, equilibrium temperature, age and atmospheric metallicity -- requires non-trivial numerical methods that solve the equations of planetary interior structure [@Miguel_rev], in addition to Bayesian methods to solve the inverse problem [@Dorn15; @Otegi20; @Acuna21; @Bloot23]. Additionally, interior models need to incorporate the effects of processes occurring in the planet atmosphere, such as radiative transfer [@Baraffe03; @Marley07], and the thermodynamical properties of the different materials, provided by the equation of state (EOS) [@HG23; @Howard23; @gastli_science].

`GASTLI` (GAS gianT modeL for Interiors) is a Python package for interior structure modelling of volatile-rich exoplanets. The package combines a FORTRAN back end with a Python-based user interface. This enables fast convergence of the interior structure to a solution that fulfills the user input, while having a flexible and user-friendly interface. `GASTLI` was specifically designed to compute a wider range of combinations of core mass and H/He content in the envelope compared to existing open-source interior models. GASTLI also incorporates state-of-the-art thermodynamical data and EOSs for rock, water and H/He. One of its key features is the coupling class, which couples the interior structure model to a grid of atmospheric temperature profiles to determine self-consistently the temperature of the outermost layer. GASTLI is released with a default atmospheric grid that assumes equilibrium chemistry and a clear atmosphere at varying internal temperatures computed with petitCODE [@Molliere17; @Molliere15]. Nonetheless, the flexible interface enables the user to use their own custom atmospheric grid, which is useful to explore the effect of the processes assumed in the upper atmosphere in the interior, such as the type of chemistry (equilibrium or disequilibrium), cloud implementations or 3D global circulation (GCM). 


# Statement of need

Interior structure codes are essential to carry out planet composition analyses and Bayesian inverse fitting - known as retrievals - of mass and radius data. MESA [@Paxton11; @Paxton13; @Paxton15; @Paxton18; @Paxton19] is a well-known stellar interior structure software that has been applied to sub-Neptunes and gas giants [@chen_rogers16; @Muller21]. The compositions of the envelope in MESA are limited to low metal contents due to the default EOSs, which is particularly limiting for planets with Neptune's mass or lower. The default MESA EOS can incorporate up to a 4x solar composition in the envelope. For reference, Saturn's upper atmosphere is approximately 10x solar [@Atreya22], while Neptune's measured envelope metallicity is 80x solar [@Neptune_met1; @Irwin21]. Therefore, the implementation of alternative EOS in MESA are necessary, which is not straightforward for the user given the complexity of the code, and the input system based on FORTRAN namelist files - known as inlists. These input namelists allow users to specify initial conditions, control settings, and adjust physical parameters within the MESA software. MESA requires at least one minute of computational time for a single forward model calculation, making it difficult to couple with atmospheric models. This coupling is important because atmospheric conditions can significantly influence the the output of the interior calculations [@LF14; @Linder18; @Poser19; @Poser24]. MAGRATHEA [@magrathea] is another open-source planet interior model that is a computationally fast ($<2$ seconds for a single forward model) and based on C++. Its layer structure (ideal gas EOS on top of liquid/ice water layer) makes its adaptation for applications to sub-Neptunes and gas giants a challenge to new users because these types of planet require the treatment of nonideal, high-pressure H/He [@saumon95; @Chabrier21] and metal (water) [@Thomas_Madhu16; @Mazevet19; @Haldemann20] envelopes. Additionally, ExoInt [@exoint] and ExoPlex [@Unterborn] are interior structure model packages that include Fe-rich cores and rock mantles dedicated to planets without volatile gas species, such as super-Earths and rocky planets ($R < 2 \ R_{\oplus}$).



# The GASTLI Python package

The interior model module is contained in the `GASTLI.int_planet` class. The interior structure equations for hydrostatic equilibrium, Gauss's theorem for gravitational acceleration, convection, and mass conservation are solved by combining integration by the trapezoidal rule with an iteration scheme over all profiles (pressure, temperature, gravity and density). In each iteration, the trapezoidal integrations are carried for each point in a 1D grid that represents the radius, from the center of the planet to the outermost boundary of the top layer. The density is calculated by evaluating the EOS for rock [@sesame; @Miguel22], water [@Mazevet19] and H/He [@Chabrier21; @HG23] at the pressure and temperature conditions of each point in the spatial grid. The iteration sequence is considered to have converged when the difference between iterations is less than a given precision that can be modified by the user. The recommended and default value for this precision is $10^{-5}$. The iteration scheme can also be stopped at a given number of iterations defined by the user, with a default of 30 maximum iterations, and a maximum possible value of 99 iterations. The boundary conditions to solve these equations include the top layer's pressure and temperature, which can be defined by the user. In addition to the profile convergence, the interior module ensures that the surface boundary conditions and the input mass are satisfied. A complete description of the interior module integration scheme can be found in @Brugger_phd_thesis, and has been widely used [@Brugger16; @Brugger17; @Mousis20; @Acuna21; @Aguichine21; @Acuna22; @Vivien22; @Acuna23]. The new feature of the interior module includes the optimization of several of its routines, such as the trapezoidal integration and the evaluation of the EOS's density tables, and the parallelization of the calculation of the interior structure profiles with OpenMP in FORTRAN. One single computation of the interior module (function `calc_radius` in the `GASTLI.int_planet` class) for the default precision and number of iterations takes 1.89 seconds in an Apple M2 CPU, with parallelization on 8 cores (four cores operating at 3.49 GHz and the other four cores operating at 2.42 GHz).

The class `atm_models_interp` encloses the function that interpolates the atmospheric grid to obtain the temperature profile and interior boundary temperature. In addition, this class also provides the functions that calculate the atmospheric thickness and the atmosphere density profile. The latter is obtained by evaluating the EOS used by the interior module for H/He, and the AQUA water EOS [@Haldemann20]. The atmospheric grid class and the interior module class are then combined in the `Coupling` class. This class self-consistently couples the interior module and the atmosphere grids to calculate the boundary temperature of the outermost layer of the interior model. See @Acuna21 for a detailed description of the coupling algorithm, which has been used in previous work [@Aguichine21; @Acuna22; @Acuna23]. Finally, the class `Thermal_evolution` contains the functions necessary to evaluate a sequence of coupled interior-atmosphere models at different internal temperatures and solve the luminosity equation to obtain the internal temperatures and radii as a function of age [@Fortney13; @Poser19]. The complete validation of GASTLI against mass-radius relationships obtained with interior models used in previous work [@Fortney13; @Muller21] can be found in @gastli_science.



![Schematic diagram of GASTLI's compositional structure, inputs and outputs. \label{fig:summary_fig}](GASTLI_summary_fig.pdf)


In an interior retrieval, we fit an interior composition model to data, particularly to mass, radius, age, atmospheric metallicity and/or internal temperature. This is usually done through Bayesian model fitting, such as Markov chain Monte Carlo [@Mosegaard] or nested sampling [@nested_sampling]. The likelihood function in interior retrievals is computed as a function of the squared residuals between the compositional model and the observed values of mass, radius, etc. The likelihood function use in interior retrievals can be found in equations 6 and 14 of @Dorn15 and @Acuna21, respectively. GASTLI can generate a grid of forward models that contains thousands of thermal evolution tracks and explores all possible parameter space for an exoplanet. For example, a grid of 4000 models can be computed in 12 hours, depending on the computer architecture and number of threads/cores available for parallelization, given that each model's computational time is 2-10 seconds. In GASTLI's documentation, we show an example of how to generate a grid and interpolate it to carry an interior retrieval. This example uses as sampler emcee [@emcee]. 

Finally, GASTLI provides a default atmospheric grid for gas giants and sub-Neptunes with equilibrium temperatures under 1000 K, and atmospheric compositions from sub-solar to 250x solar (80% metal mass fraction). Our default grid assumes clear atmospheres in chemical equilibrium. If this grid is not appropriate for a particular target, the user may choose from a variety of model grids that are publicly available. These grids include clear, self-luminous atmospheres in equilibrium [@bobcat; @Linder18; @samland17], clear, self-luminous in disequilibrium [@elfowl; @cholla], cloudy self-luminous [@diamondback; @lacy_burrows; @Jorgensen24], cloudy irradiated gas giants [@charnay18], hot Jupiters [@Molliere15; @Roth24] and temperate sub-Neptunes [@Kempton23]. GASTLI's documentation includes an example of the grid file format and input.


# Similar tools

The following table lists public software tools that provide mass-radius curves and possibly thermal evolution tracks. We distinguish between open-source interior structure packages, and neural networks together with grid interpolators. The latter may use tables and data generated with proprietary interior structure model codes.

| Name and link   | Type of software | Keywords    | Reference             | 
| :-------------- | :--------------- | :---------- | :-------------------- | 
| [MESA](https://docs.mesastar.org/en/latest/modules.html) | Interior structure model | Stars, hot Jupiters, gas giants; H/He and metals | @Paxton11, @Paxton13, @Paxton15, @Paxton18, @Paxton19 |
| [Magrathea](https://github.com/Huang-CL/Magrathea)  | Interior structure model | Super-Earths, rocky planets; liquid water| @magrathea  |
| [ExoInt](https://github.com/astro-seanwhy/ExoInt/tree/master)  | Interior structure model | Dry super-Earths, rocky planets; Fe, rock | @exoint |
| [ExoPlex](https://github.com/CaymanUnterborn/ExoPlex)  | Interior structure model  | Dry super-Earths, rocky planets; Fe, rock | @Unterborn |
| [plaNETic](https://github.com/joannegger/plaNETic)     | Trained neural network | Dry super-Earths, rocky planets; H/He, water, rock, Fe | @Egger24, @Haldemann24|
| [ExoMDN](https://github.com/philippbaumeister/exomdn)  | Trained neural network | Super-Earths, sub-Neptunes; Fe, rock, liquid water, H/He. Retrievals only. | @Baumeister |
| [Exoplanet Composition Interpolator](https://tools.emac.gsfc.nasa.gov/ECI/) | Grid interpolation (web interface)  | Sub-Neptunes; H/He, rock, Fe | @LF14 |
| [HARDCORE](https://github.com/gsuissa/hardCORE?tab=readme-ov-file)       | Grid interpolation (web interface)  | Dry super-Earths, rocky planets; Fe, rock | @Suissa18 |
| [planetsynth](https://github.com/tiny-hippo/planetsynth)| Grid interpolator | Gas giants; H/He, metals | @Muller21  | 
| [SMINT](https://github.com/cpiaulet/smint#readme)|  Grid interpolator   | Sub-Neptunes; H/He, water, rock, Fe | @Piaulet21; Tables from @LF14, @Zeng16, @Zeng19, @Aguichine21   | 
| [mr-plotter](https://github.com/castro-gzlz/mr-plotter/tree/main)| Grid interpolator | Super-Earths, sub-Neptunes; H/He, water, rock, Fe | @Amadeo23; Tables from @LF14, @Zeng16, @Zeng19, @Aguichine21, @Seager07, @Turbet20, @Haldemann24 | 
| [Mardigras](https://github.com/an0wen/MARDIGRAS)|  Grid interpolator | Sub-Neptunes; Water, rock, Fe | @mardigras; Tables from @LF14, @Zeng19, @Aguichine21 | 

# Documentation

Users can find GASTLI's documentation [here](https://gastli.readthedocs.io/en/latest/).


# Acknowledgements

We thank Yamila Miguel for sharing the reformatted EOS tables for silicate (dry sand) from the SESAME database. We acknowledge the support of the Data Science Group at the Max Planck Institute for Astronomy (MPIA) and especially Ivelina Momcheva for her invaluable assistance in optimizing the software and developing the GASTLI Python package. We thank Magali Deleuil, Olivier Mousis, Artyom Aguichine and Bastien Brugger for their contributions to the development of Marseille’s Super-Earth Interior (MSEI) model, on which GASTLI’s development is based. We thank the two JOSS reviewers, Iva Laginja and Ryan MacDonald, for their suggestions to improve the manuscript, as well as the JOSS editor, Josh Borrow. 


# References
