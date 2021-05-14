# Demory_Temperature-driven_virus_dynamics

David Demory and collaborators, 20 April 2021, School of Biological Sciences Georgia Institute of Technology

**Code for:** A thermal trade‐off between viral production and degradation drives virus‐phytoplankton population dynamics
David Demory, Joshua S. Weitz, Anne‐Claire Baudoux, Suzanne Touzeau, Natalie Simon, Sophie Rabouille, Antoine Sciandra and Olivier Bernard.
Ecology Letters, doi: 10.1111/ele.13722

This code is archived on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4445347.svg)](https://doi.org/10.5281/zenodo.4445347)

**Instructions:** all codes are written in MATLAB. To plot the figures or do the analyses, you need (1) to modify the folder path in path_setup.m. and (2) run path_setup.m before running other codes.

**Folder descriptions:**
- **data:** contains experimental data from Demory et al. 2017. (1st column: time in days and 2d column: concentration in particles/ml) and SST projection from IPCC.
- **analysis:** contains matlab analysis codes.
- **figures:** contains matlab figure codes.
- **functions:** contains additional function codes needed to plot the figures: ODE system, temperature-driven functions, etc.
- **outputs:** figure outputs in eps.
- **parameters:** hyper-parameters of the temperature-driven functions and ecological state biogeography estimates.

**Additional packages:** Figure 5 was generated using the M_Map package (Pawlowicz 2020). You need to download it at https://www.eoas.ubc.ca/~rich/map.html and add the m_map folder to your directory.

**Figure License:** CC-BY-4.0 Code License: MIT

**References:**
- Demory et aL., 2017. Temperature is a key factor in Micromonas-virus interactions, The ISME journal 11(3), 601-612.
- Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.

