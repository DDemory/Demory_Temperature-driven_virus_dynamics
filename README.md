# Demory_Temperature-driven_virus_dynamics

David Demory and collaborators, August 18, 2020 School of Biological Sciences Georgia Institute of Technology

**Code for:** A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics.

A preprint of the manuscript can be found on BioRxiv: [link](https://doi.org/10.1101/2020.08.18.256156)

This code is archived on Zenodo: [![DOI](https://zenodo.org/badge/288514967.svg)](https://zenodo.org/badge/latestdoi/288514967)

**Instructions:** all codes are written in MATLAB. To plot the figures, you need (1) to modify the folder path in path_setup.m. and (2) run path_setup.m before running other codes.

**Folder descriptions:**
- **data:** contains experimental data from Demory et al. 2017. (1st column: time in days and 2d column: concentration in particles/ml) and SST projection from IPCC.
- **figures:** contains matlab figure codes.
- **functions:** contains additional function codes needed to plot the figures: ODE system, etc., temperature-driven functions, etc.
- **outputs:** figure outputs in eps.
- **parameters:** hyper-parameters of the temperature-driven functions and ecological state biogeography estimates.

**Additional packages:** Figure 5 was generated using the M_Map package (Pawlowicz 2020). You need to download it at https://www.eoas.ubc.ca/~rich/map.html and add the m_map folder to your directory.

**Figure License:** CC-BY-4.0 Code License: MIT

**References:**
- Demory et aL., 2017. Temperature is a key factor in Micromonas-virus interactions, The ISME journal 11(3), 601-612.
- Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.

