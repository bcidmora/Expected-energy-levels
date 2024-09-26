# Expected energy levels

To run this script, the following libraries are needed:
    + pylatex
    + numpy
    + os
    + matplotlib

These scripts are organized as:

** expected_energy_levels.py ** 
    This is the main script to run. The directory at the top must be changed. This is where the Expected energy Levels from Colin are stored. Some details of other variables are:
    
    - unkwnownHads: if True, the unknown masses for a certain ensemble will be include in the calculation of energy levels. It takes the values from E250.
    
    - noEnergyCut: If False, the full list of energy levels is reduced to those that are within the chosen cutoff, otherwise it includes all possible levels.
    
    - threeParticleThreshold: if True, three particles must be chosen and then it cuts the energy two levels after this threshold. 
    
    - cmEnergy: if True the final energy levels are in terms of the E_{cm}.
    
    - plotEnergyLevels: If True, all the selected energy levels are plotted at the end of the section.
    
    - ensChoices: Ensembles to be studied. More ensembles can be added to the dictionary in "ensembles_info.py"
    
    - chosenThreshold: The cutoff value (three particles, a specific value or no cut)
    
    - hadronsInfo: This list contains all the hadrons that are known in all the ensembles in "ensembles_info.py". More hadron values can be added, but their mass values must be included for the ensemble in "ensembles_info.py". 
    
    Later, the process is the following: 
        * All energy values are computed using Colin's group theory (ENERGY_LIST_RAW). 
        * Then, a selection of levels is made according to the chosen cutoff (ENERGY_LIST_TABLES). 
        * Finally, all tables are constructed according to this selection, if the tables are too long, they are separated in several pages (CONSTRUCTING_TABLES). 
        * When plotEnergyLevels==True, the plots can include thresholds.    

** tables_latex.py **
    Builds all the tables that appear in the pdf file. These tables are written in a TeX file and compiled. It also has the plotting routine to save the plot in pdf.
    
** energy_functions.py **
    It has all the routines that calculate stuff, read from Colin's files and rewrites them, Menus, etc. 
    
** ensembles_info.py **
    This file contains all the information about the ensembles one wants to get the energy levels from. Please add more hadron masses in this file if you have more values, and share them with the group.
    
