resource-heterogeneity
======================

This repository contains code, configuration files, data, and analysis scripts for the paper "Spatial resource heterogeneity increases diversity and evolutionary potential" by Emily Dolson, Samuel Perez, Randal Olson, and Charles Ofria.

# Repository contents:

### configs:
Contains all configuration files given to Avida to run the experiments described in the paper. The most important of these configuration files are the environment files, located in the **environmentFiles** subdirectory. The name of each environment file contains the random seed of the run of Avida for which it was used. Environment files that contain task names rather than random seeds were used for the niche construction experiments. The **configs** directory also contains a copy of the Avida binary used to run these experiments. Source code for Avida is available here: https://github.com/devosoft/avida. Lastly, the **configs** direcotry contains a subdirectory called **qsub_scripts**, which contains the bash programs that were submitted to the MSU HPCC to run the experiments in this paper.

### data:
Contains all data presented in this paper.

### exploratory_analysis:
Contains analysis scripts that were not used for the final analysis presented in the paper, but nevertheless provided useful insight into the data (and data from various preliminary experiments).

### figs:
Figures presented in the paper.

### final_analysis:
Contains all of the scripts used to perform the statistical analyses discussed in the paper.

### scripts:
Scripts that were used to run the experiments. The most important of these is genHeterogenousEnvironment.py, and genHeterogenousEnvironment2.py (an updated version written to handle the randomization experiments better), which generated all of the environment configuration files. Most of the code in scripts has since been incorporated into the avidaspatial Python library: https://github.com/emilydolson/avida-spatial-tools.

### tutorial:
A tutorial on how to carry out the experiments in the first half of the paper.

