#  Benchmarking macaque brain gene expression for horizontal and vertical translation.
Authors: A.I. Luppi, Z-Q Liu,  J.Y. Hansen,  R. Cofre,  E. Kuzmin,  S. Froudist-Walsh, N. Palomero-Gallagher,  & B. Misic.

This repository provides code to reproduce the main results of Luppi et al., "Benchmarking macaque brain gene expression for horizontal and vertical translation." bioRxiv (2024) ([preprint](https://doi.org/10.1101/2024.08.18.608440)).

It was developed in MATLAB 2019a by Andrea Luppi from the the [Network Neuroscience Lab](netneurolab.github.io/) at the Montreal Neurological Institute, McGill University.

This code relies on MATLAB code from the [BrainSpace Toolbox](https://brainspace.readthedocs.io/en/latest/) for MATLAB by Vos de Wael et al. (2020) Communications Biology.
For additional plotting functionality, also include in your MATLAB path the [ENIGMA Toolbox](https://github.com/MICA-MNI/ENIGMA.git) by Lariviere et al. (2021) Nature Methods.

The study investigates how network control energy in the brain depends on the network, but also the start and destination states, and the control strategy.
In this code, we change these three elements one at a time.

## Repository Structure
### Main script
The main file is [Luppi_macaque_brain_gene_translation_code_4GitHub.m](Luppi_macaque_brain_gene_translation_code_4GitHub.m)
This script should work out of the box, if run from the parent directory. 
To run, ensure you are in the main directory of the repo.

### `data`
The [data](data/) folder contains all the data you need to make this code run. 

### `utils`
The [utils](utils/) folder contains support functions.
