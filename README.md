# Welcome to Vatpy
Vatpy (Visualisation of Arepo in the Terminal using PYthon) is a light-weight, highly customisable, visualisation tool-kit for astrophysical simulations performed using the Arepo code (Springel 2010). 
Many of its functions can be generally applied to simulations made by Arepo (as long as the output is in HDF5-format), however, more specific capabilities, 
such as creating visual maps of the gas chemistry, is at the moment only adapted to simulations run using the ArepoNoctua numerical framework (see Petersson et al. 2025).

![alt text](logo/vatpy_vertical.png)

## Installation
In this section, we provide some general guidelines to install and set up Vatpy on your own machine.

### Requirements
> `python`, `numpy`, `matplotlib`, `mpl_toolkits`, `scipy`, `h5py`, `labellines`, `pycstruct`

### Quick Installation Guide:
- Download the repository to your home directory
- Include the path to your Vatpy directory to the PYTHONPATH variable. For example, you can add the following line to your .bashrc or .bash_profile file, to ensure that this is always done when launching a new session:
  ```shell
  export PYTHONPATH=$PYTHONPATH:$HOME/vatpy
  ```
- Now, go to your Vatpy directory and make a copy of `configv-template.py` to `configv.py` (do NOT move the copied file to any other directory). This file (`configv.py`), is now your personal Vatpy configuration, and can be altered as you wish. However, it is important that before you try using Vatpy for the first time, make sure to change the variable `homedir`, to the actual home directory path on your own machine (where Vatpy should ideally be installed). 
- Voilà, that's it, you should now be ready to use Vatpy in the terminal, Jupyter Notebooks, and Python scripts. To verify that Vatpy is installed correctly, simply try:
  ```python
  import vatpy as vp
  ```
- If you have any issues installing Vatpy, please do not hesitate to contact any of the main contributors listed below, or simply open up an issue on the GitHub repository (https://github.com/AstroJPetersson/vatpy), and we will get back to you as soon as possible. 

## Usage
Here, we provide a short overview of how to use Vatpy directly in the terminal, as well as in Jupyter Notebooks.

### How to use Vatpy in the Terminal
*TODO*

#### GUI Version
*Still under development...*

#### CLI Version
*Still under development...*

### How to use Vatpy in Jupyter Notebooks
*TODO*

## Contributions
**Jonathan Petersson** - *PhD Student @ EPFL* - jonathan.petersson@epfl.ch

## License
Distributed under the MIT license (see LICENSE.md for more information).
