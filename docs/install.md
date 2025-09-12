# Installation

Here, we provide some general guidelines on how to install and set up Vatpy on your personal machine.

## Requirements
*Version requirements are not neccessarily strict. The provided versions are the ones at which Vatpy
have been tested at, and does not mean that Vatpy might run on older python/package versions.*

```
python>=3.12.7
numpy>=1.26.4
matplotlib>=3.9.2
scipy>=1.13.1
h5py>=3.11.0
matplotlib-label-lines>=0.8.1
pycstruct>=0.12.2
cmasher
```

*Only needed to deploy documentation.*

```
mkdocs>=1.6.1
mkdocs-material>=9.6.18
```

## Recommended set-up
*The recommended set-up of Vatpy is still in an early development stage. In the future, we hope to provide a complete package that can be easily installed via pip. But for now, we recommended the following set-up on your selected machine:*

* Download the repository to your home directory, e.g. via HTTPS
```
cd ~
git clone https://github.com/AstroJPetersson/vatpy.git
```
This will create a directory called **vatpy** in your home directory.

* Now, enter the **vatpy** directory and copy `template-configv.py` to `configv.py`, i.e.
```
cd vatpy
cp template-configv.py configv.py
```
The file `configv.py` contains your personal Vatpy configuration. Within this file, please make sure to change the `homedir` variable to the actual home directory on your personal machine.

* As a final step to make Vatpy work on your personal machine, you will have to export the path to your recently downloaded **vatpy** directory to your PYTHONPATH variable. For example, to ensure that this is always done when starting up a new session, you can e.g. add the following line to your `~/.bash_rc` or `~/.bash_profile` file
```
export PYTHONPATH=$PYTHONPATH:$HOME/vatpy
```

* Voilà, that's it, you should now be ready to use Vatpy via Python scripts directly in the terminal, or inside Jupyter Notebooks. To make sure that Vatpy is installed correctly, simply try:
```
import vatpy as vp
```

* If you encounter any issues when installing Vatpy, please do not hesitate to contact any of the main contributors to the code, or open up an issue on the GitHub repository (<https://github.com/AstroJPetersson/vatpy>).

<br><br>
Last updated: 2025-09-04
