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
labellines
pycstruct
cmasher
```

*Only needed to deploy documentation.*

```
mkdocs>=1.6.1
mkdocs-material
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


