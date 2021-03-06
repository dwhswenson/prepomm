# prepomm
Tools to prepare an OpenMM simulation

This is intended to be a simple package for getting OpenMM simulations up
and running; i.e., this includes some scripts to help get from a bare file
(like a .pdb or a .gro) to running MD.  My general approach is to use a
Jupyter notebook to do the setup interactively. I try to ensure that all
steps are achieved programmatically, in order to facilitate reproducibility.

The beauty of OpenMM is that, as a library, it is customizable to your
individual needs. Therefore, many people end up writing little wrappers that
facilitate their own workflows. These are my wrappes.  They are intended to
simplify some common tasks, but will not work with all setups. In many cases
you may need to write custom setup code, but that can often be based on code
in `prepomm`.

This package may rely on (and certainly does not replace) several other
useful packages for working with OpenMM, including:

* [MDTraj](http://mdtraj.org)
* [ParmEd](http://parmed.github.io/ParmEd/)
* [pdbfixer](https://github.com/pandegroup/pdbfixer)
* [openmmtools](https://openmmtools.readthedocs.io/)

Mostly, this chains things together from those packages (and from OpenMM), and
it provides a few conveniences that aren't built into any of those tools.

## Features

* `FFModel`: convenience class for working with force fields
* `max_atom_distance`: fast MDTraj-based code for identifying maximum
  distance between atoms (useful to setting box size)
* functions for creating rhombic dodecahedral box vectors
* automatic protocols for adding hydrogens/solvating, as well as running
  various equilibration simulations
* PLANNED: sanity checks on your PDB (correct termination, reasonable box
  size, etc.) These will give better error messages than the default OpenMM
  errors.

Note that if these features are included into other packages, the
implementation here may be removed. This repository isn't intended to be a
permanent project; just a central location to maintain a set of useful
tools that I couldn't find elsewhere.

## Installation

<!--This package requires OpenMM, which you should install with `conda`. The-->
<!--rest of the requirements can be installed with either `conda` (recommended)-->
<!--or with `pip`. If you are missing requirements (other than OpenMM), the-->
<!--`pip` commands below will attempt to install the requirements.-->

<!--The simplest way to install is to ask `pip` to install from source:-->

<!--```bash-->
<!--pip install git+http://github.com/dwhswenson/prepomm.git-->
<!--```-->

The recommended installation is to do a developer install. Clone it (or your
fork of it, if you think you might contribute changes), change into the
directory with `setup.py`, and run:

```bash
pip install -e .
```

Depending on your Python installation, you may need to preface that command
with `sudo`.

<!--You can also install a developer version. To do that, we recommend forking-->
<!--the project, and cloning from your own fork. Then you can install with `pip-->
<!--install -e .` from within the directory with the `setup.py` file.-->
