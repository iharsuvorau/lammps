## Overview
This package implements a *fix femocs* command which can be used in a
LAMMPS input script.  This fix can be employed to either do concurrent
coupling of MD with FE-based physics surrogates or on-the-fly
post-processing of atomic information to continuum fields.

There are example scripts for using this package in examples/USER/femocs.

This package uses an external library in lib/femocs which must be
compiled before making LAMMPS. See the instructions below how to do it.

The primary people who created this package are Mihkel Veske
(veskem at gmail.com) and Andreas Kyritsakis (akyritsos1 at gmail.com).
Contact them directly if you have questions.

## Citing
[FEMOCS](https://github.com/veskem/femocs/) is an open-source and freely available
code. The details about its algorithms are published in
[Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2018.04.031) and
[Physical Review E](https://doi.org/10.1103/PhysRevE.101.053307).
When publishing results obtained with the help of FEMOCS, please cite

    Veske, M. et al, 2020. Dynamic coupling between particle-in-cell and atomistic simulations. Physical Review E, 101(5), p.053307.

## Instructions to build LAMMPS with FEMOCS
First of all, FEMOCS needs to be built.
If you are lucky, under Ubuntu it can be done by first going to
*LAMMPS_DIR/src/USER-FEMOCS* directory and then running a command

    $ make install
    
If the installation completes without any error message, the same Makefile could
be used to compile & link FEMOCS and LAMMPS by running in the same directory

    $ make release

However, in case of no succes, both installation and compilation must be done manually.
For this, first go to *LAMMPS_DIR/lib/femocs* and see the instructions in **README.md** for building FEMOCS.
After this, come back to *LAMMPS_DIR/src* and include FEMOCS into LAMMPS build by running

    $ make yes-user-femocs
    
Most probably you plan to use EAM potential as well, therefore run also

    $ make yes-manybody
    
After this, you are ready to start compiling LAMMPS. Under Ubuntu desktop,
it can be done by running

    $ make -j4 ubuntu    # uses 4 CPU cores to speed up compilation
    
It might be needed to change *ubuntu* to something else, especially if you are 
running some other OS than Ubuntu. For instance, as currently LAMMPS-FEMOCS
could be run only with one CPU core, you might try *serial* instead of *ubuntu*.
Other options could be found either by

* running *make help* in *LAMMPS_DIR/src*,
* consulting LAMMPS manual or
* modifing file *LAMMPS_DIR/src/MAKE/Makefile.ubuntu*.

In a fresh install of Ubuntu, the build might fail due to lack of *libfftw*, 
*libjpeg* and/or *libpng* in system. You can install them by running

    $ sudo apt install libfftw3-dev libjpeg-dev libpng-dev

Under other platforms, see the instructions for *libfftw* [here](http://micro.stanford.edu/wiki/Install_FFTW3).

## Usage of FEMOCS in LAMMPS
To use FEMOCS in LAMMPS, you should add the following command to your LAMMPS input script:

    fix FIXNAME GROUPNAME femocs PATH_TO_FEMOCS_INPUT_SCRIPT VERBOSE_FLAG
   
The sample usage of FEMOCS is shown in sample LAMMPS and FEMOCS input scripts at
*LAMMPS_DIR/examples/USER/femocs* directory. In case LAMMPS was built with a flag
*ubuntu*, you can invoke a test simulation in this directory by running

    $ ../../../src/lmp_ubuntu < in.lmp
