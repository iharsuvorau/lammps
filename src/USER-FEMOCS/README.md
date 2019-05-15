## Overview
This package implements a *fix femocs* command which can be used in a
LAMMPS input script.  This fix can be employed to either do concurrent
coupling of MD with FE-based physics surrogates or on-the-fly
post-processing of atomic information to continuum fields.

There are example scripts for using this package in examples/USER/femocs.

This package uses an external library in lib/femocs which must be
compiled before making LAMMPS. See the lib/femocs/README.md for the 
instructions to build FEMOCS. See the LAMMPS manual for information 
on building LAMMPS with external libraries. The settings in the 
femocs/share/makefile.femocs file must be correct for LAMMPS to build 
correctly with this package installed.

The primary people who created this package are Mihkel Veske
(veskem at gmail.com) and Andreas Kyritsakis (akyritsos1 at gmail.com).
Contact them directly if you have questions.

## Citing
[FEMOCS](https://github.com/veskem/femocs/) is an open-source and freely available code. The details about its algorithms are published in
[Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2018.04.031). When publishing results
obtained with the help of FEMOCS, please cite

    Veske, M. et al, 2018. Dynamic coupling of a finite element solver to large-scale atomistic simulations. Journal of Computational Physics, 367, pp.279–294.

## Instructions to build LAMMPS with FEMOCS
First of all, FEMOCS needs to be built. For this, go to *LAMMPS_DIR/lib/femocs*
and see the instructions in README.md.

After this, come back to *LAMMPS_DIR/src* and include FEMOCS into LAMMPS build by running

    $ make yes-user-femocs
    
Most probably you plan to use EAM potential as well, therefore run also

    $ make yes-manybody
    
After this, you are ready to start compiling LAMMPS. Under Ubuntu desktop, it can be done by running

    $ make -j8 ubuntu    # uses 8 CPU cores to speed up compilation

It might be needed to change *ubuntu* to something else, especially if you are running some other OS than Ubuntu.
You can find better option either by

* running *make help*,
* consulting LAMMPS manual or
* modifing file *LAMMPS_DIR/src/MAKE/Makefile.ubuntu*.

In a fresh install of Ubuntu, the build might fail due to lack of *libfftw*, *libjpeg* and/or *libpng* in system.
You can install them by running

    $ sudo apt install libfftw3-dev libjpeg-dev libpng-dev

Under other platforms, see the instructions for *libfftw* [here](http://micro.stanford.edu/wiki/Install_FFTW3).

## Testing and usage of FEMOCS
To test FEMOCS installation, run
   
    $ cd LAMMPS_DIR/lib/femocs
    $ make release
    $ ./build/femocs
   
The last command invokes small FEMOCS test program. If it completes without any
error message, FEMOCS is installed properly. If not, see README.md in FEMOCS directory.

To use FEMOCS in LAMMPS, you should add the following command to your LAMMPS input script:

    fix FIXNAME GROUPNAME femocs PATH_TO_FEMOCS_INPUT_SCRIPT
   
The sample usage of FEMOCS is shown in sample LAMMPS and FEMOCS input scripts at *LAMMPS_DIR/examples/USER/femocs* directory.
In case LAMMPS was built with a flag *ubuntu*, you can invoke a test simulation in this directory by running

    $ ../../../src/lmp_ubuntu < in.lmp
   
## Instructions to build AtC library
Below are the instructions for building AtC library that is similar yet
different to FEMOCS.  
  
To use AtC, first compile it as a library. Use the same compiler as for building lammps.
Command make ubuntu is using mpic++ (see lammps/src/MAKE/MACHINES/Makefile.ubuntu). Therefore
  cd lammps/lib/atc
  make -f Makefile.mpic++
  
After this there should be libatc.a and Makefile.lammps in atc directory. The last must be altered:
  
    $ rm Makefile.lammps
    $ cp Makefile.lammps.installed Makefile.lammps

For further details, see lammps/lib/atc/README

To use libatc.a in LAMMPS, do

    $ cd lammps/src
    $ make yes-user-atc
    $ make yes-manybody
    $ make -j8 ubuntu      # uses 8 CPU cores to speed up compilation
  
## Miscellanneous notes
To get an idea of changing the forces, see src/fix_setforce.cpp/h

To use the code that is in USER-MISC, first copy *.cpp and *.h files there. After this
 
    $ cd lammps/src
    $ make yes-user-misc
    $ make -j8 ubuntu
