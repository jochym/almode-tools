Tools for alamode program to facilitate the calculations
--------------------------------------------------------

The directory structure

    <crystal_name> ---+--- opt
                      |
                      +--- harmonic
                      |
                      +--- TDEP_<volume 1> ---+--- opt
                      |                       |
                      |                       +--- harmonic (may be linked to harmonic upstairs)
                      |                       |
                      |                       +--- T_<temperature 1>
                      |                       |
                      |                       +--- T_<temperature 1>_base
                      |                       |
                      |                       +--- ...
                      |                       |
                      |
                      +--- TDEP_<volume 2>
                      |


The basic procedure:

1. Start with optimized structure in `opt` directory
2. Run the basic harmonic calculation in `harmonic` directory:
    1. Starting from the unit cell as `POSCAR` file
    2. Generate displacements pattern
    3. Generate displacements
    4. Run vasp on displacements
    5. Extract forces
    6. Calculate phonons and Mean Square Displacements (MSD)
3. Calculate TDEP cycle for several volumes in `TDEP_<volume>`


TODO
-----

create common directory in /opt/bin/ for all scripts:
* make-gen
* alm
* run-calcs
* proc-calcs
* make-tdep
* proc-dirs
* extract.py
  
 make short (one line) comments on external parameters expected in the command line