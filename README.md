Tools for alamode program to facilitate the calculations
--------------------------------------------------------

The directory structure

    <crystal_name> ---+--- opt
                      |
                      +--- harmonic
                      |
                      +--- TDEP_<volume_1>
                      |
                      +--- TDEP_<volume_2>
                      |


The basic procedure:

1. Start with optimized structure in `opt` directory
2. Run the basic harmonic calculation in `harmonic` directory
    a. Starting from the unit cell as POSCAR file
    b. Generate displacements pattern
    c. Generate displacements
    d. Run vasp on displacements
    e. Extract forces
    f. Calculate phonons and Mean Squere Displacements (MSD)
3. Calculate TDEP cycle for several volumes in `TDEP_<volume>`



