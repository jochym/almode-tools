#!/usr/bin/env python3

import numpy as np
from numpy import array, zeros, arange, savetxt
from numpy import dot
from numpy.linalg import inv
from numpy.random import choice
import ase
import ase.io
import ase.units as units
from ase.atoms import Atoms
import glob
import click
import sys

def read_potim(fn):
    pat='0123456789.'
    with open(fn) as f:
        for l in f:
            if 'POTIM' in l:
                potim=[c for c in l if c in pat]
                return float(''.join(potim))
    return None    

def normalize_traj(c, tr):
    base=Atoms(c)
    cell=tr[0].get_cell()
    base.set_cell(cell, scale_atoms=True)
    
    spos=array([a.get_scaled_positions() for a in tr])

    # Calculate unwrapped step-by-step atom displacements
    sdx=spos[1:]-spos[:-1]
    sht=(sdx < -0.5)*1 - (sdx > 0.5)*1
    sdx+=sht
    # Check if step-to-step fractional displacements are below 1/3
    assert (abs(sdx) < 1/3).all()

    # All below computed in fractional coordinates 
    
    # Calculate step-to-step CM drift and integrate it to CM trajectory
    # CM displacements
    sdcm=(sdx*base.get_masses()[None,:,None]).sum(axis=1)/base.get_masses().sum()
    # Integrate (sum) the displacements
    stcm=sdcm.cumsum(axis=0)
    # Remove CM drift wrapped for PBC
    spos[1:]-=mod(stcm[:,None,:],1)
    
    # Unwrap the initial position relative to base
    sp=spos-base.get_scaled_positions()
    spsht=(sp < -0.5)*1 - (sp > 0.5)*1
    spos+=spsht
    # Calculate mean CM position along the trajectory
    scm=(spos*base.get_masses()[None,:,None]).sum(axis=1).mean(axis=0)/base.get_masses().sum()
    # Move the mean CM position to the reference CM position
    spos+=(dot(base.get_center_of_mass(),inv(cell)) - scm)
    
    # Return carthesian positions, fractional positions, fractional cm traj, cm traj
    return dot(spos,cell), spos, sdx, stcm, dot(stcm,cell)


def read_traj(c0, dn, fn='vasprun.xml'):
    '''
    Read the trajectory from vasprun.xml or OUTCAR
    Calculate the velocities (unwrapping PBC jumps)
    Remove the CM drift from the system
    Return unwrapped positions, fractional positions, velocities and CM drift.
    The c0 is used as a reference system.
    '''
    print(f'Reading {dn}',end=':')
    tr=ase.io.read(dn+f'/{fn}',index=':')
    dt=read_potim(dn+f'/{fn}')*ase.units.fs
    print(f'{len(tr)}')
    base=Atoms(c0)
    cell=tr[0].get_cell()
    base.set_cell(cell, scale_atoms=True)
    
    pos, spos, sdx, stcm, tcm = normalize_traj(base, tr)
    
    # Calculate the central difference velocities
    v=zeros(pos.shape)
    v[:-1]=dot(sdx,cell)/dt
    v[0]=v[1]
    v[1:-1]=(v[:-2]+v[2:])/2
    
    # Store the data
    return {
        'base':base,
        'dt':dt,
        'fname': dn,
        'trj': tr,
        'pos': pos,
        'spos': spos,
        'vel': v,
        'drift':tcm,
    }



@click.command()
@click.option('-p', '--poscar', default='SPOSCAR', type=click.Path(exists=True), help='Supercell POSCAR file (SPOSCAR)')
@click.option('-t', '--traj', default='vasprun.xml', type=click.Path(exists=True), help='Trajectory file OUTCAR or vasprun.xml (vasprun.xml)')
@click.option('-n', '--number', default=50, help='Number of configurations (50)')
@click.option('-d', '--disp', default='disp.dat', help='Displacement file (disp.dat)')
@click.option('-f', '--force', default='force.dat', help='Forces file (force.dat)')
@click.option('-c', '--configs', is_flag=True, default=False, help='Write disp_xxx.POSCAR files')
def make_disp_force(poscar, traj, number, disp, force, configs):
    """Generates displacement and forces files from the MD trajectory"""
    base=ase.io.read(poscar)
    md=read_traj(base, traj)
    pos=md['pos']
    base=md['base']
    tr=md['trj']
    p0=pos.mean(axis=0)
    n=0
    print(f'#Writing {disp} and {force} files from {number} random steps out of {traj} file.')
    print('#Time steps:')
    with open(disp,'bw') as dsp, open(force, 'bw') as frc:
        for k in sorted(choice(arange(pos.shape[0]), number, replace=False)):
            savetxt(dsp, (pos[k]-p0)*units.Bohr, fmt='%24.18f')
            savetxt(frc, tr[k].get_forces()*units.Ry/units.Bohr, fmt='%24.18f')
            print(f'{k}', end=' ')
            sys.stdout.flush()
            if configs :
                a=Atoms(tr[k])
                a.set_pbc(False)
                a.set_positions(pos[k])
                ase.io.write(f'disp_{n:04d}_{k:04d}.POSCAR', a, direct=True, vasp5=True)
            n+=1    
    print()
    print(f'#Finished: {number} configs extracted')
if __name__ == '__main__':
    make_disp_force()