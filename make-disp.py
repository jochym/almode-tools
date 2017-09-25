#!/usr/bin/env python3

import sys
import ase.io
import click

from numpy.random import normal
from numpy import sqrt, array, loadtxt, savetxt, concatenate
from glob import glob
from os import path
import scipy.interpolate
from collections import OrderedDict

def read_elm_order(cr):
    d=OrderedDict(zip(cr.get_chemical_symbols(),cr.get_atomic_numbers()))
    for n, k in enumerate(d):
        d[k]=n
    return d

@click.command()
@click.option('-p', '--poscar', default='SPOSCAR', 
                type=click.Path(exists=True),
                help='Supercell POSCAR file (SPOSCAR)')
@click.option('-m', '--msdfile', default='base/QHA.msd', 
                type=click.Path(exists=True), 
                help='Base calc mean squere displacement file (base/QHA.msd)')
@click.option('-T', '--temp', default=300, help='Temperature (300)')
@click.option('-n', '--number', default=50, help='Number of configurations (50)')
@click.argument('outfile')
def disp(poscar, msdfile, temp, outfile, number):
    """Generates thermal displacement files """

    cr = ase.io.read(poscar)
    elems=read_elm_order(cr)
    msd=loadtxt(path.join(msdfile)).T
    dsp=scipy.interpolate.interp1d(msd[0],msd[1:])(temp)
    dsp=dsp.reshape((-1,3))
    basepos=cr.get_positions()
    print('Generating displacements: ', end='')
    for k in range(number):
        print(k, end=' ')
        rnd=normal(size=basepos.shape)
        np=array([sqrt(dsp[elems[z]])*r+p
                    for z, p, r in zip(cr.get_chemical_symbols(), basepos, rnd)])
        cr.set_positions(np)
        ase.io.vasp.write_vasp('%s_%03d.POSCAR' % (outfile,k), cr, direct=True, vasp5=True)
        savetxt('%s_%03d.disp' % (outfile,k), concatenate((rnd, np-basepos), axis=1),
                fmt='%8.4f', header=' '.join(['%f' % sqrt(v) for v in dsp.reshape(-1)]))
    print()
if __name__ == '__main__':
    disp()
