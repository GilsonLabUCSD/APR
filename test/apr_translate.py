#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
import re

def setup_translate(distance, ligand_residue):
    """
    Move the guest away from the host (zero distance) is used for the attachment phase.
    :param distance: distance between the guest and host
    :param host_atoms: the number of atoms in the host molecule
    :return:
    """

    print('Translating...')
    # Open PDB file and read coordinates
    cols_before_coords = []
    cols_after_coords = []
    ter_row = []
    coords = []
    ter_list = []
    total_atom = 0

    f =  open('align_z.pdb')

    for line in f:
        line = line.rstrip()
        splitdata = line.split()
        if (splitdata[0]=='ATOM')or(splitdata[0]=='HETATM'):
            total_atom += 1
            x_coord = float(line[30:38].strip())
            y_coord = float(line[38:46].strip())
            z_coord = float(line[46:54].strip())
            cols_before_coords.append(line[0:30])
            cols_after_coords.append(line[54:81])

            resname = line[17:20].strip()
            resid = line[22:26].strip()

            if resname in ligand_residue or resid in ligand_residue:
               z_coord += distance
            coords.append((x_coord, y_coord, z_coord))

        elif splitdata[0] == 'TER':
            ter_list.append(total_atom)
            ter_row.append(line)

    # Write the new pdb file
    newPDB_file = open('vac.pdb', 'w')
    j = 0
    for i in range(len(coords)):
        newPDB_file.write('%s%8.3f%8.3f%8.3f%s\n'%(cols_before_coords[i], coords[i][0], coords[i][1], coords[i][2],cols_after_coords[i]))
        if i+1 in ter_list:
            newPDB_file.write(ter_row[j]+'\n')
            j+=1 

    if total_atom not in ter_list:
        newPDB_file.write('TER\n')

    newPDB_file.write('END\n')

    f.close()
    newPDB_file.close()
