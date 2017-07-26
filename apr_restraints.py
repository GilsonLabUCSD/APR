#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys


def setup_restraints(prefix, trans_dist, rest_weight, scale_w, H1, H2, H3, G1, G2, dist_fc,
                     angle_fc, jacks_fc, jacks_dist, jacks_list, strip, d1_resid):
    """
    Add restraints between host and guest atoms.
    :param prefix: whether 'a'ttach or 'p'nattach restraints.
    :param trans_dist: the distance between the host and guest (may be zero)
    :param rest_weight: the weight of the restraint, set up to scale for the attachment phase
    :param scale_w
    :param H1: host atom 1
    :param H2: host atom 2
    :param H3: host atom 3
    :param G1: guest atom 1
    :param G2: guest atom 2
    :param dist_fc: force constant for the distance restraints
    :param angle_fc: force constant for the angle restraints
                      the target value should be in degrees but the restraint should be in radians (because AMBER)
    :param jacks_fc: force constant for the conformational restraints (optional)
    :return:
    """
    print('Adding restraints ...')
    # Get the atom serial numbers of the host, guest and dummy atoms for imposing restraints.
    idx_r1 = find_index(H1, 'dry.pdb')
    idx_r2 = find_index(H2, 'dry.pdb')
    idx_r3 = find_index(H3, 'dry.pdb')
    idx_l1 = find_index(G1, 'dry.pdb')
    idx_l2 = find_index(G2, 'dry.pdb')
    idx_d1 = find_index(':DUM@Pb','dry.pdb')
    idx_d2 = str(int(idx_d1) + 1)
    idx_d3 = str(int(idx_d2) + 1)

    # Get the residue serial numbers of the dummy atoms
    d2_resid = str(int(d1_resid) + 1)
    d3_resid = str(int(d2_resid) + 1)     

    D1 = ':'+str(d1_resid)+'@Pb'
    D2 = ':'+str(d2_resid)+'@Pb'
    D3 = ':'+str(d3_resid)+'@Pb'

    # These should scale like a lambda (from 0 to 1) and gradually turn on during the attachment phase.

    if prefix == 'r':
        dwt = dist_fc
        awt = angle_fc
    else:
        dwt = rest_weight * dist_fc * 0.01
        awt = rest_weight * angle_fc * 0.01

    jwt = rest_weight * jacks_fc * 0.01

    # Write the cpptraj input file that will print the reference values 
    pt_file = open('reference.in', 'w')
    pt_file.write('#AnchorAtoms: %9s %9s %9s %9s %9s\n' % (H1, H2, H3, G1, G2))
    pt_file.write('parm solvated.prmtop\n')
    pt_file.write('trajin solvated.inpcrd\n')
    pt_file.write('distance r0 %9s %9s noimage out reference.dat\n' % (D1, H1))
    pt_file.write('angle r1 %9s %9s %9s out reference.dat\n' % (D2, D1, H1))
    pt_file.write('dihedral r2 %9s %9s %9s %9s out reference.dat\n' % (D3, D2, D1, H1))
    pt_file.write('angle r3 %9s %9s %9s out reference.dat\n' % (D1, H1, H2))
    pt_file.write('dihedral r4 %9s %9s %9s %9s out reference.dat\n' % (D2, D1, H1, H2))
    pt_file.write('dihedral r5 %9s %9s %9s %9s out reference.dat\n' % (D1, H1, H2, H3))
    pt_file.write('distance r6 %9s %9s noimage out reference.dat\n' % (D1, G1))
    pt_file.write('angle r7 %9s %9s %9s out reference.dat\n' % (D2, D1, G1))
    pt_file.write('angle r8 %9s %9s %9s out reference.dat\n' % (D1, G1, G2))

    pt_file.close()
    p = sp.call('cpptraj -i reference.in > reference.log', shell=True)

    # Read the outfile generated by cpptraj "reference.dat"
    dat_file = open('reference.dat', 'r')
    data = dat_file.readline()
    data = dat_file.readline()
    # Obtain the values for each restraints
    reference_values = data.split()
    dat_file.close()

    guest_end = '&end #GuestTR\n'
    host_end = '&end #HostTR\n'
    jacks_end = '&end #HostJ\n'

    # Recast some floats as strings, to prevent python from complaining
    rest_weight = str(rest_weight)
    trans_dist = str(trans_dist)

    # Write those values into the disang.rest file
    rest_file = open('disang.rest', 'w')
    rest_file.write('#AnchorAtoms %9s %9s %9s %9s %9s #Type %s #Weight %s #TransDist %s #ScaleW %s\n' \
                    % (H1, H2, H3, G1, G2, prefix, rest_weight, trans_dist, str(scale_w)))
    # r0
    rest_file.write('&rst iat= %-10s' % (idx_d1+','+ idx_r1+','))
    rest_file.write('%16s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', 0.0, float(reference_values[1]), float(reference_values[1]), 999.0, dist_fc,
        dist_fc) + host_end)
    # r1
    rest_file.write('&rst iat= %-15s' % (idx_d2+','+idx_d1+','+idx_r1+','))
    rest_file.write('%11s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', 0.0, float(reference_values[2]), float(reference_values[2]), 180.0, angle_fc,
        angle_fc) + host_end)
    # r2
    rest_file.write('&rst iat= %-20s' % (idx_d3+','+idx_d2+','+idx_d1+','+idx_r1+','))
    rest_file.write('%6s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', float(reference_values[3]) - 180.0, float(reference_values[3]), float(reference_values[3]),
        float(reference_values[3]) + 180,
        angle_fc, angle_fc) + host_end)
    # r3
    rest_file.write('&rst iat= %-15s' % (idx_d1+','+idx_r1+','+idx_r2+','))
    rest_file.write('%11s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', 0.0, float(reference_values[4]), float(reference_values[4]), 180.0, angle_fc,
        angle_fc) + host_end)
    # r4
    rest_file.write('&rst iat= %-20s' % (idx_d2+','+idx_d1+','+idx_r1+','+idx_r2+','))
    rest_file.write('%6s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', float(reference_values[5]) - 180.0, float(reference_values[5]), float(reference_values[5]),
        float(reference_values[5]) + 180,
        angle_fc, angle_fc) + host_end)
    # r5
    rest_file.write('&rst iat= %-20s' % (idx_d1+','+idx_r1+','+idx_r2+','+idx_r3+','))
    rest_file.write('%6s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', float(reference_values[6]) - 180.0, float(reference_values[6]), float(reference_values[6]),
        float(reference_values[6]) + 180,
        angle_fc, angle_fc) + host_end)
    # r6
    rest_file.write('&rst iat= %-10s' % (idx_d1+','+idx_l1+','))
    rest_file.write('%16s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', 0.0, float(reference_values[7]), float(reference_values[7]), 999.0, dwt, dwt) + guest_end)
    # r7
    rest_file.write('&rst iat= %-15s' % (idx_d2+','+idx_d1+','+idx_l1+','))
    rest_file.write('%11s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', 0.0, 180.0, 180.0, 180.0, awt, awt) + guest_end)
    # r8
    rest_file.write('&rst iat= %-15s' % (idx_d1+','+idx_l1+','+idx_l2+','))
    rest_file.write('%11s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' % (
        'r1', 0.0, 180.0, 180.0, 180.0, awt, awt) + guest_end)

    # if conformational restraints are applied, write extra lines
    if jacks_fc:
        for i in range(0, len(jacks_list) / 2):
            rest_file.write('&rst iat= %-11s' % (find_index(jacks_list[i * 2])+','+find_index(jacks_list[i * 2 + 1])+','))
            rest_file.write('%16s= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, ' \
                            % ('r1', 0.0, jacks_dist, jacks_dist, 999.0, jwt, jwt) + jacks_end)

    rest_file.close()

    # Write restraints.in 
    pt2_file = open('restraints.in', 'w')
    pt2_file.write('#AnchorAtoms: %9s %9s %9s %9s %9s\n' % (H1, H2, H3, G1, G2))
    pt2_file.write('noexitonerror\n')
    if strip == 'yes':
        pt2_file.write('parm vac.prmtop\n')
    else:
        pt2_file.write('parm solvated.prmtop\n')
    for i in range(1, 21):
        pt2_file.write('trajin traj.%02d\n' % (i))

    pt2_file.write('distance d6 %9s %9s noimage out restraints.dat\n' % (D1, G1))
    pt2_file.write('hist d6,0,30,0.05,* norm out hist.dat\n')
    pt2_file.write('angle a7 %9s %9s %9s out restraints.dat\n' % (D2, D1, G1))
    pt2_file.write('hist a7,0,180,0.5,* norm out hist.dat\n')
    pt2_file.write('angle a8 %9s %9s %9s out restraints.dat\n' % (D1, G1, G2))
    pt2_file.write('hist a8,0,180,0.5,* norm out hist.dat\n')

    # if conformational restraints are applied, write extra lines
    if jacks_fc:
        for i in range(0, len(jacks_list) / 2):
            pt2_file.write('distance d%-2d %9s %9s noimage out restraints.dat\n' % (
                i + 7, jacks_list[i * 2], jacks_list[i * 2 + 1]))
    pt2_file.close()
    return 

def return_restraints_for_error_analysis(prefix, trans_dist, H1, H2, H3, G1, G2, dist_fc,
                     angle_fc, jacks_fc, jacks_dist, jacks_list, jacks):
        # Assume that theere is a reference.dat file already...
        dat_file = open('reference.dat', 'r')
        data = dat_file.readline()
        data = dat_file.readline()
        # Obtain the values for each restraints
        reference_values = data.split()
        dat_file.close()
        # These should scale like a lambda (from 0 to 1) and gradually turn on during the attachment phase.
        #if prefix == 'r':
        #    dwt = dist_fc
        #    awt = angle_fc
        #else:
        dwt = dist_fc #* 0.01
        awt = angle_fc #* 0.01
        jwt = jacks_fc #* 0.01
        #print dist_fc,dwt,awt,jwt
        if jacks == 'no':
            return reference_values[7], dwt, 180.0, awt, 0.0, 0.0, 0
        else:
            # Here we return the distance restraint target, the distance restraint weight,
            # the angle restraint target, the angle restraint weight, the jacks restraint target distance,
            # the jacks restraint weight, and the number of jacks.
            # In the future, we can make jacks_dist an array, so the jacks can have different target values.
            # Hopefully len(jacks_list) / 2 is an integer, otherwise the error analysis code will have trouble looping.
            return reference_values[7], dwt, 180.0, awt, jacks_dist, jwt, len(jacks_list) / 2

def find_index(usr_input, pdbfile):
    """
    Find the atom number in a PDB file.
    """
    atom = usr_input.split('@')[1]
    residue = usr_input.split('@')[0][1:]

    flag = 0 

    pdb_file = open(pdbfile, 'r')
     
    for line in pdb_file:
        if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM': 
            if (line[17:20].strip() == residue or line[22:26].strip() == residue) and line[12:16].strip() == atom:
                flag = 1
                break
    pdb_file.close()

    if not flag :
        print ('%s cannot be found.'%(usr_input))
        sys.exit()    
    return line[6:11].strip()
