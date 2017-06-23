#!/usr/bin/env python2

#*********************
#      zalign.py      
# written by Jane Yin 
#*********************

import sys
import os
import math
import string
import re

def align_exception():
  print ''
  print '*****************************************************************************************************************************'
  print '                                             Welcome to Zalign!\n'
  print 'Options:' 
  print '    -f: location and name of the PDB file\n'
  print '    -a1, -a2, -a3: indices of three non-linear atoms (receptor or ligand) that you think should' 
  print '     lie in the same XY plane after the Z-axis alignment.\n'
  print '    -a4: a fourth atom that is supposed to have more negative Z coordinates than the first three,'
  print '     so that the molecule will point to the desired +Z or -Z direction.\n'
  print '    -c: the new origin of the aligned structure. The selected atom should be the same as the L1'
  print '        option in the APR input file.\n'
  print '    -r: the rotation angle (in degrees) around the z axis.\n'
  print '    -cutoff: cutoff value for generating the final structure. Using a larger cutoff value (in Angstrom)'
  print '     will generate alignment with higer precision, but will of course cost more time for searching.\n'
  print 'For example: python zalign.py -f pdb/oa_cba.pdb -a1 :OCT@O3 -a2 :OCT@O7 -a3 :OCT@O8 -a4 :OCT@O11 -r 85 -cutoff 0.05 -c :MOL@N1'
  print '             python zalign.py -f pdb/hsa-ligand.pdb -a1 :545@NZ -a2 :406@CB -a3 :541@NZ -a4 :573@NZ -cutoff 0.1 -r 0 -c :MYR@C8'
  print '             python zalign.py -f pdb/b-hex-s.pdb -a1 :6@O1 -a2 :8@O1 -a3 :10@O1 -a4 :8@H62 -cutoff 0.05 -r 0 -c :HEX@C6'
  print '*******************************************************************************************************************************'
  print '\n'
  return

# Check if the user input is an integer, a float number or a string
def ismyinstance(val_type, val, str):
  # if no input
  if not val: 
        return 0
  if val_type == 'float': 
    try:
        float(val)
    except ValueError:
        print  val
        print ('Please enter a float value for %s.'%(str))
        sys.exit()
  elif val_type == 'int':
    try:
        int(val)
    except ValueError:
        return False
        sys.exit()
  elif (val_type == 'string'):
        return val
  
  if (float(val) < 0):
        print  val
        print ('Please enter a non-negative value for %s.'%(str))
        sys.exit()
   
  if (val_type == 'float'):
        return float(val)
  elif (val_type == 'int'):
        return int(val) 

def find_index(usr_input, pdbfile):
    """
    Find the atom index of an atom in a PDB file.
    :param name: atom name (in PDB file)
    :param resname: atom resname (in PDB file)
    :return:
    """
    atom = usr_input.split('@')[1]
    residue = usr_input.split('@')[0][1:]
    resname  = 'None'
    res_id = 99999

    if ismyinstance('int',residue,'N/A'):
        # residue number  was provided
        res_id = int(residue)
    else:
        resname = residue

    flag = 0 

    with open(pdbfile) as pdb_file:
        lines = (line.rstrip('\n') for line in pdb_file)
        lines = list(line for line in pdb_file) 
 
    for i in range(0, len(lines)):
        newline = lines[i].split()
        if newline[0] == 'ATOM' or newline[0] == 'HETATM':
            if (newline[3] == resname or int(lines[i][22:27].strip())==res_id) and newline[2] == atom:
                flag = 1
                break
    pdb_file.close()

    if not flag :
        print ('%s cannot be found.'%(usr_input))
        sys.exit()    

    return int(newline[1])


# Check if a file exists
def check_file(pdb_file):
    try:
      f = open(pdb_file,'r');
    except IOError:
      print 'Error: File does not exist'
      sys.exit()
    f.close()  

# Release the memory of a list
def release_list(a):
    del a[:]
    del a

################################
#      zalign starts          #
                          
if (len(sys.argv)!=17):
    align_exception()
    print 'Aborted. Check to see if you missed any flags.'
    sys.exit()

#read arguments from the command line

arg_list = [1,3,5,7,9,11,13,15]
for i in arg_list:
    if sys.argv[i] == '-f':
      pdb_file = sys.argv[i+1]
    elif sys.argv[i] == '-a1':
      R1 = sys.argv[i+1]
    elif sys.argv[i] == '-a2':
      R2 = sys.argv[i+1]
    elif sys.argv[i] == '-a3':
      R3 = sys.argv[i+1]
    elif sys.argv[i] == '-a4':
      R4 = sys.argv[i+1]
    elif sys.argv[i] == '-c':
      L1 = sys.argv[i+1]
    elif sys.argv[i] == '-cutoff':
      cutoff = ismyinstance('float', sys.argv[i+1], sys.argv[i])
    elif (sys.argv[i] == '-r'):
      rotangle = ismyinstance('float', sys.argv[i+1], sys.argv[i])  
    else:
      print 'Wrong flags!! Please only use -f, -a1, -a2, -a3, -a4, -c, -r or -cutoff.'  

check_file(pdb_file)

idx1 = find_index(R1,pdb_file)
idx2 = find_index(R2,pdb_file)
idx3 = find_index(R3,pdb_file)
idx4 = find_index(R4,pdb_file)
idx_origin = find_index(L1,pdb_file)

# read coordinates and other information from the PDB file

total_atom  = 0
coords = []
ter_list = []

cols_before_coords = []
cols_after_coords = []
ter_row = []

# Keep track of the residue numbers and atom numbers
resid = []
atom_id = []

newPDB_file = open('align_z.pdb', 'w')

with open(pdb_file) as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

for i in range(len(lines)):
    splitdata = lines[i].split()
    # skip the header lines and seperating lines  
    if (splitdata[0]=='ATOM')or(splitdata[0]=='HETATM'):
        total_atom += 1
        coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
        resid.append(int(lines[i][22:26]))
        atom_id.append(int(lines[i][6:11]))
        cols_before_coords.append(lines[i][0:30])
        cols_after_coords.append(lines[i][54:81])
 
    elif splitdata[0] == 'TER':
        ter_list.append(total_atom)
        ter_row.append(lines[i])

flag = 0
print 'Start searching'
# rotate around the X-axis, 1 degree each time 
for i in range (360):
    coords_new = []
    dx = i*math.pi/180.0
    for num1 in range (total_atom):
        tmp1 = coords[num1][0]
        tmp2 = math.cos(dx)*coords[num1][1] + (-1)*math.sin(dx)*coords[num1][2]
        tmp3 = math.sin(dx)*coords[num1][1] + math.cos(dx)*coords[num1][2]
        coords_new.append((tmp1,tmp2,tmp3))
     
    #rotate around the Y-axis, 1 degree each time
 
    for j in range (360):
        coords_new2 =[]
        sys.stdout.write('\rScanning... '+ str((i*360+j)/(36*36)) + '%')
        dy = j*math.pi/180.0
        for num2 in range (total_atom):
            tmp1 = math.cos(dy)*coords_new[num2][0] + math.sin(dy)*coords_new[num2][2] 
            tmp2 = coords_new[num2][1]
            tmp3 = (-1)*math.sin(dy)*coords_new[num2][0] + math.cos(dy)*coords_new[num2][2]
            coords_new2.append((tmp1,tmp2,tmp3))
        average = (float(coords_new2[idx1-1][2]) + float(coords_new2[idx2-1][2]) + float(coords_new2[idx3-1][2]))/3.0

        diff1 = float(coords_new2[idx1-1][2])-average
        diff2 = float(coords_new2[idx2-1][2])-average
        diff3 = float(coords_new2[idx3-1][2])-average
        if (abs(diff1) < cutoff)and(abs(diff2)<cutoff)and(abs(diff3)<cutoff)and(float(coords_new2[idx4-1][2])<average):
            flag = 1
            print '\n'
            print 'Solution found.'
            break     #break the inner loop 
  
        release_list(coords_new2)

    if (flag==1):
        break   #break the outer loop 
    release_list(coords_new)

if flag == 0:
     print '\nCannot find a solution. Please try a larger cutoff value.\n'
     sys.exit(0) 

# Translate the coordinates according to the new origin
coords_new = []
for i in range(0, total_atom):
    tmpx = float(coords_new2[i][0]) - float(coords_new2[idx_origin-1][0])
    tmpy = float(coords_new2[i][1]) - float(coords_new2[idx_origin-1][1])
    tmpz = float(coords_new2[i][2]) - float(coords_new2[idx_origin-1][2])
    coords_new.append((tmpx, tmpy, tmpz))

# Rotate around the z axis
coords_new3 = []
dz = rotangle*math.pi/180.0

for i in range(total_atom):
    tmpx = math.cos(dz)*float(coords_new[i][0]) + (-1)*math.sin(dz)*float(coords_new[i][1]);
    tmpy = math.sin(dz)*float(coords_new[i][0]) + math.cos(dz)*float(coords_new[i][1]);  
    tmpz = float(coords_new[i][2]); #  z coordinate doesn't change
    coords_new3.append((tmpx, tmpy, tmpz))

j = 0

# Write the new pdb file
if (flag == 1):
    for i in range(len(coords)):
        newPDB_file.write('%s%8.3f%8.3f%8.3f%s\n'%(cols_before_coords[i], coords_new3[i][0], coords_new3[i][1], coords_new3[i][2],cols_after_coords[i]))
        if i+1 in ter_list:
            newPDB_file.write(ter_row[j]+'\n')
            j+=1 

    if total_atom not in ter_list:
       newPDB_file.write('TER\n')

    print 'The coordinates of three dummy atoms were appended in the end.'
    pb_resid = max(resid) + 1
    pb_atom_id = max(atom_id) + 1    
 
    newPDB_file.write('%s%5s%4s%2s%3s%2s%4s%s\n'%('ATOM  ',pb_atom_id, 'Pb','','DUM', '', pb_resid,'       0.000   0.000  -6.000  1.00  0.00')) 
    newPDB_file.write('TER\n')
    newPDB_file.write('%s%5s%4s%2s%3s%2s%4s%s\n'%('ATOM  ',pb_atom_id + 1, 'Pb','','DUM', '', pb_resid + 1,'       0.000   0.000 -11.000  1.00  0.00'))
    newPDB_file.write('TER\n')
    newPDB_file.write('%s%5s%4s%2s%3s%2s%4s%s\n'%('ATOM  ',pb_atom_id + 2, 'Pb','','DUM', '', pb_resid + 2,'       0.000   3.500 -14.500  1.00  0.00'))
    newPDB_file.write('TER\n')
    newPDB_file.write('END\n')
    print 'A new pdb file, align_z.pdb has been generated.'

f_in.close()
newPDB_file.close() 


