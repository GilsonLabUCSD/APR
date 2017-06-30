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
  print '***********************************************************************'
  print '                      Welcome to Zalign!\n'
  print 'Options:' 
  print '    -f: location and name of the PDB file\n'
  print '    -a1, -a2, -a3: indices of three non-linear atoms (could be in either' 
  print '     receptor or ligand) that you think should lie in the same XY plane'
  print '     after the Z-axis alignment.\n'
  print '    -a4: a fourth atom that is supposed to have more negative Z'
  print '     coordinates than the first three, so that the molecule will'
  print '     point to the desired +Z or -Z direction.\n'
  print '    -c: the new origin of the aligned structure. The selected atom'
  print '     should be the same as the L1 option in the APR input file.\n'
  print '    -r: the rotation angle (in degrees) around the z axis.\n'
  print '    -cutoff: cutoff value for generating the final structure. Using'
  print '     a larger cutoff value (in Angstrom) will generate alignment with'
  print '     higher precision, but will of course takes more time to search.\n'
  print 'For example: python zalign.py -f pdb/oa_cba.pdb -a1 :OCT@O3 -a2 :OCT@O7'
  print '             -a3 :OCT@O8 -a4 :OCT@O11 -r 85 -cutoff 0.05 -c :MOL@N1'
  print '***********************************************************************'
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

def find_index(usr_input, PDBEntries):
    """
    Find the atom index of an atom in a PDB file.
    :param name: atom name (in PDB file)
    :param resname: atom resname (in PDB file)
    :return:
    """
    atom = usr_input.split('@')[1]
    residue = usr_input.split('@')[0][1:]

    flag = 0 
    count = 0

    for i in range(0, len(PDBEntries)):
        if PDBEntries[i][0:6].strip() == 'ATOM' or PDBEntries[i][0:6].strip() == 'HETATM':
            count += 1
            if (PDBEntries[i][17:20].strip() == residue or PDBEntries[i][22:26].strip() == residue) and PDBEntries[i][12:16].strip() == atom:
                flag = 1
                break

    if not flag :
        print ('%s cannot be found.'%(usr_input))
        sys.exit()    

    return count-1


# Check if a file exists
def checkFile(pdb_file):
    try:
      f = open(pdb_file,'r');
    except IOError:
      print 'Error: File does not exist'
      sys.exit()
    f.close()  

def readFile(pdb_file):
    content = []

    with open(pdb_file) as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line) # Non-blank lines in a list

    for i in range(len(lines)):
        splitdata = lines[i].split()
        # skip the header lines and seperating lines
        if (splitdata[0]=='ATOM')or(splitdata[0]=='HETATM')or(splitdata[0]=='TER'):
            content.append(lines[i])

    return content

# Release the memory of a list
def release_list(a):
    del a[:]
    del a

################################ zalign starts ##########################################
                          
if (len(sys.argv)!=17):
    align_exception()
    print 'Aborted. Check to see if you missed any flags.\n'
    sys.exit(1)

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

checkFile(pdb_file)
PDBatomEntries = readFile(pdb_file)

idx1 = find_index(R1,PDBatomEntries)
idx2 = find_index(R2,PDBatomEntries)
idx3 = find_index(R3,PDBatomEntries)
idx4 = find_index(R4,PDBatomEntries)
idx_origin = find_index(L1,PDBatomEntries)

# read coordinates and other information from the PDB file

total_atom  = 0
coords = []
ter_list = []
cols_before_resNumber = []
aChar_list = []
cols_after_coords = [] # Including everything after the coordinate columns
ter_row = []
resNumber_list = []

newPDB_file = open('align_z.pdb', 'w')

total_atom = 0 # not including "TER" entries

resid = 1

for i in range(len(PDBatomEntries)):
    if (PDBatomEntries[i][0:6].strip()=='ATOM')or(PDBatomEntries[i][0:6].strip()=='HETATM'):
        cols_before_resNumber.append(PDBatomEntries[i][0:22])
        total_atom += 1
        coords.append((float(PDBatomEntries[i][30:38].strip()), float(PDBatomEntries[i][38:46].strip()), float(PDBatomEntries[i][46:54].strip())))
        if i > 0:
            if PDBatomEntries[i-1][22:26]!=PDBatomEntries[i][22:26]:
                resid += 1
        resNumber_list.append(resid)
        aChar_list.append(PDBatomEntries[i][26:30])  
        cols_after_coords.append(PDBatomEntries[i][54:81])

    elif PDBatomEntries[i][0:6].strip() == 'TER':
        ter_list.append(total_atom)
        

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
        average = (float(coords_new2[idx1][2]) + float(coords_new2[idx2][2]) + float(coords_new2[idx3][2]))/3.0

        diff1 = float(coords_new2[idx1][2])-average
        diff2 = float(coords_new2[idx2][2])-average
        diff3 = float(coords_new2[idx3][2])-average
        if (abs(diff1) < cutoff)and(abs(diff2)<cutoff)and(abs(diff3)<cutoff)and(float(coords_new2[idx4][2])<average):
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
    tmpx = float(coords_new2[i][0]) - float(coords_new2[idx_origin][0])
    tmpy = float(coords_new2[i][1]) - float(coords_new2[idx_origin][1])
    tmpz = float(coords_new2[i][2]) - float(coords_new2[idx_origin][2])
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
        newPDB_file.write('%s%4d%s%8.3f%8.3f%8.3f%s\n'%(cols_before_resNumber[i], resNumber_list[i], aChar_list[i], coords_new3[i][0],\
                                                           coords_new3[i][1], coords_new3[i][2],cols_after_coords[i]))
        if i+1 in ter_list:
            newPDB_file.write('TER\n')
            j+=1 

    if total_atom not in ter_list:
       newPDB_file.write('TER\n')

    print 'The coordinates of three dummy atoms were appended in the end.'
    pb_resid = resNumber_list[i] + 1
    pb_atom_id = int(cols_before_resNumber[i][6:11].strip()) + 1    
 
    newPDB_file.write('%s%5s%4s%2s%3s%2s%4s%s\n'%('ATOM  ',pb_atom_id, 'Pb','','DUM', '', pb_resid,'       0.000   0.000  -6.000  1.00  0.00')) 
    newPDB_file.write('TER\n')
    newPDB_file.write('%s%5s%4s%2s%3s%2s%4s%s\n'%('ATOM  ',pb_atom_id + 1, 'Pb','','DUM', '', pb_resid + 1,'       0.000   0.000 -11.000  1.00  0.00'))
    newPDB_file.write('TER\n')
    newPDB_file.write('%s%5s%4s%2s%3s%2s%4s%s\n'%('ATOM  ',pb_atom_id + 2, 'Pb','','DUM', '', pb_resid + 2,'       0.000   3.500 -14.500  1.00  0.00'))
    newPDB_file.write('TER\n')
    newPDB_file.write('END\n')
    print 'A new pdb file, align_z.pdb has been generated.'

newPDB_file.close() 


