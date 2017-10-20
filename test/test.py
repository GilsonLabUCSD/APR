import os
import sys
import shutil
import glob as glob
import subprocess as sp

# This script tests the error_estimate module with a set of pre-generated restraints.dat files.

targetValUpper = -5.70
targetValLower = -5.71

# Copy over the scripts to the test folder
print('Copying over scripts ...')
for file in glob.glob('../*.py'):
    if file != '../test.py':   # Do not copy file that has the same name with the current script!
        shutil.copy2(file, './')

for file in glob.glob('../apr.in'):
    shutil.copy2(file, './')

for file in glob.glob('test.out'):
    os.remove(file)

print ('Testing apr analysis ... Sit tight :)')
sp.call('python2 apr.py analysis -i apr.in > test.out', shell = True)

# Grab the final result of binding free energy from the output file
with open('test.out') as fr:
    lines = fr.read().splitlines()
    lines = list(line for line in lines if line)

for i, line in enumerate(lines):
    line = lines[i].split()
    if line[0] == 'Binding':
        val = float(line[3])
        break

if (targetValLower < val < targetValUpper):
    print ('Test passed!')

# Clean up
for file in glob.glob('*.py'):
    if file != 'test.py':   # Do not delete itself!
        os.remove(file)

for file in glob.glob('*.pyc'):
    os.remove(file)

for file in glob.glob('apr.in'):
    os.remove(file)

for file in glob.glob('*.dat'):
    os.remove(file)
