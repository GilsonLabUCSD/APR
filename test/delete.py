import os
import sys
import glob as glob
import subprocess as sp
import shutil as shutil


dirlist = os.listdir('windows')

os.chdir('windows')

for i in range(0, len(dirlist)):
    path = dirlist[i]
    for file in glob.glob('%s/*' % (path) ):
        if (file != '%s/restraints.dat' % (path)) and (file != '%s/reference.dat' % (path)):
            os.remove(file)

# Write mock files
for i in range(0, len(dirlist)):
    file = dirlist[i] + '/eqnpt50.rst7'
    f = open(file, 'w')
    f.write('This is a mock file; it exists only for the testing purpose') 
    f.close()  
