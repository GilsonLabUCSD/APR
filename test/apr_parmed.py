#!/usr/bin/env python2
import sys
import os
from parmed.tools import changeLJSingleType
from parmed.tools import HMassRepartition
from parmed.tools import parmout
from parmed.tools import change
from parmed.amber import AmberParm
from parmed.structure import Structure

def perturb_parameters(perturb, hmr, prmtop):
    """
    Use ParmEd to make parameter perturbations
    """
    new_param = prmtop
    parm = AmberParm(prmtop+'.prmtop')
    logfile = open('parmed.log','w')

    if perturb == 'yes':
        # Parsing the input file that includes new parameters
        with open('new_params.dat') as fr:
            lines = fr.read().splitlines()
            lines = list(line for line in lines if line) # Remove blank lines

        for line in lines:
            if not line[0] == ';':
              splitline = line.split()
              if len(splitline) < 4:
                  print 'Aborted. Wrong format! Please provide atom type(e.g.,OW), parameter type(e.g.,vdw), and the new parameters\n '
                  sys.exit(1)
              atom_type = splitline[0]
              param_type = splitline[1]
              if len(splitline) == 4:
                  radius = splitline[2]
                  epsilon = splitline[3]

              if param_type.lower() == 'vdw':
                  action = changeLJSingleType(parm, "@%"+atom_type, radius, epsilon)
                  action.execute()
                  logfile.write(('%s\n' % action))
              else:
                  print 'Aborted. Only the feauture of perturbing non-bonded parameters is supported for now. (keyword: vdw)\n' 
                  sys.exit(1)
        new_param += '.perturbed'
        if os.path.isfile(new_param + '.prmtop'):
            os.remove(new_param + '.prmtop')
        Structure.save(parm, new_param + '.prmtop') 

    if hmr == 'yes':
        parm = AmberParm(new_param + '.prmtop')
        # Use the ParmEd API, which is more stable than calling a subprocess
        action = HMassRepartition(parm)
        action.execute()
        logfile.write(('%s\n' % action))
        new_param += '.hmr'
        if os.path.isfile(new_param + '.prmtop'):
       	    os.remove(new_param + '.prmtop')
        Structure.save(parm, new_param + '.prmtop')

    return new_param
     

