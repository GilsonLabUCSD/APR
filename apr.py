#!/usr/bin/env python2
import os as os
import sys as sys
import glob as glob
import shutil as shutil
import datetime as dt
import signal as signal
import subprocess as sp
import re as re
import numpy as np

import apr_mdin as apr_mdin
import apr_solvate as apr_solvate
import apr_translate as apr_translate
import apr_restraints as apr_restraints
import apr_parmed as apr_parmed
import err_estimate as err_estimate

class APR:
    def __init__(self):
        """
        Default user and system parameters specification. These are now overwritten with the values found in apr.in.
        More information on these variables can be found in the input file and the online tutorial.
        """
        self.amber16 = 'yes'
        self.action1 = 'setup'
        self.action2 = 'continue'
        self.input_file = 'apr.in'
        self.exe_path = 'pmemd.cuda'
        self.perturb = 'no'           # Whether to perturb the GAFF parameters (via ParmEd)
        self.prmtop = 'solvated'      # The intial filename of Amber topology (+.prmtop)

        # System parameters
        # These attachment restraints are scaled from 0 to 1 in the setup_restraints() function.
        self.attach_fc = []
        # Distance targets for pulling phase in Angstroms
        self.trans_dist = []
        self.md_temperature = 298.15       # Kelvin
        self.analysis_temperature = 298.15 

        # Restraint force constants
        self.dist_fc = 5.0     # kcal/mol/Angstrom**2
        self.angle_fc = 100.0  # kcal/mol/rad**2

        # Type and number of counterions in system,
        # including ions for neutralization only,
        # and extra ions to mimic the salt buffer
        self.ions = []

        # Target number of waters (plus ions) in each window
        self.waters = 5000
        self.warning = 'on'   # Provide an estimation for number of water molecules needed for full solvation

        # Total number of atoms in solute (dummy atoms included)
        self.solute_atoms = 0 

        # Water model for solvation
        self.water_model = 'TIP3P'

        # Receptor and ligand atoms for the restraints (accept AMBER style mask)
        self.lig_name = 'None'
        self.dum_resid = 99999
        self.R1 = ''
        self.R2 = ''
        self.R3 = ''
        self.L1 = ''
        self.L2 = ''
        
        # Attributes for the production phase
        self.hmr = 'no'
        self.stepsize = 2 # unit: fs
        self.steps = 1250000  
        self.ntpr = 500  # energy output frequency
        self.ntwx = 500  # trajectory output frequency          
        self.barostat = 2   # Monte Carlo barostat
        self.cutoff = 9     # vdW cutoff in angstrom
        self.strip = 'yes'  # strip water molecules and ions in the MD trajectories (to save disc space)
 
        # Attributes for the equilibration phase
       	self.eq_stepsize = 2 # unit: fs
       	self.eq_steps = 2500
       	self.eq_barostat = 2   # Monte Carlo barostat
       	self.eq_cutoff = 9	    # vdW cutoff in angstrom

        # The maximum number of iterations and error threshold, for the production phase only
        self.maxcycle = 10
        self.maxsem_a = 1.2
        self.maxsem_u = 0.1
        self.maxsem_r = 0.1

        self.jacks = 'no'

        # Force constant of conformational restraints and jacks input file; Not used if jacks = no
        self.jacks_dist = 15  # Angstrom
        self.jacks_fc = 25.0  # kcal/mol/Angstrom**2
        self.jacks_list = []

    def now(self):
        """
        Return the date (for the the log file).
        :return:
        """
        d_date = dt.datetime.now()
        return d_date.strftime("%Y-%m-%d %I:%M:%S %p")

    def check_arguments(self):
        """
        Ensure the user calls this script with -overwrite or -continue to avoid overwriting existing
        simulation data in the windows/ folder unintentionally.
        :return:
        """
        if len(sys.argv) == 1:
            help_message()
            # Return without an error
            sys.exit(0)
        elif (len(sys.argv) == 2) and ('-h' in sys.argv[1].lower()):
            help_message()
            # Return without an error
            sys.exit(0)
        elif (len(sys.argv) == 6) and ('eq' == sys.argv[1].lower()):
            print('Setting up the APR framework ...')
            self.action1 = 'setup'
        elif (len(sys.argv) == 6) and ('prod' == sys.argv[1].lower()):
            print('Running the production phase ...')
            self.action1 = 'production'
        elif (len(sys.argv) == 4) and ('analysis' == sys.argv[1].lower()):
            print('Analyzing the binding calculations ...')
            self.action1 = 'analysis'
        elif (len(sys.argv) > 4) and ('analysis' == sys.argv[1].lower()):
            help_message()
            print('Analysis does not need -s flags. Aborted.')
            sys.exit(1)
        else:
            help_message()
            print('Choose among eq, prod and analysis.')
            print('Flags -s and -i needed for eq(equilibration) and prod(production). Flag -i needed for analysis. Aborted!')
            sys.exit(1)

        if self.action1 != 'analysis':
            # Loop over the list of command line arguments and look for the `-s` and `-i` flags
            for i in [2, 4]:
                if '-s' == sys.argv[i].lower():
                    if 'overwrite' == sys.argv[i + 1]:
                        print('Overwrite mode enabled. Start freshly ...')
                        self.action2 = 'overwrite'
                    elif 'continue' == sys.argv[i + 1]:
                        print('Continue mode enabled. Checking previous simulations ...')
                        self.action2 = 'continue'
                    else:
                        help_message()
                        print('Please choose between -continue and -overwrite. Aborted!')
                        sys.exit(1)
                elif '-i' == sys.argv[i].lower():
                    self.input_file = sys.argv[i + 1]
                    # check the input file
                    if not os.path.isfile(self.input_file):
                        print('The input file %s does not exist.' % self.input_file)
                        sys.exit(1)
                else:
                    help_message()
                    print('I could not find a -s or -i flag in your input.')
                    sys.exit(1)

    def process_input_file(self):
        """
        Reading the user input file
        """
        with open(self.input_file) as f_in:
        # Delete all spaces/tabs at both ends
            lines = (line.strip(' \t\n\r') for line in f_in)
            lines = list(line for line in lines if line)  # Non-blank lines in a list

        for i in range(0, len(lines)):
            # Combine the lines that belong to the same entry
            if not lines[i][0] == ';':
                lines[i] = lines[i].split(';')[0].split('=')
                if len(lines[i]) == 1:
                    j = i
                    while True:
                        if lines[j - 1] != ';':
                            lines[j - 1][1] += lines[i][0]
                            lines[i] = ';'
                            break
                        j -= 1
               
        for i in range(0, len(lines)):
            if not lines[i][0] == ';':
                lines[i][0] = lines[i][0].strip().lower()
                lines[i][1] = lines[i][1].strip()
                if lines[i][0] == 'amber16':
                    self.amber16 = read_bool_attributes(lines[i][1],lines[i][0])
                elif lines[i][0] == 'hmr':
                    self.hmr = read_bool_attributes(lines[i][1],lines[i][0])
                elif lines[i][0] == 'perturb':
                    self.perturb = read_bool_attributes(lines[i][1],lines[i][0])
                elif lines[i][0] == 'strip_water_ions':
                    self.strip = read_bool_attributes(lines[i][1],lines[i][0])
                elif lines[i][0] == 'temperature':
                    self.analysis_temperature = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                    self.md_temp  = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'exe_path':
                    self.exe_path = lines[i][1].strip('\'\"-,.:;][')
                elif lines[i][0] == 'distance_force':
                    self.dist_fc = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'angle_force':
                    self.angle_fc = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'water_model':
                    # A known water model is not required here to enhance flexibility.
                    # If you want to use your own water model, make sure the water topology files are available and
                    # Amber can read them in correctly.
                    self.water_model = lines[i][1].upper()
                elif lines[i][0] == 'warning':
                    self.warning = read_bool_attributes(lines[i][1],lines[i][0])
                elif lines[i][0] == 'neutralizing_cation':
                    neu_cation = ismyinstance('string', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'neutralizing_anion':
                    neu_anion = ismyinstance('string', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'cations':
                    extra_cation = ismyinstance('string', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'number_cations':
                    num_cations = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'anions':
                    extra_anion = ismyinstance('string', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'number_anions':
                    num_anions = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'waters':
                    self.waters = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'attach_list':
                    newline = lines[i][1].strip('\'\"-,.:;][').split(',')
                    for j in range(0, len(newline)):
                        self.attach_fc.append(ismyinstance('float', newline[j], self.input_file, lines[i][0]))
                elif lines[i][0] == 'translate_list':
                    newline = lines[i][1].strip('\'\"-,.:;][').split(',')
                    for j in range(0, len(newline)):
                        self.trans_dist.append(ismyinstance('float', newline[j], self.input_file, lines[i][0]))
                elif lines[i][0] == 'lig':
                    self.lig_name = ismyinstance('string', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'r1':
                    self.R1 = lines[i][1]
                elif lines[i][0] == 'r2':
                    self.R2 = lines[i][1]
                elif lines[i][0] == 'r3':
                    self.R3 = lines[i][1]
                elif lines[i][0] == 'l1':
                    self.L1 = lines[i][1]
                elif lines[i][0] == 'l2':
                    self.L2 = lines[i][1]
                elif lines[i][0] == 'maxcycle':
                    self.maxcycle = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'maxsem_attach':
                    self.maxsem_a = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'maxsem_pull':
                    self.maxsem_u = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'maxsem_release':
                    self.maxsem_r = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'jacks':
                    self.jacks = read_bool_attributes(lines[i][1],lines[i][0])
                elif lines[i][0] == 'jacks_distance':
                    self.jacks_dist = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'jacks_force':
                    self.jacks_fc = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'jacks_list':
                    newline = lines[i][1].strip('\'\"-,.;][').split(',')
                    for j in range(0, len(newline)):
                        self.jacks_list.append(newline[j].strip())
                elif lines[i][0] == 'nstlim':
                    self.steps = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'eq_nstlim':
                    self.eq_steps = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'ntpr':
                    self.output_freq = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'ntwx':
                    self.trajout_freq = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'dt':
                    self.stepsize = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'eq_dt': 
                    self.eq_stepsize = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'cutoff':
                    self.cutoff = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'eq_cutoff':
                    self.eq_cutoff = ismyinstance('float', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'barostat':
                    self.barostat = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                    if self.barostat != 1 and self.barostat!=2 :
                       print ('Wrong input! The value for barostat is either 1 or 2.')
                       sys.exit(1)
                elif lines[i][0] == 'eq_barostat':
                    self.eq_barostat = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                    if self.eq_barostat != 1 and self.eq_barostat!=2 :
                       print ('Wrong input! The value for barostat is either 1 or 2.')
                       sys.exit(1)
                else:
                    print('Wrong entry name in %s! ' % self.input_file)
                    print(lines[i][0] + '\n')
                    print('Please use the same keywords as in the template input file. Aborted.')
                    sys.exit(1)

        # Update self.ions list
        self.ions = [[neu_cation, 0], [neu_anion, 0], [extra_cation, num_cations], [extra_anion, num_anions]]

        if not os.path.exists('setup/align_z.pdb'):
            print('Aborted! align_z.pdb cannot be found in the setup folder. Please use zalign.py to generate it first.\n')
            sys.exit(1)            

        # Find the total number of solute atoms and the residue serial number of the first dummy atom:
        f =  open('setup/align_z.pdb','r') 
        for line in f:
            if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM':
                self.solute_atoms += 1
        f.seek(0,0)
        for line in f:
            if  line[17:20].strip() == 'DUM':
                self.dum_resid = int(line[22:26].strip())   
                break
        f.close()
      
        if self.dum_resid == 99999:
            print 'Dummy atoms were not detected in align_z.pdb.\n'

        if self.perturb == 'yes':
            if not os.path.exists('setup/param_files/new_params.dat'):
                print('Aborted! new_params.dat cannot be found in the setup/param_files folder.')
                print ('Please provide new parameters for perturbation, or switch the perturb option from yes to no in the APR input.\n')
                sys.exit(1)


    def make_files_and_directories(self, method):
        """
        Copy the structure and input files to each subdirectory for individual APR windows.
        :param method: attachment or translation (pulling)
        :return:
        """
        if method == 'attachment':
            prefix = 'a'
            windows = len(self.attach_fc)
        elif method == 'translation':
            prefix = 'p'
            windows = len(self.trans_dist)
        elif method == 'release':
            prefix = 'r'
            windows = len(self.attach_fc)
        else:
            print('Error making files and directories.')

        if not os.path.exists('windows'):
            os.makedirs('windows')

        for window in range(0, windows):
            destination = 'windows/%s%02d' % (prefix, window)
            if not os.path.exists(destination):
                os.makedirs(destination)

            if (self.action2 == 'continue') and (os.path.isfile('%s/eqnpt50.rst7' % destination)):
                print('Skipping folder %s' % destination)

            else:
                for file in glob.glob('setup/param_files/*'):
                    shutil.copy(file, destination)
                #for file in glob.glob('setup/input_files/*'):
                    #shutil.copy(file, destination)
                shutil.copy('setup/align_z.pdb', destination)
            if self.amber16 == 'yes':
                shutil.copy('setup/input_files/tleap.in.amber16', '%s/tleap.in'%(destination))
            else:
                # Assuming Amber14 or Amber12
                shutil.copy('setup/input_files/tleap.in.amber14', '%s/tleap.in'%(destination))

    def check_executable(self):
        """
        Check user defined pmemd/sander settings and echo
        """
        flag = 0
        exe_str = self.exe_path.split()
        for i in range(0, len(exe_str)):
          if '-np' in exe_str[i]:
            flag = i
          if 'pmemd' in exe_str[i] or 'sander' in exe_str[i]:
            print ('%s is on.'%(exe_str[i]))
          else:
            print 'Neither pmemd nor sander is defined. Are you sure this is the right executable?'

        if flag < len(exe_str) and flag > 0:
          print ('Totally %d processors are used.'%(int(exe_str[flag+1])))

    def prepare_and_simulate(self, method):
        """
        Set the parameters for  translation, solvation, restraints, and thermalization of the
        host-guest complex.
        :param method: attachment or translation (pulling)
        :return:
        """

        if method == 'attachment':
            prefix = 'a'
            windows = len(self.attach_fc)
            values = self.attach_fc
        elif method == 'translation':
            prefix = 'p'
            windows = len(self.trans_dist)
            values = self.trans_dist
        elif method == 'release':
            prefix = 'r'
            windows = len(self.attach_fc)
        else:
            print('Error preparing.')
        
        for window in range(windows):
            destination = 'windows/%s%02d' % (prefix, window)

            if self.action2 == 'continue' and os.path.isfile('%s/eqnpt50.rst7' % (destination)):
                # Skip to next window
                continue

            # Change directory and call the setup functions
            os.chdir(destination)
            print('Preparing folder %s %s'%(destination, self.now()))

            sys.stdout.flush()

            # Calculate scale_w
            # Can we store the scale w somewhere?
            if method == 'attachment' or method == 'release':
                rest_weight = self.attach_fc[window]
                if window == 0:
                    if len(self.attach_fc) == 1:
                        scale_w = 0.001
                    else:
                        scale_w = 0.01 * (self.attach_fc[1]-self.attach_fc[0]) / 2
                elif window == windows - 1:
                    scale_w = 0.01 * (100 - self.attach_fc[window - 1]) / 2
                else:
                    scale_w = 0.01 * (self.attach_fc[window + 1] - self.attach_fc[window - 1]) / 2
                if method == 'attachment':
                    shutil.copyfile('align_z.pdb', 'vac.pdb')
                if method == 'release':
                    apr_translate.setup_translate(self.trans_dist[-1], self.lig_resid, self.lig_resname)

            # Set up windows to pull the guest away from the host
            if method == 'translation':
                rest_weight = 100.0
                scale_w = 1.0
                apr_translate.setup_translate(self.trans_dist[window], self.lig_resid, self.lig_resname)

            # Add the correct number of waters 
            apr_solvate.setup_solvate(self.warning, self.water_model, self.waters, self.ions, self.amber16)

            # Add restraints to keep the guest at the target distances
            if 'no' in self.jacks.lower():
                self.jacks_fc = 0

            if method == 'attachment':
                apr_restraints.setup_restraints(prefix, 0.0, rest_weight, scale_w, self.R1, self.R2, self.R3, self.L1, self.L2,
                                 self.dist_fc, self.angle_fc, self.jacks_fc, self.jacks_dist, self.jacks_list, self.strip)
            if method == 'translation':
                apr_restraints.setup_restraints(prefix, self.trans_dist[window], rest_weight, scale_w, self.R1,
                                 self.R2, self.R3, self.L1, self.L2, self.dist_fc,
                                 self.angle_fc, self.jacks_fc, self.jacks_dist, self.jacks_list, self.strip)
            if method == 'release':
                apr_restraints.setup_restraints(prefix, self.trans_dist[-1], rest_weight, scale_w, self.R1,
                                 self.R2, self.R3, self.L1, self.L2, self.dist_fc,
                                 self.angle_fc, self.jacks_fc, self.jacks_dist, self.jacks_list, self.strip)

            sys.stdout.flush()
           
            apr_mdin.write_min_in(self.dum_resid)
            apr_mdin.write_therm1_in(self.dum_resid)
            apr_mdin.write_therm2_in(self.md_temp, self.dum_resid)
            apr_mdin.write_eqnpt_in(self.md_temp, self.eq_stepsize, self.eq_steps, self.eq_barostat, self.eq_cutoff, self.dum_resid)
            if self.strip == 'yes':
                self.ntwprt = self.solute_atoms  # Only solute atoms will be included in the trajectory
            else:
                self.ntwprt = 0  # All atoms will be included in the trajectory
            apr_mdin.write_mdin(self.md_temp, self.stepsize, self.steps, self.barostat, self.cutoff, self.ntpr, self.ntwx, self.ntwprt, self.dum_resid)

            local_prmtop = apr_parmed.perturb_parameters(self.perturb, self.hmr, self.prmtop)
            local_prmtop += '.prmtop'

            for equilibration_counter in range(5):
                # Thermalize and equilibrate
                status = equilibrate(self.exe_path, local_prmtop)
                if status == 0:
                    break
                elif status == 1:
                    print('Equilibration failed in window %s%02d' % (prefix, window))
                    print('Re-running solvation step.')
                    apr_solvate.setup_solvate(self.warning, self.water_model, self.waters, self.ions, self.amber16)
                    equilibration_counter += 1

            if equilibration_counter == 5:
                print 'APR stopped at the solvation stage. Please check whether AMBER, PMEMD, CUDA, or MPI were corrected installed and called.'
                print 'Also check the solvation.log output file generated by tLeap to see if there are any warnings or error messages.'
                print 'Another possiblity is that the restraints were not set correctly or reasonably. Incorrect restraints will be reflected'
                print 'in the local disang.rest and restart (.rst7) files.\n'   
                sys.exit(1)

            sys.stdout.flush()

            # Return to the root directory
            os.chdir('../../')

    def checkfile(self, stage):
        """
        Check whether the equilibration phase has been finished in all windows.
        """
        # Check to see if windows folder exists
        if not os.path.exists('windows'):
            sys.stderr.write('Aborted! The windows folder does not exist.')
            sys.stderr.write('Please make sure the equilibration phase has been performed.\n')
            sys.exit(1)

        # Check eqnpt50.rst7 in all windows
        for i in ['a', 'p', 'r']:
            if i == 'a':
                windows = len(self.attach_fc)
            elif i == 'p':
                windows = len(self.trans_dist)
            elif i == 'r':
                if self.jacks == 'no':
                    return
                else:
                    windows = len(self.attach_fc)
            for window in range(windows):
                if not os.path.isfile('windows/%s%02d/eqnpt50.rst7' % (i, window)):
                    sys.stderr.write(
                        'Aborted! Either the %s%02d window does not exist, or the eqnpt50.rst7 '
                        'file is missing in this window.\n' % (i, window))
                    sys.stderr.write('Please make sure the equilibration phase has been completed.\n')
                    sys.exit(1)
             
                # Check restraints.dat in all windows 
                if stage == 'analysis':
                   if not os.path.isfile('windows/%s%02d/restraints.dat' % (i, window)):
                      sys.stderr.write(
                          'Aborted! File: restraints.dat'
                          ' is missing in the %s%02d window.\n' % (i, window))
                      sys.stderr.write('Please make sure the production phase has been completed.\n')
                      sys.exit(1)
     

    def simulate(self, method):
        """
        Run pmemd for the production phase. 
        """
        if method == 'attachment':
            prefix = 'a'
            windows = len(self.attach_fc)
            maxsem = self.maxsem_a
        elif method == 'translation':
            prefix = 'p'
            windows = len(self.trans_dist)
            maxsem = self.maxsem_u
        elif method == 'release':
            prefix = 'r'
            windows = len(self.attach_fc)
            maxsem = self.maxsem_r

        for window in range(windows):
            destination = 'windows/%s%02d' % (prefix, window)
            # Change directory and run simulations in each window
            os.chdir(destination)

            print('Simulating folder {:<25} {:<10}'.format(destination, self.now()))
            
            mdin = 'mdin' 
            if self.hmr == 'yes' and self.perturb == 'yes':
                prmtop = self.prmtop + '.perturbed' + '.hmr' + '.prmtop'
            elif self.perturb == 'yes':
                prmtop = self.prmtop + '.perturbed' + '.prmtop'
            elif self.hmr == 'yes':
                prmtop = self.prmtop + '.hmr' + '.prmtop'
            else:
                prmtop = self.prmtop + '.prmtop'

            # Clean up restraints.dat and restraints.log for both continue and overwrite mode
            if os.path.isfile('restraints.dat'):
                sp.call('rm restraints.dat', shell = True)
            if os.path.isfile('restraints.log'):
                sp.call('rm restraints.log', shell = True)

            if self.action2 == 'continue':
                mdoutstep = 1
                while (len(sp.Popen('grep TIMINGS mdout.%02.0f | grep -o TIMINGS'%mdoutstep,
                           stdout=sp.PIPE, stderr=sp.PIPE, shell=True).stdout.read().splitlines()) == 1):
                    mdoutstep += 1
                currcycle = mdoutstep - 1
            else:
                if os.path.isfile('hist.dat'):
                    sp.call('rm hist.dat', shell = True)
                for i in ['traj', 'rst', 'mden', 'mdinfo', 'mdout']:
                    if glob.glob('%s.*'%(i)):
                        # Delete previous trajectories
                        sp.call('rm %s.*'%(i),shell = True)
                currcycle = 0

            if not os.path.isfile('rst.00'):
              shutil.copyfile('eqnpt50.rst7','rst.00')

            while True:
                # Get scale_w ... again!
                if method == 'attachment' or method == 'release':
                    rest_weight = self.attach_fc[window]
                    if window == 0:
                        if len(self.attach_fc) == 1:
                            scale_w = 0.001
                        else:
                            scale_w = 0.01 * (self.attach_fc[1]-self.attach_fc[0]) / 2
                    elif window == windows - 1:
                        scale_w = 0.01 * (100 - self.attach_fc[window - 1]) / 2
                    else:
                        scale_w = 0.01 * (self.attach_fc[window + 1] - self.attach_fc[window - 1]) / 2

                # Set up windows to pull the guest away from the host
                if method == 'translation':
                    rest_weight = 100.0
                    scale_w = 1.0

                if method == 'attachment':
                    restraint_list = apr_restraints.return_restraints_for_error_analysis(prefix, 0.0,
                                                                                         self.R1, self.R2, self.R3,
                                                                                         self.L1, self.L2, self.dist_fc,
                                                                                         self.angle_fc, self.jacks_fc,
                                                                                         self.jacks_dist, self.jacks_list,
                                                                                         self.jacks)
                if method == 'translation':
                    restraint_list = apr_restraints.return_restraints_for_error_analysis(prefix, self.trans_dist[window],
                                                                                         self.R1,
                                                    self.R2, self.R3, self.L1, self.L2, self.dist_fc,
                                                    self.angle_fc, self.jacks_fc, self.jacks_dist, self.jacks_list, self.jacks)
                if method == 'release':
                    restraint_list = apr_restraints.return_restraints_for_error_analysis(prefix, self.trans_dist[-1],
                                                                                         self.R1,
                                                    self.R2, self.R3, self.L1, self.L2, self.dist_fc,
                                                    self.angle_fc, self.jacks_fc, self.jacks_dist, self.jacks_list, self.jacks)

                # Generate restraints.dat file

                if os.path.isfile('traj.01') and len(sp.Popen('grep TIMINGS mdout.01 | grep -o TIMINGS',stdout=sp.PIPE,
                                                      stderr=sp.PIPE, shell=True).stdout.read().splitlines()) == 1:
                    
                    sp.call('cpptraj -i restraints.in >& restraints.log', shell=True)

                if os.path.isfile('restraints.dat'):
                    error_value = err_estimate.error_estimate(scale_w, method, restraint_list)
                    print ('SEM of the forces: %6.3f kcal/mol; Threshold(maxsem): %6.3f kcal/mol'%(error_value,maxsem))
                else:
                    error_value = maxsem + 1.0
                if (currcycle == self.maxcycle) or (error_value < maxsem):
                    break
                else:
                    currcycle += 1

                sys.stdout.flush()
                 
                excstr = '%s -O -p %s -ref solvated.inpcrd -c rst.%02d\
                   -i %s -o mdout.%02d -r rst.%02d -x traj.%02d -inf mdinfo.%02d -e mden.%02d' % (
                self.exe_path, prmtop, currcycle - 1, mdin, currcycle, currcycle, currcycle, currcycle, currcycle)

                p = sp.call(excstr, shell=True)
            os.chdir('../../')

    def WRelease(self,kr,r0,kt,t0,kb,b0,T):
        """
           Do the analytical "release" of guest to standard concentration
        """
        ### SETUP
        # Example: print WRelease(5.0,24.00,100.0,180.0,100.0,180.0,298.15)
        # kr, r0: distance restraint (r) force constant, target value
        # kt, t0: angle restraint (theta) force constant, target value
        # kb, b0: angle restraint (b) force constant, target value
        # T: Temperature
        R = 1.987204118e-3 # kcal/mol-K, a.k.a. boltzman constant
        beta = 1/(T*R)
        rlb,rub,rst = [0.0,100.0,0.0001]  # r lower bound, upper bound, step size
        tlb,tub,tst = [0.0,np.pi,0.00005] # theta ",       ",            "
        blb,bub,bst = [0.0,np.pi,0.00005] # b     ",       ",            "
        def fr(val):
          return (val**2)*np.exp(-beta*kr*(val-r0)**2)
        def ft(val):
          return np.sin(val)*np.exp(-beta*kt*(val-np.radians(t0))**2)
        def fb(val):
          return np.sin(val)*np.exp(-beta*kb*(val-np.radians(b0))**2)
        ### Integrate
        rint,tint,bint = [0.0,0.0,0.0]
        intrange = np.arange(rlb,rub,rst)
        rint = np.trapz(fr(intrange),intrange)
        intrange = np.arange(tlb,tub,tst)
        tint = np.trapz(ft(intrange),intrange)
        intrange = np.arange(blb,bub,bst)
        bint = np.trapz(fb(intrange),intrange)
        return R*T*np.log(np.pi*(1.0/1660.0)*rint*tint*bint)

    def ti_analysis(self,method):
        """
            Set windows and run TI analysis
        """        
        ### Determine number of windows for this phase
        windows = 0
        if method == 'attachment':
            windows = len(self.attach_fc)
            prefix = 'a'
            end_point = 100                  # lambda percentage
        elif method == 'release':
            windows = len(self.attach_fc)
            prefix = 'r'
            end_point = 100                   # lambda percentage
        elif method == 'translation':
            windows = len(self.trans_dist)
            prefix = 'p'

        ### Initialize Forces
        # First remember, that first and last windows of translation are also part of attachment and release, respectively.
        if method == 'attachment' or method == 'release':
          total_windows = windows + 1
        else:
          total_windows = windows
        FrcMeans = np.zeros([total_windows], np.float64)
        FrcSEMs = np.zeros([total_windows], np.float64)
        RxnCrds = np.zeros([total_windows], np.float64)
    
    
        for window in range(windows):
            destination = 'windows/%s%02d' % (prefix, window)
            os.chdir(destination)
    
            ### Copied this from above.  Need to get restraints info for get_forces!
            if method == 'attachment':
                restraint_list = apr_restraints.return_restraints_for_error_analysis(prefix, 0.0,
                                                                                     self.R1, self.R2, self.R3,
                                                                                     self.L1, self.L2, self.dist_fc,
                                                                                     self.angle_fc, self.jacks_fc,
                                                                                     self.jacks_dist, self.jacks_list,
                                                                                     self.jacks)
            if method == 'translation':
                restraint_list = apr_restraints.return_restraints_for_error_analysis(prefix, self.trans_dist[window],
                                                                                     self.R1,
                                                self.R2, self.R3, self.L1, self.L2, self.dist_fc,
                                                self.angle_fc, self.jacks_fc, self.jacks_dist, self.jacks_list, self.jacks)
            if method == 'release':
                restraint_list = apr_restraints.return_restraints_for_error_analysis(prefix, self.trans_dist[-1],
                                                                                          self.R1,
                                                     self.R2, self.R3, self.L1, self.L2, self.dist_fc,
                                                     self.angle_fc, self.jacks_fc, self.jacks_dist, self.jacks_list, self.jacks)


            if method == 'attachment' or method == 'release':
                frcvals = err_estimate.get_forces(method,restraint_list)
                FrcMeans[window] = frcvals[1]    #NMH: to match old scripts, use frcvals[3] ie, block mean.
                FrcSEMs[window] = frcvals[4]
                RxnCrdLoc = self.attach_fc[window] / end_point
                RxnCrds[window] = RxnCrdLoc
                if window == windows-1 and method == 'attachment':
                    # Do something special on the last window
                    os.chdir('../../')
                    # We also need the starting and endpoint
                    destination = 'windows/p00'
                    os.chdir(destination)
                    frcvals = err_estimate.get_forces(method,restraint_list)
                    FrcMeans[window + 1] = frcvals[1]
                    FrcSEMs[window + 1] = frcvals[4]
                    RxnCrdLoc = 1.0
                    RxnCrds[window + 1] = RxnCrdLoc
                if window == windows-1 and method == 'release':
                    # Do something special on the last window
                    os.chdir('../../')
                    # We also need the starting and endpoint
                    destination = 'windows/p%02d' % (windows - 1)
                    os.chdir(destination)
                    frcvals = err_estimate.get_forces(method,restraint_list)
                    FrcMeans[window + 1] = frcvals[1]
                    FrcSEMs[window + 1] = frcvals[4]
                    RxnCrdLoc = 1.0
                    RxnCrds[window + 1] = RxnCrdLoc
    
            else:
                frcvals = err_estimate.get_forces(method,restraint_list)
                FrcMeans[window] = frcvals[1]
                FrcSEMs[window] = frcvals[4]
                RxnCrdLoc = restraint_list[0]
                RxnCrds[window] = RxnCrdLoc

            os.chdir('../../')
    
        ### Prepare to Integrate
        BootCyc = 50000  # Num of Boot strap cycles
        if method == 'attachment' or method == 'release':
          windows += 1
        Intg = np.zeros([windows, BootCyc], np.float64)  # Integration value at each window for each boot cycle
        Xspl = np.zeros([0], np.float64)  # array for x dimension spline points
        Idx = np.zeros([windows], np.int32)  # index location for each window in the spline array
        Idx[0] = 0
        for i in range(1, windows):
            Xspl = np.append(Xspl, np.linspace(RxnCrds[i - 1], RxnCrds[i], num=100,
                                               endpoint=False))  # 100 spline poins between windows
            Idx[i] = len(Xspl)
        Xspl = np.append(Xspl, RxnCrds[-1])
    
        ### Integrate with bootstrapping
        Yboot = np.zeros([windows], np.float64)
        for c in range(BootCyc):
            for i in range(windows):
                Yboot[i] = np.random.normal(FrcMeans[i], FrcSEMs[i], 1)
            #print RxnCrds, Yboot, Xspl
            Yspl = err_estimate.interpolate(RxnCrds, Yboot, Xspl)
            for i in range(windows):
                if method == 'attachment' or method == 'release':  ### Attach/Release
                    Intg[i, c] = np.trapz(Yspl[0:Idx[i]], Xspl[0:Idx[i]])
                else:  # Umbrella Translate
                    Intg[i, c] = -1 * np.trapz(Yspl[0:Idx[i]], Xspl[0:Idx[i]])
    
        ### Print Integration
        ti_file = open('TI_%s.dat' % method, 'w')
        if method == 'attachment' or method == 'release':
            ti_file.write ('#Weight of restraints (%), Accumulative work (in kcal/mol), SEM (in kcal/mol)\n')
        elif method == 'translation':
            ti_file.write ('#Translation (pulling) distance (in Angstrom), Accumulative work (in kcal/mol), SEM (in kcal/mol)\n')
        for i in range(windows):
            ti_file.write("%6.4f %12.5f %12.5f\n"%(RxnCrds[i], np.mean(Intg[i]), np.std(Intg[i])))
        ti_file.close()

        # This assumes identical force constants and target equilibrium values for both angles
        if method == 'translation':
            wrelease = self.WRelease(float(restraint_list[1]),float(restraint_list[0]),float(restraint_list[3]),
                                float(restraint_list[2]),float(restraint_list[3]),float(restraint_list[2]), self.analysis_temperature)
            return np.mean(Intg[-1]),np.std(Intg[-1]),wrelease
        else:
            return np.mean(Intg[-1]),np.std(Intg[-1])


    def run_setup(self):
        """
        This wrapper sets up the simulation and runs through equilibration.
        :return:
        """
        print('Translational and rotational restraints will be imposed.')

        if self.action2 == 'overwrite':
            if os.path.isdir('windows'):
                shutil.rmtree('windows')
                os.makedirs('windows')

        if 'yes' in self.jacks.lower():
            print('Conformational restraints will be applied.')
            self.make_files_and_directories('release')
        else:
            print('No conformational restraints applied.')
        
        if self.hmr == 'yes':
            print('Hydrogen mass repartitioning (HMR) is on.\n')

        if self.perturb == 'yes':
            print 'The Original GAFF parameters have been perturbed. Pleace check the parmed.log file'
            print 'in each umbrella sampling window to make sure the parameters were perturbed as intended.\n' 

        print ('You can use Ctrl+C to quit the program.\n')

        self.make_files_and_directories('attachment')
        self.make_files_and_directories('translation')
  
        # Check parameter files listed in tleap.in
        tleap_in =  open('windows/a00/tleap.in','r') 
        for line in tleap_in:
            splitline = line.split()
            for str in splitline:
                if '.frcmod' in str or '.mol2' in str and splitline[0][0]!='#':
                    if not os.path.exists('setup/param_files/%s'%(str)):                                        
                        print ('Aborted. %s cannot be found in the setup/param_files folder.'%(str))
                        print 'Please make sure all the mol2 and frcmod files listed in tleap.in are available in the setup/param_files directory.\n' 
                        sys.exit(1)
        tleap_in.close() 
        # find out the serial numbers of ligand atoms based on the user input 
        self.lig_resid, self.lig_resname = select_ligand_atoms(self.lig_name) 

        self.prepare_and_simulate('attachment')
        self.prepare_and_simulate('translation')
        if 'yes' in self.jacks.lower():
            self.prepare_and_simulate('release')

    def run_production(self):
        """
        This wrapper runs through the production phase.
        :return:
        """
        if self.hmr == 'yes':
            print('Hydrogen mass repartitioning (HMR) is on.\n')
        print ('You can use Ctrl+C to quit the program.\n')
        self.simulate('attachment')
        self.simulate('translation')
        if 'yes' in self.jacks.lower():
            self.simulate('release')

    def run_analysis(self):
        """
        This wrapper runs thermodynamic integration calculations
        :return:
        """
        bind_dg = 0
        bind_sem = 0
        attachment_dg,attachment_sem = self.ti_analysis('attachment')
        print 'all in kcal/mol                   Mean                SEM'
        print ('%-25s%15.6f%20.6f'%('Attachment:', attachment_dg,attachment_sem))
        bind_dg += attachment_dg
        bind_sem += attachment_sem**2
        translation_dg,translation_sem,wrelease = self.ti_analysis('translation')
        print ('%-25s%15.6f%20.6f'%('Translation (pulling):',translation_dg,translation_sem))
        bind_dg += translation_dg
        bind_sem += translation_sem**2
        if 'yes' in self.jacks.lower():
            release_dg,release_sem = self.ti_analysis('release')
            print ('%-25s%15.6f%20.6f'%('Release:',-1*release_dg,release_sem))
            bind_dg += -1*release_dg
            bind_sem += release_sem**2
        print ('%-25s%15.6f%18s'%('Release-standard:',wrelease,'N/A'))
        print ('%-25s%15.6f%20.6f'%('Binding free energy:',-1*(bind_dg + wrelease),np.sqrt(bind_sem)))

def signal_handler(*args):
    """
    This function is called when the user presses Ctrl+C and will kill the threads launched with subprocess.call().
    """
    p.kill()
    print('Killed process %s' % p)
    sys.exit(0)

def equilibrate(pmemd_option, prmtop):
    """
    Thermalize and equilibrate each window.
    :return:
    """
    print('Equilibrating ...')
    sp.call('%s -O -i min.in -p %s -c solvated.inpcrd -o min.out -r solvated.rst7 '
            '-inf mdinfo -ref solvated.inpcrd'%(pmemd_option, prmtop), shell=True)
    p = check_eq_outputs('min.out')
    if (p == 1):
        return p

    sp.call('%s -O -i therm1.in -p %s -c solvated.rst7 -o therm1.out -r therm1.rst7 '
            '-x therm1.nc -inf mdinfo -ref solvated.inpcrd'%(pmemd_option, prmtop), shell=True)

    p = check_eq_outputs('therm1.out')
    if (p == 1):
        return p

    sp.call('%s -O -i therm2.in -p %s -c therm1.rst7 -o therm2.out -r therm2.rst7 '
            '-x therm2.nc -inf mdinfo -ref solvated.inpcrd'%(pmemd_option, prmtop), shell=True)

    p = check_eq_outputs('therm2.out')
    if (p == 1):
        return p

    sp.call('%s -O -i eqnpt.in -p %s -c therm2.rst7 -o eqnpt.out -r eqnpt01.rst7 -x eqnpt.nc '
            '-inf mdinfo -ref solvated.inpcrd'%(pmemd_option, prmtop), shell=True)

    p = check_eq_outputs('eqnpt.out')
    if (p == 1):
        return p

    # Do 50 runs of short equilibration for each window...
    for run in range(1, 50):
        sp.call('%s -O -i eqnpt.in -p %s -c eqnpt%02d.rst7 -o eqnpt.out -r eqnpt%02d.rst7 '
                '-x eqnpt.nc -inf mdinfo -ref solvated.inpcrd' % (pmemd_option, prmtop, run, run + 1), shell=True)
        p = check_eq_outputs('eqnpt.out')
        if p == 1:
            return p
        os.remove('eqnpt%02d.rst7' % run)
    return 0


def read_bool_attributes(input, attribute):
    """
    Read values of the attributes that require 'yes' or 'no' as input.
    """
    if input.lower() == 'yes' or input.lower() == 'no':
        return input.lower()
    else:
        if attribute == 'amber16':       
            hint = 'whether the MD will be performed using Amber16 or an older version of Amber.\n'
        if attribute == 'perturb':
            hint = 'whether parameters need to be perturbed by ParmEd.'
            hint += ' An input file including the new values of parameters (new_params.dat) is then required in the folder setup/param_files.\n'
        elif attribute == 'hmr':
            hint = 'whether hydrogen mass repartitioning will be used.\n'
        elif attribute == 'strip_water_ions':
            hint = 'whether to strip water and ions in the MD trajectories.\n'
        elif attribute == 'jacks':
            hint = 'whether the conformational restraints are needed.\n'
        elif attribute == 'warning':
            hint = 'whether the warning messages will be printed out.\n'        
        print 'Wrong input! Please use yes or no to indicate', hint
        sys.exit(1)

def check_eq_outputs(output):
    try:
        sp.check_output(['grep -q TIMINGS {}'.format(output)], shell=True)
    except sp.CalledProcessError as e:
        return 1

def select_ligand_atoms(lig_str):

    resid_list = []
    resname_list = []

    lig_str = lig_str.strip(':, ')
    lig_str = filter(None,lig_str.split(','))
 
    for str in lig_str:
        if RepresentsInt(str):
            resid_list.append(int(str))
        elif '-' in str:
            substr = str.split('-')
            if RepresentsInt(substr[0].strip()) and RepresentsInt(substr[0].strip()):
                for i in range(int(substr[0]), int(substr[1])+1):
                    resid_list.append(i)   
        else:
            resname_list.append(str.strip())          

    return resid_list, resname_list

def RepresentsInt(val):
    try: 
        int(val)
        return True
    except ValueError:
        return False

def ismyinstance(variable_type, parameter_value, filename, parameter):
    """
    Check the data type of the user inputs
    """
    if not parameter_value:
    # If the value of entry is left as blank
         if variable_type == 'string':
             return 'None'
         elif variable_type == 'list':
             return []
         elif variable_type == 'float':
             return 0.0
         elif variable_type == 'int':
             return 0
    
    # if input was provided
    if variable_type == 'float':
        try:
            float(parameter_value)
        except ValueError:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter a float value for %s.' % parameter)
            sys.exit()

        if float(parameter_value) < 0:
            print('Wrong input in %s:' % (filename))
            print(parameter_value)
            print('Please enter a non-negative value for %s.' % parameter)
            sys.exit()
        else:
            return float(parameter_value)

    elif variable_type == 'int':
        try:
            int(parameter_value)
        except ValueError:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter an int value for %s.' % parameter)
            sys.exit()

        if int(parameter_value) < 0:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter a non-negative value for %s.' % parameter)
            sys.exit()
        else:
            return int(parameter_value)

    elif (variable_type == 'list') or (variable_type == 'string'):
        return parameter_value

def welcome_message():
    print('**********************************************************************************')
    print(' Welcome to APR: a tool for binding calculations. Version 1.1')
    print(' Written by:')
    print('   Niel M. Henriksen, Jian (Jane) Yin, David R. Slochower')
    print(' and the project leader:')
    print('   Michael K. Gilson ')
    print(' Copyright (c) 2016-2017, University of California, San Diego')
    print('**********************************************************************************')
    print('The current APR scripts may not be directly applied to protein systems, especially')
    print('for those with buried binding sites. Careful adjustments of the protocols and scripts') 
    print('will be needed, based on the requirements of every particular system.\n') 
    
def help_message():
    print('python2 apr.py')
    print('    eq            Set up the APR framework and run the equilibration')
    print('    prod          Run the production phase')
    print('    analysis      Run the analysis and print the final results\n')
    print('-s  Simulation options to be taken:')
    print('    overwrite     Start freshly from the first umbrella sampling window.')
    print('                  Be aware that all exsiting files in windows will be overwritten.\n')
    print('    continue      Pick up from where it was interrupted. You can also use this for the first time run.')
    print('                  The exsiting files (including the input files) in windows will not be overwritten.\n')
    print('-i  APR input file\n')
    print('For example:')
    print('python2 apr.py eq -i apr.in -s continue')
    print('python2 apr.py prod -i apr.in -s overwrite')
    print('python2 apr.py analysis -i apr.in (analysis does not need the -s flag)\n')
    print('For more details, please visit the APR tutorial on the AMBER website:'
          'http://ambermd.org/tutorials/advanced/tutorial29/\n')

def check_versions():
    """
    Check the versions of python and openmm
    """
    print 'Checking the version of Python installed ...'
    if sys.version_info[0]==2 and sys.version_info[1] == 7:
      print ('You are using Python 2.7, which is perfect for running the APR scripts.\n')
    else:
      print ('You are using Python %d.%d. Please install or switch to Python 2.7 to run APR scripts.'%(sys.version_info[0],sys.version_info[1]))
      sys.exit(1)

def main():
    # Register the signal handler to catch Ctrl+C
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    # Every time we call an external program, assign the subprocess p and then when we catch Ctrl+C,
    # we should be able to also run p.kill()
    global p
    welcome_message()
    check_versions()
    this = APR()
    this.check_arguments()
    if this.action1 == 'setup':
        this.process_input_file()
        this.check_executable()
        # Run setup and equilibration
        this.run_setup()
        print('Done!')
        print('Now you can issue the command python2 apr.py prod -i apr.in -s continue to run the production phase.')

    elif this.action1 == 'production':
        # Production runs
        this.process_input_file()
        this.check_executable()
        this.checkfile('production')
        this.run_production()
        print('Done!')
        print('Now you can issue the command python2 apr.py analysis -i apr.in to analyze your data.')

    elif this.action1 == 'analysis':
        this.process_input_file()
        this.checkfile('analysis')
        this.run_analysis()

        
    else:
        print('I could not understand. Please choose something to do among eq (equilibration), prod (production) or analysis.')
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
