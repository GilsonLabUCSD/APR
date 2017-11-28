#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys


def setup_solvate(warning, solvent_model, solvents, ion_list, isAmber16):
    """
    Solvate the system to a certain number of solvent molecules.
    :param target_solvents: the target number of solvent molecules in the window
    :param solvent_box: type of solvent box (must be a string that tleap recognizes)
    :param ions: type of counterions to be added for neutralization
    :return:
    """
    print('Solvating ...')
    # Copy tleap.in to vac_tleap (because this might be needed to analyze the trajectories later)
    shutil.copy('tleap.in', 'vac_tleap.in')

    # Append vac_tleap.in
    vac_in_file = open('vac_tleap.in', 'a')
    vac_in_file.write('check model\n')
    vac_in_file.write('saveamberparm model vac.prmtop vac.inpcrd\n')
    vac_in_file.write('savepdb model dry.pdb\n')
    vac_in_file.write('quit\n')
    vac_in_file.close()

    # Generate vacuum structures
    p = sp.call('tleap -s -f vac_tleap.in > vac_tleap.log', shell=True)

    # Find out the type and number of neutralizing ions needed
    # Search the phrase "WARNING: The unperturbed charge of the unit:" in vac_tleap.log
    f = open('vac_tleap.log', 'r')
    for line in f:
        if "WARNING: The unperturbed charge of the unit:" in line:
            splitline = line.split()
            if float(splitline[7]) < 0:
                ion_list[0][1] = round(float(re.sub('[+-]', '', splitline[7])))
            elif float(splitline[7]) > 0:
                ion_list[1][1] = round(float(re.sub('[+-]', '', splitline[7])))
    f.close()
    
    if solvent_model == 'SPCE':
        solvent_box = 'SPCBOX'
    else:
        solvent_box = solvent_model +'BOX'

    # Update target_solvents according to the number of ions; 
    target_solvents = solvents + int(ion_list[0][1]) + int(ion_list[1][1]) + int(ion_list[2][1]) + int(ion_list[3][1])

    # Initialize a log file to record the thickness of the solvent buffer
    # How much we should adjust the buffer in each iteration, trying to match the correct number of solvent molecules
    buffer_adjustment = 0.1  # Angstroms
    buffer_thickness = 9.0  # Angstroms

    # Call tleap for an initial estimate
    write_leap_input(solvent_model, buffer_thickness, isAmber16)
    current_solvents = run_leap_and_read_output('tmp_tleap.in', 'tmp.log')

    count = 0
    max_count = 100
    refinement_threshold = 10

    # Generate a warning message when the number of solvent molecules requested is not close to our estimation.
    if warning == 'yes' and abs(current_solvents - solvents) > 500:   
        if current_solvents - solvents > 500:
            print ('Warning: The requested number of solvent molecules added to the simulation box (%d) is probably a lot less than what is actually needed.'%(solvents))
        elif current_solvents - solvents < 500: 
            print ('Warning: The requested number of solvent molecules added to the simulation box (%d) is probably a lot more than what is actually needed.'%(solvents))
        print ('We recommend you to try something close to our estimated number of solvent molecules (%d).'%(current_solvents)) 
        print ('If you do not want to see this warning message, switch the warning attribute from ON to OFF in the APR input file.\n')        

    manual_removal = None

    while current_solvents != target_solvents:
        count += 1
        if count > max_count:
        # Try a larger refinement threshold
             refinement_threshold += 1
        # If we are really close, redefine the target_adjustment to make smaller steps.
        if abs(current_solvents - target_solvents) < 20:
            buffer_adjustment = 0.01
        # If we are really really close and we have more solvent molecules than necessary, we can manually remove a couple
        if current_solvents > target_solvents and (current_solvents - target_solvents) < refinement_threshold:
            difference = current_solvents - target_solvents
            manual_removal = [target_solvents + 1 + i for i in range(difference)]
            write_leap_input(solvent_model, buffer_thickness, isAmber16, manual_removal)
            run_leap_and_read_output('tmp_tleap.in', 'tmp.log')
            break
        # If we have more solvents than we want, reduce the buffer thickness
        if current_solvents > target_solvents:
            buffer_thickness -= buffer_adjustment
        # Otherwise, we need more solvents, increase the buffer thickness
        if current_solvents < target_solvents:
            buffer_thickness += buffer_adjustment
        write_leap_input(solvent_model, buffer_thickness,isAmber16)
        current_solvents = run_leap_and_read_output('tmp_tleap.in', 'tmp.log')    

    # Compose the final leap input file from the solvate leap input stub and the final number of solvent molecules
    shutil.copy('tleap.in', 'solvate_tleap.in')
    sol_file = open('solvate_tleap.in', 'a')
    sol_file.write('# load the solvent parameters\n')
    if isAmber16!='yes' and solvent_model!='TIP3P':
        sol_file.write('# load the solvent parameters\n')
        sol_file.write('loadamberparams frcmod.%s\n'%(solvent_model.lower()))
    elif ((isAmber16=='yes') and (solvent_model in ['TIP3P', 'TIP4PEW', 'OPC'])):
        sol_file.write('source leaprc.water.%s\n'%(solvent_model.lower()))
    elif ((isAmber16=='yes') and (solvent_model == 'SPCE')):
        sol_file.write('source leaprc.water.%s\n'%(solvent_model.lower()))
        sol_file.write('loadamberparams frcmod.%s\n'%(solvent_model.lower()))
    elif ((isAmber16=='yes') and (solvent_model in ['CHCL3','MEOH','NMA'])): 
        sol_file.write('loadOff solvents.lib\n')
        sol_file.write('loadamberparams frcmod.%s\n'%(solvent_model.lower()))

    sol_file.write('solvatebox model ' + solvent_box + ' {10.0 10.0 ' + str(buffer_thickness) + '}\n\n')
    if manual_removal is not None:
        for solvent in manual_removal:
            sol_file.write('remove model model.%s\n' % solvent)
        sol_file.write('\n')
    # Make the system neutral by putting a final '0'
    # For the tutorial, this should add 9 ions at the expense of solvent molecules
    if ((ion_list[0][1] != 0) and (solvent_model.lower() not in ['chcl3','meoh','nma'])):
        sol_file.write('addionsrand model %s 0\n' % (ion_list[0][0]))
    if ((ion_list[1][1] != 0) and (solvent_model.lower() not in ['chcl3','meoh','nma'])):
        sol_file.write('addionsrand model %s 0\n' % (ion_list[1][0]))
    if (ion_list[2][1] != 0):
        sol_file.write('addionsrand model %s %d\n' % (ion_list[2][0], ion_list[2][1]))
    if (ion_list[3][1] != 0):
        sol_file.write('addionsrand model %s %d\n' % (ion_list[3][0], ion_list[3][1]))
    sol_file.write('\n')
    #sol_file.write('desc model\n')
    sol_file.write('savepdb model solvated.pdb\n')
    sol_file.write('saveamberparm model solvated.prmtop solvated.inpcrd\n')
    sol_file.write('quit')
    sol_file.close()
    p = sp.call('tleap -s -f solvate_tleap.in > solvate_tleap.log', shell=True)
    
    f = open('solvate_tleap.log', 'r')
    added_solvents = target_solvents
    for line in f:
        if "Could not open file" in line:
           print 'WARNING!!!'
           print line
           sys.exit(1)
        if ("WARNING: The unperturbed charge of the unit:" in line) and warning == 'yes':
           print line
           print ('The system is not neutralized properly. Are you sure you want to continue?')
           print('Maybe you should check the solvated_tleap.log file, and ion types specified in the APR input file.')
        if "addIonsRand: Argument #2 is type String must be of type: [unit]" in line:
           print('Aborted. The ion types specified in the APR input file could be wrong.')
           print('Please check the solvation.log file, and the ion types specified in the APR input file.\n')
           sys.exit(1)
    f.close()

    # Clean up
    sp.call('rm tmp.log', shell=True)
    sp.call('rm tmp_tleap.in', shell=True)
    sp.call('rm leap.log', shell=True)


def write_leap_input(solvent_model, local_buff, isAmber16, manual_removal=None):
    """
    Write a tleap input file that solvates according to the solvent type and buffer size in Angstroms.
    :param solvent_box: type of solvent box (must be a string that tleap recognizes)
    :param local_buff: size of z-axis buffer in Angstroms
    :param manual_removal: can be a list to tleap to move a few, specific solvents.
    :return:
    """
    # Copy tleap.in to tmp_tleap.in
    shutil.copy('tleap.in', 'tmp_tleap.in')
    tmp_file = open('tmp_tleap.in', 'a')

    if solvent_model == 'SPCE':
        solvent_box = 'SPCBOX'
    else:
        solvent_box = solvent_model +'BOX'

    if isAmber16!='yes' and solvent_model!='TIP3P':  
        tmp_file.write('loadamberparams frcmod.%s\n'%(solvent_model.lower()))
    elif ((isAmber16=='yes') and (solvent_model.lower() in ['tip3p', 'tip4pew', 'opc'])):
        tmp_file.write('source leaprc.water.%s\n'%(solvent_model.lower()))
    elif ((isAmber16=='yes') and (solvent_model	== 'SPCE')):
       	tmp_file.write('source leaprc.water.%s\n'%(solvent_model.lower()))
        tmp_file.write('loadamberparams frcmod.%s\n'%(solvent_model.lower()))
    elif ((isAmber16=='yes') and (solvent_model.lower() in ['chcl3','meoh','nma'])):
        tmp_file.write('loadOff solvents.lib\n')
        tmp_file.write('loadamberparams frcmod.%s\n'%(solvent_model.lower()))
    elif isAmber16=='yes':
        print('Warning: your solvent model is currently not supported by APR.\n')
        print('currently supported types include: TIP3P   TIP4P-Ew   OPC   SPC/E   CHCl3   MeOH   NMA\n') 
        sys.exit(1)
    tmp_file.write('solvatebox model ' + solvent_box + ' {10.0 10.0 ' + str(local_buff) + '}\n')
    if manual_removal is not None:
        for solvent in manual_removal:
            tmp_file.write('remove model model.%s' % solvent)
            tmp_file.write('\n')
    tmp_file.write('quit')
    tmp_file.close()


def run_leap_and_read_output(tleap_in, outlog):
    """
    Run tleap and get the number of solvents for a given buffer padding.
    :return:
    """
    p = sp.call('tleap -s -f %s > %s'%(tleap_in, outlog), shell=True)
    # Search the key word "Added" in tleap log file
    current_solvents = None
    f = open(outlog, 'r')
    for line in f:
        if "Added" in line:
            current_solvents = line[line.find('Added') + len('Added'):].split()
            current_solvents = int(current_solvents[0])
    f.close()
    if current_solvents is None:
        print('It appears there was a problem running tleap.')
        print('Maybe check tleap.in and your AMBER version option in your APR input file.') 
        print('You may also need to adjust the number of solvent molecules requested.\n')
        sys.exit(1)
    return current_solvents

