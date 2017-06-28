### APR 1.1 User Manual (last updated 2017/06) ###

# Contributors to the current code #

    Niel M. Henriksen
    Jian (Jane) Yin
    David R. Slochower

and the project leader:

    Michael K. Gilson 


# Acknowledgements #

The development of APR (attach-pull-release), a binding calculation tool was made possible by support from NIH and 
Air Force Office of Scientific Research (AFOSR) Basic Research Initiative (BRI) grant.


# APR Papers #

Please use the following papers when citing APR method or scripts:

Henriksen NM, Fenley AT, Gilson MK. Computational Calorimetry: High-Precision Calculation of Host-Guest Binding Thermodynamics. 
J. Chem. Theory Comput., 2015, 11(9), 4377-4394. http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00405.

Yin, J, Henriksen, NM, Slochower, DR, Gilson, MK. The SAMPL5 Host-Guest Challenge: 
Computing Binding Free Energies and Enthalpies from Explicit Solvent Simulations by the Attach-Pull-Release (APR) Method 
J. Comput. Aided Mol. Des., 2016. http://link.springer.com/article/10.1007/s10822-016-9970-8.

Velez-Vega C, Gilson MK. Overcoming Dissipation in the Calculation of Standard Binding Free Energies by Ligand Extraction. 
J. Comput. Chem., 2013, 34(27), 2360-2371. http://onlinelibrary.wiley.com/doi/10.1002/jcc.23398/full.


# Introduction #

The APR (attach-pull-release) protocols have been shown to generate moderate to strong correlations between experimental and computational binding thermodynamics based on 
a broad testing of host-guest systems including cucurbit[7]uril (CB7), β-cyclodextrin (β-CD), octa acid (OA) and tetra-endo-methyl octa-acid (TEMOA) with guest molecules. 
For the detailed methodology and validation of APR approach, please refer to the publications listed aboved and a tutorial for Amber users:
http://ambermd.org/tutorials/advanced/tutorial29/.

The current version of APR scripts (v1.1) is a demonstration of how to use pulling approach to compute binding thermodynamics. You can apply it on host-guest complexes 
with a minimal effort of setting up your workflow. 
 
If you are planning on computing the binding affinities of protein-ligand systems, please be aware that careful adjustments of the protocols and scripts will be needed, 
based on the requirements of every particular system. In addition, be extremely cautious about using the APR approach to compute binding affinities for proteins with 
buried binding sites, as those may present convergence issues due to the significant conformational change of the protein during the pulling process. A recent advance of 
APR applications on protein systems can be found here: http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00275.


# What can APR scripts do? #

The scripts will automatically set up the APR framework: building each umbrella sampling window, pulling the guest away from the host, and imposing the necessary restraints 
for the system along the way. Equilibration, simulation, and data analysis are also carried out in a highly-automated way. 


# How does APR work? #

## Installation ##
The APR scripts were written in Python language. To execute them, you need Python version 2.7 installed on your machine. If you don't have access to Python 2.7, 
consider using the miniconda distribution which can be installed optionally in Amber16. You can invoke this Python version with $AMBERHOME/miniconda/bin/python. 
Another alternative is to install Python 2.7 through Anaconda (https://www.continuum.io/downloads).

Work stations that enable GPU acceleration of Amber are also highly recommended for running APR scripts. Please remember to set the environment variable 
"CUDA_VISIBLE_DEVICES" to prevent multiple simulations from running on a single GPU. It is possible to run the tutorial with pmemd.MPI, pmemd, and even sander, 
but it will take significantly longer to finish the calculations.

## Topology and coordinate Files ##
The current APR scripts do not include a built-in docking program. Therefore, a PDB file of the bound structure needs to be provided with all the water molecules and counterions removed,
since those will be added within the APR workflow later on. For small molecules,
mol2 files and possibly frcmod files are required. Those files need to be saved in the ./APR/setup/pdb and ./APR/setup/parameter_files directories, accordingly. They should also be indicated
in ./APR/setup/input_files/tleap.in.amber16, if using Amber16, and in ./APR/setup/input_files/tleap.in.amber14, if using Amber14. The Amber topology (.prmtop) and coordinate (.inpcrd)
files will be generated within the workflow.

## Atom and residue selections ##
The Amber-mask style syntax is adopted for atom and residue selections. The characters ":" and "@" are used for atom and residue selections, respectively. For the residue selections,
either the residue name or the sequence number can be recognized in the latest version of APR. A mixture of those two is also recognizable. Examples are:

    :MOL              Only select residues with the name MOL
    :6-10             Select residues 6 to 10
    :MOL,OCT,3,10-12  Select all residues with names MOL and OCT, as well as residues 3, 10, 11 and 12
    :MOL@C1           Select an atom with the name C1 at the residue MOL

Wildcard features are not supported at this point.  

## System alignment ## 
APR is a pulling approach. To pull the ligand out along a straight line, aligning the system with the longest dimension of the simulation box, the Z-axis will make things a lot easier. We
currently provide a stand-alone script called zalign.py to help the users with this task. The usage of zalign is introduced by issuing the command python zalign.py.
Alternatively, a series of molecular visualization programs such as VMD and Chimera can be used for system alignment.

It is also noteworthy that the residue numbers in the PDB file are subject to change after alignment. The reason is that the tleap program in Amber re-numbers residues if they do not start at "1".
Zalign adopts the same numbering convention used in tleap to keep things consistent. The residue numbers in align_z.pdb (output from zalign.py) should be used for selecting atoms to impose restraints, rather than
those in the original PDB file, if they are not the same.       
     

## Dummy atoms ##
Imposing restraints is crucial for performing binding calculations with the APR method. Distance, angle, and dihedral restraints are usually required to fix the position and orientation 
of the system (For more details see Henriksen NM, Fenley AT, Gilson MK. J. Chem. Theory Comput., 2015, 11(9), 4377-4394).
Anchor particles, i.e. dummy atoms come in handy for setting up restraints. In APR version 1.1, three dummy atoms are appended to the end of the output file
generated by zalign.py (align_z.pdb). These dummy atoms will be assigned by the atom name of "Pb" and the residue name of "DUM", but they have no charge or volume, yet an atomic mass of 220 Da. 
The coordinates of the dummy atoms were hand-coded at this moment, but they can be directly modified in align_z.pdb to adapt to the need of each unique system.    

## Command lines ##
APR is a three-step program which carries out the tasks of running equilibration, production and conducting analysis sequentially. To start the program,
you need to first indicate which step you would like to run with the corresponding keyword:

    eq            Set up the APR framework and run the equilibration
    prod          Run the production phase
    analysis      Run the analysis and print the final results

Meanwhile, two restarting modes, "overwrite" and "continue", are available for both "eq" and "prod" steps. The overwrite mode should be used whenever any setting in the APR input file 
(see below) is changed, or when the tleap.in, PDB or mol2 files have been added, replaced or modified in the setup folder. The overwrite mode will discard
previous equilibration (if used with eq) or simulation results (if used with prod) and start freshly from the first umbrella sampling window. In contrast, the continue mode will 
pick up from where the equilibration or production stops. Restarting modes are indicated with the "-s" flag. Note that both "continue" and "overwrite" can be used for the first time run. 

Examples of the command lines are:

    python2 apr.py eq -i apr.in -s continue
    python2 apr.py prod -i apr.in -s overwrite
    python2 apr.py analysis -i apr.in (analysis does not need the -s flag)


## Test case ##
Input and parameters files for a test system -- octa acid and its guest 4-cyanobenzoic acid were saved within the setup folder in the APR package. Those should be either modified or replaced 
when computing the bind free energy for your own system. The template user input file apr.in was also written based on this example.    


## How to write an APR input file ##
The APR input file can be named by the users and should be indicated by the flag "-i" in the command line. A template of the APR input file, apr.in, is available in the APR package.
Comments starting with a semicolon in this file will not be parsed. The values of the options are case sensitive (except for YES, NO, ON and OFF) while the options themselves are not. 
Not providing options or simply leaving the values blank may not cause abortion of the program, but it will likely cause unexpected consequences. Therefore, it is strongly recommended to specify
all the listed options as instructed in the template file. Options and values are separated with an equal sign (" = "). The spaces before and after the equal sign are not mandatory.
On the other hand, some options contain underscores (" _ ") to make sure that they are parsed as continuous strings. Please do not replace them with spaces or hyphens.   
      
Options are explained below in details.

### Amber16 <YES/NO> ###  
The current APR workflow is closely coupled with the Amber suite of programs. This option is to indicate the version of Amber installed in the computing environment. 
The reason to distinguish different versions is that the tleap module in Amber16 uses slightly different syntax from that in older versions of Amber.

If Amber16 was installed, put YES; otherwise put NO.

### HMR <YES/NO> ###
This option is to indicate whether the hydrogen mass repartitioning is intended to be used in accelerating MD simulations for the equilibration and the production phases. 
If HMR is used, remember to use relatively large step sizes (e.g. 4 fs) for the dt and eq_dt options (see below). The HMR feature requires the installation of 
the Amber built-in or the stand-alone version of parmEd.  

For more details for the HMR technique, please read: Hopkins, Chad W., et al. "Long-time-step molecular dynamics through hydrogen mass repartitioning." 
Journal of chemical theory and computation 11.4 (2015): 1864-1874.
    
### exe_path \<'pmemd.cuda','mpirun -np 12 pmemd.MPI', 'pmemd', 'sander'...> ###
Executables to run MD simulations, depending on your computing environment. 

### temperature \<float> ###
The same temperature will be used for equilibration, production and data analysis.

### perturb <YES/NO> ###
GAFF (general Amber force field) will be used to parameterize small molecules. This option allows the users to run MD with perturbed GAFF parameters. Currently only the feature of perturbing
LJ parameters (radius and epsilon) is available. If the value of this option is specified as YES, new parameters need to be listed in a file named new_parameters.dat and saved in the
./APR/setup/param_files directory. A template file, new_parameters.dat is available in the package. This feature also requires the installation of 
the Amber built-in or the stand-alone version of ParmEd (http://parmed.github.io/ParmEd/html/index.html).    

### attach_list \<a python list of float values> ###
The number of windows in the attach phase is equal to the number of elements in the attach_list option. The value of each element indicates the weight (%) of restraints imposed on the ligand atoms
in each umbrella sampling window, which will be printed out in the first line of the disang.rest file, after the keyword "Weight". 

### translate_list \<a python list of float values> ###
The number of windows in the pulling phase is equal to the number of elements in the translate_list option. The value of each element defines the how far the ligand will be pulled away
relative to its original position in each umbrella sampling window, which will printed out in the first line of the disang.rest file, after the keyword "TransDist". 

It is OK to write the list into multiple lines as long as it is wrapped by brackets.

### distance_force \<float> ###
Force constant for the distance restraints imposed on the receptor and ligand atoms, in the unit of kcal/mol/Angstrom**2. This is only for the translational and rotational restraints.
The force constant of conformational restraints (if any) is specified separately using the jacks_force option.

### angle_force \<float> ###
Force constant for the angle and torsion restraints imposed on the receptor and ligand atoms. The unit is kcal/mol/rad**2. This is only for the translational and rotational restraints.
The angle and torsion restraints are currently not supported for the conformational restraints (distance restraints only). 

### water_model \<string> ###
This option supports all the water models available in Amber: TIP3P, TIP4P, TIP4PEW, TIP5P, OPC etc.

### waters \<int> ###
The target value of water molecules to be added to the simulation box in each window. The tleap program in Amber does not provide the feature of solvating a system with a fixed number of 
water molecules. To add the same number of water molecules to the simulation box in each window, the APR solvation module adopts an iterative approach which will gradually reduce
the difference between the amount of added water and the target value. 

### warning <YES/NO> ###
Estimations can be made in each umbrella sampling window based on the size of the system about how many water molecules are needed to solvate the system. A warning message can be printed out
if the difference between the estimated number and target value (via the option waters) is larger than 500. Put yes to turn this warning on, and no to turn it off if 
the users are certain about their choices.

### neutralizing_cation \<string> ###     
The type of cations used for neutralization. Any ion name that can be recognized by tleap should work here. The number of the neutralizing cations will be determined automatically,
therefore there is no need to specify it here. For instance,
if the system has a net negative charge of 10, a totally of 10 neutralizing cations (assuming monovalent) will be added during the solvation process, and neutralizing anions will not be
added (see below). For positively charged systems, neutralizing cations will not be added, so it is okay to leave the default value (Na+) there, or just simply make a blank option.

### neutralizing_anion	\<string> ###
The type of anions used for neutralization. Similarly, the number of the neutralizing cations will be determined automatically. Neutralizing anions will not be added to negatively charged
systems.

### cations \<string> ###
This option is for adding extra cations on top of the neutralizing ones. If not needed, specify the number of cations (via the option number_cations) as 0.

### number_cations \<int> ###
The number of extra cations added to mimic the buffer conditions. Also see cations.

### anions \<string> ###
For adding extra anions on top of the neutralizing ones. If not needed, specify the number of anions (via the option number_anions) as 0.

### number_anions \<int> ###
The number of extra anions added to mimic the buffer conditions. Also see anions.

### LIG \<an Amber-mask style string for residue selection> ###
This option is for the users to indicate all the atoms that belong to the ligand. This information is needed for manipulating the coordinates of the ligand atoms in each 
umbrella sampling window in order to pull the ligand out. Both residue names and residue sequence numbers with the Amber mask style are acceptable (e.g. :MOL :5). 
Directly selecting atoms are currently not supported for this option.   

### R1 \<an Amber-mask style string for atom selection> ### 
One of the three receptor atoms selected for imposing restraints. Examples are: :OCT@C1, :1@N2. 

R1, R2 and R3 (see below) should not lie in the same line, and they should locate at
relatively rigid regions of the receptor in order to maintain robust restraints.

### R2 \<an Amber-mask style string for atom selection> ###
One of the three receptor atoms selected for imposing receptor restraints. Also see R1.

### R3 \<an Amber-mask style string for atom selection> ###
One of the three receptor atoms selected for imposing receptor restraints. Also see R1.

### L1 \<an Amber-mask style string for atom selection> ###
One of the two ligand atoms selected for imposing ligand restraints. L1 is also the origin of the aligned structure in align_z.pdb. Given how the restraints are currently set up,
L1 should locate between L2 (see below) and the dummy atoms. (For more detailed illustrations, see Figure S1 in Henriksen et al. http://pubs.acs.org/doi/suppl/10.1021/acs.jctc.5b00405.)

### L2 \<an Amber-mask style string for atom selection> ###
One of the two ligand atoms selected for imposing ligand restraints. 

### eq_dt \<int> ###
The step size for the NPT runs in the equilibration phase, in the unit of femtosecond (fs). Normally it is 2 fs. If hydrogen mass repartitioning is used (controlled by the option HMR), 
this value can be increased to around 4.
Note that the equilibration phase consists of four steps: a minimization step, a NVT phase of 1 ps at 10 K, a 100 ps heat-up from 10 K 
to a user-defined target temperature (via the option temperature), and lastly, 50 cycles of NPT runs. Note This option only controls the step size for the last step of NPT runs. The step sizes of other steps were 
hand coded as 2 fs.

### eq_nstlim \<int> ###     
The number of steps for the NPT runs in the equilibration phase. For example, if the step size is 2 fs (controlled by the option eq_dt), assigning it a value of 2500 will launch a NPT run of 
5 ps per cycle and 250 ps in total (the number of cycles is fixed as 50).

### eq_barostat \<1/2> ###
Barostat options for the NPT runs in the equilibration phase; 1 is Berendsen and 2 is Monte Carlo.

### eq_cutoff \<int> ###
Cutoff distance of the van der Waals for the 50 cycles of NPT runs in the equilibration phase, in the unit of angstrom. Long-range corrections will be automatically
applied on top of it (see Amber manual).

### dt \<int> ###
The step size for the production runs, in the unit of fs. Normally it is 2 fs. If hydrogen mass repartitioning is used (controlled by the option HMR),
this value can be increased to around 4.

### nstlim \<int> ###
The number of steps for the production runs. For example, if the step size is 2 fs (controlled by the option dt), assigning it a value of 2500000 will launch a 
production run of 5 ns per cycle and 100 ns in maximum (the maximum number of cycles is fixed as 20). Note that the maximum of simulation time may not always be achieved. ...

### barostat \<1/2> ###
Barostat options for the production runs; 1 is Berendsen and 2 is Monte Carlo.

### cutoff \<int> ###
Cutoff distance of the van der Waals for the production runs, in the unit of angstrom. Long-range corrections will be automatically applied on top of it (see Amber manual).

### ntpr \<int> ### 
Frequency of printing energy terms to the Amber output file (.mdout). If the step size is 2 fs (controlled by the option dt), giving it a value of 500 will allow the Amber program
to write the output file every 1 ps.

### ntwx \<int> ###
Frequency of recording frames to the trajectory (in NetCDF format). If the step size is 2 fs (controlled by the option dt), giving it a value of 500 will allow the Amber program
to write the trajectory every 1 ps.

### strip_water_ions <YES/NO> ###
This option determines whether the water molecules and counterions will be stripped in the MD trajectories. If yes, only the solute atoms (including the dummy atoms) will be saved in the
trajectories, otherwise all atoms will be stored. This option is corresponding to the ntwprt variable in Amber input files. Only saving the solute atoms is recommended for saving
the disk space, if the analysis later on does not require any information from water and counterions.



The restart files and the trajectories are in the format of NetCDF but can be converted to other formats such as mdcrd or PDB using Cpptraj. If water and ions are stripped, 
vac.prmtop (topology in the gas phase) should be used for parsing the trajectories instead of solvated.prmtop.

If necessary, more options can be modified manually through the apr_mdin.py file. 

more here...

## How to interpret the APR output files ##
To be continued ...
