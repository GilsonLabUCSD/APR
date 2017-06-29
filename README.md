# APR(attach-pull-release): a binding calculation approach

The current version of APR scripts (v1.1) is a demonstration of how to use the pulling approach to compute binding thermodynamics. 
You can apply it on most host-guest complexes with a minimal effort of setting up your system. 
If you are planning on computing the binding affinities for proteins, please be aware that
careful adjustments of the protocols and scripts will be needed, based on the requirements of every particular system. In addition,
be extremely cautious about using the APR approach to compute binding affinities for proteins with buried binding sites, as those
may present convergence issues due to the significant conformational change of the protein during the pulling process.     

The APR protocols have been used to generate moderate to strong correlations between experimental and computational binding thermodynamics based on a broad testing of host-guest systems including cucurbit[7]uril (CB7), β-cyclodextrin (β-CD), octa acid (OA) and tetra-endo-methyl octa-acid (TEMOA) with guest molecules. For the detailed theoretical framework, methodology, and validation of APR, please refer to the following publications:

Velez-Vega C, Gilson MK. Overcoming Dissipation in the Calculation of Standard Binding Free Energies by Ligand Extraction. J. Comput. Chem., 2013, 34(27), 2360-2371. http://onlinelibrary.wiley.com/doi/10.1002/jcc.23398/full

Henriksen NM, Fenley AT, Gilson MK. Computational Calorimetry: High-Precision Calculation of Host-Guest Binding Thermodynamics. J. Chem. Theory Comput., 2015, 11(9), 4377-4394. http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00405
 
Yin, J, Henriksen, NM, Slochower, DR, Gilson, MK. The SAMPL5 Host-Guest Challenge: Computing Binding Free Energies and Enthalpies from Explicit Solvent Simulations by the Attach-Pull-Release (APR) Method J. Comput. Aided Mol. Des., 2016. http://link.springer.com/article/10.1007/s10822-016-9970-8


A recent advance of APR applications on protein systems can be found here: http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00275
