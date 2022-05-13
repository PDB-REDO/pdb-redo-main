#!/bin/tcsh -f

# pdb_redo.csh: a general method for optimising macromolecular crystallographic structures.
#
# PDB_REDO installation directory. EDIT THIS!
setenv BASE /zata

# This script was created by Robbie P. Joosten, Bart van Beusekom, and Wouter Touw
# Correspondence to r.joosten@nki.nl or robbie_joosten@hotmail.com
#
# Reference: If you publish results (directly or indirectly) obtained by using this protocol, please refer to (any of)
# these publications:
# 1) Robbie P. Joosten, Gert Vriend: "PDB improvement starts with data deposition" Science, 317, p. 195-196 (2007)
# 2) Robbie P. Joosten, Thomas Womack, Gert Vriend and Gerard Bricogne: "Re-refinement fromdeposited X-ray data can
#    deliver improved models for most PDB entries"  Acta Cryst. D65, p. 176-185 (2009)
# 3) Robbie P. Joosten, Jean Salzemann, Vincent Bloch, Heinz Stockinger, Ann-Charlott Berglund, Christophe Blanchet, Erik
#    Bongcam-Rudloff, Christophe Combet, Ana L. Da Costa, Gilbert Deleage, Matteo Diarena, Roberto Fabbretti, Geraldine
#    Fettahi, Volker Flegel, Andreas Gisel, Vinod Kasam, Timo Kervinen, Eija Korpelainen, Kimmo Mattila, Marco Pagni,
#    Matthieu Reichstadt, Vincent Breton, Ian J. Tickle, Gert Vriend: "PDB_REDO: automated re-refinement of X-ray
#    structure models in the PDB" J. Appl. Cryst., 42, p. 376-384 (2009)
# 4) Robbie P. Joosten, Tim A.H. te Beek, Elmar Krieger, Maarten Hekkelman, Rob W.W. Hooft, Reinhard Schneider, Chris
#    Sander, Gert Vriend: "A series of PDB related databases for everyday needs" Nucl. Acids Res., 39, p. D411-D419 (2011)
# 5) Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend, Anastassis Perrakis: "Automatic rebuilding and
#    optimisation of crystallographic structures in the Protein Data Bank" Bioinformatics, 27, p. 3392-3398 (2011)
# 6) Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastass Perrakis: "PDB_REDO: constructive validation, more
#    than just looking for errors" Acta Cryst. D68, p. 484-496 (2012)
#

# Usage notice:

if ($#argv == 0 ) then
  echo "How to use PDB-REDO:"
  echo " "
  echo "Optimising your structure model (locally):"
  echo "'$0 --local --xyzin=a_coordinate_file --hklin=an_mmCIF_file (or --mtzin=an_mtz_file) (--tlsin=tls_file) (--homin=pdb_file) (--restin=cif_restraint_file) (--extin=external_restraints) (--seqin=fasta_sequence_file) --dirout=place_for_output (--flags) (--params=parameter_JSON_file)'"
  echo " "
  echo "Optimising your structure model (as a PDB-REDO webserver):"
  echo "'$0 --server --xyzin=a_coordinate_file --hklin=an_mmCIF_file (or --mtzin=an_mtz_file) (--tlsin=tls_file) (--homin=pdb_file) (--restin=cif_restraint_file) (--extin=external_restraints) (--seqin=fasta_sequence_file) --dirout=place_for_output (--flags) (--params=parameter_JSON_file)'"
  echo " "
  echo "Optimising an entry from the Protein Data Bank (PDB):"
  echo "'$0 [PDBid] (--download) (--tlsin=tls_file) (--homin=pdb_file) (--restin=cif_restraint_file) (--extin=external_restraints) (--dirout=place_for_output) (--flags) (--params=parameter_JSON_file)'"
  echo " "
  echo "Please note:"
  echo "-The --xyzin option takes input in PDB or mmCIF format. Validation of mmCIF format is strict." 
  echo "-The --restin option takes a single restraint file for Refmac in mmCIF format."
  echo " It overrides all the automatic restraint generation."
  echo "-The --extin option takes a single external restraint file in Refmac format."
  echo " The external restraints are added to those generated within PDB-REDO."
  echo "-The --seqin option takes a sequence file in fasta format."
  echo " The sequence is cross checked with the (ATOM and SEQRES) sequence in the pdb file."
  echo "-The --tlsin option takes an extra TLS group definition in Refmac format."
  echo " You can give the --tlsin option multiple times to test more TLS group definitions."  
  echo "-You can add a comment to the output: --comment='This is a comment'"
  echo "-The --download flag lets the required data from be downloaded from PDBe. Without"
  echo " this flag the required data is taken from a local copy of the PDB."
  echo "-The --params option takes a JSON file with the basic debug flags below. The key is" 
  echo " the flag without the '--' and the value is 1 to activate it. The tighter and looser"
  echo " keys take positive integers. Higher numbers increase the tightness (or looseness)"
  echo " of the restraints. The fitlig key takes an array of residue names and the maxres "
  echo " key takes a float. Example:"
  echo '  {"legacy": 1, "fitlig": ["SO4", "GOL"], "crossval":0, "tighter":2}'    
  echo " "
  echo "These basic debug flags are available:"
  echo "--nohyd       : do not add hydrogens (in riding postions) during refinement"
  echo "--legacy      : for legacy PDB entries. R-factor is not checked and the number of "
  echo "                refinement cycles is increased (a lot)"
  echo "--notls       : no TLS refinement is performed"
  echo "--notlsupdate : use TLS, but do not update the tensors in the final refinement"
  echo "--noncs       : no NCS restraints are applied"
  echo "--nojelly     : switch off jelly body refinement"
  echo "--notwin      : no detwinning is performed"
  echo "--newmodel    : always take an updated model from the re-refinement for the rebuilding"
  echo "                steps. Only use this option when all else fails"
  echo "--tighter     : try tighter restraints than usual (use '--tighter --tighter' for even"
  echo "                tighter restraints)"
  echo "--looser      : try looser restraints than usual (use '--looser --looser' for even"
  echo "                looser restraints)"
  echo "--fitlig=RES  : fit ligand with residue name RES (use --fitlig= multiple times to"
  echo "                search for multiple compounds)"
  echo "--noloops     : do not try to complete loops"  
  echo "--nofixdmc    : do not add missing backbone atoms"
  echo "--nopepflip   : no peptide flips are performed"
  echo "--noscbuild   : side chains will not be rebuilt"
  echo "--nocentrifuge: waters with poor density will not be deleted"
  echo "--norebuild   : all rebuilding steps are skipped"
  echo "--nosugarbuild: no (re)building of carbohydrates"
  echo "--noanomalous : ignore all anomalous data if Fmean or Imean are available"
  echo "--maxres=VAL  : cut the resolution to VAL Angstrom"
  echo "--fewrefs     : deals with very small data sets by switching off R-free set sanity checks"
  echo "--crossval    : performs (very lengthy) k-fold cross validation on the final results"
  echo "--intens      : force the use of the intensities from the reflection file"
  echo "--noocc       : do not refine occupancies"
  echo "--notruncate  : do not use truncate to convert intensities to amplitudes"
  echo "--nosigma     : do not use sigF or sigI for scaling etc."
  echo "--nometalrest : do not generate special metal restraints"
  echo "--nohomology  : do not use homology-based restraints"
  echo "--homology    : force homology-based restraints"
  echo "--hbondrest   : use hydrogen bond restraints"
  echo "--nonucrest   : do not use nucleic acid restraints"
  echo "--paired      : force paired refinement"
  echo "--isotropic   : force isotropic B-factors if there are fewer than 30 reflections/atom"
# Hidden echo "--nqa_isotropic   : force isotropic B-factors no matter what"  
  echo "--anisotropic : force anisotropic B-factors if there are more than 13 reflections/atom"
  echo "                and the resolution is better than 1.95A"
  echo " "
  echo "Debug flags only for working on local PDB entries:"
  echo "--nopdb       : any existing PDB file in the working directory is not replaced"
  echo "--nosf        : any existing reflection file in the working directory is not replaced"
  echo "--mtzdb       : the structure factor files are stored in MTZ format"
  echo " "
  echo "Debug flags for performance:"
  echo "--nproc=value : the maximum number of CPU cores you want to use (the default is 1);"
  echo "                values greater than 7 or the actual number of CPU cores on your"
  echo "                system are not recommended"
  echo "--lowmem      : reduce memory consumption for massive models that otherwise stop REFMAC"
  echo " "  
  echo "Citing PDB-REDO:"
  echo "-Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend, Anastassis Perrakis: "
  echo " Automatic rebuilding and optimization of crystallographic structures in the Protein Data Bank"
  echo " Bioinformatics, 27, p. 3392-3398 (2011)"
  echo "-Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastassis Perrakis:"
  echo " PDB_REDO: constructive validation, more than just looking for errors"
  echo " Acta Cryst. D68, p. 484-496 (2012)"
  exit(0)

endif
#
# Required software:
# - PDB-REDO       Download it from http://pdb-redo.eu
# - CCP4 package   Download it from http://www.ccp4.ac.uk
#
# Optional but highly recommended software
# - BlastP         apt install ncbi-blast+
#                  Needed for homology-based functionality
# - YASARA         Download it from http://www.yasara.org
#                  You need at least YASARA model for pictures describing atomic shifts and TLS grouping and
#                  YASARA dynamics for ligand validation
# - x3dna-dssr     Needed for nucleic acid restraints and validation.
#
####################################################### Change log #######################################################
set VERSION = '7.37' #PDB-REDO version

# Version 7.37:
# - Bugfix for twinning detection when running in legacy mode.
# - Bugfix for using PDBPeep.
# - Added hidden option to always use isotropic B-facters, no questions asked: --nqa_isotropic.
# - Bugfix for the --notwin option.
#
# Version 7.36:
# - Ligand validation output is now in JSON format.
# - Ligand fitting is slightly more agressive in terms of map fit in Coot.
# - Replaced carbonanza with a new mmCIF enabled version.
# - WHAT_CHECK checkdb files are no longer stored
# - Anomalous maps are now always made regardless of whether Refmac decides that there are anomalous scatterers.
# - Added contingency for loopwhole crashes.
# - SIGATM and SIGUIJ records are now removed.
#
# Version 7.35:
# - Added contingency for problems with modelcompare.
# - Added contingency for compounds deleted by loopwhole that are on the list for ligand validation.
#
# Version 7.34:
# - Added new nucleic acid restraints and updated base pair validation.
# - Removed torsion restraints for nucleic acids.
# - Removed KRAB.
# - Versions of PDB coordinate file and reflections data are now also stored.
# - Added support for a JSON formatted parameter file.
# - Robustness fix for dealing with poor LINK records.
# - Fixed reporting error for flipper.
# - Changed the option --isotropic to also force isotropic B-facturs at high(-ish) resolution.
# - Added the --anisotropic option to force anisotropic B-factorsdown to 13 reflections/atom.
# - Added a workaround for incomplete TER records.
#
# Version 7.33:
# - In databank mode, if a PDB file does not exist we try to convert the mmCIF file after renaming the chains.
# - Added contingency in case fixDMC fails.
# - Moved to a new version of platonyzer with mmCIF support.
#
# Version 7.32:
# - Dropped support for OSX.
# - Software versions are now written in json.
# - Stopped making the YASARA scenes.
# - Nucleic acid restraints are now used in all but the highest resolution categories.
# - Added command line switch --nonucrest.
# - Added nucleic acid specific validation. 
# - Bugfix in distel.py to return NA if there are no matched restraints.
# - Bugfix in the way seqrescopier is called.
#
# Version 7.31:
# - Updated the model for the R-factor ratio. It is now empirical based on the PDB-REDO databank.
# - Added a model for sigma value of the R-factor ratio.
# - Decision-making based on the R-free Z-score is now based on the R-factor ratio Z-score.
# - In cases of suspiciously low completeness PDBpeep is used to check whether the reported completeness is ellipsoidal 
#   rather than spherical.
#
# Version 7.30:
# - Updated picker to perform better when a choice is forced. 
# - Added additional step to check quality of TLS refinement.
#
# Version 7.29:
# - Changed the search space for geometric restraint weights.
# - Added contingency for corrupted rwcontents output.
# - Changed the 'make connectivity' setting of Refmac to Refmac's default. 
# - Added a molrep run if the R-factor is over 0.50.
# - Now using tortoize to calculate Ramachandran and torsion Z-scores.
# - Using seqrescopier to add back the sequence in all stages.
#
# Version 7.28:
# - The option --nofixdmc switches of fixDMC.
# - Stopped using PDB-care.
# - Isotropic B-factors can be forced with the '--isotropic' flag.
# - Bugfix for when the --download option is used.
#
# Version 7.27:
# - Added contingency for flipper failure.
# - Bugfix for cases where stats fails after the re-refinement. 
#
# Version 7.26:
# - The analysis of density fit changes is now done by the new tool dRSCC, removing the dependency on R.
# - Added contingency in case PDB-care does not produce output.
# - The output for stats is now given in JSON format.
# - PDB-care can now be switched of from the command line (--careless).
# - Ligand fitting is now only done when R-free < 0.3500 to avoid false positives.
# - If stats fails after the re-refinement, it is run without the restraint file.
# - If it still fails PDB-REDO stops.
#
# Version 7.25:
# - Cell dimensions in the reflection data in mmCIF format are now retained.
# - If there is a conflict in cell dimensions, the reflection data are now leading.
# - If the coordinate file does not have cell dimensions they are added from the reflection data. 
# - Removed 'CRYST1' based tests for input coordinate file.
# - Binliner now reads cell dimensions from the reflection data.
# - Added workaround for cases where metal sites are no longer detected before the final refinement.
# - Unnumbered REMARK records are now stripped out.
#
# Version 7.24:
# - Moved the WHAT_CHECK validation of the input model to after the re-refinement so that it can run in parallel.
# - Some consistency fixes of which input file is used.
# - If prepper cannot parse the input file more verbose output is given to the user and workarounds may be triggered.
#
# Version 7.23:
# - Moved to a new implementation of stripper (now called prepper) that can work with mmCIF.
# - Cleaned up file naming in the initial model editing steps.
# - Added workaround for non-valid XML output from PDB-care.
# - Removed the second pass of PDB-care.
# - Users get an explicit warning if their anomalous data is too incomplete on one of the sides.
# - Bugfix to copying back the SEQRES records to the output.
# - Made the RMSZ cut-offs for cases with RMSZ > 1 tighter.
#
# Version 7.22:
# - Removed the html output and the file compression for the pdb-redo as this output is now generated by the server.
# - Removed the making of boxplots and dRSCC bar charts as this is now also done by the pdb-redo server.
# - The data completeness check is now skipped if paired refinement is forced. 
# - Paired refinement output is now more verbose.
#
# Version 7.21:
# - Added workaround for cases where SFCHECK fails.
# - Replaced all Perl code and removed most environment variables and a overal reduction of variables.
# - The EXPDTA record is now also written out. 
# - Lowered the initial map threshold for ligand fitting.
# - Fitted ligands are now subjected to occupancy refinement.
# - Paired refinement can now be forced with the '--paired' option. 
# - In server mode the output data is no longer precomressed.
#
# Version 7.20:
# - First implementation of carbivore that adds, extends, and corrects N-glycans. This adds a dependancy on 
#   privateer-validate.  
# - Stats now uses the user-provided or auto-generated restraint files if available.
# - Bugfix in handling of corrupt user-provided restraint files.
# - Bugfix in stats stopping errors in server mode.
# - B-factors are reset to 2.00A^2 in stripper if lower.
# - Added contingency if there is no suitable model after the second refinement.
# - Made the cut-offs for ligand fitting less strict to decrease the false negative rate.
#
# Version 7.19:
# - New tool carbonanza adds missing LINKs between ASN and NAG/NDG.
# - Problems with the density plots now cause PDB-REDO to stop.
# - Replaced the Python program rnbterror.py by an improved C++ implementation.
#
# Version 7.18:
# - Bug fix to fixDMC removes false positives for OXT atoms.
# - Centrifuge, pepflip and SideAide exclude lists are now also populated with platonyzer output.
# - Loop building now has an explicit filter on RSCC, which means that the initial fitting is followed Refmac to generate 
#   refined B-factor values and new map coefficients. 
# - Made water deletion in SideAide less agressive.
# - Changed platonyzer restraints.
# - The number of deleted waters is now more reliable.
#
# Version 7.17:
# - Replaced zen by a new program platonyzer that also makes angle restraints for octahedral Na+ and Mg2+.
# - More elegant handling of LINK records.
# - EDIAm and OPIA scores are now reported for new and existing ligands.
# - Changes of EDIAm > |0.1| and OPIA > |12.5%|are now also reported.
#
# Version 7.16:
# - Ligand fitting now uses a new map validation tool that includes the EDIAm metric and removes the Refmac 0-cycle step.
# - Ligands are sorted by size and fitted in order of decreasing size.
# - Because of an update in loopwhole density statistics are also calculated for the besttls model.
#
# Version 7.15:
# - Loop building is now default behaviour. It can be switched off by giving the --noloops flag.
# - Now working with a new version of fixDMC which removes the MCfix dependency.
# - The new fixDMC can handle more complex forms of incomplete residues.
#
# Version 7.14:
# - Added basic mmCIF input support. Currently the mmCIF file is converted to PDB format.
#
# Version 7.13:
# - In server and local mode the additional restraint files are copied to the output directory.
# - Added a fallback and a warning if mmCIF conversion fails.
# - Stability fix that stops Refmac from stopping at later PDB-REDO stages if there are dodgy LINKs. 
# - In databank mode, B-factor overflow now causes a stopping error.
#
# Version 7.12:
# - Low resolution side-chain completion is enforced for side-chains deleted by loopwhole.
# - Missing backbone atoms are now added back with fixDMC.
# - Explicitly defining the radiation type for EDSTATS.
#
# Version 7.11:
# - Now also writing mmCIF coordinate files as output.
# - Dropped support for FoldX.
# - Stability fix for partially merged data sets.
# - Added a second validation step for newly fitted ligands.
#
# Version 7.10:
# - Added loopwhole to complete missing loops.
# - Made the 'omit' map generation for low completeness datasets more elegant. 
# - Ligands can be fitted with the --fitlig= keyword. 
# - The keyword can be given multiple times and fitting occurs sequentially in order of the keywords given.
#
# Version 7.09:
# - Replaced rotacompare with the program modelcompare which also takes over the fuctionality of coot_tour. 
# - New txt2json is now run as python3.
#
# Version 7.08:
# - Added restraint geration with KRAB. Use the --krab switch to activate it.
# - Changed the setting for restraint generation in Refmac to 'connectivity YES'.
# - Simplified the way the version-tracking file is made.
# - In electron diffraction mode, no anomalous maps are made.
#
# Version 7.07:
# - If the R-factor cannot be reproduced in databank mode, but a previous PDB-REDO entry exists, the 0-cycle coordinates
#   are downloaded and used to try to reproduce the R-factors.
# - The hidden --norb switch switches off rigid-body refinement, even in legacy mode.
#
# Version 7.06:
# - More elegant overwriting of databank entries. This does mean that old data in the output directory is removed.
# - The PDBe.json files are now always created.
# - The index.html files are no longer created.
# - WHYNOT and DEBUG messages for server runs go to a seperate WHYNOT file.  
#
# Version 7.05:
# - Updated data for percentile calculation and boxplots.
# - Now also calculating percentiles for coarse packing.
# - CB atoms of GLY are removed automatically.
#
# Version 7.04:
# - REMARK 350 records (biological assembly) are now retained if available.
# - Workaround for cases with alternate residues that get mangled by the rebuilding tools.
# - Workaround for cases with a rediculous amount of TLS groups.
# - Better handling of carbohydrate LINKs with alternates.
#
# Version 7.03:
# - With '--homin=homologous_pdb_file' users can add additional homologous structures that are not yet in the PDB to 
#   improve the homology restraints. 
# - This option automatically triggers homology restraints.
# - Bugfix in writing out the refmac script for low resolution cases.
# - The DNA/RNA restraints from LibG are no longer a hidden option.
# - Fix importing external restraints.
#
# Version 7.02:
# - Added symmetry support for Zen. This required generation of a temporary mtz file early in PDB-REDO.
# - More elegant handling of TLS groups not working in TLSanl. 
# - Experimental phases without HL coefficients or FOMs are supplemented with default FOMs.
# - If the NCS alignment failes due to terminal insertion codes. The residues are now automatically renumbered.
# - Bugfix for cases that have problems with ctruncate and have relatively small test sets.
#
# Version 7.01:
# - Changed the slider values for PDBe to make sure that the step from neutral to positive happens for values greater than
#   2.58.
# - More elegant handling of unstable TLS refinement.
# - More elegant handling of cases where no TLS file could be made.
# - Better handling of cases where Refmac refuses to make anomalous maps.
# - More elegant handling of cases where LibG does not make all types of restraints. 
# - Bug fix in the way Phaser errors are intercepted.
# - Fix for reading logfiles of backbone-only models. 
#
# Version 7.00:
# - rmsZ values for (homology-based) H-bond restraints are now reported. 
# - Better handling for missing high resolution data.
# - Changed the way data.txt is written to avoid race conditions.
#
# Version 6.29:
# - The output PDB files now also have SEQRES records.
# - Added flipper to do DEFY standardisation flips. This replaces the WHAT_CHECK run on the not-quite final pdb file.
# - Flipper is also run before extractor to avoid potential harmful DEFY flips by other programs.
# - Workaround for ctruncate problems if there are anomalous intensities.
#
# Version 6.28:
# - Input MTZ files with many systematic absent reflection now no longer cause infinite loops.
# - Fix for twinning detection by PHASER.
# - Specified NOHARVEST for the REFMAC jobs.
# - Fix for new version of WHAT_CHECK.
#
# Version 6.27:
# - Homology restraints are now automatically used in the third lowest resolution category.
# - Added the --nohomology flag to not use homology restraints.
# - Added H-bond satisfaction percentile for the server.
# - Writing out percentiles to data.txt.
# - Removed the link to EDS.
# - More robust treatment of alternative directory structures and compression of local PDB/reflection files.
# - The besttls files are now gzipped in the databank.
# - Missing percentiles are no longer in <em> tags on the server.
#
# Version 6.26:
# - Homology restraints are now also used in the second lowest resolution category.
# - Fix to deal with = characters in the file names.
# - Carbohydrate residues are now only renamed if they are part of N-glycans.
# - PDB-care is run differently based on the presence or absence of CONECT records.
#
# Version 6.25:
# - Added a twin test from PHASER to reduce twinning false positives.
# - Fixed bug in handeling corrupt restraint files.
# - SEGID records (not part of the official PDB format) are now deleted.
# - Added the --mtzdb flag.
#
# Version 6.24:
# - The configuration of PDB-REDO is moved to a separate file which makes updating easier.
# - Changed the geometric weight for mini-rsr to 10.00 based on feedback from Paul Emsley.
# - Users can give a sequence file.
# - Now using Refmac's default treatment of sugars.
# - Stability fix for parsing WHAT_CHECK output.
# - Better handling of incorrect external restraint files.
#
# Version 6.23:
# - Fix in the map handeling for pepflip.
# - Improved error message for atom naming problems.
# - Bug fix for PHASER output detection.
# - Kollumer now ignores anomalous data with very low completeness.
# - Homology restraints are now automatically used in the lowest resolution category.
# - If the resolution cut-off is changed the solvent mask parameters are re-optimised.
# - Bugfix in YASARA atom shift scene generation.
#
# Version 6.22:
# - A JSON version of data.txt is now also written.
# - Started using WHAT_CHECK 14.
# - Stability fix for homology restraints.
# - Better handling of rediculously high (total) B-factors after TLS refinement.
#
# Version 6.21:
# - Updated the information in the data.txt that allows regeneration of all relevant warnings in the index.html files.
# - The entry creation date is now also stored in data.txt.
# - Stabitility fix in the real-space map validation.
# - Fixed support for users without FoldX.
# - Substantial reduction in messages to debug.txt.
#
# Version 6.20:
# - Added workaround for Refmac deleting LINK records after rigid-body refinement.
# - If ctruncate fails or writes out an unusable reflection file the I to F conversion is performed by cif2cif.
# - Models with strict NCS are now also automatically run in databank mode.
#
# Version 6.19:
# - Ligand validation in YASARA is now also parallelised.
# - Fixed bug in the residue name extraction for validated ligands.
# - The wavelength is now read form the PDB header if no value is given in the reflection file.
# - Stability fix for k-fold cross validation.
# - Switched to the CCP4 version of syminfo.lib for the rebuilding tools.
#
# Version 6.18:
# - Switched back to the mini-rsr (now named coot-mini-rsr) distributed with CCP4.
# - At low resolution the B-factor model selection is slightly less conservative.
# - A debug statement is written for legacy entries with calculated R-free > 0.50.
# - PDB-REDO is stopped if the resolution of the reflection data is higher than 0.30A, which is an indication reflection
#   data problems.
# - Stability fix in the real-space validation.
#
# Version 6.17:
# - Anomalous data is now used to make anomalous maps, but not for refinement.
# - Generating anomalous maps is not compatible with twinning or phased refinement.
# - Anomalous data can be ignored with the '--noanomalous' flag.
# - Small improvements to how the MTZ file is made.
#
# Version 6.16:
# - The completeness of the data is now reported in data.txt.
# - At very low completeness (less than 50%) and omit map is created with COMIT. This map is used for water deletion,
#   side-chain rebuilding and selection of peptide flipping candidates. The map is of too low quality for mini-rsr.
# - Made centrifuge a bit less aggressive in water deletion.
# - The new program rotacompare replaces YASARA for rotamer and peptide analysis.
# - External restraint files are now checked to see whether they are not regular mmCIF restraint files.
#
# Version 6.15:
# - Change in pepflip that should improve the initial filtering of flipping candidates.
# - Changed resolution cut-off for peptide flipping to 3.3A (from 3.5A).
# - Switched to a new rotamer dictionary for SideAide, based on the Top8000/Top7200.
# - Reporting of missing wavelengths is now limited to entries from 2007 onwards.
#
# Version 6.14:
# - A JSON-formatted datafile for PDBe is now written in databank mode.
# - A debug message is written when the wavelength is not reported in the reflection file.
# - Residues named WAT are now excluded from ligand validation, UNL is now included.
# - The output from the rebuilding tools is now slightly less verbose.
# - Stability fix in making the dRSCC plots.
#
# Version 6.13:
# - The solvent percentage is now also given for structures with strict NCS.
# - Hydrogen bond restraints now also apply to side chains.
# - Added optional homology based restraints, this requires a local installation of BlastP.
# - Homology-based restraints currently only work in databank mode.
# - Activate homology-based restraints with --homology.
#
# Version 6.12:
# - LINKs are now added based on the output of PDB-care.
#
# Version 6.11:
# - CIF2CIF now keeps phase information from the input reflection file.
# - Several related stability fixes to cif2cif.
# - Started using CCP4 7.0.
#
# Version 6.10:
# - First attempt to deal with electron diffraction data. It is not very sophisticated yet.
# - The type of experiment is now written to data.txt.
# - Writing of success.txt is more elegant now.
# - Implemented a new version of detectHbonds.
#
# Version 6.09:
# - Changed to a new version of mini-rsr that should improve pepflip's performance. This temporarily adds a lot of
#   dependencies.
# - The dependency problem will be solved once CCP4 starts distributing COOT 0.8.3.
# - Added a specific file for debug messages from the rebuilding tools in PDB-REDO.
#
# Version 6.08:
# - Failure of the second TLSANL run now triggers a debug message.
#
# Version 6.07:
# - Using a new version of pepflip that uses maps in which ligands an metals are masked out. This should improve the real-
#   space refinement in mini-rsr and reduce severe cases of severe backbone distortion.
# - Pepflip now uses the mini-rsr distributed with CCP4.
# - Fix to the second run of zen.
# - If different metal sites are discoverd in the second run of zen, the number of refinement cycles is increased.
# - Directories are now deleted without prompting.
#
# Version 6.06:
# - Now in k-fold cross validation, the coordinate perturbation is now very small and the number of refinement cycles is
#   increased.
# - Having more than one CRYST1 card now triggers a fatal error.
#
# Version 6.05:
# - SideAide can now use secondary structure specific rotamers. This is currently only done for Ile.
# - Stability fix for ion restraints.
# - Changed the tolerance for the density fit score in SideAide. This reduces the bias towards new conformations a bit.
# - Side chain flipping is now skipped for structures without protein. This gives a slight speedup for those structures.
# - Now properly deals with different tiers of YASARA.
# - UNL residues are now no longer deleted in local and in server mode.
#
# Version 6.04:
# - Extra metal, hydrogen bond, and nucleic acid restraints are re-evaluated after rebuilding.
#
# Version 6.03:
# - Jelly-body refinement can be switched off with the keyword '--nojelly'.
# - Removed the '-s' flag for picker in TLS group selection for better comparison of TLS group sections consisting of the
#   same number of groups.
# - Fixed the restraint pruning for stacking restraints.
#
# Version 6.02:
# - Zn-Cys4 clusters are now treated automatically. LINKs and angle restraints are added, disfulfide bridge detection in
#   REFMAC is switched off.
# - This option can be switched off with '--nometalrest'.
# - More verbose reporting on additional external restraints.
# - Hydrogen bond restraints can be used with '--hbondrest'.
#
# Version 6.01:
# - Unmerged symmetry related reflections are averaged using SFTOOLS rather than randomly selected by CAD.
# - The B-factor resetting now has a minumum of 10A^2.
# - The detection of disulfide bridges can be switched of with the hidden keyword --noss.
# - Moved the cispeptide detection to a separate YASARA run to improve stability.
# - Strict NCS and local NCS are no longer mutually exclusive.
# - Now doing a proper geometric weight optimisation in the lowest resolution category.
#
# Version 6.00:
# - First version with full OS X support.
# - Started reporting the solvent content from the density-based estimate from RWCONTENTS. Note that it is unreliable for
#   models with deuteriums, many alternates or strictncs.
# - Harmonic restraints are used for refinements with very small data sets (< 1000 reflections). They can be switched off
#   with the keyword '--noharmonic'.
# - Fix for generating dRSCC plots for really small structures.
# - Separated the restraints from LibG from user-supplied external restraints. The user-supplied restraints have a default
#   scale of 10, but this can be overruled in the restraint file using 'external weight scale [value]'
# - A warning is given when there are unmerged relections.
#
# A complete changelog is available from the PDB-REDO website
#
echo " "
if ($1 == "--local" || $1 == "--server") then
  setenv PDBID `mktemp -u XXXX`
else
  setenv PDBID $1
  #Check the validity of the PDBid
  if ($PDBID != `echo $PDBID | cut -c 1-4`) then
    echo "This is not a valid PDB identifier. Cannot continue."
    exit(1)
  endif
endif
set D2       = `echo $PDBID | cut -c 2-3`    #Essential for certain environment values


# Check wether PDB-REDO is properly configured.
setenv TOOLS $BASE/tools                    #Directory with PDB-REDO software and data
if (! -e $TOOLS/pdb_redo.setup) then
  #PDB_REDO is not properly installed or configured.
  echo "FATAL ERROR!"
  echo "------------"
  echo "PDB-REDO cannot find its installation directory."
  echo "Please, set the correct value for 'BASE' in $0 on line 6."
  exit(1)
else
  #set additional environment parameters (define directories and files).
  source $TOOLS/pdb_redo.setup
endif

########################## Initialise important variables. Usually no need to edit this ##################################

#Add extra libraries
if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH $TOOLS/lib:$LD_LIBRARY_PATH
else
  setenv LD_LIBRARY_PATH $TOOLS/lib
endif


#Refinement parameters
set NCYCLE        = 20      #Initial number of (restrained) refinement cycles
set RBCYCLE       = 10      #Number of rigid-body refinement cycles
set BTESTCYC      = 10      #Number of restrained refinement cycles in the B-restraint weight optimisation
set TLSCYCLE      = 15      #Number of TLS-refinement cycles
set TLSCMD        =         #Empty TLS command
set TLSFILS       =         #No TLS files specified
set DOTLS         = 1       #Use TLS
set TLSUPDATE     = 1       #Update TLS tensors in final refinement
set OPTTLSG       = 'none'  #Fallback value if no TLS is used in refinement
set TIME          = `date +'%F'`
set BCMD          =         #Don't change the B-factors
set WGTTYPE       = MATRIX  #Use weight matrix by default
set BREFTYPE      = ISOT    #Do isotropic B-factor refinement by default
set NCSTYPE       =         #Empty NCS command
set RBNCS         =         #No NCS during rigid-body refinement
set NCSALIGN      =         #Only usefull when NCS is used (alignment cut-off)
set NCSNEIGH      =         #Only usefull when NCS is used (include neighbouring atoms in restraints)
set SOLVENT       = SIMP    #Use simple solvent model by default
set MASKPAR       =         #Use default masking parameters
set JELLY         =         #Do not use 'jelly body' refinement by default
set TORSION       =         #Do not use torsion restraints by default
set LOGSTEP       = 17      #Line number of Refmac log with required output (counted backwards)
set ESTRICT       =         #Picker is not extra strict
set SHIFTCO       = 0.50    #Cut-off for atomic shift in ligand validation
set EXTSCAL       = 10      #Scale for external restraints
set REFIRES       =         #Resolution cut-offs for refinement; use all data by default
set OCCCMD        =         #Occupancy refinement commands
set METALCMD      =         #Command for metal site restraints
set RESTCMD       =         #Command for external restraints
set NUCRCMD       =         #DNA/RNA restraints commands
set HBONDWGTCMD   =         #Command for weighting external hydrogen bond restraints
set HBONDCMD      =         #Command for using hydrogen bond restraints
set HBRESTRWGT    = 2       #Hydrogen bond restraint weight
set NHBONDREST    = 0       #Number of H-bond restraints
set HOMOLWGTCMD   =         #Command for weighting external homology restraints
set HOMOLCMD      =         #Command for using homology restraints
set HOMOLRESTRWGT = 2       #Homology restraint weight
set NHOMOLREST    = 0       #Number of homology restraints
set HOMINCMD      =         #Command for extra homologous structures
set HARMCMD       =         #Harmonic restraints command
set WAVELCACHE    =         #Wavelength cashed from input MTZ file
set MAXPRUNE      = 20      #Give up after removing MAXPRUNE restraints
set SCATLIN       =         #Line that overruls the default scattering factor file
set SCATTERCMD    =         #Use the default X-ray scattering factors
set PHASES        =         #No phase columns as input by default
set ANOMCOEF      =         #Input anomalous data for Refmac
set ANOMCMD       =         #Command for anomalous refinement
set FASTAIN       =         #No input fasta by default
set DICTCMD       =         #No dictionary to append by default  

#Refmac settings
set CONNECTIVITY  =         #Empty value uses Refmac's default other options are 'connectivity [YES|NO|DEFINE]'
set SUGAR         =         #Empty value uses Refmac's default other options are 'sugar [YES|NO|DEFINE]'
set RSYMM         =         #Empty value uses Refmac's default other options are 'symmetry [YES|NO]'

#Default values
set NTLS         = 0            #Zero TLS groups by default. The number of groups is reset when TLS is used
set NPRUNEM      = 0            #No metal restraints pruned yet
set NPRUNEN      = 0            #No nucleic acid restraints pruned yet
set NPRUNE       = 0            #Total number of pruned restraints
set NMETALREST   = 0            #No metal restraints to start with
set NMETALREST2  = 0            #No metal restraints to start with
set NNUCLEICREST = 0            #No nucleic acid restraints to start with
set EXPTYP       = 'X-ray'      #The type of experiment
set GOTOLD       = 0            #By default not starting with previous PDB-REDO entry 
set NUMBERING    =              #The file with the residue numbering mapping from rnbterror
set NLIG_LIST    =              #No newly added ligands
set OXTADD       =              #Add OXT atoms by default
set TRUSTSEQ     = 0            #Trust the sequence by default
set SFFILE       = 'atomsf'     #Default scattering factor file for stats
set SFTYPE       =              #Set scattering type for stats (Default: X-ray) 
set DIDMR        = 0            #No molecular replacement was done 

#Error flags
set TTEST      = 0 #No errors in the TLS group optimisation
set UTRUNCATE  = 0 #Truncate was not used

#Deal with potentially missing sigF/sigI values
set EXPSIG  = 'n'       #Do not use experimental sigmas for scaling
set WGTSIG  =           #Use experimental sigmas for scaling

#Debug flags and derivatives
set NOPDB        = 0
set CEDIT        = false  #The coordinate file is edited from the original PDB entry
set NOSF         = 0
set REDIT        = false  #The relection file is edited from the original PDB entry
set USEMTZ       = 0
set COORDCONV    = 0      #The coordinates were converted from mmCIF version of PDB entry
set DOWNLOAD     = 0
set NOHYD        = 0
set HYDROGEN     = ALL
set LEGACY       = 0
set DOTWIN       = 1
set DOHARMONIC   = 1      #Use harmonic restraints if needed
set TWIN         = 'test' #By default test for twinning
set ISTWIN       = 0      #Data not treated as twinned
set FALSETWIN    = 0      #Data not previously incrorrectly treated as twinned
set ISED         = 0      #this is an electron diffraction set
set DOTASER      = 0
set DOPEPFLIP    = 1
set DOSCBUILD    = 1
set DOCENTRIFUGE = 1
set DOREBUILD    = 1
set DOSUGARBUILD = 1 
set DOFIXDMC     = 1     #Use fixDMC by default
set DONCS        = 1
set DORB         = 1 
set RESOCHECK    = 1
set DOJELLY      = 1     #Allow jelly-body refinement
set DOMETALREST  = 1     #Use special metal restraints by default
set STRICTNCS    = 0
set NCSSTRICT    =       #Empty strict NCS keyword for Refmac
set LOCAL        = 0
set SERVER       = 0
set NPROC        = 1     #Default number of processors to use
set MAXPROC      = 20    #The maximum number of processors to use
set PDBIN        =
set XYZIN        =
set HKLIN        =
set MTZIN        =
set PARAMS       =       #No parameter file specified
set XTLS         =       #User submitted TLS files
set XHOM         =       #User submitted homologous structure models
set INREST       =       #The location of the geometric restraint file
set INEXT        =       #The location of the external restraint file
set INSEQ        =       #The location of the input sequence
set RELAX        =
set NEWMODEL     =
set LOWMEM       =
set FEWREFS      =
set CROSS        = 0
set C2CERR       = 0
set INTENS       =       #Use intensities if set to '-i'
set USTATUS      =       #Use the reflection status column. Set to '-s' to create a new R-free set
set SIGMA        =
set USIGMA       = 1     #Set to 0 to ignore sigma values
set C2CCONV      =       #Set to '-c' to convert intensities to amplitudes with cif2cif instead of ctruncate
set DOOCC        = 1     #Refine occupancies for certain residues
set DONUCR       = 1     #Use nucleic acid restraints by default
set CLEAN        = 0
set RESCATSTEP   = 0     #Do not change the resolution category
set SSBOND       = 'YES' #Autodetect disulfide bridges in Refmac
set HBONDREST    = 0
set DOHOMOLOGY   = 0     #Do not use homology restraints by default
set NOHOMOLOGY   = 0     #Modifier to switch off homology restraints
set FORCEPAIRED  = 0     #Force paired refinenement if set to 1
set SMODE        =       #Running mode for prepper
set C2CANO       = '-a'  #Run cif2cif in anomalous mode
set DOLOOPS      = 1     #Complete loops with loopwhole
set FITLIGANDS   = 0     #Do not fit ligands by default 
set FITLIGS      =       #Ligands to fit
set BISOT        = 0     #Do the Hamilton test to select the B-factor model (ISOT vs OVER and ISOT vs ANISOT)
set FISOT        = 0     #Force isotropic B-factors, no questions asked
set BANISOT      = 0     #Do the Hamilton test to select the B-factor model (ISOT vs ANISOT)
set COMPERROR    = 0     #No error in the completeness

#Other settings
set WEBGET = 'wget -q'
set ALFNUM = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

#Rebuilding settings
set NO_REBUILD =      #List of residues that must not be rebuilt
set H2O_KEEP   =      #List of waters to keep without question

#Legacy values
set OGFOLD = 'NA'     #Output scores from FoldX
set NGFOLD = 'NA'     #Output scores from FoldX
set FGFOLD = 'NA'     #Output scores from FoldX

############################################### You need not edit beyond this line #######################################

#Make sure the LOG file can be created (the :h modifier removes the filename from $LOG leaving only the directory)
mkdir -p $LOG:h

#make a temporary log file
if ($1 == "--local") then
  echo "Results from PDB-REDO $*" > $LOG
endif

#Write out header
echo " __   __   __      __   ___  __   __    __   __  __ " | tee -a $LOG
echo "|__) |  \ |__) __ |__) |__  |  \ /  \    /   __)  / " | tee -a $LOG
echo "|    |__/ |__)    |  \ |___ |__/ \__/   /  o __) /  " | tee -a $LOG
echo " "

#Font for the header
# _     __  __       __  _  __  _   _
#/ \ /|  _) __) |_| |_  |_   / (_) (_|
#\_/  | /__ __)   | __) |_) /  (_)  _|


############################################### Check for debug flags ####################################################
#Get debug flags
echo " " | tee -a $LOG
echo "****** Optimisation setup ******" | tee -a $LOG
foreach ARG ($*)
  if ($ARG == $PDBID) then
    echo "-Evaluating PDB entry $PDBID" | tee -a $LOG
  else if ($ARG == "--local") then
    set LOCAL = 1
    #Set additional flags for programs
    set RELAX = "-r"  #Relaxed mode for extractor
    set SMODE = '--server'  #Server mode for prepper
  else if ($ARG == "--server") then
    set LOCAL  = 1
    set SERVER = 1
    set STDIR  = $cwd
    #Write different debug files
    set WHYNOT = $WHYNOT.server
    set DEBUG  = $DEBUG.server
    #Make the running status file
    touch $STDIR/processRunning.txt
    #Reset the logfile
    mkdir -p $STDIR/output
    cat $LOG    >> $STDIR/output/process.log
    setenv LOG     $STDIR/output/process.log
    #Set additional flags for programs
    set RELAX = "-r"  #Relaxed mode for extractor
    set SMODE = '--server'  #Server mode for prepper
  else if (`echo $ARG | cut -c 1-8` == "--pdbin=") then
    #Get the file in two steps to expand relative paths, but also paths using '~'
    set PDBIN = `echo $ARG | cut -d '=' -f 2-`
    set PDBIN = `readlink -m $PDBIN`
    echo "-PDB-REDO will optimise the structure model in $PDBIN" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--xyzin=") then
    #Get the file in two steps to expand relative paths, but also paths using '~'
    set XYZIN = `echo $ARG | cut -d '=' -f 2-`
    set XYZIN = `readlink -m $XYZIN`
    #Choose the PDB file if two coordinate files are specified
    if ($PDBIN != "") then
      echo "-You specified two coordinate files. Using $PDBIN" | tee -a $LOG
      set XYZIN =
    else
       echo "-PDB-REDO will optimise the structure model in $XYZIN" | tee -a $LOG  
    endif    
  else if (`echo $ARG | cut -c 1-8` == "--hklin=") then
    set HKLIN = `echo $ARG | cut -d '=' -f 2-`
    set HKLIN = `readlink -m $HKLIN`
    #Choose the mmCIF file if two reflection files are specified
    if ($MTZIN != "") then
      echo "-You specified two reflection files. Using $HKLIN" | tee -a $LOG
      set MTZIN =
    else
      echo "-Using the experimental data in $HKLIN" | tee -a $LOG
    endif
  else if (`echo $ARG | cut -c 1-8` == "--mtzin=") then
    #Choose the mmCIF file if two reflection files are specified
    if ($HKLIN != "") then
      echo "-You specified two reflection files. Using $HKLIN" | tee -a $LOG
    else
      set MTZIN = `echo $ARG | cut -d '=' -f 2-`
      set MTZIN = `readlink -m $MTZIN`
      echo "-Using the experimental data in $MTZIN" | tee -a $LOG
    endif
  else if (`echo $ARG | cut -c 1-9` == "--dirout=") then
    set OUTPUT = `echo $ARG | cut -d '=' -f 2-`
    set OUTPUT = `readlink -m $OUTPUT`
  else if (`echo $ARG | cut -c 1-9` == "--params=") then
    set PARAMS = `echo $ARG | cut -d '=' -f 2-`
    set PARAMS = `readlink -m $PARAMS`
  else if (`echo $ARG | cut -c 1-8` == "--tlsin=") then
    set INTLS = `echo $ARG | cut -d '=' -f 2-`
    set INTLS = `readlink -m $INTLS`
    echo "-Testing TLS group definitions from $INTLS" | tee -a $LOG
    #Push the input TLS file on a stack
    set XTLS  = `echo $XTLS $INTLS`
  else if (`echo $ARG | cut -c 1-8` == "--homin=") then
    set INHOM = `echo $ARG | cut -d '=' -f 2-`
    set INHOM = `readlink -m $INHOM`
    echo "-Using homologous structure model $INHOM" | tee -a $LOG
    #Push the input homologous file on a stack
    set XHOM  = `echo $XHOM $INHOM`   
  else if (`echo $ARG | cut -c 1-9` == "--fitlig=") then
    set FITLIG  = `echo $ARG | cut -d '=' -f 2-`
    echo "-Trying to fit ligand $FITLIG" | tee -a $LOG
    #Add the ligand to the list
    set FITLIGS    = `echo $FITLIGS $FITLIG`
    set FITLIGANDS = 1
  else if (`echo $ARG | cut -c 1-9` == "--restin=") then
    set INREST = `echo $ARG | cut -d '=' -f 2-`
    set INREST = `readlink -m $INREST`
    echo "-Using geometric restraints from $INREST" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--extin=") then
    set INEXT = `echo $ARG | cut -d '=' -f 2-`
    set INEXT = `readlink -m $INEXT`
    echo "-Using additional restraints from $INEXT" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--seqin=") then
    set INSEQ = `echo $ARG | cut -d '=' -f 2-`
    set INSEQ = `readlink -m $INSEQ`
    echo "-Using sequence from $INSEQ" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--nproc=") then
    #Using multiple CPUs
    set NPROC = `echo $ARG | cut -d '=' -f 2-`
    #Is the value an integer?
    test $NPROC -eq $NPROC >& /dev/null
    if ($status) then
      #Not an integer, use a single core
      set NPROC = 1
      echo "-Using $NPROC CPU core" | tee -a $LOG
    else
      if ($NPROC < 1) then
        set NPROC = 1
      else if ($NPROC > $MAXPROC) then
        set NPROC = $MAXPROC
      endif
      echo "-Using a maximum of $NPROC CPU cores where possible" | tee -a $LOG
    endif
  else if ($ARG == --nopdb) then
    set NOPDB    = 1
    echo "-Keeping the PDB file in the working directory" | tee -a $LOG
    #Assume the file was edited
    set CEDIT = true
  else if ($ARG == --nosf) then
    set NOSF     = 1
    echo "-Keeping the reflection file in the working directory" | tee -a $LOG
    #Assume the file was edited
    set REDIT = true
  else if ($ARG == --mtzdb) then
    set USEMTZ   = 1
    echo "-The reflection data in the data bank will be treated as MTZ" | tee -a $LOG
    set MTZIN = $SF/r${PDBID}sf.mtz
  else if ($ARG == --download) then
    set DOWNLOAD = 1
    echo "-The required files will be downloaded" | tee -a $LOG
  else if ($ARG == --nohyd) then
    set NOHYD    = 1
    set HYDROGEN = NO
    echo "-Hydrogens will NOT be added in the riding position during refinement" | tee -a $LOG
  else if ($ARG == --legacy) then
    set LEGACY = 1
    echo "-The PDB file will be treated as a legacy PDB entry" | tee -a $LOG
  else if ($ARG == --tighter) then
    @ RESCATSTEP = ($RESCATSTEP - 1)
    echo "-Trying tighter-than-default restraints" | tee -a $LOG
  else if ($ARG == --looser) then
    @ RESCATSTEP = ($RESCATSTEP + 1)
    echo "-Trying looser-than-default restraints" | tee -a $LOG
  else if ($ARG == --notls) then
    set DOTLS = 0
    echo "-No TLS refinement will be performed" | tee -a $LOG
  else if ($ARG == --notlsupdate) then
    set TLSUPDATE = 0
    echo "-No additional TLS refinement will be performed after rebuilding" | tee -a $LOG
  else if ($ARG == --noncs) then
    set DONCS = 0
    echo "-No NCS restraints will be used" | tee -a $LOG
  else if ($ARG == --nojelly) then
    set DOJELLY = 0
    echo "-Not doing jelly-body refinement" | tee -a $LOG
  else if ($ARG == --norb) then
    set DORB = 0
    echo "-Not doing rigid-body refinement" | tee -a $LOG    
  else if ($ARG == --noharmonic) then
    set DOHARMONIC = 0
    echo "-Not using harmonic restraints for small datasets" | tee -a $LOG
  else if ($ARG == --nometalrest) then
    set DOMETALREST = 0
    echo "-Not using additional metal site restraints" | tee -a $LOG
  else if ($ARG == --notwin) then
    set TWIN   =     #No detwinning
    set DOTWIN = 0
    echo "-No detwinning will be performed" | tee -a $LOG
  else if ($ARG == --noanomalous) then
    set C2CANO = ""
    echo "-Anomalous data will be ignored" | tee -a $LOG
  else if ($ARG == --newmodel) then
    set NEWMODEL = "-f" #Force rerefinement to finish with new model
    echo "-Rerefinement will always return a new model" | tee -a $LOG
  else if ($ARG == --nocentrifuge) then
    set DOCENTRIFUGE = 0    #No water deletion
    echo "-Poor waters will not be deleted" | tee -a $LOG
  else if ($ARG == --nopepflip) then
    set DOPEPFLIP = 0    #No peptide flips
    echo "-No peptide flips will be performed" | tee -a $LOG
  else if ($ARG == --noscbuild) then
    set DOSCBUILD = 0    #No side chain rebuilding
    echo "-Side chains will not be rebuilt" | tee -a $LOG
  else if ($ARG == --norebuild) then
    set DOREBUILD = 0    #No rebuilding at all
    echo "-No rebuilding at all" | tee -a $LOG
  else if ($ARG == --nosugarbuild) then
    set DOSUGARBUILD = 0    #No sugar rebuilding at all
    echo "-No sugar (re)building at all" | tee -a $LOG       
  else if ($ARG == --nonucrest) then
    set DONUCR = 0    #Do not generate basepair and stacking restraints
    echo "-Not using DNA/RNA restraints" | tee -a $LOG
  else if ($ARG == --noocc) then
    set DOOCC = 0    #No occupancy refinement
    echo "-No occupancy refinement" | tee -a $LOG
  else if ($ARG == --notruncate) then
    set C2CCONV = '-c'    #Use cif2cif instead of ctruncate to convert intensities
    echo "-Not using truncate" | tee -a $LOG
  else if ($ARG == --lowmem) then
    #Reduce Refmac's memory consumption by not 'cleaning' the solvent mask
    set LOWMEM = "solvent process islands noremove"
    echo "-Working in limited memory mode" | tee -a $LOG
  else if ($ARG == --crossval) then
    set CROSS = 1
    echo "-Full k-fold cross validation will be performed" | tee -a $LOG
  else if ($ARG == --fewrefs) then
    set FEWREFS = "-n"
    echo "-Working with very small data set" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-9` == "--maxres=") then
    set MRESO = `echo $ARG | cut -d '=' -f 2-`
    echo "-The data will be cut to $MRESO A" | tee -a $LOG
  else if ($ARG == --intens) then
    set INTENS = "-i"
    echo "-Using intensities from the reflection file" | tee -a $LOG
  else if ($ARG == --nosigma) then
    set SIGMA  = "-g"
    set USIGMA = 0
    echo "-Not using sigF or sigI values from reflection file" | tee -a $LOG
  else if ($ARG == --paired) then
    set FORCEPAIRED = 1 
    echo "-Performing paired refinement" | tee -a $LOG
  else if ($ARG == --noss) then
    set SSBOND = 'NO'
    echo "-Switched off detection of disulfide bridges" | tee -a $LOG
  else if ($ARG == --noloops) then
    set DOLOOPS = 0
    echo "-No loops will be built" | tee -a $LOG 
  else if ($ARG == --nofixdmc) then
    set DOFIXDMC = 0
    echo "-No correction of missing backbone atoms" | tee -a $LOG     
  else if ($ARG == --hbondrest) then
    set HBONDREST = 1
    echo "-Using hydrogen bond restraints" | tee -a $LOG
  else if ($ARG == --homology) then
    set DOHOMOLOGY = 1
    echo "-Using homology-based restraints" | tee -a $LOG
  else if ($ARG == --nohomology) then
    set NOHOMOLOGY = 1
    echo "-No homology-based restraints will be used" | tee -a $LOG
  else if ($ARG == --isotropic) then
    set BISOT = 1
    echo "-Forcing isotropic B-factors, with fewer than 30 reflections/atom" | tee -a $LOG   
  else if ($ARG == --nqa_isotropic) then
    set FISOT = 1
    echo "-Forcing isotropic B-factors, no questions asked" | tee -a $LOG       
  else if ($ARG == --anisotropic) then
    set BANISOT = 1
    echo "-Forcing anisotropic B-factors" | tee -a $LOG  
  else if (`echo $ARG | cut -c 1-10` == "--comment=") then
    set COMMENT = `echo $* | grep -o -E -e '--comment='+'.*' | sed 's\--comment=\\' | sed 's\--.*\\'`
  else
    if (`echo $ARG | cut -c 1-2` == '--') then
      #It is an invalid keyword
      echo "Invalid keyword: $ARG" | tee -a $LOG
    else
      #Do nothing; it is part of the comment
    endif
  endif
end

#Parse the parameter file
if ("$PARAMS" != "") then
  echo "-Reading additional parameters from $PARAMS" | tee -a $LOG
  
  #Go through all possible parameters
  if (`jq .fitlig $PARAMS` != null) then
    foreach FITLIG (`jq -r '.fitlig[]' $PARAMS`)
      echo " o Trying to fit ligand $FITLIG" | tee -a $LOG
      #Add the ligand to the list
      set FITLIGS    = `echo $FITLIGS $FITLIG`
      set FITLIGANDS = 1
    end
  endif  
  if (`jq .nohyd $PARAMS` == 1) then
    set NOHYD    = 1
    set HYDROGEN = NO
    echo " o Hydrogens will NOT be added in the riding position during refinement" | tee -a $LOG
  endif  
  if (`jq .legacy $PARAMS` == 1) then
    set LEGACY = 1
    echo " o The PDB file will be treated as a legacy PDB entry" | tee -a $LOG
  endif  
  if (`jq .tighter $PARAMS` != null && `jq .looser $PARAMS` > 0) then
    @ RESCATSTEP = ($RESCATSTEP - `jq .tighter $PARAMS`)
    echo " o Restraint tightness set to $RESCATSTEP (lower is tighter)" | tee -a $LOG
  endif  
  if (`jq .looser $PARAMS` != null && `jq .looser $PARAMS` > 0) then
    @ RESCATSTEP = ($RESCATSTEP + `jq .looser $PARAMS`)
    echo " o Restraint tightness set to $RESCATSTEP (higher is looser)" | tee -a $LOG
  endif  
  if (`jq .notls $PARAMS` == 1) then
    set DOTLS = 0
    echo " o No TLS refinement will be performed" | tee -a $LOG
  endif  
  if (`jq .notlsupdate $PARAMS` == 1) then
    set TLSUPDATE = 0
    echo " o No additional TLS refinement will be performed after rebuilding" | tee -a $LOG
  endif
  if (`jq .noncs $PARAMS` == 1) then
    set DONCS = 0
    echo " o No NCS restraints will be used" | tee -a $LOG
  endif 
  if (`jq .nojelly $PARAMS` == 1) then
    set DOJELLY = 0
    echo " o Not doing jelly-body refinement" | tee -a $LOG
  endif
  if (`jq .norb $PARAMS` == 1) then
    set DORB = 0
    echo " o Not doing rigid-body refinement" | tee -a $LOG    
  endif
  if (`jq .noharmonic $PARAMS` == 1) then
    set DOHARMONIC = 0
    echo " o Not using harmonic restraints for small datasets" | tee -a $LOG
  endif
  if (`jq .nometalrest $PARAMS` == 1) then
    set DOMETALREST = 0
    echo " o Not using additional metal site restraints" | tee -a $LOG
  endif
  if (`jq .notwin $PARAMS` == 1) then
    set TWIN   =     #No detwinning
    set DOTWIN = 0
    echo " o No detwinning will be performed" | tee -a $LOG
  endif
  if (`jq .noanomalous $PARAMS` == 1) then
    set C2CANO = ""
    echo " o Anomalous data will be ignored" | tee -a $LOG
  endif
  if (`jq .newmodel $PARAMS` == 1) then
    set NEWMODEL = "-f" #Force rerefinement to finish with new model
    echo " o Rerefinement will always return a new model" | tee -a $LOG
  endif
  if (`jq .nocentrifuge $PARAMS` == 1) then
    set DOCENTRIFUGE = 0    #No water deletion
    echo " o Poor waters will not be deleted" | tee -a $LOG
  endif
  if (`jq .nopepflip $PARAMS` == 1) then
    set DOPEPFLIP = 0    #No peptide flips
    echo " o No peptide flips will be performed" | tee -a $LOG
  endif
  if (`jq .noscbuild $PARAMS` == 1) then
    set DOSCBUILD = 0    #No side chain rebuilding
    echo " o Side chains will not be rebuilt" | tee -a $LOG
  endif
  if (`jq .norebuild $PARAMS` == 1) then
    set DOREBUILD = 0    #No rebuilding at all
    echo " o No rebuilding at all" | tee -a $LOG
  endif
  if (`jq .nosugarbuild $PARAMS` == 1) then
    set DOSUGARBUILD = 0    #No sugar rebuilding at all
    echo " o No sugar (re)building at all" | tee -a $LOG       
  endif
  if (`jq .nonucrest $PARAMS` == 1) then
    set DONUCR = 0    #Do not generate basepair and stacking restraints
    echo " o Not using DNA/RNA restraints" | tee -a $LOG
  endif
  if (`jq .noocc $PARAMS` == 1) then
    set DOOCC = 0    #No occupancy refinement
    echo " o No occupancy refinement" | tee -a $LOG
  endif
  if (`jq .notruncate $PARAMS` == 1) then
    set C2CCONV = '-c'    #Use cif2cif instead of ctruncate to convert intensities
    echo " o Not using truncate" | tee -a $LOG
  endif
  if (`jq .crossval $PARAMS` == 1) then
    set CROSS = 1
    echo " o Full k-fold cross validation will be performed" | tee -a $LOG
  endif
  if (`jq .fewrefs $PARAMS` == 1) then
    set FEWREFS = "-n"
    echo " o Working with very small data set" | tee -a $LOG
  endif
  if (`jq .maxres $PARAMS` != null) then
    set MRESO = `jq .maxres $PARAMS`
    echo " o The data will be cut to $MRESO A" | tee -a $LOG
  endif  
  if (`jq .intens $PARAMS` == 1) then
    set INTENS = "-i"
    echo " o Using intensities from the reflection file" | tee -a $LOG
  endif  
  if (`jq .nosigma $PARAMS` == 1) then
    set SIGMA  = "-g"
    set USIGMA = 0
    echo " o Not using sigF or sigI values from reflection file" | tee -a $LOG
  endif
  if (`jq .paired $PARAMS` == 1) then
    set FORCEPAIRED = 1 
    echo " o Performing paired refinement" | tee -a $LOG
  endif
  if (`jq .noloops $PARAMS` == 1) then
    set DOLOOPS = 0
    echo " o No loops will be built" | tee -a $LOG 
  endif
  if (`jq .nofixdmc $PARAMS` == 1) then
    set DOFIXDMC = 0
    echo " o No correction of missing backbone atoms" | tee -a $LOG     
  endif
  if (`jq .hbondrest $PARAMS` == 1) then
    set HBONDREST = 1
    echo " o Using hydrogen bond restraints" | tee -a $LOG
  endif
  if (`jq .homology $PARAMS` == 1) then
    set DOHOMOLOGY = 1
    echo " o Using homology-based restraints" | tee -a $LOG
  endif
  if (`jq .nohomology $PARAMS` == 1) then
    set NOHOMOLOGY = 1
    echo " o No homology-based restraints will be used" | tee -a $LOG
  endif
  if (`jq .isotropic $PARAMS` == 1) then
    set BISOT = 1
    echo " o Forcing isotropic B-factors, with fewer than 30 reflections/atom" | tee -a $LOG  
  endif
  if (`jq .nqa_isotropic $PARAMS` == 1) then
    set FISOT = 1
    echo " o Forcing isotropic B-factors, no questions asked" | tee -a $LOG  
  endif 
  if (`jq .anisotropic $PARAMS` == 1) then
    set BANISOT = 1
    echo " o Forcing anisotropic B-factors" | tee -a $LOG  
  endif
endif


if ($DOHOMOLOGY == 0 && "$XHOM" != "") then
  set DOHOMOLOGY = 1
  echo "-Using homology-based restraints" | tee -a $LOG
endif  

if ($LOCAL == 1) then
  echo "-During the process your structure model will be called '$PDBID'" | tee -a $LOG
  echo "-The final result will be written to $OUTPUT" | tee -a $LOG
endif

############################################### Check for input files ####################################################

#For local runs make sure we have what we need
if ($LOCAL == 1) then
  #Is there a coordinate file?
  if ($PDBIN == "" && $XYZIN == "") then
    echo " " | tee -a $LOG
    echo "FATAL ERROR!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    echo "You did not provide an input coordinate file. Cannot continue." | tee -a $LOG
    if ($SERVER == 0) then
      #Give a helpfull message
      echo "Use '--xyzin=some_coordinate_file' and try again." | tee -a $LOG
    else
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif
  #Is there reflection data
  if ($HKLIN == "" && $MTZIN == "") then
    echo " " | tee -a $LOG
    echo "FATAL ERROR!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    echo "You did not provide an input reflection file file. Cannot continue." | tee -a $LOG
    if ($SERVER == 0) then
      #Give a helpfull message
      echo "Use '--hklin=some_mmCIF_file' or '--mtzin=some_MTZ_file' and try again." | tee -a $LOG
    else
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif
endif

#Delete previous results?
if ($NOPDB == 1 || $NOSF == 1) then
  cd $WORKDIR   #Go to workdir
  #Remove some files
  rm $WORKDIR/*.tls >& /dev/null
  rm -rf $WORKDIR/download >& /dev/null
  rm $WORKDIR/*.rest >& /dev/null
  rm $WORKDIR/mapval.log >& /dev/null
  rm $WORKDIR/renumber.json >& /dev/null
else
  #Delete all previous results!
  if (-e $WORKDIR) then
    rm -rf $WORKDIR
  endif
  mkdir -p $WORKDIR
  cd $WORKDIR
endif


########################################## Initialise the provenance record ##############################################

#Spawn reusable versions file or bypass the file creation
jq --arg version $VERSION --arg pdbid $PDBID --argjson cedit $CEDIT --argjson redit $REDIT -n '{"data": {"PDBID":$pdbid, "coordinates_revision_date_pdb": null, "coordinates_revision_major_mmCIF": null, "coordinates_revision_minor_mmCIF": null, "coordinates_edited": $cedit, "reflections_revision": null, "reflections_edited": $redit}, "software":{"pdb-redo": {"version":$version, "used":true}}}' > $WORKDIR/versions.json
  
#Add the other programs  
$TOOLS/versions.csh $TOOLS $WORKDIR/versions.json

############################################### Parse the input files ####################################################

#Report
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Parsing input files ******" | tee -a $LOG

#Download the files if the --download flag is given
if ($DOWNLOAD == 1) then
  echo "-Downloading" | tee -a $LOG
  #Make a download dir
  mkdir -p $WORKDIR/download
  cd $WORKDIR/download

  #Download the stuff (reflection data)
  #$WEBGET ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/structure_factors/r${PDBID}sf.ent.gz
  $WEBGET ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/structure_factors/r${PDBID}sf.ent.gz
  if($status) then
    echo " o Cannot download experimental data file" | tee -a $LOG
    exit(1)
  else
    echo " o Downloaded experimental data file" | tee -a $LOG
  endif
  gzip -df r${PDBID}sf.ent.gz
  setenv SF $WORKDIR/download

  #Download the stuff (model)
  #$WEBGET ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb${PDBID}.ent.gz
  $WEBGET ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/pdb/pdb${PDBID}.ent.gz
  if($status) then
    echo " o Cannot download PDB file" | tee -a $LOG
    echo "COMMENT: No PDB-format coordinate file available" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
    exit(1)
  else
    echo " o Downloaded PDB file" | tee -a $LOG
  endif
  gzip -df pdb${PDBID}.ent.gz
  setenv PDB $WORKDIR/download

  #Go back to the working directory
  cd $WORKDIR
endif

#PDB file
if ($NOPDB == 0) then
  if ($LOCAL == 1) then
    if ($XYZIN != "") then
      #Test the file type
      if (-e $XYZIN && `grep -c ^data_ $XYZIN` > 0) then
        #It's not a PDB file assume it is mmCIF and try to convert it.
        $TOOLS/cif2pdb $XYZIN $WORKDIR/pdb${PDBID}.pdb >& $WORKDIR/cif2pdb.log
        if (! -e $WORKDIR/pdb${PDBID}.pdb) then
          echo " " | tee -a $LOG
          echo "FATAL ERROR!" | tee -a $LOG
          echo "------------" | tee -a $LOG
          echo "Cannot use the input file $XYZIN. Please, upload a file in valid PDB or mmCIF format." | tee -a $LOG
          if ($SERVER == 1) then
            #Write out status files
            touch $STDIR/stoppingProcess.txt
            touch $STDIR/processStopped.txt
          endif
          cd $BASE
          exit(1)
        else
          mv $WORKDIR/pdb${PDBID}.pdb $WORKDIR/pdb${PDBID}.ent
        endif           
      else   
        cp $XYZIN $WORKDIR/pdb${PDBID}.ent
      endif
    else if (-e $PDBIN) then
      #Copy PDB file
      cp $PDBIN $WORKDIR/pdb${PDBID}.ent
    else
      echo " " | tee -a $LOG
      echo "FATAL ERROR!" | tee -a $LOG
      echo "------------" | tee -a $LOG
      echo "Cannot find the input file $PDBIN" | tee -a $LOG
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  else
    #Copy PDB file from different directory structures
    if (-e $PDB/pdb${PDBID}.ent) then
      cp $PDB/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.ent
    else if (-e $PDB/$D2/pdb${PDBID}.ent) then
      cp $PDB/$D2/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.ent
    else if (-e $PDB/pdb${PDBID}.ent.gz) then
      zcat $PDB/pdb${PDBID}.ent.gz > $WORKDIR/pdb${PDBID}.ent      
    else if (-e $PDB/$D2/pdb${PDBID}.ent.gz) then
      zcat $PDB/$D2/pdb${PDBID}.ent.gz > $WORKDIR/pdb${PDBID}.ent
    else
      echo "-No PDB-format coordinate file" | tee -a $LOG
      #Test for the existence of a mmCIF file
      if (-e $COORD/$D2/${PDBID}.cif.gz) then
        echo " o Switching to mmCIF-format coordinate file" | tee -a $LOG
        #Count the number of unique 'chains'
        cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software."cif-grep".used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
        
        set NASYM = `$TOOLS/cif-grep -i '_atom_site.auth_asym_id' '.' $COORD/$D2/${PDBID}.cif.gz | sort -u | wc -l`
        set NATOM = `$TOOLS/cif-grep -c -i '_atom_site.auth_asym_id' '.' $COORD/$D2/${PDBID}.cif.gz`
        if ($NASYM < 63 && $NATOM < 100000) then  #Only 62 possible chainIDs and 99999 atoms in PDB format
          echo " o Renaming chains and converting coordinate file" | tee -a $LOG
          #Construct command for mmCQL
          set MMCQLCMD = 
          set IASYM = 0
          #Loop over all single character chains and remove them as renaming options
          set POSSIBLE = $ALFNUM
          foreach ASYM (`$TOOLS/cif-grep -i '_atom_site.auth_asym_id' '.' $COORD/$D2/${PDBID}.cif.gz | sort -u | grep '^.$'`)
            set POSSIBLE = `echo $POSSIBLE | tr -d "$ASYM"`
#            echo $POSSIBLE
          end

          #Rename all chains except the ones that already have a single character
          foreach ASYM (`$TOOLS/cif-grep -i '_atom_site.auth_asym_id' '.' $COORD/$D2/${PDBID}.cif.gz | sort -u | grep '..'`)
            @ IASYM = ($IASYM + 1)
            set CHAIN = `echo $POSSIBLE | cut -c $IASYM-$IASYM`
            set MMCQLCMD = "$MMCQLCMD UPDATE pdbx_poly_seq_scheme SET pdb_strand_id = '$CHAIN' WHERE pdb_strand_id = '$ASYM'; UPDATE atom_site SET auth_asym_id = '$CHAIN' WHERE auth_asym_id = '$ASYM'; "
          end

#          echo $MMCQLCMD 
          #Run mmCQL
          cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.mmCQL.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
          $TOOLS/mmCQL -v --force \
          $COORD/$D2/${PDBID}.cif.gz \
          $WORKDIR/${PDBID}_converted.pdb \
          <<eof >& $WORKDIR/mmCQL_chains.log
            $MMCQLCMD
eof
          
          if (-e $WORKDIR/${PDBID}_converted.pdb) then
            set COORDCONV = 1
            cp $WORKDIR/${PDBID}_converted.pdb $WORKDIR/pdb${PDBID}.ent
          else
            #Report error and stop
            echo " o mmCIF coordinate file cannot be converted" | tee -a $LOG
            echo "COMMENT: No usable coordinate file available" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                              >> $WHYNOT
            cd $BASE
            exit(1)
          endif
        else  
          #Report error and stop
          echo " o mmCIF coordinate file cannot be used" | tee -a $LOG
          echo "COMMENT: No usable coordinate file available" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                              >> $WHYNOT
          cd $BASE
          exit(1)
        endif
      else
        #Report error and stop
        echo "COMMENT: No coordinate file available" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                       >> $WHYNOT
        cd $BASE
        exit(1)
      endif  
    endif
  endif
endif


#Cache the input file
cp $WORKDIR/pdb${PDBID}.ent $WORKDIR/cache.pdb

#Does the input PDB file have SEQRES records?
if (`grep -c ^SEQRES $WORKDIR/pdb${PDBID}.ent` > 0) then
  set GOTSEQRES = 1
else
  set GOTSEQRES = 0
endif


#Check the experimental method (first through the EXPDTA record then by looking for other clues)
if (`grep -c EXPDTA $WORKDIR/pdb${PDBID}.ent` == 0) then
  set ISXRAY = 1 #Assume X-ray and do not check experimental method further
else
  set ISXRAY = `grep EXPDTA $WORKDIR/pdb${PDBID}.ent | head -n 1 | grep -c X-RAY`
endif
if ($ISXRAY == 0) then #Not an X-ray structure, is it electron diffraction?
  set ISED = `grep EXPDTA $WORKDIR/pdb${PDBID}.ent | head -n 1 | grep -c 'ELECTRON CRYSTALLOGRAPHY'`
  if ($ISED == 1) then
    echo "-Switching to electron diffraction mode" | tee -a $LOG
    #Use converted X-ray to electron form factors
    set SCATTERCMD = 'source electron MB'
    set SFTYPE     = '--electron-scattering'
    set SFFILE     = 'atomsf_electron'
    set EXPTYP     = 'ED'
    #Ignore anomalous data
    set C2CANO = ""
  else
    #Create a WHY_NOT comment
    echo "-This is not an X-ray diffraction structure" | tee -a $LOG
    echo "COMMENT: Not an X-ray structure" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                 >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif


#If this is a multi-entry complex, do not re-refine, create WHY_NOT comment
if (`grep -c -E '^SPLIT' $WORKDIR/pdb${PDBID}.ent` != 0) then
  echo "-This is part of a structure divided over multiple PDB entries" | tee -a $LOG
  echo "COMMENT: Multi-entry complex" >> $WHYNOT
  echo "PDB-REDO,$PDBID"              >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#Check for implicit strict NCS
if (`grep -E '^MTRIX' $WORKDIR/pdb${PDBID}.ent | grep -c -v -E '^MTRIX'+'.{54}1'` != 0) then
  set STRICTNCS = 1 
  set NCSSTRICT = 'ncsconstraints'
endif

#If implicit atoms are used as described in REMARK 285, create WHY_NOT comment
if (`grep -c -E '^REMARK 285 THE ENTRY PRESENTED HERE DOES NOT CONTAIN THE COMPLETE' $WORKDIR/pdb${PDBID}.ent` != 0) then
  echo "-This PDB file has implicit atoms" | tee -a $LOG
  echo "COMMENT: Asymmetric unit not complete" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                       >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#If it is a multi-model refinement structure, create WHY_NOT comment
if (`grep -c -E '^ENDMDL' $WORKDIR/pdb${PDBID}.ent` != 0) then
  echo "-This is a multi-model structure" | tee -a $LOG
  echo "COMMENT: Multi-model refinement structure" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                           >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif


#Check to see if it is a CA-only structure
if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | sort -u | wc -l` == 1) then
  #There is only one type of atom. check to see if it is CA
  if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "CA"` != 0) then
    echo "-This is a C-alpha-only structure" | tee -a $LOG
    echo "COMMENT: C-alpha-only structure" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                 >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif

#Check to see if it is a CA/P-only structure
if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | sort -u | wc -l` == 2) then
  #There is only two types of atom. check to see if they are CA and P
  if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "CA"` != 0 && \
      `grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "P "` != 0) then
    echo "-This is a C-alpha/backbone-phosphorus-only structure" | tee -a $LOG
    echo "COMMENT: C-alpha/backbone-phosphorus-only structure" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                     >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif


#Is there only 1 type of residue and is it UNK?
if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 18-20 | sort -u | wc -l` == 1 && \
    `grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 18-20 | grep -c "UNK"` != 0) then
  #Do not use autoNCS (poly-UNKs cannot be alligned properly)
  echo "-The main chain consists only of UNK residues, NCS will not be used!" | tee -a $LOG
  set DONCS = 0
  #Are there ligands or hetero compounds (count the number of different residues)
  if (`grep -E '^ATOM|^HETATM' $WORKDIR/pdb${PDBID}.ent | cut -c 18-20 | sort -u | wc -l` == 1) then
    echo "-This is an UNK-only structure" | tee -a $LOG
    echo "COMMENT: UNK-only structure" >> $DEBUG
    echo "PDB-REDO,$PDBID"             >> $DEBUG
  endif
endif

#Is it a PDB file that came straight from PHASER?
if (`grep -c 'REMARK Log-Likelihood Gain' $WORKDIR/pdb${PDBID}.ent` != 0 && `grep -c 'BUSTER' $WORKDIR/pdb${PDBID}.ent` == 0) then
  #Switch off occupancy refinement
  set DOOCC = 0

  #Report and write a debug message
  echo "-PHASER output model detected. Switching off occupancy refinement." | tee -a $LOG
  echo "COMMENT: PHASER model; no occupancy refinement" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                >> $DEBUG
endif

#Remove SIGATM and SIGUIJ records
if (`grep -c '^SIGATM' $WORKDIR/pdb${PDBID}.ent` > 0 || `grep -c '^SIGATM' $WORKDIR/pdb${PDBID}.ent` > 0) then
  echo "-Removing SIGATM and SIGUIJ records." | tee -a $LOG
  cp $WORKDIR/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.bak
  cat $WORKDIR/pdb${PDBID}.bak | grep -v '^SIGATM' | grep -v '^SIGUIJ' > $WORKDIR/pdb${PDBID}.ent
endif

#Structure factors
if ($NOSF == 0) then
  #Get the user supplied reflection file...
  if ($LOCAL == 1 || $USEMTZ == 1) then
    #If supplied, analyse the input reflection file and convert it to mmCIF if needed...
    if ($MTZIN != "") then
      if (-e $MTZIN) then
        #Dump the MTZ file. PROGRAM: mtzdump
        mtzdmp $MTZIN > $WORKDIR/cifcreate.log
        if ($status) then
          #Give error message
          if ($USEMTZ == 1) then
            #Local databank-type message
            echo "-Cannot parse the input MTZ file." | tee -a $LOG
            echo "COMMENT: Cannot parse the input MTZ file" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                          >> $WHYNOT
          else
            #--local type error
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "The input reflection file $MTZIN:t is not in the MTZ format." | tee -a $LOG
            echo "Please, rerun PDB-REDO with reflections in the MTZ format." | tee -a $LOG
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
          endif
          cd $BASE
          exit(1)
        endif

        #Stop if we are dealing with an unmerged reflection file
        if (`grep -c 'Y  M/ISYM' $WORKDIR/cifcreate.log` != 0) then
          #Give error message
          if ($USEMTZ == 1) then
            #Local databank-type message
            echo "-The input MTZ file has unmerged reflections. Cannot continue." | tee -a $LOG
            echo "COMMENT: Input MTZ file has unmerged reflections" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
          else
            #--local type error
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "The MTZ file $MTZIN:t seems to contain unmerged reflections." | tee -a $LOG
            echo "Please, rerun PDB-REDO with merged reflections." | tee -a $LOG
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
          endif
          cd $BASE
          exit(1)
        endif

        #Cache the wavelength
        set WAVELCACHE = `grep -A 6 'Dataset ID, project/crystal/dataset names, cell' $WORKDIR/cifcreate.log | tail -n 1`

        #Sanity check for the wavelength cache
        if ($WAVELCACHE == 0.00000) then
          set WAVELCACHE =       #Empty value
        endif

        #Get column labels PROGRAM: kollumer
        cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.kollumer.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
        $TOOLS/kollumer -v $WORKDIR/cifcreate.log > $WORKDIR/kollumer.log
        set LABELS = `tail -n 1 $WORKDIR/kollumer.log`

        #Warn if no proper columns were found
        if (`echo $LABELS | grep -c Problem` != 0) then
          #Give error message
          if ($USEMTZ == 1) then
            #Local databank-type message
            echo "-No useful information extracted from MTZ file. Cannot continue." | tee -a $LOG
            echo "COMMENT: No useful data in input MTZ file" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                           >> $WHYNOT
          else
            #--local type error
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            #Give specific warning for anomalous data asymetry
            if (`grep -c "Anomalous intensity data too incomplete!" $WORKDIR/kollumer.log` > 0) then
              echo "Anomalous data too incomplete to use. Please, merge your Friedel pairs." | tee -a $LOG
            else
              echo "Could not extract useful reflection data from $MTZIN:t" | tee -a $LOG
            endif  
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
          endif
          cd $BASE
          exit(1)
        endif

        #Convert the MTZ file to mmCIF. PROGRAM: mtz2various
        mtz2various \
        HKLIN  $MTZIN \
        HKLOUT $WORKDIR/r${PDBID}sf.ent \
        <<eof >> $WORKDIR/cifcreate.log
          OUTPUT CIF data_$PDBID
          LABIN $LABELS
          END
eof

      else
        #Give error message (--local type)
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Cannot find the reflection file $MTZIN" | tee -a $LOG
        if ($SERVER == 1) then
          #Write out status files
         touch $STDIR/stoppingProcess.txt
         touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      endif
    else
      #...or just copy the mmCIF file
      if (-e $HKLIN) then
        cp $HKLIN $WORKDIR/r${PDBID}sf.ent
      else
        #Give error message (--local type)
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Cannot find the reflection file $HKLIN" | tee -a $LOG
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      endif
    endif
  #...or the reflection file from (a local copy of) the PDB
  else
    if (-e $SF/r${PDBID}sf.ent) then
      cp $SF/r${PDBID}sf.ent $WORKDIR/
    else if (-e $SF/$D2/r${PDBID}sf.ent.gz) then
      cp $SF/$D2/r${PDBID}sf.ent.gz $WORKDIR/
      gzip -df $WORKDIR/r${PDBID}sf.ent.gz
    else if (-e $SF/r${PDBID}sf.ent.gz) then
      cp $SF/r${PDBID}sf.ent.gz $WORKDIR/
      gzip -df $WORKDIR/r${PDBID}sf.ent.gz
    else
      echo "-No structure factors found" | tee -a $LOG
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  endif
else
  #Clean some stuff up to avoid problems later on
  if (-e $WORKDIR/$PDBID.cif) then
    rm $WORKDIR/$PDBID.cif
  endif
endif


#Strip out proprietary REMARKS, USER records, gap LINKs, and poor TER records
cat $WORKDIR/pdb$PDBID.ent | grep -v -e '^REMARK [ ,a-z,A-Z][ ,a-z,A-Z][a-z,A-Z]' | grep -v -E '^TER.$' | grep -v -E '^TER$' | grep -v -E '^USER' | grep -v -E 'LINKR.{67}gap' > $WORKDIR/pdb$PDBID.pdb

#Improve LINK compatibility (convert LINK records without distance to LINKR)
sed -i -e '/^LINK .\{69\}     /s/LINK /LINKR/g' -e '/^LINK .\{53\}$/s/LINK /LINKR/g' -e '/^LINK .\{52\}$/s/LINK /LINKR/g' $WORKDIR/pdb$PDBID.pdb

#De-Buster the file if needed
if (`grep -c 'REMARK --------------------- added by autoBUSTER -' $WORKDIR/pdb${PDBID}.pdb` > 0) then
  #Give detailed warning message
  if ($LOCAL == 1) then
    echo " " | tee -a $LOG
    echo "WARNING!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    echo "Your input coordinate file has Buster-specific format errors. Trying to compensate." | tee -a $LOG 
    echo " " | tee -a $LOG
  else
    echo "-The coordinate file has format errors. Trying to compensate." | tee -a $LOG
  endif
  #De-Buster and restart
  cp $WORKDIR/pdb${PDBID}.pdb $WORKDIR/pdb${PDBID}.bak
  set DEBUST = `grep -n 'REMARK --------------------- added by autoBUSTER -' $WORKDIR/pdb${PDBID}.bak | tail -n 1 | cut -d ':' -f 1 | awk '{print $1 +1}'`
  tail -n +$DEBUST $WORKDIR/pdb${PDBID}.bak > $WORKDIR/pdb${PDBID}.pdb
endif

#De-Phenix the file if needed
if (`grep -c 'REMARK IF THIS FILE IS FOR PDB DEPOSITION: REMOVE ALL FROM THIS LINE UP.' $WORKDIR/pdb${PDBID}.pdb` > 0) then
  #Give detailed warning message
  if ($LOCAL == 1) then
    echo " " | tee -a $LOG
    echo "WARNING!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    echo "Your input coordinate file has phenix.refine-specific format errors. Trying to compensate." | tee -a $LOG
    echo " " | tee -a $LOG
  else
    echo "-The coordinate file has format errors. Trying to compensate." | tee -a $LOG
  endif
  #De-Phenix and restart
  cp $WORKDIR/pdb${PDBID}.pdb $WORKDIR/pdb${PDBID}.bak
  grep -A 250000 'REMARK IF THIS FILE IS FOR PDB DEPOSITION:' $WORKDIR/pdb${PDBID}.bak | grep -v 'REMARK IF THIS FILE IS FOR PDB DEPOSITION:' | grep . >  $WORKDIR/pdb${PDBID}.pdb
endif


#Check the cell dimensions
#Get model cell dimensions and space group from coordinate file
if (`grep -c ^CRYST1 $WORKDIR/pdb${PDBID}.pdb` != 1) then
  set MAAXIS = NA
  set MBAXIS = NA
  set MCAXIS = NA 
  set MALPHA = NA
  set MBETA  = NA
  set MGAMMA = NA
  set MSPACE = 
else 
  set MAAXIS = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.pdb | cut -c 7-15`
  set MBAXIS = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.pdb | cut -c 16-24`
  set MCAXIS = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.pdb | cut -c 25-33`
  set MALPHA = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.pdb | cut -c 34-40`
  set MBETA  = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.pdb | cut -c 41-47`
  set MGAMMA = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.pdb | cut -c 48-54`
  set MSPACE = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.pdb | cut -c 56-66 | sed "s/P 21 21 2 A/P 21 21 2 (a)/g" | sed "s/P 1-       /P -1/g"`
endif

#Read the cell dimensions for validation
set RAAXIS = `grep _cell.length_a $WORKDIR/r${PDBID}sf.ent | head -n 1 | awk '{printf "%.3f\n", int(1000*$2 + 0.5)/1000}'`
set RBAXIS = `grep _cell.length_b $WORKDIR/r${PDBID}sf.ent | head -n 1 | awk '{printf "%.3f\n", int(1000*$2 + 0.5)/1000}'`
set RCAXIS = `grep _cell.length_c $WORKDIR/r${PDBID}sf.ent | head -n 1 | awk '{printf "%.3f\n", int(1000*$2 + 0.5)/1000}'`
set RALPHA = `grep _cell.angle_alpha $WORKDIR/r${PDBID}sf.ent | head -n 1 | awk '{printf "%.2f\n", int(100*$2 + 0.5)/100}'`
set RBETA  = `grep _cell.angle_beta $WORKDIR/r${PDBID}sf.ent  | head -n 1 | awk '{printf "%.2f\n", int(100*$2 + 0.5)/100}'`
set RGAMMA = `grep _cell.angle_gamma $WORKDIR/r${PDBID}sf.ent | head -n 1 | awk '{printf "%.2f\n", int(100*$2 + 0.5)/100}'`
set RSPACE = `grep _symmetry.space_group_name_H-M $WORKDIR/r${PDBID}sf.ent | head -n 1 | cut -d "'" -f 2 | sed "s/P 21 21 2 A/P 21 21 2 (a)/g" | sed "s/P 1-       /P -1/g"`
#Contingency for alternative quote use 
if (`echo $RSPACE | grep -c space_group` != 0) then
  set RSPACE = `grep _symmetry.space_group_name_H-M $WORKDIR/r${PDBID}sf.ent | head -n 1 | cut -d '"' -f 2 | sed "s/P 21 21 2 A/P 21 21 2 (a)/g" | sed "s/P 1-       /P -1/g"`
endif

#Give a warning if the cell dimensions do not match
if ($RAAXIS != $MAAXIS || $RBAXIS != $MBAXIS || $RCAXIS != $MCAXIS || $RALPHA != $MALPHA || $RBETA != $MBETA || $RGAMMA != $MGAMMA) then
  #The warning
  echo " " | tee -a $LOG
  echo "WARNING!" | tee -a $LOG
  echo "--------" | tee -a $LOG
  echo "The cell dimensions from the reflection data and the model do not match!" | tee -a $LOG
  echo "From reflections: $RAAXIS $RBAXIS $RCAXIS $RALPHA $RBETA $RGAMMA" | tee -a $LOG
  echo "From model      : $MAAXIS $MBAXIS $MCAXIS $MALPHA $MBETA $MGAMMA" | tee -a $LOG
  echo "Cell dimensions from the input reflection data will be used." | tee -a $LOG
  echo " " | tee -a $LOG
  
  #Write debug message
  echo "COMMENT: cell dimension mismatch" >> $DEBUG
  echo "PDB-REDO,$PDBID"                  >> $DEBUG
 
mmcqlrun:

  echo "-Running mmCQL" | tee -a $LOG
  #Reset the cell dimensions
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.mmCQL.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  $TOOLS/mmCQL -v --force \
  $WORKDIR/pdb${PDBID}.pdb \
  $WORKDIR/${PDBID}_cell.pdb \
  <<eof >& $WORKDIR/mmCQL.log
    UPDATE cell SET length_a = $RAAXIS, length_b = $RBAXIS, length_c = $RCAXIS, angle_alpha = $RALPHA, angle_beta = $RBETA, angle_gamma = $RGAMMA;
eof

  #Check the file
  if (-z $WORKDIR/${PDBID}_cell.pdb || ! -e $WORKDIR/${PDBID}_cell.pdb) then
    #Something is really wrong. Halt and give WHY_NOT comment
    echo "COMMENT: mmCQL: cannot parse structure model" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                              >> $WHYNOT
    #Give detailed help message
    if ($LOCAL == 1) then
      echo " " | tee -a $LOG
      echo "FATAL ERROR!" | tee -a $LOG
      echo "------------" | tee -a $LOG
      echo "Your input model is not valid (enough) PDB or mmCIF format." | tee -a $LOG
      echo "Diagnostic output from mmCQL:" | tee -a $LOG
      echo " " | tee -a $LOG
      grep -E -v 'Expected|Ignoring' $WORKDIR/mmCQL.log | tee -a $LOG
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
    else
      echo "   * The input model is not valid enough to parse. Exit." | tee -a $LOG	  
    endif   
    exit(1)
  endif    
else 
  cp $WORKDIR/pdb${PDBID}.pdb $WORKDIR/${PDBID}_cell.pdb 
endif


#Get versions of the input data (only when working from the PDB)
if ($LOCAL == 0) then
  #Get the version of the reflection data
  set RREVIS = `$TOOLS/cif-grep -i '_audit.revision_id' '.' $WORKDIR/r${PDBID}sf.ent | sort -n | tail -n 1`

  #Assign version
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq --arg rrevis $RREVIS '.data.reflections_revision |= $rrevis' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  #Set the version of the coordinates
  if ($COORDCONV == 1) then
    #mmCIF model is used 
    set CREVMAJOR = `$TOOLS/cif-grep -i '_pdbx_audit_revision_history.major_revision' '.' $COORD/$D2/${PDBID}.cif.gz | sort -n | tail -n 1`
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq --arg crevmajor $CREVMAJOR '.data.coordinates_revision_major_mmCIF |= $crevmajor' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    set CREVMINOR = `echo "SELECT minor_revision FROM pdbx_audit_revision_history WHERE major_revision = $CREVMAJOR;" | $TOOLS/mmCQL $COORD/$D2/${PDBID}.cif.gz >& $WORKDIR/temp.log &&  sort -n $WORKDIR/temp.log | grep -v 'minor_revision' | grep -v 'valid' | tail -n 1`
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq --arg crevminor $CREVMINOR '.data.coordinates_revision_minor_mmCIF |= $crevminor' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  else
    #PDB model is used
    set CDATE  = `grep -e '^REVDAT....  ' $WORKDIR/pdb$PDBID.ent | cut -c 8- | sort -n | cut -c 7-15 | tail -n 1`
    set CDATEC = `date -d $CDATE "+%Y-%m-%d"`
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq --arg cdate $CDATEC '.data.coordinates_revision_date_pdb |= $cdate' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  endif
endif

#Setup restraint library file for ligands and LINKs
if ("$INREST" != "") then
  echo "-Importing geometric restraints" | tee -a $LOG
  cp $INREST $WORKDIR/${PDBID}_het.cif
  #Test the integrety of the restraint file
  if (`grep -a -c '_chem_comp_bond.value_dist' $WORKDIR/${PDBID}_het.cif` == 0 && `grep -a -c '_chem_link_bond.value_dist' $WORKDIR/${PDBID}_het.cif` == 0) then
    #There are no distance restraints. Report...
    echo " o The restraint file contained no valid restraints"   | tee -a $LOG
    echo " o Using standard and on-the-fly generated restraints" | tee -a $LOG
    echo "COMMENT: Corrupt user-provided restraints" >> $DEBUG
    echo "PDB-REDO,$PDBID"                           >> $DEBUG
    # ...and use the regular restraint generation
    set INREST = ""
    rm $WORKDIR/${PDBID}_het.cif
    set LIBLIN = `echo LIBOUT $WORKDIR/${PDBID}_het.cif` #Output Refmac library file for new compounds or LINKs
  else
    #Force Refmac to use the uploaded restraint file
    set LIBLIN  = `echo LIBIN $WORKDIR/${PDBID}_het.cif LIBOUT $WORKDIR/${PDBID}_het2.cif`
    set DICTCMD = "--dict $WORKDIR/${PDBID}_het.cif"
  endif
else
  set LIBLIN = `echo LIBOUT $WORKDIR/${PDBID}_het.cif` #Output Refmac library file for new compounds or LINKs
endif

################################################# Remove unwanted atoms ##################################################

echo "-Preparing the structure model" | tee -a $LOG


renumbered:

#Perform DEFY flips before reannotating any LINKs. PROGRAM: flipper
echo " o Performing DEFY flips" | tee -a $LOG
$TOOLS/flipper -v \
$WORKDIR/${PDBID}_cell.pdb \
$DICTCMD \
-o $WORKDIR/${PDBID}_flipper.pdb >& $WORKDIR/flipper.log

#Contingency for flipper failure
if (! -e $WORKDIR/${PDBID}_flipper.pdb) then
  cp $WORKDIR/${PDBID}_cell.pdb $WORKDIR/${PDBID}_flipper.pdb
  
  #See if this is an error or note
#  if (`grep -c '
  echo "   * Error running flipper"   | tee -a $LOG
  echo "COMMENT: flipper: general error" >> $DEBUG
  echo "PDB-REDO,$PDBID"                 >> $DEBUG
endif

#Run carbonanza to LINK some disconnected sugars. PROGRAM: carbonanza
echo " o Checking carbohydrates" | tee -a $LOG
echo "   * Running carbonanza"   | tee -a $LOG
$TOOLS/carbonanza -v \
$WORKDIR/${PDBID}_flipper.pdb \
$DICTCMD \
-o $WORKDIR/${PDBID}_carbonanza.pdb  >& $WORKDIR/carbonanza.log

#Do we have a new model; then copy the new model
if (-e $WORKDIR/${PDBID}_carbonanza.pdb) then
  if (`grep -c 'Generating LINK record' $WORKDIR/carbonanza.log` > 0) then
    echo "   * Carbonanza added `grep -c 'Generating LINK record' $WORKDIR/carbonanza.log` new LINKs"   | tee -a $LOG
  endif 
else  
  #Copy the old model
  cp $WORKDIR/${PDBID}_flipper.pdb $WORKDIR/${PDBID}_carbonanza.pdb
endif

renumbered:

#Generate metal restraints
if ($DOMETALREST == 1) then
  
  #Run Platonyzer. PROGRAM: Platonyzer
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.platonyzer.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  echo " o Running Platonyzer" | tee -a $LOG
  $TOOLS/platonyzer -v \
  $WORKDIR/${PDBID}_carbonanza.pdb \
  -o $WORKDIR/${PDBID}_platonyzed.pdb \
  --delete-vdw-rest \
  >& $WORKDIR/platonyzer.log
  
  if (-e $WORKDIR/${PDBID}_platonyzed.restraints) then
    
    #Set up the use of the restraints
    if (`cat $WORKDIR/${PDBID}_platonyzed.restraints | wc -l` > 0) then
      mv ${PDBID}_platonyzed.restraints $WORKDIR/metal.rest
      set METALCMD = "@$WORKDIR/metal.rest"
    endif
  else  
    cp $WORKDIR/${PDBID}_carbonanza.pdb $WORKDIR/${PDBID}_platonyzed.pdb
  endif
else  
  cp $WORKDIR/${PDBID}_carbonanza.pdb $WORKDIR/${PDBID}_platonyzed.pdb
endif

#Delete atoms
prepperrun:
echo " o Running prepper" | tee -a $LOG

#Run prepper to remove all hydrogens, deuteriums, atoms with type X, unknown ligands (UNL), UNK atoms not described in
#the refmac library, and superfluous carbohydrate oxygen atoms. Also remove crazy LINKs. PROGRAM: prepper
$TOOLS/prepper \
$WORKDIR/${PDBID}_platonyzed.pdb \
-o $WORKDIR/${PDBID}_prepped.pdb \
-v $SMODE \
--pdb-redo-data $TOOLS/pdb-redo-data.cif \
--debug 3 \
$DICTCMD \
>& $WORKDIR/prepper.log

#Was prepper succesful?
if (! -e $WORKDIR/${PDBID}_prepped.pdb || -z $WORKDIR/${PDBID}_prepped.pdb) then
  if (`grep Error $WORKDIR/prepper.log | grep -c ${PDBID}_het.cif` != 0) then
    #The restraint file is fishy. Warn, but continue.
    echo "COMMENT: prepper: cannot parse restraint file" >> $DEBUG
    echo "PDB-REDO,$PDBID"                               >> $DEBUG
    #Give detailed warning message
    if ($LOCAL == 1) then
      echo " " | tee -a $LOG
      echo "WARNING!" | tee -a $LOG
      echo "------------" | tee -a $LOG
      echo "Your input restraint file has format errors. This may cause refinement problems." | tee -a $LOG
      echo "Diagnostic output from prepper:" | tee -a $LOG
      echo " " | tee -a $LOG
      grep -E -v 'Expected|Ignoring' $WORKDIR/prepper.log | tee -a $LOG
      echo " " | tee -a $LOG
    else
      echo "   * The restraint file has format errors. Trying to compensate." | tee -a $LOG
    endif
    #Rerun prepper without the restraint file.
    set DICTCMD =
    goto prepperrun 
  endif
  #Check for UNL problems
  if (`grep -c "File contains unspecified UNL residue" $WORKDIR/prepper.log` > 0) then
    #Input file contains an UNL that needs restraints
    echo "COMMENT: prepper: UNL without description" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                           >> $WHYNOT
    #Give detailed help message
    if ($LOCAL == 1) then
      echo " " | tee -a $LOG
      echo "FATAL ERROR!" | tee -a $LOG
      echo "------------" | tee -a $LOG
      echo "Your input model contains an UNL residue without restraints." | tee -a $LOG
      echo "Please provide a valid restraint file for UNL" | tee -a $LOG
      echo " " | tee -a $LOG
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
    endif
    exit(1)
  endif    
  #Check for other problems
  if (`grep -c 'Error trying to load file' $WORKDIR/prepper.log` != 0) then
    #The PDB file is fishy. Can the problem be solved?
    if (`grep -c 'REMARK IF THIS FILE IS FOR PDB DEPOSITION: REMOVE ALL FROM THIS LINE UP.' $WORKDIR/${PDBID}_platonyzed.pdb` > 0) then
      #Give detailed warning message
      if ($LOCAL == 1) then
        echo " " | tee -a $LOG
        echo "WARNING!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Your input coordinate file has phenix.refine-specific format errors. Trying to compensate." | tee -a $LOG
        echo "Diagnostic output from prepper:" | tee -a $LOG
        echo " " | tee -a $LOG
        grep -E -v 'Expected|Ignoring' $WORKDIR/prepper.log | tee -a $LOG
        echo " " | tee -a $LOG
      else
        echo "   * The coordinate file has format errors. Trying to compensate." | tee -a $LOG
      endif
      #De-phenix and restart
      cp $WORKDIR/${PDBID}_platonyzed.pdb $WORKDIR/${PDBID}_platonyzed.bak
      grep -A 250000 'REMARK IF THIS FILE IS FOR PDB DEPOSITION:' $WORKDIR/${PDBID}_platonyzed.bak | grep -v 'REMARK IF THIS FILE IS FOR PDB DEPOSITION:' | grep . > $WORKDIR/${PDBID}_platonyzed.pdb	    
      goto prepperrun
    else if (`grep -c 'REMARK --------------------- added by autoBUSTER -' $WORKDIR/${PDBID}_platonyzed.pdb` > 0) then
      #Give detailed warning message
      if ($LOCAL == 1) then
        echo " " | tee -a $LOG
        echo "WARNING!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Your input coordinate file has Buster-specific format errors. Trying to compensate." | tee -a $LOG
        echo "Diagnostic output from prepper:" | tee -a $LOG
        echo " " | tee -a $LOG
        grep -E -v 'Expected|Ignoring' $WORKDIR/prepper.log | tee -a $LOG
        echo " " | tee -a $LOG
      else
        echo "   * The coordinate file has format errors. Trying to compensate." | tee -a $LOG
      endif
      #De-phenix and restart
      cp $WORKDIR/${PDBID}_platonyzed.pdb $WORKDIR/${PDBID}_platonyzed.bak
      set DEBUST = `grep -n 'REMARK --------------------- added by autoBUSTER -' $WORKDIR/${PDBID}_platonyzed.pdb | tail -n 1 | cut -d ':' -f 1 | awk '{print $1 +1}'`
      tail -n +$DEBUST $WORKDIR/${PDBID}_platonyzed.bak > $WORKDIR/${PDBID}_platonyzed.pdb
      goto prepperrun
    else
      #Something is really wrong. Halt and give WHY_NOT comment
      echo "COMMENT: prepper: cannot parse structure model" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      #Give detailed help message
      if ($LOCAL == 1) then
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Your input model is not valid (enough) PDB or mmCIF format." | tee -a $LOG
        echo "Diagnostic output from prepper:" | tee -a $LOG
        echo " " | tee -a $LOG
        grep -E -v 'Expected|Ignoring' $WORKDIR/prepper.log | tee -a $LOG
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
      else
        echo "   * The input model is not valid enough to parse. Exit." | tee -a $LOG	  
      endif   
      exit(1)
    endif
  endif    
endif
  
#Report changes
grep -A 12 Report $WORKDIR/prepper.log | grep -v -E 'Report|CPU|Invalid' | xargs -d '\n' -n 1 echo "   *" | tee -a $LOG

########################## Extract essential information from structure factors and pdb file #############################

#Check structure factor file and reformat
#Label to come back to for forced use of intensities
beintens:

#Run cif2cif
echo "-Checking reflection data" | tee -a $LOG
c2cagain:

#PROGRAM: cif2cif
$TOOLS/cif2cif $C2CANO $USTATUS $FEWREFS $INTENS $SIGMA $C2CCONV \
$WORKDIR/r${PDBID}sf.ent \
$WORKDIR/$PDBID.cif \
$WAVELCACHE \
> $WORKDIR/${PDBID}c2c.log

#Success or not?
if (-e $WORKDIR/$PDBID.cif) then

  #Check for status flag column (only once)
  if ($?GOTR) then
    #GOTR is already set. Do nothing.
  else
    #GOTR is 1 if the status flag (R/Rfree) is used, 0 if not
    set GOTR = `grep -c _refln.status $WORKDIR/$PDBID.cif`
    if ($GOTR == 0) then
      #Check the reason for the missing set
      if (`grep -c 'Warning: too small R-free set.' $WORKDIR/${PDBID}c2c.log` != 0) then
        echo " o Original R-free set too small, a new one will be created" | tee -a $LOG
      else
        echo " o No R-free set defined, a new one will be created" | tee -a $LOG
      endif
    endif
  endif
else if (! $?GOTR) then
  #...warn about the R-free set and try again
  echo " o Error(s) in structure factors; not using status flag" | tee -a $LOG
  set C2CERR = 1
  set GOTR = 0
  set USTATUS = '-s'
  mv $WORKDIR/${PDBID}c2c.log $WORKDIR/${PDBID}c2cv1.log
  goto c2cagain
else
  #Something is really wrong. Halt and give WHY_NOT comment
  echo " o Cannot parse structure factor file. Exit." | tee -a $LOG
  echo "COMMENT: cif2cif: cannot parse structure factors" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
  rm core.*
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  exit(1)
endif

#Check for fishy sigF or sigI data
if (`grep -c 'using the -g switch!' $WORKDIR/${PDBID}c2c.log` != 0) then
  #Give a warning
  echo " o Suspicious sigma values detected; they will be ignored" | tee -a $LOG
  echo "COMMENT: cif2cif: suspicious sigma values" >> $DEBUG
  echo "PDB-REDO,$PDBID"                           >> $DEBUG

  #Go back to cif2cif, now ignorig the sigma column.
  mv $WORKDIR/${PDBID}c2c.log $WORKDIR/${PDBID}c2cv2.log
  mv $WORKDIR/$PDBID.cif $WORKDIR/${PDBID}_badsig.cif
  set SIGMA = '-g'
  goto c2cagain
endif

#Check for multiple wavelength data
if (`grep -c "Second wavelength dataset found!" $WORKDIR/${PDBID}c2c.log` != 0) then
  echo " o Only using reflection data from the first wavelength" | tee -a $LOG
  echo "COMMENT: cif2cif: only using data from first wavelength" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                         >> $DEBUG
endif

#Check to see whether a suitable set of experimental sigmas exists (checked now because DYEAR is needed)
if (`grep -c "sigma values" $WORKDIR/${PDBID}c2c.log` != 0) then
  set EXPSIG  =  'n' #Do not use experimental sigmas for weighting
  set WGTSIG  = "NOEX"
  set USIGMA  = 0
  if (`grep -c "No experimental sigmas found" $WORKDIR/${PDBID}c2c.log` != 0) then
    echo " o No usable experimental sigmas in reflection data"    | tee -a $LOG
    #Give debug message later
  else
    echo " o Cannot use experimental sigmas for weighting" | tee -a $LOG
    #Give debug message later
  endif
endif

#Check for phase information (only for electron diffraction)
if ($ISED == 1) then
  if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0) then
    echo " o Experimental phases in HL format will be used"     | tee -a $LOG
    set PHASES = "HLA=HLA HLB=HLB HLC=HLC HLD=HLD"
    set TWIN   =     #No detwinning
    set DOTWIN = 0
  else if (`grep -c '_refln.phase_meas' $WORKDIR/$PDBID.cif` != 0) then
    echo " o Experimental phases with figures of merit will be used" | tee -a $LOG
    set PHASES = "PHIB=PHIB FOM=FOM"
    set TWIN   =     #No detwinning
    set DOTWIN = 0
  endif
endif

#Report on anomalous data
if (`grep -c '_refln.pdbx_F_plus' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
  if (`echo $PHASES | cut -c 1` == "") then
    echo " o Anomalous data will be used to calculate anomalous maps" | tee -a $LOG
  endif
endif

#Extract information from a pdbfile en the newly created reflection file
#TLS group information is extracted from the PDB header. If no groups are defined, every chain gets its own TLS group
echo "-Extracting data from the PDB file" | tee -a $LOG

#PROGRAM: extractor
$TOOLS/extractor -f $RELAX \
$WORKDIR/${PDBID}_prepped.pdb \
$WORKDIR/${PDBID}.cif \
$CLIBD_MON/list/mon_lib_list.cif \
$TOOLS/pdb_redo.dat \
$WORKDIR/$PDBID.extracted \
$WORKDIR/$PDBID.tls > $WORKDIR/extractor.log

#Success or not?
if (-e $WORKDIR/$PDBID.extracted) then
#Do nothing
else
  #Give PDB-REDO mode-specific error messages
  if ($LOCAL == 1) then
    #Give the long error message
    echo " " | tee -a $LOG
    echo " " | tee -a $LOG
    echo "FATAL ERROR!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    if (`grep -c 'No space group in the PDB file' $WORKDIR/extractor.log` != 0) then
      echo "Space group missing in the PDB file. Please add it." | tee -a $LOG
    else if (`grep -c 'Cannot interpret the PDB file:' $WORKDIR/extractor.log` != 0) then
      echo "Cannot intrepret the input PDB file at this line:" | tee -a $LOG
      grep -A 1 'Cannot interpret the PDB file:' $WORKDIR/extractor.log | tail -n 1 | tee -a $LOG
      echo "Please, ensure that you provide a valid PDB file."
    else
      echo "Could not parse the input PDB file. Please, ensure it is a valid file." | tee -a $LOG
    endif
  else
    #Give the short error message
    echo " " | tee -a $LOG
    echo " o Cannot parse pdb file. Exit." | tee -a $LOG
  endif
  echo "COMMENT: extractor: error using PDB file" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                          >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#Are there any TLS groups from extractor?
if (-e $WORKDIR/$PDBID.tls || -e $WORKDIR/REDO.tls) then
  if (`grep -c 'TLS groups created' $WORKDIR/extractor.log` != 0) then
    echo " o `grep 'TLS groups created' $WORKDIR/extractor.log`" | tee -a $LOG
  endif
  if (`grep -c 'TLS groups extracted' $WORKDIR/extractor.log` != 0) then
    echo " o `grep 'TLS groups extracted' $WORKDIR/extractor.log`" | tee -a $LOG
  endif
else
  echo " o Extractor could not create TLS group definitions" | tee -a $LOG
  echo "COMMENT: extractor: cannot create TLS groups" >> $DEBUG
  echo "PDB-REDO,$PDBID"                              >> $DEBUG
endif

#Check for TLS parsing errors
if (`grep -c "negative selectors" $WORKDIR/extractor.log` != 0 || `grep -c "Cannot read TLS residue range" $WORKDIR/extractor.log` != 0) then
  echo " o Extractor could not parse original TLS group definitions" | tee -a $LOG
  echo "COMMENT: extractor: cannot parse TLS groups" >> $DEBUG
  echo "PDB-REDO,$PDBID"                             >> $DEBUG
endif

#Copy in extra TLS files
if ("$XTLS" != "") then
  echo "-Importing extra TLS group definitions" | tee -a $LOG
  set CNT = 0
  #Loop over files
  foreach FIL ($XTLS)
    set CNT  = `expr $CNT + 1`
    set RANK = `seq -w $CNT 99 | head -n 1`
    #Only copy the group definitions, not the tensors
    grep -E 'TLS|RANGE|^ ' $FIL > $WORKDIR/in$RANK.tls
    #Reject unusable files
    if (`grep -c 'RANGE' $WORKDIR/in$RANK.tls` == 0) then
      echo " o The file $FIL does not contain valid TLS group definitions" | tee -a $LOG
      rm $WORKDIR/in$RANK.tls
      set CNT  = `expr $CNT - 1`
    endif
  end
  echo " o Definitions imported: $CNT" | tee -a $LOG
endif

#Get important parameters from $PDBID/$PDBID.extracted
set AAXIS      = `head -n 1  $WORKDIR/$PDBID.extracted`
set BAXIS      = `head -n 2  $WORKDIR/$PDBID.extracted | tail -n 1`
set CAXIS      = `head -n 3  $WORKDIR/$PDBID.extracted | tail -n 1`
set ALPHA      = `head -n 4  $WORKDIR/$PDBID.extracted | tail -n 1`
set BETA       = `head -n 5  $WORKDIR/$PDBID.extracted | tail -n 1`
set GAMMA      = `head -n 6  $WORKDIR/$PDBID.extracted | tail -n 1`
set RESOLUTION = `head -n 7  $WORKDIR/$PDBID.extracted | tail -n 1`
set DATARESH   = `head -n 8  $WORKDIR/$PDBID.extracted | tail -n 1`
set DATARESL   = `head -n 9  $WORKDIR/$PDBID.extracted | tail -n 1`
set RFACT      = `head -n 10 $WORKDIR/$PDBID.extracted | tail -n 1`
set RFREE      = `head -n 11 $WORKDIR/$PDBID.extracted | tail -n 1`
set NO_REBUILD = `head -n 12 $WORKDIR/$PDBID.extracted | tail -n 1`
#Extend the list if needed
if (-e $WORKDIR/${PDBID}_platonyzed.skip-sideaid) then
  set NO_REBUILD = "$NO_REBUILD`cat $WORKDIR/${PDBID}_platonyzed.skip-sideaid`"
endif
set BAVER      = `head -n 13 $WORKDIR/$PDBID.extracted | tail -n 1`
#set BREF       = `head -n 14 $WORKDIR/$PDBID.extracted | tail -n 1` #Not used
set REFCNT     = `head -n 15 $WORKDIR/$PDBID.extracted | tail -n 1`
set TSTCNT     = `head -n 16 $WORKDIR/$PDBID.extracted | tail -n 1`
set TSTPRC     = `head -n 17 $WORKDIR/$PDBID.extracted | tail -n 1`
set PROG       = `head -n 18 $WORKDIR/$PDBID.extracted | tail -n 1`
set DYEAR      = `head -n 19 $WORKDIR/$PDBID.extracted | tail -n 1`
set SPACEGROUP = `head -n 20 $WORKDIR/$PDBID.extracted | tail -n 1`
set H2O_KEEP   = `head -n 21 $WORKDIR/$PDBID.extracted | tail -n 1`
#Extend the list if needed
if (-e $WORKDIR/${PDBID}_platonyzed.skip-waters) then
  set H2O_KEEP = "$H2O_KEEP`cat $WORKDIR/${PDBID}_platonyzed.skip-waters`"
endif
set BBN_KEEP   = `head -n 22 $WORKDIR/$PDBID.extracted | tail -n 1`
#Extend the list if needed
if (-e $WORKDIR/${PDBID}_platonyzed.skip-pepflipN) then
  set BBN_KEEP = "$BBN_KEEP`cat $WORKDIR/${PDBID}_platonyzed.skip-pepflipN`"
endif
set BBO_KEEP   = `head -n 23 $WORKDIR/$PDBID.extracted | tail -n 1`
#Extend the list if needed
if (-e $WORKDIR/${PDBID}_platonyzed.skip-pepflipO) then
  set BBO_KEEP = "$BBO_KEEP`cat $WORKDIR/${PDBID}_platonyzed.skip-pepflipO`"
endif
set GOT_PROT   = `head -n 24 $WORKDIR/$PDBID.extracted | tail -n 1`
set VDWPROBE   = `head -n 25 $WORKDIR/$PDBID.extracted | tail -n 1`
set IONPROBE   = `head -n 26 $WORKDIR/$PDBID.extracted | tail -n 1`
set RSHRINK    = `head -n 27 $WORKDIR/$PDBID.extracted | tail -n 1`
set COMPLETEH  = `head -n 28 $WORKDIR/$PDBID.extracted | tail -n 1`
set LIG_LIST   = `head -n 29 $WORKDIR/$PDBID.extracted | tail -n 1`
set GOT_NUC    = `head -n 30 $WORKDIR/$PDBID.extracted | tail -n 1`
set WAVELPDB   = `head -n 31 $WORKDIR/$PDBID.extracted | tail -n 1`
set TITLE      = "`tail -n 1 $WORKDIR/$PDBID.extracted`"

#Check for strict NCS and count the number of atoms
#Calculate the number of atoms in the refinement...
set ATMCNT = `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_prepped.pdb`
set NMTRIX  = "1"
#... if strict NCS is used, multiply the apparent number of atoms with the number of MTRIX records
if ($STRICTNCS == 1) then
  #Count the number of MTRIX records
  set NMTRIX = `grep -E '^MTRIX' $WORKDIR/${PDBID}_prepped.pdb | cut -c 7-10 | sort -u | wc -l`
  #Get the total
  set ATMCNT = `echo $NMTRIX $ATMCNT | awk '{print $1*$2}'`
endif


#Get the solvent content. PROGRAM: rwcontents
rwcontents \
XYZIN $WORKDIR/${PDBID}_prepped.pdb <<eof > $WORKDIR/rwcont.log
eof
if ($status) then
  #Give a debug message
  echo "COMMENT: rwcontents failed"   >> $DEBUG
  echo "PDB-REDO,$PDBID"              >> $DEBUG
  set SOLVD = 'NA'
else
  #Get the numbers
  set SOLVD = `grep 'Assuming protein density is' $WORKDIR/rwcont.log | cut -d ':' -f 2 | sed 's/\*\*\*\*\*\*\*\*\*\*\*\*/NA/g'`
endif

#Get the sequence
if ($INSEQ != "") then
  echo " o Importing amino acid sequence" | tee -a $LOG
  cp $INSEQ $WORKDIR/user.fasta
  set FASTAIN = "-fastain $WORKDIR/user.fasta"
endif

#Extract and or check the sequence. PROGRAM pdb2fasta
if ($GOT_PROT == 'T') then
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.pdb2fasta.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json

  echo " o Extracting the amino acid sequence" | tee -a $LOG
  $TOOLS/pdb2fasta \
  -v \
  $FASTAIN \
  -pdb $WORKDIR/${PDBID}_prepped.pdb \
  -fasta $WORKDIR/$PDBID.fasta \
  -tools $TOOLS \
  >& $WORKDIR/pdb2fasta.log
  #make sure there is a fasta file
  if (! -e $WORKDIR/$PDBID.fasta) then
    if ($DOHOMOLOGY == 1) then
      echo "  * Cannot extract sequence. Not using homology restraints." | tee -a $LOG
      set DOHOMOLOGY = 0
    else
      echo "   * Cannot extract sequence." | tee -a $LOG
    endif
    echo "COMMENT: Cannot make fasta file" >> $DEBUG
    echo "PDB-REDO,$PDBID"                 >> $DEBUG
  else if (-e $WORKDIR/user.fasta) then
    #There is a fasta file check if there were sequence conflicts
    if (`grep -c 'mis-matched' $WORKDIR/pdb2fasta.log` > 0) then
      #Show the sequince conflicts
      echo " " | tee -a $LOG
      echo "WARNING!" | tee -a $LOG
      echo "--------" | tee -a $LOG
      echo "Conflict(s) between the input sequence and the input model detected." | tee -a $LOG
      echo "Please, check whether these conflicts can be resolved." | tee -a $LOG
      echo " " | tee -a $LOG
      echo "Details from pdb2fasta output:" | tee -a $LOG
      grep 'mis-matched' $WORKDIR/pdb2fasta.log | cut -c 10- | tee -a $LOG
      echo "--------" | tee -a $LOG
      echo " " | tee -a $LOG
    endif
    #Check whether there was a dodgy input file
    if (`grep -c 'could not be matched to any FASTA input sequence' $WORKDIR/pdb2fasta.log` > 0) then
      set TRUSTSEQ = 1
    endif
  else  
    #Check whether there was a dodgy input file
    if (`grep -c '^SEQRES' $WORKDIR/cache.pdb` == 0) then
      set TRUSTSEQ = 1
    endif
  endif
  
  #Find PDB entries with homologous sequences. PROGRAM: BLASTp
  if ($?BLASTP) then
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.BLASTp.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
    $BLASTP \
    -evalue 0.001 \
    -num_descriptions 9999 \
    -num_alignments 9999 \
    -query $WORKDIR/$PDBID.fasta \
    -out $WORKDIR/$PDBID.blast \
    -db $TOOLS/pdbredo_seqdb.txt >& $WORKDIR/blast.log
  else
    echo " o BlastP is missing or not configured properly" | tee -a $LOG
    echo " o Cannot generate homology-based restraints"    | tee -a $LOG
    echo "homology: BlastP unavailable" >> $DEBUG
    echo "PDB-REDO,$PDBID"              >> $DEBUG
  endif
endif

#Check to see whether a suitable set of experimental sigmas exists (checked now because DYEAR is needed)
#Only give error message is the structure was deposited this century
if ($DYEAR > 2009 && `grep -c "sigma values" $WORKDIR/${PDBID}c2c.log` != 0 ) then
  if (`grep -c "No experimental sigmas found" $WORKDIR/${PDBID}c2c.log` != 0 && $USIGMA == 1) then
    echo "COMMENT: cif2cif: no expertimental sigmas found"   >> $DEBUG
    echo "PDB-REDO,$PDBID"                                   >> $DEBUG
  endif
endif

#Debug message for missing R-free set
if ($C2CERR == 1 && $DYEAR > 2009 && `grep -c useful $WORKDIR/${PDBID}c2c.log` == 0) then
  if (`grep -c useful $WORKDIR/${PDBID}c2c.log` == 0) then
    echo "COMMENT: cif2cif: cannot use _refln.status column" >> $DEBUG
    echo "PDB-REDO,$PDBID"                                   >> $DEBUG
  endif
endif

#Set mask parameter
set MASKPAR  = "solvent vdwprobe $VDWPROBE ionprobe $IONPROBE rshrink $RSHRINK"

#Check for freakishly high resolution
if (`echo $DATARESH | awk '{if ($1 < 0.30) {print "1"} else {print "0"}}'` == 1) then
  echo " o Unlikely high resolution. Reflection data may be corrupted. Cannot continue." | tee -a $LOG
  echo "COMMENT: Suspiciously high data resolution" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                            >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#Set resolution cut-offs if needed
if (`echo $RESOLUTION $DATARESH | awk '{if ($1 - $2 > 0.10) {print "1"} else {print "0"}}'` == 1) then
  set REFIRES = "RESO $RESOLUTION"
  set URESO = $RESOLUTION
else
  set URESO = $DATARESH
endif

#Check for missing high resolution data
if (`echo $DATARESH $RESOLUTION | awk '{if ($1 - $2 > 0.10) {print "1"} else {print "0"}}'` == 1) then
  if ($LOCAL == 1) then
    echo " " 
    echo "WARNING!" | tee -a $LOG
    echo "--------" | tee -a $LOG
    echo "The resolution of the model is higher than that of the data." | tee -a $LOG
    echo "Please, check whether the correct reflection data file was used." | tee -a $LOG
  else
    echo " o High resolution data is missing, calculated R-factors are unrealiable" | tee -a $LOG
  endif
  #Write debug warning
  echo "COMMENT: High resolution data missing" >> $DEBUG
  echo "PDB-REDO,$PDBID"                       >> $DEBUG
  
  #Compensate
  if ($?COMMENT) then
    set COMMENT = "$COMMENT; High resolution data is missing" 
  else
    set COMMENT = "High resolution data is missing" 
  endif
  set NEWMODEL = "-f" #Force rerefinement to finish with new model
endif

#Set legacy mode for PDB entries from the seventies and eighties and ED structures predating 1995
if ($DYEAR < 1990 && $ISXRAY == 1) then
  set LEGACY = 1
  echo " o Warning pre-1990 X-ray PDB entry. It will be treated as legacy entry." | tee -a $LOG
  echo "COMMENT: Run in legacy mode" >> $DEBUG
  echo "PDB-REDO,$PDBID"             >> $DEBUG
else if ($DYEAR < 1995 && $ISED == 1) then
  set LEGACY = 1
  echo " o Warning pre-1995 electron diffraction PDB entry. It will be treated as legacy entry." | tee -a $LOG
  echo "COMMENT: Run in legacy mode" >> $DEBUG
  echo "PDB-REDO,$PDBID"             >> $DEBUG
endif

#Stop and give WHY_NOT comment if no R-factor can be extracted from the PDB header...
if ($RFACT == 0.9990) then
  # ... unless PDB-REDO is running in legacy mode or working on a local file
  echo " o Cannot extract R-factor from PDB header. Recalculated value will be used as refinement target." | tee -a $LOG
  set LEGACY = 1
  echo " o Switching to legacy mode." | tee -a $LOG
endif

#Obtain solvent model from REMARK records
if (`grep '^REMARK   3' $WORKDIR/cache.pdb | grep 'METHOD USED' | grep -c -E 'BABINET|SWAT|BULK|MOEWS|KRETSINGER|TNT'` != 0) then
  set SOLVENT = BULK
endif
#Fall back to a simple model when the keywords 'FLAT' or 'CNS BULK' are found.
if (`grep '^REMARK   3' $WORKDIR/cache.pdb | grep 'METHOD USED' | grep -c -E 'MASK|FLAT|CNS BULK'` != 0) then
  set SOLVENT = SIMP
endif

############################# Create MTZ file for structure factor handling in Refmac ####################################

echo "-Creating MTZ file" | tee -a $LOG

mtzmaking:

#Import CIF file (using all reflections). PROGRAM: cif2mtz
cif2mtz \
HKLIN  $WORKDIR/$PDBID.cif \
HKLOUT $WORKDIR/raw.mtz \
<<eof >> $WORKDIR/mtz_creation.log
  END
eof
if ($status) then
  echo " o Error using CIF2MTZ. Cannot continue." | tee -a $LOG
  echo "COMMENT: cif2mtz: general error" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                 >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  exit(1)
endif

#Remove the phase columns if there is anomalous data and we do not need to do phased refinement
if (`grep -c '_refln.pdbx_F_plus' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
  #We have anomalous data, do we also have phases?
  if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.phase_meas' $WORKDIR/$PDBID.cif` != 0) then
    #Do we need phased refinement?
    if (`echo $PHASES | cut -c 1` == "") then
      #Make a back-up
      cp $WORKDIR/raw.mtz $WORKDIR/raw_phase.mtz

      #Which labels must be removed
      if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0 && `grep -c '_refln.phase_meas' $WORKDIR/$PDBID.cif` != 0) then
        set DELLABEL = "HLA HLB HLC HLD PHIB FOM"
      else if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0) then
        set DELLABEL = "HLA HLB HLC HLD"
      else
        set DELLABEL = "PHIB FOM"
      endif

      #Strip out the phase columns. PROGRAM: mtzutils
      mtzutils \
      HKLIN  $WORKDIR/raw_phase.mtz \
      HKLOUT $WORKDIR/raw.mtz \
<<eof >> $WORKDIR/mtz_creation.log
        EXCLUDE $DELLABEL
        END
eof

    endif
  endif
endif

#Cut the resolution if needed
if ($?MRESO) then
  echo " o Cutting the data at $MRESO A" | tee -a $LOG
  set DATARESH = $MRESO
  set URESO    = $DATARESH

  cp $WORKDIR/raw.mtz $WORKDIR/raw_uncut.mtz
  #Run MTZUTILS to cut the data.
  mtzutils \
  HKLIN  $WORKDIR/raw_uncut.mtz \
  HKLOUT $WORKDIR/raw.mtz \
  <<eof >>$WORKDIR/mtz_creation.log
    RESOLUTION $DATARESH $DATARESL
    END
eof
endif

#Add the wavelength if possible
#Returnpoint
wavelengthadd:

if (`grep -c '_diffrn_radiation_wavelength.wavelength' $WORKDIR/$PDBID.cif` != 0) then
  #Get the wavelenghth
  set WAVELENGTH = `grep '_diffrn_radiation_wavelength.wavelength' $WORKDIR/$PDBID.cif | awk '{print $2}'`
else if ($WAVELPDB != '0.00000') then
  set WAVELENGTH = $WAVELPDB
else
  set WAVELENGTH = 'NA'
  if ($DYEAR > 2006) then
    echo "COMMENT: wavelength unknown" >> $DEBUG
    echo "PDB-REDO,$PDBID"             >> $DEBUG
  endif
endif

#Add the wavelength to the MTZ file
if ($WAVELENGTH != 'NA') then
  #Make a back-up
  cp $WORKDIR/raw.mtz $WORKDIR/raw_nowavel.mtz

  #Run CAD. PROGRAM: cad
  cad \
  HKLIN1 $WORKDIR/raw_nowavel.mtz \
  HKLOUT $WORKDIR/raw.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    LABIN FILE 1  ALLIN
    DWAVELENGTH FILE_NUMBER 1 1 $WAVELENGTH
    END
eof
  #Check for a large reduction in the number of reflections. This indicates unmerged reflections.
  set PREFIN  = `grep 'and a total of' $WORKDIR/mtz_creation.log | tail -n 1 | awk '{print $5}'`
  set PREFOUT = `grep 'Final Total of Unique records to HKLOUT' $WORKDIR/mtz_creation.log | tail -n 1 | cut -d '=' -f 2`
  if (`echo $PREFIN $PREFOUT | awk '{if ($1 / $2 > 1.9) {print "1"} else {print "0"}}'` == 1) then
    #There seem to be unmerged reflections or Friedel pairs on separate lines or many systematically absent reflections
    if (`grep -c 'Systematic absent reflection rejected' $WORKDIR/mtz_creation.log` > 1000) then
      #Something seroiusly wrong with the dataset
      if ($LOCAL == 1) then
        #Give the long error message
        echo " " | tee -a $LOG
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "There are too many systematic absent reflections. Cannot continue. Please, check the input reflection data." | tee -a $LOG
      else
        echo " o There are too many systematic absent reflections. Cannot continue." | tee -a $LOG
      endif
      echo "COMMENT: Too many systematic absent reflections" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                 >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      exit(1)
    endif

    #Report only once
    if ($?MERGEREP) then
      #Problem already reported
    else
      echo " o The reflections are not completely merged" | tee -a $LOG
      echo "COMMENT: unmerged reflections" >> $DEBUG
      echo "PDB-REDO,$PDBID"               >> $DEBUG
      set MERGEREP = 1
    endif

    #First get the sigF or sigI column if it was rejected before
    if (-e $WORKDIR/${PDBID}_badsig.cif) then
      #Report
      echo "   * Recovering sigma value data before merging" | tee -a $LOG

      #Take the reflection file with the sigma columns
      mv $WORKDIR/${PDBID}_badsig.cif $WORKDIR/$PDBID.cif

      #Remake the mtz file
      goto mtzmaking
    endif

    #Delete the output from cad.
    rm $WORKDIR/raw.mtz

    #Merge the reflections. PROGRAM: sftools
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.sftools.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    sftools << eof >> $WORKDIR/mtz_creation.log
    read $WORKDIR/raw_nowavel.mtz
    reduce
    merge average
    write $WORKDIR/raw.mtz
    stop
eof

    #Add the wavelength again
    goto wavelengthadd
  endif
endif

#create a unique list of reflections for the given unit cell-symmetry-resolution. 
unique \
HKLOUT $WORKDIR/unique.mtz \
<<eof >> $WORKDIR/mtz_creation.log
  CELL $AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA
  SYMM '$SPACEGROUP'
  LABOUT F=FUNI SIGF=SIGFUNI
  RESOLUTION $DATARESH
  END
eof

refree:

#Combine data (generate R-free set if one doesn't already exist)
if ($GOTR == 0) then

  #Calculate the required R-free fraction (between 0.05 and 0.10)
  set FRAC = `echo $REFCNT | awk '{if ($1 > 20000) {print "0.05"} else if ($1 < 10000) {print "0.10"} else {printf ("%.4f\n", 1000/$1)}}'`

  #Generate R-free flag, then ....
  echo " o Generating R-free set using $FRAC fraction of all reflections" | tee -a $LOG

  #PROGRAM: freerflag
  freerflag \
  HKLIN $WORKDIR/unique.mtz \
  HKLOUT $WORKDIR/freer.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    FREERFRAC $FRAC
    END
eof

  #.... combine data
  cad \
  HKLIN2 $WORKDIR/freer.mtz \
  HKLIN1 $WORKDIR/raw.mtz \
  HKLOUT $WORKDIR/combined.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    LABIN FILE 2  E1=FreeR_flag
    LABIN FILE 1  ALLIN
    END
eof

  #Rename FreeR_flag to FREE
  mtzutils \
  HKLIN  $WORKDIR/combined.mtz \
  HKLOUT $WORKDIR/merged.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    COLUMN_LABELS FreeR_flag=FREE
    END
eof

  #Clean up extra temporary files
  rm $WORKDIR/freer.mtz

else

  #Just combine data
  cad \
  HKLIN2 $WORKDIR/unique.mtz \
  HKLIN1 $WORKDIR/raw.mtz \
  HKLOUT $WORKDIR/combined.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    LABIN FILE 2  ALLIN
    LABIN FILE 1  ALLIN
    END
eof

  #Strip out the FUNI and SIGFUNI columns
  mtzutils \
  HKLIN  $WORKDIR/combined.mtz \
  HKLOUT $WORKDIR/merged.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    EXCLUDE FUNI SIGFUNI
    END
eof
endif

#Clean up
rm $WORKDIR/combined.mtz

#Need to work with intensities?
if (`mtzdmp $WORKDIR/merged.mtz -e | grep -A 2 'Column Types' | grep -c J` != 0) then
  echo " o Using ctruncate to convert intensities to amplitudes" | tee -a $LOG
  set UTRUNCATE = 1
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.ctruncate.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  #Create a back-up
  cp $WORKDIR/merged.mtz $WORKDIR/mergedbu.mtz

  if (`grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
    #Convert I to F, also for anomalous data. PROGRAM: ctruncate
    ctruncate \
    -mtzin  $WORKDIR/merged.mtz \
    -mtzout $WORKDIR/ctruncate.mtz \
    -no-aniso \
    -freein "/*/*/[FREE]" \
    -colano "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]" \
    -colin  "/*/*/[I,SIGI]" >>& $WORKDIR/mtz_creation.log
  else
    #Convert I to F. PROGRAM: ctruncate
    ctruncate \
    -mtzin  $WORKDIR/merged.mtz \
    -mtzout $WORKDIR/ctruncate.mtz \
    -no-aniso \
    -freein "/*/*/[FREE]" \
    -colin "/*/*/[I,SIGI]" >>& $WORKDIR/mtz_creation.log
  endif

  if ($status || ! -e $WORKDIR/ctruncate.mtz) then
    #Try running without using ctruncate
    echo "   * Error using ctruncate" | tee -a $LOG
    echo "   * Using cif2cif to convert intensities" | tee -a $LOG

    #Make backup
    cp $WORKDIR/$PDBID.cif $WORKDIR/${PDBID}_notruncate.cif

    #Run cif2cif
    $TOOLS/cif2cif -c $C2CANO $FEWREFS $INTENS $USTATUS $SIGMA \
    $WORKDIR/${PDBID}_notruncate.cif \
    $WORKDIR/$PDBID.cif \
    >> $WORKDIR/mtz_creation.log

    #Check for the existence of a test set.    
    set GOTR = `grep -c _refln.status $WORKDIR/$PDBID.cif`
    if ($GOTR == 0) then
      #Check the reason for the missing set
      if (`grep -c 'Warning: too small R-free set.' $WORKDIR/mtz_creation.log` != 0) then
        echo "   * Original R-free set too small, a new one will be created" | tee -a $LOG
      else
        echo "   * No R-free set defined, a new one will be created" | tee -a $LOG
      endif
    endif
    
    #Restart making the MTZ file
    goto mtzmaking
  endif

  #Rename amplitude columns and remove intensity columns
  if (`grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
    set DELLABEL = 'DANO SIGDANO ISYM I(+) SIGI(+) I(-) SIGI(-)'
  else
    set DELLABEL = 'I SIGI'
  endif

  mtzutils \
  HKLIN  $WORKDIR/ctruncate.mtz \
  HKLOUT $WORKDIR/merged.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    COLUMN_LABELS F=FP SIGF=SIGFP
    EXCLUDE $DELLABEL
    END
eof

  #Clean up
  rm $WORKDIR/ctruncate.mtz
endif


#Clean up data and create output file (use separate log file)
freerflag \
HKLIN $WORKDIR/merged.mtz \
HKLOUT $WORKDIR/$PDBID.mtz \
<<eof >& $WORKDIR/freerflag.log
  COMPLETE FREE=FREE
  END
eof
if ($status) then
  #Problem with FREERFLAG
  echo " o Problem with test set, creating a new one." | tee -a $LOG

  #Give debug message
  echo "COMMENT: freerflag: test set problem" >> $DEBUG
  echo "PDB-REDO,$PDBID"                      >> $DEBUG

  #Strip old test set
  mtzutils \
  HKLIN  $WORKDIR/merged.mtz \
  HKLOUT $WORKDIR/nofree.mtz \
<<eof >> $WORKDIR/freerflag.log
    EXCLUDE FREE
    END
eof

  #Calculate the required R-free fraction (between 0.05 and 0.10)
  set FRAC = `echo $REFCNT | awk '{if ($1 > 20000) {print "0.05"} else if ($1 < 10000) {print "0.10"} else {printf ("%.4f\n", 1000/$1)}}'`

  #Generate R-free flag, then ....
  echo "   * Generating R-free set using $FRAC fraction of all reflections" | tee -a $LOG

  freerflag \
  HKLIN  $WORKDIR/nofree.mtz \
  HKLOUT $WORKDIR/newfree.mtz \
<<eof >> $WORKDIR/freerflag.log
    FREERFRAC $FRAC
    END
eof

  #Rename FreeR_flag to FREE
  mtzutils \
  HKLIN  $WORKDIR/newfree.mtz \
  HKLOUT $WORKDIR/$PDBID.mtz \
<<eof >> $WORKDIR/freerflag.log
    COLUMN_LABELS FreeR_flag=FREE
    END
eof

endif

#Append freerflag.log to mtz_creation.log
cat $WORKDIR/freerflag.log >> $WORKDIR/mtz_creation.log
rm  $WORKDIR/freerflag.log


#Create lowres mtz file for sfcheck and mtz2various
if (`echo $RESOLUTION $DATARESH | awk '{if ($1 - $2 > 0.10) {print "1"} else {print "0"}}'` == 1) then
  #Run MTZUTILS to cut the data
  mtzutils \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/lowres.mtz \
  <<eof >> $WORKDIR/mtz_creation.log
    RESOLUTION $RESOLUTION $DATARESL
    END
eof
  if ($status) then
    echo "   * Error using MTZUTILS. Cannot continue." | tee -a $LOG
    echo "COMMENT: mtzutils: general error" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                  >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif
else
  cp $WORKDIR/$PDBID.mtz $WORKDIR/lowres.mtz
endif

#Maximise the used CPU for sfcheck
limit cputime 1h

#Calculate Wilson B, completeness and twinning statistics. PROGRAM: SFcheck
sfcheck \
-f $WORKDIR/lowres.mtz \
-mem 150 \
>> $WORKDIR/mtz_creation.log

#Reset the limit for sfcheck
limit cputime unlimited

#See if sfcheck ran propely and truncate was used
if (`grep -c "ERR: NUMBER OF REFLS" sfcheck.log` != 0 && $UTRUNCATE == 1) then
  #Try running without using ctruncate
  echo "   * Error using ctruncate" | tee -a $LOG
  echo "   * Using cif2cif to convert intensities" | tee -a $LOG

  #Make backup
  cp $WORKDIR/$PDBID.cif $WORKDIR/${PDBID}_notruncate.cif

  #Run cif2cif
  $TOOLS/cif2cif -c $C2CANO $FEWREFS $INTENS $USTATUS $SIGMA \
  $WORKDIR/${PDBID}_notruncate.cif \
  $WORKDIR/$PDBID.cif \
  >> $WORKDIR/mtz_creation.log
  
  #Check for the existence of a test set.    
  set GOTR = `grep -c _refln.status $WORKDIR/$PDBID.cif`
  if ($GOTR == 0) then
    #Check the reason for the missing set
    if (`grep -c 'Warning: too small R-free set.' $WORKDIR/mtz_creation.log` != 0) then
      echo "   * Original R-free set too small, a new one will be created" | tee -a $LOG
    else
      echo "   * No R-free set defined, a new one will be created" | tee -a $LOG
    endif
  endif

  #Restart making the MTZ file
  goto mtzmaking
endif

#Check the SFCHECK logfile
if (`grep -c 'by Wilson' sfcheck.log` > 0 && `grep -c 'Completeness :' sfcheck.log` > 0) then

  #Mine sfcheck output
  set BWILS     = `grep 'by Wilson' sfcheck.log | awk '{print $7}'`
  set COMPLETED = `grep 'Completeness :' sfcheck.log | cut -c 16- | awk '{print $1}'`
  set TWINA     = `grep 'Alpha(twin fraction)' sfcheck.log | cut -c 38-43`

else
  #Contingency for SFCHECK probles
  echo "COMMENT: SFCHECK: general error" >> $DEBUG
  echo "PDB-REDO,$PDBID"                 >> $DEBUG
  
  set BWILS     = $BAVER
  set TWINA     = ""
  set COMPLETED = `mtzdmp $WORKDIR/lowres.mtz | grep '  FP' | awk '{print $6}'`
endif
  
#Check for twinning with Phaser if SFCHECK suggests twinning. PROGRAM: PHASER
if ($TWINA != "") then
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.PHASER.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json

  phaser \
  << eof > $WORKDIR/phaser.log
    MODE NCS
    ROOT $WORKDIR
    HKLIN $WORKDIR/lowres.mtz
    LABIN  F=FP SIGF=SIGFP
    COMPOSITION BY SOLVENT
    COMPOSITION PERCENTAGE $SOLVD
    HKLOUT OFF
    END
eof
  if (`grep -c "EXIT STATUS: SUCCESS" $WORKDIR/phaser.log` == 0) then
    #Phaser problem.
    echo "COMMENT: phaser: general error" >> $DEBUG
    echo "PDB-REDO,$PDBID"                >> $DEBUG
    set PHASTWIN = 0
  else
    #Check whether PHASER found twinning
    if (`grep -c "Warning: Intensity moments suggest significant twinning" $WORKDIR/phaser.log` == 0 && `grep -c "Warning: Intensity moments suggest possibility of twinning" $WORKDIR/phaser.log` == 0) then
      #Phaser did not find twinning
      set PHASTWIN = 0
    else
      echo " o Both SFCHECK and PHASER suggest the data are twinned" | tee -a $LOG
      set PHASTWIN = 1
    endif
  endif
else
  #Don't look for twinning in Phaser
  set PHASTWIN = 0
endif


#Get the correct number of (R-free) reflections
mtz2various \
HKLIN $WORKDIR/lowres.mtz \
HKLOUT $WORKDIR/temp.hkl \
<<eof >> $WORKDIR/mtz_creation.log
  OUTPUT CIF data_temp
  LABIN FP=FP SIGFP=SIGFP FREE=FREE
  END
eof

set NTSTCNT = `grep -cE ' f ' $WORKDIR/temp.hkl`
set NREFCNT = `grep -cE ' [of] ' $WORKDIR/temp.hkl`
set WORKCNT = `grep -cE ' o ' $WORKDIR/temp.hkl`

#Repeat process if test set was merged away
if ($?MERGEREP && $NTSTCNT < 30) then
  set GOTR = 0
  
  cp $WORKDIR/raw.mtz $WORKDIR/rawbu.mtz
  #Strip old test set
  mtzutils \
  HKLIN  $WORKDIR/rawbu.mtz \
  HKLOUT $WORKDIR/raw.mtz \
<<eof >> $WORKDIR/freerflag.log
    EXCLUDE FREE
    END
eof

  #Go back to creating the test set
  goto refree
endif

#Check for anomalous data
if (`grep -c '_refln.pdbx_F_plus' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
  set ANOMCOEF = 'F+=F(+) SIGF+=SIGF(+) F-=F(-) SIGF-=SIGF(-)'
  set ANOMCMD  = 'ANOM maponly'
endif


#Cleanup
rm sfcheck.xml
if (-e sfcheck_XXXX.ps) then
  rm sfcheck_XXXX.ps
endif  
rm $WORKDIR/temp.hkl
rm $WORKDIR/lowres.mtz


###################################### Report structure model and data details ###########################################

#Print values to screen. Compensate if R(-free) was not reported.
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Structure model and data details ******" | tee -a $LOG
echo "Cell axes            : $AAXIS $BAXIS $CAXIS" | tee -a $LOG
echo "Cell angles          : $ALPHA $BETA $GAMMA"  | tee -a $LOG
echo "Space group          : $SPACEGROUP" | tee -a $LOG
echo "Resolution (pdb)     : $RESOLUTION" | tee -a $LOG
echo "Resolution (data)    : $DATARESH"   | tee -a $LOG  #Data resolution (only reflections with F/SIGF>1.0 are used)
echo "Lowest resolution    : $DATARESL"   | tee -a $LOG #Lowest resolution in the dataset (only reflections with F/SIGF>1.0 are used)
if ($RFACT == 0.9990) then
  echo "Reported Rfactor     : none"   | tee -a $LOG
  set RFACT = "NA"
  set RHEAD = 0
else
  echo "Reported Rfactor     : $RFACT" | tee -a $LOG
  set RHEAD = 1
endif
if ($RFREE == 0.9990) then
  echo "Reported R-free      : none" | tee -a $LOG
  set RFHEAD = 0
  set RFREE = "NA"
else
  echo "Reported R-free      : $RFREE" | tee -a $LOG
  set RFHEAD = 1
endif
echo "Average B-factor     : $BAVER"    | tee -a $LOG
echo "Wilson B from SFcheck: $BWILS"    | tee -a $LOG
echo "Solvent percentage   : $SOLVD"    | tee -a $LOG
echo "Solvent model        : $SOLVENT"  | tee -a $LOG
echo "VdW probe            : $VDWPROBE" | tee -a $LOG
echo "Ion probe            : $IONPROBE" | tee -a $LOG
echo "Shrinkage            : $RSHRINK"  | tee -a $LOG
echo "Reflections (data)   : $REFCNT"   | tee -a $LOG
echo "Work set size        : $WORKCNT"  | tee -a $LOG
echo "Test set size (data) : $TSTCNT ($TSTPRC%)" | tee -a $LOG
echo "Test set size (used) : $NTSTCNT"   | tee -a $LOG
echo "Completeness (pdb)   : $COMPLETEH" | tee -a $LOG
echo "Completeness (used)  : $COMPLETED" | tee -a $LOG
echo "Wavelength           : $WAVELENGTH"| tee -a $LOG
if ($TWINA != "") then
  echo "Twin fraction alpha  : $TWINA"   | tee -a $LOG
endif
echo "Refinement tool      : $PROG"      | tee -a $LOG
echo "Deposition year      : $DYEAR"     | tee -a $LOG


############################################# Recalculation of R and R-free ##############################################
molrepped:

echo " " | tee -a $LOG
echo " " | tee -a $LOG

#Stop if the completeness is too low (less than 20% or less than two thirds of reported completeness), but not if paired refinement is forced
if ($FORCEPAIRED == 0 && `echo $COMPLETED $COMPLETEH | awk '{if ($1 < 20.0) {DIFF = ($1 - $2); if (DIFF < 0) {DIFF = -1*DIFF}; if (DIFF < 2.5) {print "0"} else {print "1"}} else if ($1 < 0.67*$2) {print "1"} else {print "0"}}'` == 1) then
  #If this is a PDB entry check PDBpeep to see if the reported completeness is the ellipsoidal value
  if ($LOCAL == 0) then
    #Get the ellipsoidal completeness from PDBpeep
    $WEBGET -q http://staraniso.globalphasing.org/PDB/${PDBID}-1.log -O $WORKDIR/PDBpeep.log
    set ELLCOMP = `grep 'All data' $WORKDIR/PDBpeep.log | head -n 1 | awk '{print 100*$14}'`
    if (`echo $ELLCOMP $COMPLETEH | awk '{DIFF = ($1 - $2); if (DIFF < 0) {DIFF = -1*DIFF}; if (DIFF < 5.0) {print "0"} else {print "1"}}'` == 1) then
      #There a real problem with the dataset completeness
      set COMPERROR = 1
    endif
  else 
    #The completeness is suspicious
    set COMPERROR = 1
  endif
endif  
    
if ($COMPERROR == 1) then    
  echo "-These data are much too incomplete to use" | tee -a $LOG
  echo " o Cannot continue" | tee -a $LOG
  echo "COMMENT: Too much missing experimental data" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                             >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

echo "****** R(-free) recalculation ******" | tee -a $LOG

#NCS settings
if ($DONCS == 1) then
  set NCSTYPE  = "ncsr local"
  set NCSALIGN = "ncsr align level 0.90 iterate N rmslevel 2.00"
  set NCSNEIGH = "ncsr neighbours exclude"
endif

#Import the additional restraints
if ("$INEXT" != "") then
  echo "-Importing additional restraints" | tee -a $LOG
  #Check whether the file is not a cif restraint file or a PDB file or a TLS definition
  if (`grep -a -c '_chem_comp_bond.value_dist' $INEXT` == 0 && `grep -a -c '_chem_link_bond.value_dist' $INEXT` == 0 && `grep -a -c '^[AH][TE][OT][MA]' $INEXT` == 0 && `grep -a -c '^RANGE' $INEXT` == 0 && `grep -c 'exte' $INEXT` != 0) then
    #Set up the external restraint file
    echo "external UndefinedAtoms ignore"  > $WORKDIR/external.rest
    echo "external weight scale $EXTSCAL" >> $WORKDIR/external.rest
    sed 's/  / /g' $INEXT >> $WORKDIR/external.rest
    #Set-up Refmac to use the uploaded additional restraints
    set RESTCMD = "@$WORKDIR/external.rest"
  else
    #This is not the right type of restraint file. Report...
    echo " o The additional restraint file is not in REFMAC format; it will be ignored" | tee -a $LOG
    echo "COMMENT: Corrupt external restraints" >> $DEBUG
    echo "PDB-REDO,$PDBID"                      >> $DEBUG
  endif
endif


#Set up scaling function
set SCALING = "lssc function a sigma $EXPSIG"

#Restarting oint with previous PDB-REDO entry
previousredo:

#TLS settings
junkedtls:
set ORITLS = 0
set ISTLS  = 'notls'
set TLSLIN =     #Do not use static TLS tensors unless...
#.... a TLS exists...
if (-e $WORKDIR/$PDBID.tls && $DOTLS == 1) then
  #...and has complete TLS tensors.
  if (`grep -c ORIGIN $WORKDIR/$PDBID.tls` != 0) then
    set ORITLS = 1
    set ISTLS  = 'tls'
    set TLSLIN = `echo TLSIN $WORKDIR/$PDBID.tls` #Use static TLS-groups for better R/R-free calculation
    echo -n "-Recalculating R-factors with original TLS tensors    " | tee -a $LOG
  else
    echo -n "-Running Refmac for 0 cycles "    | tee -a $LOG
  endif
else
  echo -n "-Running Refmac for 0 cycles " | tee -a $LOG
endif

#Run refmac for 0 cycles to check R-factors (with TLS). PROGRAM: REFMAC
refmac5 \
XYZIN  $WORKDIR/${PDBID}_prepped.pdb \
XYZOUT $WORKDIR/${PDBID}_0cyc$ISTLS.pdb \
HKLIN  $WORKDIR/$PDBID.mtz \
HKLOUT $WORKDIR/${PDBID}_0cyc$ISTLS.mtz \
$LIBLIN \
$TLSLIN \
$SCATLIN \
<<eof >& $WORKDIR/${PDBID}_0cyc$ISTLS.log
  $SCATTERCMD
  make check NONE
  make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
    ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
  refi type REST resi MLKF meth CGMAT bref MIXE
  $REFIRES
  ncyc 0
  tlsd waters exclude
  scal type $SOLVENT $SCALING
  solvent YES
  $MASKPAR
  $LOWMEM
  weight $WGTSIG MATRIX 0.5
  monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
    chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
  $NCSTYPE
  $NCSALIGN
  $NCSNEIGH
  $NCSSTRICT
  labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
  $ANOMCMD
  pdbout copy remarks 200 280 350
  pdbout copy expdta 
  NOHARVEST
  END
eof

if ($status || `grep -c 'Error: Fatal error. Cannot continue' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
  echo " " | tee -a $LOG
  #Check whether there are too many TLS goups
  if (`grep -c 'Too many tls groups. Maximum allowed is' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
    #Give error message
    echo " o Input model had more TLS groups than Refmac can handle. They will be ignored." | tee -a $LOG
    echo "COMMENT: too many TLS groups" >> $DEBUG
    echo "PDB-REDO,$PDBID"              >> $DEBUG  
    
    #Recover by ignoring te original TLS model
    mv $WORKDIR/$PDBID.tls $WORKDIR/$PDBID.notls
    goto junkedtls
  endif
  #Check to see if there is a problem with alternate residues
  if (`grep -c 'different residues have the same number' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
    if ($LOCAL == 1) then
      #Give the long error message
      echo " " | tee -a $LOG
      echo "FATAL ERROR!" | tee -a $LOG
      echo "------------" | tee -a $LOG
      echo "Refmac had problems using these residues with alternate identities:" | tee -a $LOG
      grep 'ERROR:' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 12- | tee -a $LOG
      echo "This problem may be solved by renumbering residues or by ensuring that" | tee -a $LOG
      echo "alternate atoms directly follow eachother in the PDB file." | tee -a $LOG
      echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details."  | tee -a $LOG
    else
      #Give the simple error message
      echo " o Cannot use structure with alternate residues" | tee -a $LOG
    endif
    #Write out WHY_NOT mesage
    echo "COMMENT: refmac: error with alternate residues" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
  #Check to see if there are dictionary problems
  if (! -e $WORKDIR/${PDBID}_het.cif || "$INREST" != "") then
    #Check for atoms not in the dictionary
    if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then

      #Give PDB-REDO mode-specific error messages
      if ($LOCAL == 1) then
        #Give the long error message
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Some residues have (atoms with) naming conflicts." | tee -a $LOG
        echo "Please, use standard atom and residue names in your input PDB file or upload a custom restraint file." | tee -a $LOG
        echo " " | tee -a $LOG
        echo "Residue  PDB standard description" | tee -a $LOG
        foreach HETID (`grep 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 22-37 | sort -u`)
          set HETID1 = `echo $HETID |cut -c 1-1`
          #Is the entry in LigandExpo
          wget --spider -q http://ligand-expo.rcsb.org/reports/$HETID1/$HETID
          if ($status) then
            echo "$HETID      New compound: make sure all $HETID residues are consistent within the input PDB" | tee -a $LOG
          else
            echo "$HETID      http://ligand-expo.rcsb.org/reports/$HETID1/$HETID" | tee -a $LOG
          endif
        end

        #Give server-specific extra information
        if ($SERVER == 1) then
          echo "Details from Refmac output:" | tee -a $LOG
          grep ' ERROR :' $WORKDIR/${PDBID}_0cyc$ISTLS.log | tee -a $LOG
        else
          echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details." | tee -a $LOG
        endif
      else
        #Give the short error message
        echo " o Residue or atom naming conflict. Cannot continue." | tee -a $LOG
      endif

      #Give WHY_NOT comment
      echo "COMMENT: refmac: residue or atom naming conflict" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                  >> $WHYNOT

    #Check for problems making a compound description
    else if (`grep -a -c 'is not completely connected' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then

      #Give PDB_RED mode-specific error messages
      if ($LOCAL == 1) then
        #Give the long error message
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Cannot create restraint file for residue:" | tee -a $LOG
        grep 'program will create complete description for' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -d ':' -f 2 | sort -u | tee -a $LOG
        echo "Please, supply a restraint file using '--restin=Your_restraints.cif'." | tee -a $LOG
      else
        #Give the short error message
        echo " o Cannot create restraint file. Cannot continue." | tee -a $LOG
      endif

      #Write WHY_NOT comment
      echo "COMMENT: refmac: cannot create restraint file" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                               >> $WHYNOT

    #Check for NCS alignment problems
    else if (`grep -a -c 'ncs_ncs_generate.f90' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0 && ! -e $WORKDIR/renumber.json) then
        
      #Renumber the PDB file and start again
      cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.rnbterror.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
      
      cp $WORKDIR/${PDBID}_prepped.pdb $WORKDIR/${PDBID}_prepped.bak
      $TOOLS/rnbterror -v \
      -jsonout $WORKDIR/renumber.json \
      -pdbin $WORKDIR/${PDBID}_prepped.bak \
      -pdbout $WORKDIR/${PDBID}_prepped.pdb > $WORKDIR/renumber.log
          
      #Write DEBUG message
      echo " o Problem with NCS alignment. Renumbering terminal residues and restarting." | tee -a $LOG
      echo "COMMENT: residues renumbered" >> $DEBUG
      echo "PDB-REDO,$PDBID"              >> $DEBUG  
     
     #Start again
      goto renumbered

    #All other problems
    else
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
    endif

    #Write out status files
    if ($SERVER == 1) then
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif


#Check to see whether a new ligand was encountered (only if no extra restraints were provided).
if ("$INREST" == "") then
  if(-e $WORKDIR/${PDBID}_het.cif) then
    echo " " | tee -a $LOG
    echo -n " o New ligand encountered, retrying " | tee -a $LOG

    #Backup old files
    cp $WORKDIR/${PDBID}_0cyc$ISTLS.log $WORKDIR/${PDBID}_0cycv1.log
    if(-e $WORKDIR/${PDBID}_0cyc$ISTLS.pdb) then
      cp $WORKDIR/${PDBID}_0cyc$ISTLS.pdb $WORKDIR/${PDBID}_0cycv1.pdb
    endif

    #Start using the new ligand library
    set LIBLIN = `echo LIBIN $WORKDIR/${PDBID}_het.cif LIBOUT $WORKDIR/${PDBID}_het2.cif`

    #Rerun refmac for 0 cycles
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_prepped.pdb \
    XYZOUT $WORKDIR/${PDBID}_0cyc$ISTLS.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc$ISTLS.mtz \
    $LIBLIN \
    $TLSLIN \
    $SCATLIN \
<<eof >& $WORKDIR/${PDBID}_0cyc$ISTLS.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
	ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref MIXE
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
      $ANOMCMD
      pdbout copy remarks 200 280 350
      pdbout copy expdta
      NOHARVEST
      END
eof
    if($status) then
      #Try to give a specific error message
      #Check to see if there is a problem with alternate residues
      if (`grep -c 'different residues have the same number' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
        if ($LOCAL == 1) then
          #Give the long error message
          echo " " | tee -a $LOG
          echo "FATAL ERROR!" | tee -a $LOG
          echo "------------" | tee -a $LOG
          echo "Refmac had problems using these residues with alternate identities:" | tee -a $LOG
          grep 'ERROR:' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 12- | tee -a $LOG
          echo "This problem may be solved by renumbering residues or by ensuring that" | tee -a $LOG
          echo "alternate atoms directly follow eachother in the PDB file." | tee -a $LOG
          echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details."  | tee -a $LOG
        else
          #Give the simple error message
          echo " " | tee -a $LOG
          echo " o Cannot use structure with alternate residues" | tee -a $LOG
        endif
        #Write out WHY_NOT mesage
        echo "COMMENT: refmac: error with alternate residues" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      endif

      #Check to see if there are dictionary problems
      if (! -e $WORKDIR/${PDBID}_het2.cif) then
        #Check for atoms not described in the dictionary
        if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
          if ($LOCAL == 1) then
            #Give the long error message
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Some residues have (atoms with) naming conflicts." | tee -a $LOG
            echo "Please, use standard atom and residue names in your input PDB file or upload a custom restraint file." | tee -a $LOG
            echo " " | tee -a $LOG
            echo "Residue  PDB standard description" | tee -a $LOG
            echo "---------------------------------" | tee -a $LOG
            foreach HETID (`grep 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 22-37 | sort -u`)
              set HETID1 = `echo $HETID |cut -c 1-1`
              #Is the entry in LigandExpo
              wget --spider -q http://ligand-expo.rcsb.org/reports/$HETID1/$HETID
              if ($status) then
                echo "$HETID      New compound: make sure all $HETID residues are consistent within the input PDB" | tee -a $LOG
              else
                echo "$HETID      http://ligand-expo.rcsb.org/reports/$HETID1/$HETID" | tee -a $LOG
              endif
            end
            echo " " | tee -a $LOG

            #Give server-specific extra information
            if ($SERVER == 1) then
              echo "Details from Refmac output:" | tee -a $LOG
              grep ' ERROR :' $WORKDIR/${PDBID}_0cyc$ISTLS.log | tee -a $LOG
            else
              echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details." | tee -a $LOG
            endif
          else
            #Give the short error message
            echo " " | tee -a $LOG
            echo " o Residue or atom naming conflict. Cannot continue." | tee -a $LOG
          endif
          echo "COMMENT: refmac: residue or atom naming conflict" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
        #Check for residues for which a restraint file cannot be generated
        else if (`grep -a -c 'is not completely connected' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
          if ($LOCAL == 1) then
            #Give the long error message
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Cannot create restraint file for residue:" | tee -a $LOG
            grep 'program will create complete description for' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -d ':' -f 2 | sort -u
            echo "Please, supply a restraint file using '--restin=Your_restraints.cif'." | tee -a $LOG
          else
            #Give the short error message
            echo " " | tee -a $LOG
            echo " o Cannot create restraint file. Cannot continue." | tee -a $LOG
          endif
          echo "COMMENT: refmac: cannot create restraint file" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                               >> $WHYNOT
        #Check for problems with local NCS restraints
        else if (`grep -a -c 'ncs_ncs_generate.f90' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0 && ! -e $WORKDIR/renumber.json) then
          #Renumber the PDB file and start again
          cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.rnbterror.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
          cp $WORKDIR/${PDBID}_prepped.pdb $WORKDIR/${PDBID}_prepped.bak
          
          $TOOLS/rnbterror \
          -jsonout $WORKDIR/renumber.json \
          -pdbin $WORKDIR/${PDBID}_prepped.bak \
          -pdbout $WORKDIR/${PDBID}_prepped.pdb > $WORKDIR/renumber.log
          
          #Write DEBUG message
          echo " " | tee -a $LOG
          echo " o Problem with NCS alignment. Renumbering terminal residues and restarting." | tee -a $LOG
          echo "COMMENT: residues renumbered" >> $DEBUG
          echo "PDB-REDO,$PDBID"              >> $DEBUG  
          
          #Start again
          goto renumbered
         
        else
          echo " " | tee -a $LOG
          echo " o Problem with refmac. Cannot continue." | tee -a $LOG
          echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
        endif
      endif
      #Stop the run
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  endif
endif

#Check the value for LOGSTEP and RBLS (only once)
if ($GOTOLD == 0) then
  if (`grep -a -c 'Rms ChirVolume' $WORKDIR/${PDBID}_0cyc$ISTLS.log` == 0) then
    @ LOGSTEP = ($LOGSTEP - 1)
    @ RBLS = ($LOGSTEP - 2)
  else
    @ RBLS = ($LOGSTEP - 3)
  endif
endif  

#Get calculated R-factor
set PRCAL1 = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $2}'`
set TRFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $3}'`
echo "(R-free = $TRFREE)" | tee -a $LOG

#Check R-factor without TLS tensors
if ($ORITLS == 1) then
  echo -n "-Recalculating R-factors without original TLS tensors " | tee -a $LOG
  set ISTLS = 'notls'

  #Rerun refmac for 0 cycles, again
  refmac5 \
  XYZIN  $WORKDIR/${PDBID}_prepped.pdb \
  XYZOUT $WORKDIR/${PDBID}_0cyc$ISTLS.pdb \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/${PDBID}_0cyc$ISTLS.mtz \
  $LIBLIN \
  $SCATLIN \
<<eof > $WORKDIR/${PDBID}_0cyc$ISTLS.log
    $SCATTERCMD
    make check NONE
    make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
      ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
    refi type REST resi MLKF meth CGMAT bref MIXE
    $REFIRES
    ncyc 0
    tlsd waters exclude
    scal type $SOLVENT $SCALING
    solvent YES
    $MASKPAR
    $LOWMEM
    weight $WGTSIG MATRIX 0.5
    monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
      chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
    $NCSTYPE
    $NCSALIGN
    $NCSNEIGH
    $NCSSTRICT
    labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
    $ANOMCMD
    pdbout copy remarks 200 280 350
    pdbout copy expdta
    NOHARVEST
    END
eof
  if ($status) then
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif

  #Get calculated R-factor
  set PRCAL2 = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $2}'`
  set TRFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $3}'`
  echo "(R-free = $TRFREE)" | tee -a $LOG

  #Which is better (i.e. gives lower R-factor), with or without TLS?
  if (`echo $PRCAL1 $PRCAL2 | awk '{if ($1 - $2 < 0) {print "with"} else {print "without"}}'` == 'with') then
    echo " o Recalculating the R-factor with TLS tensors gives the best result" | tee -a $LOG
    set ISTLS = 'tls'
  else
    echo " o Recalculating the R-factor without TLS tensors gives the best result" | tee -a $LOG
    set TLSLIN =   #Do not use static TLS tensors
  endif
endif

#Clean up a bit
mv $WORKDIR/${PDBID}_0cyc$ISTLS.log $WORKDIR/${PDBID}_0cyc.log
mv $WORKDIR/${PDBID}_0cyc$ISTLS.pdb $WORKDIR/${PDBID}_0cyc.pdb
mv $WORKDIR/${PDBID}_0cyc$ISTLS.mtz $WORKDIR/${PDBID}_0cyc.mtz

#Evaluate log file (dirty)
set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`


#Did SFCHECK report a possible twin
if ($DOTWIN == 1) then
  set TWIN = `echo $TWINA | awk '{if ($1 < 0.05) {print ""} else {print "test"}}'`
endif  

#Check whether de-twinning is needed to repoduce R-factors. Only when NOT using legacy mode
if ($LEGACY == 0 && $DOTWIN == 1) then
  #Check whether the R-factors can be reproduced
  if (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 0 && $TWIN == 'test') then #R-factors do not fit, try de-twinning.
    echo " o Problem reproducing the R-factors" | tee -a $LOG
    echo -n "-Trying de-twinning " | tee -a $LOG

    #Set up detwinning
    set TWIN = 'twin'

    #Create a backup
    cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv3.log

    #Run Refmac
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_prepped.pdb \
    XYZOUT $WORKDIR/${PDBID}_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc.mtz \
    $TLSLIN \
    $LIBLIN \
    $SCATLIN \
<<eof > $WORKDIR/${PDBID}_0cyc.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref MIXE
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $TWIN
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy remarks 200 280 350
      pdbout copy expdta
      NOHARVEST
      END
eof
    if ($status) then
      echo " " | tee -a $LOG
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Switch off detwinning if Refmac finds just one twin domain
    if (`grep -a 'The number of twin domains' $WORKDIR/${PDBID}_0cyc.log | tail -n 1 | awk '{print $7}'` == 1) then
      set TWIN =  #No detwinning
      echo " " | tee -a $LOG
      echo " o Refmac detected no twinning" | tee -a $LOG
    else
      #Report R-free
      set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`
      echo "(R-free = $RFCAL)" | tee -a $LOG

      #Did Refmac find higher symmetry
      set REFHSYMM = `grep -c 'twin or higher symmetry' $WORKDIR/${PDBID}_0cyc.log`

      #Evaluate R-factors
      set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
      if (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 1) then
        #Detwinning was needed to reproduce R-factors
        #Give warnings
        if ($PHASTWIN == 0 && $REFHSYMM > 0) then
          if ($LOCAL == 1) then
            #Give big warning for server and local mode
            echo " " | tee -a $LOG
            echo "WARNING!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Data seems to have treated as twinnned before, but Phaser did not find significant twinning." | tee -a $LOG
            echo "Additionally, Refmac finds a potential space group problem." | tee -a $LOG
            echo "The data will not be treated as twinned in PDB-REDO." | tee -a $LOG
            echo "Please, check the space group carefully!" | tee -a $LOG
            echo " " | tee -a $LOG

            #Do not use twinning and copy back the old log file
            cp $WORKDIR/${PDBID}_0cycv3.log $WORKDIR/${PDBID}_0cyc.log
            set TWIN =  #No detwinning

          else
            #Give small warning for databank mode
            echo " o Data seems to have treated as twinnned before."  | tee -a $LOG
            echo "   * Phaser did not find significant twinning."     | tee -a $LOG
            echo "   * Refmac finds a potential space group problem." | tee -a $LOG
            echo "   * Reluctantly treating the data as twinned."     | tee -a $LOG

            set FALSETWIN = 1

            #Write DEBUG comment
            echo "COMMENT: false twinning problem" >> $DEBUG
            echo "PDB-REDO,$PDBID"                 >> $DEBUG

            #The data are treated as twinned. Do not calculate anomalous maps
            set ANOMCOEF =
            set ANOMCMD  =
          endif
        else if ($PHASTWIN == 0) then
          if ($LOCAL == 1) then
            #Give big warning for server and local mode
            echo " " | tee -a $LOG
            echo "WARNING!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Data seems to have treated as twinnned before and Refmac and sfcheck indicate that the data are twinned, but Phaser does not. " | tee -a $LOG
            echo "The data will be treated as twinned in PDB-REDO, reluctantly." | tee -a $LOG
            echo "Please, check your data carefully!" | tee -a $LOG
            echo " " | tee -a $LOG
          else
            #Give small warning for databank mode
            echo " o Data seems to have treated as twinnned before." | tee -a $LOG
            echo "   * Phaser did not find significant twinning."    | tee -a $LOG
            echo "   * Reluctantly treating the data as twinned."    | tee -a $LOG

            set FALSETWIN = 1
            #Write DEBUG comment
            echo "COMMENT: ambiguous twinning problem" >> $DEBUG
            echo "PDB-REDO,$PDBID"                     >> $DEBUG
          endif

          #The data are treated as twinned. Do not calculate anomalous maps
          set ANOMCOEF =
          set ANOMCMD  =
        else
          #The data are twinned. Do not calculate anomalous maps
          set ANOMCOEF =
          set ANOMCMD  =
        endif
      else
        #Detwinning did not solve the problem. Only use twinning when it is detected
        if ($PHASTWIN == 1 && $REFHSYMM == 0) then
          #The data are twinned. Do not calculate anomalous maps
          set ANOMCOEF =
          set ANOMCMD  =
        else
          #Do not use twinning and copy back the old log file
          cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv4.log
          cp $WORKDIR/${PDBID}_0cycv3.log $WORKDIR/${PDBID}_0cyc.log
          set TWIN =  #No detwinning
        endif

      endif
    endif
  endif

  #Evaluate log file (dirty)
  set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
  set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`

endif


#Check the fit with the data. If it is poor, try rigid-body refinement. This is always done when in legacy mode
if ($DORB == 1 && (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 0 || $LEGACY == 1)) then #R-factors do not fit, try rigid-body refinement.

  #Create a backup
  cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv5.log

  #Do $RBCYCLE cycles of rigid body refinement
  if ($LEGACY == 1) then
    echo "-Legacy mode" | tee -a $LOG
    echo -n " o Doing rigid-body refinement " | tee -a $LOG
    set RBCYCLE = 15
  else
    echo " o Problem reproducing the R-factors"  | tee -a $LOG
    echo -n "-Trying rigid-body refinement " | tee -a $LOG
  endif

  #Set up NCS
  if ($STRICTNCS == 1) then
    set RBNCS = ncsconstraints
  endif

  #Run Refmac
  refmac5 \
  XYZIN $WORKDIR/${PDBID}_prepped.pdb \
  XYZOUT $WORKDIR/${PDBID}_refmacrb.pdb \
  HKLIN $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/${PDBID}_refmacrb.mtz \
  $SCATLIN \
<<eof > $WORKDIR/${PDBID}_rb.log
    $SCATTERCMD
    make check NONE
    make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
      ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
    refi type RIGID resi MLKF meth CGMAT
    $REFIRES
    mode rigid
    rigid ncycle $RBCYCLE
    scal type $SOLVENT $SCALING
    scale mlscale nrfr 5
    solvent YES
    $MASKPAR
    $LOWMEM
    weight $WGTSIG MATRIX 0.5
    monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
      chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
    $RBNCS
    $TWIN
    labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
    $ANOMCMD
    pdbout copy remarks 200 280 350
    pdbout copy expdta
    NOHARVEST
    END
eof
  if ($status) then
    echo " " | tee -a $LOG
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in rigid-body refinement" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                 >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif

  #Transplant the LINKs if needed
  if (`grep -c ^LINK $WORKDIR/${PDBID}_prepped.pdb` != 0 && `grep -c ^LINK $WORKDIR/${PDBID}_refmacrb.pdb` == 0) then
    cp  $WORKDIR/${PDBID}_refmacrb.pdb $WORKDIR/${PDBID}_refmacrb.bak
    grep -B 100000 ^CRYST1 $WORKDIR/${PDBID}_refmacrb.bak | grep -v ^CRYST1 > $WORKDIR/${PDBID}_refmacrb.pdb
    grep ^LINK $WORKDIR/${PDBID}_prepped.pdb >> $WORKDIR/${PDBID}_refmacrb.pdb
    grep -A 250000 ^CRYST1 $WORKDIR/${PDBID}_refmacrb.bak >> $WORKDIR/${PDBID}_refmacrb.pdb
  endif


  #Now do another restrained refinement run to include the TLS contribution (if needed)
  if ($ORITLS == 1) then
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_refmacrb.pdb \
    XYZOUT $WORKDIR/${PDBID}_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc.mtz \
    $TLSLIN \
    $LIBLIN \
    $SCATLIN \
<<eof > $WORKDIR/${PDBID}_0cyc.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref MIXE
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $TWIN
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
      $ANOMCMD
      pdbout copy remarks 200 280 350
      pdbout copy expdta
      NOHARVEST
      END
eof
    if ($status) then
      echo " " | tee -a $LOG
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Evaluate log file
    set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
    set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`
  else
    #Evaluate log file from the rigid body refinement
    set RCAL  = `tail -n $RBLS $WORKDIR/${PDBID}_rb.log | head -n 1 | awk '{print $2}'`
    set RFCAL = `tail -n $RBLS $WORKDIR/${PDBID}_rb.log | head -n 1 | awk '{print $3}'`

    #Use the MTZ file from the rigid-body refinement
    cp $WORKDIR/${PDBID}_refmacrb.mtz $WORKDIR/${PDBID}_0cyc.mtz
  endif

  #Fill the line with the R-free value
  echo "(R-free = $RFCAL)" | tee -a $LOG

  #Copy files so that the rigid-body refined model is used
  cp $WORKDIR/${PDBID}_refmacrb.pdb $WORKDIR/${PDBID}_0cyc.pdb
  
else  
  cp $WORKDIR/${PDBID}_prepped.pdb $WORKDIR/${PDBID}_refmacrb.pdb
endif

#Check the fit with the data. If it is poor, try short TLS-refinement
if (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 0) then
  #Only do something if original TLS groups were given
  if ($ORITLS == 1) then
    echo " o Problem reproducing the R-factors" | tee -a $LOG
    echo -n "-Trying TLS tensor re-evaluation " | tee -a $LOG

    #Create a backup
    cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv5.log
    cp $WORKDIR/${PDBID}_0cyc.mtz $WORKDIR/${PDBID}_0cycv5.mtz

    #Use TLS group definitions only
    grep -E 'TLS|RANGE|^ ' $WORKDIR/${PDBID}.tls > $WORKDIR/${PDBID}_0cycin.tls

    #Run Refmac with 5 TLS cycles
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_refmacrb.pdb \
    XYZOUT $WORKDIR/${PDBID}_TLS0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc.mtz \
    $LIBLIN \
    TLSIN $WORKDIR/${PDBID}_0cycin.tls TLSOUT $WORKDIR/${PDBID}_0cyc.tls \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_0cyc.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi tlsc 5
      refi type REST resi MLKF meth CGMAT bref ISOT
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $TWIN
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES  $ANOMCOEF
      $ANOMCMD
      pdbout copy remarks 200 280 350
      pdbout copy expdta
      NOHARVEST
      kill $TOOLS/pdb_redo.refmac
      END
eof
    if ($status) then
      if (`grep -a -c 'Program terminated by user' $WORKDIR/${PDBID}_0cyc.log` != 0) then
        #Problems with the TLS group definition.
        mv $WORKDIR/${PDBID}.tls $WORKDIR/${PDBID}.notls
        #Cannot use this run. Resore the backup.
        cp $WORKDIR/${PDBID}_0cycv5.log $WORKDIR/${PDBID}_0cyc.log
      else
        echo " " | tee -a $LOG
        echo " o Problem with refmac. Cannot continue." | tee -a $LOG
        echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      endif
    else if (`grep -c -E '\*{6}' $WORKDIR/${PDBID}_TLS0cyc.pdb` != 0) then
      #TLS refinement caused rediculously high B-factors. Don't use this TLS group selection.
      mv $WORKDIR/${PDBID}.tls $WORKDIR/${PDBID}.notls
      #Cannot use this run. Resore the backup.
      cp $WORKDIR/${PDBID}_0cycv5.log $WORKDIR/${PDBID}_0cyc.log
      #Report
      echo " " | tee -a $LOG
      echo " o TLS refinement was unstable. Original TLS group selection will not be used." | tee -a $LOG
      echo "COMMENT: Problem with TLS tensor re-evaluation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                                >> $DEBUG
    else
      #Strip ANISOU records and clean up
      grep -a -v -E '^ANISOU' $WORKDIR/${PDBID}_TLS0cyc.pdb > $WORKDIR/${PDBID}_0cyc.pdb
      rm $WORKDIR/${PDBID}_TLS0cyc.pdb
      rm $WORKDIR/${PDBID}_0cycin.tls
      rm $WORKDIR/${PDBID}_0cyc.tls

      #Did the R-factor improve over the rigid-body results?
      set PRCAL1 = $RCAL
      set PRCAL2 = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`

      #Discard the TLS results if they are poor.
      if (`echo $PRCAL1 $PRCAL2 | awk '{if ($1 - $2 < 0) {print "rb"} else {print "tls"}}'` == 'rb') then
        cp $WORKDIR/${PDBID}_refmacrb.pdb $WORKDIR/${PDBID}_0cyc.pdb
        cp $WORKDIR/${PDBID}_0cycv5.mtz $WORKDIR/${PDBID}_0cyc.mtz
        cp $WORKDIR/${PDBID}_0cycv5.log $WORKDIR/${PDBID}_0cyc.log
      endif

      #Evaluate log file (dirty)
      set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
      set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`

      echo "(R-free = $RFCAL)" | tee -a $LOG
    endif
  endif
endif

#Get the geometry
set RMSZB = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | sed 's/\*\*\*\*\*\*/huge/g' | awk '{print $8}'`
set RMSZA = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | sed 's/\*\*\*\*\*\*/huge/g' | awk '{print $10}'`

#Do last ditch MR attempts if it is really bad
if ($DIDMR == 0 && `$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL 0.10` == 0 && `echo $RFCAL | awk '{if ($1 > 0.500) {print 1} else {print 0}}'` == 1) then 
  #R-factors do not fit after several tries and are very high.
  echo " o Suspiciously high R-factors" | tee -a $LOG
  
  #Run Molrep (PROGRAM MOLREP)
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.molrep.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  echo "   * Running MOLREP" | tee -a $LOG
  set DIDMR = 1
  molrep -f $WORKDIR/$PDBID.mtz -m $WORKDIR/${PDBID}_prepped.pdb > $WORKDIR/molrep.log
  
  #Start over if there is usable output
  if (-e molrep.pdb) then 
    echo "   * Restarting with MOLREP output" | tee -a $LOG
    cp molrep.pdb $WORKDIR/${PDBID}_prepped.pdb
    goto molrepped
  endif  
endif  


#Check if R-factor can (very roughly) be reproduced. If not, abort. Do not do this in legacy mode.
if ($LEGACY == 1) then
  #Write a warning if R-free is really high
  if (`echo $RFCAL | awk '{if ($1 > 0.500) {print 1} else {print 0}}'` == 1) then
    echo "COMMENT: Suspiciously high calculated R-free" >> $DEBUG
    echo "PDB-REDO,$PDBID"                              >> $DEBUG
  endif
else
  set FITR = `$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL 0.10`
  if ($FITR == 0) then #R-factors do not fit after several tries.
    echo " o Cannot reproduce Rfactor within 0.10 tolerance:" | tee -a $LOG
    echo "   * Reported Rfactor : $RFACT" | tee -a $LOG
    echo "   * Calculated R     : $RCAL"  | tee -a $LOG

    if ($INTENS != "-i" && `grep -c ' F ' $WORKDIR/${PDBID}c2c.log` > 0 &&  `grep -c ' I ' $WORKDIR/${PDBID}c2c.log` > 0) then

      #Go back and use intensities
      echo "   * Retrying with reflection intensities intead of amplitudes"  | tee -a $LOG
      echo " " | tee -a $LOG

      #Create debug entry
      echo "COMMENT: Using intensities instead of amplitudes" >> $DEBUG
      echo "PDB-REDO,$PDBID"                                  >> $DEBUG

      #Set the flag for cif2cif and go back there
      set INTENS = "-i"
      set SIGMA = ''
      goto beintens
    else if ($LOCAL == 0 && $GOTOLD == 0) then
      #See if a previous PDB-REDO entry exists and use its 0-cycle file
      $WEBGET http://pdb-redo.eu/db/$PDBID/${PDBID}_0cyc.pdb.gz
      if($status) then
        echo "   * Re-refinement aborted" | tee -a $LOG
      else
        zcat $WORKDIR/${PDBID}_0cyc.pdb.gz > $WORKDIR/${PDBID}_prepped.pdb
        echo "   * Trying previous PDB-REDO entry as starting point" | tee -a $LOG
        set GOTOLD = 1
        goto previousredo
      endif
    else
      echo "   * Re-refinement aborted" | tee -a $LOG
    endif

    #Create whynot entry
    echo "COMMENT: Cannot reproduce Rfactor within 0.10 tolerance" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                         >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif

#Copy back the SEQRES records
if ($GOTSEQRES == 1 && `grep -c '^SEQRES' $WORKDIR/${PDBID}_0cyc.pdb` == 0) then
  cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_0cyc.bak
  #Run seqrescopier. PROGRAM: seqrescopier
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.seqrescopier.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  $TOOLS/seqrescopier -v \
  -pdbinw  $WORKDIR/cache.pdb \
  -pdbinwo $WORKDIR/${PDBID}_0cyc.bak \
  -pdbout  $WORKDIR/${PDBID}_0cyc.pdb > $WORKDIR/seqrescopier.log 
else
  #Do nothing
endif 


####################################### Decide on using detwinning in refinement #########################################

#If detwinning is active, test to see if Refmac detects a twin
if ($TWIN == 'test' && $PHASTWIN == 1) then
  echo "-Evaluating twinning" | tee -a $LOG

  #Set up detwinning
  set TWIN = 'twin'

  refmac5 \
  XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
  XYZOUT $WORKDIR/${PDBID}_twin.pdb \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/${PDBID}_twin.mtz \
  $LIBLIN \
  $SCATLIN \
<<eof >$WORKDIR/${PDBID}_twin.log
    $SCATTERCMD
    make check NONE
    make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
      ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
    refi type REST resi MLKF meth CGMAT bref MIXE
    $REFIRES
    ncyc 0
    scal type $SOLVENT $SCALING
    solvent YES
    $MASKPAR
    $LOWMEM
    $NCSSTRICT
    weight $WGTSIG AUTO 2.50
    monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
      chiral  10.0 bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
    $TWIN
    labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES
    pdbout copy expdta
    pdbout copy remarks 200 280 350
    END
eof
  if ($status) then
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in twin evaluation" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                           >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif

  #Switch off detwinning if Refmac finds just one twin domain
  if (`grep -a 'The number of twin domains' $WORKDIR/${PDBID}_twin.log | tail -n 1 | awk '{print $7}'` == 1) then
    set TWIN =  #No detwinning
    echo " o Refmac detected no twinning" | tee -a $LOG
  else
    #The data are twinned. Do not calculate anomalous maps.
    set ANOMCOEF =
    set ANOMCMD  =
  endif

  #Clean up
  rm $WORKDIR/${PDBID}_twin.mtz
  rm $WORKDIR/${PDBID}_twin.pdb
else if ($TWIN == 'test') then
  #Do not try detwinning because PHASER detected no twinning
  set TWIN =  #No detwinning
endif

#Set the final twin status
if ($TWIN == "twin") then
  set ISTWIN = 1
endif

############################################## Check R-factors and geometry ##############################################

#Warn for high structure RMSZ scores
if ($RMSZB == 'huge' || $RMSZA == 'huge' || `echo $RMSZB $RMSZA | awk '{if ($1 > 10.000) {print "1"} else if ($2 > 10.000) {print "1"} else {print "0"}}'` == 1 ) then
  if ($LOCAL == 1) then
    #The warning
    echo " " | tee -a $LOG
    echo "WARNING!" | tee -a $LOG
    echo "--------" | tee -a $LOG
    echo "Extremely large geometric outliers!" | tee -a $LOG
    echo "Bond length RMSZ: $RMSZB" | tee -a $LOG
    echo "Bond angle RMSZ : $RMSZA" | tee -a $LOG
    echo "This is likely a problem with the restraint generation." | tee -a $LOG
    if ("$INREST" != "") then
      echo "Make sure that all non-standard compounds and LINKs are properly described in your restraint file."   | tee -a $LOG
    else
      echo "Consider running PDB-REDO with a restraint file that describes all non-standard compounds and LINKs." | tee -a $LOG
    endif
    echo " " | tee -a $LOG
  else
    echo "-Suspiciously high Refmac bond length or bond angle RMSZ detected" | tee -a $LOG
    echo " o Bond length RMSZ: $RMSZB"                               | tee -a $LOG
    echo " o Bond angle RMSZ : $RMSZA"                               | tee -a $LOG
    echo " o This is likely a problem with the restraint generation" | tee -a $LOG
  endif

  #Write debug message
  echo "COMMENT: suspiciously high rmsZ" >> $DEBUG
  echo "PDB-REDO,$PDBID"                 >> $DEBUG
endif

set PPATM = "4" #Assumes isotropic B-factors for now

#Calculate sigma(R-free)
set SIGRFCAL = `echo $RFCAL $NTSTCNT | awk '{printf ("%.4f\n", $1/sqrt($2))}'`

#Calculate expected R-free/R ratio (Ticke model and isotropic empirical model 2020)
set RFRRAT2 = `echo $ATMCNT $WORKCNT $PPATM $URESO | awk '{X = $1/$2; Z = 1 - $3 * X/(1 + 2.5 * X); if (Z < 0) {Z = 0}; A = $3-2.5*(1-Z*Z*Z*Z*Z); X = A*X; RATIO = 1-X; if(RATIO <= 0) {RATIO = 1.010} else {RATIO = sqrt((1+X)/RATIO)};  if ($4 > 2.65 && RATIO > 1.200000) {RATIO = 1.200000}; if ($4 > 3.0 && RATIO < 1.011) {RATIO = 1.20}; if (RATIO > 1.454) {RATIO = 1.45}; printf ("%.4f\n", RATIO)}'`
set RFRRAT  = `echo $ATMCNT $WORKCNT | awk '{printf ("%.4f\n", -0.005 * $2/$1 + 1.2286)}'`
set SRFRRAT = `echo $ATMCNT $WORKCNT | awk '{printf ("%.4f\n", -0.024 * log($2/$1) + 0.1064)}'`
set ZRFRRATCAL = `echo $ATMCNT $WORKCNT $RCAL $RFCAL | awk '{printf ("%.2f\n", (((-0.005 * $2/$1 + 1.2286) - $4/$3)/(-0.024 * log($2/$1) + 0.1064)))}'`

#Calculate expected R-free
set RFCALUNB = `echo $RCAL $RFRRAT  | awk '{printf ("%.4f\n", $1*$2)}'`

#R-free Z-score (comparison of R-free with its expected value)
set RFCALZ = `echo $RFCALUNB $RFCAL $SIGRFCAL | awk '{printf ("%.2f\n", ($1-$2)/$3)}'`

#Rcalculate the R-ratio Z-score"

#Print values
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** R-factor and R-free details ******" | tee -a $LOG
echo "Calculated R      : $RCAL"       | tee -a $LOG
echo "Calculated R-free : $RFCAL"      | tee -a $LOG
echo "Expected R-free/R : $RFRRAT (new)"   | tee -a $LOG
echo "Expected R-free/R : $RFRRAT2 (old)"  | tee -a $LOG
echo "Expected R-free   : $RFCALUNB"   | tee -a $LOG
echo " " | tee -a $LOG
echo "sigma(R-free/R)   : $SRFRRAT"    | tee -a $LOG
echo "R-free/R Z-score  : $ZRFRRATCAL" | tee -a $LOG
echo " " | tee -a $LOG
echo "sigma(R-free)     : $SIGRFCAL"   | tee -a $LOG
echo "R-free Z-score    : $RFCALZ"     | tee -a $LOG 


############################################ Fix atom chiralities if needed ##############################################

#Count the chirality problems
set CHIRERR = `grep -a -A 100 "Chiral volume deviations from the" $WORKDIR/${PDBID}_0cyc.log | grep -E "^.{19}mod" | wc -l`

#Are fixes needed?
if ($CHIRERR != 0) then

  #Report
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Chirality validation ******" | tee -a $LOG
  echo "-Found $CHIRERR chirality problems" | tee -a $LOG
  echo " o Running chiron"   | tee -a $LOG

  #Do the fixes and report
  mv $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_0cyc.old

  #PROGRAM: chiron
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.chiron.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  $TOOLS/chiron -v $TOOLS/pdb_redo.dat $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cyc.old $WORKDIR/${PDBID}_0cyc.pdb > $WORKDIR/chiron.log
  set CHIFIX = `grep "ters fixed" $WORKDIR/chiron.log | awk '{print $5}'`
  echo "   * `grep 'tested' $WORKDIR/chiron.log`"         | tee -a $LOG
  echo "   * `grep 'ters fixed' $WORKDIR/chiron.log`"     | tee -a $LOG
  echo "   * `grep 'not fixed' $WORKDIR/chiron.log`"      | tee -a $LOG
  echo "   * `grep 'Unknown chiral' $WORKDIR/chiron.log`" | tee -a $LOG

  #Make a debug record for unfixable problems
  if (`grep -c Unknown $WORKDIR/chiron.log` != 0 && `grep Unknown $WORKDIR/chiron.log | awk '{print $5}'` != 0) then
    if ($LOCAL != 1) then
      echo "$PDBID :"                                           >> $CHIRALS
      grep 'Unknown residue' $WORKDIR/chiron.log | cut -c 19-61 >> $CHIRALS
    else
      echo "COMMENT: Unknown chirality errors" >> $DEBUG
      echo "PDB-REDO,$PDBID"                   >> $DEBUG
    endif
  endif

  #Make a debug record if chiron failed
  if (-e $WORKDIR/${PDBID}_0cyc.pdb) then
    #Do nothing
  else
    mv $WORKDIR/${PDBID}_0cyc.old $WORKDIR/${PDBID}_0cyc.pdb
    echo "COMMENT: CHIRON: general error" >> $DEBUG
    echo "PDB-REDO,$PDBID"                >> $DEBUG
  endif
else
  set CHIFIX = 0
endif

#############################################  Fix the backbone if needed  ###############################################


if ($GOT_PROT == 'T' && $DOFIXDMC == 1) then
  #Report
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Backbone correction ******" | tee -a $LOG
  echo "-Running fixDMC" | tee -a $LOG
  
  #Add OXT only when the sequence is trusted
  if ($TRUSTSEQ == 1) then
    set OXTADD = '-noaddoxt'
    echo " o The sequence is untrusted, no OXT atoms can be added" | tee -a $LOG
  endif

  #Run fixDMC
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.fixDMC.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  $TOOLS/fixDMC -v \
  -pdb $WORKDIR/${PDBID}_0cyc.pdb \
  -output-name $PDBID \
  -tools $TOOLS \
  $OXTADD \
  -fasta $WORKDIR/$PDBID.fasta >& $WORKDIR/fixDMC.log 

  #Report
  set NATMADD = `grep -A 2 'Total number backbone atoms added:' $WORKDIR/fixDMC.log |  awk 'BEGIN {SUM = 0}{SUM = SUM+$6} END {print SUM}'` 
  echo " o FixDMC added $NATMADD missing backbone atoms" | tee -a $LOG

  if (-e $WORKDIR/${PDBID}_fixDMC.pdb) then 
    mv $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_0cyc.old
    mv $WORKDIR/${PDBID}_fixDMC.pdb $WORKDIR/${PDBID}_0cyc.pdb
    mv $WORKDIR/${PDBID}.fasta.new $WORKDIR/${PDBID}.fasta >& /dev/null
  else
    echo " o FixDMC failed" | tee -a $LOG
    echo "COMMENT: CHIRON: general error" >> $DEBUG
    echo "PDB-REDO,$PDBID"                >> $DEBUG
  endif
endif  
  

###############################################  Set up occupancy fixing  ################################################
#Set counter
set OCCREF = 0

#Set up occupancy refinement if it is not surpressed
if ($DOOCC == 1) then
  #Create file with hetero compounds with more than two occupancies if needed
  if (`grep -E '^HETATM'+'.{10} ' $WORKDIR/${PDBID}_0cyc.pdb |  grep -v MSE | cut -c 22-27,56-60 | sort -u | cut -c 1-6 | uniq -c | awk '{if ($1 > 2) {$1 = ""; print substr($0,2,6)}}' | wc -l` != 0) then
    #Grep the hetero compounds, but leave out MSE (seleno-methionine), then cut out the chains, residue numbers and
    #occupancies and sort. This will leave only the unique occupancies per residue. The chains and residue numbers are
    #cut out and for each residue the count is given. If this count is > 2, the chain and residue number is printed to
    # a file.
    grep -E '^HETATM'+'.{10} ' $WORKDIR/${PDBID}_0cyc.pdb | grep -v MSE | cut -c 22-27,56-60 | sort -u | cut -c 1-6 | uniq -c | awk '{if ($1 > 2) {print substr($0,length($0)-5, length($0))}}' > $WORKDIR/occupancy.lst

  endif

  #Find hetero compounds with at least one atom with occupancy 0.01
  foreach RES ("`grep -E '^HETATM'+'.{10} ' $WORKDIR/${PDBID}_0cyc.pdb | grep -v MSE | cut -c 22-27,56-60 | grep ' 0.01' | cut -c 1-6 | sort -u`")
    echo "$RES" >> $WORKDIR/occupancy.lst
  end

  #Filter the residue list and write the Refmac command file
  if (-e $WORKDIR/occupancy.lst) then
    foreach RES ("`sort -u $WORKDIR/occupancy.lst`")
      @ OCCREF = ($OCCREF + 1)
      set CHNID = `echo $RES | cut -c 1-1`
      set RESN  = `echo $RES | cut -c 2-6`
      #Skip cases where the residue has an insertion code (i.e. $RESN has non-numeric characters)
      if (`echo $RESN | awk '{if ($0 ~ /^[0-9]+$/) {print "1"} else {print "0"}}'` == 1 ) then
        echo "occupancy group id $OCCREF chain $CHNID residue $RESN" >> $WORKDIR/occupancy_cmd.refmac
      else
        @ OCCREF = ($OCCREF - 1)
      endif
    end
  endif

  #Is occupancy refinement needed?
  if (-e $WORKDIR/occupancy_cmd.refmac) then

    #Append refinement command
    echo "occupancy refine" >> $WORKDIR/occupancy_cmd.refmac

    #Make Refmac use command file
    set OCCCMD = "@$WORKDIR/occupancy_cmd.refmac"
  endif
endif


############################################### Resolution-based settings ################################################
#Label to come back to after updating the resolution
resobasedsettings:

#Find the resolution type (6 categories, do not use reflections per atom for heavy strict NCS)
# Resolution categories: 0 = extremely low, 1 = very low, 2 = low, 3 = medium, 4 = high, 5 = atomic)
if ($STRICTNCS == 1 && $NMTRIX > 9) then
  set RESOTYPE = `echo $URESO | awk '{if ($1 > 4.99) {print "0"} else if ($1 > 3.49 && $1 < 5.00) {print "1"} else if ($1 > 2.79 && $1 < 3.50) {print "2"} else if ($1 < 1.21) {print "5"} else if ($1 > 1.20 && $1 < 1.71) {print "4"} else {print "3"}}'`
else
  set RESOTYPE = `echo $URESO $WORKCNT $ATMCNT | awk '{if ($1 > 4.99) {print "0"} else if ($2/$3 < 1.0) {print "0"} else if ($1 > 3.49 && $1 < 5.00) {print "1"} else if ($2/$3 < 2.5) {print "1"} else if ($1 > 2.79 && $1 < 3.50) {print "2"} else if ($1 < 1.21) {print "5"} else if ($1 > 1.20 && $1 < 1.71) {print "4"} else {print "3"}}'`
endif

#Modify the resolution category based on user input and experiment type
if ($ISED == 1) then
  @ RESOTYPE = ($RESOTYPE - 1)
endif
@ RESOTYPE = ($RESOTYPE + $RESCATSTEP)

#Set refinement parameters
if ($RESOTYPE < 1) then
  set WEIGHTS    = "1e-7 1e-6 1e-5 1e-4 5e-4 .001"
  set BWEIGHTS   = "2.50 2.00 1.50 1.20 1.00"  #Do not try looser B-factor restraints
  set JELLY      = "ridg dist sigm 0.02"       #Use jelly-body refinement
  @ NCYCLE = ($NCYCLE + 15)
  set TORSION    = "restr tors include group peptide"
  set ESTRICT    = "-e"                        #Picker is extra strict
  #Use homology restraints
  if (-e $WORKDIR/$PDBID.fasta) then
    set DOHOMOLOGY = 1
  endif
else if ($RESOTYPE == 1) then
  set WEIGHTS  = "1e-4 5e-4 .001 .002 .005 0.01"
  set BWEIGHTS = "2.00 1.50 1.20 1.00 0.80 0.50"
  set JELLY    = "ridg dist sigm 0.05" #Use jelly-body refinement
  @ NCYCLE = ($NCYCLE + 10)
  set ESTRICT  = "-e"                  #Picker is extra strict
  #Use homology restraints
  if (-e $WORKDIR/$PDBID.fasta) then
    set DOHOMOLOGY = 1
  endif
else if ($RESOTYPE == 2) then
  set WEIGHTS  = ".001 .002 .005 0.01 0.03 0.05"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30"
  set JELLY    = "ridg dist sigm 0.10" #Use jelly-body refinement (with looser restraints)
  @ NCYCLE = ($NCYCLE + 5)
  #Use homology restraints
  if (-e $WORKDIR/$PDBID.fasta) then
    set DOHOMOLOGY = 1
  endif
else if ($RESOTYPE == 3) then
  set WEIGHTS  = ".005 0.01 0.03 0.05 0.10 0.30 0.50"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30 0.10"
else if ($RESOTYPE == 4) then
  set WEIGHTS  = "0.05 0.10 0.30 0.50 0.70 1.00 1.50"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30 0.10"
else if ($RESOTYPE > 4) then
  set WEIGHTS  = "0.50 0.70 1.00 1.50 2.00 3.00 5.00"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30 0.10"
endif

#Set NCS restraints
if ($DONCS == 1) then
  #Check whether NCS is really detected
  if (`grep -c 'Ncsr group' $WORKDIR/${PDBID}_0cyc.log` == 0) then
    #Switch off NCS
    set DONCS      = 0
    set NCSTYPE    =        #Empty NCS command
    set NCSALIGN   =        #Only usefull when NCS is used (alignment cut-off)
    set NCSNEIGH   =        #Only usefull when NCS is used (include neighbouring atoms in restraints)
  endif
endif

#Use harmonic retraints for very small data sets
if ($DOHARMONIC == 1 && $NREFCNT < 1000) then
  set HARMCMD = "ridge atoms 0.05"
  @ NCYCLE = ($NCYCLE + 10)
endif

#Switch off jelly-body refinement
if ($DOJELLY == 0) then
  set JELLY = ""
endif

#Switch off homology based restraints
if ($NOHOMOLOGY == 1) then
  set DOHOMOLOGY = 0
endif

#Use the Wilson B-factor unless it is negative or the resolution is lower than 4.00A; maximise at 50A^2.
set BSET  = `echo $URESO $BWILS $BAVER | awk '{if ($1 > 3.99 || $2 < 0) {BSET = 0.5*$3} else {BSET = 0.5*$2}; if (BSET > 50) {BSET = 50.00}; if (BSET < 10) {BSET = 10.00}; printf ("%.2f\n", BSET)}'`
set TBCMD = "bfac set $BSET"


############################################# Generate external restraints ###############################################
if ($GOT_NUC == 'T' || $HBONDREST == 1 || $DOHOMOLOGY == 1) then
  echo "" | tee -a $LOG
  echo "" | tee -a $LOG
  echo "****** Structure specific restraints and targets ******" | tee -a $LOG

  #Genrate nucleic acid restraints. Not for high or atomic resolution categories. PROGRAMS: libg and bphbond
  if ($GOT_NUC == 'T' && $DONUCR == 1 && $RESOTYPE < 4) then
    echo "-Generating nucleic acid restraints" | tee -a $LOG
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.libg.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    libg -p $WORKDIR/${PDBID}_0cyc.pdb -o $WORKDIR/tstack.rest -w dsp >  $WORKDIR/nucrest.log
    touch $WORKDIR/tstack.rest
    #Clean the stacking restraints
    awk '{if (NF < 68) {print $0}}' $WORKDIR/tstack.rest | grep -v 'plan 3' > $WORKDIR/stacking.rest
    
    #Create base-pair H-bond restraints 
    python3 $TOOLS/bphbond.py $WORKDIR/${PDBID}_0cyc.pdb $TOOLS/x3dna-dssr > $WORKDIR/bphbond.tmp
    
    #Clean the restraints
    #Remove gap LINKRs
    grep -v -E 'LINKR.{67}gap' $WORKDIR/${PDBID}_0cyc.pdb > $WORKDIR/${PDBID}_0cyc.tmp
    
    #Make list of all atoms with alternates
    echo "SELECT auth_asym_id,auth_seq_id,auth_atom_id FROM atom_site WHERE label_alt_id IS NOT NULL;" | $TOOLS/mmCQL $WORKDIR/${PDBID}_0cyc.tmp >& $WORKDIR/allalt.list
    
    #Loop over all restraints
    foreach REST (`cat $WORKDIR/bphbond.tmp | tr ' ' '_'`) 
      set CHID1 = `echo $REST | cut -d '_' -f 5`
      set RESN1 = `echo $REST | cut -d '_' -f 7`
      set ATOM1 = `echo $REST | cut -d '_' -f 11`
      set CHID2 = `echo $REST | cut -d '_' -f 14`
      set RESN2 = `echo $REST | cut -d '_' -f 16`
      set ATOM2 = `echo $REST | cut -d '_' -f 20`
      
      #Only write out restraint if neither party is in the list of atoms with alternates
      if (`grep $CHID1 $WORKDIR/allalt.list | grep -w -e "$RESN1" | grep -c "$ATOM1"$` == 0 && `grep $CHID2 $WORKDIR/allalt.list | grep -w -e "$RESN2" | grep -c "$ATOM2"$` == 0) then
        echo $REST | tr '_' ' ' >> $WORKDIR/bphbond.rest
      else
        #echo $REST
      endif
    end

    #Make Refmac use command files
    if (-e $WORKDIR/stacking.rest || -e $WORKDIR/bphbond.rest) then
      echo "external weight scale 5" > $WORKDIR/nucleic.rest
      if (-e $WORKDIR/stacking.rest) then
        sed 's/  / /g' $WORKDIR/stacking.rest >> $WORKDIR/nucleic.rest
        rm $WORKDIR/stacking.rest
      endif 
      echo "external weight scale 2" >> $WORKDIR/nucleic.rest
      if (-e $WORKDIR/bphbond.rest) then
        cat $WORKDIR/bphbond.rest >> $WORKDIR/nucleic.rest
      endif
      echo "external weight scale 1" >> $WORKDIR/nucleic.rest
      set NUCRCMD = "@$WORKDIR/nucleic.rest"
    endif
  else if ($GOT_NUC == 'T') then  
    echo "-Generating basepair hydrogen bond targets" | tee -a $LOG
    python3 $TOOLS/bphbond.py $WORKDIR/${PDBID}_0cyc.pdb $TOOLS/x3dna-dssr > $WORKDIR/bphbond.rest
  endif
 
  #Generate a DSSP file
  if ($HBONDREST == 1 || $DOHOMOLOGY == 1) then
    #Define the secondary structure
    mkdssp -i $WORKDIR/${PDBID}_0cyc.pdb -o $WORKDIR/${PDBID}_0cyc.dssp >& $WORKDIR/dssp.log

    if ( ! -e $WORKDIR/${PDBID}_0cyc.dssp) then
      #DSSP failed. Cannot make restraints
      echo " o Cannot produce hydrogen bond or homology-based restrains" | tee -a $LOG
      echo "DSSP: general error in restraint generation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                             >> $DEBUG

      #Switch off hydrogen bond and homology restraints
      set HBONDREST  = 0
      set DOHOMOLOGY = 0
  endif

  #Generate hydrogen bond restraints
  if ($HBONDREST == 1) then

    echo "-Generating hydrogen bond restraints" | tee -a $LOG

    #Generate hydrogen bond restraints if a DSSP file can be made.

      #Make the hydrogen bond restraints. PROGRAM: detectHbonds
      cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.detectHbonds.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
      
      $TOOLS/detectHbonds -v \
      -pdb $WORKDIR/${PDBID}_0cyc.pdb \
      -dssp $WORKDIR/${PDBID}_0cyc.dssp \
      -output-name $PDBID \
      -tools $TOOLS > $WORKDIR/hbondrest.log
      if ( ! -e $WORKDIR/${PDBID}_hbonds.rest) then
        echo " o Cannot generate hydrogen bond restraints" | tee -a $LOG
        echo "detectHbonds: general error" >> $DEBUG
        echo "PDB-REDO,$PDBID"             >> $DEBUG
      else
        #Set up restraint commands
        sed 's/ [ ]*/ /g' $WORKDIR/${PDBID}_hbonds.rest > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
        rm $WORKDIR/${PDBID}_hbonds.rest
        set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
        set HBONDCMD    = "@$WORKDIR/hbond.rest"
      endif
    endif
  endif

  #Generate homology-based restraints
  if ($DOHOMOLOGY == 1 && $GOT_PROT == 'T' && -e $WORKDIR/$PDBID.blast) then
    echo "-Generating homology-based restraints" | tee -a $LOG

    #Make a temporary directory
    mkdir -p $WORKDIR/homol
    cd $WORKDIR/homol
    
    #Copy in extra homologous structure models
    if ("$XHOM" != "") then
      echo "-Importing homologous structure models" | tee -a $LOG
      set CNT = 0
      #Loop over files
      foreach FIL ($XHOM)
        #Only use files with exactly 1 CRYST1 card.
        if (`grep -c CRYST1 $FIL` != 1) then
          echo " o Homologous structure model $FIL does not have exactly one CRYST1 card. It will be ignored" | tee -a $LOG
        else
          set CNT  = `expr $CNT + 1`
          set RANK = `seq -w $CNT 99 | head -n 1`
          cp $FIL $WORKDIR/homol/hom$CNT.pdb
          
          #Construct command line
          set HOMINCMD = "$HOMINCMD -homol hom$CNT.pdb"
        endif
      end
      echo " o Homologous models imported: $CNT" | tee -a $LOG
    endif

    #Generate the restraints
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.HODER.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    $TOOLS/hoder -v \
    -alignout \
    -hitsummary \
    $HOMINCMD \
    -pdb $WORKDIR/${PDBID}_0cyc.pdb \
    -dssp $WORKDIR/${PDBID}_0cyc.dssp \
    -blast $WORKDIR/$PDBID.blast \
    -fasta $WORKDIR/$PDBID.fasta \
    -output-name $PDBID \
    -output-dir $WORKDIR/homol \
    -pdbdir $REDODIR \
    -edsdir $EDSDIR \
    -tools $TOOLS >> homologs.log

    #Set up the restraint commands
    #Homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/homol/${PDBID}_homologyBased.restr > $WORKDIR/homology.rest #Replace one or more spaces by a single one
    set HOMOLWGTCMD = "EXTERNAL WEIGHT SCALE $HOMOLRESTRWGT"
    set HOMOLCMD    = "@$WORKDIR/homology.rest"

    #Non-homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/homol/${PDBID}_general.restr > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
    set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
    set HBONDCMD    = "@$WORKDIR/hbond.rest"

    #Go back down
    cd $WORKDIR
  endif
  
endif


############################################### Solvent mask optimisation ################################################


#Swith to simple model
set SOLVENT = SIMP

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Solvent mask optimisation ******" | tee -a $LOG

#Label to go back to after problems with external restraints
solventtest:

#Run Refmac
echo "-Running Refmac grid search" | tee -a $LOG

refmac5 \
XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
XYZOUT $WORKDIR/${PDBID}_solvent.pdb \
HKLIN  $WORKDIR/$PDBID.mtz \
$LIBLIN \
$SCATLIN \
<<eof >$WORKDIR/${PDBID}_solvent.log
  $SCATTERCMD
  make check NONE
  make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
    ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
  make newligand continue
  refi type REST resi MLKF meth CGMAT bref MIXE
  $REFIRES
  ncyc 0
  scal type $SOLVENT $SCALING
  solvent YES
  solvent optimise
  $NCSSTRICT
  $LOWMEM
  weight $WGTSIG AUTO 2.50
  monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
   chiral 10.0   bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
  $TWIN
  $METALCMD
  $NUCRCMD
  $RESTCMD
  $HOMOLWGTCMD
  $HOMOLCMD
  $HBONDWGTCMD
  $HBONDCMD
  END
eof
if ($status) then
  echo " o Problem with refmac." | tee -a $LOG
  echo "COMMENT: refmac: error solvent optimisation" >> $DEBUG
  echo "PDB-REDO,$PDBID"                             >> $DEBUG

  #Set up fallback values
  set VDWPROBE = 'NA'
  set IONPROBE = 'NA'
  set RSHRINK  = 'NA'

else if (`grep -c 'Error: At least one of the atoms from the restraints' $WORKDIR/${PDBID}_solvent.log` > 0) then
  #There are libg or zen restraint problems report only the first time
  if ($NPRUNE == 0) then
    echo " o Problematic external restraints" | tee -a $LOG
    echo "COMMENT: refmac: external restraint problems" >> $DEBUG
    echo "PDB-REDO,$PDBID"                              >> $DEBUG
  endif

  #try pruning the restraints and try again.
  if ($NPRUNEN < $MAXPRUNE && $NPRUNEM < $MAXPRUNE) then
    cp $WORKDIR/${PDBID}_solvent.log $WORKDIR/${PDBID}_solventv$NPRUNE.log
    if (`grep -c 'Error: At least one of the atoms from the restraints' $WORKDIR/${PDBID}_solvent.log` != 0) then
      set BADREST = `grep -A 1 'Error: At least one of the atoms from the restraints' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | cut -c 6-`
    else if (`grep -c 'Plane number can be 1 or 2' $WORKDIR/${PDBID}_solvent.log` != 0) then
      set BADREST = `grep -A 1 'Plane number can be 1 or 2' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | cut -c 6-`
    else
      echo "   * Cannot locate problematic restraint. Stopping." | tee -a $LOG
      echo "COMMENT: refmac: unsolvable external restraint problems" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                         >> $WHYNOT
      exit(1)
    endif
    echo "   * Removing problematic restraint" | tee -a $LOG
    if (-e $WORKDIR/nucleic.rest) then
      cp $WORKDIR/nucleic.rest $WORKDIR/temp.rest
      grep -v "$BADREST" temp.rest > $WORKDIR/nucleic.rest
      @ NPRUNEN = ($NPRUNEN + 1)
    endif
    if (-e $WORKDIR/metal.rest) then
      cp $WORKDIR/metal.rest $WORKDIR/temp.rest
      grep -v "$BADREST" temp.rest > $WORKDIR/metal.rest
      @ NPRUNEM = ($NPRUNEM + 1)
    endif
    @ NPRUNE = ($NPRUNEM + $NPRUNEN)
    goto solventtest
  else
    #Give up for the specific restraint type
    if ($NPRUNEM == $MAXPRUNE)  then
      #Stop using metal restraints (reset the counter)
      set METALCMD =
      @ NPRUNEM = ($MAXPRUNE - 1)
    endif
    if ($NPRUNEN == $MAXPRUNE)  then
      #Stop using nucleic acid restraints (reset the counter)
      set NUCRCMD =
      @ NPRUNEN = ($MAXPRUNE - 1)
    endif
    cp $WORKDIR/${PDBID}_solvent.log $WORKDIR/${PDBID}_solventv$NPRUNE.log
    goto solventtest
  endif
else
  #Set up solvent parameters.
  set VDWPROBE = `grep -a -A 3 -e 'Minimum R f*ree' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | awk '{print $4}' | cut -c 1-3`
  set IONPROBE = `grep -a -A 4 -e 'Minimum R f*ree' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | awk '{print $4}' | cut -c 1-3`
  set RSHRINK  = `grep -a -A 5 -e 'Minimum R f*ree' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | awk '{print $3}' | cut -c 1-3`
  set MASKPAR  = "solvent vdwprobe $VDWPROBE ionprobe $IONPROBE rshrink $RSHRINK"

  #Report
  echo " o VDW probe: $VDWPROBE" | tee -a $LOG
  echo " o Ion probe: $IONPROBE" | tee -a $LOG
  echo " o Shrinkage: $RSHRINK"  | tee -a $LOG

  #Clean
  rm $WORKDIR/${PDBID}_solvent.pdb
endif


########################################## Extend the resolution #########################################################

#Extend only if the data has higher resolution than the PDB header
set RESOGAP = `echo $RESOLUTION $DATARESH | awk '{if ($1 - $2 > 0.10) {print "1"} else {print "0"}}'`
if ( ($RESOGAP == 1 || $FORCEPAIRED == 1) && $RESOCHECK == 1) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Testing resolution cut-off ******" | tee -a $LOG

  #Don't come back here
  set RESOCHECK = 0

  #Run binliner
  echo "-Running binliner" | tee -a $LOG
  #PROGRAM: binliner
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.binliner.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  if ($RESOGAP == 1) then
    $TOOLS/binliner -v \
    $WORKDIR/$PDBID.cif \
    $RESOLUTION $DATARESH > $WORKDIR/binliner.log
  else
    $TOOLS/binliner -v \
    $WORKDIR/$PDBID.cif > $WORKDIR/binliner.log
  endif  
  
  
  echo "-Testing resolution cut-offs:" `tail -n 1 $WORKDIR/binliner.log` | tee -a $LOG

  #Initialise values
  set RLOWER   = 1.00
  set RFLOWER  = 1.00
  set WRFLOWER = 1.00
  set FLLOWER  = 999999.9
  set FCCLOWER = 0.00
  set RESLOWER = $RESOLUTION
  set RESSTEP  = -1
  @ NRCYCLE = ($NCYCLE + 10)

  #Start a pretty logfile
  echo "CUTOFF RLOWER RHIGHER RFLOWER RFHIGHER WRFLOWER WRFHIGHER FLLOWER FLHIGHER FCCLOWER FCCHIGHER" > $WORKDIR/${PDBID}_resotest.log

  #Refine to extend the resolution
  foreach RESCO (`tail -n 1 $WORKDIR/binliner.log`)

    echo " o Testing resolution $RESCO" | tee -a $LOG
    #Set resolution cut-offs
    set REFIRES = "RESO $RESCO"

    #Refine against data
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
    XYZOUT $WORKDIR/${PDBID}_res$RESCO.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_restest.mtz \
    $LIBLIN \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_all$RESCO.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout YES peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      make newligand continue     
      $TBCMD
      refi type REST resi MLKF meth CGMAT bref ISOT
      $REFIRES
      ncyc $NRCYCLE
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG AUTO 2.50
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $HARMCMD
      $METALCMD
      $NUCRCMD
      $RESTCMD
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy expdta
      pdbout copy remarks 200 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      END
eof
    if ($status) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in resolution extension" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Calculate R(-free) against lower resolution
    set TESTRES = "RESO $RESLOWER"

    #Run Refmac 0-cycles
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_res$RESCO.pdb \
    XYZOUT $WORKDIR/${PDBID}_restest_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_restest_0cyc.mtz \
    $LIBLIN \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_high2low$RESCO.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout YES peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      make newligand continue
      refi type REST resi MLKF meth CGMAT bref ISOT
      $TESTRES
      ncyc 0
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG AUTO 2.50
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $HARMCMD
      $METALCMD
      $NUCRCMD
      $RESTCMD
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy expdta
      pdbout copy remarks 200 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      END
eof
    if ($status) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in resolution extension" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Accept or reject the higher resolution data
    set RHIGHER   = `tail -n $LOGSTEP $WORKDIR/${PDBID}_high2low$RESCO.log | head -n 1 | awk '{print $2}'`
    set RFHIGHER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_high2low$RESCO.log | head -n 1 | awk '{print $3}'`
    set WRFHIGHER = `grep "Free weighted R2 factor" $WORKDIR/${PDBID}_high2low$RESCO.log | awk '{print $6}'`
    set FLHIGHER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_high2low$RESCO.log | head -n 1 | awk '{print $6}'`
    set FCCHIGHER = `grep "Free correlation coefficient" $WORKDIR/${PDBID}_high2low$RESCO.log | awk '{print $5}'`

    #Continue logfile
    echo $RESCO $RLOWER $RHIGHER $RFLOWER $RFHIGHER $WRFLOWER $WRFHIGHER $FLLOWER $FLHIGHER $FCCLOWER $FCCHIGHER >> $WORKDIR/${PDBID}_resotest.log

    #Is the best resolution so far the same as the higher resolution? PROGRAM: resolute
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.resolute.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    $TOOLS/resolute -v $WORKDIR/${PDBID}_resotest.log 1 > $WORKDIR/resolute.log
    #Report deteriorated metrics, if any
    if (`grep -A 12 "Resolution $RESCO" $WORKDIR/resolute.log | grep -c deteriorated` > 0) then
      grep -A 12 "Resolution $RESCO" $WORKDIR/resolute.log | grep deteriorated | xargs -d '\n' -n 1 echo "   *" | tee -a $LOG
    endif
    
    
    #Accept new cut-off?
    if (`tail -n 1 $WORKDIR/resolute.log` == $RESCO) then
      #Yes. Current resolution is better than previous one. Continue.
      @ RESSTEP = ($RESSTEP + 1)
    else
      #No. Reject the resolution and exit the loop.
      goto gotbestreso
    endif

    #Prepare for next cycle get R(-free) from 0 cycles
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_res$RESCO.pdb \
    XYZOUT $WORKDIR/${PDBID}_restest_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_restest_0cyc.mtz \
    $LIBLIN \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_0cyc$RESCO.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout YES peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      make newligand continue  
      refi type REST resi MLKF meth CGMAT bref ISOT
      $REFIRES
      ncyc 0
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG AUTO 2.50
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $HARMCMD
      $METALCMD
      $NUCRCMD
      $RESTCMD
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy expdta
      pdbout copy remarks 200 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      END
eof
    if ($status) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in resolution extension" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Set values for next round
    set RLOWER   = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$RESCO.log | head -n 1 | awk '{print $2}'`
    set RFLOWER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$RESCO.log | head -n 1 | awk '{print $3}'`
    set WRFLOWER = `grep "Free weighted R2 factor" $WORKDIR/${PDBID}_0cyc$RESCO.log | awk '{print $6}'`
    set FLLOWER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$RESCO.log | head -n 1 | awk '{print $6}'`
    set FCCLOWER = `grep "Free correlation coefficient" $WORKDIR/${PDBID}_0cyc$RESCO.log | awk '{print $5}'`
    set URESO    = $RESCO
    set RESLOWER = $RESCO
  end

  #Update the reflection count
gotbestreso:

  #Run MTZUTILS to cut the data
  mtzutils \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/lowres.mtz \
  <<eof >>$WORKDIR/mtz_creation.log
    RESOLUTION $URESO $DATARESL
    END
eof
  if ($status) then
    echo "   * Error using MTZUTILS. Cannot continue." | tee -a $LOG
    echo "COMMENT: mtzutils: general error (2)" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                      >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif

  #Get the correct number of (R-free) reflections
  mtz2various \
  HKLIN $WORKDIR/lowres.mtz \
  HKLOUT $WORKDIR/temp.hkl \
<<eof >>$WORKDIR/mtz_creation.log
    OUTPUT CIF data_temp
    LABIN FP=FP SIGFP=SIGFP FREE=FREE
    END
eof
  set NTSTCNT = `grep -cE ' f ' $WORKDIR/temp.hkl`
  set NREFCNT = `grep -cE ' [of] ' $WORKDIR/temp.hkl`
  set WORKCNT = `grep -cE ' o ' $WORKDIR/temp.hkl`

  #Cleanup
  rm $WORKDIR/temp.hkl
  rm $WORKDIR/lowres.mtz
  rm $WORKDIR/${PDBID}_restest.mtz
  rm $WORKDIR/${PDBID}_restest_0cyc.pdb
  rm $WORKDIR/${PDBID}_restest_0cyc.mtz

  #Update the resolution based settings
  echo "-High resolution cut-off: $URESO" | tee -a $LOG
  set REFIRES = "RESO $URESO"

  #Update the number of refinement cycles only if there was a resolution gap:
  if ($RESSTEP > 0 && $RESOGAP == 1) then
    set TCYCLE = $NCYCLE
    @ NCYCLE = ($NCYCLE + $RESSTEP * 5)
    echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
  endif

  #Update R-factor statistics
  set RCAL  = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_all$URESO.log | tail -n 1 | awk '{print $2}'`
  set RFCAL = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_all$URESO.log | tail -n 1 | awk '{print $3}'`

  goto resobasedsettings
endif


####################################### Decide on type of B-factor refinement ############################################

#Try refinement with anisotropic B's and no TLS if the data parameter ratio allows it (not with heavy strict NCS)
if ($STRICTNCS == 1 && $NMTRIX > 9) then
  set BREFTYPE     = "ISOT"
else if ($BISOT == 1) then  
  #Force isotropic B-factors, within reason
  set BREFTYPE = `echo $WORKCNT $ATMCNT | awk '{if ($1/$2 > 30.0) {print "ANISOT"} else {print "ISOT"}}'`
else if ($FISOT == 1) then  
  #Force isotropic B-factors, outside reason
  set BREFTYPE = "ISOT"
else if ($BANISOT == 1) then
  #Force anisotropic B-factors, within reason
  set BREFTYPE = `echo $URESO $WORKCNT $ATMCNT | awk '{if ($1 > 1.94) {print "ISOT"} else if ($2/$3 > 13.0) {print "ANISOT"} else {print "ISOT"}}'`
else
   #Use the Hamilton test if needed
  set BREFTYPE = `echo $URESO $WORKCNT $ATMCNT | awk '{if ($1 > 1.94) {print "ISOT"} else if ($2/$3 > 30.0) {print "ANISOT"} else if ($2/$3 > 13.0) {print "TEST"} else {print "ISOT"}}'`
endif


set REFPATM  = `echo $WORKCNT $ATMCNT | awk '{printf ("%.1f\n", $1/$2)}'`

#Test the best B-factor restraint type
if ($BREFTYPE == TEST) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Testing B-factor refinement type ******" | tee -a $LOG

  #Run refmac with and without anisotropic B-factors.
  foreach TYPETEST (`echo "ANISOT ISOT"`)

    #Return label for job launching
anisooriso:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log
    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      echo "-Testing ${TYPETEST}ropic B-factors" | tee -a $LOG

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
      XYZOUT $WORKDIR/${PDBID}_${TYPETEST}ropic.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_${TYPETEST}ropic.mtz \
      $LIBLIN \
      $SCATLIN \
<<eof >$WORKDIR/${PDBID}_${TYPETEST}ropic.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        make newligand continue
        $TBCMD
        refi type REST resi MLKF meth CGMAT bref $TYPETEST
        $REFIRES
        ncyc 50
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG AUTO 2.50
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $NUCRCMD
        $RESTCMD
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy expdta
        pdbout copy remarks 200 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto anisooriso
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists; then clean up
  foreach TYPETEST (`echo "ANISOT ISOT"`)
    if (! -e $WORKDIR/${PDBID}_${TYPETEST}ropic.mtz) then
      echo "o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in first B-factor type selection" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                         >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      rm $WORKDIR/${PDBID}_${TYPETEST}ropic.mtz
      rm $WORKDIR/${PDBID}_${TYPETEST}ropic.pdb
    endif
  end

  #Compare refinement results and decide on the B-factor type. PROGRAM: bselect
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.bselect.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  $TOOLS/bselect -v $WORKDIR/${PDBID}_ANISOTropic.log $WORKDIR/${PDBID}_ISOTropic.log > $WORKDIR/bselect.log 
  set BREFTYPE = `tail -n 1 $WORKDIR/bselect.log`
  
  #Report on the individual tests
  if (`grep -c 'Taking the simplest model:' $WORKDIR/bselect.log` > 0 || `grep -A 9 'Hamilton test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $4}'` != 'NULL') then
    if (`grep -c 'R-factor gap very similar in the' $WORKDIR/bselect.log` > 0) then 
      echo " o Hamilton test overruled: ${BREFTYPE}tropic B-factors selected" | tee -a $LOG
    else
      echo " o Conclusive Hamilton test: ${BREFTYPE}tropic B-factors selected" | tee -a $LOG
    endif  
  else
    #Check the other tests (R-free Z-score)
    echo " o Inconclusive Hamilton test, checking for overfitting:" | tee -a $LOG
    set TBTYPE = `grep -A 9 'R-free/R Z-score test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $6}'`
    if ($TBTYPE == 'NULL') then
      echo "   * R-free/R Z-score test inconclusive" | tee -a $LOG
    else
      echo "   * R-free/R Z-score test: ${TBTYPE}tropic B-factors preferred" | tee -a $LOG
    endif 
    #R-factor gap
    set TBTYPE = `grep -A 5 'R-free difference test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $4}'`
    if ($TBTYPE == 'NULL') then
      echo "   * R-factor gap test inconclusive" | tee -a $LOG
    else
      echo "   * R-factor gap test: ${TBTYPE}tropic B-factors preferred" | tee -a $LOG
    endif 
    #Exchange rate
    set TBTYPE = `grep -A 6 'Exchange rate test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $4}'`
    if ($TBTYPE == 'BOTH') then
      echo "   * Exchange rate test inconclusive" | tee -a $LOG
    else
      echo "   * Exchange rate test: ${TBTYPE}tropic B-factors preferred" | tee -a $LOG
    endif 
    echo " o Used combined tests to select ${BREFTYPE}tropic B-factors" | tee -a $LOG
  endif
endif

#Setup the refinement
if ($BREFTYPE == "ISOT") then
  #Set number of parameters per atom
  set PPATM = "4"

  #Use TLS unless it is surpressed
  if ($DOTLS == 1) then
    #Set command file
    if(`ls $WORKDIR/????.tls | wc -l` != 0) then
      set TLSCMD  = `echo refi tlsc $TLSCYCLE`

      #Keep only the TLS group definitions, recalculate the origin and the tensors
      if (-e $WORKDIR/${PDBID}.tls) then
        cp $WORKDIR/${PDBID}.tls $WORKDIR/${PDBID}.tls_original
        grep -E 'TLS|RANGE|^ ' $WORKDIR/${PDBID}.tls_original > $WORKDIR/${PDBID}.tls
      endif
    endif
  endif

else
  #Set number of parameters per atom
  set PPATM = "9"

  #Do not use TLS
  set DOTLS    = 0  #No TLS at all
  cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb

  #Do extra cycles in B-restraint weight optimisation
  set BTESTCYC = 15
endif


################################################### Show refinement setting ##############################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Refinement settings ******" | tee -a $LOG

#Type of refinement
echo "-B-factor model" | tee -a $LOG
echo " o Number of atoms      : $ATMCNT"  | tee -a $LOG
echo " o Number of reflections: $WORKCNT" | tee -a $LOG
echo " o Reflections per atom : $REFPATM" | tee -a $LOG
echo " o B-factor type        : ${BREFTYPE}ropic" | tee -a $LOG
if ($DOTLS == 0) then
  echo " o TLS-models           : not used" | tee -a $LOG
else
  echo " o TLS-models           : used" | tee -a $LOG
endif

#NCS restraint type
echo "-Non-crystallographic symmetry" | tee -a $LOG
if ($DONCS == 1 && $STRICTNCS == 1) then
  echo " o Using local NCS restraints and global NCS constraints" | tee -a $LOG
else if ($DONCS == 1) then
  echo " o Using local NCS restraints" | tee -a $LOG
else if ($STRICTNCS == 1) then
  echo " o Using NCS constraints" | tee -a $LOG
else
  echo " o No NCS restraints/constraints used" | tee -a $LOG
endif

#Twinning
echo "-Twinning" | tee -a $LOG
if ($DOTWIN == 1) then
  if ($TWIN == 'twin') then
    echo " o Detwinning during refinement" | tee -a $LOG
  else
    echo " o No twinning detected" | tee -a $LOG
  endif
else
  echo " o Detwinning surpressed"  | tee -a $LOG
endif

#No R-free was reported so a new target must be defined. The expected R-free is used.
if ($RFHEAD == 0) then
  echo "-No R-free reported in PDB header" | tee -a $LOG
  echo " o R-free target is now $RFCALUNB" | tee -a $LOG
endif

#Warn for high structure RMSZ scores
if ($RMSZB == 'huge' || $RMSZA == 'huge' || `echo $RMSZB $RMSZA | awk '{if ($1 > 1.000) {print "1"} else if ($2 > 1.000) {print "1"} else {print "0"}}'` == 1) then
  echo "-High Refmac bond length or bond angle RMSZ detected" | tee -a $LOG
  echo " o Bond length RMSZ: $RMSZB"                          | tee -a $LOG
  echo " o Bond angle RMSZ : $RMSZA"                          | tee -a $LOG
  echo " o Geometric targets will be relaxed"                 | tee -a $LOG
endif

#Check for strangely low R-free values (only when the original R-free set was used)
if ($GOTR == 1) then
  #Check for high Z-score or R-free lower than a certain cut-off (0.0 or R+0.33(Rfree-R))
  if ($RHEAD == 1 && $RFHEAD == 1) then
    set ZCALERR = `echo $ZRFRRATCAL $RFCAL $RCAL $RFREE $RFACT | awk '{if ($1 > 2.6) {print "1"} else if (($2 - $3) < 0.33*($4 - $5)) {print "1"} else {print "0"}}'`
  else
    set ZCALERR = `echo $ZRFRRATCAL $RFCAL $RCAL | awk '{if ($1 > 2.6) {print "1"} else if (($2 - $3) < 0.0) {print "1"} else {print "0"}}'`
  endif
  if ($ZCALERR == 1) then
    echo "-Severe R-free/R ratio bias and possible test set problem!" | tee -a $LOG
    #Do 10 cycles extra refinement
    set TCYCLE = $NCYCLE
    @ NCYCLE = ($NCYCLE + 10)
    echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
  endif
else
  set ZCALERR = 0
endif

#Compensate for new R-free set
if ($GOTR == 0) then
  echo "-New R-free set" | tee -a $LOG
  echo " o Calculated R and expected R-free values will be used for reference" | tee -a $LOG

  #Use 0.5*Wilson B-factor unless the resolution lower than 4.00A; maximise at 50A^2.
  set BCMD = "bfac set $BSET"
  echo " o Resetting B-factors to $BSET to remove model bias" | tee -a $LOG

  #Do 10 restrained refinement cycles extra
  set TCYCLE = $NCYCLE
  @ NCYCLE = ($NCYCLE + 10)
  echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
endif

#Compensate for legacy model
if ($LEGACY == 1 ) then
  echo "-Legacy mode" | tee -a $LOG
  echo " o Using calculated R-factor ($RCAL) as refinement target" | tee -a $LOG

  #Do 10 restrained refinement cycles extra
  set TCYCLE = $NCYCLE
  @ NCYCLE = ($NCYCLE + 10)
  echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
endif

#Compensate for anisotropic B-factors
if ($BREFTYPE == "ANISOT") then
  echo "-Anisotropic B-factor refinement" | tee -a $LOG
  #Do extra cycles of refinement
  set TCYCLE = $NCYCLE
  @ NCYCLE = ($NCYCLE + 20)
  echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
  echo " o Updating R-factor and R-free details"
endif

#Warn about occupancy refinement
if ($OCCREF > 0) then
  echo "-Occupancy refinement" | tee -a $LOG
  echo " o Performing occupancy refinement on $OCCREF residues" | tee -a $LOG
endif

#warn about external restraints
if ($METALCMD != "" || $NUCRCMD != "" || $RESTCMD != "" || $HBONDCMD != "" || $HOMOLCMD != "" || `echo $JELLY | grep -c ridg` != 0 || `echo $HARMCMD | grep -c ridg` != 0) then
  echo "-Additional restraints" | tee -a $LOG
  if (`echo $JELLY | grep -c ridg` != 0) then
    echo " o Using jelly-body restraints" | tee -a $LOG
  endif
  if (`echo $HARMCMD | grep -c ridg` != 0) then
    echo " o Using harmonic restraints" | tee -a $LOG
  endif
  if (-e $WORKDIR/metal.rest) then
    set NMETALREST = `grep -c exte $WORKDIR/metal.rest`
    if ($NMETALREST > 0) then
      echo " o Using $NMETALREST metal site restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/nucleic.rest) then
    set NNUCLEICREST = `grep -c exte $WORKDIR/nucleic.rest`
    if ($NNUCLEICREST > 0) then
      echo " o Using $NNUCLEICREST DNA/RNA restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/homology.rest) then
    set NHOMOLREST = `grep -c exte $WORKDIR/homology.rest`
    if ($NHOMOLREST > 0) then
      echo " o Using $NHOMOLREST homology restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/hbond.rest) then
    set NHBONDREST = `grep -c exte $WORKDIR/hbond.rest`
    if ($NHBONDREST > 0) then
      echo " o Using $NHBONDREST hydrogen bond restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/external.rest) then
    set NEXTERNALREST = `grep -c exte $WORKDIR/external.rest`
    #Subtract 2 instoduced statements
    @ NEXTERNALREST = ($NEXTERNALREST - 2)
    if ($NEXTERNALREST > 0) then
      echo " o Using $NEXTERNALREST user-provided restraints" | tee -a $LOG
    endif
  endif
endif

#Update the R(-free) statistics if the data/parameter ratio changed
if ($BREFTYPE == "ANISOT" || $URESO != $RESOLUTION) then
  #Calculate sigma(R-free)
  set SIGRFCAL = `echo $RFCAL $NTSTCNT | awk '{printf ("%.4f\n", $1/sqrt($2))}'`
  
  #B-factor model can be isotropic or anisotropic  
  if ($BREFTYPE == "ANISOT") then
    set PPATM = 9 #Anisotropic B-factors now
    #Calculate expected R-free/R ratio (Ticke model and anisotropic empirical model 2020)
    set RFRRAT2    = `echo $ATMCNT $WORKCNT $PPATM $URESO | awk '{X = $1/$2; Z = 1 - $3 * X/(1 + 2.5 * X); if (Z < 0) {Z = 0}; A = $3-2.5*(1-Z*Z*Z*Z*Z); X = A*X; RATIO = 1-X; if(RATIO <= 0) {RATIO = 1.010} else {RATIO = sqrt((1+X)/RATIO)};  if ($4 > 2.65 && RATIO > 1.200000) {RATIO = 1.200000}; if ($4 > 3.0 && RATIO < 1.011) {RATIO = 1.20}; if (RATIO > 1.454) {RATIO = 1.45}; printf ("%.4f\n", RATIO)}'`
    set RFRRAT     = `echo $ATMCNT $WORKCNT | awk '{printf ("%.4f\n", -0.14 * log($2/$1) + 1.6674)}'`
    set SRFRRAT    = `echo $ATMCNT $WORKCNT | awk '{printf ("%.4f\n", -0.018 * log($2/$1) + 0.0979)}'`
    set ZRFRRATCAL = `echo $ATMCNT $WORKCNT $RCAL $RFCAL | awk '{printf ("%.2f\n", (((-0.0046 * $2/$1 + 1.3375) - $4/$3)/(-0.018 * log($2/$1) + 0.0979)))}'`
  else  
    #Calculate expected R-free/R ratio (Ticke model and isotropic empirical model 2020)
    set RFRRAT2    = `echo $ATMCNT $WORKCNT $PPATM $URESO | awk '{X = $1/$2; Z = 1 - $3 * X/(1 + 2.5 * X); if (Z < 0) {Z = 0}; A = $3-2.5*(1-Z*Z*Z*Z*Z); X = A*X; RATIO = 1-X; if(RATIO <= 0) {RATIO = 1.010} else {RATIO = sqrt((1+X)/RATIO)};  if ($4 > 2.65 && RATIO > 1.200000) {RATIO = 1.200000}; if ($4 > 3.0 && RATIO < 1.011) {RATIO = 1.20}; if (RATIO > 1.454) {RATIO = 1.45}; printf ("%.4f\n", RATIO)}'`
    set RFRRAT     = `echo $ATMCNT $WORKCNT | awk '{printf ("%.4f\n", -0.005 * $2/$1 + 1.2286)}'`
    set SRFRRAT    = `echo $ATMCNT $WORKCNT | awk '{printf ("%.4f\n", -0.024 * log($2/$1) + 0.1064)}'`
    set ZRFRRATCAL = `echo $ATMCNT $WORKCNT $RCAL $RFCAL | awk '{printf ("%.2f\n", (((-0.005 * $2/$1 + 1.2286) - $4/$3)/(-0.024 * log($2/$1) + 0.1064)))}'`
  endif
  
  #Calculate expected R-free
  set RFCALUNB = `echo $RCAL $RFRRAT | awk '{printf ("%.4f\n", $1*$2)}'`

  #R-free Z-score (comparison of R-free with its expected value)
  set RFCALZ = `echo $RFCALUNB $RFCAL $SIGRFCAL | awk '{printf ("%.2f\n", ($1-$2)/$3)}'`

  #Print values
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** R-factor and R-free details (updated) ******" | tee -a $LOG
  echo "Calculated R      : $RCAL"       | tee -a $LOG
  echo "Calculated R-free : $RFCAL"      | tee -a $LOG
  echo "Expected R-free/R : $RFRRAT (new)"  | tee -a $LOG
  echo "Expected R-free/R : $RFRRAT2 (old)" | tee -a $LOG
  echo "Expected R-free   : $RFCALUNB"   | tee -a $LOG
  echo " " | tee -a $LOG
  echo "sigma(R-free/R)   : $SRFRRAT"    | tee -a $LOG
  echo "R-free/R Z-score  : $ZRFRRATCAL" | tee -a $LOG
  echo " " | tee -a $LOG
  echo "sigma(R-free)     : $SIGRFCAL"   | tee -a $LOG
  echo "R-free Z-score    : $RFCALZ"     | tee -a $LOG 
endif

################################################ TLS group optimisation ##################################################

#Only do this if TLS refinement isn't surpressed
if ($DOTLS == 1) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** TLS group optimisation ******" | tee -a $LOG

  #Optimise TLS groups iff any definitions exist.
  set NTLS = `find $WORKDIR -name "????.tls" | wc -l`
  if ($NTLS != 0) then
    echo "-Testing $NTLS TLS group definition(s)" | tee -a $LOG

    #Maximise the used CPU time
    limit cputime 24h

    #Run refmac with different TLS group configurations
    foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)

      #Return label for job launching
ttestrunning:

      #Only launch new jobs when the number of cores is not exceeded
      #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
      jobs > $WORKDIR/jobs.log
      if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

        #Run refmac
        refmac5 \
        XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
        XYZOUT $WORKDIR/${PDBID}_ttest$TLSG.pdb \
        HKLIN  $WORKDIR/$PDBID.mtz \
        HKLOUT $WORKDIR/${PDBID}_ttest$TLSG.mtz \
        $LIBLIN \
        TLSIN $WORKDIR/$TLSG.tls TLSOUT $WORKDIR/${TLSG}_out.tls \
        $SCATLIN \
<<eof >& $WORKDIR/${PDBID}_ttest$TLSG.log &
          $SCATTERCMD
          make check NONE
          make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
           ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
          make newligand continue
          $TBCMD
          refi type REST resi MLKF meth CGMAT bref ISOT
          $REFIRES
          $TLSCMD
          tlsd waters exclude
          ncyc 0
          scal type $SOLVENT $SCALING
          solvent YES
          $MASKPAR
          $LOWMEM
          weight AUTO
          monitor MEDIUM -
            torsion 10.0 distance 10.0 angle 10.0 plane 10.0 chiral 10.0 -
            bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
          $NCSTYPE
          $NCSALIGN
          $NCSNEIGH
          $NCSSTRICT
          $TWIN
          blim 2.0 999.0
          labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
          $ANOMCMD
          pdbout copy expdta
          pdbout copy remarks 200 280 350
          NOHARVEST
          $HOMOLWGTCMD
          $HOMOLCMD
          $HBONDWGTCMD
          $HBONDCMD
          kill $TOOLS/pdb_redo.refmac
          END
eof
      else
        #Wait a bit to start again
        sleep 10
        goto ttestrunning
      endif
    end

    #Wait for the jobs to finish
    wait

    #Unset the CPU time limit
    limit cputime unlimited

    #Check for errors and report results
    foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)
      if (! -e $WORKDIR/${TLSG}_out.tls) then
        echo " o Problem with Refmac using $TLSG.tls" | tee -a $LOG
        mv $WORKDIR/$TLSG.tls $WORKDIR/$TLSG.notls
        @ NTLS = $NTLS - 1 
      else
        #Mine out the R-free and report
        set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$TLSG.log | head -n 1 | awk '{print $3}'`
        echo " o Tested tls groups from $TLSG.tls (R-free = $TFREE)" | tee -a $LOG
      endif
    end

    #Remove TLS group definitions that cause crazy values in the tensors or B-factors
    if ($NTLS != 0) then
      echo "-Filtering results" | tee -a $LOG
      foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)
        if (`grep -c -E '\*{6}' ${TLSG}_out.tls` != 0) then
          echo " o Problem with TLS group definition $TLSG.tls" | tee -a $LOG
          mv $WORKDIR/$TLSG.tls $WORKDIR/$TLSG.notls
          @ NTLS = $NTLS - 1
        endif
        if (-e $WORKDIR/$TLSG.tls && `grep '^[AH][TE]' ${PDBID}_ttest${TLSG}.pdb | grep -c -E '\*{6}'` != 0) then
          echo " o Problem with TLS group definition $TLSG.tls" | tee -a $LOG
          mv $WORKDIR/$TLSG.tls $WORKDIR/$TLSG.notls
          @ NTLS = $NTLS - 1
        endif
      end
    endif  

    #Create .ttest file for picker iff any TLS groups are approved
    if ($NTLS != 0) then
      #Grep out the refinement statistics, and take the first data line to obtain R(-free) targets. Use the first valid log file.
      set VALID1 = `ls -Sr ????.tls | head -n 1 | cut -c 1-4`
      set TTR  = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_ttest$VALID1.log | tail -n 1 | awk '{print $2}'`
      set TTRF = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_ttest$VALID1.log | tail -n 1 | awk '{print $3}'`

      echo "$TTR $TTRF" > $WORKDIR/${PDBID}.ttest   #Use R and R-free obtained after resetting the B-factor
      foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)
        if (-e $WORKDIR/${PDBID}_ttest$TLSG.pdb) then
          if ($TLSG == $VALID1) then
            set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$TLSG.log | head -n 1`
            echo "$TLSG $LINE" >> $WORKDIR/${PDBID}.ttest
          else
            #Do a hamilton test first
            cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.bselect.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
            
            if (`$TOOLS/bselect -t $WORKDIR/${PDBID}_ttest$VALID1.log $WORKDIR/${PDBID}_ttest$TLSG.log` != 'LOG1') then
              set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$TLSG.log | head -n 1`
              echo "$TLSG $LINE" >> $WORKDIR/${PDBID}.ttest
            else
              echo " o TLS group definition $TLSG.tls causes over-fitting" | tee -a $LOG
            endif
          endif
        else
          echo " o $WORKDIR/${PDBID}_ttest$TLSG.pdb was missing" | tee -a $LOG
          set TTEST = 1
        endif
      end

      #Pick the best TLS group configuration. PROGRAM: picker
      set OPTTLSG = `$TOOLS/picker -z $WORKDIR/${PDBID}.ttest $NTSTCNT $RFRRAT $SRFRRAT` #The geometry is ignored here
      
      #Also pick second best TLS group.
      if (`ls -Sr ????.tls | wc -l` != 1) then
        grep -v $OPTTLSG $WORKDIR/${PDBID}.ttest > $WORKDIR/${PDBID}.ttest2
        set OPTTLSG2 = `$TOOLS/picker -z $WORKDIR/${PDBID}.ttest2 $NTSTCNT $RFRRAT $SRFRRAT`
      else
        set OPTTLSG2 = 'none'
      endif

    
      #Print values
      if ($OPTTLSG == 'none') then
        echo "-TLS does not seem to work for this structure" | tee -a $LOG
        echo " o TLS refinement will not be used"  | tee -a $LOG
        #Set up input for B-weight optimisation
        cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
        set DOTLS    = 0
        set TLSCMD   =        #Empty TLS command
        set TLSFILS  =        #No TLS files specified
      else
        echo "-TLS groups in $OPTTLSG.tls may be used in refinement" | tee -a $LOG
        #Set up input for B-weight optimisation
        cp $WORKDIR/${OPTTLSG}_out.tls         $WORKDIR/${PDBID}_refin.tls
        set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_refin.tls`
        set BCMD    =  #No need to reset the B-factors anymore
        #Run TLSanl to get a proper input PDB file if an older Refmac was used
        cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.tlsanl.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
        
        tlsanl \
        XYZIN $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb \
        XYZOUT $WORKDIR/${PDBID}_refin.pdb \
<<eof >& $WORKDIR/${PDBID}_tlsanl.log
        BINPUT t
        BRESID f
        ISOOUT RESI
        NUMERICAL
        END
eof
        if($status) then
          echo " o Problem with TLSanl for TLS groups in $OPTTLSG.tls" | tee -a $LOG
          echo "COMMENT: TLSanl: general error for optimal TLS group" >> $DEBUG
          echo "PDB-REDO,$PDBID"                                      >> $DEBUG
          #Try using the second best TLS group configuration
          if ($OPTTLSG2 != 'none') then

            set OPTTLSG = $OPTTLSG2
            echo " o TLS groups in $OPTTLSG.tls may be used in refinement instead" | tee -a $LOG

            #Run TLSanl again to get a proper input PDB file
            tlsanl \
            XYZIN $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb \
            XYZOUT $WORKDIR/${PDBID}_refin.pdb \
<<eof >$WORKDIR/${PDBID}_tlsanl2.log
            BINPUT t
            BRESID f
            ISOOUT RESI
            NUMERICAL
            END
eof
            if($status) then #second failure
              echo " o Problem with TLSanl for TLS groups in $OPTTLSG.tls" | tee -a $LOG
              echo "COMMENT: TLSanl: general error for second best TLS group" >> $DEBUG
              echo "PDB-REDO,$PDBID"                                          >> $DEBUG
              cp $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb $WORKDIR/${PDBID}_refin.pdb
              #The output files will have the total B-factor instead of the residual. Circumvent by resetting the B-factors.
              set BCMD = `echo $TBCMD`
            endif

            #Set up input for B-weight optimisation (again)
            cp $WORKDIR/${OPTTLSG}_out.tls         $WORKDIR/${PDBID}_refin.tls
            set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_refin.tls`

            
          else
            rm $WORKDIR/${PDBID}.ttest2
            cp $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb $WORKDIR/${PDBID}_refin.pdb
            #The output files will have the total B-factor instead of the residual. Circumvent by resetting the B-factors.
            set BCMD = `echo $TBCMD`
          endif
        endif
        
        #Test the refinement stability with and without TLS
        echo "-Testing refinement performance" | tee -a $LOG
        
        #Run refmac with different TLS group configurations
        foreach TMODE (TLSY TLSN)
   
          if ($TMODE == TLSY) then
            set TLSLIN = "TLSIN $WORKDIR/$OPTTLSG.tls TLSOUT $WORKDIR/tmodeTLSY_out.tls"
          else
            #No TLS whatsoever
            set TLSLIN = 
            set TLSCMD =
          endif
          
          #Return label for job launching
tmoderunning:

          #Only launch new jobs when the number of cores is not exceeded
          #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
          jobs > $WORKDIR/jobs.log
          if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

            #Run refmac
            refmac5 \
            XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
            XYZOUT $WORKDIR/${PDBID}_tmode$TMODE.pdb \
            HKLIN  $WORKDIR/$PDBID.mtz \
            HKLOUT $WORKDIR/${PDBID}_tmode$TMODE.mtz \
            $LIBLIN \
            $TLSLIN \
            $SCATLIN \
<<eof >& $WORKDIR/${PDBID}_tmode$TMODE.log &
              $SCATTERCMD
              make check NONE
              make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
                ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
              make newligand continue
              $TBCMD
              refi type REST resi MLKF meth CGMAT bref ISOT
              $REFIRES
              $TLSCMD
              tlsd waters exclude
              ncyc 20
              scal type $SOLVENT $SCALING
              solvent YES
              $MASKPAR
              $LOWMEM
              weight $WGTSIG AUTO 2.50
              monitor MEDIUM -
                torsion 10.0 distance 10.0 angle 10.0 plane 10.0 chiral 10.0 -
                bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
              $NCSTYPE
              $NCSALIGN
              $NCSNEIGH
              $NCSSTRICT
              $TWIN
              blim 2.0 999.0
              labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
              $ANOMCMD
              pdbout copy expdta
              pdbout copy remarks 200 280 350
              NOHARVEST
              $HOMOLWGTCMD
              $HOMOLCMD
              $HBONDWGTCMD
              $HBONDCMD
              kill $TOOLS/pdb_redo.refmac
              END
eof
          else
            #Wait a bit to start again
            sleep 10
            goto tmoderunning
          endif
        end

        #Wait for the jobs to finish
        wait
    
        #Analyse the results benchmark against TLS refinement only
        set TMR  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$OPTTLSG.log | head -n 1 | awk '{print $2}'`
        set TMRF = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$OPTTLSG.log | head -n 1 | awk '{print $3}'`
        echo "$TMR $TMRF" > $WORKDIR/${PDBID}.tmode
        foreach TMODE (TLSY TLSN)
          set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_tmode$TMODE.log | head -n 1 | awk '{print $3}'`
          if ($TMODE == TLSY) then
            echo " o Tested refinement with TLS    (R-free = $TFREE)" | tee -a $LOG
          else
            echo " o Tested refinement without TLS (R-free = $TFREE)" | tee -a $LOG
          endif  
          set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_tmode$TMODE.log | head -n 1`
          echo "$TMODE $LINE" >> $WORKDIR/${PDBID}.tmode
        end
        set BTMODE = `$TOOLS/picker -f $WORKDIR/${PDBID}.tmode $NTSTCNT $RFRRAT $SRFRRAT $RMSZB $RMSZA`
    
        # make the call
        if ($BTMODE == TLSY) then
          echo "-Model refinement with TLS works best" | tee -a $LOG 
        else
          echo "-Model refinement without TLS works best" | tee -a $LOG
          set DOTLS   = 0  #Do not use TLS
          set TLSCMD  =    #Empty TLS command
          set TLSFILS =    #No TLS files specified
          cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
        endif
    
      endif
    else
      #All TLS definitions failed
      echo "-No usable TLS group definitions"   | tee -a $LOG
      echo " o TLS refinement will not be used" | tee -a $LOG
      set DOTLS   = 0  #Do not use TLS
      set TLSCMD  =    #Empty TLS command
      set TLSFILS =    #No TLS files specified
      cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
    endif
  else
    #No TLS groups could be made
    echo "-Cannot create any TLS group definitions" | tee -a $LOG
    echo " o TLS refinement will not be used"       | tee -a $LOG
    echo "COMMENT: Could not create any TLS groups" >> $DEBUG
    echo "PDB-REDO,$PDBID"                          >> $DEBUG
    set DOTLS   = 0  #Do not use TLS
    set TLSCMD  =    #Empty TLS command
    set TLSFILS =    #No TLS files specified
    cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
  endif

else
  #Not using TLS ensure that the input file for the B-weight optimisation exists
  cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
endif

##################################### Individual B-factors or one overall B-factor? ######################################

#Try refinement with one overall B if TLS works and there are less than 4 reflections per atom.
if ( ($STRICTNCS == 1 && $NMTRIX > 9) || $BISOT == 1) then
  set BTYPE = "ISOT"
else
  set BTYPE = `echo $WORKCNT $ATMCNT | awk '{if ($1/$2 < 4) {print "TEST"} else {print "ISOT"}}'`
endif

#Test the best B-factor restraint type
if ($BTYPE == TEST) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Testing B-factor refinement type ******" | tee -a $LOG

  #Set the B-factor restraint weight
  set TIGHTB = `echo $BWEIGHTS | cut -c 1-4`

  #Set all B-factors to a single value if no TLS is used
  if ($DOTLS == 0) then
    set BBCMD = `echo $TBCMD`
  else
    set BBCMD =  #No B-factor resetting (the B-factors were reset during the TLS optimisation)
  endif

  #Run refmac with and without individual B-factors.
  foreach TYPETEST (`echo "ISOT OVER"`)

    #Return label for job launching
isoorover:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log
    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      #Report the type of refinement
      if ($TYPETEST == "ISOT") then
        echo "-Testing isotropic B-factors"  | tee -a $LOG
      else
        echo "-Testing one overall B-factor" | tee -a $LOG
      endif

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_refin.pdb \
      XYZOUT $WORKDIR/${PDBID}_${TYPETEST}.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_${TYPETEST}.mtz \
      $LIBLIN \
      $TLSFILS \
      $SCATLIN \
<<eof >$WORKDIR/${PDBID}_${TYPETEST}.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        make newligand continue
        $BBCMD
        refi type REST resi MLKF meth CGMAT bref $TYPETEST
        $REFIRES
        ncyc $NCYCLE
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG AUTO 2.50
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $NUCRCMD
        $RESTCMD
        temp $TIGHTB
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy expdta
        pdbout copy remarks 200 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto isoorover
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists; then clean up
  foreach TYPETEST (`echo "ISOT OVER"`)
    if (! -e $WORKDIR/${PDBID}_${TYPETEST}.mtz) then
      echo "o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in second B-factor type selection" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                          >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      rm $WORKDIR/${PDBID}_${TYPETEST}.mtz
      rm $WORKDIR/${PDBID}_${TYPETEST}.pdb
    endif
  end

  #Compare refinement results and decide on the B-factor type
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.bselect.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  $TOOLS/bselect -v $WORKDIR/${PDBID}_ISOT.log $WORKDIR/${PDBID}_OVER.log > $WORKDIR/bselect.log 
  set BREFTYPE = `tail -n 1 $WORKDIR/bselect.log`
  
  #Report on the individual tests
  if (`grep -c 'Taking the simplest model:' $WORKDIR/bselect.log` > 0 || `grep -A 9 'Hamilton test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $4}'` != 'NULL') then
    if (`grep -c 'R-factor gap very similar in the' $WORKDIR/bselect.log` > 0) then 
      echo " o Hamilton test overruled: isotropic B-factors selected" | tee -a $LOG
    else if ($BREFTYPE == "ISOT") then
      echo " o Conclusive Hamilton test: isotropic B-factors selected" | tee -a $LOG
    else  
      echo " o Conclusive Hamilton test: one overall B-factor selected" | tee -a $LOG
    endif  
  else
    #Check the other tests (R-free Z-score)
    echo " o Inconclusive Hamilton test, checking for overfitting:" | tee -a $LOG
    set TBTYPE = `grep -A 9 'R-free/R Z-score test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $6}'`
    if ($TBTYPE == 'NULL') then
      echo "   * R-free/R Z-score test inconclusive" | tee -a $LOG
    else if ($TBTYPE == "ISOT") then
      echo "   * R-free/R Z-score test: isotropic B-factors preferred" | tee -a $LOG
    else  
      echo "   * R-free/R Z-score test: one overall B-factor preferred" | tee -a $LOG
    endif 
    #R-factor gap
    set TBTYPE = `grep -A 5 'R-free difference test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $4}'`
    if ($TBTYPE == 'NULL') then
      echo "   * R-factor gap test inconclusive" | tee -a $LOG
    else if ($TBTYPE == "ISOT") then
      echo "   * R-factor gap test: isotropic B-factors preferred" | tee -a $LOG
    else  
      echo "   * R-factor gap test: one overall B-factor preferred" | tee -a $LOG      
    endif 
    #Exchange rate
    set TBTYPE = `grep -A 6 'Exchange rate test:' $WORKDIR/bselect.log | grep 'Best B-factor' | awk '{print $4}'`
    if ($TBTYPE == 'BOTH') then
      echo "   * Exchange rate test inconclusive" | tee -a $LOG
    else if ($TBTYPE == "ISOT") then
      echo "   * Exchange rate test: isotropic B-factors preferred" | tee -a $LOG
    else
      echo "   * Exchange rate test: one overall B-factor preferred" | tee -a $LOG
    endif 
    #Show the conclusion
    if ($BREFTYPE == "ISOT") then
      echo " o Used combined tests to select isotropic B-factors" | tee -a $LOG
    else
      echo " o Used combined tests to select one overall B-factor" | tee -a $LOG
    endif
  endif
  
  #Setup the refinement
  if ($BREFTYPE == "OVER") then
    echo "-PDB-REDO will use one overall B-factor in refinement" | tee -a $LOG

    #Set number of parameters per atom
    set PPATM = "3"

    #Calculate sigma(R-free)
    set SIGRFCAL = `echo $RFCAL $NTSTCNT | awk '{printf ("%.4f\n", $1/sqrt($2))}'`
    
    #Set all B-factors to a single value if no TLS is used
    if ($DOTLS == 0) then
      set BCMD = `echo $TBCMD`
    endif

    #Re-calculate R-free/R ratio
    set RFRRAT2    = `echo $ATMCNT $WORKCNT $PPATM $URESO | awk '{X = $1/$2; Z = 1 - $3 * X/(1 + 2.5 * X); if (Z < 0) {Z = 0}; A = $3-2.5*(1-Z*Z*Z*Z*Z); X = A*X; RATIO = 1-X; if(RATIO <= 0) {RATIO = 1.010} else {RATIO = sqrt((1+X)/RATIO)};  if ($4 > 2.65 && RATIO > 1.200000) {RATIO = 1.200000}; if ($4 > 3.0 && RATIO < 1.011) {RATIO = 1.20}; if (RATIO > 1.454) {RATIO = 1.45}; printf ("%.4f\n", RATIO)}'`
    set RFRRAT     = 1.1700 #Simple flat model due to limited data
    set SRFRRAT    = 0.0850 #Simple flat model due to limited data
    set ZRFRRATCAL = `echo $RCAL $RFCAL | awk '{printf ("%.2f\n", (1.1170 - ($2/$1))/0.0850)}'`
    
    #Calculate expected R-free
    set RFCALUNB = `echo $RCAL $RFRRAT | awk '{printf ("%.4f\n", $1*$2)}'`

    #R-free Z-score (comparison of R-free with its expected value)
    set RFCALZ = `echo $RFCALUNB $RFCAL $SIGRFCAL | awk '{printf ("%.2f\n", ($1-$2)/$3)}'`
    
    #Print values
    echo " " | tee -a $LOG
    echo " " | tee -a $LOG
    echo "****** R-factor and R-free details (updated) ******" | tee -a $LOG
    echo "Calculated R      : $RCAL"       | tee -a $LOG
    echo "Calculated R-free : $RFCAL"      | tee -a $LOG
    echo "Expected R-free/R : $RFRRAT (new)"  | tee -a $LOG
    echo "Expected R-free/R : $RFRRAT2 (old)" | tee -a $LOG
    echo "Expected R-free   : $RFCALUNB"   | tee -a $LOG
    echo " " | tee -a $LOG
    echo "sigma(R-free/R)   : $SRFRRAT"    | tee -a $LOG
    echo "R-free/R Z-score  : $ZRFRRATCAL" | tee -a $LOG
    echo " " | tee -a $LOG
    echo "sigma(R-free)     : $SIGRFCAL"   | tee -a $LOG
    echo "R-free Z-score    : $RFCALZ"     | tee -a $LOG 

  else
    echo "-PDB-REDO will use isotropic B-factors in refinement" | tee -a $LOG
  endif
endif

################################################# B-weight optimisation ##################################################

if ($BREFTYPE == ISOT || $BREFTYPE == ANISOT) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** B-weight optimisation ******" | tee -a $LOG
  echo "-Testing B-factor restraint weights: $BWEIGHTS" | tee -a $LOG

  #Run refmac with predefined B-factor restraint weights
  foreach BWGT (`echo $BWEIGHTS`)

    #Return label for job launching
bwgtrunning:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log

    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      #Run Refmac. Errors are caught later by checking the output.
      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_refin.pdb \
      XYZOUT $WORKDIR/${PDBID}_btest$BWGT.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_btest$BWGT.mtz \
      $LIBLIN \
      $TLSFILS \
      $SCATLIN \
<<eof >$WORKDIR/${PDBID}_btest$BWGT.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        make newligand continue
        $BCMD
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        ncyc $BTESTCYC
        tlsd waters exclude
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG AUTO 2.50
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $NUCRCMD
        $RESTCMD
        temp $BWGT
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy expdta
        pdbout copy remarks 200 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto bwgtrunning
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists
  foreach BWGT (`echo $BWEIGHTS`)
    if (! -e $WORKDIR/${PDBID}_btest$BWGT.mtz) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in B-weight optimisation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                 >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  end

  #Report
  foreach BWGT (`echo $BWEIGHTS`)
    set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_btest$BWGT.log | head -n 1 | awk '{print $3}'`
    echo " o Performed B-weight test with weight $BWGT (R-free = $TFREE)" | tee -a $LOG
  end

  #Create .btest file for picker
  echo "-Selecting best B weight" | tee -a $LOG
  set BTEST = 0
  if ($GOTR == 0 || $ZCALERR == 1) then
    echo "$RCAL $RFCALUNB" > $WORKDIR/${PDBID}.btest  #Use recalculated R and expected R-free as benchmarks
  else
    echo "$RCAL $RFCAL" > $WORKDIR/${PDBID}.btest  #Use R and R-free obtained from the recalculation as benchmarks
  endif
  foreach BWGT (`echo $BWEIGHTS`)
    if (-e $WORKDIR/${PDBID}_btest$BWGT.pdb) then
      set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_btest$BWGT.log | head -n 1`
      echo "$BWGT $LINE" >> $WORKDIR/${PDBID}.btest
    else
      echo " o $WORKDIR/${PDBID}_btest$BWGT.pdb was missing." | tee -a $LOG
      set BTEST = 1
    endif
  end

  #Pick the best weight. If no weight is found 'none' is returned and the re-refinement runs with default settings.
  set BBEST = `$TOOLS/picker -s $ESTRICT $WORKDIR/${PDBID}.btest $NTSTCNT $RFRRAT $SRFRRAT $RMSZB $RMSZA`

  #Print values
  echo " o Best B weight: $BBEST" | tee -a $LOG
else
  set BBEST = overall
  set BTEST = 0
endif

######################################################### re-refinement ##################################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Structure re-refinement ******" | tee -a $LOG
echo "-Refining with geometric restraint weights: $WEIGHTS" | tee -a $LOG


#Run refmac with predefined matrix weights
foreach WGT (`echo $WEIGHTS`)

  #Return label for job launching
refirunning:

  #Only launch new jobs when the number of cores is not exceeded
  #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
  jobs > $WORKDIR/jobs.log

  if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_refin.pdb \
    XYZOUT $WORKDIR/${PDBID}_refmac${WGT}.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_refmac${WGT}.mtz \
    $LIBLIN \
    $TLSFILS \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_${WGT}.log &
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      make newligand continue
      $BCMD
      refi type REST resi MLKF meth CGMAT bref $BREFTYPE
      $REFIRES
      tlsd waters exclude
      ncyc $NCYCLE
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG $WGTTYPE $WGT
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $OCCCMD
      $HARMCMD
      $METALCMD
      $NUCRCMD
      $RESTCMD
      temp $BBEST
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
      $ANOMCMD
      pdbout copy expdta
      pdbout copy remarks 200 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      END
eof
  else
    #Wait a bit to start again
    sleep 10
    goto refirunning
  endif
end

#Wait for the jobs to finish
wait

#Check for problems by seeing if the output mtz file exists
foreach WGT (`echo $WEIGHTS`)
  if (! -e $WORKDIR/${PDBID}_refmac$WGT.mtz) then
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in re-refinement" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                         >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
end

################################################# Find best results ######################################################

#Report
foreach WGT (`echo $WEIGHTS`)
  #Mine out the R-free
  set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_${WGT}.log | head -n 1 | awk '{print $3}'`
  echo " o Performed refinement with $WGTTYPE weight $WGT (R-free = $TFREE)" | tee -a $LOG
end


echo "-Selecting best geometric restraint weight" | tee -a $LOG

#No known errors in the refinement
set TLSERR  = 0

#Create second .refi file
if ($GOTR == 0 || $ZCALERR == 1) then
  echo "$RCAL $RFCALUNB" > $WORKDIR/${PDBID}.refi
else
  echo "$RCAL $RFCAL" > $WORKDIR/${PDBID}.refi
endif
foreach WGT (`echo $WEIGHTS`)
  if (-e $WORKDIR/${PDBID}_refmac${WGT}.pdb) then
    set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_${WGT}.log | head -n 1`
    echo "$WGT $LINE" >> $WORKDIR/${PDBID}.refi
  else
    echo "$WORKDIR/${PDBID}_refmac${WGT}.pdb is missing." | tee -a $LOG
    set TLSERR = 1
  endif
end

#Pick the best re-refined structure
set TLSBEST = `$TOOLS/picker $NEWMODEL $ESTRICT $WORKDIR/${PDBID}.refi $NTSTCNT $RFRRAT $SRFRRAT $RMSZB $RMSZA`

if ($TLSBEST == none) then

  #Set R(-free)
  set RTLS  = $RCAL
  set RFTLS = $RFCAL

  #Copy files
  cp $WORKDIR/${PDBID}_refin.pdb $WORKDIR/${PDBID}_besttls.pdb
  if ($DOTLS == 1) then
    cp $WORKDIR/${PDBID}_ttest$OPTTLSG.mtz $WORKDIR/${PDBID}_besttls.mtz
    cp $WORKDIR/${PDBID}_ttest$OPTTLSG.log $WORKDIR/${PDBID}_besttls.log
  else
    cp $WORKDIR/${PDBID}_0cyc.mtz $WORKDIR/${PDBID}_besttls.mtz
    cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_besttls.log
  endif

  #No use calculating values here
  set SIGRFTLS   = 'NA'
  set RFTLSUNB   = 'NA'
  set RFTLSZ     = 'NA'
  set ZRFRRATTLS = 'NA'

else

  #Copy files
  cp $WORKDIR/${PDBID}_refmac${TLSBEST}.pdb $WORKDIR/${PDBID}_besttls.pdb
  cp $WORKDIR/${PDBID}_refmac${TLSBEST}.mtz $WORKDIR/${PDBID}_besttls.mtz
  cp $WORKDIR/${PDBID}_${TLSBEST}.log $WORKDIR/${PDBID}_besttls.log

  #Set R(-free)
  set RTLS  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_besttls.log | head -n 1 | awk '{print $2}'`
  set RFTLS = `tail -n $LOGSTEP $WORKDIR/${PDBID}_besttls.log | head -n 1 | awk '{print $3}'`

  #Validate R-free and R-free/R
  set SIGRFTLS   = `echo $RFTLS $NTSTCNT | awk '{printf ("%.4f\n", $1/sqrt($2))}'`
  set RFTLSUNB   = `echo $RTLS $RFRRAT   | awk '{printf ("%.4f\n", $1*$2)}'`
  set ZRFRRATTLS = `echo $RFRRAT $SRFRRAT $RTLS $RFTLS | awk '{printf ("%.2f\n", ($1-($4/$3))/$2)}'`
  set RFTLSZ     = `echo $RFTLSUNB $RFTLS $SIGRFTLS | awk '{printf ("%.2f\n", ($1-$2)/$3)}'`
  
  #Update RMSZ targets if they were higher than 1.00 before
  if (`echo $RMSZB $RMSZA | awk '{if ($1 > 1.000) {print "1"} else if ($2 > 1.000) {print "1"} else {print "0"}}'` == 1) then
    set RMSZB = `tail -n $LOGSTEP $WORKDIR/${PDBID}_besttls.log | head -n 1 | sed 's/\*\*\*\*\*\*/huge/g' | awk '{print $8}'`
    set RMSZA = `tail -n $LOGSTEP $WORKDIR/${PDBID}_besttls.log | head -n 1 | sed 's/\*\*\*\*\*\*/huge/g' | awk '{print $10}'`
  endif

endif
#Print values
echo " o Best geometric restraint weight: $TLSBEST" | tee -a $LOG

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Re-refinement details ******" | tee -a $LOG
echo "Best matrix weight: $TLSBEST" | tee -a $LOG
echo "Resulting R-factor: $RTLS"    | tee -a $LOG
echo "Resulting R-free  : $RFTLS"   | tee -a $LOG
echo "R-free/R Z-score  : $ZRFRRATTLS" | tee -a $LOG
echo " " | tee -a $LOG
echo "Expected R-free   : $RFTLSUNB" | tee -a $LOG
echo "sigma(R-free)     : $SIGRFTLS" | tee -a $LOG
echo "R-free Z-score    : $RFTLSZ"   | tee -a $LOG

#Copy back the SEQRES records
if ($GOTSEQRES == 1  && `grep -c '^SEQRES' $WORKDIR/${PDBID}_besttls.pdb` == 0) then
  cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_besttls.bak
  #Run seqrescopier. PROGRAM: seqrescopier
  $TOOLS/seqrescopier -v \
  -pdbinw  $WORKDIR/cache.pdb \
  -pdbinwo $WORKDIR/${PDBID}_besttls.bak \
  -pdbout  $WORKDIR/${PDBID}_besttls.pdb >> $WORKDIR/seqrescopier.log 
else
  #Do nothing
endif 


####################################################### Clean up round 1 #################################################

#Delete files...

#General
rm $WORKDIR/mtz_creation.log
rm $WORKDIR/raw.mtz
rm $WORKDIR/unique.mtz
rm $WORKDIR/merged.mtz
if (-e $WORKDIR/rawbu.mtz) then
  rm $WORKDIR/rawbu.mtz
endif
if (-e $WORKDIR/${PDBID}_0cycv1.log) then
    rm $WORKDIR/${PDBID}_0cycv1.log
    rm $WORKDIR/${PDBID}_0cycv1.pdb
endif

#TLS groups
if ($DOTLS == 1) then
  foreach TLSG (`ls ????.tls | cut -c 1-4`)
    rm $WORKDIR/${TLSG}_out.tls
    rm $WORKDIR/${PDBID}_ttest$TLSG.mtz
  end
endif

#B-weight
if ($BREFTYPE == ISOT || $BREFTYPE == ANISOT) then
  foreach BWGT (`echo $BWEIGHTS`)
    rm $WORKDIR/${PDBID}_btest$BWGT.pdb
    rm $WORKDIR/${PDBID}_btest$BWGT.mtz
    rm $WORKDIR/${PDBID}_btest$BWGT.log
  end
endif

#Re-refinement
foreach WGT (`echo $WEIGHTS`)
  rm $WORKDIR/${PDBID}_refmac${WGT}.pdb
  rm $WORKDIR/${PDBID}_refmac${WGT}.mtz
  rm $WORKDIR/${PDBID}_${WGT}.log
end

############################################   Validate structures  ######################################################

#Start reporting
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Model validation ******" | tee -a $LOG

#Initialise
set WCERR = 0


#Run WHAT_CHECK on the original model "input" and the re-refined model "reref"
foreach STAGE (`echo "0cyc besttls"`)

  #Return label for job launching
wcrunning:

  #Only launch new jobs when the number of cores is not exceeded
  #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
  jobs > $WORKDIR/jobs.log

  if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

    #Create and go to temporary running directory
    set WCWORK = $WORKDIR/wctemp_$STAGE
    mkdir -p $WCWORK
    cd $WCWORK
    cp $WORKDIR/${PDBID}_${STAGE}.pdb $WCWORK/

    #Report
    if ($STAGE == "0cyc") then
      echo "-Validating the input model with WHAT_CHECK" | tee -a $LOG
    else  
      echo "-Validating the re-refined model with WHAT_CHECK" | tee -a $LOG
    endif
  
    #Do the actual validation. PROGRAM: WHAT_CHECK
    $WC/bin/whatcheck $WCWORK/${PDBID}_$STAGE.pdb Y Y Y >& $WORKDIR/wc_$STAGE.log &
    
  else
    #Wait a bit to start again
    sleep 5
    goto wcrunning
  endif     
end  

#Wait for the jobs to finish
wait

#Check the jobs and clean up  
foreach STAGE (`echo "0cyc besttls"`)
  #Check for an output file
  set WCWORK = $WORKDIR/wctemp_$STAGE
  cd $WCWORK
  if (-e $WCWORK/pdbout.txt) then
    #Do nothing
  else
    #Give warning
    echo " o Validation failed" | tee -a $LOG
    set WCERR = 1
    echo "COMMENT: WHAT_CHECK cannot validate $STAGE model" >> $DEBUG
    echo "PDB-REDO,$PDBID"                                  >> $DEBUG
  endif

  #Create webpage
  $WC/bin/pdbout2html >>& $WORKDIR/wc_original.log

  #Check index.html completeness
  if (`grep -c "Final summary" $WCWORK/pdbout.html` == 0 && `grep -c "Summary report" $WCWORK/pdbout.html` == 0) then
    echo " o Validation failed" | tee -a $LOG
    set WCERR = 1
    echo "COMMENT: WHAT_CHECK cannot validate $STAGE model" >> $DEBUG
    echo "PDB-REDO,$PDBID"                                  >> $DEBUG
  endif

  #Setup data move and data mining
  if ($STAGE == "0cyc") then
    set WCDIR = $WORKDIR/wo
  else   
    set WCDIR = $WORKDIR/wc
  endif  
  mkdir -p $WCDIR

  if (-e $WCWORK/pdbout.txt) then
    mv $WCWORK/pdbout.txt $WCDIR/
  endif
  
  #Copy html and images only for the server
  if ($SERVER == 1) then
    mv $WCWORK/pdbout.html $WCDIR/index.html
    mv $WCWORK/*.gif $WCDIR/ >& /dev/null
  endif

  #Go to start directory
  cd $WORKDIR
  rm -rf $WCWORK
end  

#Run tortoize
echo "-Validating the input model with tortoize" | tee -a $LOG
$TOOLS/tortoize --xyzin $WORKDIR/${PDBID}_0cyc.pdb --output $WORKDIR/${PDBID}_0cyc_tortoize.json >& $WORKDIR/tortoize.log
echo "-Validating the rerefined model with tortoize" | tee -a $LOG
$TOOLS/tortoize --xyzin $WORKDIR/${PDBID}_besttls.pdb --output $WORKDIR/${PDBID}_besttls_tortoize.json >>& $WORKDIR/tortoize.log

#Extract statistics from original structure model
set PDBOUT = $WORKDIR/wo/pdbout.txt
set ONATOM = `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_0cyc.pdb`
if (-e $PDBOUT) then
  if (`grep -a -c '1st generation packing quality :' $PDBOUT` == 0) then
    set OZPAK1 = 'NA'
  else
    set OZPAK1 = `grep -a '1st generation packing quality :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c '2nd generation packing quality :' $PDBOUT` == 0) then
    set OZPAK2 = 'NA'
  else
    set OZPAK2 = `grep -a '2nd generation packing quality :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Ramachandran plot appearance   :' $PDBOUT` == 0) then
    set OZRAMA = 'NA'
  else
    set OZRAMA = `grep -a 'Ramachandran plot appearance   :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'chi-1/chi-2 rotamer normality  :' $PDBOUT` == 0) then
    set OCHI12 = 'NA'
  else
    set OCHI12 = `grep -a 'chi-1/chi-2 rotamer normality  :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Backbone conformation          :' $PDBOUT` == 0) then
    set OBCONF = 'NA'
  else
    set OBCONF = `grep -a 'Backbone conformation          :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond lengths                   :' $PDBOUT` == 0) then
    set OBRMSZ = 'NA'
  else
    set OBRMSZ = `grep -a 'Bond lengths                   :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond angles                    :' $PDBOUT` == 0) then
    set OARMSZ = 'NA'
  else
    set OARMSZ = `grep -a 'Bond angles                    :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  #Get bump scores
  if (`grep -a -c 'Total number of bumps:' $PDBOUT` == 0) then
    set OBUMPS = 0
    set OSBMPL = 0
    set OWBMPS = 0.000
  else
    set OBUMPS = `grep -a 'Total number of bumps:' $PDBOUT | cut -c 24-28`
    set OSBMPL = `grep -a 'Total squared bump value:' $PDBOUT | cut -c 27-33`
    set OWBMPS = `echo $ONATOM $OSBMPL | awk '{printf "%5.3f\n", 100*$2/$1}'`
  endif

  #Count unsatisfied donors and acceptors
  if (`grep -a -A 80 'The buried hydrogen bond donors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set OHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set OHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $PDBOUT | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  if (`grep -A 80 'The buried side-chain hydrogen bond acceptors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set OHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set OHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $PDBOUT | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  @ OHBUNS = ($OHBDON + $OHBACC)
  if (`grep -a -c 'Buried donors:' $PDBOUT` == 0) then
    set OHBSAT = 'NA'
  else
    @ OBHBDA = (`grep -a 'Buried donors:' $PDBOUT | awk '{print $3}'` + `grep -a 'Buried acceptors:' $PDBOUT | awk '{print $3}'`)
    if ($OBHBDA == 0) then
      set OHBSAT = 'NA'
    else  
      set OHBSA1 = `grep -a 'with a H-bond:' $PDBOUT | awk '{SUM += $5} END {print SUM}'`
      set OHBSA2 = `grep -a 'with a poor H-bond:' $PDBOUT | awk '{SUM += $6} END {print 0.5*SUM}'`
      set OHBSA3 = `grep -a 'with only a very poor H-bond:' $PDBOUT | awk '{SUM += $8} END {print 0.25*SUM}'`
      set OHBSA4 = `grep -a 'essentially without H-bond:' $PDBOUT | awk '{SUM += $5} END {print 0.125*SUM}'`
      set OHBSAT = `echo $OHBSA1 $OHBSA2 $OHBSA3 $OHBSA4 $OBHBDA | awk '{printf "%5.3f\n" ,($1 + $2 + $3 + $4)/$5}'`
    endif  
  endif
else
  set OZPAK1 = 'NA'
  set OZPAK2 = 'NA'
  set OZRAMA = 'NA'
  set OCHI12 = 'NA'
  set OBCONF = 'NA'
  set OBRMSZ = 'NA'
  set OARMSZ = 'NA'
  set OBUMPS = 'NA'
  set OWBMPS = 'NA'
  set OHBUNS = 'NA'
  set OHBSAT = 'NA'
endif


#Extract statistics from re-refined structure
set NNATOM = `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_besttls.pdb`
if (-e $WORKDIR/wc/pdbout.txt) then
  if (`grep -a -c '1st generation packing quality :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NZPAK1 = 'NA'
  else
    set NZPAK1 = `grep -a '1st generation packing quality :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c '2nd generation packing quality :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NZPAK2 = 'NA'
  else
    set NZPAK2 = `grep -a '2nd generation packing quality :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Ramachandran plot appearance   :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NZRAMA = 'NA'
  else
    set NZRAMA = `grep -a 'Ramachandran plot appearance   :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NCHI12 = 'NA'
  else
    set NCHI12 = `grep -a 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Backbone conformation          :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NBCONF = 'NA'
  else
    set NBCONF = `grep -a 'Backbone conformation          :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond lengths                   :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NBRMSZ = 'NA'
  else
    set NBRMSZ = `grep -a 'Bond lengths                   :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond angles                    :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NARMSZ = 'NA'
  else
    set NARMSZ = `grep -a 'Bond angles                    :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  #Get bump scores
  if (`grep -a -c 'Total number of bumps:' $WORKDIR/wc/pdbout.txt` == 0) then
    set NBUMPS = 0
    set NSBMPL = 0
    set NWBMPS = 0.000
  else
    set NBUMPS = `grep -a 'Total number of bumps:' $WORKDIR/wc/pdbout.txt | cut -c 24-28`
    set NSBMPL = `grep -a 'Total squared bump value:' $WORKDIR/wc/pdbout.txt | cut -c 27-33`
    set NWBMPS = `echo $NNATOM $NSBMPL | awk '{printf "%5.3f\n", 100*$2/$1}'`
  endif

  #Count unsatisfied donors and acceptors
  if (`grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set NHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set NHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  if (`grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set NHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set NHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  @ NHBUNS = ($NHBDON + $NHBACC)
  if (`grep -a -c 'Buried donors:' $WORKDIR/wc/pdbout.txt` == 0) then
    set NHBSAT = 'NA'
  else
    @ NBHBDA = (`grep -a 'Buried donors:' $WORKDIR/wc/pdbout.txt | awk '{print $3}'` + `grep -a 'Buried acceptors:' $WORKDIR/wc/pdbout.txt | awk '{print $3}'`)
    if ($NBHBDA == 0) then
      set NHBSAT = 'NA'
    else
      set NHBSA1 = `grep -a 'with a H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $5} END {print SUM}'`
      set NHBSA2 = `grep -a 'with a poor H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $6} END {print 0.5*SUM}'`
      set NHBSA3 = `grep -a 'with only a very poor H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $8} END {print 0.25*SUM}'`
      set NHBSA4 = `grep -a 'essentially without H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $5} END {print 0.125*SUM}'`
      set NHBSAT = `echo $NHBSA1 $NHBSA2 $NHBSA3 $NHBSA4 $NBHBDA | awk '{printf "%5.3f\n" ,($1 + $2 + $3 + $4)/$5}'`
    endif  
  endif
else
  set NZPAK1 = 'NA'
  set NZPAK2 = 'NA'
  set NZRAMA = 'NA'
  set NCHI12 = 'NA'
  set NBCONF = 'NA'
  set NBRMSZ = 'NA'
  set NARMSZ = 'NA'
  set NBUMPS = 'NA'
  set NWBMPS = 'NA'
  set NHBUNS = 'NA'
  set NHBSAT = 'NA'
endif

#Get statistics from tortoize
if (-e $WORKDIR/${PDBID}_0cyc_tortoize.json) then 
  set OZRAMA  =  `jq '."model"."1"."ramachandran-z"' $WORKDIR/${PDBID}_0cyc_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set OSZRAMA =  `jq '."model"."1"."ramachandran-jackknife-sd"' $WORKDIR/${PDBID}_0cyc_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set OCHI12  =  `jq '."model"."1"."torsion-z"' $WORKDIR/${PDBID}_0cyc_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set OSCHI12 =  `jq '."model"."1"."torsion-jackknife-sd"' $WORKDIR/${PDBID}_0cyc_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
else
  #Scream bloody murder
  echo "COMMENT: tortoize cannot validate 0cyc model" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                              >> $WHYNOT
  exit(1)
endif

#Get statistics from tortoize
if (-e $WORKDIR/${PDBID}_besttls_tortoize.json) then 
  set NZRAMA  =  `jq '."model"."1"."ramachandran-z"' $WORKDIR/${PDBID}_besttls_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set NSZRAMA =  `jq '."model"."1"."ramachandran-jackknife-sd"' $WORKDIR/${PDBID}_besttls_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set NCHI12  =  `jq '."model"."1"."torsion-z"' $WORKDIR/${PDBID}_besttls_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set NSCHI12 =  `jq '."model"."1"."torsion-jackknife-sd"' $WORKDIR/${PDBID}_besttls_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
else
  #Scream bloody murder
  echo "COMMENT: tortoize cannot validate besttls model" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                                 >> $WHYNOT
  exit(1)
endif

#Run distel if needed
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  #Consolidate the distance restraints
  if (-e $WORKDIR/homology.rest) then
    cat $WORKDIR/homology.rest > $WORKDIR/allhb.rest
  endif
  if (-e $WORKDIR/hbond.rest) then
    cat $WORKDIR/hbond.rest >> $WORKDIR/allhb.rest
  endif
  
  #Calculate rmsZ values. PROGRAM: distel.py
  set OHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  set NHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`


  #Clean up
  rm $WORKDIR/allhb.rest
else 
  set OHRMSZ = 'NA'
  set NHRMSZ = 'NA'
endif


#Nucleic acid validation
#Set fallback values
set OBPHBRMSZ = 'NA'
set NBPHBRMSZ = 'NA'
set OSHEAR    = 'NA'
set OSTRETCH  = 'NA'
set OBUCKLE   = 'NA'
set OPROPEL   = 'NA'
set NSHEAR    = 'NA'
set NSTRETCH  = 'NA'
set NBUCKLE   = 'NA'
set NPROPEL   = 'NA'
set ODNRMSD   = 'NA'
set OCONFAL   = 'NA'
set TOCONFAL  = 'NA'
set NDNRMSD   = 'NA'
set NCONFAL   = 'NA'
set TNCONFAL  = 'NA'
set OBPGRMSZ =  'NA'
set NBPGRMSZ =  'NA'

if ($GOT_NUC == 'T') then

  echo "-Validating nucleic acids" | tee -a $LOG

  #Get the basepair values if any were detected
  if (-e $WORKDIR/bphbond.rest) then
    if (`grep -c 'exte dist' $WORKDIR/bphbond.rest` > 0) then
      echo " o Validating base pairs" | tee -a $LOG
      set OBPHBRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/bphbond.rest | tail -n 1 | awk '{print $4}'`
      set NBPHBRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/bphbond.rest | tail -n 1 | awk '{print $4}'`
    
      #Run nucrmsz and get out the results
      python3 $TOOLS/nucrmsz.py $WORKDIR/${PDBID}_0cyc.pdb $TOOLS/x3dna-dssr > $WORKDIR/nucrmsz.log
      set OSHEAR   = `grep 'Shear rmsZ'     $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set OSTRETCH = `grep 'Stretch rmsZ'   $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set OBUCKLE  = `grep 'Buckle rmsZ'    $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set OPROPEL  = `grep 'Propeller rmsZ' $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set OBPGRMSZ = `grep 'bpG rmsZ'       $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
    
      #Run nucrmsz again and get out the results
      python3 $TOOLS/nucrmsz.py $WORKDIR/${PDBID}_besttls.pdb $TOOLS/x3dna-dssr > $WORKDIR/nucrmsz.log
      set NSHEAR   = `grep 'Shear rmsZ'     $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set NSTRETCH = `grep 'Stretch rmsZ'   $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set NBUCKLE  = `grep 'Buckle rmsZ'    $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set NPROPEL  = `grep 'Propeller rmsZ' $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set NBPGRMSZ = `grep 'bpG rmsZ'       $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
    endif
  endif  
  
  #Validate dinucleotides with DNATCO
  echo " o Validating dinucleotides with DNATCO" | tee -a $LOG
  $TOOLS/dnatco.py $WORKDIR/${PDBID}_0cyc.pdb >& $WORKDIR/dnatco.log
  if (-e ${PDBID}_0cyc.pdb.dnatco.json.gz) then 
    set ODNRMSD  = `zcat ${PDBID}_0cyc.pdb.dnatco.json.gz | jq .overall."average_rmsd" | tr -d '"'`
    set OCONFAL  = `zcat ${PDBID}_0cyc.pdb.dnatco.json.gz | jq .overall."confal_score" | tr -d '"'`
    set TOCONFAL = `zcat ${PDBID}_0cyc.pdb.dnatco.json.gz | jq .overall."confal_percentile" | tr -d '"'`
  else
    #Check what went wrong in DNATCO
    if (`grep -c "doesn't contain enough DNA/RNA steps" $WORKDIR/${PDBID}_0cyc.pdb.html` > 0) then
      echo "   * Not enough dinucleotides to use DNATCO" | tee -a $LOG
    else
      echo "   * Unknown DNATCO error" | tee -a $LOG
      echo "COMMENT: DNATCO cannot validate 0cyc model" >> $DEBUG
      echo "PDB-REDO,$PDBID"                            >> $DEBUG
    endif  
  endif
  
  $TOOLS/dnatco.py $WORKDIR/${PDBID}_besttls.pdb >>& $WORKDIR/dnatco.log
  if (-e ${PDBID}_besttls.pdb.dnatco.json.gz) then 
    set NDNRMSD = `zcat ${PDBID}_besttls.pdb.dnatco.json.gz | jq .overall."average_rmsd" | tr -d '"'`
    set NCONFAL  = `zcat ${PDBID}_besttls.pdb.dnatco.json.gz | jq .overall."confal_score" | tr -d '"'`
    set TNCONFAL = `zcat ${PDBID}_besttls.pdb.dnatco.json.gz | jq .overall."confal_percentile" | tr -d '"'`
  else
    #Check what went wrong in DNATCO
    if (`grep -c "doesn't contain enough DNA/RNA steps" $WORKDIR/${PDBID}_besttls.pdb.html` > 0) then
      echo "   * Not enough dinucleotides to use DNATCO" | tee -a $LOG
    else
      echo "   * Unknown DNATCO error" | tee -a $LOG
      echo "COMMENT: DNATCO cannot validate besttls model" >> $DEBUG
      echo "PDB-REDO,$PDBID"                               >> $DEBUG
    endif  
  endif
endif

#Print results
#Start reporting
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Validation details ******" | tee -a $LOG
echo '                                     Before  After' | tee -a $LOG
echo "1st generation packing quality     : $OZPAK1 $NZPAK1"  | tee -a $LOG
echo "2nd generation packing quality     : $OZPAK2 $NZPAK2"  | tee -a $LOG
if ($GOT_PROT == 'T') then
  echo "Ramachandran plot Z-score          : $OZRAMA $NZRAMA"  | tee -a $LOG
  echo "Ramachandran plot RMSD             : $OSZRAMA $NSZRAMA" | tee -a $LOG
  echo "chi-1/chi-2 rotamer Z-score        : $OCHI12 $NCHI12"  | tee -a $LOG
  echo "chi-1/chi-2 rotamer RMSD           : $OSCHI12 $NSCHI12" | tee -a $LOG
  echo "Backbone conformation              : $OBCONF $NBCONF"  | tee -a $LOG
endif  
echo " " | tee -a $LOG
echo "Bond length RMS Z-score            : $OBRMSZ $NBRMSZ"  | tee -a $LOG
echo "Bond angle RMS Z-score             : $OARMSZ $NARMSZ"  | tee -a $LOG
echo " " | tee -a $LOG
echo "Total number of bumps              : $OBUMPS $NBUMPS"  | tee -a $LOG
echo "Weighted bump severity score       : $OWBMPS $NWBMPS"  | tee -a $LOG
echo " " | tee -a $LOG
echo "Unsatisfied H-bond donors/acceptors: $OHBUNS $NHBUNS"  | tee -a $LOG
echo "H-bond satisfaction fraction       : $OHBSAT $NHBSAT"  | tee -a $LOG
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  echo "H-bond restraint RMS Z-score       : $OHRMSZ $NHRMSZ" | tee -a $LOG
endif
if ($GOT_NUC == 'T') then
  echo " " | tee -a $LOG
  if (-e $WORKDIR/bphbond.rest) then
    if (`grep -c 'exte dist' $WORKDIR/bphbond.rest` > 0) then
      echo "Basepair H-bond RMS Z-score        : $OBPHBRMSZ $NBPHBRMSZ" | tee -a $LOG
      echo "Basepair shear RMS Z-score         : $OSHEAR $NSHEAR"       | tee -a $LOG
      echo "Basepair stretch RMS Z-score       : $OSTRETCH $NSTRETCH"   | tee -a $LOG
      echo "Basepair buckle RMS Z-score        : $OBUCKLE $NBUCKLE"     | tee -a $LOG
      echo "Basepair propeller RMS Z-score     : $OPROPEL $NPROPEL"     | tee -a $LOG
      echo "Basepair geometry RMS Z-score      : $OBPGRMSZ $NBPGRMSZ"   | tee -a $LOG 
    endif
  endif  
  echo "Dinucleotide CONFAL score          : $OCONFAL  $NCONFAL"  | tee -a $LOG
  echo "Dinucleotide conformation rmsd     : $ODNRMSD $NDNRMSD"     | tee -a $LOG
endif

################################################## Rebuild structure #####################################################

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Structure rebuilding ******" | tee -a $LOG

#Temporary HETATM and LINK workaround
mv $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_besttls.old
sed '/MSE/s/HETATM/ATOM  /g' $WORKDIR/${PDBID}_besttls.old > $WORKDIR/${PDBID}_besttls.pdb
# $TOOLS/stripper -v $SMODE \
# $WORKDIR/${PDBID}_besttls.mse \
# $WORKDIR/${PDBID}_besttls.pdb \
# $TOOLS/pdb_redo.dat \
# >> $WORKDIR/stripper.log

#Set fall-back values
set NWATDEL = 0
set NBBFLIP = 0
set NSCBLT  = 0
set NSCFLIP = 0
set NCHIRFX = 0
set BUILT   = 0  #The model was not rebuilt
set SCBUILT = 0  #No side-chains were rebuilt 
set NLOOPS  = 0

if ($DOREBUILD == 0) then
  echo "-All rebuilding steps are skipped" | tee -a $LOG
  cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_built.pdb
  cp $WORKDIR/${PDBID}_besttls.mtz $WORKDIR/${PDBID}_loopwhole.mtz
else
  #Calculate BBUILD
  set BBUILD = `echo $BSET | awk '{printf ("%.2f\n", 7*log($1))}'`
  
  #Manage residue renumbering
  if (-e $WORKDIR/renumber.json) then
    set NUMBERING = "-renum $WORKDIR/renumber.json"
  endif
  
  #Append restraint dictionary (again)?
  if (-e $WORKDIR/${PDBID}_het.cif) then
    set DICTCMD = "--dict $WORKDIR/${PDBID}_het.cif"
  endif

  
  #Calculate density fit
  $TOOLS/stats \
  --use-auth-ids \
  --sampling-rate 1.5 \
  --hklin $WORKDIR/${PDBID}_besttls.mtz \
  $SFTYPE \
  --xyzin $WORKDIR/${PDBID}_besttls.pdb \
  -o $WORKDIR/${PDBID}_besttls.json \
  --output-format json \
  $DICTCMD >& $WORKDIR/loopwhole.log 
  
  #Try without dictionary if stats fails
  if (! -e $WORKDIR/${PDBID}_besttls.json && -e $WORKDIR/${PDBID}_het.cif) then
    #Problem probably is in the dictionary
    set DICTCMD = 
    
    #Calculate density fit
    $TOOLS/stats \
    --use-auth-ids \
    --sampling-rate 1.5 \
    --hklin $WORKDIR/${PDBID}_besttls.mtz \
    $SFTYPE \
    --xyzin $WORKDIR/${PDBID}_besttls.pdb \
    -o $WORKDIR/${PDBID}_besttls.json \
    --output-format json \
    $DICTCMD >>& $WORKDIR/loopwhole.log
  endif 
  
  #Stop of no stats file can be made
  if (! -e $WORKDIR/${PDBID}_besttls.json) then
    #Write error statement
    echo " o Cannot calculate density fit. Cannot continue." | tee -a $LOG
    echo "COMMENT: error in stats" >> $WHYNOT
    echo "PDB-REDO,$PDBID"         >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif   
    
    exit(1)
  endif      
        
  #Deal with TLS if needed
  if ($DOTLS == 1) then
    set LTLSCMD  = "-tls $WORKDIR/${PDBID}_refin.tls"
    set LRTLSCMD = `echo TLSIN $WORKDIR/homol/${PDBID}_refin.tls.new`
  else  
    set LTLSCMD  =
    set LRTLSCMD = 
  endif
  
  #Try to complete loops with loopwhole
  if ($DOLOOPS == 1 && $GOT_PROT == T && -e $WORKDIR/$PDBID.blast) then
    echo "-Completing loops with loopwhole" | tee -a $LOG
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.loopwhole.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software."coot-mini-rsr".used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    $TOOLS/loopwhole -v \
    -pdb $WORKDIR/${PDBID}_besttls.pdb \
    -mtz $WORKDIR/${PDBID}_besttls.mtz \
    -eds $WORKDIR/${PDBID}_besttls.json \
    -blast $WORKDIR/$PDBID.blast \
    -output-dir $WORKDIR/homol \
    -output-name $PDBID \
    -fasta $WORKDIR/$PDBID.fasta \
    -tools $TOOLS \
    -pdbdir $REDODIR \
    $LTLSCMD \
    $HOMINCMD \
    $NUMBERING >> $WORKDIR/loopwhole.log 
    
    #Count and report the new loops
    set NLOOPS  = `grep 'Number of loops added' $WORKDIR/loopwhole.log | awk '{print $5}'`
    if ($NLOOPS == '') then
      NLOOPS = 0
    endif
#    set NPLOOPS = `grep -c 'Loop is partly kept' $WORKDIR/loopwhole.log`
    
    #Update files if loops were added
    if ($NLOOPS > 0) then
    
      echo " o Found $NLOOPS potential loops" | tee -a $LOG
      echo " o Validating potential loops" | tee -a $LOG
      echo "   * Updating atomic B-factors and density map" | tee -a $LOG
      
      #Run Refmac 5 cycles to update B-factors
      refmac5 \
      XYZIN  $WORKDIR/homol/${PDBID}_loopwhole.pdb \
      XYZOUT $WORKDIR/homol/${PDBID}_loopwhole_ref.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/homol/${PDBID}_loopwhole_ref.mtz \
      $LIBLIN \
      $LRTLSCMD \
      $SCATLIN \
      <<eof >> $WORKDIR/loopwhole.log
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        make newligand continue
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        tlsd waters exclude
        tlsout addu
        ncyc 5
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG AUTO 2.50
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSSTRICT
        temp $BBEST
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy expdta
        pdbout copy remarks 200 280 350
        NOHARVEST
        END
eof
      
      #Run stats
      echo "   * Calculating density fit" | tee -a $LOG
      $TOOLS/stats \
      --use-auth-ids \
      --sampling-rate 1.5 \
      --hklin $WORKDIR/homol/${PDBID}_loopwhole_ref.mtz \
      $SFTYPE \
      --xyzin $WORKDIR/homol/${PDBID}_loopwhole_ref.pdb \
      -o $WORKDIR/${PDBID}_loopwhole_validate.json \
      --output-format json \
      $DICTCMD >>& $WORKDIR/density_loopwhole.log 
      
      #Filter out bad loops
      cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software."loopwhole-validate".used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
      
      $TOOLS/loopwhole-validate -v \
      -pdb $WORKDIR/homol/${PDBID}_loopwhole.pdb \
      -pdbold $WORKDIR/${PDBID}_besttls.pdb \
      -mtz $WORKDIR/homol/${PDBID}_loopwhole_ref.mtz \
      -output-dir $WORKDIR/homol \
      -output-name $PDBID \
      -tools $TOOLS \
      -eds $WORKDIR/${PDBID}_loopwhole_validate.json >> $WORKDIR/loopwhole.log
      
      @ NLOOPS = ($NLOOPS - `grep 'Number of loops deleted:' $WORKDIR/loopwhole.log | awk '{print $5}'`)
      
      #Report and update maps
      echo "   * $NLOOPS loops accepted" | tee -a $LOG
      echo " o Updating maps from new model" | tee -a $LOG
      
      #Run Refmac 0cyc to update the maps
      refmac5 \
      XYZIN  $WORKDIR/homol/${PDBID}_loopwhole-validate.pdb \
      XYZOUT $WORKDIR/homol/${PDBID}_loopwhole_0cyc.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_loopwhole.mtz \
      $LIBLIN \
      $LRTLSCMD \
      $SCATLIN \
      <<eof >> $WORKDIR/${PDBID}_loopwhole_0cyc.log
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        make newligand continue
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        tlsd waters exclude
        ncyc 0
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG MATRIX 0.5
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSSTRICT
        temp $BBEST
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy expdta
        pdbout copy remarks 200 280 350
        NOHARVEST
        END
eof

      
      #Copy over files
      cp $WORKDIR/homol/${PDBID}_loopwhole-validate.pdb $WORKDIR/${PDBID}_loopwhole.pdb
      if (-e $WORKDIR/homol/${PDBID}.fasta.new) then
        cp $WORKDIR/homol/${PDBID}.fasta.new $WORKDIR/${PDBID}.fasta
      endif  
      if (-e $WORKDIR/homol/${PDBID}_refin.tls.new) then
        cp $WORKDIR/homol/${PDBID}_refin.tls.new $WORKDIR/${PDBID}_refin.tls
      endif  
      
    else  
      #Report and just copy the files
      echo " o No new loops were found" | tee -a $LOG
      cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_loopwhole.pdb
      cp $WORKDIR/${PDBID}_besttls.mtz $WORKDIR/${PDBID}_loopwhole.mtz
    endif
    
    #Check for updated renumbering file
    if (-e $WORKDIR/renumber.json.new) then
      cp $WORKDIR/renumber.json $WORKDIR/renumber.json.old
      mv $WORKDIR/renumber.json.new $WORKDIR/renumber.json
    else if (-e $WORKDIR/homol/${PDBID}_renum.json) then
      cp $WORKDIR/homol/${PDBID}_renum.json $WORKDIR/renumber.json
    endif
  else  
    #Just copy the pdb file
    cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_loopwhole.pdb
    cp $WORKDIR/${PDBID}_besttls.mtz $WORKDIR/${PDBID}_loopwhole.mtz
  endif

  #Run dssp to find secondary structure elements
  if ($GOT_PROT == T) then
    echo "-Assigning secondary structure with DSSP" | tee -a $LOG
    mkdssp $WORKDIR/${PDBID}_loopwhole.pdb $WORKDIR/$PDBID.dssp >& $WORKDIR/dssp.log
    if($status) then
      echo " o DSSP failed" | tee -a $LOG
      echo "   * Not using secondary structure in rebuilding" | tee -a $LOG
      echo "COMMENT: DSSP: general error (1)" >> $DEBUG
      echo "PDB-REDO,$PDBID"                  >> $DEBUG
    endif

    #Is there usable output?
    if (-e $WORKDIR/$PDBID.dssp) then
      set DSSPFILE = $WORKDIR/$PDBID.dssp
    else
      set DSSPFILE = ""
    endif
  endif

  #Check the model completeness and create an omit map of sorts if needed
  if (`echo $COMPLETED | awk '{if ($1 < 50.0) {print "1"} else {print "0"}}'` == 1) then
    #Make a map with SER/PRO/VAL/THR side chains trimmed
    grep -v '^A[NT][IO][SM].........[OC][GD][ 12].[SVTP][EAHR][RLO]' $WORKDIR/${PDBID}_loopwhole.pdb > $WORKDIR/${PDBID}_besttls_trim.pdb

    #Run Refmac 0cyc to make the maps
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_besttls_trim.pdb \
    XYZOUT $WORKDIR/${PDBID}_scomit.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_scomit.mtz \
    $LIBLIN \
    $TLSFILS \
    $SCATLIN \
    <<eof > $WORKDIR/${PDBID}_omit.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      make newligand continue  
      refi type REST resi MLKF meth CGMAT bref $BREFTYPE
      $REFIRES
      tlsd waters exclude
      ncyc 0
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSSTRICT
      temp $BBEST
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy expdta
      pdbout copy remarks 200 280 350
      NOHARVEST
      END
eof
    #Set map coefficients for building
    set SCBUILDMTZ = $WORKDIR/${PDBID}_scomit.mtz
    set BUILDF     = 'FWT'
    set BUILDP     = 'PHWT'
  else
    #Set map coefficients for building
    set SCBUILDMTZ = $WORKDIR/${PDBID}_besttls.mtz
    set BUILDF     = 'FWT'
    set BUILDP     = 'PHWT'

  endif


  #Count the number of waters
  set NWATER = `grep -c -E '^[AH][TE][OT][MA].{13}HOH' $WORKDIR/${PDBID}_loopwhole.pdb`

  if ($DOCENTRIFUGE == 1 && $NWATER > 0) then
    #Delete poor waters
    echo "-Removing waters without density" | tee -a $LOG
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.centrifuge.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json

    $TOOLS/centrifuge \
    <<eof >& $WORKDIR/${PDBID}_centrifuge.log
      #Settings
      InputPDB = "$WORKDIR/${PDBID}_loopwhole.pdb"
      PDBOutputFilename = "$WORKDIR/${PDBID}_centrifuge.pdb"
      MessageFilename = "$WORKDIR/${PDBID}_centrifuge.msg"
      XMLOutputFilename = "$WORKDIR/${PDBID}_centrifuge.xml"

      #residue names that represent waters. Default HOH WAT H2O EAU
      WaterNames = HOH
      #Recognised waters, and dummies, with a density fit lower than this threshold will be removed
      WaterRejectionThreshold = 0.30
      #If true, everything except Dummies and recognised Waters is placed in the density map, default False
      PlaceNonSolventInDensity = False
      #Colon separated list of water ids which won't be affected by the methods of centrifuge
      UnaffectedWaterList = "$H2O_KEEP"
      #InputMap or here: input mtz
      InputMTZ = "$WORKDIR/${PDBID}_loopwhole.mtz"
      FWTLabel = $BUILDF
      PHIWTLabel = $BUILDP
      SpaceGroup = "$SPACEGROUP"
      XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"
      DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
      SymmetryFilename = "$CLIBD/syminfo.lib"

      #type of method used to determine the density fit
      CFitTarget::Type = accelerated
      UniformAtomBFactor = "$BBUILD"
      UniformAtomRadius = 0.74
      UseAtomicB = True

      #Factor with which to invert the density of placed atoms
      AntiBumpFactor = 2.0

      #Program information
      ProgramName = centrifuge
      MessageLevel = 6
      AbortLevel = 8

      KeepPDBheaderInfo = true
      StoreOriginalChainAndSegID = true
      SelectByOrigChainID = true
      KeepSideChain = true
      KeepFragmentLongerThan = 0
      TrustWaters = false

      RemoveBasedOnDensity = true
eof
      if($status) then
        #Write general debug statement
        echo " o Problem with centrifuge" | tee -a $LOG
        echo "COMMENT: error in centrifuge" >> $DEBUG
        echo "PDB-REDO,$PDBID"              >> $DEBUG
        cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_centrifuge.pdb
      else
        #Write debug messsages for building tools
        if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_centrifuge.msg` > 0) then
          grep -H WRNG $WORKDIR/${PDBID}_centrifuge.msg >> $DEBUGB
        endif
      endif

    set NWATDEL = `grep -c 'removing atom' $WORKDIR/${PDBID}_centrifuge.log`

    #If waters were deleted, consider model rebuilt
    if ($NWATDEL > 0) then
      set BUILT = 1
      #Delete external restraints refering to the waters
      if ($RESTCMD != "") then
        foreach WAT (`grep 'removing atom' $WORKDIR/${PDBID}_centrifuge.log | cut -c 16-20`)
          cp $WORKDIR/external.rest $WORKDIR/external.bak
          set CHID = `echo $WAT | cut -c 1-1`
          set RESN = `echo $WAT | cut -c 2-`
          grep -v "chain $CHID resi $RESN" $WORKDIR/external.bak > $WORKDIR/external.rest
        end
      endif
    endif

  else
    echo "-Skipping centrifuge run" | tee -a $LOG
    cp $WORKDIR/${PDBID}_loopwhole.pdb $WORKDIR/${PDBID}_centrifuge.pdb
  endif

  #Further rebuilding steps

  #Perform peptide flips if resolution is adequate and there is protein
  if ($DOPEPFLIP == 1 && $GOT_PROT == T && `echo $URESO | awk '{if ($1 < 3.30) {print "1"}}'` == 1) then
    echo "-Performing peptide flips" | tee -a $LOG

    #Run pepflip
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.pepflip.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software."coot-mini-rsr".used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    $TOOLS/pepflip \
<<eof >& $WORKDIR/${PDBID}_pepflip.log
      ProgramName = pepflip
      MessageFilename = $WORKDIR/${PDBID}_pepflip.msg
      MessageLevel = 6
      XMLOutputFilename = $WORKDIR/${PDBID}_pepflip.xml
      AbortLevel = 8

      #input pdb, and its handling
      InputPDB = $WORKDIR/${PDBID}_centrifuge.pdb
      SelectByOrigChainID =true
      KeepPDBheaderInfo =true

      #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
      InputMTZ = $WORKDIR/${PDBID}_loopwhole.mtz
      FWTLabel = $BUILDF
      PHIWTLabel = $BUILDP
      #Use this masked map
      OutputMaskedNonAAmap = $WORKDIR/masked.map

      #Output map with the pepflip version of asu limits, needed for loopfit
      #OutputExtendedMap = $WORKDIR/${PDBID}_map.ext
      #density settings for the computations on the difference density at the oxygen
      DifDensity::InputMTZ = $WORKDIR/${PDBID}_loopwhole.mtz
      DifDensity::FWTLabel = DELFWT
      DifDensity::PHIWTLabel = PHDELWT
      DifDensity::CFitTarget::Type = accelerated
      DifDensity::CFitTarget::ConvertToZscore = false

      UniformAtomRadius= 0.74
      AntiBumpFactor= 2.0
      SpaceGroup = "$SPACEGROUP"
      XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"
      UniformAtomBFactor = "$BBUILD"

      FlipMode = repair
      #Output pdb with the flipped peptides
      FlippedOutputPDB = $WORKDIR/${PDBID}_pepflip.pdb

      #Exclude lists
      DontFlipListN = "$BBN_KEEP"
      DontFlipListO = "$BBO_KEEP"

      #Optional filename (output from DSSP) to set the secondary structure of the complex
      DSSPfilename = $DSSPFILE
      #Number of residues (default 2) at the ends of secondary structure, which won't be trusted
      BufferSS=1
      #If true (default false), and the secondary structure is set, trust the middle of helices (H)
      #and don't consider those residues for flips
      TrustHelices = true
      Trust310Helices = true
      TrustPiHelices = true
      TrustBetaStrands = true

      #Refinement parameters
      Loopfit::UseDefaultLoopfit =false
      Loopfit::UseMiniRSR=true
      Loopfit::InputMap= $WORKDIR/masked.map
      Loopfit::LoopfitExe = coot-mini-rsr
      Loopfit::LoopfitLog = $WORKDIR/pepflip.fit
      Loopfit::MiniRSRtorsions = true
      Loopfit::MiniRSRrama = true
      Loopfit::MiniRSRweight = 10.00
      #temporary pdb file for loopfit/mini-rsr
      Loopfit::TmpOutputPDB = $WORKDIR/${PDBID}_fit.pdb
      #Location of the CCP4 dictionary
      Loopfit::MonomerLib = $CLIBD_MON
      WeightGeoGooF = 24

      #local dictionary files
      DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
      SymmetryFilename = "$CLIBD/syminfo.lib"

      #Methods for density computations
      CFitTarget::Type = cubic
      CFitTarget::MaskRadius = 2.0
      CFitTarget::UseObservedDensityAroundAtoms = "CA CB"

      #In repair mode - If true, place the part of the flanking residues that doesn't belong to the peptide in the map,
      #default false
      RepairPlaceFlanking = true
      #If true, place known residues in the map, before checking the density fit of the peptide - This will slow down
      #the progress significantly, default false
      CheckDensityPlaceKnown = false

      #Reject flip if the new density < average(density) + 'MinSigmaLevel'*sigma(density)
      MinSigmaLevel = -3.5
      #Ratio between the density fit of the original peptide and the flipped version, to be considered as possible flip (default 0.9)
      MinRatioOrigFlip =0.9
      #Density fit of Oxygen in the difference map, indicating a clear positive monopole (default 2.)
      DensityPosOmonopole = 2.4
      #Density fit of Oxygen in the difference map, indicating a clear negative monopole (default -3.)
      DensityNegOmonopole = -2.4
eof
    if($status) then
      echo " o Problem with pepflip. Cannot continue." | tee -a $LOG
      echo "COMMENT: error in pepflip" >> $WHYNOT
      echo "PDB-REDO,$PDBID"           >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      #Write debug messsages for building tools
      if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_pepflip.msg` > 0) then
        grep -H WRNG $WORKDIR/${PDBID}_pepflip.msg >> $DEBUGB
      endif
    endif

    #Count the number of flips
    set NBBFLIP = `grep -A 80 'List of peptides' $WORKDIR/${PDBID}_pepflip.log | grep -c -E '^.[0-9]'`

    #If peptides were flipped, consider model rebuilt
    if ($NBBFLIP > 0) then
      set BUILT = 1
    endif

  else
    echo "-Skipping pepflip run"
    cp $WORKDIR/${PDBID}_centrifuge.pdb $WORKDIR/${PDBID}_pepflip.pdb
  endif

  #Rebuild side chains if resolution is adequate and there is protein
  if ($DOSCBUILD == 1 && $GOT_PROT == T && `echo $URESO | awk '{if ($1 < 3.30) {print "1"}}'` == 1) then
    echo "-Refitting side chains" | tee -a $LOG

    #Run dssp again to find secondary structure elements if there were any peptide flips
    if ($NBBFLIP > 0) then
      mkdssp $WORKDIR/${PDBID}_centrifuge.pdb $WORKDIR/$PDBID.dssp >>& $WORKDIR/dssp.log
      if($status) then
        echo " o Problem with DSSP" | tee -a $LOG
        echo "   * Not using secondary structure in SideAide" | tee -a $LOG
        echo "COMMENT: DSSP: general error (2)" >> $DEBUG
        echo "PDB-REDO,$PDBID"                  >> $DEBUG
      endif

      #Is there usable output?
      if (-e $WORKDIR/$PDBID.dssp) then
        set DSSPFILE = $WORKDIR/$PDBID.dssp
      else
        set DSSPFILE = ""
      endif
    endif

    #Run SideAide
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.SideAide.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    $TOOLS/SideAide \
  <<eof >& $WORKDIR/${PDBID}_scbuild.log
      ProgramName = SideAide
      MessageFilename = $WORKDIR/${PDBID}_scbuild.msg
      MessageLevel = 6
      XMLOutputFilename = $WORKDIR/${PDBID}_scbuild.xml
      AbortLevel = 8

      #Input pdb, and its handling
      #IncludeChains = '-'
      InputPDB = "$WORKDIR/${PDBID}_pepflip.pdb"

      #list of residue names representing waters, default HOH H2O WAT EAU
      #WaterNames = HOH WAT H2O EAU

      #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
      InputMTZ = "$SCBUILDMTZ"
      FWTLabel = $BUILDF
      PHIWTLabel = $BUILDP
      SpaceGroup = "$SPACEGROUP"
      XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"

      #B factor used in map calculation (for every atom in every residue):
      UniformAtomBFactor = $BBUILD
      UniformAtomRadius = 0.74
      AntiBumpFactor = 1.5

      #Optional filename (output from DSSP) to set the secondary structure of the complex
      DSSPfilename = $DSSPFILE

      #Output pdb
      PDBOutputFilename = "$WORKDIR/${PDBID}_scbuild.pdb"

      #selection of the side chains to rebuild/flip
      #list of single residues:
      RebuildAll = true #If false, no residues are rebuilt (use for flipping without rebuilding)
      UseAvBforNewSCatoms = true #default
      #BFactorForNewSCatoms = 15 #ouse if UseAvBforNewSCatoms = false

      #RebuildList = ":A33 :A66 :A69 :A87 :A167 :A245 :A275 :A304 :A314 :A412 :"
      DontRebuildList = "$NO_REBUILD"
      TrustedWaterList = "$H2O_KEEP"
      #list of residue regions
      #RebuildDefinition = "A36(7)A42:AAAA90(6)AAAA95"
      #list of residues to flip
      #FlipResidueList = ":A23 :"

      #local dictionary files
      DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
      SymmetryFilename = "$CLIBD/syminfo.lib"

      # Methods for density computations
      CFitTarget::Type = accelerated
      CFitTarget::MaskRadius = 2.0
      CFitTarget::UseObservedDensityAroundAtoms = "CA"
      RefineScope::CFitTarget::Type = accelerated
      FitRotamerScope::CFitTarget::Type = correlation

      #Dummies and dubious ignored residues, i.e. possible waters, are removed when their density lies below this value
      DummyRejectionThreshold = -0.25
      RejectionScope::CFitTarget::type = accelerated #kies de fittarget die hoort by DummyRejectionThreshold.

      #Advanced parameter settings
      KeepPDBheaderInfo = true
      KeepUnknownAsIs = true

      #if true, the order of fitting rotamers and refinement is not only based on the AA type, but only on whether the
      # rotamer fit exceeds a threshold (these will be fit and refined in a second round
      OrderOnRotamerValidation = true
      #list of AA type and their validation threshold ('all' can be used to set the threshold of every aa type).
      ValidationThresholds = "all 0.4 ALA 0.3"

      #Below:  default is false, if true, the original chain ID and Seg ID are saved in the output pdb.
      #Only used for CDummyResidues and CFullRes, i.e. ignored residues
      StoreOriginalChainAndSegID = true
      KeepSideChain = true
      KeepFragmentLongerThan = 0
      #TrustWaters = false
      AminoAcidFittingOrder = "GLY ALA PRO VAL THR SER CYS TRP PHE TYR HIS ILE LEU ASP ASN GLU GLN MSE MET LYS ARG"
      #place residues outside IncludeChains in the map:
      InitTargetWithIgnored = true
      ThresholdNewRefinedSC = -0.04 #Only keep sidechains if the fit improves by this value

      #Parameters to shift the CA to try and find better rotamers, optional
      #Values to twist the N-CA-CB angle over, to try and achieve better rotamer fitting;
      #Radius of the CA shell of shifted options
      ShiftingCAshell = 0.2
      #number of CA, or closest fibonaccinumber, in the shell of CA options
      ShiftingCAnumber = 8
      #only try to find a better rotamer if the density score is lower than this value, 10 = NO THRESHOLD
      ShiftingCAthreshold = 10.
      #only keep a rotamer with a shifted CA if the density score is higher than this value
      ShiftingCAkeepThreshold = -10
      #-only keep a rotamer with a shifted CA if the relative improvement in the score is higher than this value
      ShiftingCAimprovementThreshold = 0.05

      #Neither RebuildDefinition nor ElalXML will rebuild all the side chains
      SelectByOrigChainID = true
      Overlap = 0
eof
    if($status) then
      echo " o Problem with SideAide. Cannot continue." | tee -a $LOG
      echo "COMMENT: error in SideAide" >> $WHYNOT
      echo "PDB-REDO,$PDBID"            >> $WHYNOT
      if ($SERVER == 1) then
	#Write out status files
	touch $STDIR/stoppingProcess.txt
	touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      #Write debug messsages for building tools
      if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_scbuild.msg` > 0) then
        grep -H WRNG $WORKDIR/${PDBID}_scbuild.msg >> $DEBUGB
      endif
      
      #Filthy hack for SideAide quirk! Add back waters from metal restraints.
      if (-e $WORKDIR/metal.rest) then
        #Construct a list of potential waters among the restraints
        cat $WORKDIR/metal.rest | awk '{if ($11 == "O") {printf "HOH %s%4d\n", $5, $7}; if ($20 == "O") {printf "HOH %s%4d\n", $14, $16}; if ($29 == "O") {printf "HOH %s%4d\n", $23, $25}}' | sort -u > $WORKDIR/water.list
        #Make list of actual waters
        grep --file=$WORKDIR/water.list $WORKDIR/${PDBID}_centrifuge.pdb | cut -c 18-26 > $WORKDIR/water.have
        #Make list of waters that must be replaced
        grep --file=$WORKDIR/water.list $WORKDIR/${PDBID}_scbuild.pdb | cut -c 18-26 | comm -3 $WORKDIR/water.have - > $WORKDIR/water.back
        #Replace them
        grep --file=$WORKDIR/water.back $WORKDIR/${PDBID}_centrifuge.pdb >> $WORKDIR/${PDBID}_scbuild.pdb
      endif
    endif

    #Count the number of completed side chains.
    set NSCBLT  = `grep -c 'building the complete side chain' $WORKDIR/${PDBID}_scbuild.log`

    #Some side chains must have changed, so consider the model rebuilt
    set BUILT = 1
    set SCBUILT = 1

  else
    echo "-Skipping side chain rebuilding" | tee -a $LOG
    cp $WORKDIR/${PDBID}_pepflip.pdb $WORKDIR/${PDBID}_scbuild.pdb
  endif

######################################## Flip side chains to optimise hydrogen bonding ###################################

  #Always do this if there is protein
  if ($GOT_PROT == T) then

    echo "-Validation-based rebuilding" | tee -a $LOG

    #Check for candidates using WHAT_CHECK
    #Go to temporary running directory
    mkdir -p $WORKDIR/flip
    cd $WORKDIR/flip

    #Get the PDB file
    cp $WORKDIR/${PDBID}_scbuild.pdb $WORKDIR/flip

    #Do the actual validation
    $WC/bin/whatcheck $WORKDIR/${PDBID}_scbuild.pdb Y Y Y >& $WORKDIR/flip.log

    #Check for an output file
    if (-e $WORKDIR/flip/pdbout.txt) then

      #Create a fliplist
      cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.what_todo.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
      
      $TOOLS/what_todo $WORKDIR/flip/pdbout.txt $WORKDIR/$PDBID.todo $WORKDIR/$PDBID.extracted

      #Do the flipping
      cd $WORKDIR
      set FLIPLIST = `head -n 2  $WORKDIR/$PDBID.todo | tail -n 1`
      echo " o Flipping side chains"

      #Run SideAide
      cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.SideAide.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
      
      $TOOLS/SideAide \
<<eof >& $WORKDIR/${PDBID}_scflip.log
        ProgramName = SideAide
        MessageFilename = $WORKDIR/${PDBID}_scflip.msg
        MessageLevel = 6
        XMLOutputFilename = $WORKDIR/${PDBID}_scflip.xml
        AbortLevel = 8

        #Input pdb, and its handling
        InputPDB = "$WORKDIR/${PDBID}_scbuild.pdb"

        #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
        InputMTZ = "$WORKDIR/${PDBID}_loopwhole.mtz"
        FWTLabel = $BUILDF
        PHIWTLabel = $BUILDP
        SpaceGroup = "$SPACEGROUP"
        XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"

        #B factor used in map calculation (for every atom in every residue):
        UniformAtomBFactor = $BBUILD
        UniformAtomRadius = 0.74
        AntiBumpFactor = 2.0

        #Optional filename (output from DSSP) to set the secondary structure of the complex
        DSSPfilename = $DSSPFILE

        #Output pdb
        PDBOutputFilename = "$WORKDIR/${PDBID}_built.pdb"

        #selection of the side chains to rebuild/flip
        #list of single residues:
        RebuildAll = false #If false, no residues are rebuilt (use for flipping without rebuilding)
        UseAvBforNewSCatoms = true #default
        #BFactorForNewSCatoms = 15 #ouse if UseAvBforNewSCatoms = false
        DontRebuildList = "$NO_REBUILD"
        TrustedWaterList = "$H2O_KEEP"
        #list of residues to flip
        FlipResidueList = "$FLIPLIST"

        #local dictionary files
        DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
        SymmetryFilename = "$CLIBD/syminfo.lib"

        #Methods for density computations
        CFitTarget::Type = accelerated
        CFitTarget::MaskRadius = 2.0
        CFitTarget::UseObservedDensityAroundAtoms = "CA"
        RefineScope::CFitTarget::Type = accelerated
        FitRotamerScope::CFitTarget::Type = correlation

        #Dummies and dubious ignored residues, i.e. possible waters, are removed when their density lies below this value
        DummyRejectionThreshold = -0.1
        RejectionScope::CFitTarget::type = accelerated #kies de fittarget die hoort by DummyRejectionThreshold.

        #Advanced parameter settings
        KeepPDBheaderInfo = true
        KeepUnknownAsIs = true

        #if true, the order of fitting rotamers and refinement is not only based on the AA type, but only on whether the
        # rotamer fit exceeds a threshold (these will be fit and refined in a second round
        OrderOnRotamerValidation = true
        #list of AA type and their validation threshold ('all' can be used to set the threshold of every aa type).
        ValidationThresholds = "all 0.4 ALA 0.3"

        #Below:  default is false, if true, the original chain ID and Seg ID are saved in the output pdb.
        #Only used for CDummyResidues and CFullRes, i.e. ignored residues
        StoreOriginalChainAndSegID = true
        KeepSideChain = true
        KeepFragmentLongerThan = 0
        #TrustWaters = false
        AminoAcidFittingOrder = "GLY ALA PRO VAL THR SER CYS TRP PHE TYR HIS ILE LEU ASP ASN GLU GLN MET LYS ARG"
        #place residues outside IncludeChains in the map:
        InitTargetWithIgnored = true
        ThresholdNewRefinedSC = -0.04 #Only keep sidechains if the fit improves by this value

        #Parameters to shift the CA to try and find better rotamers, optional
        #Values to twist the N-CA-CB angle over, to try and achieve better rotamer fitting;
        #Radius of the CA shell of shifted options
        ShiftingCAshell = 0.2
        #number of CA, or closest fibonaccinumber, in the shell of CA options
        ShiftingCAnumber = 8
        #only try to find a better rotamer if the density score is lower than this value, 10 = NO THRESHOLD
        ShiftingCAthreshold = 10.
        #only keep a rotamer with a shifted CA if the density score is higher than this value
        ShiftingCAkeepThreshold = -10
        #only keep a rotamer with a shifted CA if the relative improvement in the score is higher than this value
        ShiftingCAimprovementThreshold = 0.10

        #Neither RebuildDefinition nor ElalXML will rebuild all the side chains
        SelectByOrigChainID = true
        Overlap = 0
eof
      if($status) then
        echo "   * Problem with SideAide. Cannot continue." | tee -a $LOG
        echo "COMMENT: error in side chain flipping" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                       >> $WHYNOT
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      else
        #Write debug messsages for building tools
        if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_scflip.msg` > 0) then
          grep -H WRNG $WORKDIR/${PDBID}_scflip.msg >> $DEBUGB
        endif
      endif

      #Get the number of flipped residues
      set NSCFLIP = `grep -c 'Flipping residue' $WORKDIR/${PDBID}_scflip.log`

      #If sidechains were flipped, consider model rebuilt
      if ($NSCFLIP > 0) then
        set BUILT = 1
      endif

      #Force rebuilding seriously problematic side chains
      if ( (`grep -c '#Residues to be rebuilt forcefully' $WORKDIR/$PDBID.todo` != 0) || ($SCBUILT == 0 && -e $WORKDIR/homol/${PDBID}_sc_rebuild.txt) ) then
        set CHIRLIST = `grep -A 1 '#Residues to be rebuilt forcefully' $WORKDIR/$PDBID.todo | tail -n 1`
        cp $WORKDIR/${PDBID}_built.pdb $WORKDIR/${PDBID}_chirfix_in.pdb

        #Add notify the user
        if (`grep -c '#Residues to be rebuilt forcefully' $WORKDIR/$PDBID.todo` != 0) then
          echo " o Chirality fixes are needed" | tee -a $LOG
        endif
        if ($SCBUILT == 0 && -e $WORKDIR/homol/${PDBID}_sc_rebuild.txt)  then
          echo " o Completing residues trucated by loopwhole" | tee -a $LOG
          set CHIRLIST2 = `cat $WORKDIR/homol/${PDBID}_sc_rebuild.txt`
          set CHIRLIST  = `echo ${CHIRLIST}":"${CHIRLIST2}`
          
        endif

        #Run SideAide
        cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.SideAide.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
        
        $TOOLS/SideAide \
<<eof >& $WORKDIR/${PDBID}_chirfix.log
        ProgramName = SideAide
        MessageFilename = $WORKDIR/${PDBID}_chirfix.msg
        MessageLevel = 6
        XMLOutputFilename = $WORKDIR/${PDBID}_chirfix.xml
        AbortLevel = 8

        #Input pdb, and its handling
        #IncludeChains = '-'
        InputPDB = "$WORKDIR/${PDBID}_chirfix_in.pdb"

        #list of residue names representing waters, default HOH H2O WAT EAU
        #WaterNames = HOH WAT H2O EAU

        #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
        InputMTZ = "$WORKDIR/${PDBID}_loopwhole.mtz"
        FWTLabel = $BUILDF
        PHIWTLabel = $BUILDP
        SpaceGroup = "$SPACEGROUP"
        XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"

        #B factor used in map calculation (for every atom in every residue):
        UniformAtomBFactor = $BBUILD
        UniformAtomRadius = 0.74
        AntiBumpFactor = 2.0

        #Optional filename (output from DSSP) to set the secondary structure of the complex
        DSSPfilename = $DSSPFILE

        #Output pdb
        PDBOutputFilename = "$WORKDIR/${PDBID}_built.pdb"

        #selection of the side chains to rebuild/flip
        #list of single residues:
        RebuildAll = false #If false, no residues are rebuilt (use for flipping without rebuilding)
        UseAvBforNewSCatoms = true #default
        #BFactorForNewSCatoms = 15 #ouse if UseAvBforNewSCatoms = false

        RebuildList = "$CHIRLIST"
        TrustedWaterList = "$H2O_KEEP"
        #list of residue regions
        #list of residues to flip

        #local dictionary files
        DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
        SymmetryFilename = "$CLIBD/syminfo.lib"

        #Methods for density computations
        CFitTarget::Type = accelerated
        CFitTarget::MaskRadius = 2.0
        CFitTarget::UseObservedDensityAroundAtoms = "CA"
        RefineScope::CFitTarget::Type = accelerated
        FitRotamerScope::CFitTarget::Type = correlation

        #Dummies and dubious ignored residues, i.e. possible waters, are removed when their density lies below this value
        DummyRejectionThreshold = -0.1
        RejectionScope::CFitTarget::type = accelerated #kies de fittarget die hoort by DummyRejectionThreshold.

        #if true, the order of fitting rotamers and refinement is not only based on the AA type, but only on whether the
        # rotamer fit exceeds a threshold (these will be fit and refined in a second round
        OrderOnRotamerValidation = true
        #list of AA type and their validation threshold ('all' can be used to set the threshold of every aa type).
        ValidationThresholds = "all 0.4 ALA 0.3"

        #Advanced parameter settings
        KeepPDBheaderInfo = true
        KeepUnknownAsIs = true

        #Below:  default is false, if true, the original chain ID and Seg ID are saved in the output pdb.
        #Only used for CDummyResidues and CFullRes, i.e. ignored residues
        StoreOriginalChainAndSegID = true
        KeepSideChain = true
        KeepFragmentLongerThan = 0
        #TrustWaters = false
        AminoAcidFittingOrder = "GLY ALA PRO VAL THR SER CYS TRP PHE TYR HIS ILE LEU ASP ASN GLU GLN MET LYS ARG"
        #place residues outside IncludeChains in the map:
        InitTargetWithIgnored = true
        ThresholdNewRefinedSC = -1.00 #Very low value to ensure that original side chains are not kept

        #Parameters to shift the CA to try and find better rotamers, optional
        #Values to twist the N-CA-CB angle over, to try and achieve better rotamer fitting;
        #Radius of the CA shell of shifted options
        ShiftingCAshell = 0.2
        #number of CA, or closest fibonaccinumber, in the shell of CA options
        ShiftingCAnumber = 8
        #only try to find a better rotamer if the density score is lower than this value, 10 = NO THRESHOLD
        ShiftingCAthreshold = 10.
        #only keep a rotamer with a shifted CA if the density score is higher than this value
        ShiftingCAkeepThreshold = -10
        #-only keep a rotamer with a shifted CA if the relative improvement in the score is higher than this value
        ShiftingCAimprovementThreshold = 0.05

        #Neither RebuildDefinition nor ElalXML will rebuild all the side chains
        SelectByOrigChainID = true
        Overlap = 0
eof
        if($status) then
          echo "   * Problem with SideAide. Cannot continue." | tee -a $LOG
          echo "COMMENT: error in chirality fixing" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                    >> $WHYNOT
          if ($SERVER == 1) then
            #Write out status files
            touch $STDIR/stoppingProcess.txt
            touch $STDIR/processStopped.txt
          endif
          cd $BASE
          exit(1)
        else
          #Write debug messsages for building tools
          if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_chirfix.msg` > 0) then
            grep -H WRNG $WORKDIR/${PDBID}_chirfix.msg >> $DEBUGB
          endif
        endif

        #Reset the number of chirality fixes
        set NCHIRFX = `echo $CHIRLIST | awk -F ':' 'BEGIN {SUM = 0} {for (i=1; i<=NF; i++) {if ($i != "") {SUM = SUM + 1}}} END {print SUM}'`

        #If chiralities were fixed, consider model rebuilt
        if ($NCHIRFX > 0) then
          set BUILT = 1
        endif
      endif

      #Clean up
      rm -rf $WORKDIR/flip

    else
      #Give warning
      echo " o Validation-based rebuilding failed" | tee -a $LOG
      set WCERR = 1
      echo "COMMENT: WHAT_CHECK failed in rebuilding step." >> $DEBUG
      echo "PDB-REDO,$PDBID"                                >> $DEBUG
      cd $WORKDIR
      cp $WORKDIR/${PDBID}_scbuild.pdb $WORKDIR/${PDBID}_built.pdb
    endif
  else
    #Fill in some of the blanks
    set NSCFLIP = 0
    cp $WORKDIR/${PDBID}_scbuild.pdb $WORKDIR/${PDBID}_built.pdb
  endif

  #Correct the total number of deleted waters
  @ NWATDEL = (`grep '^[AH][TE][OT]' $WORKDIR/cache.pdb | grep -c HOH` - `grep '^[AH][TE][OT]' $WORKDIR/${PDBID}_built.pdb | grep -c HOH`)
  
  #Present a rebuilding summary.
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Rebuilding details ******" | tee -a $LOG
  echo "Loops added        : $NLOOPS"  | tee -a $LOG
  echo "Waters removed     : $NWATDEL" | tee -a $LOG
  echo "Peptides flipped   : $NBBFLIP" | tee -a $LOG
  echo "Side chains built  : $NSCBLT"  | tee -a $LOG
  echo "Side chains flipped: $NSCFLIP" | tee -a $LOG
  echo "Chirality fixes    : $NCHIRFX" | tee -a $LOG
endif

#Calculate the total number of chirality fixes
@ NCHIRFX = ($CHIFIX + $NCHIRFX)

################################################## (Re)build carbohydrates ###############################################
if ($DOSUGARBUILD == 1 && $GOT_PROT == 'T') then

  #Start reporting
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Carbohydrate (re)building ******" | tee -a $LOG

  #Make a backup
  cp $WORKDIR/${PDBID}_built.pdb $WORKDIR/${PDBID}_carb.pdb
  
  #Run carbivore
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.carbivore.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.COOT.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  echo "-Running carbivore" | tee -a $LOG
  
  $TOOLS/carbivore -v \
  -fasta $WORKDIR/$PDBID.fasta \
  -pdb $WORKDIR/${PDBID}_carb.pdb \
  -mtz $WORKDIR/${PDBID}_loopwhole.mtz \
  -cache-pdb $WORKDIR/cache.pdb \
  -output-name $PDBID \
  -tools $TOOLS \
  -output-dir $WORKDIR/sugar >& $WORKDIR/carbivore.log
  
  #Copy over the output and anaylse results
  if (-e $WORKDIR/sugar/${PDBID}_carbivore.pdb) then
    cp $WORKDIR/sugar/${PDBID}_carbivore.pdb $WORKDIR/${PDBID}_built.pdb
    
    #Analyse results
    echo " o Carbohydrate tree processing results" | tee -a $LOG
    echo "   * Whole-trees built: " `grep 'whole trees built:' $WORKDIR/carbivore.log | awk '{print $6}'` | tee -a $LOG
    echo "   * Trees rebuilt    : " `grep 'trees rebuilt:'     $WORKDIR/carbivore.log | awk '{print $5}'` | tee -a $LOG
    echo "   * Trees extended   : " `grep 'trees extended:'    $WORKDIR/carbivore.log | awk '{print $5}'` | tee -a $LOG
    echo " o Overall results" | tee -a $LOG
    echo "   * Carbohydrate residues added  : " `grep 'carbohydrates built:'   $WORKDIR/carbivore.log | awk '{print $5}'` | tee -a $LOG
    echo "   * Carbohydrate residues deleted: " `grep 'carbohydrates deleted:' $WORKDIR/carbivore.log | awk '{print $5}'` | tee -a $LOG
    echo "   * Overlapping waters deleted   : " `grep 'waters deleted:'        $WORKDIR/carbivore.log | awk '{print $5}'` | tee -a $LOG
    
    #Correct number of deleted waters
    @ NWATDEL = (`grep '^[AH][TE][OT]' $WORKDIR/pdb$PDBID.ent | grep -c HOH` - `grep '^[AH][TE][OT]' $WORKDIR/${PDBID}_built.pdb | grep -c HOH`)
  else
    #Report error
    echo " o Carbivore failed, no carbohydrates (rebuilt)" | tee -a $LOG
    echo "Carbivore: general error" >> $DEBUG
    echo "PDB-REDO,$PDBID"          >> $DEBUG
  endif
endif

################################################### Fit ligands ##########################################################

if ($FITLIGANDS == 1) then
  #Report
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Ligand fitting ******" | tee -a $LOG
  
  #Check whether R-free is good enough to allow ligand fitting 
  if (`echo $RFTLS | awk '{if ($1 > 0.35) {print "1"} else {print "0"}}'` == 1) then
    #Abort ligand fitting
    echo " o R-free too high, model not yet suited for ligand fitting" | tee -a $LOG
  else
    #Try to fit ligands
    set RESNUM = -999
    set CHID   = 'z'
    set NLIG_LIST = 
    set LIGSIGMA = 1.05
  
    cp $WORKDIR/${PDBID}_built.pdb $WORKDIR/${PDBID}_ligfit.pdb

    #Loop over all ligands to prepare for fitting
    foreach LIG (`echo $FITLIGS`)
  
      #Clean a bit if needed
      if (-e $WORKDIR/fitted-ligand-0-0.pdb) then
        rm $WORKDIR/fitted-ligand*.pdb
      endif
  
      #Report
      echo "-Preparing for ligand $LIG" | tee -a $LOG
      set FITCOUNT = 0
      set ACCEPTED = 0
    
      #Find the restraint file (user restraints take precedence)
      set D1 = `echo $LIG | tr "[A-Z]" "[a-z]" | cut -c 1`
      if (-e $WORKDIR/${PDBID}_het.cif) then
        if (`grep -c data_comp_$LIG $WORKDIR/${PDBID}_het.cif` == 1) then
          set LIGREST = $WORKDIR/${PDBID}_het.cif
        else if (-e $CLIBD_MON/$D1/$LIG.cif) then
          set LIGREST = $CLIBD_MON/$D1/$LIG.cif 
        else  
          set LIGREST = 'none'
        endif  
      else if (-e $CLIBD_MON/$D1/$LIG.cif) then
        set LIGREST = $CLIBD_MON/$D1/$LIG.cif 
      else   
        set LIGREST = 'none'
      endif
    
      #Continue if there is a restraint file
      if ($LIGREST != 'none') then
        echo "Making idealised coordinates for $LIG" >> $WORKDIR/ligfit.log
    
        #Make idealised model of ligand. The empty line is needed!
        cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.libcheck.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
        
        libcheck << eof >> $WORKDIR/ligfit.log
          N    
          file_l $LIGREST
          mon $LIG
          nodist y
          coor y
        
eof
      
        #Remove hydrogens from idealised coordinates
        if (! -e $WORKDIR/libcheck_$LIG.pdb) then
          echo " o Could not generate coordinates for $LIG. Skipping ligand." | tee -a $LOG
          continue  
        else  
          grep -v -E '.{76} H' $WORKDIR/libcheck_$LIG.pdb > libcheck_noH_$LIG.pdb
        endif

        #Get liggand sizes in number of heavy atoms
        echo `grep -c '^ATOM' libcheck_noH_$LIG.pdb` $LIG $LIGREST >> $WORKDIR/liglist.txt
      else
        #Give error message
        echo " o There are no restraints available for $LIG. Please, provide a restraint file." | tee -a $LOG
      endif
    end 
     
    #Fit in order of decreasing size
    foreach LIG (`sort -gr $WORKDIR/liglist.txt | awk '{print $2}'`)  
   
      #Report
      echo "-Fitting ligand $LIG" | tee -a $LOG
      set FITCOUNT = 0
      set ACCEPTED = 0
    
      #Get the correct restraint file
      set LIGREST = `grep "$LIG " $WORKDIR/liglist.txt | awk '{print $3}'`
   
      #Do an initial fit of the ligand
      cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.findligand.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
      
      echo "Trying to fit $LIG" >> $WORKDIR/ligfit.log
      findligand \
      --pdbin $WORKDIR/${PDBID}_ligfit.pdb \
      --hklin $WORKDIR/${PDBID}_loopwhole.mtz \
      --f FWT --phi PHWT \
      --flexible \
      --dictionary $LIGREST \
      --samples 100 \
      --sigma $LIGSIGMA \
      $WORKDIR/libcheck_noH_$LIG.pdb >> $WORKDIR/ligfit.log
      
      #Was a ligand added?
      if (`find $WORKDIR -name "fitted-ligand-?-?.pdb" | wc -l` > 0) then
        #Report
        echo "Found `ls $WORKDIR/fitted-ligand*.pdb | wc -l` candidates for $LIG" >> $WORKDIR/ligfit.log
        echo " o Found `ls $WORKDIR/fitted-ligand*.pdb | wc -l` candidates for $LIG" | tee -a $LOG
        echo " o Refinement and validation of candidates" | tee -a $LOG

        #Loop over candidates
        foreach FITTED (`ls $WORKDIR/fitted-ligand*.pdb`)
         
          @ FITCOUNT = ($FITCOUNT + 1)
          #Do some renumbering
          while (`grep -c "z$RESNUM" $WORKDIR/${PDBID}_ligfit.pdb` > 0)
            @ RESNUM = ($RESNUM + 1)
          end
          grep '^HETATM' $FITTED | sed "s/A   1/z$RESNUM/g" > $WORKDIR/append.pdb
         
          #Append the ligand
          cat $WORKDIR/${PDBID}_ligfit.pdb $WORKDIR/append.pdb | grep -v '^END' > $WORKDIR/newlig.pdb
          #Remove offending waters
          if ($NWATER > 0 ) then 
            echo "Removing waters for candidate $RESNUM" >> $WORKDIR/ligfit.log
            cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.centrifuge.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
            
            $TOOLS/centrifuge \
            <<eof >& $WORKDIR/centrifuge_lig.log
              #Settings
              InputPDB = "$WORKDIR/newlig.pdb"
              PDBOutputFilename = "$WORKDIR/newlig_dry.pdb"
              MessageFilename = "$WORKDIR/newlig_centrifuge.msg"
              XMLOutputFilename = "$WORKDIR/newlig_centrifuge.xml"
      
              #residue names that represent waters. Default HOH WAT H2O EAU
              WaterNames = HOH
              #Reject waters with negative density (i.e. only clashing waters)
              WaterRejectionThreshold = 0.00
              #If true, everything except Dummies and recognised Waters is placed in the density map, default False
              PlaceNonSolventInDensity = true
              #Colon separated list of water ids which won't be affected by the methods of centrifuge
              UnaffectedWaterList = "$H2O_KEEP"
              #InputMap or here: input mtz
              InputMTZ = "$WORKDIR/${PDBID}_loopwhole.mtz"
              FWTLabel = FWT
              PHIWTLabel = PHWT
              SpaceGroup = "$SPACEGROUP"
              XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"
              DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
              SymmetryFilename = "$CLIBD/syminfo.lib"

              #type of method used to determine the density fit
              CFitTarget::Type = accelerated
              UniformAtomBFactor = "$BBUILD"
              UniformAtomRadius = 0.74
              UseAtomicB = True

              #Factor with which to invert the density of placed atoms
              AntiBumpFactor = 2.0

              #Program information
              ProgramName = centrifuge
              MessageLevel = 6
              AbortLevel = 8
              KeepPDBheaderInfo = true
              StoreOriginalChainAndSegID = true
              SelectByOrigChainID = true
              KeepSideChain = true
              KeepFragmentLongerThan = 0
              TrustWaters = false
              RemoveBasedOnDensity = true
eof
            cat $WORKDIR/centrifuge_lig.log >> $WORKDIR/ligfit.log
          else
            cp $WORKDIR/newlig.pdb $WORKDIR/newlig_dry.pdb
          endif

          #Real-space refine against 2mFo-DFc
          echo "Refining candidate $RESNUM" >> $WORKDIR/ligfit.log
          cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software."coot-mini-rsr".used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
          
          echo "coot-mini-rsr --pdbin $WORKDIR/newlig_dry.pdb --hklin $WORKDIR/${PDBID}_loopwhole.mtz --f FWT --phi PHWT --pdbout $WORKDIR/newlig_refi.pdb --chain-id z --residues-around $RESNUM --weight 40 --rama --torsions --no-trans-peptide-restraints --correlations --dictin $LIGREST"  >> $WORKDIR/ligfit.log

          coot-mini-rsr \
          --pdbin $WORKDIR/newlig_dry.pdb \
          --hklin $WORKDIR/${PDBID}_loopwhole.mtz \
          --f FWT --phi PHWT \
          --pdbout $WORKDIR/newlig_refi.pdb \
          --chain-id z --residues-around $RESNUM \
          --weight 40 \
          --rama --torsions --no-trans-peptide-restraints \
          --correlations \
          --dictin $LIGREST > $WORKDIR/ligref.log
          
          #Calculate density fit on updated maps
          $TOOLS/stats \
          --use-auth-ids \
          --sampling-rate 1.5 \
          --hklin $PDBID.mtz \
          $SFTYPE \
          --recalc \
          --xyzin $WORKDIR/newlig_refi.pdb \
          -o $WORKDIR/density_$RESNUM.txt \
          $DICTCMD >>& $WORKDIR/density_ligfit.log 
        
          #Save the scores
          grep "${CHID}_${RESNUM}[[:space:]]" $WORKDIR/density_$RESNUM.txt >> $WORKDIR/ligfit.scores         

          #Decide whether to keep the ligand (remove ligands with too low correlation, too low RSCC, too low EDIAm, or too high RSR)
          if (`grep -A 40 'Residue Correlation Table:' $WORKDIR/ligref.log | grep $LIG | grep -e "$RESNUM" | awk '{if ($4 < 0.53) {print "reject"} else {print "accept"}}'` == "accept" && `grep "${CHID}_${RESNUM}[[:space:]]" $WORKDIR/density_$RESNUM.txt | awk '{if ($2 < 0.21 && (($4 > 0.76 && $6 > 0.50) || ($4 + $6 > 1.1) || ($2 < 0.185 && ($4 + $6 > 1.0)))) {print 1} else {print 0}}'` == 1) then
                
            #keep the new ligand
            set BUILT = 1
            cp $WORKDIR/newlig_refi.pdb $WORKDIR/${PDBID}_ligfit.pdb
            
            #Append it to the validation list
            echo "   * Accepting candidate $FITCOUNT" | tee -a $LOG
            if (`echo $NLIG_LIST | grep -c ':'` == 0) then
              set NLIG_LIST = ':'
            endif
            set NLIG_LIST = "${NLIG_LIST}$CHID$RESNUM :"  
          
            #Occupancy refine the new ligand
            @ OCCREF = ($OCCREF + 1)
            echo "occupancy group id $OCCREF chain $CHID residue $RESNUM" >> $WORKDIR/occupancy_cmd.refmac

            #Increase counter
            @ RESNUM   = ($RESNUM + 1)
            @ ACCEPTED = ($ACCEPTED + 1)
          
          else
            #Keep using the previous model
            echo "   * Rejecting candidate $FITCOUNT" | tee -a $LOG
          endif
          
          #Consolidate log files
          mv  $WORKDIR/ligfit.log $WORKDIR/ligfit.old
          cat $WORKDIR/ligfit.old $WORKDIR/ligref.log > $WORKDIR/ligfit.log
        end
      else
        echo " o Could not find a good place to fit $LIG" | tee -a $LOG
      endif        
     
    end

    cp $WORKDIR/${PDBID}_ligfit.pdb $WORKDIR/${PDBID}_built.pdb
  
    #Is occupancy refinement needed?
    if (-e $WORKDIR/occupancy_cmd.refmac) then

      #Append refinement command
      echo "occupancy refine" >> $WORKDIR/occupancy_cmd.refmac

      #Make Refmac use command file
      set OCCCMD = "@$WORKDIR/occupancy_cmd.refmac"
    endif
  endif
endif  

############################################### Do another refinement run ################################################

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Final refinement ******" | tee -a $LOG

#Skip the final refinement if the model was not rebuilt
unbuilt:

if ($BUILT == 0) then

  #Report
  echo "-Skipping the final refinement, because the model was not changed after re-refinement"  | tee -a $LOG

  #Copy some files
  cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_final.pdb
  cp $WORKDIR/${PDBID}_besttls.mtz $WORKDIR/${PDBID}_final.mtz
  cp $WORKDIR/${PDBID}_besttls.log $WORKDIR/${PDBID}_final.log
  if ($DOTLS == 1) then
    cp $WORKDIR/${PDBID}_refin.tls   $WORKDIR/${PDBID}_final.tls
  endif
  set BLTBEST = 'skip'

  #Copy refinement results
  set RFIN       = $RTLS
  set RFFIN      = $RFTLS
  set RFFINUNB   = $RFTLSUNB
  set SIGRFFIN   = $SIGRFTLS
  set RFFINZ     = $RFTLSZ
  set ZRFRRATFIN = $ZRFRRATTLS

else
  #Update restraints
  #Nucleic acid restraints(only if the previous run had no problems)
  if (-e $WORKDIR/nucleic.rest && $NPRUNEN == 0) then
    echo "-Updating nucleic acid restraints" | tee -a $LOG
    #Remove old restraints
    if (-e $WORKDIR/bphbond.rest) then
      rm $WORKDIR/bphbond.rest
    endif
    
    libg -p $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/tstack.rest -w dsp >> $WORKDIR/libg.log
    #Clean up stacking restraints
    awk '{if (NF < 68) {print $0}}' $WORKDIR/tstack.rest | grep -v 'plan 3' > $WORKDIR/stacking.rest
    
    #Generate base pair H-bond restraints
    python3 $TOOLS/bphbond.py $WORKDIR/${PDBID}_built.pdb $TOOLS/x3dna-dssr > $WORKDIR/bphbond.tmp

    #Clean the restraints
    #Make list of all atoms with alternates
    grep -v LINK $WORKDIR/${PDBID}_built.pdb > $WORKDIR/built.pdb
    echo "SELECT auth_asym_id,auth_seq_id,auth_atom_id FROM atom_site WHERE label_alt_id IS NOT NULL;" | $TOOLS/mmCQL $WORKDIR/built.pdb >& $WORKDIR/allalt.list
    rm $WORKDIR/built.pdb
    
    #Loop over all restraints
    foreach REST (`cat $WORKDIR/bphbond.tmp | tr ' ' '_'`) 
      set CHID1 = `echo $REST | cut -d '_' -f 5`
      set RESN1 = `echo $REST | cut -d '_' -f 7`
      set ATOM1 = `echo $REST | cut -d '_' -f 11`
      set CHID2 = `echo $REST | cut -d '_' -f 14`
      set RESN2 = `echo $REST | cut -d '_' -f 16`
      set ATOM2 = `echo $REST | cut -d '_' -f 20`
      
      #Only write out restraint if neither party is in the list of atoms with alternates
      if (`grep $CHID1 $WORKDIR/allalt.list | grep -w -e "$RESN1" | grep -c "$ATOM1"$` == 0 && `grep $CHID2 $WORKDIR/allalt.list | grep -w -e "$RESN2" | grep -c "$ATOM2"$` == 0) then
        echo $REST | tr '_' ' ' >> $WORKDIR/bphbond.rest
      else
        #echo $REST
      endif
    end
    #Make Refmac use command files
    echo "external weight scale 5" >  $WORKDIR/nucleic.rest
    if (-e $WORKDIR/stacking.rest) then
      sed 's/  / /g' $WORKDIR/stacking.rest >> $WORKDIR/nucleic.rest
      rm $WORKDIR/stacking.rest
    endif
    if (-e $WORKDIR/bphbond.rest) then
      echo "external weight scale 2" >> $WORKDIR/nucleic.rest
      cat $WORKDIR/bphbond.rest >> $WORKDIR/nucleic.rest
    endif
    echo "external weight scale 1" >> $WORKDIR/nucleic.rest
  else if (-e $WORKDIR/bphbond.rest) then
    #Only update the targets for validation
    python3 $TOOLS/bphbond.py $WORKDIR/${PDBID}_built.pdb $TOOLS/x3dna-dssr > $WORKDIR/bphbond.rest
  endif

  #Hydrogen bonds
  if (-e $WORKDIR/hbond.rest && ! -e $WORKDIR/homology.rest) then
    echo "-Updating hydrogen bond restraints" | tee -a $LOG

    #Define the secondary structure
    mkdssp -i $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/${PDBID}_built.dssp >& $WORKDIR/dssp.log

    #Generate hydrogen bond restraints if a DSSP file can be made.
    if (! -e $WORKDIR/${PDBID}_built.dssp) then
      #DSSP failed. Cannot make restraints
      echo " o Cannot produce hydrogen bond restraints" | tee -a $LOG
      echo "DSSP: general error in second H-bond restraint generation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                                           >> $DEBUG
    else
      #Make the hydrogen bond restraints. PROGRAM: detectHbonds
      $TOOLS/detectHbonds -v \
      -pdb $WORKDIR/${PDBID}_built.pdb \
      -dssp $WORKDIR/${PDBID}_built.dssp \
      -output-name $PDBID \
      -tools $TOOLS >> $WORKDIR/hbondrest.log
      if ( ! -e $WORKDIR/${PDBID}_hbonds.rest) then
        echo " o Cannot generate hydrogen bond restraints for rebuilt model" | tee -a $LOG
        echo "detectHbonds: general error" >> $DEBUG
        echo "PDB-REDO,$PDBID"             >> $DEBUG
      else
        #Set up restraint commands
        sed 's/ [ ]*/ /g' $WORKDIR/${PDBID}_hbonds.rest > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
        rm $WORKDIR/${PDBID}_hbonds.rest
        set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
        set HBONDCMD    = "@$WORKDIR/hbond.rest"
        rm $WORKDIR/${PDBID}_built.dssp
      endif
    endif
  endif

  #Homology restraints
  if (-e $WORKDIR/homology.rest) then
    echo "-Updating homology-based restraints" | tee -a $LOG

    #Define the secondary structure
    mkdssp -i $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/${PDBID}_built.dssp >& $WORKDIR/dssp.log

    #Make the restraints
    cd $WORKDIR/homol

    $TOOLS/hoder -v \
    -alignout \
    -hitsummary \
    $HOMINCMD \
    -pdb $WORKDIR/${PDBID}_built.pdb \
    -dssp $WORKDIR/${PDBID}_built.dssp \
    -blast $WORKDIR/$PDBID.blast \
    -fasta $WORKDIR/$PDBID.fasta \
    -output-name $PDBID \
    -output-dir $WORKDIR/homol \
    -pdbdir $REDODIR \
    -edsdir $EDSDIR \
    -tools $TOOLS >> homologs.log

    #Set up the restraint commands
    #Homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/homol/${PDBID}_homologyBased.restr > $WORKDIR/homology.rest #Replace one or more spaces by a single one
    set HOMOLWGTCMD = "EXTERNAL WEIGHT SCALE $HOMOLRESTRWGT"
    set HOMOLCMD    = "@$WORKDIR/homology.rest"

    #Non-homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/homol/${PDBID}_general.restr > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
    set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
    set HBONDCMD    = "@$WORKDIR/hbond.rest"

    #Go back down
    cd $WORKDIR
    rm $WORKDIR/${PDBID}_built.dssp
  endif

  #Metal restraints
  if ($DOMETALREST == 1 && $NPRUNEM == 0) then

    #Run Platonyzer. PROGRAM: Platonyzer
    echo "-Updating metal restraints" | tee -a $LOG
    if (-e $WORKDIR/metal.rest) then
      rm $WORKDIR/metal.rest
    endif  
  

    $TOOLS/platonyzer -v \
    $WORKDIR/${PDBID}_built.pdb \
    -o $WORKDIR/${PDBID}_built_platonyzed.pdb \
    --delete-vdw-rest \
    >>& $WORKDIR/platonyzer.log

  
    if (-e $WORKDIR/${PDBID}_built_platonyzed.restraints) then
      cp $WORKDIR/${PDBID}_built_platonyzed.pdb $WORKDIR/${PDBID}_built.pdb
    
      #Set up the use of the restraints
      if (`cat $WORKDIR/${PDBID}_built_platonyzed.restraints | wc -l` > 0) then
        mv $WORKDIR/${PDBID}_built_platonyzed.restraints $WORKDIR/metal.rest
        set METALCMD = "@$WORKDIR/metal.rest"
      endif
    else  
      set METALCMD = ""
    endif

    #Check to see wether there are new metal restraints
    if (-e $WORKDIR/metal.rest) then
      set NMETALREST2 = `grep -c exte $WORKDIR/metal.rest`
      if ($NMETALREST2 != $NMETALREST) then
        echo " o Found different metal sites"      | tee -a $LOG
        echo " o Performing extra refinement using $NMETALREST2 restraints" | tee -a $LOG
      endif
    endif
  endif

  #set refinement parameters
  if ($TLSBEST == none) then
    #Autoweight
    set WGTTYPE = "AUTO"
    set WGRANGE = "2.50"
    set FCYCLE  = `echo $NCYCLE`
    set BLTBCMD = `echo $TBCMD`
  else
    set FCYCLE  = 20
    set BLTBCMD = 
    #Increase the number of cycles if needed
    if ($NMETALREST2 != $NMETALREST) then
      @ FCYCLE = ($FCYCLE + 5)
    endif
    #Do short weight optimisation
    set WGTTYPE = "MATRIX"
    if ($TLSBEST == 5.00) then
      set WGRANGE = "5.00 7.00 3.00"
    else if ($TLSBEST == 3.00) then
      set WGRANGE = "3.00 4.00 2.00"
    else if ($TLSBEST == 2.00) then
      set WGRANGE = "2.00 2.50 1.50"
    else if ($TLSBEST == 1.50) then
      set WGRANGE = "1.50 1.80 1.00"
    else if ($TLSBEST == 1.00) then
      set WGRANGE = "1.00 1.30 0.70"
    else if ($TLSBEST == 0.70) then
      set WGRANGE = "0.70 0.80 0.50"
    else if ($TLSBEST == 0.50) then
      set WGRANGE = "0.50 0.60 0.30"
    else if ($TLSBEST == 0.30) then
      set WGRANGE = "0.30 0.40 0.10"
    else if ($TLSBEST == 0.10) then
      set WGRANGE = "0.10 0.15 0.05"
    else if ($TLSBEST == 0.05) then
      set WGRANGE = "0.05 0.07 0.03"
    else if ($TLSBEST == 0.03) then
      set WGRANGE = "0.03 0.05 0.02"
    else if ($TLSBEST == 0.01) then
      set WGRANGE = "0.01 0.02 .005"
    else if ($TLSBEST == .005) then
      set WGRANGE = ".005 .007 .002"
    else if ($TLSBEST == .002) then
      set WGRANGE = ".002 .005 .001"
    else if ($TLSBEST == .001) then
      set WGRANGE = ".001 .002 5e-4"
    else if ($TLSBEST == 5e-4) then
      set WGRANGE = "5e-4 .001 1e-4"
    else if ($TLSBEST == 1e-4) then
      set WGRANGE = "1e-4 2e-4 5e-5"
    else if ($TLSBEST == 1e-5) then
      set WGRANGE = "1e-5 2e-5 5e-6"
    else if ($TLSBEST == 1e-6) then
      set WGRANGE = "1e-6 2e-6 5e-7"
    else
      set WGRANGE = "1e-7 2e-7 5e-8"
    endif
  endif

  #Set up TLS and killswitch
  if ($DOTLS == 1 && $TLSUPDATE == 1) then
    set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_refin.tls TLSOUT $WORKDIR/${PDBID}_final.tls`
    #Do only a minor update of the TLS tensors (avoid instability)
    set TLSCMD     = `echo refi tlsc 5`
    set KILLSWITCH = "kill $TOOLS/pdb_redo.refmac"
  else
    set KILLSWITCH =    #Do not use the killswitch
    set TLSCMD     =    #Do not refine TLS
  endif

  #Run refmac with predefined matrix weights
  echo "-Final restraint weight optimisation with $WGTTYPE weights: $WGRANGE" | tee -a $LOG

  #Return label for running after TLS troubles
notlsblt:

  foreach WWGT (`echo $WGRANGE`)

    #Return label for job launching
bltrunning:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log

    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_built.pdb \
      XYZOUT $WORKDIR/${PDBID}_blt$WWGT.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_blt$WWGT.mtz \
      $LIBLIN \
      ${TLSFILS} \
      $SCATLIN \
<<eof >& $WORKDIR/${PDBID}_blt$WWGT.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        make newligand continue
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        $TLSCMD
        $BLTBCMD
        tlsd waters exclude
        ncyc $FCYCLE
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG $WGTTYPE $WWGT
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $OCCCMD
        $HARMCMD
        $METALCMD
        $NUCRCMD
        $RESTCMD
        temp $BBEST
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
        $ANOMCMD
        pdbout copy expdta
        pdbout copy remarks 200 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        $KILLSWITCH
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto bltrunning
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists
  foreach WWGT (`echo $WGRANGE`)
    if (! -e $WORKDIR/${PDBID}_blt$WWGT.mtz) then
      #Was the killswitch activated
      if (`grep -a -c 'Program terminated by user' $WORKDIR/${PDBID}_blt$WWGT.log` != 0) then
        #Yes, check if TLS was used.
        if ($DOTLS == 1 && $TLSUPDATE == 1) then
          #Try running without updating the TLS tensors.
          set TLSFILS    = `echo TLSIN $WORKDIR/${PDBID}_refin.tls`
          set TLSCMD     =
          set TLSUPDATE  = 0
          set KILLSWITCH =    #Do not use the killswitch

          #Give warning and write debug statement.
          echo " " | tee -a $LOG
          echo " o The refinement was unstable." | tee -a $LOG
          echo "   * Trying refinement without updating the TLS model " | tee -a $LOG

          #Rerun Refmac without updating the TLS
          goto notlsblt
        else if
          #Check for atom naming problems
          if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log` != 0) then
            #Is this caused by an alternative compound
            set RESN = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 39-42`
            set CHID = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 45-45`
            #Check whether there is more than 1 residue name for a given chain and residue number.
            if (`grep -E $CHID' {0,3}'$RESN $WORKDIR/${PDBID}_besttls.pdb | cut -c 18-20 | sort -u | wc -l` > 1) then
             
              #Write warning
              echo " " | tee -a $LOG
              echo " o The rebuilt model caused problems in Refmac" | tee -a $LOG
              echo "   * Falling back to the re-refined model"      | tee -a $LOG
              echo "COMMENT: refmac: could not use rebuilt model" >> $DEBUG
              echo "PDB-REDO,$PDBID"                              >> $DEBUG
              
              #Undo the rebuilding
              set NWATDEL = 0
              set NBBFLIP = 0
              set NSCBLT  = 0
              set NSCFLIP = 0
              set NCHIRFX = $CHIFIX
              set BUILT   = 0

              #Start again form the pre-rebuilding model
              goto unbuilt
              
            else
              #Fatal atom naming problem
              echo " " | tee -a $LOG
              echo " o Problem with restraint generation in refmac. Cannot continue." | tee -a $LOG
              echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
              echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
              if ($SERVER == 1) then
                #Write out status files
                touch $STDIR/stoppingProcess.txt
                touch $STDIR/processStopped.txt
              endif
              exit(1)
            endif
          endif
        else
          #TLS was not the problem. The error is fatal.
          echo " " | tee -a $LOG
          echo " o Problem with refmac. Cannot continue." | tee -a $LOG
          echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
          if ($SERVER == 1) then
            #Write out status files
            touch $STDIR/stoppingProcess.txt
            touch $STDIR/processStopped.txt
          endif
          exit(1)
        endif
      else
        #Check for atom naming problems
        if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log` != 0) then
          #Is this caused by an alternative compound
          set RESN = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 39-42`
          set CHID = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 45-45`
          #Check whether there is more than 1 residue name for a given chain and residue number.
          if (`grep -E $CHID' {0,3}'$RESN $WORKDIR/${PDBID}_besttls.pdb | cut -c 18-20 | sort -u | wc -l` > 1) then
           
            #Write warning
            echo " " | tee -a $LOG
            echo " o The rebuilt model caused problems in Refmac" | tee -a $LOG
            echo "   * Falling back to the re-refined model"      | tee -a $LOG
            echo "COMMENT: refmac: could not use rebuilt model" >> $DEBUG
            echo "PDB-REDO,$PDBID"                              >> $DEBUG
              
            #Undo the rebuilding
            set NWATDEL = 0
            set NBBFLIP = 0
            set NSCBLT  = 0
            set NSCFLIP = 0
            set NCHIRFX = $CHIFIX
            set BUILT   = 0

            #Start again form the pre-rebuilding model
            goto unbuilt
        
          else
            #Fatal atom naming problem
            echo " " | tee -a $LOG
            echo " o Problem with restraint generation in refmac. Cannot continue." | tee -a $LOG
            echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
            exit(1)
          endif
        else
          #Something else is wrong, give an error mesage and quit.
          echo " " | tee -a $LOG
          echo " o Problem with refmac. Cannot continue." | tee -a $LOG
          echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
          if ($SERVER == 1) then
            #Write out status files
            touch $STDIR/stoppingProcess.txt
            touch $STDIR/processStopped.txt
          endif
          cd $BASE
          exit(1)
        endif  
      endif
    endif
  end

  #Report
  foreach WWGT (`echo $WGRANGE`)
    set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | awk '{print $3}'`
    echo " o Performed ${FCYCLE}-cycle refinement with $WGTTYPE weight $WWGT (R-free = $TFREE)" | tee -a $LOG
  end

  #Give refinement output
  #No known errors in the refinement
  set BLTERR  = 0

  #If the refinement was done with automated weighting, just rename the files...
  if ($WGTTYPE == "AUTO") then
    set BLTBEST = 'auto'
    mv $WORKDIR/${PDBID}_blt$WWGT.pdb $WORKDIR/${PDBID}_final.pdb
    mv $WORKDIR/${PDBID}_blt$WWGT.mtz $WORKDIR/${PDBID}_final.mtz
    mv $WORKDIR/${PDBID}_blt$WWGT.log $WORKDIR/${PDBID}_final.log
  else
    #Create second.refi file
    if ($GOTR == 0 || $ZCALERR == 1) then
      echo "$RCAL $RFCALUNB" > $WORKDIR/${PDBID}blt.refi  #Use recalculated R and expected R-free as benchmarks
    else
      echo "$RCAL $RFCAL" > $WORKDIR/${PDBID}blt.refi  #Use R and R-free obtained from the recalculation as benchmarks
    endif
    foreach WWGT (`echo $WGRANGE`)
      if (-e $WORKDIR/${PDBID}_blt$WWGT.pdb) then
        set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1`
        echo "$WWGT $LINE" >> $WORKDIR/${PDBID}blt.refi
      else
        echo " o $WORKDIR/${PDBID}_blt$WWGT.pdb is missing" | tee -a $LOG
        set BLTERR = 1
      endif
    end

    #Pick the best re-refined structure
    set BLTBEST = `$TOOLS/picker -s -f $ESTRICT $WORKDIR/${PDBID}blt.refi $NTSTCNT $RFRRAT $SRFRRAT $RMSZB $RMSZA`
    if ($BLTBEST == 'none') then
      #Something is very wrong, give an error mesage and quit.
      echo " " | tee -a $LOG
      echo " o Problem with refmac: no suitable models. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: post-rebuilding refinement gave unsuitable models" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                                    >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
    cp $WORKDIR/${PDBID}_blt$BLTBEST.pdb $WORKDIR/${PDBID}_final.pdb
    cp $WORKDIR/${PDBID}_blt$BLTBEST.mtz $WORKDIR/${PDBID}_final.mtz
    cp $WORKDIR/${PDBID}_blt$BLTBEST.log $WORKDIR/${PDBID}_final.log
  endif

  #Set R(-free)
  set RFIN  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_final.log | head -n 1 | awk '{print $2}'`
  set RFFIN = `tail -n $LOGSTEP $WORKDIR/${PDBID}_final.log | head -n 1 | awk '{print $3}'`

  #Validate R-free
  set SIGRFFIN   = `echo $RFFIN $NTSTCNT | awk '{printf ("%.4f\n", $1/sqrt($2))}'`
  set RFFINUNB   = `echo $RFIN $RFRRAT | awk '{printf ("%.4f\n", $1*$2)}'`
  set RFFINZ     = `echo $RFFINUNB $RFFIN $SIGRFFIN | awk '{printf ("%.2f\n", ($1-$2)/$3)}'`
  set ZRFRRATFIN = `echo $RFRRAT $SRFRRAT $RFIN $RFFIN | awk '{printf ("%.2f\n", ($1-($4/$3))/$2)}'`

  echo " " | tee -a $LOG
  echo " " | tee -a $LOG

  echo "****** Final refinement details ******" | tee -a $LOG
  if ($WGTTYPE == 'MATRIX') then
    echo "Best matrix weight: $BLTBEST"  | tee -a $LOG
  endif
  echo "Resulting R-factor: $RFIN"       | tee -a $LOG
  echo "Resulting R-free  : $RFFIN"      | tee -a $LOG
  echo "R-free/R Z-score  : $ZRFRRATFIN" | tee -a $LOG
  echo " " | tee -a $LOG
  echo "Expected R-free   : $RFFINUNB"   | tee -a $LOG
  echo "sigma(R-free)     : $SIGRFFIN"   | tee -a $LOG
  echo "R-free Z-score    : $RFFINZ"     | tee -a $LOG
endif

#Copy back the SEQRES records
if ($GOTSEQRES == 1  && `grep -c '^SEQRES' $WORKDIR/${PDBID}_final.pdb` == 0) then
  cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final.bak
  #Run seqrescopier. PROGRAM: seqrescopier
  $TOOLS/seqrescopier -v \
  -pdbinw  $WORKDIR/cache.pdb \
  -pdbinwo $WORKDIR/${PDBID}_final.bak \
  -pdbout  $WORKDIR/${PDBID}_final.pdb >> $WORKDIR/seqrescopier.log 
else
  #Do nothing
endif 


#Run flipper again
grep -v '^TER' $WORKDIR/${PDBID}_final.pdb | grep -v -E '^LINKR.{67}gap' > $WORKDIR/${PDBID}_final.bak
$TOOLS/flipper -v \
$WORKDIR/${PDBID}_final.bak \
$DICTCMD \
-o $WORKDIR/${PDBID}_final.pdb >>& $WORKDIR/flipper.log

###########################################  Correlation coefficients  ###################################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG

#Extract the data
set CCWOLD = `grep 'CORRELATION COEFFICIENT FO-FC     ' $WORKDIR/${PDBID}_0cyc.pdb  | awk '{print $7}'`
set CCWFIN = `grep 'CORRELATION COEFFICIENT FO-FC     ' $WORKDIR/${PDBID}_final.pdb | awk '{print $7}'`
set CCFOLD = `grep 'CORRELATION COEFFICIENT FO-FC FREE' $WORKDIR/${PDBID}_0cyc.pdb  | awk '{print $8}'`
set CCFFIN = `grep 'CORRELATION COEFFICIENT FO-FC FREE' $WORKDIR/${PDBID}_final.pdb | awk '{print $8}'`

#Get the Z-scores
set ZCCW = `echo "$CCWOLD $CCWFIN $WORKCNT" | awk '{zcco = 0.5*(log((1+$1)/(1-$1))); zccf = 0.5*(log((1+$2)/(1-$2))); zchange = (zccf-zcco)/sqrt(2/($3-3))} END {printf "%6.2f\n", zchange}'`
set ZCCF = `echo "$CCFOLD $CCFFIN $NTSTCNT" | awk '{zcco = 0.5*(log((1+$1)/(1-$1))); zccf = 0.5*(log((1+$2)/(1-$2))); zchange = (zccf-zcco)/sqrt(2/($3-3))} END {printf "%6.2f\n", zchange}'`

#Report the results
echo "****** Reciprocal space correlation ******" | tee -a $LOG
echo '                              Before  Final  Z-score'    | tee -a $LOG
echo "Work correlation coefficient: $CCWOLD   $CCWFIN  $ZCCW"  | tee -a $LOG
echo "Free correlation coefficient: $CCFOLD   $CCFFIN  $ZCCF"  | tee -a $LOG



############################################   Full cross validation  ####################################################

#Do full cross-validation if the R-free set is too small or if asked by the user...
# ... but only if there was a successful refinement
if ( ($NTSTCNT < 500 || $CROSS == 1) && !($BLTBEST == 'skip' && $TLSBEST == 'none')) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Cross validation ******" | tee -a $LOG
  if ($NTSTCNT < 500) then
    echo "-The test set contained fewer than 500 reflections" | tee -a $LOG
    set CROSS = 1
  endif

  #Calculate K
  set KFOLD = `mtzdmp $WORKDIR/$PDBID.mtz -n 0 | grep '  FREE' | cut -c 21-22 | awk '{print $1+1}'`
  @ MAXSET  = ($KFOLD - 1)
  echo "-Performing $KFOLD-fold cross validation" | tee -a $LOG

  #Set up the refinement
  if ($DOTLS == 1) then
    set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_final.tls`
  endif
  set BCMD = "bfac set $BSET"
  @ NCYCLE = ($NCYCLE + 40)
  #Keep the final refinement settings unless there was no final refinement
  if ($BLTBEST == 'skip') then
    set CWEIGHT = `echo "$WGTSIG MATRIX $TLSBEST"`
  else if ($BLTBEST == 'auto') then
    set CWEIGHT = `echo "$WGTSIG AUTO 2.50"`
  else
    set CWEIGHT = `echo "$WGTSIG $WGTTYPE $BLTBEST"`
  endif

  #Perturb the atomic coordinates if the B-factor model is OVER
  if ($BREFTYPE == "OVER") then
    #Run pdbset. PROGRAM: pdbset
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.pdbset.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    pdbset \
    XYZIN $WORKDIR/${PDBID}_final.pdb \
    XYZOUT $WORKDIR/${PDBID}_crossin.pdb \
    <<eof > $WORKDIR/pdbset.log
      NOISE 0.05
eof
  else
    cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_crossin.pdb
  endif

  #Do the refinements
  foreach SET (`seq -s " " 0 $MAXSET`)

    #Return label for job launching
xvalrunning:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log

    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      echo " o Starting refinement with test set $SET" | tee -a $LOG

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_crossin.pdb \
      XYZOUT $WORKDIR/cross$SET.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/cross$SET.mtz \
      $LIBLIN \
      ${TLSFILS} \
      $SCATLIN \
  <<eof > $WORKDIR/${PDBID}_cross$SET.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        make newligand continue
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        free $SET
        $BCMD
        tlsd waters exclude
        ncyc $NCYCLE
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $CWEIGHT
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $NUCRCMD
        $RESTCMD
        temp $BBEST
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy expdta
        pdbout copy remarks 200 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        END
eof
     else
      #Wait a bit to start again
      sleep 10
      goto xvalrunning
    endif
  end

  #Wait for the jobs to finish
  wait

  #Make a summary only if the refinement was successful
  echo "-Analysing cross validation results" | tee -a $LOG
  foreach SET (`seq -s " " 0 $MAXSET`)
    if (! -e $WORKDIR/cross$SET.mtz) then
      echo " o Problem in cross validation using test set $SET" | tee -a $LOG
      echo "COMMENT: refmac: error in cross validation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                            >> $DEBUG
    else
      set LINE  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_cross$SET.log | head -n 1`
      set NTREF = `mtzdmp $WORKDIR/$PDBID.mtz -n -1 | grep -v '?' | grep -E "^.{4}[0123456789].{3}[0123456789]" | awk '{printf "%d\n", $4}' | grep -x -c "$SET"` #Grep -x to only get exact matches
      echo "$SET $LINE $NTREF" >> $WORKDIR/${PDBID}.rtest
    endif
  end

  #Analyse results and show them

  #R-complete
  foreach SET (`seq -s " " 0 $MAXSET`)

    #Copy back the testset flags (circumvents a Refmac quirk)
    cad \
    HKLIN1 $WORKDIR/cross${SET}.mtz \
    HKLIN2 $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/crossf${SET}.mtz \
    <<eof >> $WORKDIR/rcomplete.log
    LABIN FILE 1  E1=FP E2=FC_ALL_LS E3=PHIC_ALL_LS
    LABIN FILE 2  E1=FREE
    END
eof

    #Convert the file to mmCIF
    mtz2various \
    HKLIN  $WORKDIR/crossf${SET}.mtz \
    HKLOUT $WORKDIR/cross${SET}.cif << eof >> $WORKDIR/rcomplete.log
    OUTPUT CIF data_rcom
    FREEVAL $SET.0
    LABIN FP=FP FC=FC_ALL_LS PHIC=PHIC_ALL_LS FREE=FREE
    END
eof

    #Extract and consolidate the data
    grep ' 1 1 1' $WORKDIR/cross${SET}.cif | grep ' f ' >> $WORKDIR/rcomplete.hkl
  end

  #Get the R-complete
  set RCOMPLETE =  `cat $WORKDIR/rcomplete.hkl | awk '{if ($8 > $9) {sumdif += ($8 - $9)} else {sumdif += ($9 - $8)}; sumf += $8} END {printf "%6.4f", sumdif/sumf}'`
  echo " o R-complete = $RCOMPLETE" | tee -a $LOG

  #Averages
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.longinus.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
  
  $TOOLS/longinus -v \
  $WORKDIR/$PDBID.rtest > $WORKDIR/longinus.log

  #Did the crossvalidation work
  if (`grep -c 'Could not find any valid refinement data' $WORKDIR/longinus.log` != 0) then
    #Cross validation failed write an error message
    set CROSS = 0
    echo "-Cross validation failed" | tee -a $LOG
    echo "COMMENT: error in cross validation " >> $DEBUG
    echo "PDB-REDO,$PDBID"                     >> $DEBUG
  else
    #Crossvalidation successful. Write out details.
    echo " " | tee -a $LOG
    echo "****** Cross validation details ******" | tee -a $LOG
    tail -n 3 $WORKDIR/longinus.log | tee -a $LOG

    #Mine the values
    set CRFACT  = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 10-15`
    set CSRFACT = `tail -n 1 $WORKDIR/longinus.log | cut -c 10-15`
    set CRFREE  = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 20-25`
    set CSRFREE = `tail -n 1 $WORKDIR/longinus.log | cut -c 20-25`
    set CGAP    = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 28-33`
    set CSGAP   = `tail -n 1 $WORKDIR/longinus.log | cut -c 28-33`
    set CNTEST  = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 35-40`
    set CSNTEST = `tail -n 1 $WORKDIR/longinus.log | cut -c 35-40`

    #Clean up
    foreach SET (`seq -s " " 0 $MAXSET`)
      rm $WORKDIR/cross$SET.mtz
      rm $WORKDIR/crossf$SET.mtz
      rm $WORKDIR/cross$SET.pdb
      rm $WORKDIR/cross$SET.cif
    end
  endif
endif

############################################   Validate structure   ######################################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Model validation ******" | tee -a $LOG

echo "-Running DSSP" | tee -a $LOG
mkdssp -i $WORKDIR/${PDBID}_final.pdb -o $WORKDIR/${PDBID}_final.dssp >>& $WORKDIR/dssp.log

echo "-Running WHAT_CHECK" | tee -a $LOG


#Go to temporary running directory
setenv WCWORF $WORKDIR/wctemf
mkdir -p $WCWORF
cd $WCWORF

#Get the PDB file
cp $WORKDIR/${PDBID}_final.pdb $WCWORF

#Do the actual validation
$WC/bin/whatcheck $WCWORF/${PDBID}_final.pdb Y Y Y >& $WORKDIR/wc_final.log

#Check for an output file
if (-e $WCWORF/pdbout.txt) then
  #Do Nothing
else
  #Give warning
  echo " o Validation failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate final structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                     >> $DEBUG
endif

#Create webpage
$WC/bin/pdbout2html >>& $WORKDIR/wc_final.log

#Check index.html completeness
if (`grep -c "Final summary" $WCWORF/pdbout.html` == 0 && `grep -c "Summary report" $WCWORF/pdbout.html` == 0) then
  echo " o Validation failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate final structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                     >> $DEBUG
endif

#Clean up on aisle six
mkdir -p $WORKDIR/wf

if (-e $WCWORF/pdbout.txt) then
  mv $WCWORF/pdbout.txt $WORKDIR/wf/
endif
#server only
if ($SERVER == 1) then
  mv $WCWORF/pdbout.html $WORKDIR/wf/index.html
  mv $WCWORF/*.gif $WORKDIR/wf/ >& /dev/null
endif

#Go to start directory
cd $WORKDIR
rm -rf $WCWORF

echo "-Validating the final model with tortoize" | tee -a $LOG
$TOOLS/tortoize --xyzin $WORKDIR/${PDBID}_final.pdb --output $WORKDIR/${PDBID}_final_tortoize.json >>& $WORKDIR/tortoize.log

#Extract statistics from refitted structure
set FNATOM = `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_final.pdb`
if (-e $WORKDIR/wf/pdbout.txt) then
  if (`grep -a -c '1st generation packing quality :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FZPAK1 = 'NA'
  else
    set FZPAK1 = `grep -a '1st generation packing quality :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c '2nd generation packing quality :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FZPAK2 = 'NA'
  else
    set FZPAK2 = `grep -a '2nd generation packing quality :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Ramachandran plot appearance   :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FZRAMA = 'NA'
  else
    set FZRAMA = `grep -a 'Ramachandran plot appearance   :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FCHI12 = 'NA'
  else
    set FCHI12 = `grep -a 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Backbone conformation          :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FBCONF = 'NA'
  else
    set FBCONF = `grep -a 'Backbone conformation          :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond lengths                   :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FBRMSZ = 'NA'
  else
    set FBRMSZ = `grep -a 'Bond lengths                   :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond angles                    :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FARMSZ = 'NA'
  else
    set FARMSZ = `grep -a 'Bond angles                    :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  #Get bump scores
  if (`grep -a -c 'Total number of bumps:' $WORKDIR/wf/pdbout.txt` == 0) then
    set FBUMPS = 0
    set FSBMPL = 0
    set FWBMPS = 0.000
  else
    set FBUMPS = `grep -a 'Total number of bumps:' $WORKDIR/wf/pdbout.txt | cut -c 24-28`
    set FSBMPL = `grep -a 'Total squared bump value:' $WORKDIR/wf/pdbout.txt | cut -c 27-33`
    set FWBMPS = `echo $FNATOM $FSBMPL | awk '{printf "%5.3f\n", 100*$2/$1}'`
  endif

  #Count unsatisfied donors and acceptors
  if (`grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set FHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set FHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  if (`grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set FHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set FHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  @ FHBUNS = ($FHBDON + $FHBACC)
  if (`grep -a -c 'Buried donors:' $WORKDIR/wf/pdbout.txt` == 0) then
    set FHBSAT = 'NA'
  else
    @ FBHBDA = (`grep -a 'Buried donors:' $WORKDIR/wf/pdbout.txt | awk '{print $3}'` + `grep -a 'Buried acceptors:' $WORKDIR/wf/pdbout.txt | awk '{print $3}'`)
    if ($FBHBDA == 0) then
      set FHBSAT = 'NA'
    else
      set FHBSA1 = `grep -a 'with a H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $5} END {print SUM}'`
      set FHBSA2 = `grep -a 'with a poor H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $6} END {print 0.5*SUM}'`
      set FHBSA3 = `grep -a 'with only a very poor H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $8} END {print 0.25*SUM}'`
      set FHBSA4 = `grep -a 'essentially without H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $5} END {print 0.125*SUM}'`
      set FHBSAT = `echo $FHBSA1 $FHBSA2 $FHBSA3 $FHBSA4 $FBHBDA | awk '{printf "%5.3f\n" ,($1 + $2 + $3 + $4)/$5}'`
    endif  
  endif

else
  set FZPAK1 = 'NA'
  set FZPAK2 = 'NA'
  set FZRAMA = 'NA'
  set FCHI12 = 'NA'
  set FBCONF = 'NA'
  set FBRMSZ = 'NA'
  set FARMSZ = 'NA'
  set FBUMPS = 'NA'
  set FWBMPS = 'NA'
  set FHBUNS = 'NA'
  set FHBSAT = 'NA'
endif

#Get statistics from tortoize
if (-e $WORKDIR/${PDBID}_final_tortoize.json) then 
  set FZRAMA  =  `jq '."model"."1"."ramachandran-z"' $WORKDIR/${PDBID}_final_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set FSZRAMA =  `jq '."model"."1"."ramachandran-jackknife-sd"' $WORKDIR/${PDBID}_final_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set FCHI12  =  `jq '."model"."1"."torsion-z"' $WORKDIR/${PDBID}_final_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
  set FSCHI12 =  `jq '."model"."1"."torsion-jackknife-sd"' $WORKDIR/${PDBID}_final_tortoize.json | awk '{if ($1 == "null") {print "NA"} else {printf ("%.3f\n", $1)}}'`
else
  #Scream bloody murder
  echo "COMMENT: tortoize cannot validate final model" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                               >> $WHYNOT
  exit(1)
endif


#Run distel if needed (note that the restraints may have changed so the rmsZ values may be different)
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.distel.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json

  #Consolidate the distance restraints
  if (-e $WORKDIR/homology.rest) then
    cat $WORKDIR/homology.rest > $WORKDIR/allhb.rest
  endif
  if (-e $WORKDIR/hbond.rest) then
    cat $WORKDIR/hbond.rest >> $WORKDIR/allhb.rest
  endif
  
  #Calculate rmsZ values. 
  set OHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  set NHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  set FHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_final.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  rm $WORKDIR/allhb.rest
else 
  set OHRMSZ = 'NA'
  set NHRMSZ = 'NA'
  set FHRMSZ = 'NA'
endif

#Nucleic acid validation
set FBPHBRMSZ = 'NA'
set FSHEAR    = 'NA'
set FSTRETCH  = 'NA'
set FBUCKLE   = 'NA'
set FPROPEL   = 'NA'
set FDNRMSD   = 'NA'
set FCONFAL   = 'NA'
set TFCONFAL  = 'NA'
set FBPGRMSZ  = 'NA'

if ($GOT_NUC == 'T') then
  echo "-Validating nucleic acids" | tee -a $LOG

  #Get the basepair values if any were detected
  if (-e $WORKDIR/bphbond.rest) then
    if (`grep -c 'exte dist' $WORKDIR/bphbond.rest` > 0) then
      echo " o Validating base pairs" | tee -a $LOG
      set OBPHBRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/bphbond.rest | tail -n 1 | awk '{print $4}'`
      set NBPHBRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/bphbond.rest | tail -n 1 | awk '{print $4}'`
      set FBPHBRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_final.pdb $WORKDIR/bphbond.rest | tail -n 1 | awk '{print $4}'`
    
      #Run nucrmsz again and get out the results
      python3 $TOOLS/nucrmsz.py $WORKDIR/${PDBID}_final.pdb $TOOLS/x3dna-dssr > $WORKDIR/nucrmsz.log
      set FSHEAR   = `grep 'Shear rmsZ'     $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set FSTRETCH = `grep 'Stretch rmsZ'   $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set FBUCKLE  = `grep 'Buckle rmsZ'    $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set FPROPEL  = `grep 'Propeller rmsZ' $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
      set FBPGRMSZ = `grep 'bpG rmsZ'       $WORKDIR/nucrmsz.log | cut -d ':' -f 2 | awk '{print $1}'`
    endif
  endif
    
  #Validate dinucleotides with DNATCO
  echo " o Validating dinucleotides with DNATCO" | tee -a $LOG
  
  $TOOLS/dnatco.py $WORKDIR/${PDBID}_final.pdb >>& $WORKDIR/dnatco.log
  if (-e ${PDBID}_final.pdb.dnatco.json.gz) then 
    set FDNRMSD  = `zcat ${PDBID}_final.pdb.dnatco.json.gz | jq .overall."average_rmsd" | tr -d '"'`
    set FCONFAL  = `zcat ${PDBID}_final.pdb.dnatco.json.gz | jq .overall."confal_score" | tr -d '"'`
    set TFCONFAL = `zcat ${PDBID}_final.pdb.dnatco.json.gz | jq .overall."confal_percentile" | tr -d '"'`
  else
    #Check what went wrong in DNATCO
    if (`grep -c "doesn't contain enough DNA/RNA steps" $WORKDIR/${PDBID}_final.pdb.html` > 0) then
      echo "   * Not enough dinucleotides to use DNATCO" | tee -a $LOG
    else
      echo "   * Unknown DNATCO error" | tee -a $LOG
      echo "COMMENT: DNATCO cannot validate final model" >> $DEBUG
      echo "PDB-REDO,$PDBID"                             >> $DEBUG
    endif  
  endif
endif

#Calculate percentiles vs the PDB
if ($FZRAMA == 'NA') then
  set TFZRAMA = 'NA'
else
  set TFZRAMA = `cat $TOOLS/zrama.sort | awk -v VAL=$FZRAMA '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OZRAMA == 'NA') then
  set TOZRAMA = 'NA'
else
  set TOZRAMA = `cat $TOOLS/zrama.sort | awk -v VAL=$OZRAMA '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FCHI12 == 'NA') then
  set TFCHI12 = 'NA'
else
  set TFCHI12 = `cat $TOOLS/chi12.sort | awk -v VAL=$FCHI12 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OCHI12 == 'NA') then
  set TOCHI12 = 'NA'
else
  set TOCHI12 = `cat $TOOLS/chi12.sort | awk -v VAL=$OCHI12 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FWBMPS == 'NA') then
  set TFWBMPS = 'NA'
else
  set TFWBMPS = `cat $TOOLS/bumps.sort | awk -v VAL=$FWBMPS '{if ($1 > VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OWBMPS == 'NA') then
  set TOWBMPS = 'NA'
else
  set TOWBMPS = `cat $TOOLS/bumps.sort | awk -v VAL=$OWBMPS '{if ($1 > VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FZPAK2 == 'NA') then
  set TFZPAK2 = 'NA'
else
  set TFZPAK2 = `cat $TOOLS/zpak2.sort | awk -v VAL=$FZPAK2 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OZPAK2 == 'NA') then
  set TOZPAK2 = 'NA'
else
  set TOZPAK2 = `cat $TOOLS/zpak2.sort | awk -v VAL=$OZPAK2 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FZPAK1 == 'NA') then
  set TFZPAK1 = 'NA'
else
  set TFZPAK1 = `cat $TOOLS/zpak1.sort | awk -v VAL=$FZPAK1 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OZPAK1 == 'NA') then
  set TOZPAK1 = 'NA'
else
  set TOZPAK1 = `cat $TOOLS/zpak1.sort | awk -v VAL=$OZPAK1 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FHBSAT == 'NA') then
  set TFHBSAT = 'NA'
else
  set TFHBSAT = `cat $TOOLS/hbsat.sort | awk -v VAL=$FHBSAT '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OHBSAT == 'NA') then
  set TOHBSAT = 'NA'
else
  set TOHBSAT = `cat $TOOLS/hbsat.sort | awk -v VAL=$OHBSAT '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FBPGRMSZ == 'NA') then
  set TFBPGRMSZ = 'NA'
else
  set TFBPGRMSZ = `cat $TOOLS/rmszbpg.sort | awk -v VAL=$FBPGRMSZ '{if ($1 > VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OBPGRMSZ == 'NA') then
  set TOBPGRMSZ = 'NA'
else
  set TOBPGRMSZ = `cat $TOOLS/rmszbpg.sort | awk -v VAL=$OBPGRMSZ '{if ($1 > VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Model validation results ******" | tee -a $LOG
echo '                                     Before  Re-ref  Final'      | tee -a $LOG
echo "1st generation packing quality     : $OZPAK1 $NZPAK1 $FZPAK1"    | tee -a $LOG
echo "2nd generation packing quality     : $OZPAK2 $NZPAK2 $FZPAK2"    | tee -a $LOG
if ($GOT_PROT == 'T') then
  echo "Ramachandran plot Z-score          : $OZRAMA $NZRAMA $FZRAMA"    | tee -a $LOG
  echo "Ramachandran plot RMSD             : $OSZRAMA $NSZRAMA $FSZRAMA" | tee -a $LOG
  echo "chi-1/chi-2 rotamer Z-score        : $OCHI12 $NCHI12 $FCHI12"    | tee -a $LOG
  echo "chi-1/chi-2 rotamer RMSD           : $OSCHI12 $NSCHI12 $FSCHI12" | tee -a $LOG
  echo "Backbone conformation              : $OBCONF $NBCONF $FBCONF"    | tee -a $LOG
endif  
echo " " | tee -a $LOG
echo "Bond length RMS Z-score            : $OBRMSZ $NBRMSZ $FBRMSZ"    | tee -a $LOG
echo "Bond angle RMS Z-score             : $OARMSZ $NARMSZ $FARMSZ"    | tee -a $LOG
echo " " | tee -a $LOG
echo "Total number of bumps              : $OBUMPS $NBUMPS $FBUMPS"    | tee -a $LOG
echo "Weighted bump severity score       : $OWBMPS $NWBMPS $FWBMPS"    | tee -a $LOG
echo " " | tee -a $LOG
echo "Unsatisfied H-bond donors/acceptors: $OHBUNS $NHBUNS $FHBUNS"    | tee -a $LOG
echo "H-bond satisfaction fraction       : $OHBSAT $NHBSAT $FHBSAT"    | tee -a $LOG
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  echo "H-bond restraint RMS Z-score       : $OHRMSZ $NHRMSZ $FHRMSZ"  | tee -a $LOG
endif
if ($GOT_NUC == 'T') then
  echo " " | tee -a $LOG
  if (-e $WORKDIR/bphbond.rest) then
    if (`grep -c 'exte dist' $WORKDIR/bphbond.rest` > 0) then
      echo "Basepair H-bond RMS Z-score        : $OBPHBRMSZ $NBPHBRMSZ $FBPHBRMSZ" | tee -a $LOG
      echo "Basepair shear RMS Z-score         : $OSHEAR $NSHEAR $FSHEAR"          | tee -a $LOG
      echo "Basepair stretch RMS Z-score       : $OSTRETCH $NSTRETCH $FSTRETCH"    | tee -a $LOG
      echo "Basepair buckle RMS Z-score        : $OBUCKLE $NBUCKLE $FBUCKLE"       | tee -a $LOG
      echo "Basepair propeller RMS Z-score     : $OPROPEL $NPROPEL $FPROPEL"       | tee -a $LOG
      echo "Basepair geometry RMS Z-score      : $OBPGRMSZ $NBPGRMSZ $FBPGRMSZ"    | tee -a $LOG 
    endif
  endif  
  echo "Dinucleotide CONFAL score          : $OCONFAL  $NCONFAL  $FCONFAL"       | tee -a $LOG
  echo "Dinucleotide conformation rmsd     : $ODNRMSD $NDNRMSD $FDNRMSD"         | tee -a $LOG
endif

# if ($NMETALREST2 > 0) then
#   echo "Metal restraint RMS Z-score        : $OMRMSZ $NMRMSZ $FMRMSZ" | tee -a $LOG
# endif
echo " " | tee -a $LOG


#Give the side chain details
echo " " | tee -a $LOG
echo "****** Protein structure changes ******" | tee -a $LOG

#Set fallback values
set NDROTA = 'NA'
set HBFLIP = 'NA'
set STFLIP = 'NA'

#Get model changes
if (-e $WORKDIR/renumber.json) then
  set NUMBERING = "-renum $WORKDIR/renumber.json"
endif
  
echo "-Analysing model changes" | tee -a $LOG
$TOOLS/modelcompare -v \
-pdb1 $WORKDIR/cache.pdb \
-pdb2 $WORKDIR/${PDBID}_final.pdb \
-output-name $PDBID \
$NUMBERING \
-output-dir $WORKDIR > $WORKDIR/modelcompare.log
if ($status || ! -e $WORKDIR/${PDBID}_coot_tour.py) then
  #Give error message
  echo " o Problem with modelcompare" | tee -a $LOG
  echo "COMMENT: modelcompare: error in model comparison" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                  >> $DEBUG
  
  #Fill in the missing values
  set HBFLIP  = 'NA'
  set NDROTA  = 'NA'
  set CTFLIP  = 'NA'
  set TCFLIP  = 'NA'
  set PEPFIX  = 'NA'
  set PEPDIS  = 'NA'
  set NWATDEL = 'NA'
else  
  #Rename the COOT scripts
  mv $WORKDIR/${PDBID}_coot_tour.py  $WORKDIR/${PDBID}_final.py
  mv $WORKDIR/${PDBID}_coot_tour.scm $WORKDIR/${PDBID}_final.scm
  
  #Count things
  set HBFLIP  = `grep 'H-bond flips marked' $WORKDIR/modelcompare.log | awk '{print $1}'`
  set NDROTA  = `grep 'changed rotamers marked' $WORKDIR/modelcompare.log | awk '{print $1}'`
  set CTFLIP  = `grep 'cis-trans isomerisations marked' $WORKDIR/modelcompare.log | awk '{print $1}'`
  set TCFLIP  = `grep 'trans-cis isomerisations marked' $WORKDIR/modelcompare.log | awk '{print $1}'`
  set PEPFIX  = `grep 'removed distorted peptides marked' $WORKDIR/modelcompare.log | awk '{print $1}'`
  set PEPDIS  = `grep 'introduced distorted peptides marked' $WORKDIR/modelcompare.log | awk '{print $1}'`
  set NWATDEL = `grep 'removed waters marked' $WORKDIR/modelcompare.log | awk '{print $1}'`
endif



#Give the summary
echo " " | tee -a $LOG
echo "Changed rotamers     : $NDROTA" | tee -a $LOG
echo "Hydrogen bond flips  : $HBFLIP" | tee -a $LOG
echo "Cis-trans flips      : $CTFLIP" | tee -a $LOG
echo "Trans-cis flips      : $TCFLIP" | tee -a $LOG
echo "Peptides fixed       : $PEPFIX" | tee -a $LOG
echo "Peptides distorted   : $PEPDIS" | tee -a $LOG
echo " " | tee -a $LOG


#Give the density fit details
echo " " | tee -a $LOG
echo "****** Electron density map details ******" | tee -a $LOG

#Set fallback values
set RSRB  = 'NA'
set RSRW  = 'NA'
set RSCCB = 'NA'
set RSCCW = 'NA'

echo "-Analysing density maps" | tee -a $LOG

#Generate total B-factors (the 0-cycle PDB file already has total B-factors)
if ($DOTLS == 0) then
  cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final_tot.pdb
else
  tlsanl \
  XYZIN $WORKDIR/${PDBID}_final.pdb \
  XYZOUT $WORKDIR/${PDBID}_final_tot.pdb \
<<eof > $WORKDIR/mapval.log
    BINPUT t
    BRESID t
    ISOOUT FULL
    NUMERICAL
    END
eof
  if ($status) then
    if (`grep '^[AH][TE]' ${PDBID}_final.pdb | grep -c -E '\*{6}'` != 0) then
      echo " o Problem with the atomic B-factor overflow" | tee -a $LOG
      if ($LOCAL == 0) then
        echo "COMMENT: B-factor value overflow" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                  >> $WHYNOT
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
        exit(1)
      else
        #Use the residual B-factors (not pretty)
        cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final_tot.pdb
      endif
    else
      #TLSANL failed
      echo " o Problem with TLSANL" | tee -a $LOG
      echo "COMMENT: TLSANL cannot calculate total B-factors" >> $DEBUG
      echo "PDB-REDO,$PDBID"                                  >> $DEBUG

      #Use the residual B-factors (not pretty)
      cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final_tot.pdb
    endif  
  endif
endif

#Run stats on the original model
$TOOLS/stats \
--use-auth-ids \
--sampling-rate 1.5 \
--hklin $WORKDIR/${PDBID}_0cyc.mtz \
$SFTYPE \
--xyzin $WORKDIR/${PDBID}_0cyc.pdb \
-o $WORKDIR/${PDBID}_0cyc.json \
--output-format json \
$DICTCMD >& $WORKDIR/density_0cyc.log 

#Run stats on the final model
$TOOLS/stats \
--use-auth-ids \
--sampling-rate 1.5 \
--hklin $WORKDIR/${PDBID}_final.mtz \
$SFTYPE \
--xyzin $WORKDIR/${PDBID}_final.pdb \
-o $WORKDIR/${PDBID}_final.json \
--output-format json \
$DICTCMD >& $WORKDIR/density_final.log 

#Check for weird things
@ NANTOT = (`grep -c 'nan' $WORKDIR/${PDBID}_0cyc.json` + `grep -c 'nan' $WORKDIR/${PDBID}_final.json`)
if ($NANTOT > 12 && $SERVER == 0) then
  #Only stop if not in server mode
  echo "COMMENT: suspicious output from stats" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                       >> $WHYNOT
  exit(1)
endif
    
#Analyse changes in density fit
$TOOLS/dRSCC -v \
-eds1 $WORKDIR/${PDBID}_0cyc.json \
-eds2 $WORKDIR/${PDBID}_final.json \
-output-name mapval >>& $WORKDIR/mapval.log

if ($status) then
  #Error making the figures
  echo " o Problem analysing real-space validation data" | tee -a $LOG
  if ($SERVER == 0) then
    #Only stop if not in server mode
    echo "COMMENT: error in real-space validation" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                         >> $WHYNOT
    exit(1)
  endif
else

  #Get the significant changes
  set RSRB  = `grep 'Better RSR'   $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`
  set RSRW  = `grep 'Worse RSR'    $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`
  set RSCCB = `grep 'Better RSCC'  $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`
  set RSCCW = `grep 'Worse RSCC'   $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`
  set EDIAB = `grep 'Better EDIAm' $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`
  set EDIAW = `grep 'Worse EDIAm'  $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`
  set OPIAB = `grep 'Better OPIA'  $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`
  set OPIAW = `grep 'Worse OPIA'   $WORKDIR/mapval_eds_compare.txt | awk '{print $3}'`

  #Give the summary
  echo " " | tee -a $LOG
  echo '                       Better Worse'  | tee -a $LOG
  echo "Real-space CC        : $RSCCB $RSCCW" | tee -a $LOG
  echo "Real-space R-factor  : $RSRB $RSRW"   | tee -a $LOG
  echo "EDIAm density fit    : $EDIAB $EDIAW" | tee -a $LOG
  echo "OPIA density coverage: $OPIAB $OPIAW" | tee -a $LOG
  echo " " | tee -a $LOG

endif

################################################    Ligand validation    #################################################

#Validate existing ligands
echo " " | tee -a $LOG
echo "****** Ligand validation ******" | tee -a $LOG

#Make sure there are existing ligands
if (`echo $LIG_LIST | grep -c ':'` != 0) then

  echo "-Analysing existing ligands" | tee -a $LOG

  #Create YASARA validation files in parallel
  if ($?YASARA) then
    #Loop over ligands
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.YASARA.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    foreach LIG (`echo $LIG_LIST | sed 's/:/ /g'`)

      #Specify the residue
      set RESNUM = `echo $LIG | cut -c 2-`
      set CHID   = `echo $LIG | cut -c 1`
      set PDBLIG = `echo "$CHID $RESNUM" | awk '{printf "%s%4d\n", $1, $2}' | sed 's/ /_/g'`
      
      #Only run for existing ligands
      if (`grep '^[AH][TE][OT]' $WORKDIR/${PDBID}_final.pdb | cut -c 18-27 | sed 's/ /_/g' | grep -c $PDBLIG` > 0) then

        #Return label for job launching
ligvalrunning:

        #Only launch new jobs when the number of cores is not exceeded
        #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
        jobs > $WORKDIR/jobs.log

        if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

          echo " o Validating residue $LIG" | tee -a $LOG

          $YASARA -txt $TOOLS/ligval.mcr \
          "oripdb='$WORKDIR/${PDBID}_0cyc.pdb'" \
          "newpdb='$WORKDIR/${PDBID}_final.pdb'" \
          "resnum='$RESNUM'" \
          "chid='$CHID'" > $WORKDIR/ligval_$LIG.log
        else
          #Wait a bit to start again
          sleep 2
          goto ligvalrunning
        endif
      endif  
    end

    #Wait for the jobs to finish
    wait
  endif
endif

#Validate new ligands
#Make sure there are new ligands
if (`echo $NLIG_LIST | grep -c ':'` != 0) then

  echo "-Analysing new ligands" | tee -a $LOG

  #Create YASARA validation files in parallel
  if ($?YASARA) then
    #Loop over ligands
    cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.YASARA.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json
    
    foreach LIG (`echo $NLIG_LIST | sed 's/:/ /g'`)

      #Specify the residue
      set RESNUM = `echo $LIG | cut -c 2-`
      set CHID   = `echo $LIG | cut -c 1`

      #Return label for job launching
ligvalrunning:

      #Only launch new jobs when the number of cores is not exceeded
      #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
      jobs > $WORKDIR/jobs.log

      if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

        echo " o Validating residue $LIG" | tee -a $LOG

        $YASARA -txt $TOOLS/ligval_solo.mcr \
        "newpdb='$WORKDIR/${PDBID}_final.pdb'" \
        "resnum='$RESNUM'" \
        "chid='$CHID'" > $WORKDIR/ligval_$LIG.log
      else
        #Wait a bit to start again
        sleep 2
        goto ligvalrunning
      endif
    end

    #Wait for the jobs to finish
    wait
  endif
endif

#Accumulate the results if there are any ligands
if (`echo $LIG_LIST | grep -c ':'` != 0 || `echo $NLIG_LIST | grep -c ':'` != 0) then
  echo "-Compiling validation results" | tee -a $LOG
  python3 $TOOLS/ligval_json.py $WORKDIR > $WORKDIR/ligvaljson.log
endif

#Report the validation results

#Loop over existing ligands
foreach LIG (`echo $LIG_LIST | sed 's/:/ /g'`)

  #Specify the residue
  set TRESNUM = `echo $LIG | cut -c 2-`
  if (`echo $TRESNUM | grep -c "[[:alpha:]]"` != 0) then
    set RESNUM = `echo -n $TRESNUM | head -c -1`
    set INS    = `echo -n $TRESNUM | tail -c 1`
  else  
    set RESNUM = $TRESNUM
    set INS    = ""
  endif
  set CHID   = `echo $LIG | cut -c 1`
  set PDBLIG = `echo "$CHID $RESNUM" | awk '{printf "%s%4d\n", $1, $2}' | sed 's/ /_/g'`
    
  #Only run for existing ligands
  if (`grep '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_final.pdb | cut -c 18-27 | sed 's/ /_/g' | grep -c $PDBLIG$INS` > 0) then
    set RESID  = `grep '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_final.pdb | cut -c 18-27 | sed 's/ /_/g' | grep $PDBLIG$INS | head -n 1 | cut -c 1-3 | sed 's/_/ /g'`

    #Get real-space values
    set LIGRSRO  = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.density_fit.real_space_Rfactor" | awk '{printf "%5.3f", $1}'`
    set LIGRSRF  = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.real_space_Rfactor" | awk '{printf "%5.3f", $1}'` 
    set LIGCCO   = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.density_fit.real_space_correlation" | awk '{printf "%5.3f", $1}'`
    set LIGCCF   = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.real_space_correlation" | awk '{printf "%5.3f", $1}'` 
    set LIGEDIAO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.density_fit.EDIAm_density_fit" | awk '{printf "%5.3f", $1}'`
    set LIGEDIAF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.EDIAm_density_fit" | awk '{printf "%5.3f", $1}'` 
    set LIGOPIAO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.density_fit.OPIA_density_coverage" | awk '{printf "%5.1f", $1}'`
    set LIGOPIAF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.OPIA_density_coverage" | awk '{printf "%5.1f", $1}'`     
    
    #Extract YASARA results
    #Shifts
    set NSHIFT = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | ."\""atoms_shifted_more_than_0.5A"\""" | awk '{printf "%.0f", $1}'`
    set RMSD   = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .all_atom_rmsd_in_A" | awk '{printf "%5.3f", $1}'`
    
    #Heat of formation
    set EFORMO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.heat_of_formation.energy" | sed 's/null/NA/'`
    set EFORMF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.heat_of_formation.energy" | sed 's/null/NA/'`  
    #Make pretty
    if ($EFORMO != NA) then
      set EFORMO = `echo $EFORMO | awk '{printf "%5.1f", $1}'`
    endif  
    if ($EFORMF != NA) then
      set EFORMF = `echo $EFORMF | awk '{printf "%5.1f", $1}'`
    endif    

    #Hydrogen bonds
    set EHBOO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.H_bonds_strength" | awk '{printf "%5.1f", $1}'`
    set EHBOF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.H_bonds_strength" | awk '{printf "%5.1f", $1}'`  
    set NHBOO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.H_bonds_count" | awk '{printf "%.0f", $1}'`
    set NHBOF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.H_bonds_count" | awk '{printf "%.0f", $1}'`  

    #Bumps
    set NBUMPO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.bumps_count" | awk '{printf "%.0f", $1}'`
    set NBUMPF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.bumps_count" | awk '{printf "%.0f", $1}'`  

    #Hydrophobic contacts and strength
    set SHYPHO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.hydrophobic_strength" | awk '{printf "%5.3f", $1}'`
    set SHYPHF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.hydrophobic_strength" | awk '{printf "%5.3f", $1}'`  
    set NHYPHO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.hydrophobic_count" | awk '{printf "%.0f", $1}'`
    set NHYPHF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.hydrophobic_count" | awk '{printf "%.0f", $1}'`  

    #pi-pi contacts and strength
    set SPIPIO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.pi_pi_strength" | awk '{printf "%5.3f", $1}'`
    set SPIPIF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.pi_pi_strength" | awk '{printf "%5.3f", $1}'`  
    set NPIPIO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.pi_pi_count" | awk '{printf "%.0f", $1}'`
    set NPIPIF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.pi_pi_count" | awk '{printf "%.0f", $1}'`  

    #Cation-pi contacts and strength
    set SCATPO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.cation_pi_strength" | awk '{printf "%5.3f", $1}'`
    set SCATPF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.cation_pi_strength" | awk '{printf "%5.3f", $1}'`  
    set NCATPO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.cation_pi_count" | awk '{printf "%.0f", $1}'`
    set NCATPF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.cation_pi_count" | awk '{printf "%.0f", $1}'` 
  
    #Give summary
    echo " "                                                              | tee -a $LOG 
    echo "****** Ligand validation details ($RESID $CHID $RESNUM$INS) ******" | tee -a $LOG
    echo " "                                                              | tee -a $LOG 
    echo "                                    Before  Final"          | tee -a $LOG 
    echo "Real-space R-factor               : $LIGRSRO   $LIGRSRF"    | tee -a $LOG 
    echo "Real-space correlation            : $LIGCCO   $LIGCCF"      | tee -a $LOG
    echo "EDIAm density fit                 : $LIGEDIAO   $LIGEDIAF"  | tee -a $LOG
    echo "OPIA density coverage (%)         : $LIGOPIAO    $LIGOPIAF" | tee -a $LOG   
    echo " "                                                          | tee -a $LOG 
    if ($?YASARA) then
      echo "Energy of formation (kJ/mol)      : $EFORMO $EFORMF"   | tee -a $LOG 
      echo "Hydrogen bond energy (kJ/mol)     : $EHBOO  $EHBOF"    | tee -a $LOG 
      echo "Number of hydrogen bonds          : $NHBOO  $NHBOF"    | tee -a $LOG 
      echo "Number of bumps                   : $NBUMPO  $NBUMPF"  | tee -a $LOG 
      echo "Number of hydrophobic interactions: $NHYPHO  $NHYPHF"  | tee -a $LOG 
      echo "Hydrophobic interaction strength  : $SHYPHO  $SHYPHF"  | tee -a $LOG 
      echo "Number of Pi-Pi interactions      : $NPIPIO  $NPIPIF"  | tee -a $LOG 
      echo "Pi-Pi interaction strength        : $SPIPIO  $SPIPIF"  | tee -a $LOG 
      echo "Number of cation-Pi interactions  : $NCATPO  $NCATPF"  | tee -a $LOG 
      echo "Cation-Pi interaction strength    : $SCATPO  $SCATPF"  | tee -a $LOG 
      echo " "                                                     | tee -a $LOG 
      echo "Atoms shifted more than ${SHIFTCO}A     : $NSHIFT"     | tee -a $LOG 
      echo "All atom RMSD for residue (A)     : ${RMSD}"           | tee -a $LOG 
      echo " "                                                     | tee -a $LOG 
    endif
  endif  
end


#Report over new ligands
foreach LIG (`echo $NLIG_LIST | sed 's/:/ /g'`)

  #Specify the residue
  set TRESNUM = `echo $LIG | cut -c 2-`
  if (`echo $TRESNUM | grep -c "[[:alpha:]]"` != 0) then
    set RESNUM = `echo -n $TRESNUM | head -c -1`
    set INS    = `echo -n $TRESNUM | tail -c 1`
  else  
    set RESNUM = $TRESNUM
    set INS    = ""
  endif
  set CHID   = `echo $LIG | cut -c 1`
  set PDBLIG = `echo "$CHID $RESNUM" | awk '{printf "%s%4d\n", $1, $2}' | sed 's/ /_/g'`
  set RESID  = `grep '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_final.pdb | cut -c 18-27 | sed 's/ /_/g' | grep $PDBLIG | head -n 1 | cut -c 1-3 | sed 's/_/ /g'`

  #Get real-space values
  set LIGRSRF  = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.real_space_Rfactor" | awk '{printf "%5.3f", $1}'` 
  set LIGCCF   = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.real_space_correlation" | awk '{printf "%5.3f", $1}'` 
  set LIGEDIAF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.EDIAm_density_fit" | awk '{printf "%5.3f", $1}'` 
  set LIGOPIAF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.density_fit.OPIA_density_coverage" | awk '{printf "%5.1f", $1}'`     
    
  #Get occupancy
  set LIGOCC = `grep '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_final.pdb | grep $RESID | grep $CHID | grep -e "$RESNUM" | head -n 1 | cut -c 57-60` 

  #Extract YASARA results
  #Heat of formation
  set EFORMF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.heat_of_formation.energy" | sed 's/null/NA/'`  
  
  #Make pretty
  if ($EFORMF != NA) then
    set EFORMF = `echo $EFORMF | awk '{printf "%5.1f", $1}'`
  endif    

  #Hydrogen bonds
  set EHBOF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.H_bonds_strength" | awk '{printf "%5.1f", $1}'`  
  set NHBOF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.H_bonds_count" | awk '{printf "%.0f", $1}'`  

  #Bumps
  set NBUMPO = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .original_model.interactions.bumps_count" | awk '{printf "%.0f", $1}'`

  #Hydrophobic contacts and strength
  set SHYPHF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.hydrophobic_strength" | awk '{printf "%5.3f", $1}'`  
  set NHYPHF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.hydrophobic_count" | awk '{printf "%.0f", $1}'`  

  #pi-pi contacts and strength
  set SPIPIF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.pi_pi_strength" | awk '{printf "%5.3f", $1}'`  
  set NPIPIF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.pi_pi_count" | awk '{printf "%.0f", $1}'`  

  #Cation-pi contacts and strength
  set SCATPF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.cation_pi_strength" | awk '{printf "%5.3f", $1}'`  
  set NCATPF = `cat ${PDBID}_ligval.json | jq ".Ligand_validation_data | .[] | select(.pdb.seqNum == $RESNUM and .pdb.insCode == "\""$INS"\"" and .pdb.compID == "\""$RESID"\"" and .pdb.strandID == "\""$CHID"\"") | .pdb_redo_model.interactions.cation_pi_count" | awk '{printf "%.0f", $1}'` 
 
  #Give summary
  echo " "                                                       | tee -a $LOG
  echo "****** New ligand details ($RESID $CHID $RESNUM) ******" | tee -a $LOG
  echo " "                                                       | tee -a $LOG 
  echo "                                    Score"               | tee -a $LOG
  echo "Real-space R-factor               : $LIGRSRF"            | tee -a $LOG
  echo "Real-space correlation            : $LIGCCF"             | tee -a $LOG
  echo "EDIAm density fit                 : $LIGEDIAF"           | tee -a $LOG
  echo "OPIA density coverage (%)         : $LIGOPIAF"           | tee -a $LOG
  echo "Refined occupancy                 : $LIGOCC"             | tee -a $LOG
  echo " "                                                       | tee -a $LOG 
  if ($?YASARA) then
    echo "Energy of formation (kJ/mol)      : $EFORMF"           | tee -a $LOG 
    echo "Hydrogen bond energy (kJ/mol)     : $EHBOF"            | tee -a $LOG 
    echo "Number of hydrogen bonds          : $NHBOF"            | tee -a $LOG 
    echo "Number of bumps                   : $NBUMPF"           | tee -a $LOG 
    echo "Number of hydrophobic interactions: $NHYPHF"           | tee -a $LOG 
    echo "Hydrophobic interaction strength  : $SHYPHF"           | tee -a $LOG 
    echo "Number of Pi-Pi interactions      : $NPIPIF"           | tee -a $LOG 
    echo "Pi-Pi interaction strength        : $SPIPIF"           | tee -a $LOG 
    echo "Number of cation-Pi interactions  : $NCATPF"           | tee -a $LOG 
    echo "Cation-Pi interaction strength    : $SCATPF"           | tee -a $LOG
    echo " "                                                     | tee -a $LOG
  endif
end

################################################   Wrap up  #################################################

echo " " | tee -a $LOG
echo "****** Creating final output ******" | tee -a $LOG


#Create Refmac command file 

#Not for the data bank
if ($LOCAL == 1 || $SERVER == 1) then
  echo "-Creating Refmac command script" | tee -a $LOG

  #Copy all optimised settings
  echo "#Refmac command script from PDB-REDO $VERSION" > $WORKDIR/$PDBID.refmac
  echo "#"                                            >> $WORKDIR/$PDBID.refmac
  echo "#Use of riding hydrogens"                     >> $WORKDIR/$PDBID.refmac
  echo "make hydrogen $HYDROGEN"                      >> $WORKDIR/$PDBID.refmac
  echo "#B-factor model selection"                    >> $WORKDIR/$PDBID.refmac
  echo "refi bref $BREFTYPE"                          >> $WORKDIR/$PDBID.refmac
  if (`echo $REFIRES | grep -c REFI` != 0) then
    echo "#Resolution cutoffs"                        >> $WORKDIR/$PDBID.refmac
    echo "$REFIRES"                                   >> $WORKDIR/$PDBID.refmac
  endif
  echo "#Solvent related settings"                    >> $WORKDIR/$PDBID.refmac
  echo "scal type $SOLVENT $SCALING"                  >> $WORKDIR/$PDBID.refmac
  echo "solvent YES"                                  >> $WORKDIR/$PDBID.refmac
  echo "$MASKPAR"                                     >> $WORKDIR/$PDBID.refmac
  echo "tlsd waters exclude"                          >> $WORKDIR/$PDBID.refmac
  echo "#Restraint weights"                           >> $WORKDIR/$PDBID.refmac
  #Keep the final refinement settings unless there was no final refinement
  if ($RESOTYPE > 3) then
    if ($TLSBEST == "none") then
      set CWEIGHT = `echo "$WGTSIG AUTO 2.50"`
    else
      set CWEIGHT = `echo "$WGTSIG MATRIX $TLSBEST"`
    endif
  else if ($BLTBEST == 'auto') then
    set CWEIGHT = `echo "$WGTSIG AUTO 2.50"`
  else
    set CWEIGHT = `echo "$WGTSIG $WGTTYPE $BLTBEST"`
  endif
  echo "weight  $CWEIGHT"                             >> $WORKDIR/$PDBID.refmac
  #Only give the B-factor restraint weight if optimised
  if ($BBEST == "overall" || $BBEST == "none") then
    #Do nothing
  else
    echo "temp $BBEST"                                >> $WORKDIR/$PDBID.refmac
  endif
  if ($TWIN == "twin") then
    echo "#Twinning"                                  >> $WORKDIR/$PDBID.refmac
    echo "$TWIN"                                      >> $WORKDIR/$PDBID.refmac
  endif
  if (`echo $NCSTYPE | grep -c ncs` != 0) then
    echo "#NCS handling"                              >> $WORKDIR/$PDBID.refmac
    echo "$NCSTYPE"                                   >> $WORKDIR/$PDBID.refmac
    echo "$NCSALIGN"                                  >> $WORKDIR/$PDBID.refmac
    echo "$NCSNEIGH"                                  >> $WORKDIR/$PDBID.refmac
    echo "$NCSSTRICT"                                 >> $WORKDIR/$PDBID.refmac
  endif
  if (`echo $JELLY | grep -c ridg` != 0) then
    echo "#Other restraints"                          >> $WORKDIR/$PDBID.refmac
    echo "$JELLY"                                     >> $WORKDIR/$PDBID.refmac
    echo "$TORSION"                                   >> $WORKDIR/$PDBID.refmac
  endif
endif

###################################################   Clean up round 2   #################################################

#Clean out the mtz files (only for the databank)
if ($LOCAL == 0) then

  #Delete more data?
  set DELLABEL =
  if (`echo $ANOMCOEF | cut -c 1` != "") then
    set DELLABEL = "DELFAN PHDELAN"
  endif

  #Loop over all mtz files
  foreach STAGE (0cyc besttls final)

    #Make a copy
    cp $WORKDIR/${PDBID}_$STAGE.mtz $WORKDIR/allcolumn.mtz

    #Remove the FC PHIC FC_ALL_LS and PHIC_ALL_LS columns
    mtzutils \
    HKLIN  $WORKDIR/allcolumn.mtz \
    HKLOUT $WORKDIR/${PDBID}_$STAGE.mtz \
<<eof >> $WORKDIR/mtzcleanup.log
      EXCLUDE FC PHIC FC_ALL_LS PHIC_ALL_LS $DELLABEL
      END
eof
  end
endif

#################################################    Copy over files    ###################################################
echo "-Consolidating results" | tee -a $LOG


#Create directories and copy files
gzip -f $WORKDIR/${PDBID}_0cyc.pdb
gzip -f $WORKDIR/${PDBID}_0cyc.mtz
gzip -f $WORKDIR/${PDBID}_0cyc.json
gzip -f $WORKDIR/${PDBID}_besttls.pdb
gzip -f $WORKDIR/${PDBID}_besttls.mtz

#Create mmcif file

#Merge data from input model
cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software."cif-merge".used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json

$TOOLS/cif-merge -v \
-i $WORKDIR/${PDBID}_final_tot.pdb \
-o $WORKDIR/${PDBID}_final.cif \
$DICTCMD \
--donor $WORKDIR/cache.pdb >& $WORKDIR/cif-merge.log

#Fallback direct conversion
if (! -e $WORKDIR/${PDBID}_final.cif) then
  cp $WORKDIR/versions.json $WORKDIR/versions.json.bak && jq '.software.pdb2cif.used |= true' $WORKDIR/versions.json.bak > $WORKDIR/versions.json

  $TOOLS/pdb2cif -v \
  $WORKDIR/${PDBID}_final_tot.pdb \
  $WORKDIR/${PDBID}_final.cif \
  $DICTCMD >>& $WORKDIR/cif-merge.log
endif

if (! -e $WORKDIR/${PDBID}_final.cif) then
  #Give error message
  echo " o Problem with mmCIF conversion" | tee -a $LOG
  echo "COMMENT: cannot make mmCIF file" >> $DEBUG
  echo "PDB-REDO,$PDBID"                 >> $DEBUG
endif  


set TOUTPUT = "$OUTPUT.tmp"
mkdir -p $TOUTPUT
cp $WORKDIR/${PDBID}_0cyc.pdb.gz          $TOUTPUT/
cp $WORKDIR/${PDBID}_0cyc.mtz.gz          $TOUTPUT/
cp $WORKDIR/${PDBID}_0cyc.json.gz         $TOUTPUT/
cp $WORKDIR/${PDBID}_besttls.pdb.gz       $TOUTPUT/
cp $WORKDIR/${PDBID}_besttls.mtz.gz       $TOUTPUT/
cp $WORKDIR/${PDBID}_final.pdb            $TOUTPUT/
cp $WORKDIR/${PDBID}_final.cif            $TOUTPUT/
cp $WORKDIR/${PDBID}_final_tot.pdb        $TOUTPUT/
cp $WORKDIR/${PDBID}_final.mtz            $TOUTPUT/
cp $WORKDIR/${PDBID}_final.json           $TOUTPUT/
if (-e $WORKDIR/${PDBID}_final.scm) then
  cp $WORKDIR/${PDBID}_final.scm          $TOUTPUT/
  cp $WORKDIR/${PDBID}_final.py           $TOUTPUT/
endif  
cp $WORKDIR/versions.json                 $TOUTPUT/
if (-e $WORKDIR/${PDBID}_ligval.json) then
  cp $WORKDIR/${PDBID}_ligval.json        $TOUTPUT/
endif
if (-e $WORKDIR/${PDBID}.rtest) then
  cp $WORKDIR/${PDBID}.rtest              $TOUTPUT/
endif
if (-e $WORKDIR/${PDBID}_het.cif) then
  cp $WORKDIR/${PDBID}_het.cif            $TOUTPUT/
  if ($LOCAL == 0) then
    mkdir -p $RESTOUT
    cp $WORKDIR/${PDBID}_het.cif          $RESTOUT/
  endif
endif
if (-e $WORKDIR/renumber.json) then
  cp $WORKDIR/renumber.json               $TOUTPUT/
endif
if (-e $WORKDIR/${PDBID}_final.dssp) then
  cp $WORKDIR/${PDBID}_final.dssp         $TOUTPUT/
endif

if ($WCERR == 0) then
  mkdir -p $TOUTPUT/wo
  mkdir -p $TOUTPUT/wc
  mkdir -p $TOUTPUT/wf
  cp $WORKDIR/wo/* $TOUTPUT/wo/
  cp $WORKDIR/wc/* $TOUTPUT/wc/
  cp $WORKDIR/wf/* $TOUTPUT/wf/
endif


######################################### Create data files for PDBe sliders ##############################################

#Calculate the quality change scores
set DIV = 3

#Fit to the data
if ($SIGRFFIN == 'NA') then
    set SIGRFU = $SIGRFCAL
else
    set SIGRFU = $SIGRFFIN
endif
if ($TSTCNT < 1 || $GOTR == 0 || $ZCALERR == 1) then
    set ZDFREE = `echo "$RFCALUNB $RFFIN $SIGRFU" | awk '{print ($1 - $2)/$3}'`
else
    set ZDFREE = `echo "$RFCAL $RFFIN $SIGRFU" | awk '{print ($1 - $2)/$3}'`
endif

#Start the JSON file
echo '{'                                 > $WORKDIR/pdbe.json
echo "  "\""pdbid"\"": "\""$PDBID"\""," >> $WORKDIR/pdbe.json
echo '  "ddatafit": {'                  >> $WORKDIR/pdbe.json
echo "    "\""zdfree"\"": $ZDFREE,"     >> $WORKDIR/pdbe.json
echo '    "range-lower": -12.9,'        >> $WORKDIR/pdbe.json
echo '    "range-upper": 12.9'          >> $WORKDIR/pdbe.json
echo -n '    }'                         >> $WORKDIR/pdbe.json #Allow for a comma on this line

#Add protein scores if there is protein
if ($GOT_PROT == 'T') then
  #Geometric quality
  if ($OZPAK2 == 'NA' || $FZPAK2 == 'NA') then
    set DZPAK2 = 0.000
    @ DIV = ($DIV - 1)
  else
    set DZPAK2 = `echo "$OZPAK2 $FZPAK2" | awk '{print ($2 - $1)}'`
  endif
  if ($OZRAMA == 'NA' || $FZRAMA == 'NA') then
    set DZRAMA = 0.000
    @ DIV = ($DIV - 1)
  else
    set DZRAMA = `echo "$OZRAMA $FZRAMA" | awk '{print ($2 - $1)}'`
  endif
  if ($OCHI12 == 'NA' || $FCHI12 == 'NA') then
    set DCHI12 = 0.000
    @ DIV = ($DIV - 1)
  else
    set DCHI12 = `echo "$OCHI12 $FCHI12" | awk '{print ($2 - $1)}'`
  endif

  if ($DIV == 0) then
    set DZSCORE = 'null'
  else
    set DZSCORE = `echo "$DZPAK2 $DZRAMA $DCHI12 $DIV" | awk '{print (($1 + $2 + $3)/$4)}'`
  endif
  
  #Update the JSON file (not using jq as it is annoying with dashes)
  echo ','                                >> $WORKDIR/pdbe.json
  echo '  "geometry": {'                  >> $WORKDIR/pdbe.json
  echo "    "\""dzscore"\"": $DZSCORE,"   >> $WORKDIR/pdbe.json
  echo '    "range-lower": -1.17,'        >> $WORKDIR/pdbe.json
  echo '    "range-upper": 1.17'          >> $WORKDIR/pdbe.json
  echo -n '    }'                         >> $WORKDIR/pdbe.json #Allow for a comma on this line
endif  

#Add nucleic acid scores if there is nucleic acid
if ($GOT_NUC == 'T') then
  if ($FBPGRMSZ == 'NA' || $OBPGRMSZ == 'NA') then
    set DBPGRMSZ = 'null'
  else  
    set DBPGRMSZ = `echo "$OBPGRMSZ $FBPGRMSZ" | awk '{print ($1 - $2)}'` 
  endif
  
  #Update the JSON file (not using jq as it is annoying with dashes)
  echo ','                              >> $WORKDIR/pdbe.json
  echo '  "base-pairs": {'              >> $WORKDIR/pdbe.json
  echo "    "\""drmsz"\"": $DBPGRMSZ,"  >> $WORKDIR/pdbe.json
  echo '    "range-lower": -2.496,'     >> $WORKDIR/pdbe.json
  echo '    "range-upper": 2.496'       >> $WORKDIR/pdbe.json
  echo -n '    }'                       >> $WORKDIR/pdbe.json #Allow for a comma on this line
endif  
  

#Close the JSON file
echo ' ' >> $WORKDIR/pdbe.json
echo '}' >> $WORKDIR/pdbe.json


#Copy over the file
cp $WORKDIR/pdbe.json $TOUTPUT/

############################################## Final administrative steps ################################################
#Create summary statistics
echo -n "$PDBID $VERSION $RFACT $RFREE $RCAL $RFCAL $SIGRFCAL $RFCALUNB $RFCALZ $RTLS $RFTLS $SIGRFTLS $RFTLSUNB $RFTLSZ $RFIN $RFFIN $SIGRFFIN $RFFINUNB $RFFINZ $RFRRAT $BBEST $TLSBEST $BLTBEST " > $WORKDIR/data.txt
echo -n "$NWATDEL $NBBFLIP $NSCBLT $NDROTA $HBFLIP $STFLIP $NCHIRFX $OZPAK1 $NZPAK1 $FZPAK1 $OZPAK2 $NZPAK2 $FZPAK2 $OZRAMA $NZRAMA $FZRAMA $OCHI12 $NCHI12 $FCHI12 $OBCONF $NBCONF $FBCONF $OBRMSZ $NBRMSZ $FBRMSZ $OARMSZ $NARMSZ $FARMSZ " >> $WORKDIR/data.txt
echo -n "$OBUMPS $NBUMPS $FBUMPS $OHBUNS $NHBUNS $FHBUNS $OGFOLD $NGFOLD $FGFOLD $PROG $DYEAR $RESOLUTION $DATARESH $DATARESL $NREFCNT $TSTCNT $TSTPRC $NTSTCNT $REFPATM $AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA $BAVER $BWILS $BREFTYPE " >> $WORKDIR/data.txt
echo -n "$SOLVENT $VDWPROBE $IONPROBE $RSHRINK $DOTLS $NTLS $OPTTLSG $ORITLS $LEGACY '$SPACEGROUP' $RSCCB $RSCCW $RSRB $RSRW $OWBMPS $NWBMPS $FWBMPS $OHBSAT $NHBSAT $FHBSAT $URESO $CCWOLD $CCWFIN $ZCCW $CCFOLD $CCFFIN $ZCCF " >> $WORKDIR/data.txt
echo -n "$WAVELENGTH $ISTWIN $SOLVD $EXPTYP $COMPLETED $NOPDB $NOSF $USIGMA $ZCALERR $TIME $RESOTYPE $FALSETWIN $TOZRAMA $TFZRAMA $TOCHI12 $TFCHI12 $TOZPAK2 $TFZPAK2 $TOWBMPS $TFWBMPS $TOHBSAT $TFHBSAT $OHRMSZ $NHRMSZ $FHRMSZ " >> $WORKDIR/data.txt
echo -n "$TOZPAK1 $TFZPAK1 $NLOOPS $NMETALREST2 $OSZRAMA $NSZRAMA $FSZRAMA $OSCHI12 $NSCHI12 $FSCHI12 $SRFRRAT $ZRFRRATCAL $ZRFRRATTLS $ZRFRRATFIN $GOT_PROT $GOT_NUC " >> $WORKDIR/data.txt
echo -n "$NNUCLEICREST $OBPHBRMSZ $NBPHBRMSZ $FBPHBRMSZ $OSHEAR $NSHEAR $FSHEAR $OSTRETCH $NSTRETCH $FSTRETCH $OBUCKLE $NBUCKLE $FBUCKLE $OPROPEL $NPROPEL $FPROPEL $OCONFAL $NCONFAL $FCONFAL $TOCONFAL $TNCONFAL $TFCONFAL " >> $WORKDIR/data.txt
echo    "$ODNRMSD $NDNRMSD $FDNRMSD $OBPGRMSZ $NBPGRMSZ $FBPGRMSZ $TOBPGRMSZ $TFBPGRMSZ" >> $WORKDIR/data.txt
cp $WORKDIR/data.txt $TOUTPUT/
if ($?COMMENT) then
  python3 $TOOLS/txt2json.py -i $WORKDIR/data.txt -o $TOUTPUT/tdata.json -c "$COMMENT"
else
  python3 $TOOLS/txt2json.py -i $WORKDIR/data.txt -o $TOUTPUT/tdata.json
endif

#Add the title to data.json
cat $TOUTPUT/tdata.json | jq --arg title "$TITLE" '.properties += {TITLE: $title}' > $TOUTPUT/data.json
rm $TOUTPUT/tdata.json
# #Keep a long-term copies for FAIRness
# if (-e $OUTPUT/old_versions) then
#   cp -r $OUTPUT/old_versions $TOUTPUT/
# else
#   mkdir -p $TOUTPUT/old_versions
# endif  
# gzip -c $WORKDIR/${PDBID}_final.cif > $TOUTPUT/old_versions/${PDBID}_v${VERSION}_final.cif.gz
# gzip -c $WORKDIR/versions.txt       > $TOUTPUT/old_versions/${PDBID}_v${VERSION}_versions.txt.gz
# gzip -c $TOUTPUT/data.json          > $TOUTPUT/old_versions/${PDBID}_v${VERSION}_data.json.gz

#Give final results (usefull for batch jobs and quick result interpretation)
if ($SERVER == 1) then
  cat $WORKDIR/data.txt >> $BASE/success.txt
endif

#Finish the logfiles and copy them over
if ($LOCAL == 1 || $SERVER == 1) then
  #Clean and copy the log file
  grep -v '^\[' $LOG > $TOUTPUT/$PDBID.log
endif


#Copy files over for running jobs in CCP4i or deposition
if ($LOCAL == 1 || $SERVER == 1) then
  cp $WORKDIR/$PDBID.refmac      $TOUTPUT/
  cp $WORKDIR/${PDBID}_final.log $TOUTPUT/
  cp $WORKDIR/*.rest             $TOUTPUT/ >& /dev/null
  cp $WORKDIR/resolute.log       $TOUTPUT/ >& /dev/null
  if ($DOTLS == 1) then
    cp $WORKDIR/${PDBID}_final.tls $TOUTPUT/
  endif
endif


#Swap out directories
if (-e $OUTPUT) then
  mv $OUTPUT $OUTPUT.old 
  mv $TOUTPUT $OUTPUT  
  if ($SERVER == 0 && $LOCAL == 1) then
    #Do not remove old directory
  else
    if ($SERVER == 1) then
      cp $OUTPUT.old/process.log $OUTPUT
    endif  
    rm -r $OUTPUT.old
  endif
else
  mv $TOUTPUT $OUTPUT
endif

#Finish things for the server
if ($SERVER == 1) then
  #Compress the lot and add a link
  cd $STDIR/output
  mv data.txt  $STDIR/data.txt

  #Write out status file
  touch $STDIR/processEnded.txt
endif

#Remove tempdir
if ($CLEAN == 1) then
  rm -rf $WORKDIR
endif

echo "-Done. That was a lot of work!" | tee -a $LOG
