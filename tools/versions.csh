#!/bin/tcsh -f

# versions.csh: tries to find the versions of all PDB_REDO related programs
# Call:
# ./versions.csh TOOLS_Directory
#
# Changelog:
# - Added workaround for all the temporary file from WHAT_CHECK
# - Added support for not having FoldX.
# - Reduced the number of variables by recycling $VERSION.
# - Stopped listing programs in order of appearance. They are now showed in 
#   order of addition.
# - Now working in json.

#Set up the locations of the programs 
set TOOLS = $1
set VFILE = $2
set PDBID = temp
set D2    = XX
set TMPFIL = `mktemp`
source $TOOLS/pdb_redo.setup

#Programs
#  1) mtzdmp, always
set VERSION = `mtzdump -i | grep Program | awk '{print $4}'`
if ($VERSION == "") then
  set VERSION = `mtzdump -i | grep 'patch' | awk '{print $6}'`
endif
cp $VFILE $VFILE.tmp && jq '.software +={"mtzdmp":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

#  2) kollumer, optional
set VERSION = `$TOOLS/kollumer | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"kollumer":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

#  3) mtz2various, always
set VERSION = `mtz2various -i | grep Program | awk '{print $4}'`
if ($VERSION == "") then
  set VERSION = `mtz2various -i | grep 'patch' | awk '{print $6}'`
endif
cp $VFILE $VFILE.tmp && jq '.software +={"mtz2various":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

#  4) WHAT_CHECK, always
mkdir wcjunk && cd wcjunk
  $WC/bin/whatcheck \
  << eof >& versions.log
eof
set VERSION = `grep version versions.log | cut -c 34-37`
cp $VFILE $VFILE.tmp && jq '.software +={"WHAT_CHECK":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE
cd .. && rm -r wcjunk

#  6) prepper, always
set VERSION = `$TOOLS/prepper --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"prepper":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE


#  8) cif2cif, always
set VERSION = `$TOOLS/cif2cif | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"cif2cif":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

#  9) extractor, always
set VERSION = `$TOOLS/extractor | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"extractor":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 10) cif2mtz, always
set VERSION = `cif2mtz -i | grep Program | awk '{print $3}'`
if ($VERSION == "") then
  set VERSION = `cif2mtz -i | grep 'patch' | awk '{print $6}'`
endif
cp $VFILE $VFILE.tmp && jq '.software +={"cif2mtz":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 11) mtzutils, always
set VERSION = `mtzutils -i | grep Program | awk '{print $3}'`
if ($VERSION == "") then
  set VERSION = `mtzutils -i | grep 'patch' | awk '{print $6}'`
endif
cp $VFILE $VFILE.tmp && jq '.software +={"mtzutils":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 12) cad, always
set VERSION = `cad -i | grep Program | awk '{print $3}'`
if ($VERSION == "") then
  set VERSION = `cad -i | grep 'patch' | awk '{print $6}'`
endif
cp $VFILE $VFILE.tmp && jq '.software +={"cad":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 13) ctruncate, optional 
set VERSION = `ctruncate -i | grep 'Program:' | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"ctruncate":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 14) unique, always
set VERSION = `unique -i | grep Program | awk '{print $3}'`
if ($VERSION == "") then
  set VERSION = `unique -i | grep 'patch' | awk '{print $6}'`
endif
cp $VFILE $VFILE.tmp && jq '.software +={"unique":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 15) freerflag, always
set VERSION = `freerflag -i | grep Program | awk '{print $4}'`
if ($VERSION == "") then
  set VERSION = `freerflag -i | grep 'patch' | awk '{print $6}'`
endif
cp $VFILE $VFILE.tmp && jq '.software +={"freerflag":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 16) sfcheck, always
set VERSION = `sfcheck -h | grep Vers | cut -d ';' -f 1 | awk '{print $6}'`
cp $VFILE $VFILE.tmp && jq '.software +={"sfcheck":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 17) refmac, always
set VERSION = `refmac5 -i | grep Program | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"refmac":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 18) fitr, always
set VERSION = `$TOOLS/fitr | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"fitr":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 19) chiron, optional
set VERSION = `$TOOLS/chiron | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"chiron":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 20) libg, optional
cp $VFILE $VFILE.tmp && jq '.software +={"libg":{ "version":null, "used":false }}' $VFILE.tmp > $VFILE

# 21) binliner, optional
set VERSION = `$TOOLS/binliner | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"binliner":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 22) resolute, optional
set VERSION = `$TOOLS/resolute | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"resolute":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 23) bselect, optional
set VERSION = `$TOOLS/bselect | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"bselect":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 24) picker, always
set VERSION = `$TOOLS/picker | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"picker":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 25) tlsanl, optional
tlsanl \
<< eof >& versions.log
eof
set VERSION = `grep VERSION versions.log | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"tlsanl":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 26) libcheck, optional
set VERSION = `libcheck -i | grep Vers | awk '{print $5}'`
cp $VFILE $VFILE.tmp && jq '.software +={"libcheck":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 27) centrifuge, optional
set VERSION = `$TOOLS/centrifuge --version | grep Version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"centrifuge":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 28) DSSP, always
set VERSION = `mkdssp --version | head -n 1| awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"DSSP":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 29) pepflip, optional
set VERSION = `$TOOLS/pepflip --version | grep Version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"pepflip":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 30) mini-rsr, optional (with pepflip, or ligand fitting)
set VERSION = `coot-mini-rsr --version | grep version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"coot-mini-rsr":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 31) SideAide, optional
set VERSION = `$TOOLS/SideAide --version | grep Version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"SideAide":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 32) what_todo, optional
set VERSION = `$TOOLS/what_todo | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"what_todo":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 33) longinus, optional
set VERSION = `$TOOLS/longinus | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"longinus":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 34) YASARA, optional
if ($?YASARA) then
  $YASARA -txt \
  << eof >& versions.log
  Exit
eof
  set VERSION = `grep Version versions.log | awk '{print $2}'`
  cp $VFILE $VFILE.tmp && jq '.software +={"YASARA":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE
endif

# 36) stats, always
set VERSION = `$TOOLS/stats --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"stats":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 38) sftools, optional
cp $VFILE $VFILE.tmp && jq '.software +={"sftools":{ "version":null, "used":false }}' $VFILE.tmp > $VFILE

# 39) platonyzer, optional
set VERSION = `$TOOLS/platonyzer --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"platonyzer":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 40) detectHbonds, optional
set VERSION = `$TOOLS/detectHbonds | grep version | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"detectHbonds":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 41) modelcompare, always
set VERSION = `$TOOLS/modelcompare | grep version | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"modelcompare":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 42) hoder, optional
set VERSION = `$TOOLS/hoder | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"HODER":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 43) PHASER, optional
set VERSION = `phaser --version | head -n 1 | cut -c 8-  | awk '{print $1}'`
cp $VFILE $VFILE.tmp && jq '.software +={"PHASER":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 44) txt2json.py, always
set VERSION = `$TOOLS/txt2json.py -h |  grep version | awk '{print $3}' | tr -d ')'`
cp $VFILE $VFILE.tmp && jq '.software +={"txt2json":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 45) seqrescopier, optional
set VERSION = `$TOOLS/seqrescopier | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"seqrescopier":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 46) flipper, always
set VERSION = `$TOOLS/flipper --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"flipper":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 47) distel, optional
set VERSION = `$TOOLS/distel.py --version | awk '{print $3}' | sed 's/)//'`
cp $VFILE $VFILE.tmp && jq '.software +={"distel":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 48) rnbterror, optional
set VERSION = `$TOOLS/rnbterror | grep version | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"rnbterror":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 49) rwcontents, always
set VERSION = `rwcontents -i | grep 'patch' | awk '{print $6}'`
cp $VFILE $VFILE.tmp && jq '.software +={"rwcontents":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 50) loopwhole, optional
set VERSION = `$TOOLS/loopwhole | grep version | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"loopwhole":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 51) loopwhole-validate
set VERSION = `$TOOLS/loopwhole-validate | grep version | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"loopwhole-validate":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 52) findligand, optional 
cp $VFILE $VFILE.tmp && jq '.software +={"findligand":{ "version":null, "used":false }}' $VFILE.tmp > $VFILE

# 53) cif-merge, always
set VERSION = `$TOOLS/cif-merge --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"cif-merge":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 54) fixDMC, optional
set VERSION = `$TOOLS/fixDMC | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"fixDMC":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 55) pdb2cif
set VERSION = `$TOOLS/pdb2cif --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"pdb2cif":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 56) carbonanza, always
set VERSION = `$TOOLS/carbonanza --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"carbonanza":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 57) pdb2fasta
set VERSION = `$TOOLS/pdb2fasta | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"pdb2fasta":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 58) blastp 
if ($?BLASTP) then
  set VERSION = `$BLASTP -version | grep 'blastp' | awk '{print $2}'`
  cp $VFILE $VFILE.tmp && jq '.software +={"BLASTp":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE
endif

# 59) carbivore
set VERSION = `$TOOLS/carbivore | head -n 1 | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"carbivore":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 61) pdbset
set VERSION = `pdbset -i | grep 'patch' | awk '{print $6}'`
cp $VFILE $VFILE.tmp && jq '.software +={"pdbset":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 62) mmCQL
set VERSION = `$TOOLS/mmCQL --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"mmCQL":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 63) tortoize
set VERSION = `$TOOLS/tortoize --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"tortoize":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 64) molrep
set VERSION = `molrep -h | grep Vers | awk '{print $3}' | tr -d ';'`
cp $VFILE $VFILE.tmp && jq '.software +={"molrep":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 65) coot proper
set VERSION = `coot --version | head -n 1`
cp $VFILE $VFILE.tmp && jq '.software +={"COOT":{ "version":$version, "used":false }}' --arg version "$VERSION" $VFILE.tmp > $VFILE

# 66) cif-grep
set VERSION = `$TOOLS/cif-grep --version | awk '{print $3}'`
cp $VFILE $VFILE.tmp && jq '.software +={"cif-grep":{ "version":$version, "used":false }}' --arg version $VERSION $VFILE.tmp > $VFILE

# 67) dRSCC, always
set VERSION = `$TOOLS/dRSCC | grep version | awk '{print $4}'`
cp $VFILE $VFILE.tmp && jq '.software +={"dRSCC":{ "version":$version, "used":true }}' --arg version $VERSION $VFILE.tmp > $VFILE

#Cleanup
rm versions.log
rm sfcheck.log
rm $VFILE.tmp

