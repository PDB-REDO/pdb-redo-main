      PROGRAM EXTRACTOR
C=======================================================================
C  Version 6.03 2023-01-31
C  Extracts a number of stats out of a PDB and a SF file.
C
C  Usage: extractor (flags) PDB_IN CIF_IN MON_LIST REDO_DATA DATA_OUT 
C         (TLS_OUT)
C
C  A TLS file in REFMAC5 format is created if TLS_OUT is specified. TLS
C  parameters are extracted from the PDB's REMARK cards
C  If the -f flag is specified a TLS file will be produced even if there
C  is no TLS information in the file header. In this case 1 TLS group 
C  per chain is created.
C  In '-v (verbose)' mode resolutions are printed for every reflection.
C
C  Warning: this program assumes that the mmCIF files are in the format
C  recommended bij the PDB. If in doubt, preparse the mmCIF file with 
C  cif2cif
C
C  Written by Robbie Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  EXTRACTOR. Please, refer to (one of) these references:
C  - Robbie P. Joosten, Gert Vriend: "PDB improvement starts with data 
C    deposition" Science, 317, p. 195-196 (2007)
C  - Robbie P. Joosten, Thomas Womack, Gert Vriend and Gerard Bricogne:
C    "Re-refinement fromdeposited X-ray data can deliver improved models
C    for most PDB entries"  Acta Cryst. D65, p. 176-185 (2009)
C  - Robbie P. Joosten, Jean Salzemann, Vincent Bloch, Heinz Stockinger,
C    Ann-Charlott Berglund, Christophe Blanchet, Erik Bongcam-Rudloff, 
C    Christophe Combet, Ana L. Da Costa, Gilbert Deleage, Matteo 
C    Diarena, Roberto Fabbretti, Geraldine Fettahi, Volker Flegel, 
C    Andreas Gisel, Vinod Kasam, Timo Kervinen, Eija Korpelainen, Kimmo
C    Mattila, Marco Pagni, Matthieu Reichstadt, Vincent Breton, Ian J. 
C    Tickle, Gert Vriend: "PDB_REDO: automated re-refinement of X-ray
C    structure models in the PDB" J. Appl. Cryst., 42, p. 376-384 (2009)
C  - Robbie P. Joosten, Tim A.H. te Beek, Elmar Krieger, Maarten 
C    Hekkelman, Rob W.W. Hooft, Reinhard Schneider, Chris Sander, Gert 
C    Vriend: "A series of PDB related databases for everyday needs" 
C    Nucl. Acids Res., 39, p. D411-D419 (2011)
C  - Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend, 
C    Anastassis Perrakis: "Automatic rebuilding and optimization of 
C    crystallographic structures in the Protein Data Bank" 
C    Bioinformatics, 27, p. 3392-3398 (2011)
C  - Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastassis
C    Perrakis: "PDB_REDO: constructive validation, more than just 
C    looking for errors" Acta Cryst. D68, p. 484-496 (2012)
C
C  This program contains code from SFTOOLS:
C    Copyright (C) 1999 Bart Hazes
C
C Change log
C Version 6:
C 6.04 Added fix for dealing with _refln.scale_group_code
C 6.03 Added fixes for quotes in program names
C 6.02 Title records are now also mined
C 6.01 Some cleanup and tollerance for more compounds
C 6.00 The cell dimensions are now taken from the reflection data rather
C      than from the model.
C
C Version 5:
C 5.09 Fix for long Phenix type TLS selections.
C 5.08 Fix for REFMAC TLS groups that comprise all atoms.
C 5.07 An attempt is made to read the wavelength from REMARK 200.
C 5.06 The mon_lib_list.cif file is now read more like mmCIF.
C 5.05 Peptide or amino acid ligands in the same chain are now  
C      added to the parent TLS group in one-group-per-chain TLS.
C 5.04 Fixed bugs in TLS group extraction.
C 5.03 Fixed a bug in the TLS group appending of compounds.Fixed bugs in
C      TLS group extraction.
C 5.02 Stopped writing out empty TLS files if no groups were created
C 5.01 Added support for EM structure models
C 5.00 Residues that make more than one LINK can be added to a TLS group
C      if they are on an 'approved' list.
C
C Version 4:
C 4.23 Improved reading of multiline TLS group selections.
C 4.22 Added support for more chains.
C 4.21 Improved reading of multiline TLS group selections.
C 4.20 Checks whether there are nucleic acid residues. Improvement to
C      TLSNOT routine.
C 4.19 Refmac style TLS group rages are now read in free format.
C 4.18 Format fix for number of reflections. Fix for reading TLS groups
C      in Phenix format.
C 4.17 Added some padding behind the TLS selections to improve parsing.
C 4.16 Added an error routine for corrupted atom records.
C 4.15 Added explicit error message for a missing space group.
C      Added a workaround for ligands withou a CHAINid. 
C 4.14 Made the MINCLR subroutine more general purpose. Bugfix for 
C      Phenix-style TLS groups.
C 4.13 Rewrote the routine that writes exclusion lists to allow negative
C      residue numbers. 
C 4.12 More TLS parsing updates. TLS groups that only mention the chain 
C      now run from residue -10 to 9999 (instead of 1 to 9999).
C 4.11 Moved code that parses negative selections to the subroutine
C      TLSNOT. Negative selctors for Hydrogens and waters are removed,
C      other negative selectors stop the parseing of TLS selections. 
C 4.10 Extended TLS parsing logic. Removed a stupid Buster-related bug.
C 4.09 Added another format for reading the R-factor
C 4.08 TLS groups created from the coordinates will be written to the 
C      file 'REDO.tls'.
C 4.07 TLS groups with negative selectors are ignored. Modified amino 
C      acids are now also taken into account for TLS groups. Added
C      another format for reading the R-factor.
C 4.06 Fixed a bug in the creation of a list for ligand validation. 
C 4.05 Fixed a bug to deal with residue names shorter than 3 characters.
C 4.04 Modified nucleic acid residues are now considered when making TLS
C      group selections. This required a substantial rewrite of the atom
C      parser.
C 4.03 Fixed the Buster format parser.
C 4.02 Increased the read-ahead for phenix-type TLS selections.
C 4.02 Added another fix for phenix-type selections.
C 4.01 Bugfix in the way solvent mask parameters are read.
C 4.00 All unconnected ligands are now detected and added to a list for 
C      (future) ligand rebuilding. Potential ligands are selected from 
C      the NON-POLYMERs in CCP4's mon_lib_list.cif. Ligands to be 
C      excluded are extracted from pdb_redo.dat.
C
C Version 3:
C 3.13 MASION is now set to the same value as MASSHR. I.e. The ion probe
C      for the solvent mask is set to the shrinkage is no other value is
C      found.19 This is similar to the Refamc default.
C 3.12 Fixed a bug in the way the selection is passed to the PHENIX-type
C      parser.
C 3.11 Split up reading of BUSTER-type TLS groups: clean groups directly
C      from BUSTER are read neatly, groups from the PDB as junk.
C 3.11 'NULL' selections are now intercepted.
C 3.10 Fixed a bug in the PHENIX-type TLS group parser
C 3.09 Fixed a bug for reading the ion probe radius
C 3.08 Updated the TLS group extraction for another Buster format.
C 3.07 TLS groups should now contain at least 5 residues.
C 3.06 Added a basic help function which starts if no command line 
C      arguments are given.
C 3.06 Now supports longer program names. 
C 3.05 Now tries to extract the solvent mask parameters.
C 3.04 Waters are now ignored when creating TLS ranges.
C 3.03 Bugfix for phenix TLS groups spanning more than one line.
C 3.02 TLS groups are now only generated for chains that contain protein
C      or nucleic acid.
C 3.01 Substantial rewrite for the TLS group extraction from phenix.
C 3.00 When the force keyword is specified a simple TLS group definition
C      is written out even if there is one in the header. 
C
C Version 2:
C 2.14 Fixed the TLS group extraction for files coming directly from 
C      Phenix.refine.
C 2.14 Added support for comma-separated TLS ranges from Buster.
C 2.13 Fixed the TLS group extraction.
C 2.13 If all else fails, the data resolution is now taken in relaxed 
C      mode.
C 2.12 The data completeness is now extracted.
C 2.12 Deals with TLS groups consisting of all atoms by ignoring them.
C 2.11 Added support for TLS group selections from BUSTER.
C 2.11 Changed the output format for the TLS files. To solve reading 
C      problems for the 'ORIGIN' cards in Refmac.
C 2.11 Changed the output format for the TLS files. To solve reading 
C      problems for extreme tensor values in Refmac.
C 2.11 Fixed problem with implicit ranges in PHENIX.
C 2.10 Fixed the PHENIX group parser to deal with cases like in 3kdv.
C 2.09 Added support for new PHENIX TLS group definitions
C 2.08 Writes out a list of residues that should not be refitted because
C      the backbone atoms are involved in LINKs. Trivial links are not
C      taken into account.
C 2.07 Increased the maximum number of residues in the 'keep' lists. 
C 2.07 Added the -r flag which makes extractor slightly more relaxed 
C      about missing values.
C 2.06 Writes out a list of waters that should be kept because they are
C      involved in LINKs.
C 2.05 Writes out a list of residues that should not be refitted because
C      the side chains are involved in LINKs.
C 2.04 Increased the maximum number of chains to 60. PDBid: 3hf9
C 2.03 Changed the CHOP function.
C 2.02 Changed bug in the resolution extraction caused by a change in  
C      the format of the REMARK 2 records.
C 2.01 When the Wilson B-factor is negative, it is reset to 0.00 (that
C      is 'not reported/NULL'). PDBid: 1usz
C======================================================================C
      IMPLICIT NONE
C-----Declare basic variables and parameters
      INTEGER   MAXDAT, MAXCIF, MAXCHN, REFBIN, I, J, K, STATUS, ARGS,
     +          MAXRES, MAXHET, MAXLNK
      CHARACTER RESIDUES*83, BACKBONE*24, VERS*4
      PARAMETER (VERS='6.04')
C-----MAXDAT is the maximum allowed lines in a PDB file
      PARAMETER (MAXDAT=999999)
C-----MAXCIF is the maximum allowed lines in a CIF file
      PARAMETER (MAXCIF=11000000)
C-----MAXCHN is the maximum allowed chains in a PDB file
      PARAMETER (MAXCHN=99)
C-----The number of reflection bins used for calculation
      PARAMETER (REFBIN=20)
C-----The number of residues in the NO_BLD arrays
      PARAMETER (MAXRES=2000)
C-----The maximum number of hetero compounds
      PARAMETER (MAXHET=40000)
C-----The maximum number of cached LINKs
      PARAMETER (MAXLNK=2000)


C-----Declare the variables and parameters
      CHARACTER LINE*80,   CHOP*80,    NOSPCE*40,  RTEMP*8,  SYMMHM*13,
     +          PROG*40,   ARG5*10,    CJUNK*1,    C2JUNK*2,  TMPCHN1*1,
     +          TMPCHN2*1, OLDCHN*1,   TMPINS1*1,  TMPINS2*1, TMPATM1*4,
     +          TMPATM2*4, TMPRESI1*3, TMPRESI2*3, TMPALT*1,  JRES*3
C     Character variables:
C     TMPRESI1  A temporary residue name
C     Files:
C     Name      Handle  Description
C     INPDB     7       PDB file from which data is extracted
C     INCIF     8       Reflection file in mmCIF format
C     OUTDAT    9       Output file for extracted data
C     TLSOUT    10      First TLS output file
C     TLSOUT2   10      Second TLS output file
C     MON_LIST  11      List with all monomers in CCP4's dictionary
C     REDO_DATA 12      Data file from PDB_REDO (pdb_redo.dat)
      CHARACTER INPDB*255,   INCIF*255,    OUTDAT*255, TLSOUT*255, 
     +          TLSOUT2*255, MON_LIST*255, REDO_DATA*255, TITLE*255,
     +          TTITLE*70
      REAL      RESO,   AAXIS, BAXIS,  CAXIS,  ALPHA,  BETA,  GAMMA,
     +          RFACT,  RFREE,  HKLRES, HIRES,  LOWRES, TSTPRC,BTOT,
     +          BAVER,  BTEMP,  DIFR,   GTREAL, WAVEL, GETVAL
      INTEGER   ATMCNT, LABCNT, HASH,   JUNK1,  JUNK2, JUNK3, REFCNT, 
     +          TSTCNT, 
     +          EXTRA,  CHAINS, YEAR,   ANISOU, HATMCT, FORTLS, GOTTLS,
     +          TMPRESN1, TMPRESN2, NO_BLD_CNT, H2O_KEEP_CNT, OLDRESN,
     +          BBN_KEEP_CNT, BBO_KEEP_CNT, LIG_CNT, LENSTR, NHET, NNUC,
     +          NAMINO, NALL, NTLSRES, NLINK
      LOGICAL   GTRESO, USETLS, GTTEST, VERBOS, RELAX, PROT,
     +          MACRES, NUCLIC, EMMOD
C PROT   The structure contains protein
C NUCLIC The structure contains nucleic acids    

C-----Declare arrays
      CHARACTER CHAINI(MAXCHN)*1,       TSTVAL(MAXCIF)*1, 
     +          NO_BLD_CHN(MAXRES)*1,   NO_BLD_INS(MAXRES)*1, 
     +          H2O_KEEP_CHN(MAXRES)*1, H2O_KEEP_INS(MAXRES)*1,
     +          BBN_KEEP_CHN(MAXRES)*1, BBN_KEEP_INS(MAXRES)*1,
     +          BBO_KEEP_CHN(MAXRES)*1, BBO_KEEP_INS(MAXRES)*1,
     +          LIG_CHN(MAXRES)*1,      LIG_INS(MAXRES)*1,
     +          LINK_CACHE(MAXLNK)*80
C     Lists of compounds to be used
      CHARACTER HETERO(MAXHET)*3, NUCLEIC(MAXHET)*3, AMINOA(MAXHET)*3,
     +          ALLRES(MAXHET)*3, TLSRES(MAXHET)*3
C     Numeric values for reflections
      REAL      REFRES(MAXCIF),   SF(MAXCIF),       SIGSF(MAXCIF), 
     +          RATIO(MAXCIF)
      INTEGER   INDEXH(MAXCIF),   INDEXK(MAXCIF),   INDEXL(MAXCIF),
     +          CFIRST(MAXCHN),   CLAST(MAXCHN),    NO_BLD_RES(MAXRES),
     +          H2O_KEEP_RES(MAXRES), BBN_KEEP_RES(MAXRES),
     +          BBO_KEEP_RES(MAXRES), LIG_RES(MAXRES)
      LOGICAL   MACMOL(MAXCHN)
      DOUBLE PRECISION COEFS(6), COEFFS(4)


C========================= Header data block ==========================C

      COMMON /HEADER/   RFAC_ALL  ,RFAC_OBS  ,RFAC_CUT,
     +                  RFREE_ALL ,RFREE_OBS ,RFREE_CUT,
     +                  BWILS     ,RESO3     ,COMPLET,
     +                  MASPRO    ,MASION    ,MASSHR,
     +                  RPROG
      REAL    RFAC_ALL
      REAL    RFAC_OBS
      REAL    RFAC_CUT
      REAL    RFREE_ALL
      REAL    RFREE_OBS
      REAL    RFREE_CUT
      REAL    BWILS
      REAL    RESO3
      REAL    COMPLET
      REAL    MASPRO
      REAL    MASION
      REAL    MASSHR
      CHARACTER RPROG*40

C========================End of declarations===========================C

C=========================PDB header section===========================C
C-----Start with the help function:
      ARGS=IARGC() 
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   EXTRACTOR version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Extractor is the basic data mining tool for PDB_REDO.'
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '          
      WRITE(6,*)'Usage:'
      WRITE(6,*)'extractor (flags) PDB_IN CIF_IN MON_LIST REDO_DATA D'//
     +          'ATA_OUT TLS_OUT'
      WRITE(6,*)' '  
      WRITE(6,*)'PDB_IN is the PDB file from which data must be mined.'
      WRITE(6,*)'CIF_IN is the reflection file (in mmCIF format).'
      WRITE(6,*)'MON_LIST is CCP4''s mon_lib_list.cif file.'
      WRITE(6,*)'REDO_DATA is PDB_REDO''s pdb_redo.dat file.'
      WRITE(6,*)'DATA_OUT contains the extracted data.'
      WRITE(6,*)'TLS_OUT is an output TLS group definition in Refmac '//
     +          'format.'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode. The resolution for each reflectio'//
     +          'n is printed.'
      WRITE(6,*)'-f Forced TLS. A TLS group definition is always crea'//
     +          'ted, even if no definition was given in the header '//
     +          'of PDB_IN. If PDB_IN contained a TLS group definiti'//
     +          'on, a second TLS file is written out.'
      WRITE(6,*)'-r Relaxed mode. Extractor is more forgiving about '//
     +          'data format and missing values.'
      WRITE(6,*)'-e EM mode. Extractor works does not use reflection '//
     +          'data.'     
      WRITE(6,*)'    '
      WRITE(6,*)'Citing EXTRACTOR:'
      WRITE(6,*)'Joosten et al. ''PDB_REDO: automated re-refinement o'//
     +          'f X-ray structure models in the PDB'' J. Appl. Cryst'//
     +          '., 42, 376-384 (2009)'
        GO TO 999
      END IF

C-----Initialise values
100   EXTRA =0
      RESO  =0.0
      AAXIS =0.0
      BAXIS =0.0
      CAXIS =0.0
      ALPHA =0.0
      BETA  =0.0
      GAMMA =0.0
      DIFR  =0.999
      RFACT =0.999
      RFREE =0.999
      VERBOS=.FALSE.
      GTRESO=.FALSE.
      USETLS=.FALSE.
      RELAX =.FALSE.
      EMMOD =.FALSE.
      PROT  =.FALSE.
      NUCLIC=.FALSE.
      FORTLS=0
      GOTTLS=0
      WAVEL    = 0.0000
      RFAC_ALL = 0.999
      RFAC_CUT = 0.999
      RFAC_OBS = 0.999
      RFREE_CUT = 0.999
      RFREE_OBS = 0.999
      BWILS     = 0.0
      RESO3     = 0.0
      COMPLET   = 0.0
      MASPRO    = 1.20
      MASION    = 0.80
      MASSHR    = 0.80
      RPROG  = '?'
      NO_BLD_CNT   = 0
      H2O_KEEP_CNT = 0
      BBN_KEEP_CNT = 0
      BBO_KEEP_CNT = 0
      LIG_CNT      = 0
      NHET         = 0
      NNUC         = 0
      NAMINO       = 0
      NTLSRES      = 0
      NLINK        = 0
      OLDRESN      = -9999
      OLDCHN       = '?' 
      TITLE        = ''

C-----Really ugly!!!
      RESIDUES = 'GLY,ALA,VAL,THR,SER,CYS,TRP,PHE,TYR,HIS,ILE,LEU,'//
     +           'ASP,ASN,GLU,GLN,MET,ARG,LYS,PRO,MSE'
C      CARB     = 'NAG,BMA,MAN,FUC'
      BACKBONE = ' N  , CA , C  , O  , OXT'
      

C-----Parse debug flags
      DO 10, I=1, 5
        CALL GETARG(1+EXTRA, C2JUNK)
        IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
          WRITE(6,*) 'Verbose mode'
          VERBOS=.TRUE.
          EXTRA = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-f').OR.(C2JUNK.EQ.'-F')) THEN
          WRITE(6,*) 'Creating TLS parameters for 1 group/chain'
          FORTLS=1
          EXTRA = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-r').OR.(C2JUNK.EQ.'-R'))THEN
          WRITE(6,*) 'Relaxed mode' 
          RELAX =.TRUE.
          EXTRA = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-e').OR.(C2JUNK.EQ.'-E'))THEN
          WRITE(6,*) 'EM mode' 
          EMMOD =.TRUE.
          EXTRA = EXTRA+1          
        ELSE IF (C2JUNK(1:1).EQ.'-') THEN
          WRITE(6,*) 'Invalid option: ', C2JUNK 
          WRITE(6,*) 'It will be ignored.'
          EXTRA = EXTRA+1
        END IF
10    CONTINUE

C-----Extract TLS groups?
C-----Count number of command line arguments 6th argument should be 
C-----TLS output file
      IF (ARGS.GE.(6+EXTRA)) THEN
        USETLS=.TRUE.
        WRITE(6,*) 'TLS parameters will be used'
        CALL GETARG((6+EXTRA), TLSOUT)
      END IF

C-----Generate a name for the coordinate-based TLS file
      IF (FORTLS.EQ.1) THEN
        DO 110, I=LENSTR(TLSOUT),1, -1
          IF (TLSOUT(I:I).EQ.'/') THEN
            TLSOUT2 = TLSOUT(1:I)//'REDO.tls'
            GO TO 120
          END IF
          TLSOUT2 = 'REDO.tls'
110     CONTINUE
      END IF


C-----Get input PDB file
120   STATUS=0
      CALL GETARG((1+EXTRA), INPDB)
      OPEN(UNIT=7, FILE=INPDB, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) INPDB, 'PDB file cannot be opened!'
        GO TO 999
      END IF

C-----Get the ligand dictionary
      STATUS=0
      CALL GETARG((3+EXTRA), MON_LIST)
      OPEN(UNIT=11, FILE=MON_LIST, STATUS='OLD', IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) MON_LIST, 'Ligand dictionary cannot be opened!'
        GO TO 999
      END IF

C-----Get the PDB_REDO data file
      STATUS=0
      CALL GETARG((4+EXTRA), REDO_DATA)
      OPEN(UNIT=12, FILE=REDO_DATA, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) REDO_DATA, 'PDB_REDO data cannot be opened!'
        GO TO 999
      END IF

C-----Populate the compound lists
      CALL GETCMP(HETERO, MAXHET, NHET, NUCLEIC, NNUC, AMINOA, 
     +            NAMINO, ALLRES, NALL, TLSRES, NTLSRES) 

C-----Initialise year for relaxed mode
      IF ((RELAX.EQV..TRUE.).OR.(EMMOD.EQV..TRUE.)) THEN
        YEAR = 2161
      END IF

C-----Read PDB file line by line
      DO 200, J=1,MAXDAT
        READ(UNIT=7, FMT=920, END=300) LINE
!        PRINT*, LINE
     
C-------Read deposition year
        IF (LINE(1:6).EQ.'HEADER') THEN
          READ(UNIT=LINE(58:59), FMT=915, ERR=200) YEAR
          IF (YEAR.GT.70) THEN
            YEAR = YEAR + 1900
          ELSE  
            YEAR = YEAR + 2000
          END IF
          
C-------Read the title
        ELSE IF (LINE(1:6).EQ.'TITLE ') THEN
          READ(UNIT=LINE(11:80), FMT='(A70)', ERR=200) TITLE
C         Read ahead to get the full title          
          DO 210, K=1, 20
            READ(UNIT=7, FMT=920, END=300) LINE
            IF (LINE(1:6).EQ.'TITLE ') THEN
              TITLE = TITLE(1:LENSTR(TITLE))// LINE(11:80)
            ELSE
              BACKSPACE(7)
            END IF
210       CONTINUE          
C         Clean up the line
          CALL DUBSPA(TITLE)
!          PRINT*, TITLE

C-------...or resolution
        ELSE IF (LINE(1:21).EQ.'REMARK   2 RESOLUTION') THEN
C          WRITE (6,*) LINE
C         Check for empty values
          IF (INDEX(LINE, 'NOT') .GT. 0) GO TO 200 
          READ(UNIT=LINE(23:30), FMT='(A)') RTEMP
C          WRITE (6,*) RTEMP
          RTEMP = CHOP (RTEMP,8)
C          WRITE (6,*) RTEMP
          RESO  = GTREAL (RTEMP)
C          WRITE (6,*) RESO
          GTRESO=.TRUE.
          GO TO 200

C-------...or wavelength
CREMARK 200  WAVELENGTH OR RANGE        (A) : 1.28  
        ELSE IF (LINE(1:31).EQ.'REMARK 200  WAVELENGTH OR RANGE') THEN
C         Throw a way lines with multiple values
          IF (INDEX(LINE(45:80),',').NE.0) GO TO 200 
          IF (INDEX(LINE(45:80),';').NE.0) GO TO 200 
          
C         Extract the wavelength          
          WAVEL = GETVAL(LINE(45:80))
          GO TO 200          
  
C-------... or extract TLS information and write to file (if required)
        ELSE IF ((LINE(1:23).EQ.'REMARK   3  TLS DETAILS').AND.
     +(USETLS.EQV..TRUE.)) THEN
C-------Go to TLS extraction mode
          CALL GETTLS(TLSOUT,GOTTLS)

C-------... or find residues involved in LINKs
        ELSE IF (LINE(1:4).EQ.'LINK') THEN
C         Add the LINK to the cache 
          NLINK = NLINK+1
          LINK_CACHE(NLINK) = LINE
         
C---------Avoid overflow
          IF ((NO_BLD_CNT.GE.(MAXRES-1)).OR.
     +      (H2O_KEEP_CNT.GE.(MAXRES-1))) GO TO 997

          READ(UNIT=LINE(13:57), FMT=945, ERR=996) TMPATM1, TMPRESI1, 
     +      TMPCHN1, TMPRESN1, TMPINS1, TMPATM2, TMPRESI2, TMPCHN2, 
     +      TMPRESN2, TMPINS2
C---------Check for side chains and main chains
          IF (INDEX(RESIDUES, TMPRESI1).NE.0) THEN
C           If NOT backbone atom...
            IF (INDEX(BACKBONE, TMPATM1).EQ.0) THEN
              NO_BLD_CNT = NO_BLD_CNT+1
              NO_BLD_CHN(NO_BLD_CNT) = TMPCHN1
              NO_BLD_RES(NO_BLD_CNT) = TMPRESN1
              NO_BLD_INS(NO_BLD_CNT) = TMPINS1
C            WRITE(6,*) TMPCHN1, TMPRESN1, TMPINS1
C           If backbone atom
            ELSE
              IF (TMPATM1.EQ.' N  ') THEN
                IF (TMPATM2.NE.' C  ') THEN
                  BBN_KEEP_CNT = BBN_KEEP_CNT+1
                  BBN_KEEP_CHN(BBN_KEEP_CNT) = TMPCHN1
                  BBN_KEEP_RES(BBN_KEEP_CNT) = TMPRESN1
                  BBN_KEEP_INS(BBN_KEEP_CNT) = TMPINS1
                END IF
              ELSE IF (TMPATM1.EQ.' C  ') THEN
C               Do nothing                 
              ELSE
                BBO_KEEP_CNT = BBO_KEEP_CNT+1
                BBO_KEEP_CHN(BBO_KEEP_CNT) = TMPCHN1
                BBO_KEEP_RES(BBO_KEEP_CNT) = TMPRESN1
                BBO_KEEP_INS(BBO_KEEP_CNT) = TMPINS1
              END IF
            END IF            
          END IF
          IF (INDEX(RESIDUES, TMPRESI2).NE.0) THEN
C           If NOT backbone atom...
            IF (INDEX(BACKBONE, TMPATM2).EQ.0) THEN
              NO_BLD_CNT = NO_BLD_CNT+1
              NO_BLD_CHN(NO_BLD_CNT) = TMPCHN2
              NO_BLD_RES(NO_BLD_CNT) = TMPRESN2
              NO_BLD_INS(NO_BLD_CNT) = TMPINS2
C            WRITE(6,*) TMPCHN2, TMPRESN2, TMPINS2
C           If backbone atom
            ELSE
              IF (TMPATM2.EQ.' N  ') THEN
                IF (TMPATM1.NE.' C  ') THEN
                  BBN_KEEP_CNT = BBN_KEEP_CNT+1
                  BBN_KEEP_CHN(BBN_KEEP_CNT) = TMPCHN2
                  BBN_KEEP_RES(BBN_KEEP_CNT) = TMPRESN2
                  BBN_KEEP_INS(BBN_KEEP_CNT) = TMPINS2
                END IF
              ELSE IF (TMPATM2.EQ.' C  ') THEN
C               Do nothing  
              ELSE
                BBO_KEEP_CNT = BBO_KEEP_CNT+1
                BBO_KEEP_CHN(BBO_KEEP_CNT) = TMPCHN2
                BBO_KEEP_RES(BBO_KEEP_CNT) = TMPRESN2
                BBO_KEEP_INS(BBO_KEEP_CNT) = TMPINS2
              END IF
            END IF            
          END IF
C---------Check for waters
          IF (TMPRESI1.EQ.'HOH') THEN
             H2O_KEEP_CNT = H2O_KEEP_CNT+1
             H2O_KEEP_CHN(H2O_KEEP_CNT) = TMPCHN1
             H2O_KEEP_RES(H2O_KEEP_CNT) = TMPRESN1
             H2O_KEEP_INS(H2O_KEEP_CNT) = TMPINS1
          END IF
          IF (TMPRESI2.EQ.'HOH') THEN
             H2O_KEEP_CNT = H2O_KEEP_CNT+1
             H2O_KEEP_CHN(H2O_KEEP_CNT) = TMPCHN2
             H2O_KEEP_RES(H2O_KEEP_CNT) = TMPRESN2
             H2O_KEEP_INS(H2O_KEEP_CNT) = TMPINS2
          END IF
C          WRITE(6,*) TMPATM1, TMPRESI1, TMPCHN1,
C     +TMPRESN1, TMPINS1, TMPATM2, TMPRESI2, TMPCHN2, TMPRESN2, TMPINS2

C-------CRYST1 marks the end of the header          
        ELSE IF (LINE(1:6).EQ.'CRYST1') THEN
          GO TO 300     
        
C-------... or find other information
        ELSE
          CALL RDREMARK(LINE)
        END IF
200   CONTINUE

C-----Be relaxed about resolution?
300   IF ((GTRESO.EQV..FALSE.).AND.(RELAX.EQV..TRUE.).AND.
     +    (RESO3.GT.0.0)) THEN
        RESO   = RESO3
        GTRESO = .TRUE.
      END IF 

C-----Use which R-factor and R-free?
C      WRITE(6,*) RFAC_OBS, RFAC_CUT, RFAC_ALL
      IF ((RFAC_OBS.LT.0.99).AND.(RFAC_OBS.GT.0.00)) THEN
        RFACT = RFAC_OBS
      ELSE IF ((RFAC_CUT.LT.0.99).AND.(RFAC_CUT.GT.0.00)) THEN
        RFACT = RFAC_CUT
      ELSE
        RFACT = RFAC_ALL
      END IF

      IF ((RFREE_OBS.LT.0.99).AND.(RFREE_OBS.GT.0.00)) THEN
        RFREE = RFREE_OBS
      ELSE 
        RFREE = RFREE_CUT
      END IF
      DIFR = RFREE - RFACT

C-----Which refinement program was used?
      IF (RPROG(1:4).EQ.'NULL') THEN
        RPROG = 'UNKNOWN'
      END IF

      IF (RPROG.NE.'?') THEN
        PROG = NOSPCE(RPROG)
      ELSE
        PROG = 'UNKNOWN'
      END IF
      
C      PRINT*, LINK_CACHE
C====================End of PDB header section=========================C

C===========================Coordinate section=========================C

C-----Initialise
      ATMCNT   = 0
      HATMCT   = 0
      CHAINS   = 0
      BTOT     = 0.00
      BAVER    = 10.00
      ANISOU   = 0

C-----Read first atom record
      DO 351, I=1, MAXDAT
        READ(UNIT=7, FMT=920, END=400) LINE
        IF (LINE(1:6).EQ.'ANISOU') THEN
          ANISOU = ANISOU+1
        ELSE IF((LINE(1:6).EQ.'HETATM').OR.(LINE(1:6).EQ.'ATOM  ')) THEN
          ATMCNT = ATMCNT+1
          READ(UNIT=LINE(61:65), FMT=940) BTEMP
          BTOT = BTOT + BTEMP
          CHAINS = CHAINS+1
          READ(UNIT=LINE(22:22), FMT=970, ERR=351) CHAINI(CHAINS)
	  IF (CHAINI(CHAINS).EQ.' ') THEN
	    CHAINI(CHAINS) = 'X'
	  END IF
	  READ(UNIT=LINE(23:26), FMT=965, ERR=351) CFIRST(CHAINS)
	  READ(UNIT=LINE(23:26), FMT=965, ERR=351) CLAST(CHAINS)
	  GO TO 355
	END IF
351   CONTINUE

C-----Read other atoms
355   DO 360, I=1, MAXDAT
        READ(UNIT=7, FMT=920, END=380) LINE
	IF (LINE(1:6).EQ.'ENDMDL') THEN
	  WRITE(6,*) 'Multi-model refinement detected. ',
     +'Only evaluating first model.'
          GO TO 380
        ELSE IF (LINE(1:6).EQ.'ANISOU') THEN
          ANISOU = ANISOU+1
C         Iterate the loop immediately
          GO TO 360
        ELSE IF ((LINE(1:6).EQ.'HETATM').OR.(LINE(1:6).EQ.'ATOM  '))THEN
          MACRES = .FALSE.
          ATMCNT = ATMCNT+1
          IF (LINE(1:6).EQ.'HETATM') THEN
            HATMCT = HATMCT+1
          END IF
          READ(UNIT=LINE(61:65), FMT=940) BTEMP
          BTOT = BTOT + BTEMP
C         Get the residue ID and the residue number and insertion code
C         and the alternate specifier
          READ(UNIT=LINE, FMT=989, ERR=996) TMPALT, TMPRESI1, TMPCHN1,
     +      TMPRESN1, TMPINS1
C         Make sure the residue name is left aligned
          READ(UNIT=LINE(18:20),FMT=*, END=996) TMPRESI1
C         Set empty chains to 'X'
          IF (TMPCHN1.EQ.' ') THEN
            TMPCHN1 = 'X'
          END IF      
C         Is it a usable residue
          IF (((TMPRESN1.NE.OLDRESN).OR.(TMPCHN1.NE.OLDCHN))
     +        .AND.(TMPALT.EQ.' ')) THEN
C           Mark compounds of the non-polymer type as ligands
            IF (TMPRESI1.EQ.'WAT') GO TO 363
            DO 361, J=1, NHET
              IF (TMPRESI1.EQ.HETERO(J)) THEN 
C               Add it to the list
                OLDRESN = TMPRESN1
                OLDCHN  = TMPCHN1
                LIG_CNT = LIG_CNT+1
                LIG_CHN(LIG_CNT) = TMPCHN1
                LIG_RES(LIG_CNT) = TMPRESN1
                LIG_INS(LIG_CNT) = TMPINS1
C               It's a ligand so don't add it to the range for TLS 
C               assignment 
                GO TO 360
              END IF
361         CONTINUE
C           Also mark compounds not yet in the CCP4 dictionary
            DO 362, J=1, NALL
              IF (TMPRESI1.EQ.ALLRES(J)) THEN 
C               It's a known compound, so ignore it
                GO TO 363
              END IF 
362         CONTINUE
C           It is an unknown compound, treat it as a ligand
            OLDRESN = TMPRESN1
            OLDCHN  = TMPCHN1
            LIG_CNT = LIG_CNT+1
            LIG_CHN(LIG_CNT) = TMPCHN1
            LIG_RES(LIG_CNT) = TMPRESN1
            LIG_INS(LIG_CNT) = TMPINS1  
            GO TO 360       
          END IF

C---------Add residues to TLS ranges if it is a macromolecule residue
363       IF ((TMPRESI1.EQ.'HOH') .OR. (TMPRESI1.EQ.'WAT')) GO TO 360
C         Use a quick check for residues (works most of the time)
	  IF (INDEX(RESIDUES, LINE(18:20)).NE.0) THEN
	    PROT   =.TRUE.
            MACRES =.TRUE.
	  ELSE 
C           Check for nucleic acids or modified residues. Note that the 
C           residue names in NUCLEIC and AMINOA are left aligned whereas
C           names in the PDB file are right aligned. 
C           Nucleic acids
	    DO 364, J=1, NNUC
              READ(UNIT=LINE(18:20),FMT=*) JRES   
	      IF (JRES.EQ.NUCLEIC(J)) THEN
                MACRES=.TRUE.
                NUCLIC=.TRUE.
	        GO TO 366
	      END IF 
364         CONTINUE
C           Amino acids
            DO 365, J=1, NAMINO
              READ(UNIT=LINE(18:20),FMT=*) JRES   
	      IF (JRES.EQ.AMINOA(J)) THEN
                MACRES=.TRUE.
	        GO TO 366
	      END IF 
365         CONTINUE
          END IF

C---------Force a non-empty chain ID
366       IF (TMPCHN1.EQ.' ') THEN
	    TMPCHN1 = 'X'
	  END IF

C---------Is it a new chain (ignore waters)?
	  IF (MACRES.EQV..TRUE.) THEN
C           Loop over all known chains (IN INVERSE ORDER!!!!!!!!!)
            DO 367, J=CHAINS, 1, -1
	      IF (TMPCHN1.EQ.CHAINI(J)) THEN
	        MACMOL(J)=.TRUE.
C-------------  New range?
	        IF (TMPRESN1.LT.CFIRST(J)) THEN
	          CFIRST(J) = TMPRESN1
	        ELSE IF (TMPRESN1.GT.CLAST(J)) THEN
	          CLAST(J) = TMPRESN1
	        END IF
C               No thing else to do         
	        GO TO 360
              END IF
367         CONTINUE            
C           This is a new chain              
            CHAINS=CHAINS+1
            IF (CHAINS.EQ.MAXCHN) GOTO 995
	    MACMOL(CHAINS)=.TRUE.
	    CHAINI(CHAINS)=TMPCHN1
	    READ(UNIT=LINE(23:26), FMT=965, ERR=360) CFIRST(CHAINS)
	    READ(UNIT=LINE(23:26), FMT=965, ERR=360) CLAST(CHAINS)       
	  END IF
	END IF
360   CONTINUE 

C-----Calculate average B-factor
380   BAVER = BTOT / ATMCNT

C-----Write TLS file
      IF (FORTLS.EQ.0) THEN
        GO TO 400
      ELSE
C        PRINT*, NLINK, LINK_CACHE
        CALL MAKTLS(CHAINS, CHAINI, CFIRST, CLAST, MACMOL, TLSOUT2,
     +       LINK_CACHE, NLINK, TLSRES, NTLSRES)
      END IF

C======================End of coordinate section=======================C

C============================CIF section===============================C
C-----Skip this bit in EM mode
400   IF (EMMOD.EQV..TRUE.) GO TO 600

C-----Get input CIF file
      STATUS=0
      WRITE(6,*) ' '                                        
      CALL GETARG((2+EXTRA), INCIF)
      OPEN(UNIT=8, FILE=INCIF, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) INCIF, 'CIF file cannot be opened!'
        GO TO 999
      END IF 
      REWIND(8)   
      
C-----Initialise
      LABCNT=0
      GTTEST=.FALSE.   
      
C-----Skip CIF header
C-----This is done because the ciffile comes straight from cif2cif!!!!
      DO 410, I=1,MAXCIF
        READ(UNIT=8, FMT=920, END=998) LINE
        IF (LINE(1:1).EQ.'#') GO TO 410 
        IF (LINE(1:30).EQ.'_symmetry.space_group_name_H-M') THEN
          READ(LINE(31:),*) SYMMHM
C---------Special trick for space group P 1{bar}
          IF (SYMMHM.EQ.'P 1-        ') THEN
            SYMMHM = 'P -1        '
          ELSE IF (SYMMHM.EQ.'P 21 21 2 A  ') THEN
            SYMMHM = 'P 21 21 2 (a)  '
          END IF
        END IF   
C       Read the cell dimensions
        IF (LINE(1:15).EQ.'_cell.length_a ') THEN
          READ(LINE(16:),*) AAXIS
        END IF     
        IF (LINE(1:15).EQ.'_cell.length_b ') THEN
          READ(LINE(16:),*) BAXIS
        END IF 
        IF (LINE(1:15).EQ.'_cell.length_c ') THEN
          READ(LINE(16:),*) CAXIS
        END IF  
        IF (LINE(1:18).EQ.'_cell.angle_alpha ') THEN
          READ(LINE(19:),*) ALPHA
        END IF 
        IF (LINE(1:17).EQ.'_cell.angle_beta ') THEN
          READ(LINE(18:),*) BETA
        END IF      
        IF (LINE(1:18).EQ.'_cell.angle_gamma ') THEN
          READ(LINE(19:),*) GAMMA
        END IF 
        
        IF (LINE(1:5).EQ.'loop_') GO TO 410
        LINE = CHOP(LINE, 2)
        IF((LINE(1:7).EQ.'_refln.').OR.(LINE(1:7).EQ.'_refln_'))THEN
          LABCNT=LABCNT+1
          IF(LINE(1:13).EQ.'_refln.status')THEN
	    GTTEST=.TRUE.
	  END IF
        ELSE
C-------Reached the end of the header
          IF (LABCNT .GT. 1) THEN
	    BACKSPACE(8)
	    GO TO 500
	  END IF  
        END IF
410   CONTINUE

C-----Work on reflection data
C-----Initialise
500   HIRES=99.0
      LOWRES=0.0     
      REFCNT=0
      TSTCNT=0   
      
C-----Prepare for resolution calculation
      CALL RESCNST(AAXIS, BAXIS, CAXIS, ALPHA, BETA, GAMMA, COEFS)

C-----Load all reflections
      DO 510, I=1, MAXCIF
        READ(UNIT=8, FMT=920, END=520) LINE
C-------Ignore junk lines
        IF (LINE(1:20).EQ.'                    ') GO TO 510
	HASH = INDEX(LINE, '#')
        IF (HASH.NE.0) GO TO 510
C-------READ HKL FROM LINE
        IF (GTTEST.EQV..TRUE.) THEN
          REFCNT = REFCNT+1
          READ(UNIT=LINE, FMT=*, ERR=510) JUNK1, JUNK2, JUNK3, 
     +     INDEXH(REFCNT),
     +    INDEXK(REFCNT), INDEXL(REFCNT), SF(REFCNT), SIGSF(REFCNT), 
     +    TSTVAL(REFCNT)
          IF (TSTVAL(REFCNT).EQ.'f') THEN
	    TSTCNT = TSTCNT+1
	  END IF
	ELSE
	  REFCNT = REFCNT+1
          READ(UNIT=LINE, FMT=*, ERR=510) JUNK1, JUNK2, JUNK3, 
     +    INDEXH(REFCNT),
     +    INDEXK(REFCNT), INDEXL(REFCNT), SF(REFCNT), SIGSF(REFCNT)
	  TSTVAL(REFCNT)='o'
	END IF  
C-------Calculate resolution and SF/SIGSF
        REFRES(REFCNT) = HKLRES(INDEXH(REFCNT), INDEXK(REFCNT),
     +  INDEXL(REFCNT), COEFS)   
	RATIO(REFCNT)  = SF(REFCNT)/SIGSF(REFCNT)
	
C-------Verbose mode
        IF (VERBOS.EQV..TRUE.) THEN  
          WRITE(*,*) INDEXH(REFCNT), INDEXK(REFCNT), INDEXL(REFCNT),
     +    REFRES(REFCNT),  RATIO(REFCNT)
        END IF


C-------Skip values if the reflection is flagged as 'x' or F/SIGF < 0.00
        IF ((TSTVAL(REFCNT).NE.'x').AND.(RATIO(REFCNT).GE.1.0)) THEN
          IF (REFRES(REFCNT).LE.HIRES) THEN
            HIRES = REFRES(REFCNT)
	  END IF 
	  IF (REFRES(REFCNT).GE.LOWRES) THEN
            LOWRES = REFRES(REFCNT)
	  END IF
	END IF 
510   CONTINUE
 
C-----Calculate test set fraction.
520   TSTPRC = 0.0
      TSTPRC = REAL(TSTCNT)/REAL(REFCNT) * 100

C-----Get resolution from the data is in relaxed mode
      IF (GTRESO.EQV..FALSE.) THEN
        RESO   = HIRES
      END IF 

C========================End of CIF section============================C

C============================Output section============================C
C-----Get output file
600   STATUS=0                                     
      CALL GETARG((5+EXTRA), OUTDAT)
      OPEN(UNIT=9, FILE=OUTDAT, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) OUTDAT, 'cannot be created!'
        GO TO 999
      END IF
      

C-----Write to file
      WRITE(UNIT=9, FMT=930) AAXIS
      WRITE(UNIT=9, FMT=930) BAXIS
      WRITE(UNIT=9, FMT=930) CAXIS
      WRITE(UNIT=9, FMT=940) ALPHA
      WRITE(UNIT=9, FMT=940) BETA
      WRITE(UNIT=9, FMT=940) GAMMA
      WRITE(UNIT=9, FMT=940) RESO
      IF (EMMOD.EQV..FALSE.) THEN
        WRITE(UNIT=9, FMT=940) HIRES
        WRITE(UNIT=9, FMT=940) LOWRES  
      END IF   
      WRITE(UNIT=9, FMT=955) RFACT
      WRITE(UNIT=9, FMT=955) RFREE
C      WRITE(UNIT=9, FMT=940) BWILS
C-----Write NO_BLD list
      CALL RES_LIST (NO_BLD_CNT, NO_BLD_CHN, NO_BLD_RES, NO_BLD_INS)
      WRITE(UNIT=9, FMT=940) BAVER
      WRITE(UNIT=9, FMT=950) 'NULL'
      IF (EMMOD.EQV..FALSE.) THEN
        WRITE(UNIT=9, FMT=975) REFCNT
        WRITE(UNIT=9, FMT=975) TSTCNT
        WRITE(UNIT=9, FMT=980) TSTPRC
      END IF  
      WRITE(UNIT=9, FMT=925) PROG
      WRITE(UNIT=9, FMT=965) YEAR
      WRITE(UNIT=9, FMT=961) SYMMHM
C-----Write H2O_KEEP list
      CALL RES_LIST (H2O_KEEP_CNT, H2O_KEEP_CHN, H2O_KEEP_RES, 
     + H2O_KEEP_INS)
C-----Write BBN_KEEP list
      CALL RES_LIST (BBN_KEEP_CNT, BBN_KEEP_CHN, BBN_KEEP_RES, 
     + BBN_KEEP_INS)
C-----Write BBO_KEEP list
      CALL RES_LIST (BBO_KEEP_CNT, BBO_KEEP_CHN, BBO_KEEP_RES, 
     + BBO_KEEP_INS)
       WRITE(UNIT=9, FMT=985) PROT    
       WRITE(UNIT=9, FMT=905) MASPRO
       WRITE(UNIT=9, FMT=905) MASION
       WRITE(UNIT=9, FMT=905) MASSHR 
       IF (EMMOD.EQV..FALSE.) THEN
         WRITE(UNIT=9, FMT=990) COMPLET  
       END IF  
       CALL RES_LIST (LIG_CNT, LIG_CHN, LIG_RES, LIG_INS)
       WRITE(UNIT=9, FMT=985) NUCLIC
       WRITE(UNIT=9, FMT=956) WAVEL
       WRITE(UNIT=9, FMT=925) TITLE(1:LENSTR(TITLE))
C      WRITE(UNIT=9, FMT=975) ANISOU  
C      WRITE(UNIT=9, FMT=955) DIFR  
      
C-----Skip error messages
      GO TO 999
      
C===========================Formats====================================C
905   FORMAT(F4.2)
910   FORMAT(3F9.3,3F7.2)
915   FORMAT(I2)
920   FORMAT(A80)
925   FORMAT(A)
930   FORMAT(F9.3)
940   FORMAT(F6.2)
945   FORMAT(A4,1X,A3,1X,A1,I4,A1,15X,A4,1X,A3,1X,A1,I4,A1)
950   FORMAT(A4)
955   FORMAT(F6.4)
956   FORMAT(F7.5)
960   FORMAT(A11)
961   FORMAT(A13)
965   FORMAT(I4)
970   FORMAT(A1)
975   FORMAT(I8)
980   FORMAT(F4.1)
985   FORMAT(L1)
989   FORMAT(16X,A1,A3,1X,A1,I4,A1) 
990   FORMAT(F5.1) 
C===========================Error messages=============================C
995   WRITE(6,*) 'Too many chains. Increase MAXCHN and recompile.'
      GO TO 999
996   WRITE(6,*) 'Cannot interpret the PDB file:'
      WRITE(6,*) LINE
      GO TO 999
997   WRITE(6,*) 'Maximum number of residues exceeded!'
      WRITE(6,*) LINE
      GO TO 999
998   WRITE(6,*) 'Unexpected end of mmCIF file!'

C===========================End of the line============================C
999   CLOSE (7)
      CLOSE (8)
      CLOSE (9)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12) 
      END
  
C======================================================================C
C   SUBROUTINES AND FUNCTIONS                                          C
C======================================================================C

C-----------------------------------------------------------------------
C  This subroutine extracts TLS parameters from the PDB REMARK records
C  and writes them out to TLSOUT in REFMAC5 format
C
C12345678901234567890123456789012345678901234567890123456789012345678
C Refmac:
CREMARK   3  TLS DETAILS
CREMARK   3   NUMBER OF TLS GROUPS  : 16
CREMARK   3   NUMBER OF TLS GROUPS  : 1  
CREMARK   3
CREMARK   3   TLS GROUP : 1
CREMARK   3    NUMBER OF COMPONENTS GROUP : 1
CREMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI
CREMARK   3    RESIDUE RANGE :   A    -6        A   447
CREMARK   3    RESIDUE RANGE :   B       -8        B      220
CREMARK   3    ORIGIN FOR THE GROUP (A):  25.7178  19.0390  50.5421
CREMARK   3    T TENSOR
CREMARK   3      T11:   0.1788 T22:   0.4514
CREMARK   3      T33:   0.1143 T12:   0.0152
CREMARK   3      T13:   0.0228 T23:   0.0598
CREMARK   3    L TENSOR
CREMARK   3      L11:   0.4276 L22:   3.4374
CREMARK   3      L33:   4.1965 L12:   1.1072
CREMARK   3      L13:  -1.3016 L23:  -3.0033
CREMARK   3    S TENSOR
CREMARK   3      S11:   0.1322 S12:  -0.0064 S13:  -0.1383
CREMARK   3      S21:   0.0416 S22:   0.0605 S23:   0.0054
CREMARK   3      S31:   0.0448 S32:  -0.2763 S33:  -0.1928
CREMARK   3
C
C PHENIX
CREMARK   3   NUMBER OF TLS GROUPS  : 1                                 
CREMARK   3   TLS GROUP : 1                                             
CREMARK   3    SELECTION: CHAIN D AND RESSEQ 600-676                    
CREMARK   3    ORIGIN FOR THE GROUP (A):  19.1254   2.7088 121.6270
CREMARK   3    T TENSOR                                                 
CREMARK   3      T11:   1.5731 T22:   0.9883                            
CREMARK   3      T33:   0.5809 T12:  -0.0630                            
CREMARK   3      T13:   0.0873 T23:  -0.1817                 
CREMARK   3    L TENSOR                                                 
CREMARK   3      L11:  -2.2257 L22:  -0.6547                            
CREMARK   3      L33:  -2.9061 L12:   0.8561                            
CREMARK   3      L13:  -1.3039 L23:   0.9459                            
CREMARK   3    S TENSOR                                                 
CREMARK   3      S11:   0.3667 S12:   0.3662 S13:  -0.1024              
CREMARK   3      S21:  -0.3536 S22:  -0.3367 S23:   0.1040              
CREMARK   3      S31:  -0.0890 S32:  -0.0513 S33:   0.0680    
C BUSTER
CREMARK   3    ORIGIN FOR THE GROUP (A):   77.0443   19.4729   -5.2350
CREMARK   3    T TENSOR
CREMARK   3     T11:    0.0000 T22:    0.0000
CREMARK   3     T33:    0.0000 T12:    0.0000
CREMARK   3     T13:    0.0000 T23:    0.0000
CREMARK   3    L TENSOR
CREMARK   3     L11:    0.1106 L22:    1.0233
CREMARK   3     L33:    0.6541 L12:   -0.4855
CREMARK   3     L13:    0.0884 L23:   -0.7511
CREMARK   3    S TENSOR
CREMARK   3     S11:   -0.0061 S12:   -0.0010 S13:    0.0041
CREMARK   3     S21:    0.0098 S22:    0.0044 S23:   -0.0029
CREMARK   3     S31:   -0.0016 S32:    0.0018 S33:    0.0017
C123456789012345678901234567890123456789012345678901234567890
CTLS    Chain A
CRANGE  'A   1.' 'A 300.' ALL
CORIGIN   -3.680   0.700  11.039
CT    -0.0161 -0.1041 -0.0362  0.0121  0.0497  0.0194
CL     1.9572  2.9920  2.5163  0.3999 -0.3219  0.1108
CS     0.0066 -0.0199  0.0299 -0.0444 -0.1246 -0.0684  0.0286  0.0062
C Buster:
C 3MUH
CREMARK   3    SELECTION:   L  502    L  503                
C 2WZK
CREMARK   3    SELECTION: (CHAIN A AND RESID 12:18)
C 2XF6
CREMARK   3    SELECTION: (CHAIN A:2-18)
C 2XF7
CREMARK   3    SELECTION: (A3-18)
C 2XIK
CREMARK   3    SELECTION: (A3 - A22)    
C 2XJY
CREMARK   3    SELECTION: (CHAIN  A)  
CREMARK   3    SELECTION : (CHAIN N) 
C 2XC8
CREMARK   3    SELECTION: (CHAIN A 1-19) 
C 2X7F
CREMARK   3    SELECTION: ( A:15 - A:109 )                      
C 2X53
CREMARK   3    SELECTION: (A 2 - A 141)                        
CREMARK   3    SELECTION: (1 316 - 1 372) 
C 2X9B
CREMARK   3    SELECTION: ( A1 - A62 )
C 2XRG
CREMARK   3    SELECTION: (A56-A395,A402-A457, A470-A556, A563-A569)
C 3ZV0 
CREMARK   3    SELECTION: CHAIN D AND (RESID 19:42, 48:60, 258:377)
C 2ygc
CREMARK   3    SELECTION : (CHAIN C AND RESIDUES 15 - 51)
C LOCAL 
CREMARK   3    SET : { A|555 - A|792 } 
CREMARK   3    SET : { A|* }
CREMARK   3    SET : { A|885 - A|1009 }   
CREMARK   3    SET : { A|1010 - A|1017 }     
CREMARK   3    SET : { B|469 - B|549 B|553 - B|592 B|598 - B|683 }  
C-----------------------------------------------------------------------
      SUBROUTINE GETTLS(TLSOUT,GOTTLS)    
     
      IMPLICIT NONE
C-----Declare variables and constants
      INTEGER   MAXTLS, STATUS, LENSTR    
C-----MAXTLS is the maximum allowed TLS lines in a PDB header
      PARAMETER (MAXTLS=299)      
      
      CHARACTER TLSOUT*255, LINE*80, LINE2*80, CHAIN1*1, CHAIN2*1,
     +          LONG*650,   TLONG*210, ABC*52, CHOP*650, ID*1, TLINE*80
      PARAMETER (ABC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'//
     +               'abcdefghijklmnopqrstuvwxyz')
      INTEGER   GOTTLS, COLON, COMMA,  POS, DASH, CURLY, HSEL
      INTEGER   I, J, K, L, M, TLSCNT, TLSGRP, RESID1, RESID2
      REAL      T11, T12, T13, T22, T23, T33, ORI(3) 
      REAL      L11, L12, L13, L22, L23, L33
      REAL      S11, S12, S13, S21, S22, S23, S31, S32, S33
      REAL      S22MS11, S11MS33
      LOGICAL   USETLS

C-----Initialise values
      TLSCNT = 0 
      TLSGRP = 0
      ID     = '!'
      USETLS = .TRUE.
      
C-----Read the PDB file (Unit 7) line by line    
      DO 10, I=1,MAXTLS
        READ(UNIT=7, FMT=90, END=99) LINE
!        PRINT*, LINE
C-------Find the number of TLS groups    
        IF ((LINE(1:33).EQ.'REMARK   3   NUMBER OF TLS GROUPS') .OR.
     +      (LINE(1:32).EQ.'REMARK   3  NUMBER OF TLS GROUPS')) THEN
C---------Are there really TLS groups?
          IF (INDEX(LINE,'NULL').NE.0) GO TO 99
C---------Read number of groups
          COLON = INDEX(LINE,':')
	  READ(UNIT=LINE(COLON+1:55), FMT=*, ERR=98) TLSGRP
	  IF (TLSGRP.LE.0) THEN 
	    GO TO 99
	  ELSE 
	    GO TO 15
	  END IF
	END IF 
10    CONTINUE

15    GOTTLS=1
C-----Open output file
      OPEN(UNIT=10, FILE=TLSOUT, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) TLSOUT, 'cannot be created!'
        GO TO 99
      END IF
      
C-----Write TLS parameters to file
16    DO 20, J=1,MAXTLS
        IF (TLSCNT.EQ.TLSGRP) GO TO 99
        READ(UNIT=7, FMT=90, END=99) LINE
	IF (LINE(1:24).EQ.'REMARK   3   TLS GROUP :') THEN
	  TLSCNT = TLSCNT+1
	  WRITE(UNIT=10, FMT=92)'TLS    Group', TLSCNT
	  DO 25, K=1, 15
	    READ(UNIT=7, FMT=90, END=99) LINE
C-----------Read and write residue range   
C Format example: 
C123456789012345678901234567890123456789012345678901234567890     
CREMARK   3    RESIDUE RANGE :   A    -6        A   447     
CREMARK   3    RESIDUE RANGE :        -1             -1
	    IF (LINE(1:29).EQ.'REMARK   3    RESIDUE RANGE :') THEN
C             Catch cases of 'all residues in TLS'
              IF (LINE(30:54).EQ.'        -1             -1') THEN
                GOTTLS = 0
                CLOSE(10,STATUS='DELETE')
                GO TO 100
              END IF    
              READ(UNIT=LINE(33:),FMT=*,ERR=97,END=27) CHAIN1, RESID1,
     +           CHAIN2, RESID2     
              WRITE(UNIT=10, FMT=94)CHAIN1, RESID1, CHAIN2, RESID2
              GO TO 25
C             Fall back label in case of reading errors
27            READ(UNIT=LINE(33:), FMT=*, ERR=97) CHAIN1, RESID1,
     +          RESID2  
              CHAIN2 = CHAIN1
              WRITE(UNIT=10, FMT=94)CHAIN1, RESID1, CHAIN2, RESID2  
CREMARK   3    SET : { A|555 - A|792 } 
CREMARK   3    SET : { A|885 - A|1009 }   
CREMARK   3    SET : { A|1010 - A|1017 } 
CREMARK   3    SET : { B|469 - B|549 B|553 - B|592 B|598 - B|683 }
CREMARK   3    SELECTION: {A|-1 - 51 A|59 - 68}
CREMARK   3    SELECTION: { A|504 - A|504 C|503 - C|503 E|502 -   
CREMARK   3               E|503 D|503 - D|503 F|503 - F|503 }         
            ELSE IF ((INDEX(LINE,'{').NE.0).OR.
     +               (INDEX(LINE,'SET').NE.0))  THEN
              LONG = LINE
              DO 28, L=1, 3
                READ(UNIT=7, FMT=90, END=99) LINE2
                IF (INDEX(LINE2,'ORIGIN').EQ.0) THEN
                  LONG = LONG(1:LENSTR(LONG))//LINE2(11:80)
                ELSE
                  BACKSPACE(7)
                  GO TO 29
                END IF
28            CONTINUE              
29            COLON = INDEX(LINE,':')
C              PRINT*, LONG(CURLY:LENSTR(LONG)
              CALL CBUSTER(LONG(COLON+1:))
CREMARK   3    SELECTION: CHAIN D AND RESSEQ 600-676 
            ELSE IF (LINE(1:23).EQ.'REMARK   3    SELECTION') THEN
C-------------PHENIX or dirty BUSTER format (might be buggy)
              COLON = INDEX(LINE,':')
              LONG = LINE(COLON+1:)
C-------------Read ahead to find entire selection
              DO 30, L=1, 10
                READ(UNIT=7, FMT=90, END=99) LINE2
                IF (INDEX(LINE2,'ORIGIN').EQ.0) THEN
                  LONG = LONG(1:LENSTR(LONG))//LINE2(COLON+1:80)
                ELSE
                  BACKSPACE(7)
                  GO TO 31
                END IF
30            CONTINUE
C-------------Check for selections comprising all atoms. In that case 
C             ignore the TLS group.
C 2xqt:
CREMARK   3    SELECTION: ALL
CREMARK   3    SELECTION: NULL
31            IF ((INDEX(LONG,'ALL').NE.0).OR.
     +            (INDEX(LONG,'NULL').NE.0)) THEN
                GOTTLS = 0
                CLOSE(10,STATUS='DELETE')
                GO TO 100
              END IF

C-------------Check for negative selectors
C 4bgz:
CREMARK   3    SELECTION: chain G and not (resseq 138:230 or resseq 366:
C 3mna:
CREMARK   3    SELECTION: chain A and not element H
              IF ((INDEX(LONG,'not').NE.0).OR.
     +            (INDEX(LONG,'NOT').NE.0)) THEN
C               See if we can deal with this
                CALL TLSNOT(LONG,USETLS) 
C               Stop parsing TLS if the selction cannot be used
                IF (USETLS.EQV..FALSE.) THEN
                   WRITE(6,*) 'Cannot parse groups with negative',
     +             ' selectors'
                  GOTTLS = 0
                  CLOSE(10,STATUS='DELETE')
                  GO TO 100
                END IF             
              END IF
C 2xrg:
CREMARK   3    SELECTION: (A56-A395,A402-A457, A470-A556, A563-A569)
C 3zv0:
CREMARK   3    SELECTION: CHAIN D AND (RESID 19:42, 48:60, 258:377)
C 4f88:
CREMARK   3    SELECTION:   2    2    2  462
C
C-------------If there are commas on the line, split and digest the
C             sections separately
C             Cache the line
              COMMA = INDEX(LONG,',')
              IF (COMMA.NE.0) THEN
                TLONG = LONG(1:COMMA-1)
C               Use the first non-blank character after 'CHAIN'
                POS = INDEX(TLONG,'CHAIN')
                ID  = CHOP(TLONG(POS+5:LENSTR(TLONG)),2)(1:1) 
              END IF

              DO 32, I=0, 20
                COMMA = INDEX(LONG,',')
                IF (COMMA.NE.0) THEN
                  TLONG = LONG(1:COMMA-1)
C                 Test the substring for chain ids
                  DO 33, M=1, LENSTR(TLONG)
                    IF (INDEX(ABC,LONG(M:M)).NE.0) THEN
                      CALL PHENIX (TLONG//'      ')
                      GO TO 34
                    END IF
33                CONTINUE
C                 Did not find a chain ID so get it from the cache:
C                 Use the first non-blank character after 'CHAIN'
                  TLONG = 'CHAIN '//ID//' RESID '//LONG(1:COMMA-1) 
                  CALL PHENIX (TLONG)
34                LONG = LONG(COMMA+1:LENSTR(LONG))
                ELSE                  
C                 Test the substring for chain ids
                  DO 35, M=1, LENSTR(LONG)
                    IF (INDEX(ABC,LONG(M:M)).NE.0) THEN
                      CALL PHENIX (LONG(1:LENSTR(LONG))//'      ')
C                     End of the selection
                      GO TO 25
                    END IF
35                CONTINUE
C                 No chainID, get it from the cache
                  IF (ID.NE.'!') THEN
                    TLONG ='CHAIN '//ID//' RESID '//LONG(1:LENSTR(LONG))
                    CALL PHENIX(TLONG//'      ')
                  ELSE
C                   The cache is empty, treat this as dirty buster data
                    CALL BUSTER(LONG(1:LENSTR(LONG))//'      ')
                  END IF  
C                 End of the selection
                  GO TO 25
                END IF
32            CONTINUE
C---------Read and write origin and tensors (if available) use fallback
C         for BUSTER format
            ELSE IF (LINE(1:29).EQ.'REMARK   3    ORIGIN FOR THE ') THEN
              IF (LINE(44:44).EQ.'.') THEN
	        READ(UNIT=LINE(40:66),FMT=89,ERR=96)ORI(1),ORI(2),ORI(3)
              ELSE
                READ(UNIT=LINE(40:69),FMT=83,ERR=96)ORI(1),ORI(2),ORI(3)
              END IF  
C             Read T tensor  
42            READ(UNIT=7, FMT=90, END=99,ERR=96) LINE
	      READ(UNIT=7, FMT=88, END=99,ERR=43) T11, T22
              GO TO 44
43            BACKSPACE(7)
              READ(UNIT=7, FMT=82, END=99,ERR=96) T11, T22
44            READ(UNIT=7, FMT=88, END=99,ERR=45) T33, T12
              GO TO 46
45            BACKSPACE(7)
              READ(UNIT=7, FMT=82, END=99,ERR=96) T33, T12
46            READ(UNIT=7, FMT=88, END=99,ERR=47) T13, T23
              GO TO 48
47            BACKSPACE(7)
              READ(UNIT=7, FMT=82, END=99,ERR=96) T13, T23
C             Read L tensor  
48            READ(UNIT=7, FMT=90, END=99,ERR=96) LINE
	      READ(UNIT=7, FMT=88, END=99,ERR=49) L11, L22
              GO TO 50
49            BACKSPACE(7)
              READ(UNIT=7, FMT=82, END=99,ERR=96) L11, L22
50            READ(UNIT=7, FMT=88, END=99,ERR=51) L33, L12
              GO TO 52
51            BACKSPACE(7)
              READ(UNIT=7, FMT=82, END=99,ERR=96) L33, L12
52            READ(UNIT=7, FMT=88, END=99,ERR=53) L13, L23
              GO TO 54
53            BACKSPACE(7)
              READ(UNIT=7, FMT=82, END=99,ERR=96) L13, L23
C             Read S tensor  
54            READ(UNIT=7, FMT=90, END=99,ERR=96) LINE
	      READ(UNIT=7, FMT=87, END=99,ERR=55) S11, S12, S13
              GO TO 56
55            BACKSPACE(7)
              READ(UNIT=7, FMT=81, END=99,ERR=96) S11, S12, S13
56            READ(UNIT=7, FMT=87, END=99,ERR=57) S21, S22, S23
              GO TO 58
57            BACKSPACE(7)
              READ(UNIT=7, FMT=81, END=99,ERR=96) S21, S22, S23
58            READ(UNIT=7, FMT=87, END=99,ERR=59) S31, S32, S33
              GO TO 60
59            BACKSPACE(7)
              READ(UNIT=7, FMT=81, END=99,ERR=96) S31, S32, S33
60            S22MS11 = S22 - S11
	      S11MS33 = S11 - S33
C-----------  Write the lot    
	      WRITE(UNIT=10, FMT=86) ORI(1), ORI(2), ORI(3)
	      WRITE(UNIT=10, FMT=85) 'T', T11, T22, T33, T12, T13, T23
	      WRITE(UNIT=10, FMT=85) 'L', L11, L22, L33, L12, L13, L23
	      WRITE(UNIT=10, FMT=84) 'S', S22MS11, S11MS33, S12, S13,
     +                             S23, S21, S31, S32
C-----------  This is the end of the TLS group, so break the loop
              GO TO 26
            END IF
25        CONTINUE    
26        WRITE(10,*) ' ' 
        END IF
20    CONTINUE
      
C-----Skip error messages 
      GO TO 99
C-----Format statements
78    FORMAT(I4) 
79    FORMAT(A,I3)
80    FORMAT(A20)
81    FORMAT(19X,F10.4,5X,F10.4,5X,F10.4)
82    FORMAT(19X,F10.4,5X,F10.4)
83    FORMAT(3F10.4)   
84    FORMAT(A1,2X,8F10.4)
85    FORMAT(A1,2X,6F10.4)
86    FORMAT('ORIGIN ',3(F9.4,1X))
87    FORMAT(20X,F9.4,5X,F9.4,5X,F9.4)
88    FORMAT(20X,F9.4,5X,F9.4)
89    FORMAT(3F9.4)
90    FORMAT(A80)
91    FORMAT(I5)
92    FORMAT(A12,1X,I3)
93    FORMAT(A1)
94    FORMAT('RANGE  ''',A1,I4,'.'' ''',A1,I4,'.'' ALL')
      
C-----Error messages
96    WRITE(6,*) 'Cannot read TLS tensor'
      WRITE(10,*) ' '
      GO TO 16
97    WRITE(6,*) 'Cannot read TLS residue range'
      GO TO 99
98    WRITE(6,*) 'Cannot read number of TLS groups'           

99    WRITE(6,79) 'TLS groups extracted: ',TLSCNT
      CLOSE(10)
100   RETURN
      END      

C-----------------------------------------------------------------------
C  This subroutine reads a PHENIX style TLS group selection and writes 
C  it out in Refmac format. Currently very ugly.      
C
C Format examples: 
C123456789012345678901234567890123456789012345678901234567890         
CREMARK   3    SELECTION: CHAIN D AND RESSEQ 600-676 
CREMARK   3    SELECTION: (CHAIN A AND RESID 2:20)  
CREMARK   3    SELECTION: (CHAIN A OR CHAIN E) AND (RESID 201:435 OR RE 
CREMARK   3     
C 3kdv:
C3    SELECTION: ((CHAIN A AND RESID 179:181) OR (CHAIN B AND RESID     
C3               179:185) OR (CHAIN C AND RESID 179:189) OR (CHAIN D    
C3               AND RESID 179:185) OR (CHAIN E AND RESID 179:184)) 
C 2x4n:
C     SELECTION: (CHAIN A AND RESID :180)
C     SELECTION: (CHAIN B AND RESID 0:)
C 2xig:
CREMARK   3    SELECTION: CHAIN A OR CHAIN B  
C 3csk:
CREMARK   3    SELECTION: CHAIN A AND RESID 2-310 
C PoR:
CREMARK   3    SELECTION: (chain B and resid 152:182)
C 3prx:
CREMARK   3    SELECTION: (RESID 918:969 OR RESID 1270:1333) AND CHAIN B
CREMARK   3    SELECTION: (RESID 567:606 OR RESID 760:821) AND CHAIN A 
C 4h57:
CREMARK   3    SELECTION: chain 'A' and (resid 26 through 122 )
C 4bgz:
CREMARK   3    SELECTION: chain G and not (resseq 138:230 or resseq 366:
C-----------------------------------------------------------------------
      SUBROUTINE PHENIX(LINE)    
     
      IMPLICIT NONE
C-----Declare variables and constants
      INTEGER   MAXSEL, LOW, HIGH
C-----MAXSET is the maximum allowed selection groups
      PARAMETER (MAXSEL=40)  
      PARAMETER (LOW   =-10)       
      PARAMETER (HIGH  =9999) 
      
C-----ONETEN are all possible numbers
      CHARACTER ONETEN*10
      PARAMETER (ONETEN='1234567890')   
      
      CHARACTER LINE*(*), TYP*1, CHOP*650
C-----TYP is the type of value. 'R' for range, 'C' for chain. 
      CHARACTER CHID(MAXSEL)*1
      INTEGER   I, POS, POSCHN, POSRES, J, CHCNT, RANCNT, RESID1,LVL,
     +          RESID2, ORE, COLON, LENSTR, LENGTH, HAAKJE, SELCNT,
     +          DATCHK, SPACE, INHAAK
      INTEGER   RESIDF(MAXSEL), RESIDL(MAXSEL), CHLEV(MAXSEL), 
     +          RANLEV(MAXSEL), COMPL(MAXSEL)
      
C-----Initialise
      POS    = 0
      SPACE  = 0
      CHCNT  = 0
      RANCNT = 0
      SELCNT = 0
      LVL    = 1
      LENGTH = LENSTR(LINE)
      INHAAK = 0
      DO 1, I=1, MAXSEL
        CHLEV(I)  = 0
        RANLEV(I) = 0
        COMPL(I)  = 0
1     CONTINUE
      TYP = ' ' 
      
C-----Clean up the data line
      CALL UPCASE(LINE)

C-----Make the line 'compact' and reset the length
C     PRINT*, 'Upcase', LINE
      CALL DUBSPA(LINE)
      LINE   = CHOP(LINE,1)
      LENGTH = LENSTR(LINE)
      POS    = 0
C      PRINT*, 'DUBSPA ', LINE     

C-----Minusses (used as a range indicator)
      DO 2, I=1, MAXSEL
        POS = INDEX(LINE,'-')
        IF (POS.EQ.0) GO TO 3
C-------Make sure this isn't an actual minus sign
        CALL MINFIX(LINE)
C      PRINT*, 'Minfix', LINE
2     CONTINUE

C-----Replace 'TO' as a range indicator
3     POS = INDEX(LINE,' TO ')
      IF (POS.NE.0) THEN
        LINE = LINE(1:POS-1)//':'//LINE(POS+4:)
C        PRINT*, LINE
      END IF

C-----Remove space in lines of the type ': 1234'
      IF (INDEX(LINE,': ').NE.0) THEN
        CALL SPAFIX(LINE)
      END IF

C-----RESI (instead of RESID)
      DO 4, I=1, MAXSEL
        POS = INDEX(LINE,'RESI ')
        IF (POS.EQ.0) GO TO 5
        LINE = LINE(1:POS-1)//'RESID '//LINE(POS+5:)
C        PRINT*, LINE
4     CONTINUE

C-----RESIDUE (instead of RESID)
5     DO 6, I=1, MAXSEL
        POS = INDEX(LINE,'RESIDUES')
        IF (POS.EQ.0) GO TO 8
        LINE = LINE(1:POS-1)//'RESID'//LINE(POS+8:)
C        PRINT*, LINE
6     CONTINUE

C-----RESIDUE (instead of RESID)
8     DO 9, I=1, MAXSEL
        POS = INDEX(LINE,'RESIDUE')
        IF (POS.EQ.0) GO TO 10
        LINE = LINE(1:POS-1)//'RESID'//LINE(POS+7:)
C        PRINT*, LINE
9     CONTINUE

C-----RESSEQ
10    DO 11, I=1, MAXSEL
        POS = INDEX(LINE,'RESSEQ')
        IF (POS.EQ.0) GO TO 12
        LINE = LINE(1:POS-1)//'RESID'//LINE(POS+6:)
C        PRINT*, LINE
11    CONTINUE

C-----SEGID
12    DO 13, I=1, MAXSEL
        POS = INDEX(LINE,'SEGID')
        IF (POS.EQ.0) GO TO 16
        LINE = LINE(1:POS-1)//'CHAIN'//LINE(POS+5:)
C        PRINT*, LINE
13    CONTINUE

C-----THROUGH
16    DO 17, I=1, MAXSEL
        POS = INDEX(LINE,' THROUGH ')
        IF (POS.EQ.0) GO TO 21
        LINE = LINE(1:POS-1)//':'//LINE(POS+9:)
17    CONTINUE

C-----Clear values
      POS = 1
21    DO 22, I=1, MAXSEL
        LENGTH = LENSTR(LINE)
        HAAKJE = INDEX(LINE(POS:LENGTH),')')+POS-1
        IF (HAAKJE.EQ.POS-1) GO TO 23
        IF (LINE(HAAKJE-1: HAAKJE-1).EQ.' ') THEN
          POS = HAAKJE+1
        ELSE
          LINE = LINE(1:HAAKJE-1)//' )'//LINE(HAAKJE+1:LENGTH)
          POS = HAAKJE+2
        ENDIF  
22    CONTINUE


C-----Remove double quotes
23    DO 24, I=1, MAXSEL
        POS = INDEX(LINE,'"')
        IF (POS.EQ.0) GO TO 26
        LINE = LINE(1:POS-1)//LINE(POS+1:LENGTH)
24    CONTINUE

C-----Remove single quotes
26    DO 27, I=1, MAXSEL
        POS = INDEX(LINE,'''')
        IF (POS.EQ.0) GO TO 28
        LINE = LINE(1:POS-1)//LINE(POS+1:LENGTH)
27    CONTINUE

C-----Clean additional double spaces and reestablish the line length
28    CALL DUBSPA(LINE)
      LENGTH = LENSTR(LINE)
C      PRINT*, LINE
C-----Switch to BUSTER mode if the 'RESID' keyword is not given, but 
C     there are numbers on the line.
      IF (INDEX(LINE, 'RESID').EQ.0) THEN
        DO 29, I=1, LENGTH
          IF (INDEX(ONETEN,LINE(I:I)).NE.0) THEN
C            PRINT*, 'Going to BUSTER with ', LINE
            CALL BUSTER(LINE) 
            GO TO 999
          END IF
29      CONTINUE       
      END IF

C      PRINT*, LINE

C-----Digest the selection  
      LENGTH = LENSTR(LINE)
      DO 100, I=1, LENGTH
        POS = POS+1
C        PRINT*, LINE(POS:LENGTH)
        IF (POS.GT.LENGTH) GO TO 200
C       Look for keywords
        IF (INDEX('RC(OA)',LINE(POS:POS)).EQ.0) THEN
          GO TO 100
        ELSE IF (LINE(POS:POS).EQ.'(') THEN 
          INHAAK = INHAAK+1
        ELSE IF (LINE(POS:POS).EQ.')') THEN 
          INHAAK = INHAAK-1
        ELSE IF (LINE(POS:POS).EQ.'R') THEN
          IF (LINE(POS:POS+4).EQ.'RESID') THEN
            POS = POS+6
            RANCNT = RANCNT+1
C           Update chain level?
            IF (TYP.NE.'R') THEN
              IF (RANCNT.EQ.1) THEN
                CHLEV(RANCNT) = 1  
              ELSE
                CHLEV(RANCNT) = CHLEV(RANCNT-1)+1
              END IF
            ELSE
              CHLEV(RANCNT) = CHLEV(RANCNT-1)
            END IF
            TYP = 'R'
C           We can now have 1) '123'
C                           2) '123 stuff:more_stuff' 
C                           3) ':456'
C                           4) '123:'
C                           5) '123:456'
C           Find the first ':'
            COLON = INDEX(LINE(POS:LENGTH),':')+POS-1
            SPACE = INDEX(LINE(POS:LENGTH+1),' ')+POS-1
            IF (COLON.EQ.POS-1) THEN
C             Stop with empty line
              IF (DATCHK(LINE(POS:LENGTH)).EQ.0) GO TO 997
              READ (UNIT=LINE(POS:SPACE), FMT=*, ERR=997) RESIDF(RANCNT)
              RESIDL(RANCNT) = RESIDF(RANCNT)
              POS = SPACE
            ELSE IF (SPACE.LT.COLON) THEN
              READ (UNIT=LINE(POS:SPACE), FMT=*, ERR=997) RESIDF(RANCNT)
              RESIDL(RANCNT) = RESIDF(RANCNT)
              POS = SPACE
            ELSE IF (COLON.EQ.POS) THEN
              READ (UNIT=LINE(POS+1:SPACE),FMT=*,ERR=997) RESIDL(RANCNT)
              RESIDF(RANCNT) = 1
              POS = SPACE
            ELSE IF (COLON+1.EQ.SPACE) THEN
              READ (UNIT=LINE(POS:COLON-1),FMT=*,ERR=997) RESIDF(RANCNT)
              RESIDL(RANCNT) = 9999
              POS = SPACE
            ELSE
              READ(UNIT=LINE(POS:COLON-1),FMT=*,ERR=997)RESIDF(RANCNT)
              READ(UNIT=LINE(COLON+1:SPACE),FMT=*,ERR=997)RESIDL(RANCNT)
              POS = SPACE
            END IF
          ELSE
            GO TO 997
          END IF
        ELSE IF (LINE(POS:POS).EQ.'C') THEN
          IF (LINE(POS:POS+4).EQ.'CHAIN') THEN
C            PRINT*, LINE(POS:)
            POS   = POS+6
C-----------Stop if we reached the end of the line
            IF (POS.GT.LENGTH) GO TO 997
            CHCNT = CHCNT+1
C           Update range level?
            IF (TYP.NE.'C') THEN
              IF (CHCNT.EQ.1) THEN
                RANLEV(CHCNT) = 1
              ELSE
                RANLEV(CHCNT) = RANLEV(CHCNT-1)+1
              END IF
            ELSE
              RANLEV(CHCNT) = RANLEV(CHCNT-1)
            END IF
            TYP = 'C'
            READ(UNIT=LINE(POS:POS),FMT=*,ERR=997)CHID(CHCNT)
            POS = POS+1
C            PRINT*, CHID(CHCNT)
          ELSE IF (LINE(POS:POS+1).EQ.'C ') THEN
C 4bnc
CREMARK   3    SELECTION: CHAIN B AND C 
C  Assign chain ID 'C'
            CHCNT = CHCNT+1
C           Update range level?
            IF (TYP.NE.'C') THEN
              IF (CHCNT.EQ.1) THEN
                RANLEV(CHCNT) = 1
              ELSE
                RANLEV(CHCNT) = RANLEV(CHCNT-1)+1
              END IF
            ELSE
              RANLEV(CHCNT) = RANLEV(CHCNT-1)
            END IF
            TYP = 'C'
            CHID(CHCNT) = 'C'
            POS   = POS+1         
          ELSE
            GO TO 997
          END IF 
        ELSE IF (LINE(POS:POS).EQ.'O') THEN
          IF (LINE(POS:POS+1).EQ.'OR') THEN
            POS = POS+2
          ELSE
            GO TO 997
          END IF 
        ELSE IF (LINE(POS:POS).EQ.'A') THEN
          IF (LINE(POS:POS+2).EQ.'AND') THEN
            POS = POS+3
          ELSE
            GO TO 997
          END IF 
        END IF
100   CONTINUE


C-----Write out the chains and ranges
200   DO 201, I=1, CHCNT
        SELCNT = 0 
        DO 211, J=1, RANCNT
          IF (RANLEV(I).EQ.CHLEV(J))THEN
            WRITE(UNIT=10, FMT=994)CHID(I), RESIDF(J), 
     +      CHID(I), RESIDL(J)
            SELCNT = SELCNT +1
          END IF
211     CONTINUE 
        IF (SELCNT.EQ.0) THEN
          WRITE(UNIT=10, FMT=994)CHID(I), LOW, CHID(I), HIGH
        END IF    
201   CONTINUE
      
C-----Skip error messages 
      GO TO 999
C-----Format statements
993   FORMAT(A1)
994   FORMAT('RANGE  ''',A1,I4,'.'' ''',A1,I4,'.'' ALL')         
      
C-----Error messages
996   WRITE(6,*) 'Fuck'
997   WRITE(6,*) 'Cannot read TLS residue range:'
      WRITE(6,*) LINE

999   RETURN
      END      

C-----------------------------------------------------------------------
C  This subroutine reads a dirty BUSTER style TLS group selection and 
C  writes it out in Refmac format. Currently very ugly.      
C
C Format examples: 
C123456789012345678901234567890123456789012345678901234567890         
C 3MUH
CREMARK   3    SELECTION:   L  502    L  503                           
C 2WZK
CREMARK   3    SELECTION: (CHAIN A AND RESID 12:18)
C 2XF6
CREMARK   3    SELECTION: (CHAIN A:2-18)
C 2XF7
CREMARK   3    SELECTION: (A3-18)
C 2XIK
CREMARK   3    SELECTION: (A3 - A22)    
C 2XJY
CREMARK   3    SELECTION: (CHAIN  A)      
C 2XC8
CREMARK   3    SELECTION: (CHAIN A 1-19) 
C 2X7F
CREMARK   3    SELECTION: ( A:15 - A:109 )                      
C 2X53
CREMARK   3    SELECTION: (A 2 - A 141)                        
CREMARK   3    SELECTION: (1 316 - 1 372)  
C 2YAD
CREMARK   3    SELECTION : (A90 - A149 OR A181 - A197)
CREMARK   3    SET : (B89 - B149 OR B180 - B197)
C From CBUSTER:
C B AND RESEQ 1270 
C-----------------------------------------------------------------------
      SUBROUTINE BUSTER(LINE)    
     
      IMPLICIT NONE
C-----Declare variables and constants    
      CHARACTER LINE*(*), T1*5, T2*5, T3*5, T4*5, CHAIN1*1, CHAIN2*1,
     +          CLINE*200
      INTEGER   DATCHK, VALUES, LENSTR, RESID1, RESID2, POS, I,K

C-----Initialise
C      PRINT*, 'Dirty Buster ', LINE 
C-----Clean up the data line

C-----Make the line 'compact', and write to the cache
      CALL DUBSPA(LINE)
      CLINE =  LINE(:LENSTR(LINE))
      
C-----Split the lines by ORs
      DO 1, K=1, 10
        IF (LENSTR(CLINE).GT.0) THEN
          POS = INDEX(CLINE,'OR')
          IF (POS.NE.0) THEN
            LINE  = CLINE(1:POS-1)
            CLINE = CLINE(POS+2:)
          ELSE
            LINE  = CLINE
            CLINE = ''
          ENDIF  
C          PRINT*, CLINE
        ELSE  
          GO TO 99
        END IF

C-------Remove '|'
        DO 2, I=1, 8
          POS = INDEX(LINE,'|')
          IF (POS.EQ.0) GO TO 3
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
2       CONTINUE

C-------Remove 'RESEQ'
3       DO 4, I=1, 8
          POS = INDEX(LINE,'RESEQ')
          IF (POS.EQ.0) GO TO 5
          LINE = LINE(1:POS-1)//' '//LINE(POS+5:LENSTR(LINE))
4       CONTINUE

C-----  Remove 'AND'
5       DO 6, I=1, 8
          POS = INDEX(LINE,'AND')
          IF (POS.EQ.0) GO TO 11
          LINE = LINE(1:POS-1)//' '//LINE(POS+3:LENSTR(LINE))
6       CONTINUE

C-----  Remove 'CHAIN'
11      DO 15, I=1, 2
          POS = INDEX(LINE,'CHAIN')
          IF (POS.EQ.0) GO TO 16
          LINE = LINE(1:POS-1)//' '//LINE(POS+5:LENSTR(LINE))
15      CONTINUE

C-----  Remove '('
16      DO 17, I=1, 2
          POS = INDEX(LINE,'(')
          IF (POS.EQ.0) GO TO 18
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
17      CONTINUE

C-----  Remove ')'
18      DO 19, I=1, 2
          POS = INDEX(LINE,')')
          IF (POS.EQ.0) GO TO 21
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
19      CONTINUE

C-----  Remove ':'
21      DO 25, I=1, 2
          POS = INDEX(LINE,':')
          IF (POS.EQ.0) GO TO 26
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
25      CONTINUE

C-----  Remove '- '
26      DO 30, I=1, 4
          POS = INDEX(LINE,'- ')
          IF (POS.EQ.0) GO TO 31
          LINE = LINE(1:POS-1)//' '//LINE(POS+2:LENSTR(LINE))
30      CONTINUE

C-----  Remove '{'
31      DO 35, I=1, 2
          POS = INDEX(LINE,'{')
          IF (POS.EQ.0) GO TO 36
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
35      CONTINUE

C-----  Remove '}'
36      DO 40, I=1, 2
          POS = INDEX(LINE,'}')
          IF (POS.EQ.0) GO TO 41
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
40      CONTINUE


41      VALUES = DATCHK(LINE)
C        PRINT*, VALUES, LINE

C-----  Read values first, then check them
        IF (VALUES.EQ.1) THEN
          READ(UNIT=LINE, FMT=*, ERR=96) T1
          IF (LENSTR(T1).EQ.1) THEN
            CHAIN1 = T1(1:1)
C---------  These values are implied
            CHAIN2 = CHAIN1
            RESID1 = -10
            RESID2 = 9999
          ELSE
C---------  ChainID and residue number must be fused
            CHAIN1 = T1(1:1)
            READ(UNIT=T1(2:5), FMT='(I4)', ERR=96) RESID1
C---------  These values are implied
            CHAIN2 = CHAIN1
            RESID2 = RESID1
          END IF
        ELSE IF (VALUES.EQ.2) THEN
          READ(UNIT=LINE, FMT=*, ERR=96) T1, T2
C-------  Sanity check and clean read
          IF (LENSTR(T1).EQ.1) THEN
            CHAIN1 = T1(1:1)
            READ(UNIT=T2, FMT='(I5)', ERR=96) RESID1
C---------  These values are implied
            CHAIN2 = CHAIN1
            RESID2 = RESID1
          ELSE 
C---------  ChainID and residue number must be fused
            CHAIN1 = T1(1:1)
            READ(UNIT=T1(2:5), FMT='(I4)', ERR=96) RESID1
C---------  T2 should now be a residue# or a fused chainID and residue#
            READ(UNIT=T2, FMT='(I5)', ERR=50) RESID2
            CHAIN2 = CHAIN1
            GO TO 80
C---------  Got a reading error: fused values
50          CHAIN2 = T2(1:1)
            READ(UNIT=T2(2:5), FMT='(I4)', ERR=96) RESID2
          END IF
        ELSE IF (VALUES.EQ.3) THEN
          READ(UNIT=LINE, FMT=*, ERR=96) T1, T2, T3
C-------  Sanity check and clean read
          IF (LENSTR(T1).EQ.1) THEN
            CHAIN1 = T1(1:1)
          ELSE 
            GO TO 96
          END IF
          READ(UNIT=T2, FMT='(I5)', ERR=96) RESID1
          CHAIN2 = CHAIN1
          READ(UNIT=T3, FMT='(I5)', ERR=96) RESID2
        ELSE IF (VALUES.EQ.4) THEN
          READ(UNIT=LINE, FMT=*, ERR=96) T1, T2, T3, T4
C-------  Sanity check and clean read
          IF (LENSTR(T1).EQ.1) THEN
            CHAIN1 = T1(1:1)
          ELSE 
            GO TO 96
          END IF
          READ(UNIT=T2, FMT='(I5)', ERR=96) RESID1
          IF (LENSTR(T3).EQ.1) THEN
            CHAIN2 = T3(1:1)
          ELSE 
            GO TO 96
          END IF
          READ(UNIT=T4, FMT='(I5)', ERR=96) RESID2
        ELSE
          GO TO 96
        END IF   


C-------Write out what we have
80      WRITE(UNIT=10, FMT=94)CHAIN1, RESID1, CHAIN2, RESID2
1     CONTINUE
      
C-----Skip error messages 
      GO TO 99
C-----Format statements
93    FORMAT(A1)
94    FORMAT('RANGE  ''',A1,I4,'.'' ''',A1,I4,'.'' ALL')         
      
C-----Error messages
96    WRITE(6,*) 'Dirty BUSTER error'
97    WRITE(6,*) 'Cannot read TLS residue range'

99    RETURN
      END    

C-----------------------------------------------------------------------
C  This subroutine reads a clean BUSTER style TLS group selection and 
C  writes it out in Refmac format.    
C
C Format examples: 
C123456789012345678901234567890123456789012345678901234567890         
C LOCAL   
CREMARK   3    SET : { A|555 - A|792 } 
CREMARK   3    SET : { A|* }
CREMARK   3    SET : { A|885 - A|1009 }   
CREMARK   3    SET : { A|1010 - A|1017 }  
CREMARK   3    SET : { B|469 - B|549 B|553 - B|592 B|598 - B|683 }
CREMARK   3    SELECTION: {P|74:234 298:333}
C Not clean (4cod):
CREMARK   3    SELECTION: {B AND RESEQ 1270}, {D AND RESEQ 1270},                
CREMARK   3               {F AND RESEQ 1271}, {H AND RESEQ 1270}  
C-----------------------------------------------------------------------
      SUBROUTINE CBUSTER(LINE)    
     
      IMPLICIT NONE
C-----Declare variables and constants    
      CHARACTER LINE*(*), T1*6, T2*6, T3*6, T4*6, CHAIN1*1, CHAIN2*1,
     +          CLINE*200, SLINE*30, CHAINC*1
      INTEGER   DATCHK, VALUES, LENSTR, RESID1, RESID2, POS, POS2, I, J,
     +          K, PSTAR, TMP

C-----Initialise
      CHAINC = '*'
      
C-----Clean up the data line
C      PRINT*, LINE
C-----Make the line 'compact', and write to the cache
      CALL DUBSPA(LINE)
      CLINE =  LINE(:LENSTR(LINE))
      
C-----Go to dirty buster if there are no curly brackets
      IF (INDEX(LINE,'{').EQ.0) THEN
        CALL PHENIX(LINE)
        GO TO 999
      END IF 
      
C-----Split the lines by comma's
      DO 1, K=1, 10
        IF (LENSTR(CLINE).GT.0) THEN
          POS = INDEX(CLINE,',')
          IF (POS.NE.0) THEN
            LINE  = CLINE(1:POS-1)
            CLINE = CLINE(POS+1:)
          ELSE
            LINE  = CLINE
            CLINE = ''
          ENDIF  
C          PRINT*, CLINE
        ELSE  
          GO TO 999
        END IF
      
C-----  Remove '{'
        DO 10, I=1, 5
          POS = INDEX(LINE,'{')
          IF (POS.EQ.0) GO TO 11
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
10      CONTINUE

C-----  Remove '}'
11      DO 12, I=1, 5
          POS = INDEX(LINE,'}')
          IF (POS.EQ.0) GO TO 13
          LINE = LINE(1:POS-1)//' '//LINE(POS+1:LENSTR(LINE))
12      CONTINUE

C-----  Remove single quotes
13      DO 14, I=1, 10
          POS = INDEX(LINE,'''')
          IF (POS.EQ.0) GO TO 15
          LINE = LINE(1:POS-1)//LINE(POS+1:LENSTR(LINE))
14      CONTINUE

C-----  Fix newly introduced spaces
15      DO 16, I=1, 2
          POS = INDEX(LINE,'| ')
          IF (POS.EQ.0) GO TO 17
          LINE = LINE(1:POS-1)//'|'//LINE(POS+2:LENSTR(LINE))
16      CONTINUE      

C-----  Remove colons as range separators, but stop when used atom 
C       selectors
17      DO 18, I=1, 5
          POS = INDEX(LINE,':')
          IF (POS.EQ.0) GO TO 21
C         If the colon is an atom selctor, the next character will not
C         be a digit (this is not fully robust!)
          IF (INDEX('0123456789',LINE(POS+1:POS+1)).GT.0) THEN 
            LINE = LINE(1:POS-1)//'-'//LINE(POS+1:LENSTR(LINE))
          ELSE
            GO TO 996
          ENDIF
18      CONTINUE     

C-----  Insert spaces where needed
C     PRINT*, LINE
21      CALL MINCLR(LINE, 'D', 'F')  
        CALL MINCLR(LINE, 'L', 'P')
C        PRINT*, LINE

C-----  Check the amount of data on the line (should be n*3)
100     VALUES = DATCHK(LINE)
C        PRINT*, 'Val ', VALUES
        IF (MOD (VALUES,3).NE. 0) THEN
C         Check whether a full chain is selected
          PSTAR = INDEX (LINE, '|*')
          IF (PSTAR.NE.0) THEN
C           The whole chain is selected
            CHAIN1 = LINE(PSTAR-1:PSTAR-1)
            RESID1 = -10
            CHAIN2 = CHAIN1
            RESID2 = 9999
C           Write the range 
            WRITE(UNIT=10, FMT=994)CHAIN1, RESID1, CHAIN2, RESID2
C           Strip this part of the line
            LINE  = LINE(1:PSTAR-2)//LINE(PSTAR+2:LENSTR(LINE))
C           Loop back
            GO TO 100           
          ELSE
C           Try dirty Buster
            CALL BUSTER(LINE)
            GO TO 1
          END IF
        END IF

C-----  Loop over sets of 3 values. These should be ranges
        POS2 = 1
        DO 110, I=1, VALUES/3
C         Read 3 values and concatenate them
          READ(UNIT=LINE(POS2:LENSTR(LINE)), FMT=*, ERR=996) T1, T2, T3
          SLINE = T1//' '//T2//' '//T3
C          PRINT*, SLINE          
C         Calculate new value of POS (twice for cases that a range only
C         has a single residue)
          POS2 = INDEX(LINE, T3(1:LENSTR(T3))) + LENSTR(T3)
          TMP = INDEX(LINE(POS2:),T3(1:LENSTR(T3)))+LENSTR(T3) + POS2-1
C         TMP = the hit while searching the line from the end of the 
C             previous hit + the actual length of the hit + the bit of
C             the line up to the end of the first hit
          IF (TMP .NE. (LENSTR(T3) + POS2-1)) THEN
            POS2 = TMP
          END IF
C         Remove '|'
          DO 111, J=1, 3
            POS = INDEX(SLINE,'|')
            IF (POS.EQ.0) GO TO 112
            SLINE = SLINE(1:POS-1)//' '//SLINE(POS+1:LENSTR(SLINE))
111       CONTINUE

C       Remove '- '
112       DO 113, J=1, 2
C          PRINT*, SLINE
            POS = INDEX(SLINE,'- ')
            IF (POS.EQ.0) GO TO 120
            SLINE = SLINE(1:POS-1)//' '//SLINE(POS+2:LENSTR(SLINE))
113       CONTINUE  

C       Get the range (either in 2, 3 or in 4-value mode)
120       VALUES = DATCHK(SLINE)
C        PRINT*, '#'//SLINE//'#', VALUES
          IF (VALUES .EQ. 2) THEN
C         Use the chain cache if available
            IF (CHAINC .EQ. '*') THEN 
              GO TO 996
            ELSE
              READ(UNIT=SLINE, FMT=*, ERR=996) RESID1, RESID2
              CHAIN1 = CHAINC
              CHAIN2 = CHAIN1       
            END IF           
          ELSE IF (VALUES .EQ. 3) THEN
            READ(UNIT=SLINE, FMT=*, ERR=996) T1, T2, T3
C-------Sanity check and clean read
            IF (LENSTR(T1).EQ.1) THEN
              CHAIN1 = T1(1:1)
              CHAIN2 = CHAIN1
            ELSE 
             GO TO 996
            END IF
            READ(UNIT=T2, FMT='(I5)', ERR=996) RESID1
            READ(UNIT=T3, FMT='(I5)', ERR=996) RESID2
          ELSE IF (VALUES .EQ. 4) THEN
            READ(UNIT=SLINE, FMT=*, ERR=996) T1, T2, T3, T4
C-------Sanity check and clean read
            IF (LENSTR(T1).EQ.1) THEN
              CHAIN1 = T1(1:1)
            ELSE 
             GO TO 996
            END IF
            READ(UNIT=T2, FMT='(I5)', ERR=996) RESID1
            IF (LENSTR(T3).EQ.1) THEN
              CHAIN2 = T3(1:1)
            ELSE 
             GO TO 996
            END IF
            READ(UNIT=T4, FMT='(I5)', ERR=996) RESID2
          ELSE
            GO TO 996
          END IF

C         Write the range 
          WRITE(UNIT=10, FMT=994)CHAIN1, RESID1, CHAIN2, RESID2

C         Set the cache
          CHAINC = CHAIN1      

110     CONTINUE
1     CONTINUE

      
C-----Skip error messages 
      GO TO 999
C-----Format statements
993   FORMAT(A1)
994   FORMAT('RANGE  ''',A1,I4,'.'' ''',A1,I4,'.'' ALL')         
      
C-----Error messages
996   WRITE(6,*) 'Clean BUSTER error'
997   WRITE(6,*) 'Cannot read TLS residue range'

999   RETURN
      END    

C-----------------------------------------------------------------------
C  This subroutine extracts TLS parameters from the PDB coordinates
C  and writes them out to TLSOUT in REFMAC5 format      
C-----------------------------------------------------------------------
      SUBROUTINE MAKTLS(CHAINS, CHAINI, CFIRST, CLAST, MACMOL, TLSOUT,
     +           LINK_CACHE, NLINK, TLSRES, NTLSRES) 
     
      IMPLICIT NONE
C-----Declare variables and constants
      INTEGER   CHAINS, CFIRST(CHAINS), CLAST(CHAINS), STATUS, NTLSRES,
     +          BONDS 
      LOGICAL   MACMOL(CHAINS)
      CHARACTER TLSOUT*255, CHAINI(CHAINS)*1, LINK_CACHE(NLINK)*80,
     +          TLSRES(NTLSRES)*3, TMPRES*3, CHOP*3, ADDRES(20)*5  
      INTEGER   I,J,K,L, GROUPS, NLINK, NADD

C-----Initialise
      GROUPS = 0
      BONDS  = 0
      NADD   = 1
      ADDRES(1) = '     '
C      PRINT*, CHAINI, CFIRST, CLAST
C      PRINT*, CHAINS, MACMOL
C      PRINT*, NLINK
C      PRINT*, LINK_CACHE
C      PRINT*, CHAINI(1),CFIRST(1),CHAINI(1),CLAST(1), MACMOL(1)

C-----Open output file
      OPEN(UNIT=10, FILE=TLSOUT, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) TLSOUT, 'cannot be created!'
        GO TO 99
      END IF
C-----Write TLS parameters to file
      DO 20, J=1,CHAINS
        IF ((MACMOL(J).EQV..TRUE.).AND.((CLAST(J)-CFIRST(J)).GT.4)) THEN
          GROUPS = GROUPS+1
          WRITE(UNIT=10, FMT=92)'TLS    Group', GROUPS    
C123456789012345678901234567890123456789012345678901234567890     
CREMARK   3    RESIDUE RANGE :   A    -6        A   447
          WRITE(UNIT=10, FMT=94)CHAINI(J),CFIRST(J),CHAINI(J),CLAST(J)
C          PRINT*, CHAINI(J),CFIRST(J),CHAINI(J),CLAST(J)
          
C         Append compounds exclusively LINKed to this chain 
          NADD = 0
          DO 25, I=1, NLINK
C           Get all residues involved with this chain 
            IF (LINK_CACHE(I)(22:22).EQ.CHAINI(J)) THEN
C             Skip residues with insertion codes
              IF (LINK_CACHE(I)(57:57).NE.' ') GO TO 25
C             Skip residues we already have
              DO 26, K=1, NADD  
                IF (LINK_CACHE(I)(52:56).EQ.ADDRES(K)) GO TO 25
26            CONTINUE  
C             Only residues of the right type   
              TMPRES = LINK_CACHE(I)(48:50)
              TMPRES = CHOP(TMPRES, 2)
              DO 30, K=1, NTLSRES
C                PRINT*,TMPRES
                IF (TMPRES.EQ.TLSRES(K)) THEN 
                  BONDS = 1
C                 Is the residue exclusively LINKed to one chain
                  DO 40, L=1,NLINK
                    IF(LINK_CACHE(L)(52:57).EQ.LINK_CACHE(I)(52:57))THEN
                      IF (LINK_CACHE(L)(22:22).EQ.
     +                    LINK_CACHE(I)(22:22)) THEN
C                       Additional LINK to the same chain
                        BONDS = BONDS + 1
                      ELSE
C                       LINK to a different chain. Don't add the residue
                        GO TO 25
                      END IF  
               ELSE IF(LINK_CACHE(L)(22:27).EQ.LINK_CACHE(I)(52:57))THEN
                      IF (LINK_CACHE(L)(52:52).EQ.
     +                    LINK_CACHE(I)(22:22)) THEN
C                       Additional LINK to the same chain
                        BONDS = BONDS + 1
                      ELSE
C                       LINK to a different chain. Don't add the residue
                        GO TO 25
                      END IF   
                    END IF
40                CONTINUE       
C                 Get the residues LINKed more than once
                  IF (BONDS .GE. 2) THEN
                    NADD = NADD +1
                    ADDRES(NADD) = LINK_CACHE(I)(52:56)   
                  END IF        
                END IF  
30            CONTINUE
            END IF
C           Now do swapped LINKs
            IF (LINK_CACHE(I)(52:52).EQ.CHAINI(J)) THEN
C             Skip residues with insertion codes
              IF (LINK_CACHE(I)(27:27).NE.' ') GO TO 25
C             Skip residues we already have
              DO 27, K=1, NADD
                IF (LINK_CACHE(I)(22:26).EQ.ADDRES(K)) GO TO 25
27            CONTINUE
C             Only residues of the right type
              TMPRES = LINK_CACHE(I)(18:20)
              TMPRES = CHOP(TMPRES, 2)
              DO 50, K=1, NTLSRES
                IF (TMPRES.EQ.TLSRES(K)) THEN
                  BONDS = 1
C                 Is the residue exclusively LINKed to one chain
                  DO 60, L=1,NLINK
                    IF(LINK_CACHE(L)(22:27).EQ.LINK_CACHE(I)(22:27))THEN
                      IF (LINK_CACHE(L)(52:52).EQ.
     +                    LINK_CACHE(I)(52:52)) THEN
C                       Additional LINK to the same chain
                        BONDS = BONDS + 1
                      ELSE
C                       LINK to a different chain. Don't add the residue
                        GO TO 25
                      END IF
               ELSE IF(LINK_CACHE(L)(52:57).EQ.LINK_CACHE(I)(22:27))THEN
                      IF (LINK_CACHE(L)(22:22).EQ.
     +                    LINK_CACHE(I)(52:52)) THEN
C                       Additional LINK to the same chain
                        BONDS = BONDS + 1
                      ELSE
C                       LINK to a different chain. Don't add the residue
                        GO TO 25
                      END IF
                    END IF
60                CONTINUE
C                 Get the residues LINKed more than once
                  IF (BONDS .GE. 2) THEN
                    NADD = NADD +1
                    ADDRES(NADD) = LINK_CACHE(I)(22:26)
                  END IF
                END IF
50            CONTINUE
 
                
            END IF
25        CONTINUE  
C         Add the LINKed residues
          DO 85, K=1, NADD
            WRITE(UNIT=10, FMT=95) ADDRES(K), ADDRES(K) 
85        CONTINUE  
C         Close the TLS group
          WRITE(10,*) ' '
        END IF 
20    CONTINUE
      
C-----Skip error messages 
      GO TO 99
C-----Format statements
91    FORMAT(A,I3)
92    FORMAT(A12,1X,I3)
94    FORMAT('RANGE  ''',A1,I4,'.'' ''',A1,I4,'.'' ALL')
95    FORMAT('RANGE  ''',A5,'.'' ''',A5,'.'' ALL')
      
C-----Wrap-up
99    WRITE(6,91) 'TLS groups created  : ', GROUPS
      IF (GROUPS.EQ.0) THEN
        CLOSE(10,STATUS='DELETE')
      ELSE  
        CLOSE(10)
      END IF  
      RETURN
      END      

C-----------------------------------------------------------------------
C  This subroutine checks negative selectors in TLS groups      
C-----------------------------------------------------------------------
      SUBROUTINE TLSNOT(LONG,USETLS) 

      IMPLICIT NONE
C-----Declare variables and constants
      LOGICAL   USETLS
      CHARACTER LONG*400
      INTEGER   PROBE

C     Set USETLS to false and make it TRUE if a fix is applied
      USETLS = .FALSE.
  
C     Remove negative sections for hydrogens or waters
      PROBE = INDEX(LONG,' and not element H')
      IF (PROBE.NE.0) THEN
        LONG   = LONG(1:PROBE-1)//' '//LONG(PROBE+18:)
        USETLS = .TRUE. 
      END IF
      
      PROBE = INDEX(LONG,' AND NOT (ELEMENT H)')
      IF (PROBE.NE.0) THEN
        LONG   = LONG(1:PROBE-1)//' '//LONG(PROBE+18:)
        USETLS = .TRUE. 
      END IF      

      PROBE = INDEX(LONG,' AND NOT ELEMENT H')
      IF (PROBE.NE.0) THEN
        LONG   = LONG(1:PROBE-1)//' '//LONG(PROBE+18:)
        USETLS = .TRUE. 
      END IF

      PROBE = INDEX(LONG,' and not resname HOH')
      IF (PROBE.NE.0) THEN
        LONG   = LONG(1:PROBE-1)//' '//LONG(PROBE+20:)
        USETLS = .TRUE. 
      END IF

      PROBE = INDEX(LONG,' AND NOT RESNAME HOH')
      IF (PROBE.NE.0) THEN
        LONG   = LONG(1:PROBE-1)//' '//LONG(PROBE+20:)
        USETLS = .TRUE. 
      END IF

      PROBE = INDEX(LONG,' AND NOT WATER')
      IF (PROBE.NE.0) THEN
        LONG   = LONG(1:PROBE-1)//' '//LONG(PROBE+14:)
        USETLS = .TRUE. 
      END IF

      RETURN
      END      
     
      
C ----------------------------------------------------------------------
C     This subroutine extracts data from the from PDB header.
C     CODE TAKEN FROM SFCHECK and substantially adapted to deal with
C     missing values and value precedence
C ----------------------------------------------------------------------
      SUBROUTINE RDREMARK(LINE)

C-----Declare the variables and parameters       
      CHARACTER LINE*80, CH36*36,CH16*16,NUMSPA*11
      REAL      GETVAL, TEMP
      PARAMETER (NUMSPA='0123456789 ')
      INTEGER   I, J, NUMCHK, NULL

      COMMON /HEADER/   RFAC_ALL  ,RFAC_OBS  ,RFAC_CUT,
     +                  RFREE_ALL ,RFREE_OBS ,RFREE_CUT,
     +                  BWILS     ,RESO3     ,COMPLET,
     +                  MASPRO    ,MASION    ,MASSHR,
     +                  RPROG
      REAL    RFAC_ALL
      REAL    RFAC_OBS
      REAL    RFAC_CUT
      REAL    RFREE_ALL
      REAL    RFREE_OBS
      REAL    RFREE_CUT
      REAL    BWILS
      REAL    RESO3
      REAL    COMPLET
      REAL    MASPRO
      REAL    MASION
      REAL    MASSHR
      CHARACTER RPROG*40


C-------------End of variable declaration------------------------------C

C-----Initialise
      TEMP = 0.00

C-----Skip line when there is no remark number or the line contains a 
C     'NULL' value
      NULL   = INDEX(LINE, 'NULL')
      IF (NULL.NE.0) THEN
C       WRITE(6,*) LINE, 'Skipped'
        GO TO 100
      END IF
      DO 10, J=8,10
        NUMCHK = INDEX(NUMSPA,LINE(J:J))
	IF (NUMCHK.EQ.0) THEN
C          WRITE(6,*) LINE, 'Skipped'
	  GO TO 100
	END IF
10    CONTINUE

      IF(LINE(47:47).EQ.':'.OR.LINE(48:48).EQ.':'.OR.
     +   LINE(55:55).EQ.':'.OR.LINE(58:58).EQ.':'.OR.
     +   LINE(44:44).EQ.':'.OR.LINE(56:56).EQ.':'.OR.
     +   LINE(26:26).EQ.':'.OR.LINE(45:45).EQ.':'.OR.
     +   LINE(33:33).EQ.':') THEN
        GO TO 200
      END IF
      READ(UNIT=LINE,FMT=999,ERR=100) I,CH36
C      WRITE(UNIT=6,FMT='(A36)') CH36
C123456789 123456789 123456789 123456789 123456789 123456789 1
CREMARK   3   R VALUE     (WORKING + TEST SET, NO CUTOFF) :
CREMARK   3   R VALUE     (WORKING + TEST SET)  : 0.163
C1234567123 123456789 123456789 123456789 123456 
CREMARK   3   R VALUE                    0.203 
CREMARK   3   FREE R VALUE               0.31
CREMARK   3   PROGRAM                    X-PLOR 

      IF(I.EQ.3.AND.CH36(1:10).EQ.'  R VALUE ') THEN
        CH16=CH36(29:34)
        RFAC_OBS=GETVAL(CH16)
      ELSE IF(I.EQ.3.AND.CH36(1:16).EQ.'  FREE R VALUE  ') THEN
        CH16=CH36(29:34)
        RFREE_OBS=GETVAL(CH16)
      ELSE IF(I.EQ.3.AND.CH36(1:10).EQ.'  PROGRAM ') THEN
        RPROG   = LINE(41:80)
      ENDIF
      RETURN

 200  CONTINUE

C123456789 123456789 123456789 123456789 123456789 123456789
C1234567123 123456789 123456789 12345678
C
C Refinement using X-PLOR
C
cREMARK   3   PROGRAM     : X-PLOR X.X
cREMARK   3   R VALUE            (WORKING SET) :
cREMARK   3   FREE R VALUE                     :
c
c Refinement using NUCLSQ
c
cREMARK   3   PROGRAM     : NUCLSQ
cREMARK   3   R VALUE     (WORKING + TEST SET) :
cREMARK   3   R VALUE            (WORKING SET) :
cREMARK   3   R VALUE   (WORKING + TEST SET, NO CUTOFF) :
cREMARK   3   R VALUE          (WORKING SET, NO CUTOFF) :
cREMARK   3   FREE R VALUE                  (NO CUTOFF) :
C
C Refinement using PROLSQ
C
CREMARK   3   PROGRAM     : PROLSQ
CREMARK   3   R VALUE     (WORKING + TEST SET) :
CREMARK   3   R VALUE            (WORKING SET) :
CREMARK   3   FREE R VALUE                     :
CREMARK   3   R VALUE   (WORKING + TEST SET, NO CUTOFF) :
CREMARK   3   R VALUE          (WORKING SET, NO CUTOFF) :
CREMARK   3   FREE R VALUE                  (NO CUTOFF) :
C
C Refinement using SHELXL-96
C
CREMARK   3  PROGRAM     : SHELXL-XX
CREMARK   3  R VALUE   (WORKING + TEST SET, NO CUTOFF) :
CREMARK   3  R VALUE          (WORKING SET, NO CUTOFF) :
CREMARK   3  FREE R VALUE                  (NO CUTOFF) :
C
C Refinement using TNT
C
CREMARK   3   PROGRAM     : TNT XXX
CREMARK   3   R VALUE     (WORKING + TEST SET) :
CREMARK   3   R VALUE            (WORKING SET) :
CREMARK   3   FREE R VALUE                     :
CREMARK   3   R VALUE   (WORKING + TEST SET, NO CUTOFF) :
CREMARK   3   R VALUE          (WORKING SET, NO CUTOFF) :
CREMARK   3   FREE R VALUE                  (NO CUTOFF) :
C
C EM data
CREMARK   3 REFINEMENT.                                                          
CREMARK   3   SOFTWARE PACKAGES      : FREALIGN                                  
CREMARK   3   RECONSTRUCTION SCHEMA  : NULL                                      
CREMARK   3                                                                      
CREMARK   3 EM MAP-MODEL FITTING AND REFINEMENT                                  
CREMARK   3   PDB ENTRY                    : NULL                                
CREMARK   3   REFINEMENT SPACE             : NULL                                
CREMARK   3   REFINEMENT PROTOCOL          : NULL                                
CREMARK   3   REFINEMENT TARGET            : NULL                                
CREMARK   3   OVERALL ANISOTROPIC B VALUE  : NULL                                
CREMARK   3                                                                      
CREMARK   3 FITTING PROCEDURE : NULL                                             
CREMARK   3                                                                      
CREMARK   3 EM IMAGE RECONSTRUCTION STATISTICS                                   
CREMARK   3   NOMINAL PIXEL SIZE (ANGSTROMS)    : 0.637                          
CREMARK   3   ACTUAL PIXEL SIZE  (ANGSTROMS)    : 0.637                          
CREMARK   3   EFFECTIVE RESOLUTION (ANGSTROMS)  : 2.2                            
CREMARK   3   NUMBER OF PARTICLES               : 41123                          
cREMARK   3   CTF CORRECTION METHOD             : EACH PARTICLE   
C
C Missing values
C
CREMARK   3   R VALUE          (WORKING SET, F>4SIG(F)) : 0.114       
CREMARK   3   FREE R VALUE                  (F>4SIG(F)) : 0.190
CREMARK   3  R VALUE   (WORKING + TEST SET, NO CUTOFF) : 0.1277              
CREMARK   3  R VALUE          (WORKING SET, NO CUTOFF) : NULL     
CREMARK   3   R VALUE          (WORKING SET, NO CUTOFF) : 0.143
CREMARK   3   R VALUE      (WORKING + TEST SET) : 0.1934   
CREMARK   3   R VALUE     (WORKING + TEST SET)  : 0.163                      
CREMARK   3   R VALUE             (WORKING SET) : 0.1904   
CREMARK   3   FREE R VALUE                      : 0.2456 
CREMARK   3   R VALUE            (WORKING SET)  : 0.193  
CREMARK   3   R VALUE           (WORKING SET)   : 0.1771    
C 
C Resolution
CREMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :   2.00
      IF(LINE(14:48).EQ.
     +'RESOLUTION RANGE HIGH (ANGSTROMS) :') THEN
        CH16=LINE(49:64)
        RESO3=GETVAL(CH16) 
        
      ELSE IF(LINE(14:48).EQ.
     +'EFFECTIVE RESOLUTION (ANGSTROMS)  :') THEN
        CH16=LINE(49:64)
        RESO3=GETVAL(CH16)  

      ELSE IF(LINE(14:47).EQ.
     +'R VALUE                          :') THEN
        CH16=LINE(48:63)
        RFAC_OBS=GETVAL(CH16)  

      ELSE IF(LINE(14:47).EQ.      
     +'R VALUE     (WORKING + TEST SET) :') THEN
        CH16=LINE(48:63)
        RFAC_ALL=GETVAL(CH16)

      ELSE IF(LINE(14:47).EQ.
     +'R VALUE            (WORKING SET) :') THEN
        CH16=LINE(48:63)
C-------Do not assign a new value if one already exists
        IF ((RFAC_OBS.GT.0.00).AND.(RFAC_OBS.LT.0.999))GO TO 100
        RFAC_OBS=GETVAL(CH16)

      ELSE IF(LINE(14:48).EQ.      
     +'R VALUE      (WORKING + TEST SET) :') THEN
        CH16=LINE(49:63)
        RFAC_ALL=GETVAL(CH16)

      ELSE IF(LINE(14:48).EQ.
     +'R VALUE     (WORKING + TEST SET)  :') THEN
        CH16=LINE(49:63)
        RFAC_ALL=GETVAL(CH16)

      ELSE IF(LINE(14:48).EQ.
     +'R VALUE             (WORKING SET) :') THEN
        CH16=LINE(49:63)
        RFAC_OBS=GETVAL(CH16)

      ELSE IF(LINE(14:48).EQ.
     +'R VALUE            (WORKING SET)  :') THEN
        CH16=LINE(49:63)
        RFAC_OBS=GETVAL(CH16)
      ELSE IF(LINE(14:48).EQ.
     +'R VALUE           (WORKING SET)   :') THEN
        CH16=LINE(49:63)
        RFAC_OBS=GETVAL(CH16)

      ELSE IF(LINE(14:48).EQ.
     +'FREE R VALUE                      :') THEN
        CH16=LINE(49:63)
        RFREE_OBS=GETVAL(CH16)

      ELSE IF(LINE(14:47).EQ.
     +'FREE R VALUE                     :') THEN
        CH16=LINE(48:63)
C-------Do not assign a new value if one already exists
        IF ((RFREE_OBS.GT.0.00).AND.(RFREE_OBS.LT.0.999))GO TO 100
        RFREE_OBS=GETVAL(CH16)
      ELSE IF(LINE(14:47).EQ.
     +'FREE R VALUE TEST SET SELECTION  :') THEN
        CH16=LINE(48:63)
C-------Do not assign a new value if one already exists
	IF ((RFREE_OBS.GT.0.00).AND.(RFREE_OBS.LT.0.999))GO TO 100
        RFREE_OBS=GETVAL(CH16)

      ELSE IF(LINE(14:56).EQ.
     +'R VALUE          (WORKING SET, NO CUTOFF) :') THEN
        CH16=LINE(57:69)
        RFAC_OBS=GETVAL(CH16)
	
      ELSE IF(LINE(14:56).EQ.
     +'R VALUE          (WORKING SET, F>4SIG(F)) :') THEN
        CH16=LINE(57:69)
C-------Do not assign a new value if one already exists
	IF ((RFAC_CUT.GT.0.00).AND.(RFAC_CUT.LT.0.999))GO TO 100
	RFAC_CUT=GETVAL(CH16)

      ELSE IF(LINE(14:56).EQ.
     +'R VALUE   (WORKING + TEST SET, NO CUTOFF) :') THEN
        CH16=LINE(57:69)
C-------Make sure a sensible value is used
        TEMP=GETVAL(CH16)
        IF ((TEMP.GT.0.00).AND.(TEMP.LT.0.999)) THEN
	  RFAC_ALL=TEMP
        END IF 

      ELSE IF(LINE(14:56).EQ.
     +'R VALUE   (WORKING + TEST SET, F>4SIG(F)) :') THEN
        CH16=LINE(57:69)
C-------Do not assign a new value if one already exists
	IF ((RFAC_ALL.GT.0.00).AND.(RFAC_ALL.LT.0.999))GO TO 100
        RFAC_ALL=GETVAL(CH16)

      ELSE IF(LINE(13:55).EQ.
     +'R VALUE   (WORKING + TEST SET, NO CUTOFF) :') THEN
        CH16=LINE(56:71)
        RFAC_ALL=GETVAL(CH16)

      ELSE IF(LINE(13:55).EQ.
     +'R VALUE          (WORKING SET, NO CUTOFF) :') THEN
        CH16=LINE(56:71)
        RFAC_OBS=GETVAL(CH16)

      ELSE IF(LINE(13:55).EQ.
     +'FREE R VALUE                  (NO CUTOFF) :') THEN
        CH16=LINE(56:71)
        RFREE_OBS=GETVAL(CH16)

      ELSE IF(LINE(14:56).EQ.
     +'FREE R VALUE                  (NO CUTOFF) :') THEN
        CH16=LINE(57:72)
        RFREE_OBS=GETVAL(CH16)

      ELSE IF(LINE(14:56).EQ.
     +'FREE R VALUE                  (F>4SIG(F)) :') THEN
        CH16=LINE(57:72)
C-------Do not assign a new value if one already exists
	IF ((RFREE_CUT.GT.0.00).AND.(RFREE_CUT.LT.0.999))GO TO 100
	RFREE_CUT=GETVAL(CH16)

      ELSE IF(LINE(14:26).EQ.
     +'PROGRAM     :') THEN
        RPROG   = LINE(28:67)

C-----Extract the Wilson B-factor
C123456789 123456789 123456789 123456789 123456789 123456789
CREMARK   3  WILSON B VALUE (FROM FCALC, A**2) : 2.300
CREMARK   3   FROM WILSON PLOT           (A**2) : 2.40

      ELSE IF(LINE(13:47).EQ.
     +'WILSON B VALUE (FROM FCALC, A**2) :') THEN
        CH16=LINE(48:63)
        BWILS   = GETVAL(CH16)
	IF (BWILS.LT.0.0) THEN
	  BWILS = 0.0
	END IF    

      ELSE IF(LINE(14:48).EQ.
     +'FROM WILSON PLOT           (A**2) :') THEN
        CH16=LINE(49:64)
        BWILS   = GETVAL(CH16)
	IF (BWILS.LT.0.0) THEN
	  BWILS = 0.0
	END IF  

C-----Extract the completeness
C123456789 123456789 123456789 123456789 123456789 123456789
CREMARK   3   COMPLETENESS (WORKING+TEST)   (%) : 95.3
CREMARK   3   COMPLETENESS FOR RANGE        (%) : 99.8  
      ELSE IF(LINE(14:48).EQ.
     +'COMPLETENESS (WORKING+TEST)   (%) :') THEN
        CH16=LINE(49:64)
        COMPLET   = GETVAL(CH16)
	IF (COMPLET.LT.0.0) THEN
	  COMPLET = 0.0
	END IF    

      ELSE IF(LINE(14:48).EQ.
     +'COMPLETENESS FOR RANGE        (%) :') THEN
        CH16=LINE(49:64)
        COMPLET   = GETVAL(CH16)
	IF (COMPLET.LT.0.0) THEN
	  COMPLET = 0.0
	END IF  
 
C-----Extract the mask parameters
C Refmac
CREMARK   3   VDW PROBE RADIUS   : 1.40                                          
CREMARK   3   ION PROBE RADIUS   : 0.80                                          
CREMARK   3   SHRINKAGE RADIUS   : 0.80  
CREMARK   3   VDW PROBE RADIUS   :   1.10
CREMARK   3   ION PROBE RADIUS   :   0.70
CREMARK   3   SHRINKAGE RADIUS   :   0.70
C Phenix
CREMARK   3   SOLVENT RADIUS     : 1.11                                          
CREMARK   3   SHRINKAGE RADIUS   : 0.90 
      ELSE IF(LINE(14:33).EQ.
     +'VDW PROBE RADIUS   :') THEN
        CH16=LINE(34:44)
        MASPRO   = GETVAL(CH16)

      ELSE IF(LINE(14:33).EQ.
     +'ION PROBE RADIUS   :') THEN
        CH16=LINE(34:44)
        MASION   = GETVAL(CH16)
 
      ELSE IF(LINE(14:33).EQ.
     +'SOLVENT RADIUS     :') THEN
        CH16=LINE(34:44)
        MASPRO   = GETVAL(CH16)
        MASION   = MASPRO 

      ELSE IF(LINE(14:33).EQ.
     +'SHRINKAGE RADIUS   :') THEN
        CH16=LINE(34:44)
        MASSHR   = GETVAL(CH16)

      ENDIF

  100 CONTINUE
  
C-----Format statements
999   FORMAT(7X,I3,1X,A36)

      RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine prepares coefficients for calculating reflection 
C  resolution
C  This subroutine was adapted from code in the SFTOOLS program, written
C  by Bart Hazes, which was based on an earlier implementation taken 
C  from the Groningen BIOMOL protein crystallography software package.
C-----------------------------------------------------------------------
      SUBROUTINE RESCNST(AAXIS, BAXIS, CAXIS, ALPHA, BETA, GAMMA, COEFS)
C///
C     THE RESOLUTION CAN BE CALCULATED AS FOLLOWS: 
C
C     DSTAR-SQUARED = IH*IH*A11 + IH*IK*A21 + IH*IL*A31
C                   + IK*IK*A22 + IK*IL*A32 + IL*IL*A33
C     SO:       RES = 1.0/SQRT(DSTAR_SQUARED)
C
C     IN THIS SUBROUTINE THE COEFFICIENTS A11, A21, A31, A22, A32, A33
C     ARE CALCULATED
C
C     this subroutine is a modified version of dcalc in cfblib
C
C\\\

      IMPLICIT NONE
C-----Declare variables
      DOUBLE PRECISION ca,sa,cb,sb,cg,sg,cast,cbst,cgst,sast,sbst,sgst
      DOUBLE PRECISION ast,bst,cst
      DOUBLE PRECISION COEFS(6), DEG2RAD
      REAL AAXIS, BAXIS, CAXIS, ALPHA, BETA, GAMMA
      
C-----Initialise      
      DEG2RAD = ATAN(1.0)/45.0
      
      ca   = COS(DEG2RAD*ALPHA)
      sa   = SIN(DEG2RAD*ALPHA)
      cb   = COS(DEG2RAD*BETA)
      sb   = SIN(DEG2RAD*BETA)
      cg   = COS(DEG2RAD*GAMMA)
      sg   = SIN(DEG2RAD*GAMMA)
      cast = (cb*cg - ca)/ (sb*sg)
      cbst = (cg*ca - cb)/ (sg*sa)
      cgst = (ca*cb - cg)/ (sa*sb)
      sast = SQRT(1.0 - cast*cast)
      sbst = SQRT(1.0 - cbst*cbst)
      sgst = SQRT(1.0 - cgst*cgst)
      ast  = 1.0/(AAXIS*sb*sgst)
      bst  = 1.0/(BAXIS*sg*sast)
      cst  = 1.0/(CAXIS*sa*sbst)
      COEFS(1)  = ast*ast
      COEFS(2)  = 2.0*ast*bst*cgst
      COEFS(3)  = 2.0*ast*cst*cbst
      COEFS(4)  = bst*bst
      COEFS(5)  = 2.0*bst*cst*cast
      COEFS(6)  = cst*cst

      RETURN
      END
C-----------------------------------------------------------------------
C  This subroutine writes out a residue list for PDB_REDO tools
C-----------------------------------------------------------------------
      SUBROUTINE RES_LIST(CNT, CHN, RES, INS)

      IMPLICIT NONE
C-----Declare variables
      INTEGER   CNT, RES(CNT), POS, I, LENSTR
      CHARACTER CHN(CNT)*1, INS(CNT)*1
      CHARACTER OUTPUT*(7*CNT+1), DP*1, NUMBER*15, CHOP*15

      
C-----Initialise      
      DP = ':'
      OUTPUT = DP
      POS = 2

C-----Write empty line and finish if there are no residues in the list
      IF (CNT.EQ.0) THEN
        WRITE (UNIT=9, FMT='(A)')''
        GO TO 101
      END IF

C-----Compile output line (depends on residue number character count)
      DO 10, I=1, CNT
C       Convert the residue number to text and then write it out at 
C       the appropriate length 
        WRITE(UNIT=NUMBER,FMT=*) RES(I)
        NUMBER = CHOP(NUMBER,15)
        WRITE(UNIT=OUTPUT(POS:),FMT=99)CHN(I), NUMBER(:LENSTR(NUMBER)), 
     +    INS(I), DP
        POS = POS + LENSTR(NUMBER) + 3
10    CONTINUE

C-----Write out the results
      WRITE(UNIT=9, FMT='(A)') OUTPUT(1:POS-1)

C-----Format statements
99    FORMAT (A1,A,A1,A1)

101   RETURN
      END
C-----------------------------------------------------------------------
C  This function gives the length of a string, accounting for spaces
C  at the end of the line
C-----------------------------------------------------------------------
      INTEGER FUNCTION LENSTR(STRING)
      IMPLICIT NONE
      INTEGER   I,N
      CHARACTER STRING*(*)
C -------------------------------------
      LENSTR=LEN(STRING)
      DO 10 I=LENSTR,1,-1
        IF(STRING(LENSTR:LENSTR).NE.' ') GO TO 20      
        LENSTR=LENSTR-1
   10 CONTINUE
   20 RETURN
      END
      
C-----------------------------------------------------------------------
C  This function extracts a real number for a string. 0.0 is returned
C  in case of an empty or 'NULL' string
C-----------------------------------------------------------------------
      REAL FUNCTION GETVAL(LINE)
      IMPLICIT NONE
      
C-----Declare variables
      REAL      VAL,NUMBER
      INTEGER   LENSTR
      CHARACTER LINE*(*)
      
C-----Initialise
      GETVAL = 0.0
      NUMBER = 0.0
      
C-----Return if the lenght of the string is 0 or if the string 
C     contains 'NULL'     
      IF(LENSTR(LINE).LE.0) RETURN
      IF(INDEX(LINE, 'NULL').NE.0) RETURN

      READ(UNIT=LINE, FMT=*, ERR= 15) NUMBER
      GETVAL = NUMBER

15    RETURN
      END 
C-----------------------------------------------------------------------
C  This function calculates the resolution from HKL indices and
C  precalculated coefficients
C  This function was adapted from code in the SFTOOLS program, written
C  by Bart Hazes, which was based on an earlier implementation taken 
C  from the Groningen BIOMOL protein crystallography software package.
C-----------------------------------------------------------------------
      REAL FUNCTION HKLRES(INDEXH, INDEXK, INDEXL, COEFS)
      IMPLICIT NONE
      
C-----Declare varibles and constants
      DOUBLE PRECISION COEFS(6), TMPRES
      INTEGER INDEXH, INDEXK, INDEXL

C-----Calculate resolution
      TMPRES = INDEXH*INDEXH*COEFS(1) + INDEXH*INDEXK*COEFS(2) +
     +         INDEXH*INDEXL*COEFS(3) + INDEXK*INDEXK*COEFS(4) +
     +         INDEXK*INDEXL*COEFS(5) + INDEXL*INDEXL*COEFS(6)
      
C----Convert to Angstroms     
      HKLRES = (1.0/SQRT(TMPRES))

      RETURN
      END
 
      
C-----------------------------------------------------------------------
C  This function moves starting blanks from the head to the tail of a 
C  line 
C-----------------------------------------------------------------------
      CHARACTER*(*) FUNCTION CHOP(LINE, MAXCUT)
      IMPLICIT NONE

C-----Variables:
C     LINE    Text to be chopped
C     I       Just loop iterator
C     MAXCUT  The maximum number of leading spaces to be CHOP-ped
C     LANSTR  The return of the function LENSTR
C     LINLEN  The length of 'LINE' without the trailing spaces
      CHARACTER LINE*(*)
      INTEGER   I, MAXCUT, LENSTR, LINLEN

      LINLEN = LENSTR(LINE)
      DO 10, I=1, MAXCUT
        IF (LINE(1:1).NE.' ') THEN
          GO TO 20
        ELSE
          LINE(1:(LINLEN-1)) = LINE(2:LINLEN)
          LINE(LINLEN:LINLEN) = ' '
        END IF
10    CONTINUE
20    CHOP = LINE
      RETURN 
      END      
      
C-----------------------------------------------------------------------
C  This function replaces spaces and quotes by underscores but not the 
C  'fillers' at the end of a string
C-----------------------------------------------------------------------
      CHARACTER*(*) FUNCTION NOSPCE(LINE)
      IMPLICIT NONE
      CHARACTER LINE*(*)
      INTEGER   I, SPACE, LENSTR, LINLEN

      LINLEN = LENSTR(LINE)
      DO 10, I=1, LINLEN
        SPACE = INDEX(LINE, ' ')
        
        IF ((SPACE.EQ.0).OR.(SPACE.GT.LINLEN)) THEN
          GO TO 11
        ELSE
          LINE(SPACE:SPACE) = '_'
        END IF
10    CONTINUE
11    DO 20, I=1, LINLEN
        SPACE = INDEX(LINE, "'")
        
        IF ((SPACE.EQ.0).OR.(SPACE.GT.LINLEN)) THEN
          GO TO 30
        ELSE
          LINE(SPACE:SPACE) = '_'
        END IF
20    CONTINUE
30    NOSPCE=LINE     
      RETURN 
      END
C-----------------------------------------------------------------------
C  This subroutine replaces multiple spaces by single spaces 
C-----------------------------------------------------------------------
      SUBROUTINE DUBSPA(LINE)
      IMPLICIT NONE
      CHARACTER LINE*(*)
      INTEGER   I, DSPACE, LENSTR, LINLEN

      LINLEN = LENSTR(LINE)
      DO 10, I=1, LINLEN
        DSPACE = INDEX(LINE, '  ')
        IF ((DSPACE.EQ.0).OR.(DSPACE.GT.LINLEN)) THEN
          GO TO 20
        ELSE
          LINE = LINE(1:DSPACE)//LINE(DSPACE+2:)
        END IF
10    CONTINUE
     
20    RETURN 
      END
C-----------------------------------------------------------------------
C  Converts several keywords to uppercase.
C-----------------------------------------------------------------------
      SUBROUTINE UPCASE(LINE)
      IMPLICIT NONE
      CHARACTER LINE*(*)
      INTEGER   I, POS, LENSTR, LINLEN, FIXED

      LINLEN = LENSTR(LINE)
      DO 10, I=1, LINLEN
        FIXED = 0
        POS = INDEX(LINE, 'resseq')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'RESSEQ'//LINE(POS+6:)
          FIXED = FIXED+1
        END IF
        POS = INDEX(LINE, 'chain')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'CHAIN'//LINE(POS+5:)
          FIXED = FIXED+1
        END IF
        POS = INDEX(LINE, 'Chain')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'CHAIN'//LINE(POS+5:)
          FIXED = FIXED+1
        END IF
        POS = INDEX(LINE, 'resid')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'RESID'//LINE(POS+5:)
          FIXED = FIXED+1
        END IF
        POS = INDEX(LINE, 'resi ')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'RESID '//LINE(POS+5:)
          FIXED = FIXED+1
        END IF
        POS = INDEX(LINE, 'and')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'AND'//LINE(POS+3:)
          FIXED = FIXED+1
        END IF
        POS = INDEX(LINE, 'or')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'OR'//LINE(POS+2:)
          FIXED = FIXED+1
        END IF
        POS = INDEX(LINE, 'through')
        IF (POS.GT.0) THEN
          LINE  = LINE(1:POS-1)//'THROUGH'//LINE(POS+7:)
          FIXED = FIXED+1
        END IF
        IF (FIXED.EQ.0) GO TO 20
10    CONTINUE
     
20    RETURN 
      END
C-----------------------------------------------------------------------
C  This function extracts a real number from a string
C-----------------------------------------------------------------------
      REAL FUNCTION GTREAL(TLINE)

      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*(*) TLINE
      INTEGER      SPACE
      REAL VALUE

C-----Initialise
      SPACE = INDEX(TLINE, ' ')
      IF (SPACE.NE.0) THEN
        READ(UNIT=TLINE(1:SPACE), FMT=*, ERR=15) VALUE
      ELSE
        READ(UNIT=TLINE, FMT=*, ERR=15) VALUE
      END IF

C-----Skip fallback value
      GO TO 99
C-----In case of error use fall-back value
15    VALUE = 1.000

C-----Assign and return
99    GTREAL = VALUE
      RETURN
      END
C-----------------------------------------------------------------------
C  This function returns the number of data entries in a line
C-----------------------------------------------------------------------
      INTEGER FUNCTION DATCHK(LINE)
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER   I, DATCNT, NXTSPA, LENSTR
      LOGICAL   HAVDAT
      
C-----Intitialise
      DATCNT = 0
      HAVDAT = .FALSE.

C-----Analise line character by character
      DO 10, I=1, LENSTR(LINE)
        IF (LINE(I:I).NE.' ') THEN
C---------If not new value
	  IF (HAVDAT.EQV..TRUE.) GO TO 10
C---------New value
	  DATCNT = DATCNT+1
          HAVDAT = .TRUE.
	ELSE
	  HAVDAT = .FALSE.
        END IF
10    CONTINUE   
       
      DATCHK = DATCNT      
      RETURN 
      END      

C-----------------------------------------------------------------------
C  This subroutine replaces a minus character if it is used as a range 
C  separator
C-----------------------------------------------------------------------
      SUBROUTINE MINFIX(LINE)
      IMPLICIT  NONE
      CHARACTER LINE*(*), ONETEN*10, PROBE*2
      INTEGER   I, MINUS, LENSTR, LENGTH
      PARAMETER (ONETEN='1234567890')   
      
C-----Initialise
      MINUS = 0
      LENGTH = LENSTR(LINE)

C-----Check for occurences of 1-, 2-,...,9-
      DO 10, I=1, 10
        PROBE = ONETEN(I:I)//'-'
        MINUS = INDEX(LINE,PROBE)
        IF (MINUS.NE.0) THEN
          WRITE(6,*) 'Minus used as range separator'
          WRITE(6,*) LINE
C---------Replace the minus sign not the preceding number
          MINUS = MINUS+1 
          LINE = LINE(1:MINUS-1)//':'//LINE(MINUS+1:LENGTH)
        END IF 
10    CONTINUE

C-----Specific fix for cleared minus
!       MINUS = INDEX(LINE,' - ')
!       IF (MINUS.NE.0) THEN
!         LINE = LINE(1:MINUS-1)//':'//LINE(MINUS+3:LENGTH)
!       END IF
      
30    RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine clears a minus character if it directly precedes or 
C  follows a digit or letter.
C  CTYPE:     L = probe for letters; D = digits
C  POSITION: P = minus precedes the character; F = minus follows
C  Note: 'LINE' should have ample buffer spaces at the end
C-----------------------------------------------------------------------
      SUBROUTINE MINCLR(LINE, CTYPE, POSITION)
      IMPLICIT  NONE
      CHARACTER LINE*(*), DIGIT*10, PROBE*2, CTYPE*1, POSITION*1,
     +          LETTER*52
      INTEGER   I, MINUS, LENSTR, LENGTH, NFIX, LOOPC

C-----Initialise
      LETTER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
      DIGIT  = '1234567890'
      MINUS  = 0
      NFIX   = 0
      LENGTH = LENSTR(LINE)

C-----Set up the looping
      IF (CTYPE.EQ.'L') THEN
        LOOPC = 52
      ELSE IF (CTYPE.EQ.'D') THEN
        LOOPC = 10
      ELSE
        WRITE(6,*) 'Invalid character type selected!'
        GO TO 30
      END IF

C-----Check for occurences of probe
1     DO 10, I=1, LOOPC
        IF ((CTYPE.EQ.'D').AND.(POSITION.EQ.'F')) THEN
          PROBE = DIGIT(I:I)//'-'
        ELSE IF ((CTYPE.EQ.'D').AND.(POSITION.EQ.'P')) THEN
          PROBE = '-'//DIGIT(I:I)
        ELSE IF ((CTYPE.EQ.'L').AND.(POSITION.EQ.'F')) THEN
          PROBE = LETTER(I:I)//'-'
        ELSE IF ((CTYPE.EQ.'L').AND.(POSITION.EQ.'P')) THEN
          PROBE = '-'//LETTER(I:I)
        ELSE
          WRITE(6,*) 'Invalid minus position selected!'
        END IF
        MINUS = INDEX(LINE,PROBE)
        IF (MINUS.NE.0) THEN
C          PRINT*, 'Uncleared range minus found:'
C          PRINT*, LINE
C---------Replace the minus sign not the preceding number
          IF (POSITION.EQ.'F') THEN
            MINUS = MINUS+1
          END IF
          LINE  = LINE(1:MINUS-1)//' - '//LINE(MINUS+1:LENGTH)
C---------Update the length and iterate
          LENGTH = LENSTR(LINE)
          GO TO 1
        END IF
10    CONTINUE

30    RETURN
      END


C-----------------------------------------------------------------------
C  This subroutine removes a space character in line like '123: 456'
C-----------------------------------------------------------------------
      SUBROUTINE SPAFIX(LINE)
      IMPLICIT  NONE
      CHARACTER LINE*(*), ONENIN*9, PROBE*3
      INTEGER   I, GAP, LENSTR, LENGTH
      PARAMETER (ONENIN='123456789')   
      
C-----Initialise
      GAP = 0
      LENGTH = LENSTR(LINE)

C-----Check for occurences of ': 1', ': 2',..., ': 9'
      DO 10, I=1, 9
        PROBE = ': '//ONENIN(I:I)
        GAP   = INDEX(LINE,PROBE)
        IF (GAP.NE.0) THEN
          WRITE(6,*) 'Fixing space in line:'
          WRITE(6,*) LINE
C---------Remove the space character 
          LINE = LINE(1:GAP)//LINE(GAP+2:LENGTH)
        END IF 
10    CONTINUE
      
30    RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine loads all compound lists
C-----------------------------------------------------------------------
      SUBROUTINE GETCMP(HETERO, MAXHET, NHET, NUCLEIC, NNUC, AMINOA,
     +                  NAMINO, ALLRES, NALL, TLSRES, NTLSRES)
      IMPLICIT  NONE
C-----Variables:
C     MAXHET  The maximum number of hetero compounds in the list   
C     MAXIGN  The maximum number of residues on the ignore list
C     I       Just a loop iterator
C     J       Just a loop iterator
C     K       Just a loop iterator
C     MAXLIN  The maximum number of lines read form a file
C     NHET    The actual number of usable hetero compounds
C     NNUC    The actual number of linking nucleic acids
C     NAMINO  The actual number of linking amino acids
C     NIGN    The actual number of residues to be ignored
C     NTLSRES The actual number of residues that can be added to TLS 
C             groups
C     CHOP    The return of the CHOP function     
C     TMPRES  Temporary residue name
C     HETTYP  The type of residue
C     LINE    Temporary line of text 
      INTEGER   MAXHET, MAXIGN, NHET, NNUC, NAMINO, NIGN, NALL, 
     +          I, J, K, MAXLIN, NTLSRES
      PARAMETER (MAXIGN=200)
      PARAMETER (MAXLIN=90000)
      CHARACTER CHOP*250, TMPRES*3, RES, HETTYP*11, LINE*250, TLINE*250,
     +          GTTEXT*200
C-----Arrays:
C     HETERO  Usable hetero compounds 
C     NUCLEIC Polymeric nucleic acid residues
C     AMINOA  Amino acids
C     HETIGN  Hetero compounds to ignore    
      CHARACTER HETERO(MAXHET)*3, NUCLEIC(MAXHET)*3, HETIGN(MAXIGN)*3,
     +          AMINOA(MAXHET)*3, ALLRES(MAXHET)*3, TLSRES(MAXHET)*3

C-----Initialise values
      NIGN = 0

C-----First fill the ignore list and the TLS list
      DO 10, I=1, MAXLIN
C       Read until a data loop is found     
        READ(UNIT=12, FMT=998, END=20) LINE
        IF (LINE(1:5).EQ.'loop_') THEN
C         In a data loop; is it the right one?
          READ(UNIT=12, FMT=998, END=20) LINE
          IF (LINE(1:23).EQ.'_res_hetero_ignore_name') THEN
C           Loop over all entries
            DO 15, J=1, MAXLIN
              READ(UNIT=12, FMT=998, END=20) LINE
C             Stop if a new loop is reached
              IF(LINE(1:5).EQ.'loop_') GO TO 10
C             Skip comments
              IF(LINE(1:1).EQ.'#') GO TO 15
C             Read the data
              NIGN = NIGN+1
              READ(UNIT=LINE(1:3), FMT=*, ERR=989)HETIGN(NIGN)
15          CONTINUE
          ELSE IF (LINE(1:22).EQ.'_res_hetero_TLS_append') THEN
C           Loop over all entries
            DO 16, J=1, MAXLIN
              READ(UNIT=12, FMT=998, END=20) LINE
C             Stop if a new loop is reached
              IF(LINE(1:5).EQ.'loop_') GO TO 10
C             Skip comments
              IF(LINE(1:1).EQ.'#') GO TO 16
C             Read the data
              NTLSRES = NTLSRES+1
              READ(UNIT=LINE, FMT=997, ERR=989)TLSRES(NTLSRES)
16          CONTINUE               
          END IF
        END IF
10    CONTINUE


C-----Report the stats
20    WRITE(6,*)'Hetero compounds:'

C-----Read the hetero compounds
      DO 30, I=1, MAXLIN
C       Read lines to find the proper data loop
        READ(UNIT=11, FMT=998, END=990) LINE
        IF (LINE(1:5).EQ.'loop_') THEN
C         In a data loop; is it the right one?
          READ(UNIT=11, FMT=998, END=990) LINE
C         Remove leading spaces
          LINE = CHOP(LINE,6) 
          IF (LINE(1:13).EQ.'_chem_comp.id') THEN
C           Loop over all entries
            DO 33, J=1, MAXLIN
              READ(UNIT=11, FMT=998, END=990) LINE
C             Skip lines with lables or comment lines
              IF(LINE(1:1).EQ.'_') GO TO 33
              IF(LINE(1:1).EQ.'#') GO TO 33
C             Stop at the start of a new block
              IF(LINE(1:5).EQ.'data_') GO TO 40
C             Read the data from a copy of LINE to bypass side effects
              TLINE  = LINE
              TMPRES = GTTEXT(TLINE,1)
              TLINE  = LINE
              HETTYP = GTTEXT(TLINE,4)
C             Add the data if it can be used 
              NALL = NALL+1
              ALLRES(NALL) = TMPRES
              IF ((HETTYP.EQ.'non-polymer').OR. 
     +            (HETTYP.EQ.'NON-POLYMER'))THEN
C               Exclude compounds in the HETIGN list 
                DO 36, K=1, NIGN
                  IF (TMPRES.EQ.HETIGN(K)) GO TO 33                
36              CONTINUE
                NHET = NHET+1
                HETERO(NHET) = TMPRES 
C                WRITE(6,*) HETERO(NHET)
              ELSE IF (HETTYP(2:4).EQ.'NA ') THEN
                NNUC = NNUC+1
                NUCLEIC(NNUC) = TMPRES 
              ELSE IF ((HETTYP(1:7).EQ.'peptide').OR.
     +                 (HETTYP(1:9).EQ.'M-peptide')) THEN
                NAMINO = NAMINO+1
                AMINOA(NAMINO) = TMPRES 
              END IF 
33          CONTINUE            
          END IF
        END IF   
30    CONTINUE     

C-----Report
40    WRITE(UNIT=6,FMT=996) 'Loaded residues        : ', NALL
      WRITE(UNIT=6,FMT=996) 'Residues to be ignored : ', NIGN
      WRITE(UNIT=6,FMT=996) 'Usable residues loaded : ', NHET
      WRITE(UNIT=6,FMT=996) 'Nucleic acid residues  : ', NNUC
      WRITE(UNIT=6,FMT=996) 'Amino acid residues    : ', NAMINO
      WRITE(UNIT=6,FMT=996) 'Residues to add to TLS : ', NTLSRES
C      WRITE(UNIT=6,FMT='(2200A3)') (ALLRES(i),i=1,NALL)

C-----Skip error messages
      GO TO 999

C-----Error messages
989   WRITE(6,*) 'Cannot read residue'
      GO TO 999
990   WRITE(6,*) 'Unexpected end of file'
      GO TO 999
      

C-----Formats
995   FORMAT(A3,49X,A11)
996   FORMAT(A,I5)
997   FORMAT(A3)
998   FORMAT(A250)

999   RETURN
      END
      
C-----------------------------------------------------------------------
C  This function extracts a string from column CLMNR from a text 
C  line. WARNING: this modifies LINE so send in a copy!
C-----------------------------------------------------------------------
      CHARACTER*(*) FUNCTION GTTEXT(LINE, CLMNR)
      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*(*) LINE
      CHARACTER QUOTE*1, CHOP*1024
      INTEGER   CLMNR, J, SPACE


C-----Initialise
      SPACE = 1
      GTTEXT = '?' 
C      PRINT*, LINE
      
C-----Extract right column
      DO 10, J=1, CLMNR
        LINE = LINE(SPACE:LEN(LINE))
	LINE = CHOP(LINE, LEN(LINE))
C        PRINT*, TLINE

C-------Special case when the value is a quote
        READ(UNIT=LINE(1:1), FMT='(A1)') QUOTE
        IF (QUOTE.EQ."'") THEN
C         Read ahead to the next quote          
          SPACE = INDEX(LINE, "' ") + 1
C          PRINT*, LINE(1:SPACE)
          IF (SPACE.EQ.1) GO TO 99
          READ(UNIT=LINE(1:SPACE-1), FMT='(A)', END=15, ERR=15) GTTEXT
          GO TO 10
        END IF 
        
        READ(UNIT=LINE(1:1), FMT='(A1)') QUOTE
        IF (QUOTE.EQ.'"') THEN
C         Read ahead to the next quote          
          SPACE = INDEX(LINE, '" ') + 1
          IF (SPACE.EQ.1) GO TO 99
          READ(UNIT=LINE(1:SPACE-1), FMT='(A)', END=15, ERR=15) GTTEXT
          GO TO 10
        END IF        

C-------Read normally
        SPACE = INDEX(LINE, ' ')
        IF (SPACE.EQ.0) GO TO 99
        READ(UNIT=LINE(1:SPACE), FMT=*, END=15, ERR=15) GTTEXT
	GO TO 10 
C-----In case of error use fall-back value
15      GTTEXT = '?'
10    CONTINUE
        

C-----Assign and return     
C      WRITE(UNIT=6, FMT='(F14.4)') GTREAL 
99    RETURN
      END
