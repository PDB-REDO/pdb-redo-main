      PROGRAM CHIRON
C=======================================================================
C  Version 4.00 2014-10-31
C  Fixes administrative chirality problems by swapping atoms in certain 
C  compounds or by renaming the compounds.
C
C  Usage: CHIRON (-v) DF_IN LOG_IN PDB_IN PDB_OUT
C  -v        Verbose output
C  DF_IN      is the pdb_redo.dat datafile (unit 10)
C  LOG_IN     is a Refmac log file (unit 7)
C  PDB_IN     is the original PDB file (unit 8)
C  PDB_OUT    is the output PDB file (unit 9)
C
C  Written by Robbie P. Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  CHIRON. Please, refer to (one of) these references:
C  - Robbie P. Joosten, Gert Vriend: "PDB improvement starts with data 
C    deposition" Science, 317, p. 195-196 (2007)
C  - Robbie P. Joosten, Thomas Womack, Gert Vriend and Gerard Bricogne:
C    "Re-refinement fromdeposited X-ray data can deliver improved models
C    for most PDB entries"  Acta Cryst. D65, p. 176-185 (2009)
C  - Robbie P. Joosten, Jean Salzemann, Vincent Bloch, Heinz Stockinger,
C    Ann-Charlott Berglund, Christophe Blanchet, Erik Bongcam-Rudloff, 
C    Christophe Combet, Ana L. Da Costa, Gilbert Deleage, Matteo 
C    Diarena, Roberto Fabbretti, Géraldine Fettahi, Volker Flegel, 
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
C  Change log
C  Version 4.00
C  - Some chiral volume errors are now resolved by ranaming a compound 
C    to its epimer.
C 
C  Version 3.01
C  - LINK(R) records without explicit atoms do no longer trigger an
C    error message.
C
C  Version 3.00
C  - Chiral atoms in a chiral branch are now considered mutally 
C    exclusive. If an error occurs in both atoms, fixing the atom top 
C    level atom is preferred. 
C  - The second atom can be fixed after another cycle of chirality 
C    checking and chiron.
C
C  Version 2.02
C  - Added a subroutine to deal with a few specific residues.
C  - Completely flat groups are now handled better.
C
C  Version 2.01
C  - Atoms are now also swapt in LINK records.
C  - Updated the references.
C  - Increased MAXFIX.
C
C  Version 2.00
C  - Better handling of chirality errors in loops with two chiral atoms.
C  - This also enables fixing more than one atom per residue.
C
C  Version 1.03
C  - Added a basic help function.
C
C  Version 1.2
C  - Now ignores chirality deviations that do not involve a change of 
C    sign.
C  - Up to 8 (previously 6) equivalent atoms are supported.
C  - Added support for alternates.
C
C  Version 1.1
C  - The program now takes a separate file with fixable atoms.
C  - Before further evaluating a chiral center, the deviation is tested.
C    If the outlier is too small, it will be ignored.
C
C  Version 1.0
C  - The size of the chirality deviation is now taken into account.
C  - ANISOU records are now also fixed.
C
C  Version 0.0
C  - First attempt.
C
C=======================================================================
      IMPLICIT NONE
C-----Declare the basic variables and parameters
      INTEGER   I, J, K, STATUS, MAXLIN, MAXFIX, MAXDB, ARGS, EXTRA
      CHARACTER VERS*4
      PARAMETER (VERS='4.00')
C-----MAXLIN is the maximum number of lines in the logfile
      PARAMETER (MAXLIN=200000)
C-----MAXFIX is the maximum number chirality errors
      PARAMETER (MAXFIX=100)
C-----MAXDB is the maximal number of entries in the data file
      PARAMETER (MAXDB=2000)

C-----Single values
      CHARACTER LOGIN*255, PDBIN*255, PDBOUT*255, DFIN*255, C2JUNK*2, 
     +          LINE*100,  L80*80
      LOGICAL   VERBOS 
      INTEGER   FIXABL, POS, NNAMFIX
C     NNAMFIX Number of residues that must be renamed          

C-----Arrays
      CHARACTER FIXID(MAXFIX)*7, FIXATM(MAXFIX)*4, RESATM(MAXDB)*7,
     +          SWAPS(MAXDB)*64, NONFIX(MAXDB)*3,  COUPLES(MAXDB)*11,
     +          EXCLUDE(MAXDB)*15, EPIMER(MAXDB)*10, NAMFIX(MAXFIX)*9
C     FIXID  is a combination of the chain ID, residue number, insertion
C            code and alternate ID if a residue
C     FIXATM is the atom ID of a chiral center
C     RESATM is the residue ID followed by the atom name of the chiral
C            center
C     EXLUDE contains a residue name and the nems of to atoms that are 
C            mutually exclusive for fixing
C     EPIMER contains a list of epimers 
C            format: 1-4 epimer center, 5-7 resid1, 8-10 resid2
C     NAMFIX list of residues that must be renamed
C            format: 1 chain ID, 2-5 residue number, 6 insertion code,
C                    7-9 new residue name 
      INTEGER   NTRY(5), FIXSTA(MAXFIX)
C     NTRY(1) is the number of chiral atoms that are fixable 
C     NTRY(2) is the number of residues that are non-fixable 
C     NTRY(3) is the number of 'coupled' chiral atom pairs
C     NTRY(4) is the number of mutually exclusive atom pairs
C     NTRY(5) is the numebr of epimer pairs
C     FIXSTA  is the status of a chiral atom: 1 = fix, 0 = not fix, 
C             -1 = unchecked                  
                                                
C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   CHIRON version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Chiron fixes chirality problems by swapping chemical'//
     +          'ly equivalent atoms and renaming residues' 
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'chiron (flags) DF_IN LOG_IN PDB_IN PDB_OUT'
      WRITE(6,*)' '
      WRITE(6,*)'DF_IN   the pdb_redo.dat file.'
      WRITE(6,*)'LOG_IN  a log file from Refmac.'
      WRITE(6,*)'PDB_IN  the PDB file that needs chirality fixes.'
      WRITE(6,*)'PDB_OUT the fixed PDB file.'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode. Gives more info about which atoms a'//
     +          're considered.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing CHIRON:'
      WRITE(6,*)'Joosten et al. ''PDB_REDO: constructive validation, '//
     +          'more than just looking for errors'' Acta Cryst., D68'//
     +          ', 484-496 (2012)'
        GO TO 999
      END IF

C--------------------------- Main program -----------------------------C

C-----Initialise
      EXTRA    = 0
      FIXABL   = 0
      NTRY(1)  = 0
      NTRY(2)  = 0
      NTRY(3)  = 0
      NTRY(4)  = 0
      NTRY(5)  = 0
      VERBOS   = .FALSE.

C-----Check for (debug) flags
      DO 10, I=1, 2
        CALL GETARG(1+EXTRA, C2JUNK)
C-------Are there any flags?
        IF (C2JUNK(1:1).NE.'-') GO TO 20
C       Is it a valid flag 
        IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
          VERBOS = .TRUE.
          EXTRA  = EXTRA+1
        ELSE
	  EXTRA = EXTRA+1
        END IF
10    CONTINUE  

C-----Get input files
20    STATUS=0
C     Data file
      CALL GETARG(1+EXTRA, DFIN)
      OPEN(UNIT=10, FILE=DFIN, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) 'Cannot open ', DFIN 
          GO TO 999
        END IF
      REWIND (10)
C     First log file
      CALL GETARG(2+EXTRA, LOGIN)
      OPEN(UNIT=7, FILE=LOGIN, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) 'Cannot open ', LOGIN 
          GO TO 999
        END IF
      REWIND (7)
C     Second log file    
      CALL GETARG(3+EXTRA, PDBIN)
      OPEN(UNIT=8, FILE=PDBIN, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) 'Cannot open ', PDBIN
          GO TO 999
        END IF
      REWIND (8)

C-----Get output file
      STATUS=0                                     
      CALL GETARG((4+EXTRA), PDBOUT)
      OPEN(UNIT=9, FILE=PDBOUT, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) PDBOUT, 'cannot be created!'
        GO TO 999
      END IF
C-------------------------- Parse data file----------------------------C

C-----Loop over data file to find (im)possible fixes
      DO 25, I=1, MAXLIN
        READ(UNIT=10, FMT=901, END=29) LINE
C       Skip comment lines
        IF (LINE(1:1).EQ.'#') GO TO 25
        IF (LINE(1:5).EQ.'loop_') THEN
          READ(UNIT=10, FMT=901, END=29) LINE 
          IF (LINE(1:9).EQ.'_res_atom') THEN
            CALL LDSWAPS(RESATM, SWAPS, NTRY, MAXDB)
          ELSE IF (LINE(1:18).EQ.'_res_nonfix_chiron') THEN
            CALL LDNONFIX(NONFIX, NTRY, MAXDB)
          ELSE IF (LINE(1:12).EQ.'_res_couples') THEN
            CALL LDCOUPLES(COUPLES, NTRY, MAXDB)
          ELSE IF (LINE(1:12).EQ.'_atom_exclus') THEN
            CALL LDEXCLUDE(EXCLUDE, NTRY, MAXDB)
          ELSE IF (LINE(1:19).EQ.'_carb_epimer_center') THEN
            CALL LDEPIMER(EPIMER, NTRY, MAXDB)
          END IF 
        END IF
25    CONTINUE

C-------------------------- Parse log file ----------------------------C

C-----Loop over log file to find the chirality error block
29    DO 30, I=1, MAXLIN
        READ(UNIT=7, FMT=900, END=30) LINE
        IF (LINE(1:24).EQ.'Chiral volume deviations') THEN
          CALL GETERR(FIXABL,FIXID,FIXATM,VERBOS,MAXFIX,RESATM,NONFIX,
     +      COUPLES,NTRY,FIXSTA, EXCLUDE, EPIMER, NAMFIX, NNAMFIX)
C         Bail out if there was a problem
          IF (FIXABL.LT.0) GO TO 996
          GO TO 100
        END IF
30    CONTINUE

C----------------------- Apply the fixes ------------------------------C

C-----Loop over the PDB file and fix the atoms
100   DO 101, I=1, MAXLIN
        READ(UNIT=8, FMT=901, END=999) L80
C       Only look at coordinates..
        IF ((L80(1:6).EQ.'ATOM  ').OR.(L80(1:6).EQ.'HETATM').OR.
     +      (L80(1:6).EQ.'ANISOU')) THEN
C         Loop over all residues to fix by atom swapping
          DO 110, J=1, FIXABL
C           This is a fixable residue, and it really must be fixed
            IF ((L80(22:27).EQ.FIXID(J)(1:6)).AND.
     +          (L80(17:17).EQ.FIXID(J)(7:7)).AND.
     +          (FIXSTA(J).EQ.1)) THEN
C             Find the right atoms to swap              
              DO 115, K=1, NTRY(1)
C               Is it the right atom and residue type 
                IF ((FIXATM(J).EQ.RESATM(K)(4:7)).AND.
     +              (L80(18:20).EQ.RESATM(K)(1:3))) THEN
C                 Is the atom in the swaplist
                  POS = INDEX(SWAPS(K),L80(13:16))
C                  PRINT*, L80, SWAPS(K), POS
                  IF (POS.NE.0) THEN
C                   Swap the atoms
                    IF(MOD((POS-1),8).EQ.0) THEN
                      L80(13:16) = SWAPS(K)(POS+4:POS+7)
                    ELSE IF (MOD((POS-1),8).EQ.4)THEN
                      L80(13:16) = SWAPS(K)(POS-4:POS-1)
                    ELSE
C                      PRINT*, '1ST loop'
                      GO TO 997
                    END IF
                    GO TO 199
                  END IF
                END IF
115           CONTINUE     
            END IF
110       CONTINUE 
C         Loop over all residues to be fixed by renaming
          DO 120, J=1, NNAMFIX
C           Should the residue be renamed
            IF (L80(22:27).EQ.NAMFIX(J)(1:6)) THEN
              L80(18:20) = NAMFIX(J)(7:9)
              GO TO 199
            END IF         
120       CONTINUE            
        END IF

C       Check the LINK records
        IF (L80(1:4).EQ.'LINK') THEN
C         Loop over all residues to fix by swapping atoms
          DO 130, J=1, FIXABL
C           This is a fixable residue, on the left of the LINK
            IF ((L80(22:27).EQ.FIXID(J)(1:6)).AND.(FIXSTA(J).EQ.1)) THEN
C             Find the right atoms to swap              
              DO 131, K=1, NTRY(1)
C               Find the right atom swap by residue ID
                IF (L80(18:20).EQ.RESATM(K)(1:3)) THEN
C                 Move on if there is no atom
                  IF (L80(13:16).EQ.'    ') GO TO 132
C                 Is the atom in the swaplist
                  POS = INDEX(SWAPS(K),L80(13:16))
                  IF (POS.NE.0) THEN
C                   Swap the atoms
                    IF(MOD((POS-1),8).EQ.0) THEN
                      L80(13:16) = SWAPS(K)(POS+4:POS+7)
                    ELSE IF (MOD((POS-1),8).EQ.4)THEN
                      L80(13:16) = SWAPS(K)(POS-4:POS-1)
                    ELSE
C                      PRINT*, '2nd loop'
                      GO TO 997
                    END IF
C                   Continue on the right side of the LINK
                    GO TO 132
                  END IF
                END IF
131           CONTINUE
            END IF 
C           This is a fixable residue, on the right side of the LINK    
132         IF ((L80(52:57).EQ.FIXID(J)(1:6)).AND.(FIXSTA(J).EQ.1)) THEN
C             Find the right atoms to swap  
C              PRINT*, L80            
              DO 133, K=1, NTRY(1)
C               Find the right atom swap by residue ID
                IF (L80(48:50).EQ.RESATM(K)(1:3)) THEN
C                 Move on if there is no atom
                  IF (L80(43:46).EQ.'    ') GO TO 135
C                 Is the atom in the swaplist
                  POS = INDEX(SWAPS(K),L80(43:46))
                  IF (POS.NE.0) THEN
C                   Swap the atoms
                    IF(MOD((POS-1),8).EQ.0) THEN
                      L80(43:46) = SWAPS(K)(POS+4:POS+7)
                    ELSE IF (MOD((POS-1),8).EQ.4)THEN
                      L80(43:46) = SWAPS(K)(POS-4:POS-1)
                    ELSE
                      PRINT*, '3rd loop'
                      GO TO 997
                    END IF
C                   Go on to renaming residues      
                    GO TO 199
                  END IF
                END IF
133           CONTINUE 
            END IF
130       CONTINUE

C         Check for residues to rename (on both sides at once)
135       DO 136, J=1, NNAMFIX
            IF (L80(22:27).EQ.NAMFIX(J)(1:6)) THEN
              L80(18:20) = NAMFIX(J)(7:9) 
            ELSE IF (L80(52:57).EQ.NAMFIX(J)(1:6)) THEN
              L80(48:50) = NAMFIX(J)(7:9)
            END IF
136       CONTINUE
             
        END IF
199     WRITE(UNIT=9,FMT=901) L80
101   CONTINUE

C-----Skip the error messages
      GO TO 999

C-----Formats
900   FORMAT(A100)
901   FORMAT(A80)

      

C-----Error messages
996   WRITE(6,*) 'Problem reading the log file.'
      CLOSE(9,STATUS='DELETE')
      GO TO 999  
997   WRITE(6,*) 'Problem swapping atoms!'
      CLOSE(9,STATUS='DELETE')
      GO TO 999  
998   WRITE(6,*) 'Cannot read the number of atoms'        

C-----End of the line
999   CLOSE (7)
      CLOSE (8)
      CLOSE (9)
      CLOSE (10)
      END

C-------------------------- Subroutines and functions -----------------C
C-----------------------------------------------------------------------
C  This subroutine loads all the fixable chiral centres and the swaps 
C  required to make the fixes.
C-----------------------------------------------------------------------
      SUBROUTINE LDSWAPS(RESATM, SWAPS, NTRY, MAXDB)
      IMPLICIT  NONE
      INTEGER   MAXLIN, MAXDB, NTRY(*), I
      CHARACTER RESATM(MAXDB)*7, SWAPS(MAXDB)*64, LINE*100
      PARAMETER (MAXLIN=10000)

C-----Loop over the data block until the first # is found
      DO 10, I=1, MAXLIN
        READ(UNIT=10, FMT=998, END=990) LINE
C-------Skip lines with lables
        IF (LINE(1:1).EQ.'_') GO TO 10
C-------Stop at the end of the block
        IF (LINE(1:1).EQ.'#') GO TO 20
C-------Read the data
        NTRY(1) = NTRY(1)+1
        READ(UNIT=LINE, FMT=997, ERR=989)RESATM(NTRY(1)),SWAPS(NTRY(1))
10    CONTINUE

C-----Report the stats
20    WRITE(UNIT=6,FMT=996) 'Fixable chiral centres read  : ', NTRY(1) 
      GO TO 999

C-----Error messages
989   WRITE(6,*) 'Cannot read atom swap'
      GO TO 999
990   WRITE(6,*) 'Unexpected end of file'
      GO TO 999
      

C-----Formats
996   FORMAT(A31,I4)
997   FORMAT(5X,A7,4X,A64)
998   FORMAT(A100)

999   RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine loads all the non-fixable residues.
C-----------------------------------------------------------------------
      SUBROUTINE LDNONFIX(NONFIX, NTRY, MAXDB)
      IMPLICIT  NONE
      INTEGER   MAXLIN, MAXDB, NTRY(*), I
      CHARACTER NONFIX(MAXDB)*3, LINE*80
      PARAMETER (MAXLIN=1000)

C-----Loop over the data block until the first # is found
      DO 10, I=1, MAXLIN
        READ(UNIT=10, FMT=998, END=990) LINE
C-------Skip lines with lables
        IF(LINE(1:1).EQ.'_') GO TO 10
C-------Stop at the end of the block
        IF(LINE(1:1).EQ.'#') GO TO 20
C-------Read the data
        NTRY(2) = NTRY(2)+1
        READ(UNIT=LINE, FMT=997, ERR=989)NONFIX(NTRY(2))
10    CONTINUE

C-----Report the stats
20    WRITE(UNIT=6,FMT=996) 'Non-fixable residues read    : ', NTRY(2)
      GO TO 999

C-----Error messages
989   WRITE(6,*) 'Cannot read residue'
      GO TO 999
990   WRITE(6,*) 'Unexpected end of file'
      GO TO 999
      

C-----Formats
996   FORMAT(A31,I4)
997   FORMAT(A3)
998   FORMAT(A80)

999   RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine loads all the coupled chiral atoms.
C-----------------------------------------------------------------------
      SUBROUTINE LDCOUPLES(COUPLES, NTRY, MAXDB)
      IMPLICIT  NONE
      INTEGER   MAXLIN, MAXDB, NTRY(*), I
      CHARACTER COUPLES(MAXDB)*11, LINE*80
      PARAMETER (MAXLIN=1000)

C-----Loop over the data block until the first # is found
      DO 10, I=1, MAXLIN
        READ(UNIT=10, FMT=998, END=990) LINE
C-------Skip lines with lables
        IF(LINE(1:1).EQ.'_') GO TO 10
C-------Stop at the end of the block
        IF(LINE(1:1).EQ.'#') GO TO 20
C-------Read the data
        NTRY(3) = NTRY(3)+1
        READ(UNIT=LINE, FMT=997, ERR=989)COUPLES(NTRY(3))
10    CONTINUE

C-----Report the stats
20    WRITE(UNIT=6,FMT=996) 'Coupled chiral atoms read    : ', NTRY(3)
      GO TO 999

C-----Error messages
989   WRITE(6,*) 'Cannot read coupled atoms'
      GO TO 999
990   WRITE(6,*) 'Unexpected end of file'
      GO TO 999


C-----Formats
996   FORMAT(A31,I4)
997   FORMAT(5X,A11)
998   FORMAT(A80)

999   RETURN
      END
C-----------------------------------------------------------------------
C  This subroutine loads all the coupled chiral atoms.
C-----------------------------------------------------------------------
      SUBROUTINE LDEXCLUDE(EXCLUDE, NTRY, MAXDB)
      IMPLICIT  NONE
      INTEGER   MAXLIN, MAXDB, NTRY(*), I
      CHARACTER EXCLUDE(MAXDB)*15, LINE*80, TRES*3, TROOT*4, TBRANCH*8
      PARAMETER (MAXLIN=1000)

C-----Loop over the data block until the first # is found
      DO 10, I=1, MAXLIN
        READ(UNIT=10, FMT=998, END=990) LINE
C-------Skip lines with lables
        IF(LINE(1:1).EQ.'_') GO TO 10
C-------Stop at the end of the block
        IF(LINE(1:1).EQ.'#') GO TO 20
C-------Read the data
        NTRY(4) = NTRY(4)+1
        READ(UNIT=LINE, FMT=997, ERR=989) TRES, TROOT, TBRANCH
        EXCLUDE(NTRY(4)) = TRES//TROOT//TBRANCH
C        PRINT*, EXCLUDE(NTRY(4))
10    CONTINUE

C-----Report the stats
20    WRITE(UNIT=6,FMT=996) 'Mutually exclusive atoms read: ', NTRY(4)
      GO TO 999

C-----Error messages
989   WRITE(6,*) 'Cannot read mutually exclusive atoms'
      GO TO 999
990   WRITE(6,*) 'Unexpected end of file'
      GO TO 999


C-----Formats
996   FORMAT(A31,I4)
997   FORMAT(A3,2X,A4,3X,A8)
998   FORMAT(A80)

999   RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine loads all the epimers.
C-----------------------------------------------------------------------
      SUBROUTINE LDEPIMER(EPIMER, NTRY, MAXDB)
      IMPLICIT  NONE
      INTEGER   MAXLIN, MAXDB, NTRY(*), I
      CHARACTER EPIMER(MAXDB)*10, LINE*80
      PARAMETER (MAXLIN=1000)

C-----Loop over the data block until the first # is found
      DO 10, I=1, MAXLIN
        READ(UNIT=10, FMT=998, END=20) LINE
C-------Skip lines with labels, comments or empty lines
        IF((LINE(1:1).EQ.'_').OR.(LINE(1:1).EQ.'#').OR.
     +     (LINE(1:5).EQ.'     ')) GO TO 10
C-------End on a new loop
        IF (LINE(1:5).EQ.'loop_') THEN
          GO TO 20
        END IF  
C-------Read the data
        NTRY(5) = NTRY(5)+1
        READ(UNIT=LINE, FMT=997, ERR=989) EPIMER(NTRY(5))(1:4),
     +    EPIMER(NTRY(5))(5:10)
C        PRINT*, EPIMER(NTRY(5))
10    CONTINUE

C-----Report the stats
20    WRITE(UNIT=6,FMT=996) 'Epimer sets read             : ', NTRY(5)
      BACKSPACE(10)
      GO TO 999

C-----Error messages
989   WRITE(6,*) 'Cannot read epimer set'
      GO TO 999

C-----Formats
996   FORMAT(A31,I4)
997   FORMAT(1X,A4,2X,A6)
998   FORMAT(A80)

999   RETURN
      END
C-----------------------------------------------------------------------
C  This subroutine hacks the chirality data line for certain compounds
C-----------------------------------------------------------------------
      SUBROUTINE CMPFIX(LINE)
      IMPLICIT  NONE
      CHARACTER LINE*100, FIXRES*7
      PARAMETER (FIXRES='CAS BEF')

C-----Example:
CB    318 CAS AS    mod.=   0.00 id.=   6.80 dev=  6.801 sig.=  0.200  
CA    201 BEF BE    mod.=   3.71 id.=   0.00 dev= -3.711 sig.=  0.200
C1234567890123456
  
C-----Do we have a residue that needs hacking
      IF (INDEX(FIXRES,LINE(10:12)).NE.0) THEN
        IF ((LINE(10:12).EQ.'CAS').AND.(LINE(13:16).EQ.' AS ')) THEN 
          LINE(13:16) = 'AS  '
        ELSE IF((LINE(10:12).EQ.'BEF').AND.(LINE(13:16).EQ.' BE ')) THEN 
          LINE(13:16) = 'BE  '    
        END IF
      END IF 

99    RETURN
      END
C-----------------------------------------------------------------------
C  This subroutine reads a table of chiral outliers and checks whether 
C  they are known and fixable 
C-----------------------------------------------------------------------
      SUBROUTINE GETERR(FIXABL,FIXID,FIXATM,VERBOS,MAXFIX,RESATM,
     +                  NONFIX,COUPLES,NTRY,FIXSTA, EXCLUDE, EPIMER,
     +                  NAMFIX, NNAMFIX)
      IMPLICIT  NONE
      INTEGER   FIXABL, CHIATM, UNKNOW, UNFIX, MAXLIN, MAXFIX, EXFIX,
     +          I, J, K
      INTEGER   NONRES, TOFIX, NTRY(*), FIXSTA(MAXFIX), CNT, ATSWP,
     +          NNAMFIX
C     ATSWP   Chirality problems solved by atom swapping
C     NNAMFIX Numbere of chiral centers fixed by renaming residues
      REAL      DCHIR, CHIORI, CHITAR
      CHARACTER LINE*100, TEMPID*7, FIXID(MAXFIX)*7,  FIXATM(MAXFIX)*4,
     +          RESATM(NTRY(1))*7 , NONFIX(NTRY(2))*3,
     +          COUPLES(NTRY(3))*11,EXCLUDE(NTRY(4))*15,TRES(MAXFIX)*3,
     +          EPIMER(NTRY(5))*10, NAMFIX(MAXFIX)*9 
      LOGICAL   VERBOS
      PARAMETER (MAXLIN=100)

C-----Initialise
      CNT     = 0
      UNKNOW  = 0
      UNFIX   = 0
      TOFIX   = 0
      EXFIX   = 0
      ATSWP   = 0
      NNAMFIX = 0
      DCHIR   = 0.0

      DO 10, I=1, MAXLIN
        READ(UNIT=7, FMT=998, END=990) LINE
C-------Stop at the end of the table
        IF(INDEX(LINE,'****').NE.0) GO TO 40
        IF(INDEX(LINE,'Limits').NE.0) GO TO 40
C-------If there is an outlier check if it can be fixed
        IF(LINE(20:25).EQ.'mod.=') THEN
          CNT = CNT + 1          
C         Should it be fixed?
C         Is the chirality deviation big enough
          READ(UNIT=LINE(49:55),FMT=996,ERR=989) DCHIR
          IF (ABS(DCHIR).LT.3.5) THEN
            IF (VERBOS.EQV..TRUE.) THEN
              PRINT*,'Need not fix residue: ', LINE
            END IF
            UNFIX = UNFIX+1
            GO TO 10
          END IF
C         Is it really a chirality inversion
          READ(UNIT=LINE(26:31),FMT=995,ERR=989) CHIORI
          READ(UNIT=LINE(38:43),FMT=995,ERR=989) CHITAR
C         Do the chiralities have the same sign?
          IF (((CHIORI.GE.0).AND.(CHITAR.GT.0)).OR.
     +        ((CHIORI.LE.0).AND.(CHITAR.LT.0))) THEN
            IF (VERBOS.EQV..TRUE.) THEN
              PRINT*,'Need not fix residue: ', LINE
            END IF
            UNFIX = UNFIX+1
            GO TO 10
          END IF
C         Fix some compounds          
          CALL CMPFIX(LINE)
C         Is the atom in the list
          DO 20, J=1, NTRY(1)
C           Yes, so put it on the list
            IF (LINE(10:16).EQ.RESATM(J)) THEN
              TEMPID = LINE(1:1)//LINE(5:9)//LINE(18:18)  
              IF (VERBOS.EQV..TRUE.) THEN
                PRINT*,'Found fixable residue: ', LINE
              END IF
              FIXABL = FIXABL+1
C             Exit if FIXABL > MAXFIX (this is an array overflow) 
              IF (FIXABL.GT.MAXFIX) GO TO 988
              TRES(FIXABL)  = LINE(10:12)
              FIXID(FIXABL) = TEMPID
              FIXATM(FIXABL)= LINE(13:16)
              FIXSTA(FIXABL)= -1
              GO TO 10
            END IF            
20        CONTINUE  
C         No, can the problem be fixed by renaming the residue
          DO 25, J=1, NTRY(5)
C           Check whether the atom and residue are on epimer the list  
            IF ((LINE(13:16).EQ.EPIMER(J)(1:4)).AND.
     +          (LINE(10:12).EQ.EPIMER(J)(5:7))) THEN
C             Add the residue to the list and jump formward
              TOFIX   = TOFIX+1 
              NNAMFIX = NNAMFIX+1 
              NAMFIX(NNAMFIX) = LINE(1:1)//LINE(5:9)//EPIMER(J)(8:10)
              GO TO 10
            ELSE IF ((LINE(13:16).EQ.EPIMER(J)(1:4)).AND.
     +               (LINE(10:12).EQ.EPIMER(J)(8:10))) THEN
C             Add the residue to the list and jump formward
              TOFIX   = TOFIX+1 
              NNAMFIX = NNAMFIX+1 
              NAMFIX(NNAMFIX) = LINE(1:1)//LINE(5:9)//EPIMER(J)(5:7)
              GO TO 10
            END IF
25        CONTINUE
C         No, so is it on the unfixable list?
          DO 30, J=1, NTRY(2)
C           Yes, so count it
            IF (LINE(10:12).EQ.NONFIX(J)) THEN
              IF (VERBOS.EQV..TRUE.) THEN
                PRINT*, 'Non-fixable residue: ', LINE
              END IF
              UNFIX = UNFIX+1
              GO TO 10
            END IF            
30        CONTINUE 
C         This residue is unknown
          IF (VERBOS.EQV..TRUE.) THEN
            PRINT*, 'Unknown residue: ', LINE
          END IF
          UNKNOW = UNKNOW+1
        END IF
10    CONTINUE

C-----First step of summary
40    WRITE(UNIT=6,FMT=997) 'Chirality errors tested : ', CNT

C-----Filter out the coupled atoms

C     Loop over all fixable atoms
      DO 50, I=1, FIXABL
C       Skip the atoms already matched
        IF (FIXSTA(I).NE. -1) GO TO 50
C       Is the atom in the couples list
        DO 55, J=1, NTRY(3)
          IF(TRES(I).EQ.COUPLES(J)(1:3)) THEN 
            IF (FIXATM(I).EQ.COUPLES(J)(4:7)) THEN
C             Loop over the remaining atoms to find its partner
              DO 60, K=I+1, FIXABL
                IF((FIXID(I).EQ.FIXID(K)).AND.
     +             (FIXATM(K).EQ.COUPLES(J)(8:11))) THEN
                  FIXSTA(I) = 1 
                  FIXSTA(K) = 0 
                  TOFIX     = TOFIX + 2
                  GO TO 50
                END IF
60            CONTINUE
              FIXSTA(I) = 0
              UNFIX     = UNFIX + 1
              GO TO 50
            ELSE IF (FIXATM(I).EQ.COUPLES(J)(8:11))THEN
C             Loop over the remaining atoms to find its partner
              DO 65, K=I+1, FIXABL
                IF((FIXID(I).EQ.FIXID(K)).AND.
     +             (FIXATM(K).EQ.COUPLES(J)(4:7))) THEN
                  FIXSTA(I) = 1
                  FIXSTA(K) = 0
                  TOFIX     = TOFIX + 2
                  GO TO 50
                END IF
65            CONTINUE
              FIXSTA(I) = 0       
              UNFIX     = UNFIX + 1
              GO TO 50
            END IF
          END IF               
55      CONTINUE
C       The atom is not in the couples list, accept it
        FIXSTA(I) = 1
        TOFIX     = TOFIX + 1
50    CONTINUE   

C-----Filter out the mutually exclusive atoms
C     Loop over all fixable atoms
      DO 70, I=1, FIXABL
C       Only use the atoms on the list to fix
        IF (FIXSTA(I).NE.1) GO TO 70
C       Is the atom in the exclude list
        DO 75, J=1, NTRY(4)
C         Is it the right residue and the root atom
          IF((TRES(I).EQ.EXCLUDE(J)(1:3)).AND.
     +       (FIXATM(I).EQ.EXCLUDE(J)(4:7))) THEN 
C           Loop over all atoms to find branch atoms
            DO 80, K=1, FIXABL
              IF ((FIXID(I).EQ.FIXID(K)).AND.
     +            ((FIXATM(K).EQ.EXCLUDE(J)(8:11)).OR.
     +             (FIXATM(K).EQ.EXCLUDE(J)(12:15)))) THEN
                FIXSTA(K) = 0 
                TOFIX     = TOFIX - 1
                EXFIX     = EXFIX + 1
              END IF
80          CONTINUE
          END IF               
75      CONTINUE
70    CONTINUE         


! 
!           IF (VERBOS.EQV..TRUE.) THEN
!            DO 99, I=1, NNAMFIX
!               PRINT*, NAMFIX(I)
! 99          CONTINUE
!           END IF



C-----Report the stats
      WRITE(UNIT=6,FMT=997) 'Chiral centers fixed    : ', TOFIX
      WRITE(UNIT=6,FMT=997) '- By atom swapping      : ', ATSWP
      WRITE(UNIT=6,FMT=997) '- By residue renaming   : ', NNAMFIX
      WRITE(UNIT=6,FMT=997) 'Chiral centers not fixed: ', UNFIX
      WRITE(UNIT=6,FMT=997) 'Chiral centers excluded : ', EXFIX
      WRITE(UNIT=6,FMT=997) 'Unknown chiral centers  : ', UNKNOW
      GO TO 999

C-----Error messages
988   WRITE(6,*) 'Too many problematic atoms. Increase MAXFIX.'
      FIXABL = -1
      GO TO 999
989   WRITE(6,*) 'Error reading chirality deviation'
      FIXABL = -1
      GO TO 999
990   WRITE(6,*) 'Unexpected end of file'
      FIXABL = -1
      GO TO 999
      

C-----Formats
995   FORMAT(F6.2)
996   FORMAT(F7.3)
997   FORMAT(A26,I2)
998   FORMAT(A100)

999   RETURN
      END
