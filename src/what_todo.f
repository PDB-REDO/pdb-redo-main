      PROGRAM WHAT_TODO
C=======================================================================
C  Version 2.00 2016-06-03
C  Parses pdbout.txt files from WHAT_CHECK and generates lists of 
C  residues that require attention.
C
C  Usage: what_todo (flags) PDBOUT_IN LISTFILE_OUT (EXTRACTOR_OUT)
C
C  PDBOUT_IN is the pdbout.txt file from WHAT_CHECK (V.10.0 or later)
C  LISTFILE_OUT contains lists of residues:
C  1. Resues for which the terminal Chi angle should be rotated 180 
C     degrees to improve hydrogen bonding or comply with the PDB 
C     standard.
C  2. Residues to be rebuilt forcefully.
C  EXTRACTOR_OUT is the output from the extractor tool (optional)
C
C  The '-v' option switches on verbose mode.
C
C  Written by Robbie P. Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  WHAT_TODO. Please, refer to (one of) these references:
C  - Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend, 
C    Anastassis Perrakis: "Automatic rebuilding and optimization of 
C    crystallographic structures in the Protein Data Bank" 
C    Bioinformatics, 27, p. 3392-3398 (2011)
C  - Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastassis
C    Perrakis: "PDB_REDO: constructive validation, more than just 
C    looking for errors" Acta Cryst. D68, p. 484-496 (2012)
C
C    Changelog:
C    Version 2.00:
C    - Added support for WHAT_CHECK 14.x.    
C    Version 1.04:
C    - If an EXTRACTOR_OUT file is given. Any residue marked on the 
C      skiplist will not be flagged for flipping.
C    Version 1.03:
C    - Added support for old versions of WHAT_CHECK. 
C    Version 1.02:
C    - Now compatible with WHAT_CHECK 10.x
C    Version 1.01:
C    - Added a simple help function.
C======================================================================C
      IMPLICIT NONE
C-----Declare basic variables and parameters
      INTEGER   MAXLIN, I, J, K, STATUS, ARGS, MAXRES
      CHARACTER RESIDUES*83, BACKBONE*24
      CHARACTER VERS*4
      PARAMETER (VERS='2.00')
C-----MAXLIN is the maximum allowed lines in a PDBOUT file
      PARAMETER (MAXLIN=999999)
C-----The number of residues in the NO_BLD arrays
      PARAMETER (MAXRES=150)

C-----Declare the variables and parameters
      CHARACTER LINE*80, PDBOUT*255, LISTFIL*255, EXTRACT*255, C2JUNK*2,
     +          TCHAIN*2, CHOP*2     
      INTEGER   EXTRA, FLIP_CNT, CHKGRP, GROUP, KILL_CNT, SKIP_CNT
      LOGICAL   VERBOS, PARSE

C-----Declare arrays
      CHARACTER FLIP_CHN(MAXRES)*1,  FLIP_INS(MAXRES)*1,
     +          KILL_CHN(MAXRES)*1,  KILL_INS(MAXRES)*1,
     +          SKIP_CHN(MAXRES)*1,  SKIP_INS(MAXRES)*1
      INTEGER   FLIP_RES(MAXRES), KILL_RES(MAXRES), SKIP_RES(MAXRES)

C========================End of declarations===========================C
C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   WHAT_TODO version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'what_todo parses a pdbout.txt file from WHAT_CHECK a'//
     +          'nd and creates lists of residues to be treated by Si'//
     +          'deAide.' 
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'what_todo (flags) PDBOUT_IN LISTFILE_OUT '//
     +          '(EXTRACTOR_OUT)'
      WRITE(6,*)' '
      WRITE(6,*)'PDBOUT_IN is a pdbout.txt from WHAT_CHECK.' 
      WRITE(6,*)'LISTFILE_OUT contains lists of residues that must be'//
     +          ' treated by SideAide.' 
      WRITE(6,*)'EXTRACTOR_OUT containt the output from extractor.' 
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose output. Some extra information is printed'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing WHAT_TODO:'
      WRITE(6,*)'Joosten et al. ''PDB_REDO: constructive validation, '//
     +          'more than just looking for errors'' Acta Cryst., D68'//
     +          ', 484-496 (2012)'
        GO TO 999
      END IF


C======================Work through pdbout.txt=========================C
C-----Initialise values
100   EXTRA    = 0
      FLIP_CNT = 0
      KILL_CNT = 0
      SKIP_CNT = 0
      VERBOS   = .FALSE.
      PARSE    = .FALSE.
      GROUP    = 0

C-----Really ugly!!!
      RESIDUES = 'GLY,ALA,VAL,THR,SER,CYS,TRP,PHE,TYR,HIS,ILE,LEU,ASP,AS
     +N,GLU,GLN,MET,ARG,LYS,PRO,MSE'
      BACKBONE = ' N  , CA , C  , O  , OXT'

C-----Use verbose?
      CALL GETARG(1, C2JUNK)
      IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
        WRITE(6,*) 'Verbose mode'
        VERBOS=.TRUE.
        EXTRA =1
      ELSE IF (C2JUNK(1:1).EQ.'-') THEN
        WRITE(6,*) 'Invalid option: ', C2JUNK 
        WRITE(6,*) 'It will be ignored.'
        EXTRA =1
      END IF

C-----Use extractor output?
      IF ((ARGS - EXTRA).GT.2) THEN
        STATUS=0
        CALL GETARG((3+EXTRA), EXTRACT)
        OPEN (UNIT=9, FILE=EXTRACT, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) EXTRACT, 'extractor file cannot be opened!'
          GO TO 999
        END IF

C-------Populate the skip list
        CALL MINESKIP (SKIP_CHN, SKIP_RES, SKIP_INS, SKIP_CNT) 
      ENDIF

C-----Get input PDBOUT file
      STATUS=0
      CALL GETARG((1+EXTRA), PDBOUT)
      OPEN(UNIT=7, FILE=PDBOUT, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) PDBOUT, 'PDBOUT file cannot be opened!'
        GO TO 999
      END IF

C-----Read PDBOUT file line by line
      DO 200, J=1,MAXLIN
        READ(UNIT=7, FMT=920, END=300) LINE

C-------Read until you find a title
        IF (LINE(1:1).EQ.'#') THEN
C---------Reset the group number
          GROUP = CHKGRP(LINE)
        END IF

C-------Use information in this section if usable
        IF (GROUP.EQ.0) THEN
C         Keep walking! Nothing to see here!          
        ELSE IF (GROUP.EQ.1) THEN
C         This section contains side chains to be flipped
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*) LINE
          END IF
C          IF (FLIP_CNT.EQ.MAXRES) GO TO 997
C         Only use data lines
          IF (LINE(18:18).EQ.')') THEN          
            FLIP_CNT = FLIP_CNT+1
            READ(UNIT=LINE(13:21), FMT=980) FLIP_RES(FLIP_CNT), 
     +       FLIP_INS(FLIP_CNT), TCHAIN
C           Left allign the chainID            
            TCHAIN = CHOP(TCHAIN,1)
            FLIP_CHN(FLIP_CNT) = TCHAIN(1:1)
            IF (FLIP_INS(FLIP_CNT).EQ.'-') THEN
              FLIP_INS(FLIP_CNT) = ' ' 
            END IF
          END IF
C         Remove entry if it is on the skiplist
          IF (SKIP_CNT.LT.1) GO TO 200 
          DO 210, K=1, SKIP_CNT
            IF ((FLIP_RES(FLIP_CNT).EQ.SKIP_RES(K)).AND.
     +          (FLIP_CHN(FLIP_CNT).EQ.SKIP_CHN(K)).AND.
     +          (FLIP_INS(FLIP_CNT).EQ.SKIP_INS(K))) THEN
              FLIP_CNT = FLIP_CNT-1
            END IF
210       CONTINUE
           
        ELSE IF (GROUP.EQ.2) THEN
C         This section contains side chains to be killed
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*) LINE
          END IF
C         Only use data lines
          IF (LINE(18:18).EQ.')') THEN          
            KILL_CNT = KILL_CNT+1
            READ(UNIT=LINE(13:21), FMT=980) KILL_RES(KILL_CNT), 
     +       KILL_INS(KILL_CNT), TCHAIN
C           Left allign the chainID      
            TCHAIN = CHOP(TCHAIN,1)
            KILL_CHN(KILL_CNT) = TCHAIN(1:1)
            IF (KILL_INS(KILL_CNT).EQ.'-') THEN
              KILL_INS(KILL_CNT) = ' ' 
            END IF
          END IF
        END IF
200   CONTINUE


C============================Output section============================C
C-----Get output file
300   STATUS=0                                     
      CALL GETARG((2+EXTRA), LISTFIL)
      OPEN(UNIT=8, FILE=LISTFIL, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) LISTFIL, 'cannot be created!'
        GO TO 999
      END IF
      

C-----Write out lists
      IF (VERBOS.EQV..TRUE.) THEN
        WRITE(6,*) FLIP_CNT
      END IF
      WRITE(UNIT=8, FMT='(A,I4)') '#Residues to be flipped: ',FLIP_CNT 
      CALL RES_LIST (FLIP_CNT, FLIP_CHN, FLIP_RES, FLIP_INS)
      IF (KILL_CNT.GT.0) THEN
        WRITE(UNIT=8, FMT='(A,I4)') 
     +  '#Residues to be rebuilt forcefully: ',KILL_CNT  
        CALL RES_LIST (KILL_CNT, KILL_CHN, KILL_RES, KILL_INS)
      END IF 


C-----Skip error messages
      GO TO 999
      
C===========================Formats====================================C
910   FORMAT(3F9.3,3F7.2)
915   FORMAT(I2)
920   FORMAT(A80)
925   FORMAT(A21)
930   FORMAT(F9.3)
940   FORMAT(F6.2)
945   FORMAT(A4,1X,A3,1X,A1,I4,A1,15X,A4,1X,A3,1X,A1,I4,A1)
950   FORMAT(A4)
955   FORMAT(F6.4)
960   FORMAT(A11)
961   FORMAT(A13)
965   FORMAT(I4)
970   FORMAT(A1)
975   FORMAT(I7)
980   FORMAT(I4,A1,2X,A2)
 
C===========================Error messages=============================C
997   WRITE(6,*) 'Maximum number of residues exceeded!'
      GO TO 999
C===========================End of the line============================C
999   CLOSE (7)
      CLOSE (8)
      CLOSE (9)
      END
  
C======================================================================C
C   SUBROUTINES AND FUNCTIONS                                          C
C======================================================================C

C-----------------------------------------------------------------------
C  Assigns a group number to a section of a pdbout.txt file.
C  Group numbers:
C  0 - ignore
C  1 - flip residue side-chain
C  2 - explicitly rebuild the side chain 
C-----------------------------------------------------------------------
      INTEGER FUNCTION CHKGRP(STRING)
      IMPLICIT NONE
      INTEGER   I,HASH
      CHARACTER STRING*80, STRIN2*80
C -------------------------------------

C-----Default value
      CHKGRP = 0    

C-----Strip the section number
      STRIN2 = STRING(2:80)//' '
      HASH = INDEX(STRIN2, '#')
      STRIN2 = STRIN2(HASH+2:80)
C      WRITE(6,*) STRING

C-----Check for relevant groups
      IF (STRIN2(1:38).EQ.'Warning: Arginine nomenclature problem') THEN
        CHKGRP = 1
      ELSE IF (STRIN2(1:34).EQ.'Warning: Tyrosine convention probl')THEN
        CHKGRP = 1
      ELSE IF (STRIN2(1:34).EQ.'Warning: Phenylalanine convention ')THEN
        CHKGRP = 1
      ELSE IF (STRIN2(1:34).EQ.'Warning: Aspartic acid convention ')THEN
        CHKGRP = 1
      ELSE IF (STRIN2(1:34).EQ.'Warning: Glutamic acid convention ')THEN
        CHKGRP = 1
      ELSE IF (STRIN2(1:34).EQ.'Error: His, Asn, Gln side chain fl')THEN
        CHKGRP = 1
      ELSE IF (STRIN2(1:34).EQ.'Error: HIS, ASN, GLN side chain fl')THEN
        CHKGRP = 1
      ELSE IF (STRIN2(1:34).EQ.'Error: Threonine nomenclature prob')THEN
        CHKGRP = 2
      ELSE IF (STRIN2(1:34).EQ.'Error: Isoleucine nomenclature pro')THEN
        CHKGRP = 2
      ELSE IF (STRIN2(1:34).EQ.'Warning: Leucine nomenclature prob')THEN
        CHKGRP = 2
      ELSE IF (STRIN2(1:34).EQ.'Warning: Valine nomenclature probl')THEN
        CHKGRP = 2
      ELSE
        CHKGRP = 0 
      END IF

   10 CONTINUE
   20 RETURN
      END


C-----------------------------------------------------------------------
C  This subroutine writes out a residue list 
C-----------------------------------------------------------------------
      SUBROUTINE RES_LIST(CNT, CHN, RES, INS)

      IMPLICIT NONE
C-----Declare variables
      INTEGER   CNT, RES(CNT), POS, I, LENSTR
      CHARACTER CHN(CNT)*1, INS(CNT)*1
      CHARACTER OUTPUT*1051, DP*1

      
C-----Initialise      
      DP = ':'
      READ(UNIT=DP, FMT=100) OUTPUT
      POS = 2

C-----Write empty line and finish if there are no residues in LINKs
      IF (CNT.EQ.0) THEN
        WRITE (UNIT=8, FMT='(A)')''
        GO TO 101
      END IF

C-----Compile output line
      DO 10, I=1, CNT
        IF (RES(I).GE.1000) THEN
          WRITE(UNIT=OUTPUT(POS:POS+7),FMT=99)CHN(I), RES(I), INS(I), DP
          POS = POS+7
C          WRITE(6,*) OUTPUT 
        ELSE IF (RES(I).GE.100) THEN
          WRITE(UNIT=OUTPUT(POS:POS+6),FMT=98)CHN(I), RES(I), INS(I), DP
          POS = POS+6
C          WRITE(6,*) OUTPUT
        ELSE IF (RES(I).GE.10) THEN
          WRITE(UNIT=OUTPUT(POS:POS+5),FMT=97)CHN(I), RES(I), INS(I), DP
          POS = POS+5
C          WRITE(6,*) OUTPUT
        ELSE
          WRITE(UNIT=OUTPUT(POS:POS+4),FMT=96)CHN(I), RES(I), INS(I), DP
          POS = POS+4
C          WRITE(6,*) OUTPUT
        END IF   
10    CONTINUE

C-----Write out the results
      WRITE(UNIT=8, FMT='(A)') OUTPUT(1:POS)

C-----Format statements
96    FORMAT (A1,I1,A1,A1)
97    FORMAT (A1,I2,A1,A1)
98    FORMAT (A1,I3,A1,A1)
99    FORMAT (A1,I4,A1,A1)
100   FORMAT (A281)

101   RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine reads a residue list 
C-----------------------------------------------------------------------
      SUBROUTINE MINESKIP (CHN, RES, INS, CNT)

      IMPLICIT NONE
C-----Declare variables
      INTEGER   CNT, RES(150), POS, I, LENSTR
      CHARACTER CHN(150)*1, INS(150)*1
      CHARACTER SKIPIN*1051, DP*1

      
C-----Initialise (we need the 12th line form the extractor output)
      REWIND(9)
      DO 1, I=1, 12     
        READ(UNIT=9, FMT=100) SKIPIN
1     CONTINUE
C      PRINT*, SKIPIN
      POS = 0

C-----Digest the line
      DO 10, I=1, 151
        POS = INDEX(SKIPIN, ':')
        SKIPIN = SKIPIN(POS+1:LENSTR(SKIPIN)) 
        POS = INDEX(SKIPIN, ':')
        IF (POS.EQ.0) THEN
          GO TO 101
        ELSE IF (POS.EQ.4) THEN
          CNT = CNT+1
          READ(UNIT=SKIPIN(1:POS-1),FMT=96) CHN(CNT), RES(CNT), INS(CNT)
        ELSE IF (POS.EQ.5) THEN
          CNT = CNT+1
          READ(UNIT=SKIPIN(1:POS-1),FMT=97) CHN(CNT), RES(CNT), INS(CNT)
        ELSE IF (POS.EQ.6) THEN
          CNT = CNT+1
          READ(UNIT=SKIPIN(1:POS-1),FMT=98) CHN(CNT), RES(CNT), INS(CNT)
        ELSE IF (POS.EQ.7) THEN
          CNT = CNT+1
          READ(UNIT=SKIPIN(1:POS-1),FMT=99) CHN(CNT), RES(CNT), INS(CNT)
        ELSE
          WRITE(6,*) 'Unexpected format!'
        END IF
C        PRINT*, CNT
10    CONTINUE      


C-----Format statements
96    FORMAT (A1,I1,A1)
97    FORMAT (A1,I2,A1)
98    FORMAT (A1,I3,A1)
99    FORMAT (A1,I4,A1)
100   FORMAT (A1051)

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