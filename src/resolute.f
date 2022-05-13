      PROGRAM RESOLUTE
C=======================================================================
C  Version 0.03 2018-12-20
C  Returns the best resolution for a set of refinements at increasing 
C  resolution.
C
C  Usage: resolute (-v) FILE_IN (NCRIT)
C  -v        Verbose output
C  FILENAME  is a file with refinement results
C  NCRIT     the number of criteria that is allowed to worsen before a
C            resolution in rejected (optional, default = 0)      
C
C  Written by Robbie Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  PICKER. Please, refer to (one of) these references:
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
C Changelog:
C Version 0.03
C - Updated format of verbose output.
C Version 0.02
C - Resolute now compensates if no values for the free weighted R-factor
C   are given.
C Version 0.01
C - Initial version 
C=======================================================================
      IMPLICIT NONE
C-----Declare the basic variables and parameters
      INTEGER   I, J, K, STATUS, MAXLIN
      CHARACTER VERS*4
      PARAMETER (VERS='0.03')
C-----MAXLIN is the maximum number of resolution settings that can be
C     evaluated
      PARAMETER (MAXLIN=40)
      CHARACTER INFIL*255, INPUT*10, LINE*120, C2JUNK*2
      LOGICAL   VERBOS     
      REAL      RESO,   RLOW,    RHIGH,  RFLOW, RFHIGH, WRFLOW, WRFHIGH,
     +          FLLLOW, FLLHIGH, FCCLOW, FCCHIGH, BRESO
      INTEGER   NCOND,  ARGS,    EXTRA,  NCRIT, GETINT, IJUNK, REJECT,
     +          DATCHK
C-----Data arrays

C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   Resolute version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Resolute returns the best resolution for a set of re'//
     +          'finements at increasing resolution.' 
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'resolute (flags) FILE_IN (NCRIT)'
      WRITE(6,*)' '
      WRITE(6,*)'FILE_in data file with refinement results from Refmac.'
      WRITE(6,*)'NCRIT Optional number criteria that is allowed to wo'//
     +          'rsen before a resolution is rejected (default = 0).'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode. Gives more info.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing RESOLUTE:'
      WRITE(6,*)'Joosten et al. ''PDB_REDO: constructive validation, '//
     +          'more than just looking for errors'' Acta Cryst., D68'//
     +          ', 484-496 (2012)'
        GO TO 999
      END IF
C--------------------------- Main program -----------------------------C

C-----Initialise
      EXTRA  = 0
      BRESO  = 9.99
      NCRIT  = 0
      REJECT = 0
      VERBOS= .FALSE.


C-----Check for (debug) flags
      DO 10, I=1, 2
        CALL GETARG(1+EXTRA, C2JUNK)
C-------Are there any flags?
        IF (C2JUNK(1:1).NE.'-') GO TO 20

        IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
          VERBOS = .TRUE.
          EXTRA = EXTRA+1
        ELSE
	  EXTRA = EXTRA+1
        END IF
10    CONTINUE  

C-----Get input file
20    STATUS=0
      CALL GETARG(1+EXTRA, INFIL)
      OPEN(UNIT=7, FILE=INFIL, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(UNIT=*,FMT=*) 'crap'
          GO TO 999
        END IF
      REWIND (7)

C-----Get NCRIT if it is given
      IF (ARGS.GT.1+EXTRA) THEN
        CALL GETARG(2+EXTRA, INPUT)
        IJUNK  = GETINT(INPUT)
        IF (IJUNK.GE.0) THEN
          NCRIT = IJUNK
        END IF
      END IF


C-----Read log file
      NCOND = 0
      DO 100, I=1,MAXLIN
C-------Stop if there are no more lines
        READ(UNIT=7, FMT=998, END=300, ERR=901) LINE
C-------Check whether all data is given
        IF (DATCHK(LINE).EQ.11) THEN        
          READ(UNIT=LINE, FMT=*, END=300, ERR=100) RESO, RLOW, RHIGH, 
     +    RFLOW,RFHIGH,WRFLOW,WRFHIGH,FLLLOW, FLLHIGH, FCCLOW, FCCHIGH
        ELSE IF (DATCHK(LINE).EQ.10) THEN
          READ(UNIT=LINE, FMT=*, END=300, ERR=100)RESO,RLOW,RHIGH,RFLOW,
     +    RFHIGH, WRFLOW, FLLLOW, FLLHIGH, FCCLOW, FCCHIGH
          WRFHIGH = -1.00
        ELSE IF (DATCHK(LINE).EQ.9) THEN
          READ(UNIT=LINE, FMT=*, END=300, ERR=100)RESO,RLOW,RHIGH,RFLOW,
     +    RFHIGH, FLLLOW, FLLHIGH, FCCLOW, FCCHIGH
          WRFHIGH = -1.00
          WRFLOW  = -1.00
        ELSE
C         Datafile is corrupt
          GO TO 902
        END IF

        NCOND  = NCOND + 1
C       Reset the rejection score
        REJECT = 0

C-------Check the results
        IF (VERBOS) THEN
          WRITE(6,*) ' ' 
          WRITE(UNIT=6, FMT=997) 'Resolution ', RESO
          WRITE(6,998) '===============' 
          WRITE(UNIT=6, FMT=998) '                   Lower     Higher'
          WRITE(UNIT=6, FMT=996) 'R-factor        :  ',RLOW,  RHIGH
          WRITE(UNIT=6, FMT=996) 'R-free          :  ',RFLOW, RFHIGH
          WRITE(UNIT=6, FMT=996) 'Weighted R-free :  ',WRFLOW, WRFHIGH
          WRITE(UNIT=6, FMT=994) 'LL-free         : ',FLLLOW, FLLHIGH
          WRITE(UNIT=6, FMT=996) 'free corr. coef.:  ',FCCLOW, FCCHIGH
        END IF
      
        IF (RFHIGH.GT.RFLOW) THEN
          REJECT = REJECT + 1
          IF (VERBOS) THEN
            WRITE(6,'(A)') 'R-free deteriorated'
          END IF
        END IF

        IF (WRFHIGH.GT.WRFLOW) THEN
          REJECT = REJECT + 1
          IF (VERBOS) THEN
            WRITE(6,'(A)') 'Weighted R-free deteriorated'
          END IF
        END IF

        IF (FLLHIGH.GT.FLLLOW) THEN
          REJECT = REJECT + 1
          IF (VERBOS) THEN
            WRITE(6,'(A)') 'LL-free deteriorated'
          END IF
        END IF

        IF (FCCHIGH.LT.FCCLOW) THEN
          REJECT = REJECT + 1
          IF (VERBOS) THEN
            WRITE(6,'(A)') 'Free correlation coefficient deteriorated'
          END IF
        END IF

C-------Check the final score
        IF (REJECT.GT.NCRIT) THEN
C         Stop. This resolution is worse than the previous one. 
          GO TO 300
        ELSE
C         Keep this resolution and continue.
          BRESO = RESO
        END IF

100   CONTINUE




                  
C-----Write  optimal condition
300   WRITE(UNIT=6,FMT=995) BRESO
      GO TO 999

C-----Error messages
901   WRITE(6,*) 'Cannot read data file! Stopping.' 
      GO TO 999
902   WRITE(6,*) 'Data file is not complete! Stopping.'

C-----Format statements
994   FORMAT(A,F8.1,X,F8.1)
995   FORMAT(F4.2)
996   FORMAT(A,F7.4,2X,F7.4)
997   FORMAT(A,F4.2)
998   FORMAT(A)

C-----End of the line
999   CLOSE (7)
      END

C-------------------------- Subroutines and functions -----------------C
C-----------------------------------------------------------------------
C  This function extracts a real number from a string
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GTREAL(TLINE)

      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*10 TLINE
      INTEGER      SPACE
      DOUBLE PRECISION VALUE

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
C  This function extracts an integer from column CLMNR from a text line
C-----------------------------------------------------------------------
      INTEGER FUNCTION GETINT(LINE)
      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*10 LINE
      INTEGER   SPACE, VALUE

C-----Initialise
      VALUE = 1
      SPACE = INDEX(LINE, ' ')
      IF (SPACE.NE.0) THEN
        READ(UNIT=LINE(1:SPACE), FMT=*, ERR=15) VALUE
      ELSE
        READ(UNIT=LINE, FMT=*, ERR=15) VALUE
      END IF

C-----Skip fallback value
      GO TO 99

C-----In case of error use fall-back value	
15    VALUE = -9999

C-----Assign and return
99    GETINT = VALUE
      RETURN
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
C  This function returns the number of data entries in a line
C-----------------------------------------------------------------------
      INTEGER FUNCTION DATCHK(LINE)
      IMPLICIT NONE
      CHARACTER LINE*(*)
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

