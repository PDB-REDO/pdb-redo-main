      PROGRAM LONGINUS
C=======================================================================
C  Version 1.02 2013-07-12
C  Returns outlier-filtered averages for R(-free) based on k-fold cross
C  validation results. Dedicated PDB_REDO software, not for general use!
C
C  Usage: longinus (-flags) FILENAME 
C  -v        Verbose output
C  -f        Use the 'figure skating rule' when rejecting outliers
C  FILENAME  file with lines from final refinement results from Refmac,
C            preceeded by the test set number and followd by the number
C            of reflections in each set
C
C  Written by Robbie Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  LONGINUS. Please, refer to (one of) these references:
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
C Changelog:
C Version 1.02
C - Changed the outlier cut-off to 2.6 sigma (p<0.01)
C Version 1.01
C - Changed the outlier cut-off to 2.5 sigma
C
C Version 1.00
C - Feature complete version. Added a column for the R-factor gap.
C
C Version 0.02
C - The calculation of the average and the standard deviation can now
C   use the 'figure skating rule', i.e. leave out the lowest and
C   highhest value.  
C
C Version 0.01
C - First attempt
C=======================================================================
      IMPLICIT NONE
C-----Declare the basic variables and parameters
      INTEGER   I, J, K, STATUS, MAXLIN, EXTRA
      CHARACTER VERS*4
      PARAMETER (VERS='1.02')
C-----MAXLIN is the maximum number of test sets that can be evaluated
      PARAMETER (MAXLIN=100)
      CHARACTER INTAB*255, LINE*80, C2JUNK*2
      LOGICAL   VERBOS, FIG
      INTEGER   ARGS, NSET, NCYC, SDN, IAVERN
      REAL      FOM,    LLG,   LLFREE, BOND,  BONDZ,  ANGLE, ANGLEZ, 
     +          CHIRAL, AVER,  STDEV,  AVERR, AVERRF, AVERN, SDR, SDRF,
     +          SDG, AVERG, CUTOFF
      PARAMETER (CUTOFF=2.6)
C-----Data arrays
      INTEGER   SET(MAXLIN), GOOD(MAXLIN)
      REAL      RFREE(MAXLIN), RFACT(MAXLIN), NTEST(MAXLIN), DIF(MAXLIN)
C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   Longinus version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Longines returns outlier-filtered averages for R, R-'//
     +          'free, and the number of test set reflections.'  
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'longinus (flags) FILENAME'
      WRITE(6,*)' '
      WRITE(6,*)'FILENAME file with lines from final refinement resul'//
     +          'ts preceeded by the test set number and followd by t'//
     +          'he number of reflections in each set.'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode. Gives more info about the filtering'//
     +          ' process.'
      WRITE(6,*)'-f Use the ''Figure skating rule'' when filtering ou'//
     +          'tliers. I.e. use the truncated mean (and equivalent '//
     +          'standard deviation) to filter outliers.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing LONGINUS:'
      WRITE(6,*)'Joosten et al. ''Re-refinement from deposited X-ray '//
     +          'data can deliver improved models for most PDB entrie'//
     +          's'' Acta Cryst., D65, 176-185 (2009)'
        GO TO 999
      END IF
C--------------------------- Main program -----------------------------C

C-----Initialise
      VERBOS = .FALSE.
      FIG    = .FALSE.
      NSET   = 0
      EXTRA  = 0

C-----Check for (debug) flags
      DO 10, I=1, 6
        CALL GETARG(1+EXTRA, C2JUNK)
C-------Are there any flags?
        IF (C2JUNK(1:1).NE.'-') GO TO 20

        IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
          VERBOS = .TRUE.
          EXTRA = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-f').OR.(C2JUNK.EQ.'-F')) THEN
          FIG = .TRUE.
          EXTRA = EXTRA+1
        ELSE
	  EXTRA = EXTRA+1
        END IF
10    CONTINUE  

C-----Get input file
20    STATUS=0
      CALL GETARG(1+EXTRA, INTAB)
      OPEN(UNIT=7, FILE=INTAB, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(UNIT=*,FMT=*) 'Cannot open ', INTAB
          GO TO 999
        END IF
      REWIND (7)

C-----Read data file file
      DO 100, I=1,MAXLIN
C-------Test line. Line is skipped on error!
        READ(UNIT=7, FMT=*, END=200, ERR=110)SET(I), NCYC,
     +  RFACT(I), RFREE(I), FOM, LLG, LLFREE, BOND,
     +  BONDZ, ANGLE, ANGLEZ, CHIRAL, NTEST(I)

C-------Only get here if no errors!
        NSET = NSET + 1   
        BACKSPACE(7)
        READ(UNIT=7, FMT=*) SET(NSET), NCYC , RFACT(NSET), RFREE(NSET), 
     +  FOM, LLG, LLFREE ,  BOND ,     BONDZ, ANGLE      , ANGLEZ,
     +  CHIRAL  , NTEST(NSET)
        DIF(NSET) = RFREE(NSET) - RFACT(NSET)
        GOOD(NSET) = 1
        IF (VERBOS.EQV..TRUE.) THEN     
          WRITE (6,982) SET(NSET), RFACT(NSET), RFREE(NSET), DIF(NSET),
     +    NTEST(NSET)
        END IF
        GO TO 100

C-------Only in case of error!
110     IF (VERBOS.EQV..TRUE.) THEN
          WRITE (6,*) 'Error reading dataline ', I 
        END IF

100   CONTINUE

C-----If there are no data lines, go to the end of the program 
200   IF (NSET.EQ.0) GO TO 998

C-----Filter outliers 
      CALL FILTER(NTEST, GOOD, NSET, CUTOFF, VERBOS, FIG)
      CALL FILTER(RFACT, GOOD, NSET, CUTOFF, VERBOS, FIG)
      CALL FILTER(RFREE, GOOD, NSET, CUTOFF, VERBOS, FIG)
      CALL FILTER(DIF  , GOOD, NSET, CUTOFF, VERBOS, FIG)
            
C-----Calculate final values
      AVERR  = AVER(RFACT, GOOD, NSET, .FALSE.)
      AVERRF = AVER(RFREE, GOOD, NSET, .FALSE.)
      AVERG  = AVER(DIF  , GOOD, NSET, .FALSE.)
      AVERN  = AVER(NTEST, GOOD, NSET, .FALSE.)
      SDR    = STDEV(RFACT, GOOD, NSET, AVERR , .FALSE.)
      SDRF   = STDEV(RFREE, GOOD, NSET, AVERRF, .FALSE.)
      SDG    = STDEV(DIF  , GOOD, NSET, AVERG , .FALSE.)
      SDN    = NINT(STDEV(NTEST, GOOD, NSET, AVERN, .FALSE.))
      IAVERN = NINT(AVERN)

C-----Write out the values
      WRITE(UNIT=6, FMT=980)'         R-factor  R-free  R-gap  '//
     +                      'Testset size' 
      WRITE(UNIT=6, FMT=981)'Average: ',AVERR, AVERRF, AVERG,IAVERN
      WRITE(UNIT=6, FMT=981)'Sigma  : ',SDR,   SDRF,   SDG,  SDN
C-----Format statements
980   FORMAT(A)
981   FORMAT(A,F6.4,4X,F6.4,2X,F6.4,1X,I6) 
982   FORMAT(I3,1X,2(F6.4,2X), F7.4,1X,F6.0) 

C-----Skip Error messages
      GO TO 999
      
C-----Error messages
998   WRITE(6,*) 'Could not find any valid refinement data'

C-----End of the line
999   CLOSE (7)
      END

C-------------------------- Subroutines and functions -----------------C

C-----------------------------------------------------------------------
C  This function ves the standard deviation of an array, filterd using a
C  second array
C-----------------------------------------------------------------------
      REAL FUNCTION STDEV(ARRAY, FLAGS, CNT, AVER, FRULE)
      IMPLICIT NONE
      INTEGER  I, CNT, TERMS, FLAGS(CNT)
      LOGICAL  FRULE
      REAL     SOM, ARRAY(CNT), AVER, HI, LO
C -------------------------------------
      SOM   = 0.000
      TERMS = 0
      HI    = -1000
      LO    = 10000000

      DO 10 I=1,CNT
        IF(FLAGS(I).EQ.1) THEN
          SOM = SOM + ((ARRAY(I)-AVER)*(ARRAY(I)-AVER))
          TERMS = TERMS + 1
          IF (ARRAY(I).GT.HI) THEN
            HI = ARRAY(I)
          END IF
          IF (ARRAY(I).LT.LO) THEN
            LO = ARRAY(I)
          END IF
        END IF
   10 CONTINUE
C     Apply the figure skating rule
      IF (FRULE.EQV..TRUE.) THEN
        TERMS = TERMS - 2
        SOM   = SOM - ((HI-AVER)*(HI-AVER)) - ((LO-AVER)*(LO-AVER))
      END IF

      STDEV = SQRT(SOM/(TERMS-1))

   20 RETURN
      END      
      
C-----------------------------------------------------------------------
C  This function gives the average of an array, filtered using a second 
C  array
C-----------------------------------------------------------------------
      REAL FUNCTION AVER(ARRAY, FLAGS, CNT, FRULE)
      IMPLICIT NONE
      LOGICAL  FRULE
      INTEGER  I, CNT, TERMS, FLAGS(CNT)
      REAL     SOM, ARRAY(CNT), HI, LO
C -------------------------------------
      SOM   = 0.000
      TERMS = 0
      HI    = -1000
      LO    = 10000000

      DO 10 I=1,CNT
        IF(FLAGS(I).EQ.1) THEN
          SOM = SOM + ARRAY(I)
          TERMS = TERMS + 1
          IF (ARRAY(I).GT.HI) THEN
            HI = ARRAY(I)
          END IF
          IF (ARRAY(I).LT.LO) THEN
            LO = ARRAY(I)
          END IF
        END IF
   10 CONTINUE

C     Apply the figure skating rule
      IF (FRULE.EQV..TRUE.) THEN
        TERMS = TERMS - 2
        SOM   = SOM - HI - LO 
      END IF

      AVER = SOM/TERMS

   20 RETURN
      END      

C-----------------------------------------------------------------------
C  This subroutine filters outliers from an arrey by setting a flag in
C  another array. The cut-off (in standard deviations) can be specified
C-----------------------------------------------------------------------
      SUBROUTINE FILTER(ARRAY, FLAGS, CNT, CUT, VERB, FRULE)
      IMPLICIT  NONE
      INTEGER   I, CNT, FLAGS(CNT)
      REAL      Z, CUT, LAVER, AVER, LSTDEV, STDEV, ARRAY(CNT)
      LOGICAL   VERB, FRULE

C-----Calculate the average and standard deviation
      LAVER  = AVER(ARRAY, FLAGS, CNT, FRULE)
      LSTDEV = STDEV(ARRAY, FLAGS, CNT, LAVER, FRULE)
      IF (VERB.EQV..TRUE.) THEN
        WRITE(6,99) 'Average: ', LAVER, 'Standard deviation: ', LSTDEV 
      END IF
      
      DO 10 I=1,CNT
        Z = (ARRAY(I)-LAVER)/LSTDEV
        IF(ABS(Z).GT.CUT) THEN
          FLAGS(I) = 0
          IF (VERB.EQV..TRUE.) THEN
            WRITE(6,99) 'Rejecting value: ', ARRAY(I), 'Z-score = ', Z
          END IF
        ENDIF
   10 CONTINUE     

99    FORMAT(2(A,F9.4,2X))
      RETURN
      END      
