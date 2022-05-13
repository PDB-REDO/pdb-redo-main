      PROGRAM PICKER
C=======================================================================
C  Version 7.00 2020-06-01
C  Returns the refinement setting for the best structure based on a
C  typical input file. Dedicated software, not for general use!
C
C  Usage: picker (-f) FILENAME NTEST RFRRAT SRFRRAT (RMSZ_BOND 
C                     RMSZ_ANGLE)
C  -f        forces a selection even if not all criteria are met, 'none'
C            will not be returned
C  -s        A significant drop (5 nats) in LLFREE is required. Also an
C            R-free drop of 0.5*sigma(R-free) is required. Make sure 
C            that the refinement results are sorted by model simplicity.
C  -e        Extra strict mode. The difference between R and R-free may
C            not exceed twice the original difference.
C  -v        Verbose output.
C  FILENAME  is a file in the local directory
C  NTEST     is the numebr of test set reflections
C  RFRRAT    is the expected ratio R-free/R
C  SRFRRAT   is the sigma of the expected ratio R-free/R
C  RMSZ_BOND and RMSZ_ANGLE are custom cut-offs for geometric quality.
C            They are optional and should be used together
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
C Version 7.00
C - Model rejection is now directly on the R-free/R ratio instead of on
C   the R-free Z-score.
C - More elegant resolution of conflicts between the best R-free and the
C   best LLFREE.
C
C Version 6.02
C - Bugfix for forced model selection.
C              
C Version 6.01
C - RmsZ values for bonds and angles are now maximised at 1.5.
C - RmsZ values are completely ignored if the '-z' flag is given.  
C              
C Version 6.00
C - Rather than only looking at cases without errors, we look at cases
C   with the fewest errors.              
C
C Version 5.00
C - RMSZ values greater than 1.00 are downweighted such that the part 
C   over 1.00 is halved.
C
C Version 4.11
C - In extra strict model R-factor gaps of up to 2 percentage points 
C   are accepted even when they are more than double the intial value.
C
C Version 4.10
C - Changed definition of sigma(R-free) and adjusted the related 
C   cut-offs to 2.6 sigma.
C - Changed the output in verbose mode.
C
C Version 4.09
C - Fixed a bug that occured when two refinements have exactly the same 
C   R-free.
C - Fixed a bug for cases where LLFREE imporves but R-free doesn't.
C
C Version 4.08
C - The -s option now also requires a drop in R-free of 
C   0.5*sigma(R-free).
C
C Version 4.07
C - Added a help function.
C - Updated the reference for picker. 
C
C Version 4.6
C - Added an extra strict mode where the difference between R and Rfree 
C   is limited to twice the original difference.
C Version 4.5
C - Fixed a bug for dealing with R-free values.
C
C Version 4.4
C - When the LLFREE and R-free conflict, the R-free Z-score is leading.
C 
C Version 4.3
C - The maximum allowed R-free/R difference is now 6% unless the initial
C   difference was huge.
C - Records with R-free = 0.0 are now rejected.
C
C Version 4.2
C - When there is a conflict between the optimal R-free and the optimal
C   LLfree and only LLfree shows a significant difference, the R-free Z
C   is used to make the final decision iff R-free Z for optimal LLfree
C   is less than 0.0.
C=======================================================================
      IMPLICIT NONE
C-----Declare the basic variables and parameters
      INTEGER   I, J, K, STATUS, MAXLIN
      CHARACTER VERS*4
      PARAMETER (VERS='7.00')
C-----MAXLIN is the maximum number of refinement settings that can be
C     evaluated
      PARAMETER (MAXLIN=40)
      CHARACTER INLIST*255, INPUT*10, LINE*80, C2JUNK*2
      LOGICAL   IMPROV, NEWFRE, GDANGL, BETR, BETFRE, GDDIF, FORCE,
     +          VERBOS, SIGNIF, ADVANC, ESTRIC, NORMSZ, LENIENT      
      REAL      RFACTR, RFREER, LLFRER, RFRRAT, SRFRRAT, RMSZB, RMSZA,
     +          OLDRAT, DIFFR,  RFZR,   RFREEO, ODIFF, MAXRAT, RRAT,
     +          ERRORR, LLFRER2, RFREER2
      INTEGER   GETINT, NCOND,  ARGS,   NTEST,  EXTRA,  BLLFR, BRFREE, 
     +          BDIFF,  BRFZ,   BRAT, BLLFR2, BRFREE2 
      DOUBLE PRECISION  GTREAL
C-----Data arrays
      CHARACTER COND(MAXLIN)*4, OPTCON*4
      INTEGER   NCYC(MAXLIN)
      REAL      RFREE(MAXLIN), RFACT(MAXLIN), FOM(MAXLIN),
     +          DIFF(MAXLIN),  ANGLE(MAXLIN), CHIRAL(MAXLIN),
     +          BOND(MAXLIN),  LLG(MAXLIN),   LLFREE(MAXLIN),
     +          BONDZ(MAXLIN), ANGLEZ(MAXLIN),RFMAX(MAXLIN),
     +          RFERR(MAXLIN), RAT(MAXLIN),   ERROR(MAXLIN)
C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   Picker version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Picker selects the best refinement from a set.' 
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'picker (flags) SET_IN NTEST RFRRAT SRFRRAT (RMSZ_BON'//
     +          'D RMSZ_ANGLE)'
      WRITE(6,*)' '
      WRITE(6,*)'SET_IN data file with refinement results from Refmac.'
      WRITE(6,*)'NTEST the number to test set reflections'
      WRITE(6,*)'RFRRAT the expected ratio of R-free/R.'
      WRITE(6,*)'SRFRRAT the sigma of the expected ratio of R-free/R.'
      WRITE(6,*)'RMSZ_BOND and RMSZ_ANGLE optional cut-offs for bond '//
     +          'and angle rms Z-scores. Always use them together. Th'//
     +          'e defaults are 1.000. Higher than 1.0 values are '//
     +          'downweighted.'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode. Gives more info about the picking p'//
     +          'rocess.'
      WRITE(6,*)'-f Forced selection. Always pick a refinement even i'//
     +          'f not all quality criteria are met.'
      WRITE(6,*)'-s Significant drop. Tests for a significant drop (5'//
     +          ' nats) in LLFREE. And a 0.5*sigma(R-free) drop in R-'//
     +          'free. The refinement results must be sorted by incre'//
     +          'asing complexity or decreasing restraint weight.'
      WRITE(6,*)'-e Extra strict. The difference between R and R-free'//
     +          ' may not exceed twice the original difference.'
      WRITE(6,*)'-l Lenient. The difference between R and R-free'//
     +          ' may not exceed twice the original difference.'     
      WRITE(6,*)'-z No rmsZ. The rmsZ values for bonds and angles are'//
     +          ' not considered for model selection.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing PICKER:'
      WRITE(6,*)'Joosten et al. ''PDB_REDO: constructive validation, '//
     +          'more than just looking for errors'' Acta Cryst., D68'//
     +          ', 484-496 (2012)'
        GO TO 999
      END IF
C--------------------------- Main program -----------------------------C

C-----Initialise
      RMSZB   = 1.000
      RMSZA   = 1.000
      EXTRA   = 0
      BLLFR   = 0
      BLLFR2  = 0
      BRFREE  = 0
      BRFREE2 = 0
      BDIFF   = 0
      BRFZ    = 0
      BRAT    = 0
      ERRORR  = 0
      FORCE   = .FALSE.
      VERBOS  = .FALSE.
      SIGNIF  = .FALSE.
      IMPROV  = .FALSE.
      ADVANC  = .FALSE.
      ESTRIC  = .FALSE.
      NORMSZ  = .FALSE.
      LENIENT = .FALSE.

C-----Check for (debug) flags
      DO 10, I=1, 6
        CALL GETARG(1+EXTRA, C2JUNK)
C-------Are there any flags?
        IF (C2JUNK(1:1).NE.'-') GO TO 20

        IF ((C2JUNK.EQ.'-f').OR.(C2JUNK.EQ.'-F')) THEN
          FORCE = .TRUE.
          EXTRA = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-s').OR.(C2JUNK.EQ.'-S')) THEN
          SIGNIF = .TRUE.
          EXTRA  = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
          VERBOS = .TRUE.
          EXTRA  = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-e').OR.(C2JUNK.EQ.'-E')) THEN
          ESTRIC = .TRUE.
          EXTRA  = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-z').OR.(C2JUNK.EQ.'-Z')) THEN
          NORMSZ = .TRUE.
          EXTRA  = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-l').OR.(C2JUNK.EQ.'-L')) THEN
          LENIENT = .TRUE.
          EXTRA  = EXTRA+1  
        ELSE
	  EXTRA = EXTRA+1
        END IF
10    CONTINUE  

C-----Get input file
20    STATUS=0
      CALL GETARG(1+EXTRA, INLIST)
      OPEN(UNIT=7, FILE=INLIST, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(UNIT=*,FMT=*) 'crap'
          GO TO 999
        END IF
      REWIND (7)

C-----Get testset count, R-free/R ratio and RMSZ values
      IF (ARGS.LT.3+EXTRA) THEN
        WRITE(UNIT=*,FMT=*) 'crap'
        GO TO 999
      ELSE
        NTEST  = 1
        RFRRAT = 1.000
        CALL GETARG(2+EXTRA, INPUT)
        NTEST  = GETINT(INPUT)
        CALL GETARG(3+EXTRA, INPUT)
        RFRRAT = GTREAL(INPUT)
        CALL GETARG(4+EXTRA, INPUT)
        SRFRRAT = GTREAL(INPUT)
        IF (ARGS.GE.6+EXTRA) THEN
          CALL GETARG(5+EXTRA, INPUT)
          RMSZB = GTREAL(INPUT)
          RMSZB = (RMSZB - 1.00)/2.0 +1.00
          IF (RMSZB.LT.1.000) THEN
            RMSZB = 1.000
          ELSE IF (RMSZB.GT.1.500) THEN
            RMSZB = 1.500
          END IF
          CALL GETARG(6+EXTRA, INPUT)
          RMSZA = GTREAL(INPUT)
          RMSZA = (RMSZA - 1.00)/2.0 +1.00
          IF (RMSZA.LT.1.000) THEN
            RMSZA = 1.000
          ELSE IF (RMSZA.GT.1.500) THEN
            RMSZA = 1.500            
          END IF
        END IF
      END IF


C-----Get reference values
      RFACTR = 0.99
      RFREER = 0.99
      LLFRER = 9999999.9
      READ (UNIT=7, FMT=*) RFACTR, RFREEO
      OLDRAT = RFREEO/RFACTR
      MAXRAT = MAX(RFRRAT + 2.6*SRFRRAT, OLDRAT)
      ODIFF  = RFREEO - RFACTR      
      
      IF (VERBOS.EQV..TRUE.) THEN
        WRITE(6,*) 'Old R-free  : ', RFREEO 
        WRITE(6,*) 'Old R-free/R: ', OLDRAT
        WRITE(6,*) 'Max R-free/R: ', MAXRAT
      END IF


C-----Read log file
      NCOND = 0
      DO 100, I=1,MAXLIN
C-------Test line. Line is skipped on error!
        READ(UNIT=7, FMT=*, END=200, ERR=110)COND(I), NCYC(I),
     +  RFACT(I), RFREE(I), FOM(I), LLG(I), LLFREE(I), BOND(I),
     +  BONDZ(I), ANGLE(I), ANGLEZ(I), CHIRAL(I)

C-------Only get here if no errors!
        NCOND = NCOND + 1
        BACKSPACE(7)
        READ(UNIT=7, FMT=*) COND(NCOND), NCYC(NCOND),
     +  RFACT(NCOND), RFREE(NCOND), FOM(NCOND)  , LLG(NCOND),
     +  LLFREE(NCOND),BOND(NCOND) , BONDZ(NCOND), ANGLE(NCOND),
     +  ANGLEZ(NCOND),CHIRAL(NCOND)
        DIFF(NCOND)   = RFREE(NCOND) - RFACT(NCOND)
        RFERR(NCOND)  = RFREE(NCOND)/SQRT(1.0*NTEST) 
        IF (LENIENT.EQV..TRUE.) THEN
          RFREE(NCOND) = RFREE(NCOND) - RFERR(NCOND)
        END IF
        ERROR(NCOND)= 0
        RAT(NCOND)=RFREE(NCOND)/RFACT(NCOND)
        GO TO 100

C-------Only in case of error!
110     CALL WERROR (INLIST)

100   CONTINUE

C-----If there are no data lines, go to the end of the program 
200   IF (NCOND.EQ.0) GO TO 300

C-----Set (ridiculous) reference values
      RFREER = RFREEO
      RFREER2 = RFREEO + 0.5
      OPTCON = '0.00'
      LLFRER = 1.0E10
      LLFRER2 = 2.0E10
      DIFFR  = 1.000
      RRAT   = 10.0
      IF (FORCE.EQV..TRUE.) THEN
        RFREER = RFREE(1) 
        RFREER2 = RFREE(1) + 0.5
        BRFREE = 1
        BLLFR  = 1
        BDIFF  = 1
        BRAT   = 1
        ERRORR = 10
	    RFACTR = RFACT(1)
	    IMPROV = .TRUE.
        OPTCON = COND(1)
      END IF
      
C-----Loop over all entries entries
      DO 201, J=1,NCOND

C-------Reject all the problem conditions
        IF (NORMSZ.EQV..FALSE.) THEN
          IF (ANGLEZ(J).GT.RMSZA) THEN
            IF (VERBOS.EQV..TRUE.) THEN
              WRITE(6,*) COND(J), " Bad bond angles"
            END IF
            ERROR(J) = ERROR(J)+1
          ELSE IF (ANGLEZ(J).GT.1.000) THEN
            IF (VERBOS.EQV..TRUE.) THEN
              WRITE(6,*) COND(J), " Poor bond angles"
            END IF
            ERROR(J) = ERROR(J)+0.5
          END IF
          IF (BONDZ(J).GT.RMSZB) THEN
            IF (VERBOS.EQV..TRUE.) THEN
              WRITE(6,*) COND(J), " Bad bond lengths"
            END IF
            ERROR(J) = ERROR(J)+1
          ELSE IF (BONDZ(J).GT.1.000) THEN  
            IF (VERBOS.EQV..TRUE.) THEN
              WRITE(6,*) COND(J), " Poor bond lengths"
            END IF
            ERROR(J) = ERROR(J)+0.5
          END IF
        END IF  
	    IF (RAT(J).GT.MAXRAT) THEN
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*) COND(J), " Bad R-free/R: ", RAT(J)   
          END IF
          ERROR(J) = ERROR(J)+1
        END IF
        IF (RFREE(J).LT.0.001) THEN
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*) COND(J), " R-free corrupt"
          END IF
          ERROR(J) = ERROR(J)+1
        END IF
        IF (RFREE(J)-RFREEO.GE.0.00005) THEN
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*) COND(J), " R-free went up"
          END IF
          ERROR(J) = ERROR(J)+1
        END IF
        IF(ESTRIC.EQV..TRUE.)THEN
          IF ((DIFF(J).GT.2*ODIFF).AND.(DIFF(J).GT.0.02)) THEN
            IF (VERBOS.EQV..TRUE.) THEN
              WRITE(6,*) COND(J), " R-difference too great"
            END IF
            ERROR(J) = ERROR(J)+1
          END IF
        END IF

C       Take the entry with fewest errors
        IF (ERROR(J).LT.ERRORR) THEN
          ERRORR = ERROR(J)
          IF (FORCE.EQV..TRUE.) THEN
            BRFREE = J
            BLLFR  = J
            BDIFF  = J
            BRAT   = J
          END IF          
        END IF
        
        IF (VERBOS.EQV..TRUE.) THEN
          WRITE (6,'(A, 4(X,F7.4),X,F3.1)') COND(J), RFREE(J),  
     +    DIFF(J), RFERR(J), RAT(J), ERROR(J)
        END IF
        
201   CONTINUE
 


C-----Loop over all entries with the lowest number of errors
      DO 202, J=1,NCOND

C        PRINT*, J, ERROR(J)
C-------Is any of the scores better than the reference        
        IF ((ERROR(J).LE.ERRORR).AND.(LLFREE(J).LT.LLFRER)) THEN
C         Best cases        
          IF ((SIGNIF.EQV..TRUE.).AND.(LLFRER-LLFREE(J).GE.5.0)) THEN
            BLLFR2 = BLLFR
            BLLFR  = J
            LLFRER2 = LLFRER
            LLFRER = LLFREE(J)
            IMPROV = .TRUE.
          ELSE IF (SIGNIF.EQV..FALSE.) THEN
            BLLFR2 = BLLFR
            BLLFR  = J
            LLFRER2 = LLFRER
            LLFRER = LLFREE(J)
            IMPROV = .TRUE.
          END IF
        ELSE IF ((ERROR(J).LE.ERRORR).AND.(LLFREE(J).LT.LLFRER2)) THEN
C         Possible second best cases  
          BLLFR2 = J
          LLFRER2 = LLFREE(J)       
        END IF
C-------Does R-free improve? 
        IF ((ERROR(J).LE.ERRORR).AND.(RFREER-RFREE(J).GE.0.0)) THEN
          IF ((SIGNIF.EQV..TRUE.).AND.
     +        (RFREER-RFREE(J).GT.0.5*RFERR(J))) THEN
            BRFREE2 = BRFREE
            BRFREE = J
            RFREER = RFREE(J)
            IMPROV = .TRUE.
          ELSE IF (SIGNIF.EQV..FALSE.) THEN
            BRFREE2 = BRFREE
            BRFREE = J
            RFREER = RFREE(J)
            IMPROV = .TRUE.
          END IF
C         Special case to select the first time R-free is equal to the 
C         original R-free (RFREEO)
        ELSE IF ((ERROR(J).LE.ERRORR).AND.(RFREEO-RFREE(J).GE.0.0).AND.
     +    (BRFREE.EQ.0)) THEN
          IF ((SIGNIF.EQV..TRUE.).AND.
     +        (RFREER-RFREE(J).GT.0.5*RFERR(J))) THEN
            BRFREE2 = BRFREE
            BRFREE = J
            RFREER = RFREE(J)
            IMPROV = .TRUE.
          ELSE IF (SIGNIF.EQV..FALSE.) THEN
            BRFREE2 = BRFREE
            BRFREE = J
            RFREER = RFREE(J)
            IMPROV = .TRUE.
          END IF
        ELSE IF((ERROR(J).LE.ERRORR).AND.(RFREER2-RFREE(J).GE.0.0)) THEN
C         Second best cases        
          BRFREE2 = J
          RFREER2 = RFREE(J)
        END IF
        
        IF ((ERROR(J).LE.ERRORR).AND.(DIFF(J).LT.DIFFR)) THEN
          BDIFF = J
          DIFFR = DIFF(J)
        END IF

        IF ((ERROR(J).LE.ERRORR).AND.(RAT(J).LT.RRAT)) THEN
          BRAT = J
          RRAT = RAT(J)
        END IF
202   CONTINUE 

C-----Write out the best conditions
      IF (VERBOS.EQV..TRUE.) THEN
        WRITE(6,*) 'Best LLFREE    : ', COND(BLLFR)
        WRITE(6,*) 'Best R-free    : ', COND(BRFREE)
        WRITE(6,*) 'Best difference: ', COND(BDIFF)
        WRITE(6,*) 'Best R-free/R  : ', COND(BRAT)
      END IF


C-----Figure out the best condition
      IF (IMPROV.EQV..FALSE.) THEN
C------ Do nothing: refinement did not help
      ELSE IF (BLLFR.EQ.BRFREE) THEN
C-------Obvious case where LLFREE and R-free agree
        OPTCON=COND(BLLFR)
      ELSE IF (BRFREE.EQ.0) THEN
C       There is only improvement when looking at LLFREE
        OPTCON=COND(BLLFR)
      ELSE
C       Check for a compromise (second bests agree)
        IF(VERBOS.EQV..TRUE.)THEN
          WRITE(6,*) '2nd best LLFREE: ', COND(BLLFR2)
          WRITE(6,*) '2nd best R-free: ', COND(BRFREE2)
          WRITE(6,*) 'R-free/R (best LLFREE): ', RAT(BLLFR)
          WRITE(6,*) 'R-free/R (best R-free): ', RAT(BRFREE)
        END IF    
        IF ((BLLFR2.GT.0).AND.(BLLFR2.EQ.BRFREE2)) THEN
          IF(VERBOS.EQV..TRUE.)THEN
            WRITE(6,*) 'Compromise: take 2nd best settings'
          END IF  
          OPTCON=COND(BLLFR2)
        ELSE 
          IF(VERBOS.EQV..TRUE.)THEN
            WRITE(6,*) 'Take setting with best R-free/R'
          END IF 
          IF(RAT(BLLFR).LT.RAT(BRFREE)) THEN
            OPTCON=COND(BLLFR)
          ELSE
            OPTCON=COND(BRFREE)
          END IF
        END IF  
      END IF

                  
C-----Write  optimal condition
300   IF (IMPROV.EQV..FALSE.) THEN
        WRITE(UNIT=6,FMT='(A4)') 'none'
      ELSE IF (OPTCON.EQ.'tl10') THEN
	WRITE(UNIT=6,FMT='(A5)') 'tls10'
      ELSE
        WRITE(UNIT=6,FMT='(A4)') OPTCON
      END IF

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
15    VALUE = 1

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
C  This subroutine writes a line to the error file picker.err
C-----------------------------------------------------------------------
      SUBROUTINE WERROR(INLIST)
      IMPLICIT  NONE
      CHARACTER INLIST*255, JUNK*20
      INTEGER   I, STATUT, LENSTR

C-----Open the file and find end            
      OPEN(UNIT=8, FILE='picker.err', IOSTAT=STATUT)
      DO 10, I=1, 9999
        READ(UNIT=8, FMT='(A20)', END=20)JUNK
10    CONTINUE
      
20    BACKSPACE(8)
      WRITE(UNIT=8,FMT=*)'Error found in ',INLIST(1:LENSTR(INLIST))
      CLOSE(8)
      RETURN
      END      
