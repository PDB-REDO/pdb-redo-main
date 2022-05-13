      PROGRAM BSELECT
C=======================================================================
C  Version 3.00 2020-06-08
C  Compares two different types of B-factor refinement and returns the 
C  best one. Dedicated software, not for general use!
C
C  Written by Robbie Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  BSELECT. Please, refer to (one of) these references:
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
C  Version history:
C  3.00:
C  - Moved to using a new empirical model for the R-free/R ratio.
C  - The R-free Z-score is now replaced by the R-free/R Z-score in 
C    decision making.
C  - Always using the weighted R-factor now for the Hamilton test.
C  2.11:
C  - More tweaking of the R-factor gap test.
C  2.10:
C  - The overrule score in the R-factor gap test is now slightly more 
C    increased for the choice between OVER and ISOT. 
C  2.09: 
C  - Changed definition of sigma R-free.
C  - Lowered the R-free Z-score cut-off to -2.6 accordingly (which is 
C    actually more liberal). It's a parameter now.
C  - Small robustness fix.
C  2.08:
C  - If the weighted R-factors are insanely high, they are no longer
C    used.
C  2.07: 
C  - Changed cut-offs of hamilton percentages from 30% to 15% for always
C    reject and 95% to 90% for always accept.
C  2.06:
C  - The Hamilton test is overruled if the exchange rate > 1 AND the R-
C    free difference is smaller for the complex model AND the drop in R-
C    free is > 0.5*sigma(R-free)
C  2.05:
C  - The minimum values of the weights are now always 0.0. This 
C    obsoletes the '-z' switch.
C  - Used more fancy printing.
C  - Separated the R-free Z-score test and the R-free difference test.
C  2.04:
C  - If both methods seem overfitted in the Z-score test, the exchange 
C    rate is decisive.
C  - Explicit handling of cases wher R increases in the complex model,
C    but R-free drops.
C  2.03:
C  - Using weighted R as fallback if weighted R-free is not available.
C  2.02:
C  - Now only cases with a diffrenet numebr of parameters are evaluated.
C  2.01:
C  - Now explicitly deals with equivalent refinements in TLS mode.
C  - Bugfix for refinements with the same number of TLS groups, but with
C    different group selections.
C  - Updated the references.
C  2.00:
C  - Implemented limit case Hamilton test for TLS groups.
C  1.06:
C  - Negative degrees of freedom are rejected outright.
C  - Crushed a bug cause by the new feature in version 1.05.
C  1.05:
C  - If the regular R-free goes up, the simplest model is chosen 
C    immediately.
C  1.04:
C  - Switched to the proper weighted R-factor.
C  1.03:
C  - Added a simple help function.
C  1.2:
C  - Increased the number of lines read for the the restraint table.
C  1.1:
C  - Added the -z switch which sets W1MIN and W2MIN to zero.
C  1.0:
C  - First version added to PDB_REDO.
C  0.2:
C  - The R-free Z-score is now also used for non-obvious cases.
C  - Added the exchange rate test to avoid overfitting somewhat.
C  0.1: 
C  - First atempt. The Hamilton test is used to select the best B-factor
C    refinement type.
C
C  Usage: bselect (-v -t) LOG1 LOG2
C  -v        Verbose output
C  -t        TLS mode: compare different TLS models 
C  LOG1      is the first Refmac log file
C  LOG2      is the second Refmac log file 
C=======================================================================
      IMPLICIT NONE
C-----Declare the basic variables and parameters
      INTEGER   I, J, K, T, N, STATUS, MAXLIN, STEPS
      CHARACTER VERS*4
      PARAMETER (VERS='3.00')
C-----MAXLIN is the maximum number of lines in the logfile
      PARAMETER (MAXLIN=10000)
      CHARACTER LOG1*255, LOG2*255, LINE*100, C2JUNK*2, BTYPEF*6
      LOGICAL   VERBOS  , TLSMOD, TLSREF  
C-----Single values
      INTEGER   ARGS,   EXTRA,  DPAR,  DREST, PASSED,   RTABLE, SCORE,
     +          OVERRL
C     ARGS  : Number of command line arguments
C     EXTRA : Number of debug flags detected
      REAL      RATIO,  RATIO2, W1MIN,  W1MAX,  W2MIN, W2MAX, HAMILT,
     +          HAMILN, ACCEPT, RFRRAT, SRFRRAT, DLIMIT, EXCHR, ZCUTOFF
C-----Arrays
      REAL      RFACT(2), RFREE(2), RFACTW(2), RFREEW(2), DF(2), 
     +          RRAT(2),  SRRAT(2), RRATZ(2),  SIGFR(2)  
      INTEGER   PARPA(2), NATOM(2), NPAR(2),   NREF(2),   NFREE(2),
     +          NREST(2), NTLS(2)
C     PARPA: Parameters per atom
C     NATOM: Number of atoms
C     NPAR : Number of model parameters
C     NREF : Number of reflections in the work set
C     NFREE: Number of reflections in the test set
C     NREST: number of restraints
C     NTLS : Number of TLS groups
      CHARACTER BTYPE(2)*6, BTEST(4)*6

C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   BSELECT version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Bselect finds the best B-factor model based on two R'//
     +          'efmac log files.' 
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'bselect (flags) LOG_IN1 LOG_IN2'
      WRITE(6,*)' '
      WRITE(6,*)'LOG_IN1 is a Refmac log file from a refinement with '//
     +          'a specific B-factor model.'
      WRITE(6,*)'LOG_IN2 a log file from a second refinement with an '//
     +          'alternative B-factor model.'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode. Gives output for all the model sele'//
     +          'ction tests.'
      WRITE(6,*)'-t TLS mode. Compares different TLS models. Returns'//
     +          ' ''LOG1'' for the model in LOG_IN1 and ''LOG2'' for '//
     +          ' the other model.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing BSELECT:'
      WRITE(6,*)'Joosten et al. ''PDB_REDO: constructive validation, '//
     +          'more than just looking for errors'' Acta Cryst., D68'//
     +          ', 484-496 (2012)'
        GO TO 999
      END IF
C--------------------------- Main program -----------------------------C

C-----Initialise
      EXTRA    = 0
      NREST(1) = 0
      NREST(2) = 0
      PARPA(1) = -1
      PARPA(2) = -1
      NTLS(1)  = 0
      NTLS(2)  = 0
      VERBOS   = .FALSE.
      TLSMOD   = .FALSE.
      TLSREF   = .FALSE.
      BTYPEF   = 'NULL  '
      SCORE    = 0
      OVERRL   = 0
      ZCUTOFF  = -2.6


C-----Check for (debug) flags
      DO 10, I=1, 3
        CALL GETARG(1+EXTRA, C2JUNK)
C-------Are there any flags?
        IF (C2JUNK(1:1).NE.'-') GO TO 20
C       Is it a valid flag 
        IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
          VERBOS = .TRUE.
          EXTRA  = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-t').OR.(C2JUNK.EQ.'-T')) THEN
          TLSMOD = .TRUE.
          EXTRA  = EXTRA+1
        ELSE
          EXTRA = EXTRA+1
        END IF
10    CONTINUE  

C-----Get input files
20    STATUS=0
C     First log file
      CALL GETARG(1+EXTRA, LOG1)
      OPEN(UNIT=7, FILE=LOG1, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) 'Cannot open ', LOG1 
          GO TO 999
        END IF
      REWIND (7)
C     Second log file    
      CALL GETARG(2+EXTRA, LOG2)
      OPEN(UNIT=8, FILE=LOG2, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) 'Cannot open ', LOG2
          GO TO 999
        END IF
      REWIND (8)

C-------------------------- Parse log files ----------------------------

C-----Loop over log files
      DO 30, I=1, 2
        RTABLE = 0
C       Loop over lines in each log file
        DO 40, J=1, MAXLIN
          READ(UNIT=6+I, FMT=900, END=30) LINE
          IF (TLSMOD.EQV..TRUE.) THEN
C           Just set BTYPE to logfile number
            IF (I.EQ.1) THEN
              BTYPE(I) = 'LOG1  '
            ELSE
              BTYPE(I) = 'LOG2  '
            END IF   
          ELSE 
C           Get B-factor refinement type
C Posibilities:
C  Refinement of individual isotropic Bfactors
C  Refinement of individual anisotropic Bfactors
C  Refinement of one overall Bfactor
CNumber of used reflections           =      19463
C             Restraint type 
CBond angles  : refined atoms
CResolution limits
C
C12345678901234567890123456789012345678901234567890
            IF(INDEX(LINE, ' Refinement of').NE.0) THEN
              IF(LINE(17:30).EQ.'individual iso') THEN
                PARPA(I) = 4
                BTYPE(I) = 'ISOT  '
              ELSE IF(LINE(17:30).EQ.'individual ani') THEN
                PARPA(I) = 9
                BTYPE(I) = 'ANISOT'
              ELSE IF(LINE(17:30).EQ.'one overall Bf') THEN
                PARPA(I) = 3
                BTYPE(I) = 'OVER  '
              ELSE
                WRITE(6,*) 'Cannot find B-factor refinement type'
                GO TO 999
              END IF
            END IF
          ENDIF

C---------Number of TLS groups
          IF (LINE(1:8).EQ.'  Group:') THEN
            READ(LINE(9:12), FMT=*,ERR=996) NTLS(I)
          END IF

C---------Number of atoms
          IF(LINE(1:22).EQ.'  Number of atoms    :') THEN
            READ(LINE(23:30), FMT=*, ERR=998) NATOM(I)
          END IF


C---------Number of restraints (read the second full table)
          IF (TLSMOD.EQV..FALSE.) THEN
            IF (LINE(2:23).EQ.'WRITTEN OUTPUT MTZ FIL') THEN
C             Done with refinement, the next table is the right one
              RTABLE = 1
            END IF
            IF ((LINE(1:27).EQ.'             Restraint type').AND.
     +        (RTABLE.EQ.1))THEN
              CALL RESCNT(6+I, NREST(I))
              IF(NREST(I).LT.0) THEN
                WRITE(6,*) 'Cannot read the number of restraints'
                GO TO 999
              END IF
            END IF 
          END IF  

C---------R-free reflections (read every table, keep the last values) 
          IF (LINE(1:7).EQ.'NR_free') THEN
C           Read ahead one line
	    READ(UNIT=6+I, FMT=900, END=30) LINE
	    CALL GETNFR(6+I, NFREE(I))
	    IF(NFREE(I).LT.0) THEN
	      WRITE(6,*) 'Cannot read number of test set reflections'
	      GO TO 999
	    END IF  
          END IF 

C---------R values and number of reflections (read every table, keep the
C         last values)
          IF (LINE(1:17).EQ.'Resolution limits') THEN
	    CALL GETRS(6+I, RFACT(I), RFREE(I), RFACTW(I), RFREEW(I),
     +                 NREF(I))
            IF(RFACT(I).LT.0) THEN
              WRITE(6,*) 'Cannot read the R-values'
             GO TO 999
            END IF  
          END IF 

C---------Are TLS tensors refined?
          IF (LINE(5:24).EQ.'TLS refinement cycle') THEN   
            TLSREF = .TRUE.
          END IF     
        
40      CONTINUE
30    CONTINUE

C-----Number of parameters  
      IF (TLSMOD.EQV..TRUE.) THEN
       NPAR(1) = NTLS(1)*20
       NPAR(2) = NTLS(2)*20
      ELSE 
C      Only count the TLS parameters if they were refined here
       IF (TLSREF.EQV..FALSE.) THEN
         NTLS(1) = 0
         NTLS(2) = 0
       END IF
       NPAR(1) = NATOM(1)*PARPA(1)+NTLS(1)*20
       NPAR(2) = NATOM(2)*PARPA(2)+NTLS(2)*20
      END IF

C-----Give a summary if verbos
      IF(VERBOS.EQV..TRUE.)THEN
        WRITE(UNIT=6,FMT=905) 'Summary:'
        WRITE(UNIT=6,FMT=905) '--------'
        IF (TLSMOD.EQV..TRUE.) THEN
          WRITE(UNIT=6,FMT=906) 'TLS groups in log1 :', NTLS(1)
          WRITE(UNIT=6,FMT=906) 'TLS groups in log2 :', NTLS(2)  
        ELSE
          WRITE(UNIT=6,FMT=906) 'Atoms in log1      :', NATOM(1)
          WRITE(UNIT=6,FMT=906) 'Atoms in log2      :', NATOM(2)
        END IF
        WRITE(UNIT=6,FMT=906) 'Parameters in log1 :', NPAR(1)
        WRITE(UNIT=6,FMT=906) 'Parameters in log2 :', NPAR(2)
        WRITE(UNIT=6,FMT=906) 'Work refs. in log1 :', NREF(1)
        WRITE(UNIT=6,FMT=906) 'Work refs. in log2 :', NREF(2)
        WRITE(UNIT=6,FMT=906) 'Test refs. in log1 :', NFREE(1)
        WRITE(UNIT=6,FMT=906) 'Test refs. in log2 :', NFREE(2)
        IF (TLSMOD.EQV..FALSE.) THEN
          WRITE(UNIT=6,FMT=906) 'Restraints in log1 :', NREST(1)
          WRITE(UNIT=6,FMT=906) 'Restraints in log2 :', NREST(2)
        END IF
        WRITE(UNIT=6,FMT=904) 'Overal R in log1   :', RFACT(1)
        WRITE(UNIT=6,FMT=904) 'Overal R in log2   :', RFACT(2)
        WRITE(UNIT=6,FMT=904) 'Free R in log1     :', RFREE(1)
        WRITE(UNIT=6,FMT=904) 'Free R in log2     :', RFREE(2)
        WRITE(UNIT=6,FMT=904) 'Weighted R in log1 :', RFACTW(1)
        WRITE(UNIT=6,FMT=904) 'Weighted R in log2 :', RFACTW(2)
        WRITE(UNIT=6,FMT=904) 'Weighted Free R l1 :', RFREEW(1)
        WRITE(UNIT=6,FMT=904) 'Weighted Free R l2 :', RFREEW(2)  
        WRITE(6,*) ' ' 
      END IF

C---------------------------- Hamilton test ----------------------------

C-----Decide which is which. N is the refinement with most parameters
      IF(NPAR(1).GT.NPAR(2)) THEN
        T = 2
        N = 1
      ELSE IF (NPAR(1).EQ.NPAR(2)) THEN
C       Special case. Here the Hamilton test won't work.
        BTYPEF = 'BOTH  '
        IF (VERBOS.EQV..TRUE.) THEN
         WRITE(6,*) ' '
         WRITE(6,*) 'The models have the same number of parameters.'
         WRITE(6,*) 'Cowardly refusing to select the best model.'
        END IF
        GO TO 899
      ELSE 
        T = 1
        N = 2
      END IF

C-----Get the ratio 
      IF ((ABS(RFACTW(T)).GE.1.00).OR.(ABS(RFACTW(N)).GE.1.00)) THEN
          RATIO = -1
      ELSE
          RATIO = RFACTW(T)/RFACTW(N)
      END IF
      RATIO2 = RATIO*RATIO
      IF(VERBOS.EQV..TRUE.)THEN
        WRITE(UNIT=6,FMT=905) 'Hamilton test:'
        WRITE(UNIT=6,FMT=905) '--------------'
        IF (RATIO .LT. 0) THEN  
          WRITE(UNIT=6,FMT=905) 'Warning: strange weigthed R-factors!'
          WRITE(UNIT=6,FMT=905) 'Hamilton test skipped.'
        ELSE
          WRITE(UNIT=6,FMT=904) 'Weighted R ratio  :', RATIO
          WRITE(UNIT=6,FMT=904) 'Squared  R ratio  :', RATIO2
        END IF
      END IF

C-----Simple case 1: the Hamilton test cannot be performed
      IF (RATIO.LT.0) THEN
        BTEST(1) = 'NULL  '
        ACCEPT   = -0.999
        GO TO 80
      END IF

C-----Simple case 2: the ratio < 1 --> adding parameters does not work 
C     Choose the refinment with the fewest parameters   
      IF(RATIO.LT.1.000) THEN
        BTYPEF = BTYPE(T)
        IF(VERBOS.EQV..TRUE.)THEN
          WRITE(6,*) ' '
          WRITE(6,*) 'The weighted R went up.'
          WRITE(6,*) 'Taking the simplest model:'
        END IF
        GO TO 899 
      END IF

C-----Simple case 3: the regular R-free goes up
      IF(RFREE(N).GT.RFREE(T)) THEN
        BTYPEF = BTYPE(T)
        IF (VERBOS.EQV..TRUE.) THEN
          WRITE(6,*) ' '
          WRITE(6,*) 'R-free went up.'
          WRITE(6,*) 'Taking the simplest model:'
        END IF
        GO TO 899 
      END IF


C-----Calculate the ranges of the weights
      DPAR  = NPAR(N) - NPAR(T)
      DREST = 0

      IF (TLSMOD.EQV..TRUE.) THEN
        W1MIN = 1
        W1MAX = 1
        W2MIN = 1
        W2MAX = 1
      ELSE
C       DF1 = NREF + W1*NREST - NPAR
        W1MIN = 0.0
        W1MAX = 1.0*NPAR(T)/NREST(T)
        IF (W1MAX.GT.1.00) THEN
          W1MAX = 1.0
        END IF

C       DF2 = NREF + W1*NREST(old) - NPAR(old) + w2*DREST - DPAR
        DREST = NREST(N) - NREST(T)
        W2MIN = 0.0
        W2MAX = 1.0*DPAR/DREST
        IF (W2MAX.GT.1.00) THEN
          W2MAX = 1.0
        END IF
      END IF

C     Report
      IF (VERBOS.EQV..TRUE.) THEN 
        WRITE(UNIT=6,FMT=906) 'Extra parameters  :', DPAR
        WRITE(UNIT=6,FMT=906) 'Extra restraints  :', DREST
        IF (TLSMOD.EQV..FALSE.) THEN
          WRITE(UNIT=6,FMT=907) 'W1 range    :', W1MIN, W1MAX
          WRITE(UNIT=6,FMT=907) 'W2 range    :', W2MIN, W2MAX
        END IF
      END IF

C     

C-----Check the possible outcome of the Hamilton test for each w1 and w2
      STEPS  = 0
      PASSED = 0
      DO 50, I=0, 100
C       Break loop when all acceptable values of W1 are tested
        IF((W1MIN+I*0.01).GT.W1MAX) GO TO 70
        DO 60, J=0, 100
C         Break loop when all acceptable values of W2 are tested
          IF((W2MIN+J*0.01).GT.W2MAX) GO TO 50
          STEPS  = STEPS+1
          HAMILT = NREF(T)-NPAR(T)+(W1MIN+I*0.01)*NREST(T) 
          HAMILN = HAMILT -DPAR + (W2MIN+J*0.01)*DREST
          IF((RATIO2.GT.(HAMILT/HAMILN)).AND.
     +       (HAMILT.GT.0.0).AND.(HAMILN.GT.0.0)) THEN
            PASSED = PASSED+1    
          END IF
C          IF(VERBOS.EQV..TRUE.)THEN         
C            WRITE(UNIT=6,FMT=903) 'STEP:', STEPS, (W1MIN+I*0.01), 
C     +      (W2MIN+J*0.01), HAMILT/HAMILN, HAMILT, HAMILN
C          END IF
60      CONTINUE
50    CONTINUE

70    ACCEPT = 1.0*PASSED/STEPS
      IF (ACCEPT.LT.0.15) THEN
        BTEST(1) = BTYPE(T)
      ELSE IF (ACCEPT.GT.0.90) THEN
        BTEST(1) = BTYPE(N)
      ELSE
        BTEST(1) = 'NULL  '
      END IF

80    IF(VERBOS.EQV..TRUE.)THEN  
        IF (TLSMOD.EQV..TRUE.) THEN
          WRITE(UNIT=6,FMT=908) 'Best TLS model in : ', BTEST(1)
        ELSE     
          WRITE(UNIT=6,FMT=902) 'Weights acceptable: ', ACCEPT*100,'%'
          WRITE(UNIT=6,FMT=908) 'Best B-factor type:  ', BTEST(1) 
        END IF
      END IF

C--------------- R-free/R ratio test -----------------------------------

      IF (TLSMOD.EQV..FALSE.) THEN
C-----Calculate R-free/R ratio, sigma(R-free) and Z(R-free)
        RRAT(N)   = RFRRAT(NREF(N), NATOM(N), PARPA(N))
        RRAT(T)   = RFRRAT(NREF(T), NATOM(T), PARPA(T))
        SRRAT(N)  = SRFRRAT(NREF(N), NATOM(N), PARPA(N))
        SRRAT(T)  = SRFRRAT(NREF(T), NATOM(T), PARPA(T))
        RRATZ(N) = (RRAT(N)-(RFREE(N)/RFACT(N)))/SRRAT(N)
        RRATZ(T) = (RRAT(T)-(RFREE(T)/RFACT(T)))/SRRAT(T)

C-------Is the complex model (N) overfitted?
C       Reject the test if the R-free/R ratio cannot be calculated      
        IF ((RRAT(N).GE.5.0).OR.(RRAT(T).GE.5.0)) THEN
          BTEST(2) = 'NULL  '
        ELSE
          IF (RRATZ(N).GE.ZCUTOFF) THEN
C           Complex model NOT overfitted
            BTEST(2) = BTYPE(N)
          ELSE IF (RRATZ(N).GE.RRATZ(T)) THEN
C           Complex model less overfitted than simple model
            BTEST(2) = BTYPE(N)
          ELSE
            BTEST(2) = BTYPE(T)
          END IF
        END IF

C-------Report
        IF(VERBOS.EQV..TRUE.)THEN
          WRITE(6,*) ' '
          WRITE(UNIT=6,FMT=905) 'R-free/R Z-score test:'
          WRITE(UNIT=6,FMT=905) '--------------------'
          IF((RRAT(1).GT.5.0).OR.(RRAT(2).GT.5.0)) THEN
            WRITE(UNIT=6,FMT=905)'Expected R-free/R in log1:   N/A'
            WRITE(UNIT=6,FMT=905)'Expected R-free/R in log2:   N/A'
            WRITE(UNIT=6,FMT=904)'sigma(R-free/R) in log1  :   N/A'
            WRITE(UNIT=6,FMT=904)'sigma(R-free/R) in log2  :   N/A'
            WRITE(UNIT=6,FMT=905)'R-free/R Z-score in log1 :   N/A'
            WRITE(UNIT=6,FMT=905)'R-free/R Z-score in log2 :   N/A'
          ELSE
           WRITE(UNIT=6,FMT=904)'Expected R-free/R in log1:', RRAT(1)
           WRITE(UNIT=6,FMT=904)'Expected R-free/R in log2:', RRAT(2)
           WRITE(UNIT=6,FMT=904)'sigma(R-free/R) in log1  :', SRRAT(1)
           WRITE(UNIT=6,FMT=904)'sigma(R-free/R) in log2  :', SRRAT(2)
           WRITE(UNIT=6,FMT=904)'R-free/R Z-score in log1 :', RRATZ(1)
           WRITE(UNIT=6,FMT=904)'R-free/R Z-score in log2 :', RRATZ(2)
          END IF       
          WRITE(UNIT=6,FMT=909) 'R-free Z-score cut-off   :', ZCUTOFF
          WRITE(UNIT=6,FMT=908) 'Best B-factor type is    : ', BTEST(2)
        END IF


C-------------- R-free - R difference test -----------------------------

C-----Let the maximal R-free - R difference depend on the B-factor type
        IF (BTYPE(N).EQ.'ANISOT') THEN
          DLIMIT = 0.04
        ELSE
          DLIMIT = 0.06
        END IF

C-----Test the difference
        IF ((RFREE(N)-RFACT(N)).GT.DLIMIT) THEN
C         Difference too large; see if the difference is better in 
C         the simple model
          IF ((RFREE(T)-RFACT(T)).GT.(RFREE(N)-RFACT(N))) THEN
            BTEST(3) = 'NULL  '
          ELSE
            BTEST(3) = BTYPE(T)
          END IF
        ELSE
          BTEST(3) = BTYPE(N)
        END IF

C-------Increase the overrule score by one if the complex model has a 
C       similar R-factor gap
        IF ((RFREE(N)-RFACT(N)).LT.0.5*DLIMIT) THEN
          OVERRL = OVERRL + 1
        ELSE IF (BTYPE(N).EQ.'ANISOT') THEN
          IF ((RFREE(T)-RFACT(T)).GT.(RFREE(N)-RFACT(N))) THEN
            OVERRL = OVERRL + 1
          END IF
        ELSE
          IF ((RFREE(T)-RFACT(T)).GT.0.95*(RFREE(N)-RFACT(N))) THEN
            OVERRL = OVERRL + 1
          END IF
        END IF  
          

C-------Report
        IF(VERBOS.EQV..TRUE.)THEN
          WRITE(6,*) ' '
          WRITE(UNIT=6,FMT=905) 'R-free difference test:'
          WRITE(UNIT=6,FMT=905) '--------------------'
          WRITE(UNIT=6,FMT=904) 'Difference in log1:',RFREE(1)-RFACT(1)
          WRITE(UNIT=6,FMT=904) 'Difference in log2:',RFREE(2)-RFACT(2)
          WRITE(UNIT=6,FMT=904) 'Difference cut-off:', DLIMIT
          WRITE(UNIT=6,FMT=908) 'Best B-factor type: ', BTEST(3)
        END IF
      END IF

C------------------------- R-free drop test ----------------------------

C-----Calculate exchange rate
      EXCHR = (RFREE(T)-RFREE(N))/ABS(RFACT(T)-RFACT(N))
      SIGFR(N) = RFREE(N)/SQRT(1.0*NFREE(N))
      SIGFR(T) = RFREE(T)/SQRT(1.0*NFREE(T))


C-----Increase the overrule score by one if the exchange rate > 1 AND 
C     the drop in R-free > 0.5*sigma(R-free)
      IF (BTYPE(N).EQ.'ANISOT') THEN
        IF ((EXCHR.GT.1).AND.((RFREE(T)-RFREE(N)).GT.0.5*SIGFR(N))) THEN
          OVERRL = OVERRL + 1
        END IF  
      ELSE  
       IF((EXCHR.GT.0.85).AND.((RFREE(T)-RFREE(N)).GT.0.5*SIGFR(N)))THEN
          OVERRL = OVERRL + 1
        END IF  
      END IF  
        
C-----Is the exchange rate acceptable      
      IF (EXCHR.LT.0.5) THEN
C       R-free went down much less than R
        BTEST(4) = BTYPE(T)
      ELSE IF ((RFREE(T).EQ.RFREE(N)).AND.(RFACT(T).EQ.RFACT(N))) THEN
        BTEST(4) = 'BOTH  '
      ELSE
        BTEST(4) = BTYPE(N)
      END IF



      IF(VERBOS.EQV..TRUE.)THEN
        WRITE(6,*) ' '
        WRITE(UNIT=6,FMT=905) 'Exchange rate test:'
        WRITE(UNIT=6,FMT=905) '-------------------'
        WRITE(UNIT=6,FMT=904) 'Drop in R-factor  :', RFACT(T)-RFACT(N)
        WRITE(UNIT=6,FMT=904) 'Drop in R-free    :', RFREE(T)-RFREE(N)
        WRITE(UNIT=6,FMT=904) 'Exchange rate     :', EXCHR
        WRITE(UNIT=6,FMT=904) 'Rate cut-off      :', 0.5000
        IF (TLSMOD.EQV..TRUE.) THEN
          WRITE(UNIT=6,FMT=908) 'Best TLS model in : ', BTEST(4)
        ELSE
          WRITE(UNIT=6,FMT=908) 'Best B-factor type: ', BTEST(4)
        END IF
      END IF


C------------------- Decide which refinement is best -------------------

C-----Use the conclusion from the Hamilton test if it is conclusive or 
100   IF(BTEST(1).NE.'NULL  ') THEN
        BTYPEF = BTEST(1)
C       Overrule a conclusive Hamilton test if the overrule score = 2
        IF (OVERRL.GE.2) THEN
          BTYPEF = BTYPE(N)
          IF(VERBOS.EQV..TRUE.)THEN
            WRITE(6,*) ' '  
            WRITE(UNIT=6,FMT=905) 'R-factor gap very similar in the '//
     +    'complex model. Hamilton test overruled!'
          END IF
        END IF
      ELSE IF (BTEST(2).NE.'NULL  ') THEN
C       take the composite score       
        IF (BTEST(2).EQ.BTYPE(N)) THEN
          SCORE = SCORE+1
        END IF
        IF (BTEST(3).EQ.BTYPE(N)) THEN
          SCORE = SCORE+1
        END IF
        IF (BTEST(4).EQ.BTYPE(N)) THEN
          SCORE = SCORE+1
        END IF
        IF (SCORE.GE.2) THEN
          BTYPEF = BTYPE(N)
        ELSE
          BTYPEF = BTYPE(T)
        END IF
      ELSE IF (BTEST(3).NE.'NULL  ') THEN
C       or let the difference test lead
        BTYPEF = BTEST(3)
      ELSE IF (BTEST(4).NE.'BOTH  ') THEN
        BTYPEF = BTEST(4)
      ELSE
C       or take the refinement with the fewest parameters
        BTYPEF = BTYPE(T)
      END IF



C--------------------------- Final steps -------------------------------
899   WRITE(UNIT=6, FMT= 901) BTYPEF  

C-----Skip the error messages
      GO TO 999
C-----Formats
900   FORMAT(A100)
901   FORMAT(A6)
902   FORMAT(A,1X,F5.1,A)
903   FORMAT(A5,1X,I5,2(F5.2,1X),F7.3,2(1X,F9.1))  
904   FORMAT(A,1X,F7.4) 
905   FORMAT(A) 
906   FORMAT(A,1X,I7) 
907   FORMAT(A,2(2X,F5.3))
908   FORMAT(A,1X,A)
909   FORMAT(A,1X,F4.1)

C-----Error messages
996   WRITE(6,*) 'Cannot read the TLS group number'
      GO TO 999
997   WRITE(6,*) 'Cannot read the number of reflections'
      GO TO 999  
998   WRITE(6,*) 'Cannot read the number of atoms'        

C-----End of the line
999   CLOSE (7)
      CLOSE (8)
      END

C-------------------------- Subroutines and functions -----------------C
C-----------------------------------------------------------------------
C  This subroutine reads a table to extract the number of restraints 
C-----------------------------------------------------------------------
      SUBROUTINE RESCNT(FILE,NREST)
      IMPLICIT  NONE
      CHARACTER LINE*100
      INTEGER   FILE, NREST, TREST, I, MAXLIN
      PARAMETER (MAXLIN=2100)

      NREST = 0
      DO 10, I=1, MAXLIN
        READ(UNIT=FILE, FMT=99, END=20) LINE
C-------Stop at the end of the table
        IF(LINE(1:10).EQ.'----------') GO TO 40
C-------Get the number of restraints
        READ(UNIT=LINE(43:51), FMT=*, ERR=21) TREST     
        NREST = NREST + TREST
10    CONTINUE

C-----Error messages
20    WRITE(6,*) 'Unexpected end of file'
      GO TO 30
21    WRITE(6,*) 'Cannot read number of restraints'
30    NREST = -1  
      
40    RETURN

C-----Formats
99    FORMAT(A100)
      END

C-----------------------------------------------------------------------
C  This subroutine reads a table to extract the size of the R-free set 
C-----------------------------------------------------------------------
      SUBROUTINE GETNFR (FILE,NFREE)
      IMPLICIT  NONE
      CHARACTER LINE*100
      INTEGER   FILE, NFREE, TFREE, I, MAXLIN, LENSTR
      PARAMETER (MAXLIN=30)

      NFREE = 0
      DO 10, I=1, MAXLIN
        READ(UNIT=FILE, FMT=99, END=20) LINE
C-------Stop at the end of the table
        IF(LINE(1:2).EQ.'$$') GO TO 40
C-------Skip empty(-ish) line
        IF (LINE(50:57).EQ.'        ') GO TO 10
C-------Get the number of restraints
        READ(UNIT=LINE(50:57), FMT=*, ERR=10) TFREE     
        NFREE = NFREE + TFREE
10    CONTINUE

C-----Error messages
20    WRITE(6,*) 'Unexpected end of file'
      GO TO 30
21    WRITE(6,*) 'Cannot read number of free reflections'
30    NFREE = -1  
      
40    RETURN

C-----Formats
99    FORMAT(A100)
      END

C-----------------------------------------------------------------------
C  This subroutine reads a table to extract the R values
C-----------------------------------------------------------------------
      SUBROUTINE GETRS (FILE, R, RFREE, WR, WRFREE, NREF)
      IMPLICIT  NONE
      CHARACTER LINE*100
      INTEGER   FILE, I, MAXLIN, NREF
      REAL      R, RFREE, WR, WRFREE
      PARAMETER (MAXLIN=50)

C-----Examples
COverall R factor                     =     0.1755
CFree R factor                        =     0.1979
COverall weighted R factor            =     0.1708
CFree weighted R factor               =     0.1880
C12345678901234567890123456789012345678901234567890

C-----Initialise
      R      = -1.0
      RFREE  = -1.0
      WR     = -1.0
      WRFREE = -1.0
      NREF   = 0

      DO 10, I=1, MAXLIN
        READ(UNIT=FILE, FMT=99, END=20) LINE
C-------Stop at the end of the table
        IF(LINE(1:10).EQ.'----------') GO TO 40
C-------Get the values
        IF(LINE(1:16).EQ.'Overall R factor') THEN
          READ(UNIT=LINE(44:49), FMT=*, ERR=21) R     
        ELSE IF (LINE(1:13).EQ.'Free R factor') THEN
          READ(UNIT=LINE(44:49), FMT=*, ERR=21) RFREE
        ELSE IF (LINE(1:26).EQ.'Overall weighted R2 factor') THEN
          READ(UNIT=LINE(44:49), FMT=*, ERR=21) WR
        ELSE IF (LINE(1:23).EQ.'Free weighted R2 factor') THEN
          READ(UNIT=LINE(44:49), FMT=*, ERR=21) WRFREE
        ELSE IF(LINE(1:26).EQ.'Number of used reflections') THEN
          READ(LINE(39:49), FMT=*, ERR=21) NREF
        END IF
10    CONTINUE

C-----Error messages
20    WRITE(6,*) 'Unexpected end of file'
      GO TO 30
21    WRITE(6,*) 'Cannot read value'
30    R = -1.0  
      
40    RETURN

C-----Formats
99    FORMAT(A100)
      END
C-----------------------------------------------------------------------
C  This function returns the expected Rfree/R ratio based on an 
C  empirical model.
C-----------------------------------------------------------------------
      REAL FUNCTION RFRRAT (NREF, NATOM, PARPA)
      IMPLICIT  NONE
      INTEGER NREF, NATOM, PARPA 
C Isotropic  : -0.0050 * reflections/atom + 1.2286
C Anisotropic: -0.0046 * reflections/atom + 1.3375
C Overall    :  1.1700 #Simple flat model due to limited data
      
      IF (PARPA.EQ.3) THEN
        RFRRAT = 1.1700
      ELSE IF (PARPA.EQ.4) THEN 
        RFRRAT = -0.0050 * (1.0*NREF/NATOM) + 1.2286
      ELSE IF (PARPA.EQ.9) THEN
        RFRRAT = -0.0046 * (1.0*NREF/NATOM) + 1.3375
      ELSE
C       Return rediculous value (it will not be used)
        RFRRAT = 10.0
      END IF


      RETURN
      END

C-----------------------------------------------------------------------
C  This function returns the expected Rfree/R ratio based on an 
C  empirical model.
C-----------------------------------------------------------------------
      REAL FUNCTION SRFRRAT (NREF, NATOM, PARPA)
      IMPLICIT  NONE
      INTEGER NREF, NATOM, PARPA 
C Isotropic  : -0.024 * ln(reflections/atom) + 0.1064
C Anisotropic: -0.018 * reflections/atom + 0.0979
C Overall    :  0.0850 #Simple flat model due to limited data
      
      IF (PARPA.EQ.3) THEN
        SRFRRAT = 0.085
      ELSE IF (PARPA.EQ.4) THEN  
        SRFRRAT = -0.024 * LOG(1.0*NREF/NATOM) + 0.1064
      ELSE IF (PARPA.EQ.9) THEN
        SRFRRAT = -0.018 * LOG(1.0*NREF/NATOM) + 0.0979
      ELSE
C       Return rediculous value (it will not be used)
        SRFRRAT = 10.0
      END IF


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
