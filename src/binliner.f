      PROGRAM BINLINER
C=======================================================================
C  Version 2.00 2019-02-28
C  Returns the resolution bins for testing the high resolution cut-off.
C
C  Usage: binliner (-v) CIF_IN (CUT-OFF)(MAX-RES)
C  -v        Verbose output
C  CIF_IN    Is a reflection file in mmCIF format.
C  CUT_OFF   Optional user-defined lower resolution cut-off.
C  MAX_RES   Optional user-defined higher resolution cut-off.
C
C  Written by Robbie P. Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  BINLINER. Please, refer to (one of) these references:
C  - Robbie P. Joosten, Fei Long, Garib N. Murshudov, Anastassis
C    Perrakis: "The PDB_REDO server for macromolecular structure model
C    optimization" IUCrJ, 1(4) p. 213-220 (2014)
C
C Changelog:
C Version 2.00
C - Cell dimensions are now teken directly from the reflection file.
C Version 1.04
C - Bugfix in the way high F/SigF datasates are dealt with.
C - Changed the fall back from 20% to 10% highest resolution relfections
C Version 1.03
C - Fixed format bug for very low resolution cases.
C - Updated reference.
C Version 1.02
C - Updated the outlier rejection so that the resolution gap between two
C   reflections is not greater than 0.2A.
C Version 1.01
C - The target of number of reflections per bin is lowered to 2700. At
C   this number of reflection even a correlation coefficient of 0.05 is
C   still significant at the 1% level.
C Version 1.00
C - Added the optional user-provided maximum  resolution.
C Version 0.04
C - The per-bin F/sigF is now printed in verbose mode.
C Version 0.03
C - The minumum number of bins is now 1.
C - Added some subtle outlier rejection to get rid of crazy high 
C   resolution reflections.
C Version 0.02
C - Added the optional user-provided cut-off.
C Version 0.01
C - First attempt
C=======================================================================
      IMPLICIT NONE
C-----Declare the basic variables and parameters
      INTEGER   I, J, K, STATUS, MAXDAT, ARGS, EXTRA, MAXLAB, MAXCOL,
     +          MAXRES
      CHARACTER VERS*4
      PARAMETER (VERS='1.04')
C-----MAXDAT is the maximum number of reflections
      PARAMETER (MAXDAT=11000000)
C-----MAXLAB is the maximum number column labels
      PARAMETER (MAXLAB=35)
C-----MAXCOL is the number of data columns that is considered for
C-----output: h, k, l, F, sigmaF, I, sigmaI
      PARAMETER (MAXCOL=20)
C-----MAXRES  The maximum number of resolutions in the bins
      PARAMETER (MAXRES=10)

      CHARACTER INCIF*255, LINE*170, CHOP*170, C2JUNK*2, TXTARG*15
C     VERBOS    Give verbose output
C     ISTXT     Currently reading a text block in the cif file
C     GOTSIG    The intensities have sigma values
      LOGICAL   VERBOS, ISTXT,  GOTSIG  
      REAL      CUTOFF, HKLRES, SUMRAT, RATAVE, UCUT, TRESO, HRES
      INTEGER   GETINT, LABCNT, REFLEC, LOWREF, VCOL, SCOL, NRES, NBINS,
     +          MINREF, RSKIP
      DOUBLE PRECISION  GTREAL

C-----Data arrays
C     COLUMN    Assigmnent which mmCIF column contains which data
C               COLUMN(1) = H, COLUMN(2) = K, COLUMN(3) = L
C               COLUMN(4) = F, COLUMN(5) = sigF
C               COLUMN(6) = I, COLUMN(7) = sigI
      CHARACTER LABELS(MAXLAB)*25
      INTEGER   REFH(MAXDAT), REFK(MAXDAT),  REFL(MAXDAT),
     +          COLUMN(MAXCOL)
      REAL      IORF(MAXDAT), SIORF(MAXDAT), REFRES(MAXDAT), 
     +          RATIO(MAXDAT)
C     CELDIM Cell dimensions a, b, c, alpha, beta, gamma
      REAL      CELDIM(6), RESOS(MAXRES), RATS(MAXRES)
C     COEFS coefficients for resolution calculation
      DOUBLE PRECISION COEFS(6)
C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   BINLINER version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Gives bin edges for the high resolution cut-off test.'
      WRITE(6,*)' ' 
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'binliner (flags) A B C ALPHA BETA GAMMA CIF_IN (CUT-'//
     +          'OFF) (MAX-RES)'
      WRITE(6,*)' '
      WRITE(6,*)'CIF_IN    reflection file in mmCIF format.'
      WRITE(6,*)'A, B etc. cell dimensions.'
      WRITE(6,*)'CUT-OFF   optional user-defined low resolution cutoff.'
      WRITE(6,*)'MAX-RES   optional user-defined high resolution '//
     +          'cutoff (use only with CUT-OFF).'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing BINLINER:'
      WRITE(6,*)'Joosten et al. ''The PDB_REDO server for macromolecu'//
     +          'lar structure model optimization'' IUCrJ, 1(4), '//
     +          ', 213-220 (2014)'
        GO TO 999
      END IF
C--------------------------- Main program -----------------------------C

C-----Initialise
      EXTRA    = 0
      REFLEC   = 0
      LOWREF   = 0
      NRES     = 0
      RSKIP    = 0
      SUMRAT   = 0.000
      VERBOS   = .FALSE.
      ISTXT    = .FALSE.
      GOTSIG   = .TRUE.
      CELDIM(1)= -1000.0
      CELDIM(2)= -1000.0
      CELDIM(3)= -1000.0
      CELDIM(4)= -1000.0
      CELDIM(5)= -1000.0
      CELDIM(6)= -1000.0
      UCUT     = 30.00
      HRES     = 0.00

C-----Check for (debug) flags
      DO 10, I=1, 3
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


C-----Get the reflection file
20    CALL GETARG(1+EXTRA, INCIF)
      OPEN(UNIT=7, FILE=INCIF, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*) INCIF, ' cannot be opened!'
        GO TO 999
      END IF
      REWIND(7)

C-----Get the optional resolution cut-off
      IF (ARGS.GT.1+EXTRA) THEN
        CALL GETARG(2+EXTRA, TXTARG)
        UCUT = GTREAL(TXTARG, 1)
      END IF

C-----Get the optional resolution cut-off
      IF (ARGS.GT.2+EXTRA) THEN
        CALL GETARG(3+EXTRA, TXTARG)
        HRES = GTREAL(TXTARG, 1)
      END IF

C-----------------------Work on header---------------------------------C
C-----Read mmCIF header
100   LABCNT = 0
      DO 110, I=1,MAXDAT
        READ(UNIT=7, FMT=951, END=997) LINE
C-------Ignore all text between ;-characters
        IF (ISTXT.EQV..TRUE.) THEN
          LINE = CHOP(LINE,170)
          IF (LINE(1:1).EQ.';') THEN
            ISTXT = .FALSE.
            GO TO 110
          ELSE
            GO TO 110
          END IF
        ELSE
          IF (LINE(1:1).EQ.'#') GO TO 110
          IF (LINE(1:1).EQ.';') THEN
            ISTXT = .TRUE.
            GO TO 110
          END IF
          IF (LINE(1:15).EQ.'_cell.length_a ') THEN
            READ(LINE(16:),*) CELDIM(1)
          END IF     
          IF (LINE(1:15).EQ.'_cell.length_b ') THEN
            READ(LINE(16:),*) CELDIM(2)
          END IF 
          IF (LINE(1:15).EQ.'_cell.length_c ') THEN
            READ(LINE(16:),*) CELDIM(3)
          END IF  
          IF (LINE(1:18).EQ.'_cell.angle_alpha ') THEN
            READ(LINE(19:),*) CELDIM(4)
          END IF 
          IF (LINE(1:17).EQ.'_cell.angle_beta ') THEN
            READ(LINE(18:),*) CELDIM(5)
          END IF      
          IF (LINE(1:18).EQ.'_cell.angle_gamma ') THEN
            READ(LINE(19:),*) CELDIM(6)
          END IF 
          IF (LINE(1:5).EQ.'loop_') GO TO 100
          LINE = CHOP(LINE, 2)
          IF((LINE(1:7).EQ.'_refln.').OR.(LINE(1:7).EQ.'_refln_'))THEN
            LABCNT = LABCNT+1
            IF (LABCNT.GT.MAXLAB) GO TO 996
            LABELS(LABCNT) = LINE(1:25)
          ELSE
            IF (LABCNT .GT. 1) THEN
              IF (VERBOS) WRITE(6,*)'Found ', LABCNT, ' data labels, s',
     +         'electing data columns'
              BACKSPACE(7)
              GO TO 200
            END IF
          END IF
        END IF
110   CONTINUE

C-----Assign column numbers and check to see whether data is complete
C-----Initialise
200   DO 210, I=1, MAXCOL
        COLUMN(I) = 0
210   CONTINUE


C-----Assign
      CALL COLASS(LABELS, LABCNT, COLUMN, MAXCOL, VERBOS)

C-----Is there enough data?
C     Stop if H, K, or L is missing.
      IF ((COLUMN(1).EQ.0).OR.(COLUMN(2).EQ.0).OR.(COLUMN(3).EQ.0))
     +   GO TO 995

C     Use intensities if available
      IF (COLUMN(6).NE.0) THEN
        VCOL   = COLUMN(6)
        SCOL   = COLUMN(7)
        CUTOFF = 2.00
        IF (SCOL.EQ.0) GOTSIG = .FALSE.
      ELSE IF (COLUMN(4).NE.0) THEN
        VCOL   = COLUMN(4)
        SCOL   = COLUMN(5) 
        CUTOFF = 4.00 
      ELSE
        GO TO 994
      ENDIF
      IF (SCOL.EQ.0) GOTSIG = .FALSE. 

      IF (VERBOS.EQV..TRUE.) THEN
        WRITE(UNIT=6, FMT=950) 'Cell dimensions: ', CELDIM
      END IF      
      
C-----Prepare for resolution calculation
      CALL RESCNST(CELDIM, COEFS)  
   
C-------------------Read in reflection data----------------------------C

C-----Read in line by line
300   DO 310, K=1, MAXDAT
        READ(UNIT=7, FMT=951, END=320, ERR=310) LINE

C-------Skip junk lines and data lines with missing or bad reflections
        IF (INDEX(LINE, 'nan').NE.0) GO TO 993
        IF (INDEX(LINE, '*'  ).NE.0) GO TO 993
        IF (INDEX(LINE, '#'  ).NE.0) THEN
          IF (REFLEC.GT.100) THEN
            GO TO 320
          ELSE
            GO TO 994
          END IF
        END IF

C-------Ignore empty lines
        IF (LINE(1:20).EQ.'                    ') GO TO 993

C-------Intercept new dataset problem
        IF (INDEX(LINE, '_').NE.0) GO TO 992

C-------Extract data from lines
        REFLEC = REFLEC+1

C-------HKL indices
        REFH(REFLEC) = GETINT(LINE, COLUMN(1))
C       Reject line if no sensible data is found
        IF ((REFH(REFLEC).EQ.-999999).OR.(REFH(REFLEC).GE.9999)) THEN
          REFLEC=REFLEC-1
          GO TO 993
        END IF

        REFK(REFLEC) = GETINT(LINE, COLUMN(2))
C-----Reject line if no sensible data is found
        IF ((REFK(REFLEC).EQ.-999999).OR.(REFK(REFLEC).GE.9999)) THEN
          REFLEC=REFLEC-1
C          PRINT*, 'K'
          GO TO 993
        END IF

        REFL(REFLEC) = GETINT(LINE, COLUMN(3))
C-----Reject line if no sensible data is found
        IF ((REFL(REFLEC).EQ.-999999).OR.(REFL(REFLEC).GE.9999)) THEN
          REFLEC=REFLEC-1
C          PRINT*, 'L'
          GO TO 993
        END IF

C-----Calculate the resolution
        REFRES(REFLEC) = HKLRES(REFH(REFLEC), REFK(REFLEC),
     +  REFL(REFLEC), COEFS)
C       Reject the reflection if the resolution is too low
        IF (REFRES(REFLEC).LT.HRES) THEN
          REFLEC = REFLEC-1
          RSKIP  = RSKIP+1
          GO TO 310
        END IF

C------Read F or I
        IORF(REFLEC) = GTREAL(LINE, VCOL)
C       Reject line if no sensible data is found
        IF ((IORF(REFLEC).LE.0.000).OR.(IORF(REFLEC).EQ.999999.9)) THEN
          REFLEC=REFLEC-1
C          PRINT*, 'IORF'
          GO TO 993
        END IF
C-------Read the sigma value and calculate the ratio
        IF (GOTSIG.EQV..TRUE.) THEN
          SIORF(REFLEC) = GTREAL(LINE, SCOL) 
C         Reject line if no sensible data is found
          IF ((SIORF(REFLEC).LE.0.000).OR.
     +        (SIORF(REFLEC).EQ.999999.9)) THEN
            REFLEC=REFLEC-1
C            PRINT*, 'Sigma'
            GO TO 993
          END IF
C         Calculate the signal to noise           
          RATIO(REFLEC) = IORF(REFLEC)/SIORF(REFLEC)
        ELSE
          SIORF(REFLEC) = 999
          RATIO(REFLEC) = 999
        END IF

C        PRINT*, REFH(REFLEC), REFK(REFLEC), REFL(REFLEC), 
C     +          IORF(REFLEC), SIORF(REFLEC),REFRES(REFLEC),
C     +          RATIO(REFLEC)

310   CONTINUE

320   IF (VERBOS) WRITE(6,*) REFLEC, 'reflections loaded'
      IF (VERBOS) WRITE(6,*) RSKIP,'high resolution reflections ignored'

C-----Sort the resolution and ratio arrays
      CALL DSORT (REFLEC, REFRES, RATIO)

C-----Throw away outliers in the first 100 reflections     
      TRESO = REFRES(1)
      MINREF = 1
      DO 325, I=1, 100
        IF (REFRES(I)-TRESO.GT.0.20) THEN
          MINREF = I
        END IF
        TRESO = REFRES(I)
325   CONTINUE  

!       DO 330, I=1, MINREF
!           WRITE(6,*) I, REFRES(I), RATIO(I)
! 330   CONTINUE    
!       WRITE(6,*) TRESO

C------------------- Get resolution cut-offs ---------------------------

C-----Low resolution cut-off from user
      IF ((UCUT.GT.0.00).AND.(UCUT.LT.29.00)) THEN
        NRES = 1
        RESOS(1) = UCUT  
        DO 405, I=1, REFLEC
          IF (REFRES(I).GT.UCUT) THEN
            LOWREF = I
            GO TO 425
          ENDIF
405     CONTINUE

C-----Low resolution cut-off from sigma
      ELSE IF (GOTSIG.EQV..TRUE.) THEN
        NRES = 1
C-------Get the average F/sigF for the first 1% of reflections 
        DO 410, I=1, REFLEC/100
          SUMRAT = SUMRAT + RATIO(I)
410     CONTINUE
        RATAVE = SUMRAT/(REFLEC/100)
        NRES = 1
!        PRINT*, RATAVE
        IF (RATAVE.GE.CUTOFF) THEN
C---------Fall back to the 10% highest resoltion reflections        
          RESOS(1) = REFRES(REFLEC/10)  
          LOWREF   = REFLEC/10
        ELSE
C---------Get the moving average stop when it is above the cut-off
          DO 420, I=1, REFLEC-(REFLEC/100)
            IF (RATAVE.GE.CUTOFF) GO TO 425
            SUMRAT = SUMRAT - RATIO(I) + RATIO(I+(REFLEC/100))
            RATAVE = SUMRAT/(REFLEC/100)
            RESOS(1) = REFRES(I+(REFLEC/100))
            LOWREF   = I+(REFLEC/100)
420       CONTINUE    
        END IF
      ELSE  
C-------Take the 10% highest resolution reflections     
        NRES = 1
        RESOS(1) = REFRES(REFLEC/10)  
        LOWREF   = REFLEC/10
      END IF

425   IF (VERBOS) THEN
        WRITE(UNIT=6,FMT=952) 'Binning', LOWREF, 'reflections in reso'//
     +    'lution range', REFRES(1),'-',RESOS(1),'.' 
      END IF   

C-----Decide on the number of bins
      IF (NINT(LOWREF/2700.0).LT.
     +                    NINT((RESOS(1)-REFRES(MINREF))/0.10)) THEN
        NBINS = NINT(LOWREF/2700.0)
      ELSE
        NBINS = NINT((RESOS(1)-REFRES(MINREF))/0.10)
      END IF
      IF (NBINS.GT.(MAXRES-1)) THEN
        NBINS = MAXRES-1  
      END IF  
      IF (NBINS.EQ.0) THEN
        NBINS = 1  
      END IF  

C-----Assign the bins
      DO 430, I=1, NBINS-1 
        NRES       = NRES+1
        RESOS(I+1) = REFRES(LOWREF-I*(LOWREF/NBINS))   
430   CONTINUE
      NRES = NRES+1
      RESOS(NRES) = REFRES(MINREF)

C-----Calculate bin F/sigF
      DO 440, I=1, NRES
        DO 445, J=1+(NRES-I)*(LOWREF/NBINS), (NRES-I+1)*(LOWREF/NBINS)   
          RATS(I) = RATS(I)+RATIO(J)
445     CONTINUE
        RATS(I) = RATS(I)/(LOWREF/NBINS)
440   CONTINUE        

      IF (VERBOS.EQV..TRUE.) THEN
        WRITE(UNIT=6, FMT='(A)') 'F/sigF:'
        WRITE(UNIT=6,FMT=953) (RATS(I),I=1,NRES)
      END IF
 

C-----Write the bins
      WRITE(UNIT=6,FMT=953) (RESOS(I),I=1,NRES)  


C-----Format statements
950   FORMAT(A,3F9.3,3F7.2)
951   FORMAT(A170)
952   FORMAT(A,1X,I7,1X,A,1X,F5.2,A1,F5.2,A)
953   FORMAT(100(F8.2,1X))

C-----Error messages (skip normally)
      GO TO 999
992   WRITE(6,*)'Second dataset detected. This dataset will not be used'
      GO TO 320
993   WRITE(6,*)'Cannot use data, skipping line:'
      WRITE(6,*)LINE(1:60)
      GO TO 300
994   WRITE(6,*)'No structure factors or intensities, aborting!'
      GO TO 999
995   WRITE(6,*)'H, K, L columns not found, aborting!'
      GO TO 999
996   WRITE(6,*)'More labels found than expected. Increase MAXLAB and',
     +' recomplile. Program aborting.'
      GO TO 999
997   WRITE(6,*)'No labels found! Is this realy an mmCIF file?'
      GO TO 999
998   WRITE(6,*) 'Fatal error. Cannot read cell dimensions.'

C-----End of the line
999   CLOSE (7)
      END

C-------------------------- Subroutines and functions -----------------C
C-----------------------------------------------------------------------
C  This subroutine assigns the right column numbers to the data types
C-----------------------------------------------------------------------
      SUBROUTINE COLASS(LABELS, LABCNT, COLUMN, MAXCOL, VERBOS)
      IMPLICIT  NONE
C-----Declare variables and constants
      INTEGER   LABCNT, MAXCOL, COLUMN(MAXCOL), NU(MAXCOL),NUCNT, I
      CHARACTER LABELS(LABCNT)*(*)
      LOGICAL   VERBOS

C-----Initialise
      NUCNT = 0

C-----Loopt through column labels
      IF (VERBOS.EQV..TRUE.) THEN
        WRITE(6,*) ''
        WRITE(6,*) 'Data    Column'
      END IF 

      DO 10, I=1, LABCNT
        IF (LABELS(I)(1:14).EQ.'_refln.index_h') THEN
          COLUMN(1) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'H      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln_index_h') THEN
          COLUMN(1) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'H      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln.index_k') THEN
          COLUMN(2) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'K      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln_index_k') THEN
          COLUMN(2) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'K      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln.index_l') THEN
          COLUMN(3) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'L      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln_index_l') THEN
          COLUMN(3) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'L      ',I
        ELSE IF (LABELS(I)(1:16).EQ.'_refln.F_meas_au') THEN
          COLUMN(4) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'F      ',I
        ELSE IF (LABELS(I)(1:16).EQ.'_refln.F_meas   ') THEN
          COLUMN(4) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'F      ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.F_meas_sigma_au') THEN
          COLUMN(5) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'SIGF   ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.F_meas_sigma   ') THEN
          COLUMN(5) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'SIGF   ',I
        ELSE IF (LABELS(I)(1:21).EQ.'_refln.intensity_meas') THEN
          COLUMN(6) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'I      ',I
        ELSE IF (LABELS(I)(1:21).EQ.'_refln.F_squared_meas') THEN
          COLUMN(6) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'I      ',I
        ELSE IF (LABELS(I)(1:21).EQ.'_refln_F_squared_meas') THEN
          COLUMN(6) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'I      ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.intensity_sigma') THEN
          COLUMN(7) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'SIGI   ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.F_squared_sigma') THEN
          COLUMN(7) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'SIGI   ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln_F_squared_sigma') THEN
          COLUMN(7) = I
          IF (VERBOS) WRITE(UNIT=6,FMT=20) 'SIGI   ',I
        ELSE
          NUCNT = NUCNT+1
          NU(NUCNT)=I
        END IF
10    CONTINUE

C-----Report unused columns:
      IF (VERBOS.EQV..TRUE.) THEN
        IF (NUCNT.GT.0) THEN
          WRITE(6,*) ' '
          WRITE(6,*) 'Unused data columns:'
          DO 15, I=1, NUCNT
            WRITE(6,*) LABELS(NU(I))
15        CONTINUE
          WRITE(6,*) ' '
        END IF
      END IF

20    FORMAT(X,A7,X,I2)
      RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine prepares coefficients for calculating reflection
C  resolution
C  This subroutine was adapted from code in the SFTOOLS program, written
C  by Bart Hazes, which was based on an earlier implementation taken
C  from the Groningen BIOMOL protein crystallography software package.
C-----------------------------------------------------------------------
      SUBROUTINE RESCNST(CELDIM, COEFS)
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
      REAL   CELDIM(6)

C-----Initialise
      DEG2RAD = ATAN(1.0)/45.0

      ca   = COS(DEG2RAD*CELDIM(4))
      sa   = SIN(DEG2RAD*CELDIM(4))
      cb   = COS(DEG2RAD*CELDIM(5))
      sb   = SIN(DEG2RAD*CELDIM(5))
      cg   = COS(DEG2RAD*CELDIM(6))
      sg   = SIN(DEG2RAD*CELDIM(6))
      cast = (cb*cg - ca)/ (sb*sg)
      cbst = (cg*ca - cb)/ (sg*sa)
      cgst = (ca*cb - cg)/ (sa*sb)
      sast = SQRT(1.0 - cast*cast)
      sbst = SQRT(1.0 - cbst*cbst)
      sgst = SQRT(1.0 - cgst*cgst)
      ast  = 1.0/(CELDIM(1)*sb*sgst)
      bst  = 1.0/(CELDIM(2)*sg*sast)
      cst  = 1.0/(CELDIM(3)*sa*sbst)
      COEFS(1)  = ast*ast
      COEFS(2)  = 2.0*ast*bst*cgst
      COEFS(3)  = 2.0*ast*cst*cbst
      COEFS(4)  = bst*bst
      COEFS(5)  = 2.0*bst*cst*cast
      COEFS(6)  = cst*cst

      RETURN
      END

C-----------------------------------------------------------------------
C  This sub routine sorts two arrays 'SORTA' and 'MATCH' based on the
C  values of 'SORTA'
C-----------------------------------------------------------------------
      SUBROUTINE DSORT(N, SORTA, MATCH)
      IMPLICIT NONE      
      INTEGER  N, I, J
      REAL     SORTA(N), MATCH(N), TSORT, TMATCH

C-----Do the sort    
      DO 30 I=2, N
        TSORT  = SORTA(I)
        TMATCH = MATCH(I)
        J = I   
C       This is an implicit loop
10      J = J-1
        IF(J.EQ.0 .OR. SORTA(J).LE.TSORT) GO TO 20
        SORTA(J+1) = SORTA(J) 
        MATCH(J+1) = MATCH(J)
        GO TO 10  
20      SORTA(J+1) = TSORT
        MATCH(J+1) = TMATCH  
30    CONTINUE

      RETURN
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

C-----Convert to Angstroms
      HKLRES = (1.0/SQRT(TMPRES))

      RETURN
      END

C-----------------------------------------------------------------------
C  This function extracts an integer from column CLMNR from a text line
C-----------------------------------------------------------------------
      INTEGER FUNCTION GETINT(LINE, CLMNR)
      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*170 LINE, TLINE, CHOP
      INTEGER   CLMNR, J, SPACE, VALUE

C-----Initialise
      VALUE = -999999
      SPACE = 1

C-----Extract right column
      TLINE = LINE
      DO 10, J=1, CLMNR
        TLINE = TLINE(SPACE:170)
        TLINE = CHOP(TLINE, 170)
        SPACE = INDEX(TLINE, ' ')
        READ(UNIT=TLINE(1:SPACE), FMT=*, END=15, ERR=15) VALUE
        GO TO 10
C-----In case of error use fall-back value
15      VALUE = -999999
10    CONTINUE

C-----Assign and return
99    GETINT = VALUE
      RETURN
      END

C-----------------------------------------------------------------------
C  This function extracts a real number from column CLMNR from a text
C  line.
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GTREAL(LINE, CLMNR)
      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*(*) LINE
      CHARACTER*170 TLINE, CHOP
      CHARACTER*1   QMARK
      INTEGER   CLMNR, J, SPACE
      DOUBLE PRECISION VALUE

C-----Initialise
      SPACE = 1
      VALUE = -999999.9

C-----Extract right column
      TLINE = LINE
      DO 10, J=1, CLMNR
        TLINE = TLINE(SPACE:170)
C        PRINT*, TLINE
        TLINE = CHOP(TLINE, 170)
C        PRINT*, TLINE

C-------Special case when the value is a questionmark
        READ(UNIT=TLINE(1:1), FMT='(A1)') QMARK
        IF (QMARK.EQ.'?') THEN
          VALUE = 999999.9
          SPACE = INDEX(TLINE, ' ')
          IF (SPACE.EQ.0) GO TO 99
          GO TO 10
        END IF

C-------Read normally
        SPACE = INDEX(TLINE, ' ')
        IF (SPACE.EQ.0) GO TO 99
        READ(UNIT=TLINE(1:SPACE), FMT=*, END=15, ERR=15) VALUE
        GO TO 10
C-----In case of error use fall-back value
15      VALUE = -999999.9
10    CONTINUE


C-----Assign and return
99    GTREAL = VALUE
C      WRITE(UNIT=6, FMT='(F14.4)') GTREAL
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

