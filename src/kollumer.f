      PROGRAM KOLLUMER
C======================================================================C
C  Version 2.00 2022-07-11
C  Assigns column labels for mtz2various based on an mtzdmp log file.
C
C  Usage: kollumer (-v) LOG_IN
C
C  LOG_IN is a log file from mtzdmp
C  -v runs the program in verbose mode.
C
C  Written by Robbie P. Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  KOLLUMER. Please, refer to (one of) these references:
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
C  - Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend, 
C    Anastassis Perrakis: "Automatic rebuilding and optimization of 
C    crystallographic structures in the Protein Data Bank" 
C    Bioinformatics, 27, p. 3392-3398 (2011)
C  - Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastassis
C    Perrakis: "PDB_REDO: constructive validation, more than just 
C    looking for errors" Acta Cryst. D68, p. 484-496 (2012)
C
C  Changelog
C  Version 2.00:
C  - If available, phase columns are now written out.
C
C  Version 1.16:
C  - Anomalous intensities are now only written if there are sigma data
C
C  Version 1.15:
C  - Bugfix to better deal with sigma data.

C  Version 1.14:
C  - Better handling for cases without useful data.
C
C  Version 1.13:
C  - FEM map files from Phenix are now recognised.
C
C  Version 1.12:
C  - Fix for anomalous map coefficients.
C
C  Version 1.11:
C  - Fix for figure of merit.
C
C  Version 1.10:
C  - Fix for anomalous map coefficients.
C  - Low completeness anomalous data is ignored. 
C
C  Version 1.09:
C  - Now deals with mtz files without any useful data.
C
C  Version 1.08:
C  - Some standard column names from BUSTER are filtered out.
C
C  Version 1.07:
C  - Some standard column names from PHENIX are filtered out.
C
C  Version 1.06:
C  - Fixed a serious bug in the column assignament logic. Strangely 
C    few assignments were affected.
C
C  Version 1.05:
C  - Filters out anomalous map coefficients from Phenix 
C  - Fixed a bug for handling anomalous amplitude data
C
C  Version 1.04: 
C  - Added support for anomalous amplitude data
C
C  Version 1.03:
C  - Better filtering of 2mFo-DFc columns
C  - Output fix to allow more columns
C
C  Version 1.02: 
C  - Added support for anomalous intensity data
C
C  Version 1.01:
C  - Bugfix that avoids selecting the FC column
C 
C  Version 1.00:
C  - Added a simple help function
C
C  Version 0.2:
C  - Bugfix.
C
C  Version 0.1:
C  - First attempt
C
C======================================================================C
      IMPLICIT NONE
C-----Declare basic variables and parameters
      INTEGER   MAXDAT, MAXLAB, I, J, K, STATUS, ARGS
      CHARACTER VERS*4
      PARAMETER (VERS='2.00J')
C-----MAXDAT is the maximum allowed lines in the log file
      PARAMETER (MAXDAT=2000)
C-----MAXLAB is the maximum allowed different labels
      PARAMETER (MAXLAB=100)

C-----Declare the variables and parameters
      CHARACTER LINE*150, INLOG*255, OUTLAB*500, C2JUNK*2
      INTEGER   EXTRA, NLABEL, LENSTR, LABLEN, BESTFJ, BESTQ, BESTI, 
     +          BESTA, BESTP,  
     +          IDF, IDJ, IDKP, IDKM, IDGP, IDGM, IDI, ID, IDMP, IDMM,
     +          IDLP, IDLM, IDA, IDP
C-----IDF is the number of the most likely amplitude column label
      LOGICAL   VERBOS, GTDATA

C-----Declare arrays
      CHARACTER LABEL(MAXLAB)*30, ULABEL(MAXLAB)*30, CTYPE(MAXLAB)*1
      REAL      COMPL(MAXLAB)   , RESHI(MAXLAB), RESLO(MAXLAB)
      INTEGER   NTYPE(10)  
C-----NTYPE is the number of labels of a certain type 
C           in the order FJQIKMGLPA

C========================End of declarations===========================C
C
C-----Start with the help function:
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   KOLLUMER version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Kollumer uses a log file from mtzdmp to assign the '//
     +          'correct reflection data columns for PDB_REDO optimiz'//
     +          'ation of a structure model.'
      WRITE(6,*)'Written by Robbie Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'kollumer (flags) LOG_IN'
      WRITE(6,*)' '
      WRITE(6,*)'LOG_IN is a file containing the output of mtzdmp.'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-v Verbose mode. Shows the data columns that may con'//
     +          'tain useful data.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing KOLLUMER:'
      WRITE(6,*)'Joosten et al. ''The PDB_REDO server for macromolecu'//
     +          'lar structure model optimization'' IUCrJ, 1(4), '//
     +          ', 213-220 (2014)'      
        GO TO 999
      END IF



C========================    Get the input   ==========================C
C-----Initialise values
      VERBOS  = .FALSE.
      GTDATA  = .FALSE.
      EXTRA   = 0
      NLABEL  = 0
      NTYPE(1)= 0
      NTYPE(2)= 0
      NTYPE(3)= 0
      NTYPE(4)= 0
      NTYPE(5)= 0
      NTYPE(6)= 0
      NTYPE(7)= 0
      NTYPE(8)= 0
      NTYPE(9)= 0
      NTYPE(10)= 0
      IDF     = 0
      IDJ     = 0
      IDI     = 0
      IDKP    = 0
      IDMP    = 0
      IDKM    = 0
      IDMM    = 0
      OUTLAB  = ''
      LABLEN  = LENSTR(OUTLAB)

C-----Check for (debug) flags
      DO 1, I=1, 2
        CALL GETARG(1+EXTRA, C2JUNK)
C-------Are there any flags?
        IF (C2JUNK(1:1).NE.'-') GO TO 10
C       Is it a valid flag 
        IF ((C2JUNK.EQ.'-v').OR.(C2JUNK.EQ.'-V')) THEN
          VERBOS = .TRUE.
          EXTRA  = EXTRA+1
        ELSE
          EXTRA = EXTRA+1
        END IF
1     CONTINUE  

C-----Get log file
10    STATUS=0
C     First log file
      CALL GETARG(1+EXTRA, INLOG)
      OPEN(UNIT=7, FILE=INLOG, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) 'Cannot open ', INLOG 
          GO TO 999
        END IF
      REWIND (7)

C========================  Read the log file ==========================C

C-----Read log file line by line
      DO 100, I=1, MAXDAT 
        READ(UNIT=7, FMT=980, END=150) LINE
        IF (LINE(1:9).EQ.' Col Sort') THEN 
C         Skip a line
          READ(UNIT=7, FMT=980, END=150) LINE  
C         Read the Table
          DO 110, J=1, MAXDAT
            READ(UNIT=7, FMT=980, END=150) LINE
            IF (LINE(1:8).EQ.'        ') THEN
              GO TO 110
            ELSE IF (LINE(1:8).EQ.' No. of ') THEN
C             End of the table
              GO TO 150   
            ELSE IF (INDEX('FJQIKMGLAP',LINE(74:74)).NE.0) THEN
C             This is a supported data column type 
              NLABEL = NLABEL+1
              READ(UNIT=LINE, FMT=979,ERR=105,END=110) COMPL(NLABEL),
     +        RESLO(NLABEL), RESHI(NLABEL),CTYPE(NLABEL),LABEL(NLABEL)
              GO TO 110
105           READ(UNIT=LINE, FMT=981,ERR=998,END=110) COMPL(NLABEL),
     +        RESLO(NLABEL), RESHI(NLABEL),CTYPE(NLABEL),LABEL(NLABEL) 
            END IF
110       CONTINUE         
        END IF 
100   CONTINUE

C-----Also store uppercase labels
150   DO 155, I=1, NLABEL
        ULABEL(I) = LABEL(I)
        CALL UPPER (ULABEL(I))
155   CONTINUE

C-----Report
      IF(VERBOS.EQV..TRUE.) THEN
        WRITE(6,*)'Useful columns:'
        DO 160, I=1, NLABEL
          WRITE(UNIT=6, FMT=978) LABEL(I), CTYPE(I), RESHI(I), 
     +      RESLO(I), COMPL(I)
160     CONTINUE
      END IF

C========================= Assign columns =============================C
C-----Count the types
      DO 210, I=1, NLABEL
        IF (CTYPE(I).EQ.'F') THEN
          NTYPE(1) = NTYPE(1) + 1
        ELSE IF(CTYPE(I).EQ.'J') THEN
          NTYPE(2) = NTYPE(2) + 1
        ELSE IF(CTYPE(I).EQ.'Q') THEN
          NTYPE(3) = NTYPE(3) + 1    
        ELSE IF(CTYPE(I).EQ.'I') THEN
          NTYPE(4) = NTYPE(4) + 1
        ELSE IF(CTYPE(I).EQ.'K') THEN
          NTYPE(5) = NTYPE(5) + 1
        ELSE IF(CTYPE(I).EQ.'M') THEN
          NTYPE(6) = NTYPE(6) + 1
        ELSE IF(CTYPE(I).EQ.'G') THEN
          NTYPE(7) = NTYPE(7) + 1
        ELSE IF(CTYPE(I).EQ.'L') THEN
          NTYPE(8) = NTYPE(8) + 1
        ELSE IF(CTYPE(I).EQ.'P') THEN
          NTYPE(9) = NTYPE(9) + 1   
        ELSE IF(CTYPE(I).EQ.'A') THEN
          NTYPE(10)= NTYPE(10) + 1             
        ELSE
          GO TO 997       
        END IF
210   CONTINUE
      IF (VERBOS.EQV..TRUE.) THEN
        WRITE(6,*)'Intensity columns          :', NTYPE(2)
        WRITE(6,*)'Amplitude columns          :', NTYPE(1)
        WRITE(6,*)'Anomalous intensity columns:', NTYPE(5)
        WRITE(6,*)'Anomalous amplitude columns:', NTYPE(7)
        WRITE(6,*)'Sigma columns              :', NTYPE(3)
        WRITE(6,*)'Anomalous sigma I columns  :', NTYPE(6)
        WRITE(6,*)'Anomalous sigma F columns  :', NTYPE(8)
        WRITE(6,*)'R-free columns             :', NTYPE(4)
        WRITE(6,*)'Phase columns (phi)        :', NTYPE(9)
        WRITE(6,*)'Phase columns (HL coefs)   :', NTYPE(10)      
      END IF

C-----Check to see whether there is any useful data
      IF (NTYPE(2) + NTYPE(1) + NTYPE(5) + NTYPE(7) .LT. 1) GO TO 995

C-----Amplitude columns
      IF (NTYPE(1).GT.0) THEN 
C       Select the best amplitude
        IDF = BESTFJ(ULABEL, CTYPE, RESHI, RESLO, COMPL, NLABEL,'F ')
C       Exit if there are no alternative data columns
        IF (IDF.NE.0) THEN
          GTDATA = .TRUE.
          OUTLAB = OUTLAB(1:LABLEN)//'FP='//LABEL(IDF)
          LABLEN = LENSTR(OUTLAB)
        ENDIF  
      END IF

C-----Intensity columns
      IF (NTYPE(2).GT.0) THEN 
C       Select the best intensity
        IDJ = BESTFJ(ULABEL, CTYPE, RESHI, RESLO, COMPL, NLABEL,'J ')
        IF (IDJ.NE.0) THEN
          GTDATA = .TRUE.
          OUTLAB = OUTLAB(1:LABLEN)//' I='//LABEL(IDJ)
          LABLEN = LENSTR(OUTLAB)
        END IF  
      END IF

C-----Anomalous intensity columns
      IF (NTYPE(5).GT.1) THEN
C       Select the best anomalous intensity columns
        IDKP = BESTFJ(ULABEL, CTYPE, RESHI, RESLO, COMPL, NLABEL,'KP')
        IF (IDKP.EQ.0) GO TO 996
        IF (IDKP.EQ.-1) THEN 
C         Unusable anomalous data        
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*)'Anomalous intensity data too incomplete!'
          END IF
          GO TO 214
        ENDIF
        
        IDKM = BESTFJ(ULABEL, CTYPE, RESHI, RESLO, COMPL, NLABEL,'KM')
        IF (IDKM.EQ.0) GO TO 996
        IF (IDKM.EQ.-1) THEN 
C         Unusable anomalous data        
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*)'Anomalous intensity data too incomplete!'
          END IF
          GO TO 214
        ENDIF
C       Write labels later
      END IF

C-----Anomalous amplitude columns
214   IF (NTYPE(7).GT.1) THEN
C       Select the best anomalous amplitude columns
        IDGP = BESTFJ(ULABEL, CTYPE, RESHI, RESLO, COMPL, NLABEL,'GP')
        IF (IDGP.EQ.0) GO TO 996
        IF (IDGP.EQ.-1) THEN 
C         Unusable anomalous data        
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*)'Anomalous amplitude data too incomplete!'
          END IF
          GO TO 215
        ENDIF
        
        IDGM = BESTFJ(ULABEL, CTYPE, RESHI, RESLO, COMPL, NLABEL,'GM')
        IF (IDGM.EQ.0) GO TO 996
        IF (IDGM.EQ.-1) THEN 
C         Unusable anomalous data        
          IF (VERBOS.EQV..TRUE.) THEN
            WRITE(6,*)'Anomalous amplitude data too incomplete!'
          END IF
          GO TO 215
        ENDIF
C       Write labels
        GTDATA = .TRUE.
        OUTLAB = OUTLAB(1:LABLEN)//' F(+)='//LABEL(IDGP)
        LABLEN = LENSTR(OUTLAB)
        OUTLAB = OUTLAB(1:LABLEN)//' F(-)='//LABEL(IDGM)
        LABLEN = LENSTR(OUTLAB)
      END IF

C-----Sigma columns
215   IF ((NTYPE(3).GT.0).AND.(IDF.NE.0)) THEN 
C       Select the best amplitude sigma 
        ID = BESTQ(ULABEL, CTYPE, IDF, NLABEL)
        IF (ID.EQ.0) GO TO 220
        OUTLAB = OUTLAB(1:LABLEN)//' SIGFP='//LABEL(ID)
        LABLEN = LENSTR(OUTLAB)
      END IF
220   IF ((NTYPE(3).GT.0).AND.(IDJ.NE.0)) THEN 
C       Select the best intensity sigma 
        ID = BESTQ(ULABEL, CTYPE, IDJ, NLABEL)
        IF (ID.EQ.0) GO TO 230
        OUTLAB = OUTLAB(1:LABLEN)//' SIGI='//LABEL(ID)
        LABLEN = LENSTR(OUTLAB)
      END IF

C-----Anomalous sigma intensity columns
      IF ((NTYPE(6).GT.1).AND.(IDKP.GT.0).AND.(IDKM.GT.0)) THEN
C       Select the best anomalous intensity sigma columns
        IDMP = BESTQ(ULABEL, CTYPE, IDKP, NLABEL)
        IF (IDMP.EQ.0) GO TO 225
        IDMM = BESTQ(ULABEL, CTYPE, IDKM, NLABEL)
        IF (IDMM.EQ.0) GO TO 225
        
C       We now have anomalous data and their sigmas 
        GTDATA = .TRUE.
        OUTLAB = OUTLAB(1:LABLEN)//' I(+)='//LABEL(IDKP)
        LABLEN = LENSTR(OUTLAB)
        OUTLAB = OUTLAB(1:LABLEN)//' I(-)='//LABEL(IDKM)
        LABLEN = LENSTR(OUTLAB)
        OUTLAB = OUTLAB(1:LABLEN)//' SIGI(+)='//LABEL(IDMP)
        LABLEN = LENSTR(OUTLAB)
        OUTLAB = OUTLAB(1:LABLEN)//' SIGI(-)='//LABEL(IDMM)
        LABLEN = LENSTR(OUTLAB)
      END IF

C-----Anomalous sigma amplitude columns
225   IF ((NTYPE(8).GT.1).AND.(IDGP.GT.0).AND.(IDGM.GT.0)) THEN
C       Select the best anomalous intensity sigma columns
        IDLP = BESTQ(ULABEL, CTYPE, IDGP, NLABEL)
        IF (IDLP.EQ.0) GO TO 230
        OUTLAB = OUTLAB(1:LABLEN)//' SIGF(+)='//LABEL(IDLP)
        LABLEN = LENSTR(OUTLAB)

        IDLM = BESTQ(ULABEL, CTYPE, IDGM, NLABEL)
        IF (IDLM.EQ.0) GO TO 230
        OUTLAB = OUTLAB(1:LABLEN)//' SIGF(-)='//LABEL(IDLM)
        LABLEN = LENSTR(OUTLAB)
      END IF

C-----Free columns
230   IF (NTYPE(4).GT.0) THEN 
C       Select the best R-free column
        IDI = BESTI(ULABEL, CTYPE, COMPL, NLABEL)
        IF (IDI.EQ.0) GO TO 235
        OUTLAB = OUTLAB(1:LABLEN)//' FREE='//LABEL(IDI)
        LABLEN = LENSTR(OUTLAB)
      END IF
      
C-----Phase columns (prefer Hendrickson-Lattman format
235   IF (NTYPE(10).GE.4) THEN
C       Get the right HL columns
        IDA = BESTA(ULABEL, CTYPE, COMPL, NLABEL,'A')
        OUTLAB = OUTLAB(1:LABLEN)//' HLA='//LABEL(IDA)
        LABLEN = LENSTR(OUTLAB)
        IDA = BESTA(ULABEL, CTYPE, COMPL, NLABEL,'B')
        OUTLAB = OUTLAB(1:LABLEN)//' HLB='//LABEL(IDA)
        LABLEN = LENSTR(OUTLAB)
        IDA = BESTA(ULABEL, CTYPE, COMPL, NLABEL,'C')
        OUTLAB = OUTLAB(1:LABLEN)//' HLC='//LABEL(IDA)
        LABLEN = LENSTR(OUTLAB)
        IDA = BESTA(ULABEL, CTYPE, COMPL, NLABEL,'D')
        OUTLAB = OUTLAB(1:LABLEN)//' HLD='//LABEL(IDA)
        LABLEN = LENSTR(OUTLAB)
      ELSE IF (NTYPE(9).GE.1) THEN
C       Use the regular phase instead   
        IDP = BESTP(ULABEL, CTYPE, COMPL, NLABEL)
        IF (IDP.GT.0) THEN
          OUTLAB = OUTLAB(1:LABLEN)//' PHIB='//LABEL(IDP)
          LABLEN = LENSTR(OUTLAB)
        ENDIF  
      ENDIF

C-----Write final label assignment
      IF (GTDATA) THEN
        WRITE(UNIT=6,FMT=977)OUTLAB(1:LABLEN)
      ELSE
        GO TO 995
      ENDIF

C-----Skip error messages
      GO TO 999
     
      
C===========================Formats====================================C
977   FORMAT(A)
978   FORMAT(A30,X,A1,3(X,F6.2))
979   FORMAT(32X,F6.2,19X,F6.2,X,F6.2,3X,A1,2X,A30)
980   FORMAT(A150)
981   FORMAT(34X,F6.2,19X,F6.2,X,F6.2,1X,A1,2X,A30)
 
C===========================Error messages=============================C
995   WRITE(6,*) 'Problem: no useful data found. Cannot continue!'
      GO TO 999
996   WRITE(6,*) 'Problem extracting the columns. Cannot continue!'
      GO TO 999
997   WRITE(6,*) 'Unrecognised data type ', CTYPE(I),' Cannot continue!'
      GO TO 999
998   WRITE(6,*) 'Problem reading line: ', LINE, ' Cannot continue!'

C===========================End of the line============================C
999   CLOSE (7)
      END
  
C======================================================================C
C   SUBROUTINES AND FUNCTIONS                                          C
C======================================================================C

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
C  This function returns the number of the most likely amplitude or 
C  intensity column. Other return values are: 
C  0 --> no suitable column found
C -1 --> found column has low completeness
C-----------------------------------------------------------------------
      INTEGER FUNCTION BESTFJ(ULABEL,CTYPE,RESHI,RESLO,COMPL,NLABEL,TYP)
      IMPLICIT NONE
      INTEGER   I,NLABEL, BCOL
      REAL      BHRES, BLRES, BCOMPL
      CHARACTER TYP*2
      CHARACTER ULABEL(NLABEL)*30, CTYPE(NLABEL)*1
      REAL      COMPL(NLABEL)   , RESHI(NLABEL), RESLO(NLABEL)
C---------------------------------------------------------------
 
C-----Initialise
      BESTFJ = 0
      BCOL   = 0
      BHRES  = 999.99
      BLRES  = 0.00
      BCOMPL = 0.00

      DO 10 I=1,NLABEL
C       Make sure the column is of the right type
C       PRINT*,ULABEL(I)
        IF (CTYPE(I).EQ.TYP(1:1)) THEN
C         Do anomalous data first
C       PRINT*,TYP, ULABEL(I)
          IF ((TYP(2:2).EQ.'P').AND.((INDEX(ULABEL(I),'(+)').NE.0).OR.
     +        (INDEX(ULABEL(I),'PLUS').NE.0))) THEN
            IF (COMPL(I) .LT. 10.0) THEN
              BCOL = -1
            ELSE
              BCOL = I
            END IF  
            GO TO 20
          ELSE IF((TYP(2:2).EQ.'M').AND.
     +            ((INDEX(ULABEL(I),'(-)').NE.0).OR.
     +             (INDEX(ULABEL(I),'MINUS').NE.0))) THEN
            IF (COMPL(I) .LT. 10.0) THEN
              BCOL = -1
            ELSE
              BCOL = I
            END IF  
            GO TO 20     
C         Look for a 'native' flag
          ELSE IF ((TYP(2:2).EQ.' ').AND.
     +            (INDEX(ULABEL(I),'NAT').NE.0)) THEN
            BCOL = I
            GO TO 20
          ELSE IF ((ULABEL(I)(1:5).EQ.'F    ').OR.
     +             (ULABEL(I)(1:5).EQ.'FOBS ').OR.
     +             (ULABEL(I)(1:5).EQ.'FP   ')) THEN
            BCOL = I
            GO TO 20
          ELSE IF ((ULABEL(I)(1:3).EQ.'FC_').OR.
     +             (ULABEL(I)(1:3).EQ.'FC ').OR.
     +             (ULABEL(I)(1:7).EQ.'F-MODEL')) THEN
C           This is an Fcalc column, skip
            GO TO 10
          ELSE IF ((ULABEL(I)(1:3).EQ.'FWT').OR.
     +             (ULABEL(I)(1:6).EQ.'DELFWT').OR.
     +             (ULABEL(I)(1:3).EQ.'FAN').OR.
     +             (ULABEL(I)(1:6).EQ.'DELFAN')) THEN
C           This is a map coefficient column from Refmac, skip
            GO TO 10       
          ELSE IF ((ULABEL(I)(1:7).EQ.'2FOFCWT').OR.
     +             (ULABEL(I)(1:6).EQ.'FOFCWT' ).OR.
     +             (ULABEL(I)(1:4).EQ.'FCTR' )) THEN
C           This is a map coefficient column from BUSTER, skip
            GO TO 10  
          ELSE IF ((ULABEL(I)(1:4).EQ.'ANOM').OR.
     +             (ULABEL(I)(1:3).EQ.'FEM').OR.
     +             (ULABEL(I)(1:8).EQ.'2MFO-DFC')) THEN
C           These are map coefficients column from Phenix, skip
            GO TO 10 
          ELSE IF (ULABEL(I)(1:3).EQ.'FDM') THEN
C           This is density modified map coefficient column, skip
            GO TO 10
          ELSE IF (ULABEL(I)(1:3).EQ.'FOM') THEN
C           This is the figure of merit, skip
            GO TO 10
          ELSE IF (RESHI(I).LT.BHRES) THEN
            BHRES  = RESHI(I)
            BLRES  = RESLO(I)
            BCOMPL = COMPL(I)
            BCOL   = I
C            PRINT*, "RESH"
          ELSE IF ((RESHI(I).EQ.BHRES).AND.(RESLO(I).GT.BLRES)) THEN
            BHRES  = RESHI(I)
            BLRES  = RESLO(I)
            BCOMPL = COMPL(I)
            BCOL   = I
C            PRINT*, "RESHL"
          ELSE
            IF ((RESHI(I).EQ.BHRES).AND.(RESLO(I).EQ.BLRES).AND.
     + (COMPL(I).GT.BCOMPL)) THEN
              BHRES  = RESHI(I)
              BLRES  = RESLO(I)
              BCOMPL = COMPL(I)
              BCOL   = I  
C              PRINT*, "RESHL C"
            END IF
          END IF
        END IF
C        PRINT*, BCOL, BHRES, BLRES, BCOMPL, ULABEL(BCOL) 
   10 CONTINUE
   20 BESTFJ = BCOL
      RETURN
      END

C-----------------------------------------------------------------------
C  This function returns the number of the most likely sigma column for 
C  a given amplitude or intensity column
C-----------------------------------------------------------------------
      INTEGER FUNCTION BESTQ(ULABEL,CTYPE,ID,NLABEL)
      IMPLICIT NONE
      INTEGER   I,NLABEL, BCOL, ID, EIND, LENSTR, SIG
      CHARACTER ULABEL(NLABEL)*30, CTYPE(NLABEL)*1, MODLABEL*30 
C---------------------------------------------------------------
 
C-----Initialise
      BESTQ = 0
      BCOL  = 0
      EIND  = LENSTR(ULABEL(ID))

      DO 10 I=1,NLABEL
C        PRINT*, ULABEL(I),ULABEL(ID)(1:EIND)
C       Make sure the column is of the right type
        IF ((CTYPE(I).EQ.'Q').OR.(CTYPE(I).EQ.'M').OR.
     +      (CTYPE(I).EQ.'L')) THEN
C         Look for a clearly matching name
          IF (INDEX(ULABEL(I),ULABEL(ID)(1:EIND)).NE.0) THEN
            BCOL = I
            GO TO 20
          ELSE
C           Remove SIG form the lable and see if it matches
            SIG = INDEX(ULABEL(I),'SIG')
            IF (SIG .EQ. 0) GOTO 10
            MODLABEL =  ULABEL(I)(1:SIG -1)//ULABEL(I)(SIG +3:)
            IF (MODLABEL.EQ.ULABEL(ID)) THEN
              BCOL = I
            END IF
          END IF
        END IF
   10 CONTINUE
   20 BESTQ = BCOL
      RETURN
      END

C-----------------------------------------------------------------------
C  This function returns the number of the most likely R-free column
C-----------------------------------------------------------------------
      INTEGER FUNCTION BESTI(ULABEL,CTYPE,COMPL,NLABEL)
      IMPLICIT NONE
      INTEGER   I,NLABEL, BCOL
      REAL      BHRES, BLRES, BCOMPL
      CHARACTER TYP*1
      CHARACTER ULABEL(NLABEL)*30, CTYPE(NLABEL)*1
      REAL      COMPL(NLABEL) 
C---------------------------------------------------------------
 
      BCOMPL = 0.00
C-----Initialise
      BESTI  = 0
      BCOL   = 0

      DO 10 I=1,NLABEL
C       Make sure the column is of the right type
        IF (CTYPE(I).EQ.'I') THEN
C         Look for a 'free' flag
          IF (INDEX(ULABEL(I),'FREE').NE.0) THEN
            BCOL = I
            GO TO 20 
          ELSE 
            IF (COMPL(I).GT.BCOMPL) THEN
              BCOMPL = COMPL(I)
              BCOL   = I 
            END IF 
          END IF
        END IF
   10 CONTINUE
   20 BESTI = BCOL
      RETURN
      END

C-----------------------------------------------------------------------
C  This function returns the number of the most likely HL coef column
C-----------------------------------------------------------------------
      INTEGER FUNCTION BESTA(ULABEL,CTYPE,COMPL,NLABEL,COEF)
      IMPLICIT NONE
      INTEGER   I,NLABEL, BCOL
      REAL      BHRES, BLRES, BCOMPL
      CHARACTER COEF*1
      CHARACTER ULABEL(NLABEL)*30, CTYPE(NLABEL)*1
      REAL      COMPL(NLABEL) 
C---------------------------------------------------------------
 
      BCOMPL = 0.00
C-----Initialise
      BESTA  = 0
      BCOL   = 0

      DO 10 I=1,NLABEL
C       Make sure the column is of the right type
        IF (CTYPE(I).EQ.'A') THEN
C         Look for a 'free' flag
          IF (INDEX(ULABEL(I),'HL'//COEF).NE.0) THEN
            BCOL = I
            GO TO 20 
          ELSE 
            IF (COMPL(I).GT.BCOMPL) THEN
              BCOMPL = COMPL(I)
              BCOL   = I 
            END IF 
          END IF
        END IF
   10 CONTINUE
   20 BESTA = BCOL
      RETURN
      END 
      
C-----------------------------------------------------------------------
C  This function returns the number of the most likely phase column
C-----------------------------------------------------------------------
      INTEGER FUNCTION BESTP(ULABEL,CTYPE,COMPL,NLABEL)
      IMPLICIT NONE
      INTEGER   I,NLABEL, BCOL
      REAL      BHRES, BLRES, BCOMPL
      CHARACTER ULABEL(NLABEL)*30, CTYPE(NLABEL)*1
      REAL      COMPL(NLABEL) 
C---------------------------------------------------------------
 
      BCOMPL = 0.00
C-----Initialise
      BESTP  = 0
      BCOL   = 0

      DO 10 I=1,NLABEL
C       Make sure the column is of the right type
        IF (CTYPE(I).EQ.'P') THEN
C         Look for a 'free' flag
          IF (INDEX(ULABEL(I),'PHIF').NE.0) THEN
            BCOL = I
            GO TO 20 
          ELSE IF ((ULABEL(I)(1:4).EQ.'PHIC').OR.
     +             (ULABEL(I)(1:4).EQ.'PHWT').OR.
     +             (ULABEL(I)(1:7).EQ.'PHDELWT')) THEN
C           This is a phase calc column, skip  
            GO TO 10
          ELSE 
            IF (COMPL(I).GT.BCOMPL) THEN
              BCOMPL = COMPL(I)
              BCOL   = I 
            END IF 
          END IF
        END IF
   10 CONTINUE
   20 BESTP = BCOL
      RETURN
      END       
      
C-----------------------------------------------------------------------
C  This subroutine converts strings to uppercase
C-----------------------------------------------------------------------
      SUBROUTINE UPPER(STRING)
      IMPLICIT NONE
C     JOHN MAHAFFY  3/8/96
C
      CHARACTER STRING*(*)
      INTEGER   LC, ICDIFF, I

      LC=LEN(STRING)
C
      ICDIFF=ICHAR('A')-ICHAR('a')
      DO 10 I=1,LC
C       Skip the uppercase letters
        IF(STRING(I:I).LT.'a'.OR.STRING(I:I).GT.'z') GO TO 10
C       Conver the lowercase letters
        STRING(I:I)=CHAR(ICHAR(STRING(I:I))+ICDIFF)
10    CONTINUE
      RETURN
      END




            
