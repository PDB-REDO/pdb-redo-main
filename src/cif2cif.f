      PROGRAM CIF2CIF
C=======================================================================
C  Version 9.02 2023-08-01
C  Cleans up mmCIF files to ensure that they can be used in an automated
C  fashion. It works on all reasonably valid mmCIF reflection files. One
C  weakpoint is that it cannot always handle files in which the data
C  columns run together.
C
C  Usage: cif2cif (flags) FILENAME_IN FILENAME_OUT (wavelength)
C  FILENAME_IN is an mmCIF file straight from the PDB, FILENAME_OUT is
C  the output mmCIF file. Keep in mind that not all data colums are
C  kept. If the 'status' column exists, it is transferred to the output
C  file unless the 'ignore' flag is used when starting cif2cif.
C  The -i option forces the use of intensities even when amplitudes are
C  present in the input file.
C  The -c option forces the conversion of intensities to amplitudes
C  using F=sqrt(I) and sigF=sigI/(2F)
C
C  CIF2CIF can intercept the following problems:
C  1. Lack of sigma values (defaults are used)
C  2. Sigma values all set to 0.0 (defaults, 0.01, are used to avoid
C     refmac problems)
C  3. Sigma values that are all the same (no information content)
C  4. Individual reflections with sigma equal to 0.0 These reflections
C     are not really infinitely reliable. Therefore, sigma is set to
C     the highest encountered sigma value + 0.01*
C  5. Certain wrong label names
C  6. Values that are not properly space delimited (to a certain extent)
C  7. Alternative use of the refln.status flag (normal = 1, R-free = -1)
C  8. More than one dataset in a single file (only first set is used)
C  9. Status flag used for reflection bins. 0 is assumed to be the
C     R-free set (CCP4 default). The status flags are rejected if more
C     than 100 bins are used.
C 10. Status flag column without any information content (no output
C     file, use the '-i' flag to create an mmCIF file without a status 
C     column)
C 11. Non-positive values for _refln.F_meas_au*
C 12. Incomplete and strange end-of-file reflections
C 13. R-free sets that are too small to be used
C
C  *Not checked when only intensities and no amplitudes are encountered
C   or when the '-i' flag is active. Use the CCP4 program CTRUNCATE or
C   similar to convert the intensities to amplitudes and perform
C   corrections.
C
C  Written by Robbie P. Joosten
C  E-mail:    r.joosten@nki.nl, robbie_joosten@hotmail.com
C
C  If you publish results (directly or indirectly) obtained by using
C  CIF2CIF. Please, refer to (one of) these references:
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
C  Changelog
C  Version 9.02:
C  - Using double precision for cell dimensions.
C  Version 9.01:
C  - Added contingency for one of the anomalous I columns containing 
C    no information.
C  Version 9.00:
C  - Now writing a fake scale_group_code value.
C  Version 8.00:
C  - Cell dimensions and space group are retained.
C  - Added fallback for missing values.
C  - Added fix for rare space group names. 
C  Version 7.08:
C  - Fix for dealing with multiple datasets.
C  Version 7.07:
C  - If there is no FOM column when experimental phases are given,
C    default valies of 0.71 are added.
C  Version 7.05:
C  - Anomalous intensities are also converted when the -c flag is given.
C  Version 7.04:
C  - Rediculously high-index reflections are deleted.
C  Version 7.03:
C  - Finished the anomalous mode.
C  - Fixed writing HL coefficients for missing data.
C  - If there is no anomalous difference, the columns are not written.
C  Version 7.02:
C  - Made the line buffers longer. Replaced the CHOP function to be 
C    length independent. Removed the now obsolete CHOP2 function.
C  Version 7.01:
C  - The number of data per line is now tested for each reflection, not
C    just the first one. This is slower, but more reliable. The DATCHK 
C    function is updated to be faster.
C  Version 7.00:
C  - Experimental phase information is now also properly parsed.
C
C  Version 6.01:
C  - Changed cut-off for number of reflections with poor F/sigF to 45%.
C  Version 6.00:
C  - Anomalous data columns are now written out in anomalous mode.

C  Version 5.13:
C  - Added additional sanity check for F/sigF or I/sigI.
C  Version 5.12:
C  - Added a sanity check for F/sigF or I/sigI.
C  Version 5.11:
C  - Added the '-d' switch to ignore sigF and sigI values.
C  Version 5.10:
C  - Added support for another mmCIF label describing the wavelength.
C  Version 5.09:
C  - If given, the (first) wavelength is now read from the mmCIF file 
C    and copied to the output.
C  - A wavelength to be added to the output can be passed to cif2cif.
C  - Dropped support for the 'ignore' keyword. Use '-s' flag instead. 
C  - Fixed a bug that did not allow merging anomalous intensities and 
C    converting to amplitudes at the same time.
C  Version 5.08:
C  - Fixed a bug in the merging of reflection. This process is now moved
C    to separate fuctions. 
C  - Amplitudes are converted to intensies, merged, and converted back.
C  Version 5.07:
C  - Now using the variance instead of the standard deviation to merge
C    anomalous reflections.
C  - Test sets smaller than 1% of the reflections are now flagged as too
C    small because they upset the program FREERFLAG.
C  Version 5.06:
C  - Fixed a bug for dealing with question marks.
C  Version 5.05:
C  - Fixed the detection and removal of tab characters.
C  Version 5.04:
C  - Increased the maximum nuber of reflections to 11M.
C  Version 5.03:
C  - Added the option to not do an R-free sanity check.
C  Version 5.02:
C  - Added a basic help function
C  Version 5.1:
C  - Amplitudes, intensities, and their sigma values are now read with
C    DOUBLE PRECISION.
C  Version 5.0:
C  - Now also supports anomalous data. It is converted to mean data:
C    Imean = (I+/SIGI+  +  I-/SIGI-) / (1/SIGI+  +  1/SIGI-)
C  - Small format fix for very large reflection file.
C  - Any R-free set that has fewer than 500 reflections and is smaller 
C    than 10% of the total number or reflections is rejected.
C
C  Version 4.9:
C  - Now also supports multi-line text blocks that are not properly 
C    closed.
C  - Any R-free set that has fewer than 500 reflections and is smaller 
C    than 5% of the total number or reflections is rejected.
C  Version 4.8:
C  - Included support for more data labels.
C   
C=======================================================================
      IMPLICIT NONE
C-----Declare the variables and parameters
      INTEGER   MAXDAT, MAXLAB, MAXCOL, I, J, K, L, STATUS
      CHARACTER FILLER*4
      CHARACTER VERS*4
      PARAMETER (VERS='9.02')
C-----MAXDAT is the maximum number of reflections in the file. This
C-----should be enough for almost all reflection files.
      PARAMETER (MAXDAT=11000000)
C-----MAXLAB is the maximum number column labels
      PARAMETER (MAXLAB=35)
C-----MAXCOL is the number of data columns that is considered for
C-----output: h, k, l, F, sigmaF, I, sigmaI, status, wavelength_id,
C             F+, sigF+, F-, sigF-, I+, sigI+, I-, sigI-, wavelength,
C             phase, FOM, HLA, HLB, HLC, HLD (in this order)
      PARAMETER (MAXCOL=24) 
C-----Filler is just some text used to comply with the mmCIF standard
      PARAMETER (FILLER='1 1 ')
      CHARACTER INCIF*255,  OUTCIF*255, GTFREE*1,  WAVEIN*10,  C2JUNK*2,
     +          SYMMHM*13, SYMMID*13, CELLID*13
      CHARACTER LABELS(MAXLAB)*39,  RFREE(MAXDAT)*1
      CHARACTER*250    LINE,   LINE1, DEJUNK, CHOP
      INTEGER   LABCNT, HASH,   EMPTY,  GOTMIN, REFLEC, GETINT, NAN,
     +          MINONE, UNDRSC, ARGS,   TABCHA, STAR,   DATCHK, NDAT1,
     +          NDAT2,  NDAT3,  EXTRA,  CCP4RF, FRCNT,  ISAD,   IGNCNT,
     +          LOWS2N, POS, LENSTR
C  ISAD Status of anomalous data. 0 = not anomalous, 1 = anomalous F,
C       2 = Anomalous I
      INTEGER   REFH(MAXDAT),   REFK(MAXDAT),   REFL(MAXDAT), 
     +          WAVLID(MAXDAT), COLUMN(MAXCOL)
      DOUBLE PRECISION SF(MAXDAT), SIGMAF(MAXDAT), INTENS(MAXDAT),
     +          SIGMAI(MAXDAT),    FPLUS(MAXDAT),  FMIN(MAXDAT),
     +          SFPLUS(MAXDAT),    SFMIN(MAXDAT),  IPLUS(MAXDAT),  
     +          SIPLUS(MAXDAT),    IMIN(MAXDAT),   SIMIN(MAXDAT),
     +          HL(4,MAXDAT),      PHASE(MAXDAT),  FOM(MAXDAT)
      DOUBLE PRECISION GTREAL, APLUS, AMIN, SAPLUS, SAMIN,T1, T2, T3,T4,
     +                 WPLUS, WMIN, IMERGE, SIMERGE, RATSUM, SIGNAL,
     +                 SDANO,
     +                 AAXIS,  BAXIS,  CAXIS,  ALPHA,  BETA,   GAMMA   
      REAL      TMPVAL, FRFRAC, BLAAT,  BL2,    WAVEL,  LOWS2NFRAC
      LOGICAL   IGNORE, USEI,   DEFSIG, USESI,  HVFREE, FIXNEG, ISCCP4,
     +          FRINFO, ISTWIN, ISTRIP, FORCEI, FORI2F, ISTXT, FRIDLM, 
     +          USEWVL, SINGL,  INSANE, WAVELC, USES, USEA, USEHL,DEFFOM
C  SINGL is true when only one of two anomalous reflections is recorded
C  USEA  is true, then only anomalous data is written out
C  USEHL is true when all four Hendrickson Latman coefficients are found
C=======================================================================
C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   CIF2CIF version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Cif2cif standardises reflection file for PDB_REDO an'//
     +          'd performs a few sanity checks on the data.' 
      WRITE(6,*)'Written by Robbie P. Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'cif2cif (flags) CIF_IN CIF_OUT'
      WRITE(6,*)' '
      WRITE(6,*)'CIF_IN is the original (mmCIF format) reflection file.' 
      WRITE(6,*)'CIF_OUT is the standardised reflection file.'
      WRITE(6,*)'    '
      WRITE(6,*)'Possible flags:'
      WRITE(6,*)'-i Force intensities. Use reflection intensities ins'//
     +          'tead of amplitudes if possible.'
      WRITE(6,*)'-c Convert intensities. Intensities are converted to'//
     +          ' amplitudes using F=sqrt(I) and sigF=sigI/(2F).'
      WRITE(6,*)'-s Ignore status. The reflection status column is ig'//
     +          'nored altogether.'
      WRITE(6,*)'-n Not sane modus. There will be no sanity checks on'//
     +          ' the R-free set.'
      WRITE(6,*)'-g No sigma values. The sigI and sigF columns will b'//
     +          'e ignored.'
      WRITE(6,*)'-a Anomalous modus. Write out anomalous data columns.'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing CIF2CIF:'
      WRITE(6,*)'Joosten et al. ''Re-refinement from deposited X-ray '//
     +          'data can deliver improved models for most PDB entrie'//
     +          's'' Acta Cryst., D65, 176-185 (2009)'
        GO TO 999
      END IF

C=========================== Main program =============================C
C-----Initialise
      EXTRA  = 0
      STATUS = 0
      IGNORE = .FALSE.
      FORCEI = .FALSE.
      FORI2F = .FALSE.
      ISTXT  = .FALSE.
      INSANE = .FALSE.
      WAVELC = .FALSE.
      USEA   = .FALSE.
      FRIDLM = .FALSE.
      USES   = .TRUE.
      WAVEL  = -10.0000
      RATSUM = 0.0
      IGNCNT = 0
      LOWS2N = 1
      SYMMID = 'redo'
      CELLID = 'redo'
      WRITE(6,*) 'Welcome to CIF2CIF version ', VERS
      WRITE(6,*) ' '
      
      
C-----Check (debug) flags
      DO 10, I=1, 6
        CALL GETARG(1+EXTRA, C2JUNK)
C-------Are there any flags?
        IF (C2JUNK(1:1).NE.'-') GO TO 20

        IF ((C2JUNK.EQ.'-i').OR.(C2JUNK.EQ.'-I')) THEN
          WRITE(6,*) 'Using intensities when possible'
          FORCEI= .TRUE.
          EXTRA = EXTRA+1
        ELSE IF ((C2JUNK.EQ.'-c').OR.(C2JUNK.EQ.'-C')) THEN
          WRITE(6,*) 'Converting intensities to amplitudes'   
          FORI2F= .TRUE.
          EXTRA = EXTRA+1 
        ELSE IF ((C2JUNK.EQ.'-s').OR.(C2JUNK.EQ.'-S')) THEN
          WRITE(6,*) 'Status column will be ignored'   
          IGNORE= .TRUE.
          EXTRA = EXTRA+1 
        ELSE IF ((C2JUNK.EQ.'-n').OR.(C2JUNK.EQ.'-N')) THEN
          WRITE(6,*) 'There will be no checks of R-free set sanity'   
          INSANE= .TRUE.
          EXTRA = EXTRA+1  
        ELSE IF ((C2JUNK.EQ.'-g').OR.(C2JUNK.EQ.'-G')) THEN
          WRITE(6,*) 'Reflection sigma values will be ignored'   
          USES  = .FALSE.
          EXTRA = EXTRA+1  
        ELSE IF ((C2JUNK.EQ.'-a').OR.(C2JUNK.EQ.'-A')) THEN
          WRITE(6,*) 'Writing out anomalous data (only)'   
          USEA  = .TRUE.
          EXTRA = EXTRA+1      
        ELSE
          WRITE(6,*) 'Invalid option: ', C2JUNK
          WRITE(6,*) 'It will be ignored'
	  EXTRA = EXTRA+1
        END IF
10    CONTINUE  
   
C-----Get input file
20    CALL GETARG(1+EXTRA, INCIF)
      OPEN(UNIT=7, FILE=INCIF, IOSTAT=STATUS)
        IF(STATUS .NE. 0) THEN
          WRITE(6,*) INCIF, ' cannot be opened!'
          GO TO 999
        END IF
      REWIND (7)
      
C-----Read in a wavelength if needed     
      IF (ARGS.GE.3+EXTRA) THEN
        CALL GETARG(3+EXTRA, WAVEIN)
        READ(UNIT=WAVEIN, FMT=*, ERR=100) WAVEL
        IF (WAVEL.GT.0) THEN
          WAVELC = .TRUE.
          WRITE(6,*) 'Using wavelength', WAVEL,'A' 
        ELSE
          WAVEL = -10
	END IF
      END IF

      
C-----------------------Work on header---------------------------------C


C-----Read mmCIF header
C     Initialise
100   LABCNT=0
      DO 105, I=1, MAXCOL
        COLUMN(I)=0
105   CONTINUE

C     Read the header
      DO 110, I=1,MAXDAT
        READ(UNIT=7, FMT=902, END=998) LINE
        LINE = CHOP(LINE, 2)
C-------Ignore all text between ;-characters
        IF (ISTXT.EQV..TRUE.) THEN
          LINE = CHOP(LINE,250)
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
C         Read wavelength without loop
          IF (((LINE(1:35).EQ.'_diffrn_radiation_wavelength.wavele').OR.
     +         (LINE(1:33).EQ.'_diffrn_radiation.pdbx_wavelength')).AND.
     +         (WAVELC.EQV..FALSE.)) THEN
            CALL SECFIX(LINE)
            CALL QUOFIX(LINE)
            WAVEL = GTREAL(LINE,2)
            IF (WAVEL.GT.0) THEN
              WRITE(6,*) 'Using wavelength', WAVEL,'A' 
            END IF
            GO TO 110
          END IF
C         Check for data loop
          IF (LINE(1:5).EQ.'loop_') THEN
C           Read labels in the loop
            DO 120, J=1,MAXLAB+10
              READ(UNIT=7, FMT=902, END=998) LINE
              LINE = CHOP(LINE, 2)
              IF (LINE(1:1).EQ.'#') GO TO 120
              IF ((LINE(1:7).EQ.'_refln.').OR.(LINE(1:7).EQ.'_refln_')
     +          .OR.(LINE(1:17).EQ.'_diffrn_radiation')) THEN
                LABCNT=LABCNT+1
	        IF (LABCNT.GT.MAXLAB) GO TO 993
	        LABELS(LABCNT)=LINE(1:35)
	      ELSE 
                IF (LABCNT .GT. 1) THEN
C                 We are in a usable data loop
                  WRITE(6,*)'Found ', LABCNT, ' data labels, selecting',
     +            ' data columns'

C                 Assign data columns
                  CALL COLASS(LABELS, LABCNT, COLUMN, MAXCOL)

C                 Check in which loop we are
                  IF (COLUMN(18).NE.0) THEN
C                   Read wavelength (only the first one)
                    CALL SECFIX(LINE)
                    WAVEL = GTREAL(LINE,COLUMN(18))
                    IF (WAVEL.GT.0) THEN
                      WRITE(6,*) 'Using wavelength', WAVEL,'A' 
                    END IF
                    GO TO 100
                  ELSE
C                   We are in the reflection data loop  
                    BACKSPACE(7)
                    GO TO 200
                  END IF
                ELSE 
C                 We are in some other loop
                  GO TO 110                 
                END IF
              END IF
120         CONTINUE              
          END IF
C         Read the spacegroup
          IF (LINE(1:18).EQ.'_symmetry.entry_id') THEN
            READ(LINE(19:),*) SYMMID
          END IF
          IF (LINE(1:30).EQ.'_symmetry.space_group_name_H-M') THEN
            READ(LINE(31:),*) SYMMHM
C           Correct the space group
            IF (SYMMHM.EQ.'P 21 21 2 A ') THEN
              SYMMHM = 'P 21 21 2 (a)'
            END IF  
          END IF   
C         Read the cell dimensions
          IF (LINE(1:14).EQ.'_cell.entry_id') THEN
            READ(LINE(15:),*) CELLID
          END IF
          IF (LINE(1:15).EQ.'_cell.length_a ') THEN
            READ(LINE(16:),*) AAXIS
            PRINT*, AAXIS
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
        END IF
110   CONTINUE


C-----Check to see whether diffraction data is complete
200   USEI   = .FALSE.
      DEFSIG = .FALSE.
      DEFFOM = .FALSE.
      USESI  = .FALSE.
      HVFREE = .FALSE.
      ISAD   = 0
      USEWVL = .FALSE.
      USEHL  = .FALSE.
      
C-----Is there enough data?
C     Stop if H, K, or L is missing.
      IF ((COLUMN(1).EQ.0).OR.(COLUMN(2).EQ.0).OR.(COLUMN(3).EQ.0))
     +   GO TO 997

C     Stop if there are no (anomalous) I's or F's.
      IF ((COLUMN(4).EQ.0).AND.(COLUMN(6).EQ.0).AND.(COLUMN(10).EQ.0)
     +    .AND.(COLUMN(12).EQ.0).AND.(COLUMN(14).EQ.0)
     +    .AND.(COLUMN(16).EQ.0)) GO TO 996

C     Reset sigma value columns if needed.
      IF (USES.EQV..FALSE.) THEN
        COLUMN(5)  = 0
        COLUMN(7)  = 0
        COLUMN(11) = 0
        COLUMN(13) = 0
        COLUMN(15) = 0
        COLUMN(17) = 0
      END IF

C     Warn if we have no sigma values.
      IF ((COLUMN(5).EQ.0).AND.(COLUMN(7).EQ.0).AND.(COLUMN(11).EQ.0)
     +    .AND.(COLUMN(13).EQ.0).AND.(COLUMN(15).EQ.0)
     +    .AND.(COLUMN(17).EQ.0)) GO TO 995 
      

C-----Look for status flag (if it should not be ignored)
298   IF ((IGNORE.EQV..FALSE.) .AND. (COLUMN(8).NE.0)) THEN
        HVFREE=.TRUE.
      END IF

C-----Look for wavelength id
299   IF (COLUMN(9).NE.0) THEN
        USEWVL=.TRUE.
      END IF
      
C-----Look for HL coefficients
      IF ((COLUMN(21).NE.0).AND.(COLUMN(22).NE.0).AND.
     +    (COLUMN(23).NE.0).AND.(COLUMN(24).NE.0)) THEN
        USEHL=.TRUE.
      END IF  

C-----Use default FOM if there are pahses, but no FOM or HLs     
      IF ((COLUMN(19).NE.0).AND.(COLUMN(20).EQ.0).AND.
     +    (USEHL.EQV..FALSE.)) THEN
        DEFFOM = .TRUE.
        WRITE(6,*)'No phase error estimates. Default FOMs used.'
      ENDIF

C-----Decide which data will be used.
      IF ((COLUMN(4).NE.0).AND.(COLUMN(5).NE.0).AND.
     +    (FORCEI.EQV..FALSE.)) THEN
C       Use F and sigF
      ELSE IF ((COLUMN(6).NE.0).AND.(COLUMN(7).NE.0))THEN
C       Use I and sigI
        USEI  = .TRUE.
        USESI = .TRUE.    
      ELSE IF ((COLUMN(10).NE.0).AND.(COLUMN(12).NE.0).AND.
     +         (COLUMN(11).NE.0).AND.(COLUMN(13).NE.0).AND.
     +         (FORCEI.EQV..FALSE.))THEN
C       Use anomalous amplitudes
        ISAD   = 1
        FRIDLM = .TRUE.
      ELSE IF ((COLUMN(14).NE.0).AND.(COLUMN(16).NE.0).AND.
     +         (COLUMN(15).NE.0).AND.(COLUMN(17).NE.0))THEN
C       Use anomalous intensities
        ISAD  = 2
        USEI  = .TRUE.
        USESI = .TRUE. 
        FRIDLM = .TRUE.
      ELSE IF (COLUMN(4).NE.0) THEN
C       Use F without sigF
        DEFSIG = .TRUE.
      ELSE IF (COLUMN(6).NE.0) THEN
C       Use I and sigI
        USEI  = .TRUE. 
        DEFSIG= .TRUE.
      ELSE IF ((COLUMN(10).NE.0).AND.(COLUMN(12).NE.0))THEN
C       Use anomalous amplitudes
        ISAD   = 1
        DEFSIG = .TRUE.
        FRIDLM = .TRUE.
      ELSE IF ((COLUMN(14).NE.0).AND.(COLUMN(16).NE.0))THEN
C       Use anomalous intensities
        ISAD = 2
        USEI  = .TRUE.
        DEFSIG= .TRUE. 
        FRIDLM = .TRUE.
      ELSE
C       Unexpected situation
        GO TO 988      
      END IF        

C-----Write out anomalous data
      IF (USEA.EQV..TRUE.) THEN
        IF ((COLUMN(10).NE.0).AND.(COLUMN(12).NE.0).AND.
     +   (FORCEI.EQV..FALSE.)) THEN
          ISAD = 1
        ELSE IF ((COLUMN(14).NE.0).AND.(COLUMN(16).NE.0)) THEN
          ISAD = 2
        ELSE
          USEA = .FALSE.
        END IF
      END IF  
        
C-----Finish header work
      WRITE(6,*)'All required data columns found.'
      WRITE(6,*)' '
C-------------------------Done parsing header--------------------------C

C-------------------Read in reflection data----------------------------C
C-----Initialise
399   REFLEC=0

C-----Read in line by line
400   DO 410, K=1, MAXDAT
        READ(UNIT=7, FMT=902, END=499, ERR=410) LINE
C       First check for comments
        IF (LINE(1:1).EQ.'#') THEN
          IF (REFLEC.GT.100) THEN
            GO TO 499
          ELSE
            GO TO 994
          END IF
        END IF       
C       Check for new loops or data sets
        IF ((LINE(1:5).EQ.'data_').OR.(LINE(1:5).EQ.'loop_')) GO TO 991
        
C       Check data items on the line 
        DO 401, L=1, 3
          NDAT1 = DATCHK(LINE)
          IF (NDAT1.EQ.LABCNT) THEN
C           Exactly enough data, continue
            GO TO 411
          ELSE IF (NDAT1.GT.LABCNT) THEN
C           This is a corrupt line 
            WRITE(6,*) LINE
            GO TO 989
          ELSE 
C           Not enough data, read another line        
	    READ(UNIT=7, FMT=902, END=499, ERR=410) LINE1
	    LINE = DEJUNK(LINE, LINE1)
	  END IF
401     CONTINUE
C-----Debug 
        WRITE(6,*) LINE

C-----Skip junk lines and data lines with missing or bad reflections
411     NAN = INDEX(LINE, 'nan')
        IF (NAN.NE.0) GO TO 994
        STAR = INDEX(LINE, '*')
        IF (STAR.NE.0) GO TO 994

C-----Fix 'tab-problem'
        TABCHA = INDEX(LINE, CHAR(9))
        IF (TABCHA.NE.0) THEN
          CALL TABFIX(LINE)
        END IF

C-------Ignore empty lines
        IF (LINE(1:20).EQ.'                    ') GO TO 994

C-----Fix 'minus-problem'
        GOTMIN = INDEX(LINE, '-')
        IF (GOTMIN.NE.0) THEN
          CALL MINFIX(LINE)
        END IF

C-----Extract data from lines
        REFLEC=REFLEC+1

C-------Wavelength ID
        IF (USEWVL.EQV..TRUE.) THEN
          WAVLID(REFLEC) = GETINT(LINE, COLUMN(9))
C---------Use '1' as fallback value
          IF (WAVLID(REFLEC).EQ.-999999) THEN
	    WAVLID(REFLEC) = 1
	  END IF
        ELSE
C---------Set default
          WAVLID(REFLEC) = 1
        END IF
	
C-------HKL indices
        REFH(REFLEC) = GETINT(LINE, COLUMN(1))
C-----Reject line if no sensible data is found
        IF ((REFH(REFLEC).EQ.-999999).OR.(REFH(REFLEC).GE.320)) THEN
          REFLEC=REFLEC-1
          GO TO 994
        END IF

        REFK(REFLEC) = GETINT(LINE, COLUMN(2))
C-----Reject line if no sensible data is found
	IF ((REFK(REFLEC).EQ.-999999).OR.(REFK(REFLEC).GE.320)) THEN
	  REFLEC=REFLEC-1
	  GO TO 994
	END IF

        REFL(REFLEC) = GETINT(LINE, COLUMN(3))
C-----Reject line if no sensible data is found
        IF ((REFL(REFLEC).EQ.-999999).OR.(REFL(REFLEC).GE.320)) THEN
	  REFLEC=REFLEC-1
	  GO TO 994
	END IF

C-----Reject F(0,0,0) because it crashes ctruncate
        IF ((REFH(REFLEC).EQ.0).AND.(REFK(REFLEC).EQ.0).AND.
     +      (REFL(REFLEC).EQ.0)) THEN
          REFLEC=REFLEC-1
          GO TO 994   
	END IF

C-------Intensities and structure factors
        SINGL = .FALSE.
        IF (USEI.EQV..FALSE.) THEN
C         Read amplitudes
          IF (ISAD.EQ.1) THEN
C           Read F+ and F-
            APLUS = GTREAL(LINE, COLUMN(10))
            AMIN  = GTREAL(LINE, COLUMN(12))
C           Deal with missing values
	    IF (((APLUS.EQ.-999999.9).OR.(APLUS.EQ.999999.9)).AND.
     +          ((AMIN .EQ.-999999.9).OR.(AMIN.EQ.999999.9))) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
            ELSE IF ((APLUS.EQ.-999999.9).OR.(APLUS.EQ.999999.9)) THEN
              APLUS = AMIN
              SINGL = .TRUE.
            ELSE IF ((AMIN .EQ.-999999.9).OR.(AMIN .EQ.999999.9)) THEN
              AMIN  = APLUS
              SINGL = .TRUE.
            END IF
C           Read sigF+ and sigF-
            IF (DEFSIG.EQV..TRUE.) THEN
              SAPLUS = 1.0
              SAMIN  = 1.0
            ELSE
              SAPLUS = GTREAL(LINE, COLUMN(11))
              SAMIN  = GTREAL(LINE, COLUMN(13))
            END IF
C           Deal with missing values
	    IF (((SAPLUS.LT.0.0).OR.(SAPLUS.EQ.999999.9)).AND.
     +          ((SAMIN .LT.0.0).OR.(SAMIN.EQ.999999.9))) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
            ELSE IF ((SAPLUS.LT.0.0).OR.(SAPLUS.EQ.999999.9)) THEN
              SAPLUS = SAMIN
              APLUS  = AMIN
              SINGL  = .TRUE.
            ELSE IF ((SAMIN .LT.0.0).OR.(SAMIN.EQ.999999.9)) THEN
              SAMIN  = SAPLUS
              AMIN   = APLUS
              SINGL  = .TRUE.
            END IF
C           Read normal amplitudes or merge the Friedel pair          
            IF (FRIDLM.EQV..TRUE.)THEN
C             Calculate the mean values
              IF (SINGL.EQV..FALSE.) THEN
                T1 = APLUS*APLUS
                T2 = AMIN*AMIN
                T3 = SAPLUS*2*APLUS
                T4 = SAMIN*2*AMIN
                SF(REFLEC) = SQRT(IMERGE(T1, T2, T3, T4))
                SIGMAF(REFLEC) = SIMERGE(T1, T2, T3, T4)/(2* SF(REFLEC))
              ELSE
                SF(REFLEC) = APLUS
                SIGMAF(REFLEC) = SAPLUS
              END IF
            ELSE  
C             Just read F and sigF
              SF(REFLEC) = GTREAL(LINE, COLUMN(4))
C             Reject line if no sensible data is found
	      IF ((SF(REFLEC).EQ.-999999.9).OR.
     +            (SF(REFLEC).EQ.999999.9)) THEN
	        REFLEC=REFLEC-1
	        GO TO 994
	      END IF
              IF (DEFSIG.EQV..TRUE.) THEN
                SIGMAF(REFLEC) = 0.0000
              ELSE
                SIGMAF(REFLEC) = GTREAL(LINE, COLUMN(5))
C               Reject line if no sensible data is found
                IF (SIGMAF(REFLEC).EQ.-999999.9) THEN
	          REFLEC=REFLEC-1
	          GO TO 994
                END IF
C               Set value to 0.0000 for '?' (see GTREAL)
	        IF (SIGMAF(REFLEC).EQ.999999.9) THEN
	          SIGMAF(REFLEC) = 0.0000
                END IF
	      END IF
	    END IF  
	      
C           Reject negative SF
            IF (SF(REFLEC).LE.0.0000) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
            END IF
            
C           Populate the anomalous columns if needed
            IF (USEA.EQV..TRUE.) THEN
              FPLUS(REFLEC)  = APLUS
              FMIN(REFLEC)   = AMIN
              SFPLUS(REFLEC) = SAPLUS
              SFMIN(REFLEC)  = SAMIN
            END IF
            
          ELSE
C           Read F
            SF(REFLEC) = GTREAL(LINE, COLUMN(4))
C            WRITE(6,*) LINE, SF(REFLEC)
C           Reject line if no sensible data is found
	    IF ((SF(REFLEC).EQ.-999999.9).OR.
     +          (SF(REFLEC).EQ.999999.9)) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
	    END IF
C           Reject negative SF
            IF (SF(REFLEC).LE.0.0000) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
	    END IF
C           Read sigF
            IF (DEFSIG.EQV..TRUE.) THEN
              SIGMAF(REFLEC) = 0.0000
            ELSE
              SIGMAF(REFLEC) = GTREAL(LINE, COLUMN(5))
C             Reject line if no sensible data is found
              IF (SIGMAF(REFLEC).EQ.-999999.9) THEN
	        REFLEC=REFLEC-1
	        GO TO 994
              END IF
C             Set value to 0.0000 for '?' (see GTREAL)
	      IF (SIGMAF(REFLEC).EQ.999999.9) THEN
	        SIGMAF(REFLEC) = 0.0000
              END IF
	    END IF
          END IF
        ELSE
C         Read intensities
          IF (ISAD.EQ.2) THEN
C           Read I+ and I-
            APLUS = GTREAL(LINE, COLUMN(14))
            AMIN  = GTREAL(LINE, COLUMN(16))
C           Deal with missing values
	    IF (((APLUS.EQ.-999999.9).OR.(APLUS.EQ.999999.9)).AND.
     +          ((AMIN .EQ.-999999.9).OR.(AMIN.EQ.999999.9))) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
            ELSE IF ((APLUS.EQ.-999999.9).OR.(APLUS.EQ.999999.9)) THEN
              APLUS = AMIN
              SINGL = .TRUE.
            ELSE IF ((AMIN .EQ.-999999.9).OR.(AMIN .EQ.999999.9)) THEN
              AMIN  = APLUS
              SINGL = .TRUE.
            END IF  
C           Read sigI+ and sigI-
            IF (DEFSIG.EQV..TRUE.) THEN
              SAPLUS = 1.0
              SAMIN  = 1.0
            ELSE
              SAPLUS = GTREAL(LINE, COLUMN(15))
              SAMIN  = GTREAL(LINE, COLUMN(17))
            END IF
C           Deal with missing values
	    IF (((SAPLUS.LT.0.0).OR.(SAPLUS.EQ.999999.9)).AND.
     +          ((SAMIN .LT.0.0).OR.(SAMIN.EQ.999999.9))) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
            ELSE IF ((SAPLUS.LT.0.0).OR.(SAPLUS.EQ.999999.9)) THEN
              SAPLUS = SAMIN
              APLUS  = AMIN
              SINGL  = .TRUE.
            ELSE IF ((SAMIN .LT.0.0).OR.(SAMIN.EQ.999999.9)) THEN
              SAMIN  = SAPLUS
              AMIN   = APLUS
              SINGL  = .TRUE.
            END IF
            
C           Read normal intensities or merge the Friedel pair          
            IF (FRIDLM.EQV..TRUE.)THEN
C             Calculate the mean values
              IF (SINGL.EQV..FALSE.) THEN
                T1 = APLUS*APLUS
                T2 = AMIN*AMIN
                T3 = SAPLUS*2*APLUS
                T4 = SAMIN*2*AMIN
                INTENS(REFLEC) = IMERGE(APLUS, AMIN, SAPLUS, SAMIN)
                SIGMAI(REFLEC) = SIMERGE(APLUS, AMIN, SAPLUS, SAMIN)
              ELSE
                INTENS(REFLEC) = APLUS
                SIGMAI(REFLEC) = SAPLUS
              END IF                    
            ELSE
C             Read intensity
              INTENS(REFLEC) = GTREAL(LINE, COLUMN(6))
C             Reject line if no sensible data is found
	      IF ((INTENS(REFLEC).EQ.-999999.9).OR.
     +          (INTENS(REFLEC).EQ.999999.9)) THEN
	        REFLEC=REFLEC-1
	        GO TO 994
	      END IF
C             Read sigI
              IF (DEFSIG.EQV..TRUE.) THEN
                SIGMAI(REFLEC) = 0.0100
              ELSE
                SIGMAI(REFLEC) = GTREAL(LINE, COLUMN(7))
	        IF (SIGMAI(REFLEC).EQ.-999999.9) THEN
	          REFLEC=REFLEC-1
	          GO TO 994
                END IF
C               Set value to 0.0000 for '?' (see GTREAL)
	        IF (SIGMAI(REFLEC).EQ.999999.9) THEN
	          SIGMAI(REFLEC) = 0.0000             
                END IF
              END IF
            END IF  
              
C           Populate the anomalous columns if needed
            IF (USEA.EQV..TRUE.) THEN
              IPLUS(REFLEC)  = APLUS
              IMIN(REFLEC)   = AMIN
              SIPLUS(REFLEC) = SAPLUS
              SIMIN(REFLEC)  = SAMIN              
            END IF 
          ELSE  
C           Just read intensities
            INTENS(REFLEC) = GTREAL(LINE, COLUMN(6))
C           Reject line if no sensible data is found
	    IF ((INTENS(REFLEC).EQ.-999999.9).OR.
     +          (INTENS(REFLEC).EQ.999999.9)) THEN
	      REFLEC=REFLEC-1
	      GO TO 994
	    END IF
C           Read sigI
            IF (DEFSIG.EQV..TRUE.) THEN
              SIGMAI(REFLEC) = 0.0100
            ELSE
              SIGMAI(REFLEC) = GTREAL(LINE, COLUMN(7))
              IF (SIGMAI(REFLEC).EQ.-999999.9) THEN
                REFLEC=REFLEC-1
                GO TO 994
              END IF
C             Set value to 0.0000 for '?' (see GTREAL)
              IF (SIGMAI(REFLEC).EQ.999999.9) THEN
               SIGMAI(REFLEC) = 0.0000             
              END IF
            END IF
              
          END IF 
C         Convert I to F
          IF (FORI2F.EQV..TRUE.) THEN
            IF (INTENS(REFLEC).LE.0.00) THEN
              REFLEC=REFLEC-1
	      GO TO 994
            END IF
            SF(REFLEC) = SQRT(INTENS(REFLEC))
            SIGMAF(REFLEC) = SIGMAI(REFLEC)/(2*SF(REFLEC))  
            IF (USEA.EQV..TRUE.) THEN
              IF ((IPLUS(REFLEC).LT.0.0).OR.(IMIN(REFLEC).LT.0.0)) THEN
                FPLUS(REFLEC) = SF(REFLEC)
                FMIN(REFLEC)  = SF(REFLEC)
                SFPLUS(REFLEC)= SIGMAF(REFLEC)
                SFMIN(REFLEC) = SIGMAF(REFLEC)
              ELSE
                FPLUS(REFLEC)  = SQRT(IPLUS(REFLEC))
                FMIN(REFLEC)   = SQRT(IMIN(REFLEC))
                SFPLUS(REFLEC) = SIPLUS(REFLEC)/(2*FPLUS(REFLEC))
                SFMIN(REFLEC)  = SIMIN(REFLEC)/(2*FMIN(REFLEC))
              ENDIF  
C              PRINT*, REFH(REFLEC), REFK(REFLEC), REFL(REFLEC),
C     +       FPLUS(REFLEC),FMIN(REFLEC),SFPLUS(REFLEC),SFMIN(REFLEC)
            END IF  
          END IF       
        END IF
    
C-------Phases
        IF (COLUMN(19).NE.0) THEN
          PHASE(REFLEC) = GTREAL(LINE, COLUMN(19))
        END IF  
        IF (COLUMN(20).NE.0) THEN
          FOM(REFLEC) = GTREAL(LINE, COLUMN(20))
        END IF      
        IF (USEHL.EQV..TRUE.) THEN
          HL(1,REFLEC) = GTREAL(LINE, COLUMN(21))
          HL(2,REFLEC) = GTREAL(LINE, COLUMN(22))
          HL(3,REFLEC) = GTREAL(LINE, COLUMN(23))
          HL(4,REFLEC) = GTREAL(LINE, COLUMN(24))          
        END IF 
        IF (DEFFOM.EQV..TRUE.) THEN
          FOM(REFLEC) = 0.71
        ENDIF 
        
C-------Status flag
        IF (HVFREE.EQV..TRUE.) THEN
          RFREE(REFLEC) = GTFREE(LINE, COLUMN(8))
        END IF

C-----Halt on unknown status flag
        IF (RFREE(REFLEC).EQ.'a') GO TO 992

C-----Remove reflection with I(F) and sigmaI(F) zero
        IF (USEI.EQV..FALSE.) THEN
          IF (SF(REFLEC).EQ.0.000.AND.SIGMAF(REFLEC).EQ.0.000) THEN
            REFLEC=REFLEC-1
            GO TO 994
          END IF
        ELSE
          IF (INTENS(REFLEC).EQ.0.000.AND.SIGMAI(REFLEC).EQ.0.000) THEN
            REFLEC=REFLEC-1
            GO TO 994
          END IF
        END IF

410   CONTINUE
499   WRITE(6,*) REFLEC,' reflections loaded'
      WRITE(6,*) ' '
C----------------Done loading reflections------------------------------C


C------------------Check for fixable errors----------------------------C
500   IF (WAVEL.GE.9) THEN
        WRITE(6,*) 'Unlikely diffraction wavelength (', WAVEL, 'A)!'
        WRITE(6,*) 'It will not be used!'
        WRITE(6,*) ' '
        WAVEL = -10
      END IF

C-----If I to F conversion is forced, switch off USEI and USESI
      IF (FORI2F.EQV..TRUE.) THEN
        USEI  = .FALSE.
        USESI = .FALSE.
        IF (ISAD.EQ.2) THEN
          ISAD = 1
        END IF  
      END IF

C-----Check sigma values for errors only when using amplitudes
      IF (USESI.EQV..FALSE.) THEN
        CALL SIGFIX(SIGMAF, REFLEC)
      END IF

C-----Check for status flag problems
      FIXNEG = .FALSE.
      FRINFO = .FALSE.
      ISCCP4 = .FALSE.
      MINONE = 0
      CCP4RF = 0

C-----Check for multiple wavelength datasets
      DO 505, J=1, REFLEC
        IF (WAVLID(J).NE.WAVLID(1)) THEN
C-----Second wavelength id found only use reflections up to here
          REFLEC = J-1
	  WRITE(6,*) 'Second wavelength dataset found!'
          WRITE(6,*) 'Only using the first ', REFLEC, 'reflections!'
	  WRITE(6,*) ' ' 
          GO TO 506
	END IF
505   CONTINUE

C-----If no Rfree-column exists, skip validation steps
506      IF (HVFREE.EQV..FALSE.) GO TO 540

      DO 510, J=1, REFLEC
C-----Detect alternate status flag usage
        IF (RFREE(J).EQ.'b') THEN
          FIXNEG = .TRUE.
          MINONE = MINONE+1
C-----Detect and fix CCP4-style reflection binning
	ELSE IF (RFREE(J).EQ.'c') THEN
	  CCP4RF = CCP4RF+1
	END IF
510   CONTINUE

C-----Check for information content
        DO 520, J=1, REFLEC
          IF (RFREE(J).NE.RFREE(1)) THEN
C-----Information found, skip to next test
            FRINFO = .TRUE.
	    GO TO 540
	  END IF
520     CONTINUE
C-----Halt if no information content
      IF (FRINFO.EQV..FALSE.) GO TO 990
      

C-----Can alternate status flag problem be fixed?
540   IF (FIXNEG.EQV..TRUE.) THEN
        IF((REFLEC/MINONE).GE.5) THEN
	  IF (MINONE.LT.5) THEN
	    WRITE(6,*) 'Warning! ',MINONE,' Strange reflections found.'
	    WRITE(6,*) 'They will be flagged as unmeasured.'
	    WRITE(6,*) ' '
	    CALL REFKIL(RFREE, REFLEC, 'b')
	  ELSE
	    WRITE(6,*) 'Warning! Alternate status flag usage found.'
	    WRITE(6,*) ' '
            CALL NEGFIX(RFREE, REFLEC)
	  END IF
        ELSE
          GO TO 992
        END IF
      END IF
      
C-----Must CCP4 style Rfree set be fixed?
      IF (CCP4RF.EQ.0) GO TO 550
      IF((REFLEC/CCP4RF).LE.2) THEN    
          WRITE(6,*) 'Warning! CCP4-style reflection bins detected. ',
     +'Reflections flaged 0 are assumed to make up the R-free set.'
	  CALL CCPFIX(RFREE, REFLEC)
      ELSE
        WRITE(6,*) 'Warning! ',CCP4RF,' Strange reflections found.'
	WRITE(6,*) 'They will be flagged as unmeasured.'
	WRITE(6,*) ' ' 
	CALL REFKIL(RFREE, REFLEC, 'c')
      END IF
      
C-----Do R-free sanity check  
550   FRCNT  = 0
      IF ((HVFREE.EQV..TRUE.).AND.(INSANE.EQV..FALSE.)) THEN
        DO 555, J=1, REFLEC
          IF (RFREE(J).EQ.'f') THEN
	    FRCNT = FRCNT+1
	  END IF
555     CONTINUE
        FRFRAC = REAL(FRCNT)/REAL(REFLEC)*100.0
        WRITE(6,*) 'R-free sanity check:'
	WRITE(UNIT=6, FMT='(A,I7)')   ' Reflections: ', REFLEC
	WRITE(UNIT=6, FMT='(A,I7)')   ' R-free set : ', FRCNT
        WRITE(UNIT=6, FMT='(A,F4.1,A)') ' Percentage :   ', FRFRAC,'%'
	WRITE(6,*) ' '
	IF (FRCNT.EQ.0) THEN
	  WRITE(6,*) 'Warning: R-free set empty!' 
	  WRITE(6,*) 'Not writing out status column.'
	  HVFREE = .FALSE.
	ELSE IF (((FRCNT.LT.500).AND.(FRFRAC.LT.10.0)).OR.
     +  (FRFRAC.LT.1.0)) THEN
	  WRITE(6,*) 'Warning: too small R-free set.' 
	  WRITE(6,*) 'Not writing out status column.'
	  HVFREE = .FALSE.        
	ELSE IF (FRCNT.LT.1000) THEN
	  WRITE(6,*) 'Warning: very small R-free set.'  
        ELSE IF (FRFRAC.GT.75.0) THEN
	  WRITE(6,*) 'Warning: work and test set seem to be swapped.'
	  WRITE(6,*) 'They will be un-swapped.'
	  WRITE(6,*) ' '
	  CALL CCPFIX(RFREE, REFLEC)
	  GO TO 550
	ELSE IF ((FRFRAC.LE.75.0).AND.(FRFRAC.GE.25.0)) THEN
C---------There is something wrong 
	  GO TO 992
	END IF
      END IF

C-----Deal with suspicious sigF and sigI data
      IF (DEFSIG.EQV..FALSE.) THEN
C       Calculate average I/sigI or F/sigF and number of low signal-to-
C       noise reflections
        DO 560, I=1, REFLEC
          IF (USESI.EQV..TRUE.) THEN
            IF (SIGMAI(I).GT.0.0) THEN
              RATSUM = RATSUM + INTENS(I)/SIGMAI(I)
            ELSE
              IGNCNT = IGNCNT+1
            END IF
            IF (SIGMAI(I).GT.INTENS(I)) THEN
              LOWS2N = LOWS2N+1
            END IF
          ELSE 
            RATSUM = RATSUM + SF(I)/SIGMAF(I)
            IF ((2.0*SIGMAF(I)).GT.SF(I)) THEN
              LOWS2N = LOWS2N+1
            END IF
          ENDIF
560     CONTINUE

C       Calculate signal to noise
        SIGNAL     = RATSUM/(REFLEC - IGNCNT)
        LOWS2NFRAC = 1.0*LOWS2N/REFLEC

C       Report
        WRITE(6,*) 'Signal-to-noise sanity check:' 
        IF (USESI.EQV..TRUE.) THEN
          WRITE(6, FMT=905) ' Average I/sigI           :', SIGNAL
        ELSE
          WRITE(6, FMT=905) ' Average F/sigF           :', SIGNAL 
        END IF   
        WRITE(6, FMT=905) ' Poor reflections fraction:', LOWS2NFRAC
C       For intensities... 
        IF ((USESI.EQV..TRUE.).AND.
     +      ((SIGNAL.LT.0.85).OR.(LOWS2NFRAC.GT.0.55))) THEN
          WRITE(6,*) 'Suspiciously poor I/sigI!'
          WRITE(6,*) 'Consider ignoring sigI using the -g switch!' 
C       For amplitudes...
        ELSE IF ((USESI.EQV..FALSE.).AND.
     +           ((SIGNAL.LT.1.70).OR.(LOWS2NFRAC.GT.0.55))) THEN
          WRITE(6,*) 'Suspiciously poor F/sigF!'
          WRITE(6,*) 'Consider ignoring sigI using the -g switch!' 
        END IF 
       WRITE(6,*) ' '
      END IF
      
C-----Check the anomalous data for information content
      IF (USEA.EQV..TRUE.) THEN
C       Sum the anomalous difference 
        SDANO = 0.0
        IF (ISAD.EQ.1) THEN
          DO 570, I=1, REFLEC
            SDANO = SDANO + ABS(FPLUS(I)-FMIN(I))
570       CONTINUE
        ELSE 
          DO 575, I=1, REFLEC
            IF (IPLUS(I).EQ.0.0 .OR. IMIN(I).EQ.0.0) THEN
C             One of the obeservations has no data. Ignore!
            ELSE
              SDANO = SDANO + ABS(IPLUS(I)-IMIN(I))
            ENDIF  
575       CONTINUE          
        ENDIF 

C       Report the results
        WRITE(6,*) 'Anomalous sanity check:'
        WRITE(UNIT=6,FMT=912) ' Cumulative anomalous difference: ',SDANO
        IF (SDANO.GT.0.1) THEN
          WRITE(6,*) 'There is some anomalous signal'
        ELSE
          WRITE(6,*) 'There is NO anomalous signal!'
          WRITE(6,*) 'Anomalous data will be ignored'
          USEA = .FALSE.
        ENDIF
        WRITE(6,*) ' '
      ENDIF
      
      
C------------------------Output section--------------------------------C
C-----Open output file
600   STATUS=0
      CALL GETARG(2+EXTRA, OUTCIF)
      OPEN(UNIT=8, FILE=OUTCIF, IOSTAT=STATUS)
      IF(STATUS .NE. 0) THEN
        WRITE(6,*)'Output file cannot be created!'
        GO TO 999
      END IF

C-----Write out new mmCIF file
      WRITE(6,*)'Writing data to file.'

C-----Header
      WRITE(8, FMT='(A)') 'data_rredosf'
      WRITE(8, FMT='(A)') '#'
      IF (ISCCP4.EQV..TRUE.) THEN
        WRITE(8, FMT='(A)')'#CCP4-style reflection bins in input file.'
        WRITE(8, FMT='(A)')'#R-free set is based on educated guessing.'
        WRITE(8, FMT='(A)')'#'
      END IF
      WRITE(8, FMT='(A,A)') '_cell.entry_id ', CELLID
      WRITE(8, FMT='(A,F9.3)') '_cell.length_a  ', AAXIS
      WRITE(8, FMT='(A,F9.3)') '_cell.length_b  ', BAXIS
      WRITE(8, FMT='(A,F9.3)') '_cell.length_c  ', CAXIS
      WRITE(8, FMT='(A,F7.3)') '_cell.angle_alpha ', ALPHA
      WRITE(8, FMT='(A,F7.3)') '_cell.angle_beta  ', BETA
      WRITE(8, FMT='(A,F7.3)') '_cell.angle_gamma ', GAMMA
      WRITE(8, FMT='(A)') ' '
      WRITE(8, FMT='(A,A)') '_symmetry.entry_id              ', SYMMID
      WRITE(8, FMT='(A,A,A,A)') '_symmetry.space_group_name_H-M ',
     +               "'", SYMMHM(1:LENSTR(SYMMHM)), "'"
      WRITE(8, FMT='(A)') ' '     
      IF (WAVEL.GT.0) THEN
        WRITE(8, FMT='(A)') '_diffrn_radiation_wavelength.id         1'
        WRITE(8, FMT=903)'_diffrn_radiation_wavelength.wavelength',WAVEL
        WRITE(8, FMT='(A)')'#'
      END IF
      WRITE(8, FMT='(A)') '_reflns_scale.group_code   1'
      WRITE(8, FMT='(A)')'#'   
      WRITE(8, FMT='(A)') 'loop_'
      WRITE(8, FMT='(A)') '_refln.wavelength_id'
      WRITE(8, FMT='(A)') '_refln.crystal_id'
      WRITE(8, FMT='(A)') '_refln.scale_group_code'
      WRITE(8, FMT='(A)') '_refln.index_h'
      WRITE(8, FMT='(A)') '_refln.index_k'
      WRITE(8, FMT='(A)') '_refln.index_l'
      IF (USEI.EQV..TRUE.) THEN
        WRITE(8, FMT='(A)') '_refln.intensity_meas'
        WRITE(8, FMT='(A)') '_refln.intensity_sigma'
      ELSE
        WRITE(8, FMT='(A)') '_refln.F_meas_au'
        WRITE(8, FMT='(A)') '_refln.F_meas_sigma_au'
      END IF
      IF (HVFREE.EQV..TRUE.) THEN
        WRITE(8, FMT='(A)') '_refln.status'
      END IF
      IF (USEA.EQV..TRUE.) THEN
        IF (ISAD.EQ.1) THEN
          WRITE(8, FMT='(A)') '_refln.pdbx_F_plus'
          WRITE(8, FMT='(A)') '_refln.pdbx_F_plus_sigma'
          WRITE(8, FMT='(A)') '_refln.pdbx_F_minus'
          WRITE(8, FMT='(A)') '_refln.pdbx_F_minus_sigma'
        ELSE
          WRITE(8, FMT='(A)') '_refln.pdbx_I_plus'
          WRITE(8, FMT='(A)') '_refln.pdbx_I_plus_sigma'
          WRITE(8, FMT='(A)') '_refln.pdbx_I_minus'
          WRITE(8, FMT='(A)') '_refln.pdbx_I_minus_sigma'
        END IF
      END IF  
      IF (COLUMN(19).NE.0) THEN
        WRITE(8, FMT='(A)') '_refln.phase_meas'
C       Only write figure-of-merit with phases         
        IF ((COLUMN(20).NE.0).OR.(DEFFOM.EQV..TRUE.)) THEN
          WRITE(8, FMT='(A)') '_refln.fom'
        END IF  
      END IF     
      IF (USEHL.EQV..TRUE.) THEN
        WRITE(8, FMT='(A)') '_refln.pdbx_HL_A_iso'
        WRITE(8, FMT='(A)') '_refln.pdbx_HL_B_iso'
        WRITE(8, FMT='(A)') '_refln.pdbx_HL_C_iso'
        WRITE(8, FMT='(A)') '_refln.pdbx_HL_D_iso'        
      END IF 

C-----Reflections
      DO 610, L=1,REFLEC
C       Fill the data line 
        WRITE (LINE(1:25),FMT=900) WAVLID(L), FILLER, REFH(L), REFK(L),
     +    REFL(L)
        IF (USEI.EQV..FALSE.) THEN
          WRITE(LINE(26:55),FMT=901) SF(L), SIGMAF(L)
        ELSE 
          WRITE(LINE(26:55),FMT=901) INTENS(L), SIGMAI(L)
        END IF
        POS = 55
        IF (HVFREE.EQV..TRUE.) THEN
          WRITE(LINE(POS+1:POS+2),FMT=906) RFREE(L)
          POS = POS+2
        END IF  
        IF ((USEA.EQV..TRUE.).AND.(ISAD.EQ.1)) THEN
          WRITE(LINE(POS+1:POS+60),FMT=910) FPLUS(L), SFPLUS(L),
     +     FMIN(L), SFMIN(L)
          POS = POS+60
        END IF
        IF ((USEA.EQV..TRUE.).AND.(ISAD.EQ.2)) THEN
          WRITE(LINE(POS+1:POS+60),FMT=910) IPLUS(L), SIPLUS(L),
     +     IMIN(L),SIMIN(L)
          POS = POS+60
        END IF       
        IF (COLUMN(19).NE.0) THEN
          WRITE(LINE(POS+1:POS+9),FMT=907) PHASE(L)
          POS = POS+9
          IF ((COLUMN(20).NE.0).OR.(DEFFOM.EQV..TRUE.)) THEN
            WRITE(LINE(POS+1:POS+6),FMT=908) FOM(L)
            POS = POS+6
          END IF
        END IF  
        IF (USEHL.EQV..TRUE.) THEN
C         Compensate for missing values
          IF ((HL(3,L).GT.900000).OR.(HL(4,L).GT.900000)) THEN
            IF ((HL(1,L).GT.900000).OR.(HL(2,L).GT.900000)) THEN
              WRITE(LINE(POS+1:POS+36),FMT='(A)') '       ?        ? '//
     +        '       ?        ? '             
            ELSE
              WRITE(LINE(POS+1:POS+36),FMT=911) HL(1,L), HL(2,L),
     +        '       ?        ? '        
            END IF
          ELSE
            WRITE(LINE(POS+1:POS+36),FMT=909) HL(1,L), HL(2,L), HL(3,L),
     +      HL(4,L)
          ENDIF 
          POS = POS+36
        END IF  
C       Write the data line
        WRITE(8, FMT='(A)') LINE(1:POS)
        
610   CONTINUE      


C-----Write footer
      WRITE(8, FMT='(A)') '#END OF REFLECTIONS'


C-----Give success statement.
      WRITE(6,*)'Success!'
      GO TO 999

C-----Formats
900   FORMAT(I2,1X,A,3(I5,1X))
901   FORMAT(2(F14.4,1X))
902   FORMAT(A250)
903   FORMAT(A,1X,F7.5)
904   FORMAT(A,F5.2,A)
905   FORMAT(A,F7.2)
906   FORMAT(A1,1X)
907   FORMAT(F8.3,1X)
908   FORMAT(F5.3,1X)
909   FORMAT(4(F8.3,1X))
910   FORMAT(4(F14.4,1X))
911   FORMAT(2(F8.3,1X),A)
912   FORMAT(A,F10.2)

C-----Error messages
988   WRITE(6,*)'WARNING! Weird data set. CIF2CIF aborting!'
      GO TO 999
989   WRITE(6,*)'WARNING! Number of columns does not equal number of ',
     +'column labels.'
      GO TO 399
990   WRITE(6,*)'Status flag column has no useful information. ',
     +'CIF2CIF aborting. Use the "-s" flag to avoid status flag ',
     +'problems.'
      GO TO 999
991   WRITE(6,*)'Second dataset detected. This dataset will not be used'
      GO TO 499
992   WRITE(6,*)'There is something wrong with the status flags. ',
     +'CIF2CIF aborting!'
      GO TO 999
993   WRITE(6,*)'More labels found than expected. Increase MAXLAB and',
     +' recomplile. Program aborting.'
      GO TO 999
994   WRITE(6,*)'Cannot use data, skipping line:'
      WRITE(6,*)LINE(1:60)
      GO TO 400
995   WRITE(6,*)'No experimental sigmas found. Default sigmas used.'
      GO TO 298
996   WRITE(6,*)'No structure factors or intensities, aborting!'
      GO TO 999
997   WRITE(6,*)'H, K, L columns not found, aborting!'
      GO TO 999
998   WRITE(6,*)'No labels found! Is this realy an mmCIF file?'

C-----End of the line
999   CLOSE (7)
      CLOSE (8)
      END

C-----------------------------------------------------------------------
C  This subroutine assigns the right column numbers to the data types
C-----------------------------------------------------------------------
      SUBROUTINE COLASS(LABELS, LABCNT, COLUMN, MAXCOL)
      IMPLICIT  NONE
C-----Declare variables and constants
      INTEGER   LABCNT, MAXCOL, COLUMN(MAXCOL), NU(MAXCOL),NUCNT, I
      CHARACTER LABELS(LABCNT)*(*)

C-----Initialise
      NUCNT = 0

C-----Loopt through column labels
      WRITE(6,*) ''
      WRITE(6,*) 'Data    Column'


      DO 10, I=1, LABCNT
        IF (LABELS(I)(1:14).EQ.'_refln.index_h') THEN 
          COLUMN(1) = I
          WRITE(UNIT=6,FMT=20) 'H      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln_index_h') THEN
          COLUMN(1) = I
          WRITE(UNIT=6,FMT=20) 'H      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln.index_k') THEN
          COLUMN(2) = I
          WRITE(UNIT=6,FMT=20) 'K      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln_index_k') THEN
          COLUMN(2) = I
          WRITE(UNIT=6,FMT=20) 'K      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln.index_l') THEN
          COLUMN(3) = I
          WRITE(UNIT=6,FMT=20) 'L      ',I
        ELSE IF (LABELS(I)(1:14).EQ.'_refln_index_l') THEN
          COLUMN(3) = I
          WRITE(UNIT=6,FMT=20) 'L      ',I
        ELSE IF (LABELS(I)(1:16).EQ.'_refln.F_meas_au') THEN
          COLUMN(4) = I
          WRITE(UNIT=6,FMT=20) 'F      ',I
        ELSE IF (LABELS(I)(1:16).EQ.'_refln.F_meas   ') THEN
          COLUMN(4) = I
          WRITE(UNIT=6,FMT=20) 'F      ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.F_meas_sigma_au') THEN
          COLUMN(5) = I
          WRITE(UNIT=6,FMT=20) 'SIGF   ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.F_meas_sigma   ') THEN
          COLUMN(5) = I
          WRITE(UNIT=6,FMT=20) 'SIGF   ',I
        ELSE IF (LABELS(I)(1:21).EQ.'_refln.intensity_meas') THEN
          COLUMN(6) = I
          WRITE(UNIT=6,FMT=20) 'I      ',I
        ELSE IF (LABELS(I)(1:21).EQ.'_refln.F_squared_meas') THEN
          COLUMN(6) = I
          WRITE(UNIT=6,FMT=20) 'I      ',I
        ELSE IF (LABELS(I)(1:21).EQ.'_refln_F_squared_meas') THEN
          COLUMN(6) = I
          WRITE(UNIT=6,FMT=20) 'I      ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.intensity_sigma') THEN
          COLUMN(7) = I
          WRITE(UNIT=6,FMT=20) 'SIGI   ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.F_squared_sigma') THEN
          COLUMN(7) = I
          WRITE(UNIT=6,FMT=20) 'SIGI   ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln_F_squared_sigma') THEN
          COLUMN(7) = I
          WRITE(UNIT=6,FMT=20) 'SIGI   ',I
        ELSE IF (LABELS(I)(1:16).EQ.'_refln.status') THEN
          COLUMN(8) = I
          WRITE(UNIT=6,FMT=20) 'STATUS ',I
        ELSE IF (LABELS(I)(1:22).EQ.'_refln.observed_status') THEN
          COLUMN(8) = I
          WRITE(UNIT=6,FMT=20) 'STATUS ',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.wavelength_id') THEN
          COLUMN(9) = I
          WRITE(UNIT=6,FMT=20) 'WAVE ID',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_F_plus  ') THEN
          COLUMN(10)= I
          WRITE(UNIT=6,FMT=20) 'F(+)   ',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.ccp4_F_plus  ') THEN
          COLUMN(10)= I
          WRITE(UNIT=6,FMT=20) 'F(+)   ',I
        ELSE IF (LABELS(I)(1:24).EQ.'_refln.pdbx_F_plus_sigma') THEN
          COLUMN(11)= I
          WRITE(UNIT=6,FMT=20) 'SIGF(+)',I
        ELSE IF (LABELS(I)(1:24).EQ.'_refln.ccp4_F_plus_sigma') THEN
          COLUMN(11)= I
          WRITE(UNIT=6,FMT=20) 'SIGF(+)',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_F_minus ') THEN
          COLUMN(12)= I
          WRITE(UNIT=6,FMT=20) 'F(-)   ',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.ccp4_F_minus ') THEN
          COLUMN(12)= I
          WRITE(UNIT=6,FMT=20) 'F(-)   ',I
        ELSE IF (LABELS(I)(1:25).EQ.'_refln.pdbx_F_minus_sigma') THEN
          COLUMN(13)= I
          WRITE(UNIT=6,FMT=20) 'SIGF(-)',I
        ELSE IF (LABELS(I)(1:25).EQ.'_refln.ccp4_F_minus_sigma') THEN
          COLUMN(13)= I
          WRITE(UNIT=6,FMT=20) 'SIGF(-)',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_I_plus  ') THEN
          COLUMN(14)= I
          WRITE(UNIT=6,FMT=20) 'I(+)   ',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.ccp4_I_plus  ') THEN
          COLUMN(14)= I
          WRITE(UNIT=6,FMT=20) 'I(+)   ',I
        ELSE IF (LABELS(I)(1:24).EQ.'_refln.pdbx_I_plus_sigma') THEN
          COLUMN(15)= I
          WRITE(UNIT=6,FMT=20) 'SIGI(+)',I
        ELSE IF (LABELS(I)(1:24).EQ.'_refln.ccp4_I_plus_sigma') THEN
          COLUMN(15)= I
          WRITE(UNIT=6,FMT=20) 'SIGI(+)',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_I_minus ') THEN
          COLUMN(16)= I
          WRITE(UNIT=6,FMT=20) 'I(-)   ',I
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.ccp4_I_minus ') THEN
          COLUMN(16)= I
          WRITE(UNIT=6,FMT=20) 'I(-)   ',I
        ELSE IF (LABELS(I)(1:25).EQ.'_refln.pdbx_I_minus_sigma') THEN
          COLUMN(17)= I
          WRITE(UNIT=6,FMT=20) 'SIGI(-)',I
        ELSE IF (LABELS(I)(1:25).EQ.'_refln.ccp4_I_minus_sigma') THEN
          COLUMN(17)= I
          WRITE(UNIT=6,FMT=20) 'SIGI(-)',I
        ELSE IF (LABELS(I)(8:35).EQ.'_radiation_wavelength.wavele') THEN
          COLUMN(18)= I
          WRITE(UNIT=6,FMT=21) 'Wavelength',I   
        ELSE IF (LABELS(I)(8:33).EQ.'_radiation.pdbx_wavelength') THEN
          COLUMN(18)= I
          WRITE(UNIT=6,FMT=21) 'Wavelength',I  
        ELSE IF (LABELS(I)(1:17).EQ.'_refln.phase_meas') THEN
          COLUMN(19)= I
          WRITE(UNIT=6,FMT=21) 'Phase  ',I
        ELSE IF (LABELS(I)(1:10).EQ.'_refln.fom') THEN
          COLUMN(20)= I
          WRITE(UNIT=6,FMT=21) 'FOM    ',I    
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_HL_A_iso') THEN
          COLUMN(21)= I
          WRITE(UNIT=6,FMT=21) 'HL_A   ',I     
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_HL_B_iso') THEN
          COLUMN(22)= I
          WRITE(UNIT=6,FMT=21) 'HL_B   ',I       
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_HL_C_iso') THEN
          COLUMN(23)= I
          WRITE(UNIT=6,FMT=21) 'HL_C   ',I     
        ELSE IF (LABELS(I)(1:20).EQ.'_refln.pdbx_HL_D_iso') THEN
          COLUMN(24)= I
          WRITE(UNIT=6,FMT=21) 'HL_D   ',I                
        ELSE
          NUCNT = NUCNT+1
          NU(NUCNT)=I
        END IF
10    CONTINUE

C-----Report unused columns:
      IF (NUCNT.GT.0) THEN
        WRITE(6,*) ' '
        WRITE(6,*) 'Unused data columns:'
        DO 15, I=1, NUCNT
          WRITE(6,*) LABELS(NU(I))
15      CONTINUE
        WRITE(6,*) ' '
      END IF

20    FORMAT(X,A7,X,I2)
21    FORMAT(X,A,X,I2)
      RETURN
      END

C-----------------------------------------------------------------------
C  This subroutine inserts a space between values that are not space 
C  delimited
C-----------------------------------------------------------------------
      SUBROUTINE MINFIX(LINE)
      IMPLICIT  NONE
      CHARACTER LINE*170, ONETEN*10, PROBE*2
      INTEGER   I, MINUS
      PARAMETER (ONETEN='1234567890')
      
      MINUS = 0
      DO 10, I=1, 10
        PROBE = ONETEN(I:I)//'-'
        MINUS = INDEX(LINE,PROBE)
        IF (MINUS.NE.0) THEN
          WRITE(6,*) 'Columns not properly separated:'
          WRITE(6,*) LINE
          LINE = LINE(1:MINUS)//' '//LINE(MINUS+1:169)
        END IF 
10    CONTINUE
      
30    RETURN
      END
      
C-----------------------------------------------------------------------
C  This subroutine replaces tab-characters with spaces
C-----------------------------------------------------------------------
      SUBROUTINE TABFIX(LINE)
      IMPLICIT  NONE
      CHARACTER LINE*170
      INTEGER   I, TABS
      
      DO 10, I=1, 10
        TABS = INDEX(LINE,CHAR(9))
        IF (TABS.NE.0) THEN
          WRITE(6,*) 'Tab-character detected'
          LINE = LINE(1:TABS-1)//' '//LINE(TABS+1:170)
	ELSE
	  GO TO 30
        END IF 
10    CONTINUE
      
30    RETURN
      END      

C-----------------------------------------------------------------------
C  This subroutine replaces ;-characters with spaces
C-----------------------------------------------------------------------
      SUBROUTINE SECFIX(LINE)
      IMPLICIT  NONE
      CHARACTER LINE*170
      INTEGER   I, SEMCOL
      
      DO 10, I=1, 10
        SEMCOL = INDEX(LINE,';')
        IF (SEMCOL.NE.0) THEN
          WRITE(6,*) 'Semi-colon detected'
          LINE = LINE(1:SEMCOL-1)//' '//LINE(SEMCOL+1:170)
	ELSE
	  GO TO 30
        END IF 
10    CONTINUE
      
30    RETURN
      END   

C-----------------------------------------------------------------------
C  This subroutine replaces '-characters with spaces
C-----------------------------------------------------------------------
      SUBROUTINE QUOFIX(LINE)
      IMPLICIT  NONE
      CHARACTER LINE*170
      INTEGER   I, QUOTE
      
      DO 10, I=1, 10
        QUOTE = INDEX(LINE,"'")
        IF (QUOTE.NE.0) THEN
          WRITE(6,*) 'Single quote detected'
          LINE = LINE(1:QUOTE-1)//' '//LINE(QUOTE+1:170)
	ELSE
	  GO TO 30
        END IF 
10    CONTINUE
      
30    RETURN
      END   
C-----------------------------------------------------------------------
C  This subroutine looks for sigma values of 0.00. They are set to the 
C  highest sigma value + 0.01   
C-----------------------------------------------------------------------
      SUBROUTINE SIGFIX(SIGMAF, REFLEC)
      IMPLICIT   NONE
      INTEGER          REFLEC, I, J, K
      DOUBLE PRECISION SIGMAF(REFLEC), HIVAL

C----Check to see whether there is any information in the sigma values  
      HIVAL = SIGMAF(1)
      DO 10, I=2, REFLEC 
	IF ((SIGMAF(I).NE.HIVAL).AND.(SIGMAF(I).NE.0.0000)) THEN
C-------There is information, skip warning and find highest value 	
	  GO TO 20
	END IF
10    CONTINUE

C-----No-information warning
      WRITE(6,*) 'All sigma values equal, no information content!'

C-----Find the highest sigma value
20    DO 30, J=2, REFLEC
        HIVAL = MAX(SIGMAF(J), HIVAL)
30    CONTINUE

C-----Set all sigmas with value 0.0000 (or lower) to HIVAL+0.01
      DO 40, K=1, REFLEC
        IF (SIGMAF(K).LE.0.0000) THEN
          SIGMAF(K) = HIVAL + 0.01
        END IF
40    CONTINUE
      
      RETURN
      END
      
C-----------------------------------------------------------------------
C  This subroutine fixes the alternative status flag usage
C-----------------------------------------------------------------------
      SUBROUTINE NEGFIX(RFREE, REFLEC)
      IMPLICIT   NONE
      INTEGER    REFLEC, I
      CHARACTER  RFREE(REFLEC)*1
      
      DO 10, I=1, REFLEC
        IF (RFREE(I).EQ.'f') THEN
          RFREE(I)='o'
        ELSE 
          IF (RFREE(I).EQ.'b') THEN
            RFREE(I)='f'
          END IF
        END IF 
10    CONTINUE
      
30    RETURN
      END
      
C-----------------------------------------------------------------------
C  This subroutine marks all reflections with a certain flag as 
C  umeasured
C-----------------------------------------------------------------------
      SUBROUTINE REFKIL(RFREE, REFLEC, FLAG)
      IMPLICIT   NONE
      INTEGER    REFLEC, I
      CHARACTER  RFREE(REFLEC)*1, FLAG
      
      DO 10, I=1, REFLEC
        IF (RFREE(I).EQ.FLAG) THEN
          RFREE(I)='x'
        END IF 
10    CONTINUE
      
30    RETURN
      END      

C-----------------------------------------------------------------------
C  This subroutine fixes the CCP4-style status flag usage
C-----------------------------------------------------------------------
      SUBROUTINE CCPFIX(RFREE, REFLEC)
      IMPLICIT   NONE
      INTEGER    REFLEC, I
      CHARACTER  RFREE(REFLEC)*1
      
      DO 10, I=1, REFLEC
        IF ((RFREE(I).EQ.'f').OR.(RFREE(I).EQ.'c')) THEN
          RFREE(I)='o'
        ELSE 
          IF (RFREE(I).EQ.'o') THEN
            RFREE(I)='f'
          END IF
        END IF 
10    CONTINUE
      
30    RETURN
      END

C-----------------------------------------------------------------------
C  This function extracts an integer from column CLMNR from a text line
C-----------------------------------------------------------------------
      INTEGER FUNCTION GETINT(LINE, CLMNR)
      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*250 LINE, TLINE, CHOP
      INTEGER   CLMNR, J, SPACE, VALUE

C-----Initialise
      VALUE = -999999
      SPACE = 1   
      
C-----Extract right column
      TLINE = LINE
      DO 10, J=1, CLMNR
        TLINE = TLINE(SPACE:250)     
        TLINE = CHOP(TLINE, 250)
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
      CHARACTER*250 LINE, TLINE, CHOP
      CHARACTER*1   QMARK
      INTEGER   CLMNR, J, SPACE
      DOUBLE PRECISION VALUE

C-----Initialise
      SPACE = 1
      VALUE = -999999.9         
      
C-----Extract right column
      TLINE = LINE
      DO 10, J=1, CLMNR
        TLINE = TLINE(SPACE:250)
C       PRINT*, TLINE
	TLINE = CHOP(TLINE, 250)
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
C  This function merges anomalous intensities.
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION IMERGE(AP, AM, SAP, SAM)
      IMPLICIT NONE
C-----Declare variables and constants
      DOUBLE PRECISION AP, AM, SAP, SAM, WP, WM
C AP   I(+)
C AM   I(-)
C SAP  sigmaI(+)
C SAM  sigmaI(-)
C WP   1/sigmaI(+)**2
C WM   1/sigmaI(1)**2

C-----Contingency for missing values
      IF (AP.EQ.0.0 .AND. SAP.EQ.0.0) THEN
        IMERGE = AM
        RETURN
      ELSE IF (AM.EQ.0.0 .AND. SAM.EQ.0.0) THEN
        IMERGE = AP
        RETURN
      ENDIF

C-----Calculate weights
      WP = 1.0/(SAP*SAP)
      WM = 1.0/(SAM*SAM)

C-----Calculate merged I
      IMERGE = (WP*AP + WM*AM)/(WP + WM)      

C-----Return     
      RETURN
      END

C-----------------------------------------------------------------------
C  This function merges anomalous intensitie sigma values.
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION SIMERGE(AP, AM, SAP, SAM)
      IMPLICIT NONE
C-----Declare variables and constants
      DOUBLE PRECISION AP, AM, SAP, SAM, WP, WM, T1, T2
C AP   I(+)
C AM   I(-)
C SAP  sigmaI(+)
C SAM  sigmaI(-)
C WP   1/sigmaI(+)**2
C WM   1/sigmaI(1)**2

C-----Contingency for missing values
      IF (AP.EQ.0.0 .AND. SAP.EQ.0.0) THEN
        SIMERGE = SAM
        RETURN
      ELSE IF (AM.EQ.0.0 .AND. SAM.EQ.0.0) THEN
        SIMERGE = SAP
        RETURN
      ENDIF

C-----Calculate weights
      WP = 1.0/(SAP*SAP)
      WM = 1.0/(SAM*SAM)

C-----Calculate temporary values T1 and T2
      T1 = (WP*SAP*SAP+WM*SAM*SAM)/(WP+WM)
      T2 = (WP*WM*(AP-AM)*(AP-AM))/((WP+WM)*(WP+WM))

C-----Calculate merged sigma I
      SIMERGE = SQRT(T1 + T2)      

C-----Return     
      RETURN
      END

C-----------------------------------------------------------------------
C  This function extracts a string from column CLMNR from a text 
C  line. This string is then converted to a status flag. 
C  Non default flags are used when certain errors are encountered:
C   a: Label unknown or no label at all. This is a fatal problem.
C   b: Possible use of alternative status flags. Problem will be solved 
C      in the main program if possible.
C   c: CCP4-style reflection binning. Problem will be solved if fewer 
C      than 34 bins are used. The bin flagged '0' is assumed to be the
C      R-free set (CCP4 default)
C-----------------------------------------------------------------------
      CHARACTER*1 FUNCTION GTFREE(LINE, CLMNR)
      IMPLICIT   NONE
C-----Declare variables and constants
      CHARACTER LINE*250, TLINE*253, CHOP*253, FREE*2
      INTEGER    CLMNR, I, J, SPACE
      REAL       CCP4
      
C-----Append some text to avoid crashes caused by missing columns
      TLINE = LINE//' aa'

C-----Extract right column
      DO 10, J=1, CLMNR
        TLINE = CHOP(TLINE, 253)
        SPACE = INDEX(TLINE, ' ')
        READ(UNIT=TLINE(1:2), FMT=*, ERR=15) FREE
15      TLINE = TLINE(SPACE:252)
10    CONTINUE

C-----Assign and return
      IF ((FREE.EQ.'1 ').OR.(FREE.EQ.'f ').OR.(FREE.EQ.'1.').OR.
     +(FREE.EQ.'F ')) THEN   
        GTFREE = 'f'
      ELSE IF ((FREE.EQ.'0 ').OR.(FREE.EQ.'o ').OR.(FREE.EQ.'0.').OR.
     +(FREE.EQ.'O ')) THEN
        GTFREE = 'o'
C-----Warn if no column is found
      ELSE IF (FREE.EQ.'aa') THEN
        WRITE(6,*) 'Warning! Status flag column does not really exist.'
        GTFREE = 'a'     
C-----Flag for alternate use
      ELSE IF (FREE.EQ.'-1') THEN
        GTFREE = 'b'
C-----Special reflections
      ELSE IF (FREE.EQ.'l ') THEN
        GTFREE = 'l'
      ELSE IF (FREE.EQ.'h ') THEN
        GTFREE = 'h'
      ELSE IF (FREE.EQ.'< ') THEN
        GTFREE = '<'
      ELSE IF (FREE.EQ.'x ') THEN
        GTFREE = 'x'
      ELSE IF (FREE.EQ.'- ') THEN
        GTFREE = '-'
C-----Unmeasured reflections 
      ELSE IF ((FREE.EQ.'. ').OR.(FREE.EQ.'? ')) THEN
        GTFREE = 'x'
C-----CCP4 style R-free flags?
      ELSE 
        READ(UNIT=FREE, FMT=*, ERR=20) CCP4
        GO TO 30 
20      CCP4=99.0
C-----Flag other reflections as strange
30      IF ((CCP4.GE.2.0).AND.(CCP4.LE.99.0)) THEN
          GTFREE = 'c'
        ELSE
          WRITE(6,*) 'Strange status flag for reflection:', LINE(1:60)
          GTFREE = 'a'
        END IF
      END IF

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
C  This function returns the number of data entries in a line
C-----------------------------------------------------------------------
      INTEGER FUNCTION DATCHK(LINE)
      IMPLICIT NONE
      CHARACTER LINE*250, LINE1*251
      INTEGER   I, DATCNT
      
C-----Intitialise
      DATCNT = 0
      LINE1 = ' '//LINE

C-----Analise line character by character
      DO 10, I=2, 250
        IF ((LINE1(I:I).NE.' ').AND.(LINE1(I-1:I-1).EQ.' ')) THEN
	  DATCNT = DATCNT+1
        END IF
10    CONTINUE   
       
C      PRINT*, DATCNT, LINE 
      DATCHK = DATCNT      
      RETURN 
      END      

C-----------------------------------------------------------------------
C  This function reduces two lines to one
C-----------------------------------------------------------------------
      CHARACTER*250 FUNCTION DEJUNK(LINE1, LINE2)
      IMPLICIT NONE
      CHARACTER LINE1*250, LINE2*250, TMPLIN*501
      INTEGER   I, LEN1, LEN2, LEN3, MAXCUT, KILL, DUBSPA, LENSTR

C-----Find lengths
      LEN1 = LENSTR(LINE1)
      LEN2 = LENSTR(LINE2)
      TMPLIN = LINE1(1:LEN1)//' '//LINE2(1:LEN2)
      LEN3 = LENSTR(TMPLIN)
      
      IF (LEN3.LT.250) THEN
        GO TO 20
      ELSE
C-------Remove double spaces
        DO 10, I=1, LEN3
          DUBSPA = INDEX(TMPLIN, '  ')
	  IF (DUBSPA.EQ.0) THEN
	    WRITE(6,*) 'This is bound to go wrong.'
            GO TO 20
          ELSE
            TMPLIN = TMPLIN(1:DUBSPA)//TMPLIN(DUBSPA+2:500)
            TMPLIN(501:501) = ' '
          END IF
10      CONTINUE
      END IF

20    DEJUNK = TMPLIN(1:250)
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
