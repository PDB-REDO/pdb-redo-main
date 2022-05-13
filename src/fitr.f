      PROGRAM FITR
C=======================================================================
C  Version 1.01 2022-04-04
C  Checks to see wheter one set of R-factors fits the other
C
C  Usage: fitr R1 Rfree1 R2 Rfree2 (cut-off)
C  (cut-off) is an optional cut-off for the difference between the two
C  sets of R-factors. The default value is 0.05
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
C  Version 1.00:
C  - Added the help function.
C  Version 1.01:
C  - Increases the tollerance for lower than expected R-factor
C=======================================================================
      IMPLICIT NONE
C-----Declare the basic variables and parameters
      INTEGER   I
      CHARACTER VERS*4
      PARAMETER (VERS='1.00')

C-----MAXLIN is the maximum number of refinement settings that can be
C     evaluated
      REAL      GTREAL, R1, RF1, R2, RF2, CUTOFF, DIFFR, DIFFRF
      CHARACTER ARGIN*10
      INTEGER   ARGS, ISFIT

C=========================== Help function=============================C
      ARGS=IARGC()
      IF (ARGS.EQ.0) THEN
      WRITE(6,*)'*****   FITR version: ',VERS,'   *****'
      WRITE(6,*) ' '
      WRITE(6,*)'Fitr checks wheter one set of R-factors fits the oth'//
     +          'er. It returns ''1'' if it does and ''0'' if it does'//
     +          ' not.' 
      WRITE(6,*)'Written by Robbie Joosten'
      WRITE(6,*)'E-mail: r.joosten@nki.nl, robbie_joosten@hotmail.com '
      WRITE(6,*) ' '
      WRITE(6,*)'Usage:'
      WRITE(6,*)'fitr R1 RFREE1 R2 RFREE2 (CUT_OFF)'
      WRITE(6,*)' '
      WRITE(6,*)'R1 is the first R-factor as a fractional number.'
      WRITE(6,*)'RFREE1 is the first free R-factor.'
      WRITE(6,*)'RFREE1 is the second R-factor.'
      WRITE(6,*)'RFREE2 is the second free R-factor.'
      WRITE(6,*)'CUT_OFF is an optional cut-off value (default = 0.05).'
      WRITE(6,*)'    '
      WRITE(6,*)'Citing FITR:'
      WRITE(6,*)'Joosten et al. ''PDB_REDO: automated re-refinement o'//
     +          'f X-ray structure models in the PDB'' J. Appl. Cryst'//
     +          '., 42, 376-384 (2009)'
        GO TO 999
      END IF

C--------------------------- Main program -----------------------------C
C-----Get input values                                  
      CALL GETARG(1, ARGIN)
      R1 = GTREAL(ARGIN)
      CALL GETARG(2, ARGIN)
      RF1 = GTREAL(ARGIN)
      CALL GETARG(3, ARGIN)
      R2 = GTREAL(ARGIN)
      CALL GETARG(4, ARGIN)
      RF2 = GTREAL(ARGIN)                 

C-----Use custom cut-off?
      IF (ARGS.GE.5) THEN
        CALL GETARG(5, ARGIN)
        CUTOFF=GTREAL(ARGIN)
      ELSE       
        CUTOFF=0.05
      END IF
      
      
C-----Check values
      DIFFR  = ABS(R1-R2)
      DIFFRF = ABS(RF1-RF2)
      IF ((DIFFR.LE.CUTOFF).OR.(DIFFRF.LE.CUTOFF)) THEN
        ISFIT = 1    
      ELSE
        ISFIT = 0
      END IF
	
C-----Write  out result
      WRITE(UNIT=6,FMT='(I1)') ISFIT

C-----End of the line
999   END

C-------------------------- Subroutines and functions -----------------C
C-----------------------------------------------------------------------
C  This function extracts a real number from a string
C-----------------------------------------------------------------------
      REAL FUNCTION GTREAL(TLINE)

      IMPLICIT NONE
C-----Declare variables and constants
      CHARACTER*10 TLINE
      INTEGER      SPACE
      REAL         VALUE

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
15    VALUE = 0.00



C-----Assign and return
99    GTREAL = VALUE
      RETURN
      END
            
