      SUBROUTINE GFLUXI(NLL,TLEV,NW,DW,DTAU,TAUCUM,W0,COSBAR,UBARI,
     *                  RSF,BTOP,BSURF,FTOPUP,FMIDP,FMIDM)

C  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
C  FOR THE INFRARED FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
C  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
C  MEASURED FROM THE TOP OF EACH LAYER.  THE TOP OF EACH LAYER HAS  
C  OPTICAL DEPTH ZERO.  IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
C THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY 
C VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
C
C NLL            = NUMBER OF LEVELS (NLAYERS + 1) MUST BE LESS THAT NL (101)
C TLEV(L_LEVELS) = ARRAY OF TEMPERATURES AT GCM LEVELS
C WAVEN          = WAVELENGTH FOR THE COMPUTATION
C DW             = WAVENUMBER INTERVAL
C DTAU(NLAYER)   = ARRAY OPTICAL DEPTH OF THE LAYERS
C W0(NLEVEL)     = SINGLE SCATTERING ALBEDO
C COSBAR(NLEVEL) = ASYMMETRY FACTORS, 0=ISOTROPIC
C UBARI          = AVERAGE ANGLE, MUST BE EQUAL TO 0.5 IN IR
C RSF            = SURFACE REFLECTANCE
C BTOP           = UPPER BOUNDARY CONDITION ON IR INTENSITY (NOT FLUX)
C BSURF          = SURFACE EMISSION = (1-RSFI)*PLANCK, INTENSITY (NOT FLUX)
C FP(NLEVEL)     = UPWARD FLUX AT LEVELS
C FM(NLEVEL)     = DOWNWARD FLUX AT LEVELS
C FMIDP(NLAYER)  = UPWARD FLUX AT LAYER MIDPOINTS
C FMIDM(NLAYER)  = DOWNWARD FLUX AT LAYER MIDPOINTS
C
C----------------------------------------------------------------------C

      use radinc_h
      use radcommon_h, only: planckir

      implicit none

#include "comcstfi.h"

      INTEGER NLP
      PARAMETER (NLP=101) ! MUST BE LARGER THAN NLEVEL 

      INTEGER NLL, NLAYER, L, NW, NT, NT2
      REAL*8  TERM, CPMID, CMMID
      REAL*8  PLANCK
      REAL*8  EM,EP
      REAL*8  COSBAR(L_NLAYRAD), W0(L_NLAYRAD), DTAU(L_NLAYRAD)
      REAL*8  TAUCUM(L_LEVELS), DTAUK
      REAL*8  TLEV(L_LEVELS)
      REAL*8  WAVEN, DW, UBARI, RSF
      REAL*8  BTOP, BSURF, FMIDP(L_NLAYRAD), FMIDM(L_NLAYRAD)
      REAL*8  B0(NLP),B1(NLP),ALPHA(NLP),LAMDA(NLP),XK1(NLP),XK2(NLP)
      REAL*8  GAMA(NLP),CP(NLP),CM(NLP),CPM1(NLP),CMM1(NLP),E1(NLP)
      REAL*8  E2(NLP),E3(NLP),E4(NLP)

      REAL*8  FTOPUP, FLUXUP, FLUXDN

      real*8 :: TAUMAX = L_TAUMAX 

C======================================================================C
 
C     WE GO WITH THE HEMISPHERIC CONSTANT APPROACH IN THE INFRARED
 

      IF (NLL .GT. NLP) STOP 'PARAMETER NL TOO SMALL IN GFLUXI'

      NLAYER = L_NLAYRAD

      DO L=1,L_NLAYRAD-1

!----------------------------------------------------
! There is a problem when W0 = 1
!         open(888,file='W0')
!           if ((W0(L).eq.0.).or.(W0(L).eq.1.)) then
!             write(888,*) W0(L), L, 'gfluxi'
!           endif
! Prevent this with an if statement:
        if (W0(L).eq.1.D0) then
           W0(L) = 0.99999D0
        endif
!----------------------------------------------------

        ALPHA(L) = SQRT( (1.0D0-W0(L))/(1.0D0-W0(L)*COSBAR(L)) )
        LAMDA(L) = ALPHA(L)*(1.0D0-W0(L)*COSBAR(L))/UBARI

        NT    = int(TLEV(2*L)*NTfac)   - NTstar+1
        NT2   = int(TLEV(2*L+2)*NTfac) - NTstar+1

        B1(L) = (PLANCKIR(NW,NT2)-PLANCKIR(NW,NT))/DTAU(L)
        B0(L) = PLANCKIR(NW,NT)
      END DO

C     Take care of special lower layer

      L        = L_NLAYRAD

      if (W0(L).eq.1.) then
          W0(L) = 0.99999D0
      end if

      ALPHA(L) = SQRT( (1.0D0-W0(L))/(1.0D0-W0(L)*COSBAR(L)) )
      LAMDA(L) = ALPHA(L)*(1.0D0-W0(L)*COSBAR(L))/UBARI

      ! Tsurf is used for 1st layer source function
      ! -- same results for most thin atmospheres
      ! -- and stabilizes integrations
      NT    = int(TLEV(2*L+1)*NTfac) - NTstar+1
      !! For deep, opaque, thick first layers (e.g. Saturn)
      !! what is below works much better, not unstable, ...
      !! ... and actually fully accurate because 1st layer temp (JL) 
      !NT    = int(TLEV(2*L)*NTfac) - NTstar+1
      !! (or this one yields same results
      !NT    = int( (TLEV(2*L)+TLEV(2*L+1))*0.5*NTfac ) - NTstar+1

      NT2   = int(TLEV(2*L)*NTfac)   - NTstar+1
      B1(L) = (PLANCKIR(NW,NT)-PLANCKIR(NW,NT2))/DTAU(L)
      B0(L) = PLANCKIR(NW,NT2)

      DO L=1,L_NLAYRAD
        GAMA(L) = (1.0D0-ALPHA(L))/(1.0D0+ALPHA(L))
        TERM    = UBARI/(1.0D0-W0(L)*COSBAR(L))

C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C       BOTTOM OF THE LAYER.  THAT IS AT DTAU OPTICAL DEPTH

        CP(L) = B0(L)+B1(L)*DTAU(L) +B1(L)*TERM
        CM(L) = B0(L)+B1(L)*DTAU(L) -B1(L)*TERM

C       CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C       AT THE TOP OF THE LAYER, THAT IS ZERO OPTICAL DEPTH

        CPM1(L) = B0(L)+B1(L)*TERM
        CMM1(L) = B0(L)-B1(L)*TERM
      END DO
 
C     NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C     FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C     WARNING IF DTAU(J) IS GREATER THAN ABOUT 35 (VAX)
C     WE CLIP IT TO AVOID OVERFLOW.

      DO L=1,L_NLAYRAD

C       CLIP THE EXPONENTIAL HERE.

        EP    = EXP( MIN((LAMDA(L)*DTAU(L)),TAUMAX))
        EM    = 1.0D0/EP
        E1(L) = EP+GAMA(L)*EM
        E2(L) = EP-GAMA(L)*EM
        E3(L) = GAMA(L)*EP+EM
        E4(L) = GAMA(L)*EP-EM
      END DO
 
c     B81=BTOP  ! RENAME BEFORE CALLING DSOLVER - used to be to set
c     B82=BSURF ! them to real*8 - but now everything is real*8
c     R81=RSF   ! so this may not be necessary

C     Double precision tridiagonal solver

      CALL DSOLVER(NLAYER,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,
     *             BSURF,RSF,XK1,XK2)
 

C     NOW WE CALCULATE THE FLUXES AT THE MIDPOINTS OF THE LAYERS.
 
      DO L=1,L_NLAYRAD-1
        DTAUK = TAUCUM(2*L+1)-TAUCUM(2*L)
        EP    = EXP(MIN(LAMDA(L)*DTAUK,TAUMAX)) ! CLIPPED EXPONENTIAL 
        EM    = 1.0D0/EP
        TERM  = UBARI/(1.D0-W0(L)*COSBAR(L))

C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C       BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

        CPMID    = B0(L)+B1(L)*DTAUK +B1(L)*TERM
        CMMID    = B0(L)+B1(L)*DTAUK -B1(L)*TERM
        FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
        FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID

C       FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

        FMIDP(L) = FMIDP(L)*PI
        FMIDM(L) = FMIDM(L)*PI
      END DO
 
C     And now, for the special bottom layer

      L    = L_NLAYRAD

      EP   = EXP(MIN((LAMDA(L)*DTAU(L)),TAUMAX)) ! CLIPPED EXPONENTIAL 
      EM   = 1.0D0/EP
      TERM = UBARI/(1.D0-W0(L)*COSBAR(L))

C     CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C     BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

      CPMID    = B0(L)+B1(L)*DTAU(L) +B1(L)*TERM
      CMMID    = B0(L)+B1(L)*DTAU(L) -B1(L)*TERM
      FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
      FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID
 
C     FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

      FMIDP(L) = FMIDP(L)*PI
      FMIDM(L) = FMIDM(L)*PI

C     FLUX AT THE PTOP LEVEL

      EP   = 1.0D0
      EM   = 1.0D0
      TERM = UBARI/(1.0D0-W0(1)*COSBAR(1))

C     CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C     BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

      CPMID  = B0(1)+B1(1)*TERM
      CMMID  = B0(1)-B1(1)*TERM

      FLUXUP = XK1(1)*EP + GAMA(1)*XK2(1)*EM + CPMID
      FLUXDN = XK1(1)*EP*GAMA(1) + XK2(1)*EM + CMMID
 
C     FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

      FTOPUP = (FLUXUP-FLUXDN)*PI


      RETURN
      END
