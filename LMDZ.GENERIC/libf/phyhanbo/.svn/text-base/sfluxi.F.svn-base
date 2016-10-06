      SUBROUTINE SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,RSFI,WNOI,DWNI,
     *                  COSBI,WBARI,GWEIGHT,NFLUXTOPI,NFLUXTOPI_nu,
     *                  FMNETI,fluxupi,fluxdni,fluxupi_nu,
     *                  FZEROI,TAUGSURF)

      use radinc_h
      use radcommon_h, only: planckir, tlimit,sigma

      implicit none

#include "comcstfi.h"

      integer NLEVRAD, L, NW, NG, NTS, NTT

      real*8 TLEV(L_LEVELS), PLEV(L_LEVELS)
      real*8 TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8 FMNETI(L_NLAYRAD)
      real*8 WNOI(L_NSPECTI), DWNI(L_NSPECTI)
      real*8 DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 FMUPI(L_NLEVRAD), FMDI(L_NLEVRAD)
      real*8 COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 GWEIGHT(L_NGAUSS), NFLUXTOPI
      real*8 NFLUXTOPI_nu(L_NSPECTI)
      real*8 fluxupi_nu(L_NLAYRAD,L_NSPECTI)
      real*8 FTOPUP

      real*8 UBARI, RSFI, TSURF, BSURF, TTOP, BTOP, TAUTOP
      real*8 PLANCK, PLTOP
      real*8 fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD)
      real*8 FZEROI(L_NSPECTI)
      real*8 taugsurf(L_NSPECTI,L_NGAUSS-1), fzero

      real*8 fup_tmp(L_NSPECTI),fdn_tmp(L_NSPECTI)
      real*8 PLANCKSUM,PLANCKREF


C======================================================================C
 
      NLEVRAD = L_NLEVRAD
 

C     ZERO THE NET FLUXES
    
      NFLUXTOPI = 0.0D0

      DO NW=1,L_NSPECTI
        NFLUXTOPI_nu(NW) = 0.0D0
        DO L=1,L_NLAYRAD
           FLUXUPI_nu(L,NW) = 0.0D0

              fup_tmp(nw)=0.0D0
              fdn_tmp(nw)=0.0D0

        END DO
      END DO

      DO L=1,L_NLAYRAD
        FMNETI(L)  = 0.0D0
        FLUXUPI(L) = 0.0D0
        FLUXDNI(L) = 0.0D0
      END DO
 
C     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE INFRARED
C     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      TTOP  = TLEV(2)  ! JL12 why not (1) ???
      TSURF = TLEV(L_LEVELS)

      NTS   = int(TSURF*NTfac)-NTstar+1
      NTT   = int(TTOP *NTfac)-NTstar+1

!JL12 corrects the surface planck function so that its integral is equal to sigma Tsurf^4
!JL12   this ensure that no flux is lost due to:
!JL12          -truncation of the planck function at high/low wavenumber
!JL12          -numerical error during first spectral integration
!JL12          -discrepancy between Tsurf and NTS/NTfac
      PLANCKSUM=0.d0
      PLANCKREF=TSURF*TSURF
      PLANCKREF=sigma*PLANCKREF*PLANCKREF
      DO NW=1,L_NSPECTI
         PLANCKSUM=PLANCKSUM+PLANCKIR(NW,NTS)*DWNI(NW)
      ENDDO
      PLANCKSUM=PLANCKREF/(PLANCKSUM*Pi)
!JL12

      DO 501 NW=1,L_NSPECTI

C       SURFACE EMISSIONS - INDEPENDENT OF GAUSS POINTS
        BSURF = (1.-RSFI)*PLANCKIR(NW,NTS)*PLANCKSUM !JL12 plancksum see above
        PLTOP = PLANCKIR(NW,NTT)

C  If FZEROI(NW) = 1, then the k-coefficients are zero - skip to the
C  special Gauss point at the end.
 
        FZERO = FZEROI(NW)
        IF(FZERO.ge.0.99) goto 40
 
        DO NG=1,L_NGAUSS-1
         
          if(TAUGSURF(NW,NG).lt. TLIMIT) then
            fzero = fzero + (1.0D0-FZEROI(NW))*GWEIGHT(NG)
            goto 30
          end if

C         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE IR
C         CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
C         OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY
 
!          TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
          TAUTOP = TAUCUMI(2,NW,NG)
          BTOP   = (1.0D0-EXP(-TAUTOP/UBARI))*PLTOP
 
C         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 
          
          CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),
     *                TAUCUMI(1,NW,NG),
     *                WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,
     *                BSURF,FTOPUP,FMUPI,FMDI)



C         NOW CALCULATE THE CUMULATIVE IR NET FLUX

          NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*GWEIGHT(NG)*
     *                           (1.0D0-FZEROI(NW))

c         and same thing by spectral band... (RDW)
          NFLUXTOPI_nu(NW) = NFLUXTOPI_nu(NW)
     *      +FTOPUP*DWNI(NW)*GWEIGHT(NG)*(1.0D0-FZEROI(NW))


          DO L=1,L_NLEVRAD-1

C           CORRECT FOR THE WAVENUMBER INTERVALS

            FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*
     *                              GWEIGHT(NG)*(1.0D0-FZEROI(NW))
            FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*
     *                                (1.0D0-FZEROI(NW))
            FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*GWEIGHT(NG)*
     *                                (1.0D0-FZEROI(NW))

c         and same thing by spectral band... (RW)
            FLUXUPI_nu(L,NW) = FLUXUPI_nu(L,NW) + 
     *                FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*(1.0D0-FZEROI(NW))

          END DO

   30     CONTINUE

       END DO       !End NGAUSS LOOP

   40  CONTINUE

C      SPECIAL 17th Gauss point

       NG     = L_NGAUSS

!       TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
       TAUTOP = TAUCUMI(2,NW,NG)
       BTOP   = (1.0D0-EXP(-TAUTOP/UBARI))*PLTOP

C      WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C      CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C      WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 


       CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),
     *                TAUCUMI(1,NW,NG),
     *                WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,
     *                BSURF,FTOPUP,FMUPI,FMDI)
 
C      NOW CALCULATE THE CUMULATIVE IR NET FLUX

       NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*FZERO

c         and same thing by spectral band... (RW)
          NFLUXTOPI_nu(NW) = NFLUXTOPI_nu(NW)
     *      +FTOPUP*DWNI(NW)*FZERO

       DO L=1,L_NLEVRAD-1

C        CORRECT FOR THE WAVENUMBER INTERVALS

         FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*FZERO
         FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*FZERO
         FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*FZERO

c         and same thing by spectral band... (RW)
         FLUXUPI_nu(L,NW) = FLUXUPI_nu(L,NW) + FMUPI(L)*DWNI(NW)*FZERO

       END DO

  501 CONTINUE      !End Spectral Interval LOOP

C *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****

      RETURN
      END
