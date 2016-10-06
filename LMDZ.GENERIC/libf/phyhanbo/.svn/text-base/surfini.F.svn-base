      SUBROUTINE surfini(ngrid,nq,qsurf,psolaralb)

      USE surfdat_h, only: albedodat, albedice
      USE tracer_h, only: igcm_co2_ice
      use comgeomfi_h, only: lati
      use planetwide_mod, only: planetwide_maxval, planetwide_minval

      IMPLICIT NONE
c=======================================================================
c
c   creation des calottes pour l'etat initial
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------
!#include "dimensions.h"
!#include "dimphys.h"
#include "callkeys.h"
c
      INTEGER,INTENT(IN) :: ngrid
      INTEGER,INTENT(IN) :: nq
      REAL,INTENT(OUT) :: psolaralb(ngrid)
      REAL,INTENT(IN) :: qsurf(ngrid,nq) !tracer on surface (kg/m2)

      INTEGER :: ig,icap
      REAL :: min_albedo,max_albedo
c
c=======================================================================


      DO ig=1,ngrid
         psolaralb(ig)=albedodat(ig)
      ENDDO

      call planetwide_minval(albedodat,min_albedo)
      call planetwide_maxval(albedodat,max_albedo)
      write(*,*) 'surfini: minimum corrected albedo',min_albedo
      write(*,*) 'surfini: maximum corrected albedo',max_albedo

      if (igcm_co2_ice.ne.0) then
! Change Albedo if there is CO2 ice on the surface
        DO ig=1,ngrid
          IF (qsurf(ig,igcm_co2_ice) .GT. 0.) THEN
             IF(lati(ig).LT.0.) THEN
                icap=2 ! Southern hemisphere
             ELSE
                icap=1 ! Northern hemisphere
             ENDIF
             psolaralb(ig) = albedice(icap)
          END IF
        ENDDO ! of DO ig=1,ngrid     
      else
        write(*,*) "surfini: No CO2 ice tracer on surface  ..."
        write(*,*) "         and therefore no albedo change."
      endif

      call planetwide_minval(psolaralb,min_albedo)
      call planetwide_maxval(psolaralb,max_albedo)
      write(*,*) 'surfini: minimum corrected albedo',min_albedo
      write(*,*) 'surfini: maximum corrected albedo',max_albedo

      END
