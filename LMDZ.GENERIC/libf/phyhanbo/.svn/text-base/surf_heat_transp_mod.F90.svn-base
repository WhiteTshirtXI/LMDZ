!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ocean_slab_mod.F90,v 1.3 2008-02-04 16:24:28 fairhead Exp $
!
MODULE surf_heat_transp_mod


CONTAINS

      SUBROUTINE divgrad_phy(ngrid,nlevs,temp,delta)



      IMPLICIT NONE

#include "dimensions.h"
!#include "dimphys.h"
#include "paramet.h"
#include "comcstfi.h"
#include "comgeom.h"
#include "comhdiff.h"
      
      INTEGER,INTENT(IN) :: ngrid, nlevs
      REAL,INTENT(IN) :: temp(ngrid,nlevs)
      REAL,INTENT(OUT) :: delta(ngrid,nlevs)
      REAL delta_2d(ip1jmp1,nlevs)
      INTEGER :: ll
      REAL ghx(ip1jmp1,nlevs), ghy(ip1jm,nlevs)


      CALL gr_fi_dyn(nlevs,ngrid,iip1,jjp1,temp,delta_2d)
      CALL grad(nlevs,delta_2d,ghx,ghy)
      DO ll=1,nlevs
          ghx(:,ll)=ghx(:,ll)*zmasqu
! pas de diffusion zonale  
!          ghx(:,ll)=0.
          ghy(:,ll)=ghy(:,ll)*zmasqv
      END DO
      CALL diverg(nlevs,ghx,ghy,delta_2d)
      CALL gr_dyn_fi(nlevs,iip1,jjp1,ngrid,delta_2d,delta)


  END SUBROUTINE divgrad_phy



      SUBROUTINE init_masquv(ngrid,zmasq)
     
      IMPLICIT NONE

#include "dimensions.h"
!#include "dimphys.h"
#include "paramet.h"
#include "comcstfi.h"   
#include "comgeom.h"
#include "comhdiff.h"


      INTEGER,INTENT(IN) :: ngrid
      REAL zmasq(ngrid)
      REAL zmasq_2d(ip1jmp1)
      REAL ff(ip1jm)
      REAL eps
      INTEGER i


! Masques u,v
      zmasqu=1.
      zmasqv=1.

      CALL gr_fi_dyn(1,ngrid,iip1,jjp1,zmasq,zmasq_2d)

      DO i=1,ip1jmp1-1
        IF (zmasq_2d(i).GT.1e-5 .OR. zmasq_2d(i+1).GT.1e-5) THEN
                zmasqu(i)=0.
        ENDIF
      END DO
      DO i=iip1,ip1jmp1,iip1
        zmasqu(i)=zmasqu(i-iim)
      END DO
      DO i=1,ip1jm
        IF (zmasq_2d(i).GT.1e-5 .OR. zmasq_2d(i+iip1).GT.1e-5) THEN  
                zmasqv(i)=0.
        END IF
      END DO


! Coriolis (pour Ekman transp.)
      eps=1e-5
!      CALL getin('slab_eps',eps)
!      print *,'epsilon=',eps
      ff=fext*unsairez       
      DO i=1,ip1jm
        unsev(i)=eps/(ff(i)*ff(i)+eps**2)
        unsfv(i)=ff(i)/(ff(i)*ff(i)+eps**2)
      ENDDO
      CALL gr_v_scal(1,unsfv,unsfu)
      CALL gr_v_scal(1,unsev,unseu)
! Alpha variable?
!      alpha_var=.FALSE.
!      CALL getin('slab_alphav',alpha_var)


      
  END SUBROUTINE init_masquv



      SUBROUTINE slab_ekman2(ngrid,tx_phy,ty_phy,ts_phy,dt_phy)
  
          use slab_ice_h 

      IMPLICIT NONE
      
#include "dimensions.h"
!#include "dimphys.h"
#include "paramet.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "comgeom.h"
#include "comhdiff.h"

      INTEGER,INTENT(IN) :: ngrid
      INTEGER ij
      REAL txv(ip1jm),fluxm(ip1jm),tyv(ip1jm)
      REAL fluxtm(ip1jm,noceanmx),fluxtz(ip1jmp1,noceanmx)
      REAL tyu(ip1jmp1),txu(ip1jmp1),fluxz(ip1jmp1),fluxv(ip1jmp1)
      REAL dt(ip1jmp1,noceanmx),ts(ip1jmp1,noceanmx)
      REAL tx_phy(ngrid),ty_phy(ngrid)
      REAL dt_phy(ngrid,noceanmx),ts_phy(ngrid,noceanmx)




! Passage taux,y sur grilles 2d
      CALL gr_fi_dyn(1,ngrid,iip1,jjp1,tx_phy,txu)
      CALL gr_fi_dyn(1,ngrid,iip1,jjp1,ty_phy,tyu)
! Passage grille u,v
! Multiplication par f ou eps.
      CALL gr_v_scal(1,txu,txv)
      CALL gr_v_scal(1,tyu,tyv)
      fluxm=tyv*unsev-txv*unsfv
!      fluxm=-txv*unsfv
      CALL gr_u_scal(1,txu,txu)
      CALL gr_u_scal(1,tyu,tyu)
      fluxz=tyu*unsfu+txu*unseu
!      fluxz=tyu*unsfu
            
! Calcul de T, Tdeep
      CALL gr_fi_dyn(2,ngrid,iip1,jjp1,ts_phy,ts)
       
! Flux de masse
      fluxm=fluxm*cv*cuvsurcv*zmasqv
      fluxz=fluxz*cu*cvusurcu*zmasqu
! Flux de masse vertical
      DO ij=iip2,ip1jm
        fluxv(ij)=fluxz(ij)-fluxz(ij-1)-fluxm(ij)+fluxm(ij-iip1)
      ENDDO
      DO ij=iip1,ip1jmi1,iip1
         fluxv(ij+1)=fluxv(ij+iip1)
      END DO
! Poles
      fluxv(1)=-SUM(fluxm(1:iim))     
      fluxv(ip1jmp1)=SUM(fluxm(ip1jm-iim:ip1jm-1))
      fluxv=fluxv

!Calcul flux de chaleur méridiens
      DO ij=1,ip1jm
          fluxtm(ij,1)=fluxm(ij)*(ts(ij+iip1,1)+ts(ij,1))/2.
          fluxtm(ij,2)=-fluxm(ij)*(ts(ij+iip1,2)+ts(ij,2))/2.
      END DO

!Calcul flux chaleur zonaux
      DO ij=iip2,ip1jm
      IF (fluxz(ij).GE.0.) THEN
             fluxtz(ij,1)=fluxz(ij)*ts(ij,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij+1,2)
      ELSE
             fluxtz(ij,1)=fluxz(ij)*ts(ij+1,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij,2)
      ENDIF
      END DO
      DO ij=iip1*2,ip1jmp1,iip1
             fluxtz(ij,:)=fluxtz(ij-iim,:)
      END DO
                   
! Calcul de dT
      DO ij=iip2,ip1jm
         dt(ij,:)=fluxtz(ij-1,:)-fluxtz(ij,:)   &
                  +fluxtm(ij,:)-fluxtm(ij-iip1,:)
         IF (fluxv(ij).GT.0.) THEN
           dt(ij,1)=dt(ij,1)+fluxv(ij)*ts(ij,2)
           dt(ij,2)=dt(ij,2)-fluxv(ij)*ts(ij,2)
         ELSE
           dt(ij,1)=dt(ij,1)+fluxv(ij)*ts(ij,1)
           dt(ij,2)=dt(ij,2)-fluxv(ij)*ts(ij,1)
         ENDIF
         dt(ij,:)=dt(ij,:)/aire(ij)
      END DO
      DO ij=iip1,ip1jmi1,iip1
         dt(ij+1,:)=dt(ij+iip1,:) 
      END DO
! Pôles
      dt(1,:)=SUM(fluxtm(1:iim,:),dim=1)
        IF (fluxv(1).GT.0.) THEN
          dt(1,1)=dt(1,1)+fluxv(1)*ts(1,2)
          dt(1,2)=dt(1,2)-fluxv(1)*ts(1,2)
        ELSE
          dt(1,1)=dt(1,1)+fluxv(1)*ts(1,1)
          dt(1,2)=dt(1,2)-fluxv(1)*ts(1,1)
        ENDIF
      dt(1,:)=dt(1,:)/apoln
      dt(ip1jmp1,:)=-SUM(fluxtm(ip1jm-iim:ip1jm-1,:),dim=1)
       IF (fluxv(ip1jmp1).GT.0.) THEN
         dt(ip1jmp1,1)=dt(ip1jmp1,1)+fluxv(ip1jmp1)*ts(ip1jmp1,2)
         dt(ip1jmp1,2)=dt(ip1jmp1,2)-fluxv(ip1jmp1)*ts(ip1jmp1,2)
       ELSE
         dt(ip1jmp1,1)=dt(ip1jmp1,1)+fluxv(ip1jmp1)*ts(ip1jmp1,1)
         dt(ip1jmp1,2)=dt(ip1jmp1,2)-fluxv(ip1jmp1)*ts(ip1jmp1,1)
       ENDIF
      dt(ip1jmp1,:)=dt(ip1jmp1,:)/apols
      dt(2:iip1,1)=dt(1,1)
      dt(2:iip1,2)=dt(1,2)
      dt(ip1jm+1:ip1jmp1,1)=dt(ip1jmp1,1)
      dt(ip1jm+1:ip1jmp1,2)=dt(ip1jmp1,2)

! Retour grille physique
      CALL gr_dyn_fi(2,iip1,jjp1,ngrid,dt,dt_phy)


  END SUBROUTINE slab_ekman2

 
END MODULE surf_heat_transp_mod



