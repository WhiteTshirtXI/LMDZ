! ADAPTATION GCM POUR CP(T)
c======================================================================
c S. Lebonnois, 10/2007:
c
c VENUS: Cp(T) = cpp*(T/T0)^nu 
c avec T0=460. et nu=0.35
c cpp=RCPD=cp0 = 1000.
c R/RCPD = RKAPPA
c
c La fonction d'Exner reste pk = RCPD*(play/pref)**RKAPPA
c 
c T et teta (temperature potentielle) sont liees par:
c 
c   integrale[teta a T](cp/T dT) = integrale[pref a p](R/p dp)
c
c Dans le cas de l'expression pour Venus, ca donne:
c
c   teta**nu = T**nu - nu * T0**nu * ln[ (p/pref)**RKAPPA ]
c ou
c   teta**nu = T**nu - nu * T0**nu * ln[pk/RCPD]
c
c On passe de T a teta par t2tpot(t,teta,pk)
c On passe de teta a T par tpot2t(teta,t,pk)
c
c Pour DT <-> Dteta, on utilise: dteta = dT *(T/teta)**(nu-1)
c -> routine dt2dtpot(dt,dteta,t,teta) 
c (utilisee seulement pour le contregradient)
c
c======================================================================

      FUNCTION cpdet(t)
      IMPLICIT none
#include "cpdet.h"

      real cpdet,t

      cpdet = cp0*(t/t0)**nu

      return
      end
      
c======================================================================
c======================================================================

      SUBROUTINE t2tpot(npoints,yt, yteta, ypk)
      IMPLICIT none
c======================================================================
c Arguments:
c
c yt   --------input-R- Temperature
c yteta-------output-R- Temperature potentielle
c ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
c
c======================================================================
#include "cpdet.h"

      integer npoints
      REAL    yt(npoints), yteta(npoints), ypk(npoints)
      
      yteta = yt**nu - nu * t0**nu * log(ypk/cp0)
      yteta = yteta**(1./nu)
        

      !yteta = yt/ypk


      return
      end

c======================================================================
c======================================================================

      SUBROUTINE tpot2t(npoints,yteta, yt, ypk)
      IMPLICIT none
c======================================================================
c Arguments:
c
c yteta--------input-R- Temperature potentielle
c yt   -------output-R- Temperature
c ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
c
c======================================================================
#include "cpdet.h"

      integer npoints
      REAL yt(npoints), yteta(npoints), ypk(npoints)

      yt = yteta**nu + nu * t0**nu * log(ypk/cp0)
      yt = yt**(1./nu)

      !yt = yteta*ypk

      return
      end

c$$$c======================================================================
c$$$c======================================================================
c$$$
c$$$      SUBROUTINE dt2dtpot(npoints,ydt, ydteta, yt, yteta, ypk, ydpk)
c$$$      IMPLICIT none
c$$$c======================================================================
c$$$c Arguments:
c$$$c
c$$$c ydt  --------input-R- D Temperature
c$$$c ydteta------output-R- D Temperature potentielle
c$$$c yt   --------input-R- Temperature
c$$$c yteta--------input-R- Temperature potentielle
c$$$c ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
c$$$c ydpk  -------input-R- D Fonction d'Exner
c$$$c
c$$$c======================================================================
c$$$#include "cpdet.h"
c$$$
c$$$      integer npoints
c$$$      REAL    ydt(npoints), ydteta(npoints), ydpk(npoints)
c$$$      REAL    yt(npoints), yteta(npoints), ypk(npoints)
c$$$      
c$$$      ydteta =  yteta**(1.-nu)* (ydt * yt**(nu-1.) 
c$$$     .                           - t0**nu * ydpk / ypk)
c$$$        
c$$$      return
c$$$      end

