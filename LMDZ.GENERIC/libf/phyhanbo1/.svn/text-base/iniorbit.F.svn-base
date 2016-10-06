      SUBROUTINE iniorbit
     $     (papoastr,pperiastr,pyear_day,pperi_day,pobliq)
      
      USE planete_mod, only: apoastr, periastr, year_day, obliquit,
     &                       peri_day, e_elips, p_elips, timeperi
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:
c   -------
c     Frederic Hourdin      22 Fevrier 1991
c
c   Objet:
c   ------
c    Initialisation du sous programme orbite qui calcule
c    a une date donnee de l'annee de duree year_day commencant
c    a l'equinoxe de printemps et dont le periastre se situe
c    a la date peri_day, la distance au soleil et la declinaison.
c
c   Interface:
c   ----------
c   - Doit etre appele avant d'utiliser orbite.
c   - initialise une partie du common planete_mod
c
c   Arguments:
c   ----------
c
c   Input:
c   ------
c   apoastr       \   apoastron and periastron of the orbit in AU
c   periastr      /   
c
c=======================================================================

c-----------------------------------------------------------------------
c   Declarations:
c   -------------

!#include "planete.h"
#include "comcstfi.h"

c   Arguments:
c   ----------

      REAL papoastr,pperiastr,pyear_day,pperi_day,pobliq

c   Local:
c   ------

      REAL zxref,zanom,zz,zx0,zdx
      INTEGER iter

c-----------------------------------------------------------------------

      pi=2.*asin(1.)

      apoastr =papoastr
      periastr=pperiastr
      year_day=pyear_day
      obliquit=pobliq
      peri_day=pperi_day


      !!!! SPARADRAP TEMPORAIRE !!!!
      !!!! SPARADRAP TEMPORAIRE !!!!
      !!!! We hope that all cases are above 25 Mkm [OK with Gliese 581d]
      IF ( apoastr .gt. 25.) THEN
        PRINT*,'!!!!! WARNING !!!!!'
        PRINT*,'!!!!! YOU ARE ABOUT TO WITNESS A DIRT HACK !!!!!'
        PRINT*,'This must be an old start file.'
        PRINT*,'The code changed 19/03/2012: we now use AU.'
        PRINT*,'So I am assuming units are in Mkm here'
        PRINT*,'and I am performing a conversion towards AU.'
        periastr = periastr / 149.598 ! Mkm to AU
        apoastr = apoastr / 149.598 ! Mkm to AU
      ENDIF
      !!!! SPARADRAP TEMPORAIRE !!!!
      !!!! SPARADRAP TEMPORAIRE !!!!

 
      PRINT*,'Periastron in AU  ',periastr
      PRINT*,'Apoastron in AU  ',apoastr 
      PRINT*,'Obliquity in degrees  :',obliquit


      e_elips=(apoastr-periastr)/(periastr+apoastr)
      p_elips=0.5*(periastr+apoastr)*(1-e_elips*e_elips)

      print*,'e_elips',e_elips
      print*,'p_elips',p_elips

c-----------------------------------------------------------------------
c calcul de l'angle polaire et de la distance au soleil :
c -------------------------------------------------------

c  calcul de l'zanomalie moyenne

      zz=(year_day-pperi_day)/year_day
      zanom=2.*pi*(zz-nint(zz))
      zxref=abs(zanom)
      PRINT*,'zanom  ',zanom

c  resolution de l'equation horaire  zx0 - e * sin (zx0) = zxref
c  methode de Newton

      zx0=zxref+e_elips*sin(zxref)
      DO 110 iter=1,100
         zdx=-(zx0-e_elips*sin(zx0)-zxref)/(1.-e_elips*cos(zx0))
         if(abs(zdx).le.(1.e-12)) goto 120
         zx0=zx0+zdx
110   continue
120   continue
      zx0=zx0+zdx
      if(zanom.lt.0.) zx0=-zx0
      PRINT*,'zx0   ',zx0

c zteta est la longitude solaire

      timeperi=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
      PRINT*,'Solar longitude of periastron timeperi = ',timeperi

      RETURN
      END
