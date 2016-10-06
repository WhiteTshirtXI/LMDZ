      SUBROUTINE defrun_new( tapedef, etatinit )
c
c-----------------------------------------------------------------------
c     Auteurs :   L. Fairhead , P. Le Van  .
c      Modif C. Hourdin F. Forget VERSION MARTIENNE
c
c
c  -------------------------------------------------------------------
c
c                    MODIF JUIN 2000 (zoom)
c       .........     Version  du 29/04/97       ..........
c
c   Nouveaux parametres nitergdiv,nitergrot,niterh,tetagdiv,tetagrot,
c   tetatemp   ajoutes  pour la dissipation   .
c
c   Autre parametre ajoute en fin de liste : ** fxyhypb ** 
c
c   Si fxyhypb = .TRUE. , choix de la fonction a derivee tangente hyperb.
c   Sinon , choix de fxynew  , a derivee sinusoidale  ..
c
c   ......  etatinit = . TRUE. si defrun_new  est appele dans NEWSTART
c   ETAT0_LMD  ou  LIMIT_LMD  pour l'initialisation de start.dat (dic) et
c   de limit.dat (dic)  ...........
c   Sinon  etatinit = . FALSE .
c
c   Donc etatinit = .F.  si on veut comparer les valeurs de  alphax ,
c   alphay,clon,clat, fxyhypb  lues sur  le fichier  start  avec
c   celles passees  par run.def ,  au debut du gcm, apres l'appel a 
c   lectba .  
c   Ces parametres definissant entre autres la grille et doivent etre
c   pareils et coherents , sinon il y aura  divergence du gcm .
c
c
c-----------------------------------------------------------------------
c   Declarations :
c   --------------
! to use  'getin'
      USE ioipsl_getincom
      use sponge_mod,only: callsponge,nsponge,mode_sponge,tetasponge
      use control_mod,only: nday, day_step, iperiod, anneeref,
     &                      iconser, idissip, iphysiq, ecritphy
      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"
!#include "control.h"
#include "logic.h"
#include "serre.h"
#include "comdissnew.h"
c
c   arguments:
c   ---------
      LOGICAL  etatinit ! should be .false. for a call from gcm.F
                        ! and .true. for a call from newstart.F
      INTEGER  tapedef  ! unit number to assign to 'run.def' file
c
c   local variables:
c   ---------------

      CHARACTER ch1*72,ch2*72,ch3*72,ch4*8 ! to store various strings
      INTEGER tapeout ! unit numbers for (standard) outputs
      parameter (tapeout=6)
      integer tapeerr ! unit number for error message
      parameter (tapeerr=0)

c     REAL clonn,clatt,alphaxx,alphayy
c     LOGICAL  fxyhypbb
      INTEGER ierr
      REAL clonn,clatt,grossismxx,grossismyy
      REAL dzoomxx,dzoomyy,tauxx,tauyy,temp
      LOGICAL  fxyhypbb, ysinuss


c   initialisations:
c   ----------------
 
!      tapeout=6

c-----------------------------------------------------------------------
c  Parametres de controle du run:
c-----------------------------------------------------------------------


!Initialisation des parametres "terrestres", qui ne concernent pas
!le modele martien et ne sont donc plus lues dans "run.def"

        anneeref=0
        ! Note: anneref is a common in 'control.h'

      OPEN(tapedef,file='run.def',status='old',form='formatted'
     .     ,iostat=ierr)
      CLOSE(tapedef) ! first call to getin will open the file

      IF(ierr.EQ.0) THEN ! if file run.def is found
        WRITE(tapeout,*) "DEFRUN_NEW: reading file run.def"
        
        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "Number of days to run:"
        nday=1 ! default value
        call getin("nday",nday)
        WRITE(tapeout,*)" nday = ",nday

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "Number of dynamical steps per day:",
     & "(should be a multiple of iperiod)"
        day_step=960 ! default value
        call getin("day_step",day_step)
        WRITE(tapeout,*)" day_step = ",day_step

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "periode pour le pas Matsuno (en pas)"
        iperiod=5 ! default value
        call getin("iperiod",iperiod)
        WRITE(tapeout,*)" iperiod = ",iperiod

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "periode de sortie des variables de ",
     &  "controle (en pas)"
        iconser=120 ! default value
        call getin("iconser",iconser)
        WRITE(tapeout,*)" iconser = ",iconser

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "periode de la dissipation (en pas)"
        idissip=5 ! default value
        call getin("idissip",idissip)
        WRITE(tapeout,*)" idissip = ",idissip

ccc  ....   P. Le Van , modif le 29/04/97 .pour la dissipation  ...
ccc
        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "choix de l'operateur de dissipation ",
     & "(star ou  non star )"
        lstardis=.true. ! default value
        call getin("lstardis",lstardis)
        WRITE(tapeout,*)" lstardis = ",lstardis

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "avec ou sans coordonnee hybrides"
        hybrid=.true. ! default value
        call getin("hybrid",hybrid)
        WRITE(tapeout,*)" hybrid = ",hybrid

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "nombre d'iterations de l'operateur de ",
     & "dissipation   gradiv "
        nitergdiv=1 ! default value
        call getin("nitergdiv",nitergdiv)
        WRITE(tapeout,*)" nitergdiv = ",nitergdiv

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "nombre d'iterations de l'operateur de ",
     & "dissipation  nxgradrot"
        nitergrot=2 ! default value
        call getin("nitergrot",nitergrot)
        WRITE(tapeout,*)" nitergrot = ",nitergrot

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "nombre d'iterations de l'operateur de ",
     & "dissipation  divgrad"
        niterh=2 ! default value
        call getin("niterh",niterh)
        WRITE(tapeout,*)" niterh = ",niterh

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "temps de dissipation des plus petites ",
     & "long.d ondes pour u,v (gradiv)"
        tetagdiv=4000. ! default value
        call getin("tetagdiv",tetagdiv)
        WRITE(tapeout,*)" tetagdiv = ",tetagdiv

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "temps de dissipation des plus petites ",
     & "long.d ondes pour u,v(nxgradrot)"
        tetagrot=5000. ! default value
        call getin("tetagrot",tetagrot)
        WRITE(tapeout,*)" tetagrot = ",tetagrot

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "temps de dissipation des plus petites ",
     & "long.d ondes pour  h ( divgrad)"
        tetatemp=5000. ! default value
        call getin("tetatemp",tetatemp)
        WRITE(tapeout,*)" tetatemp = ",tetatemp

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "coefficient pour gamdissip"
        coefdis=0. ! default value
        call getin("coefdis",coefdis)
        WRITE(tapeout,*)" coefdis = ",coefdis
c
c    ...............................................................

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "choix du shema d'integration temporelle ",
     & "(true==Matsuno ou false==Matsuno-leapfrog)"
        purmats=.false. ! default value
        call getin("purmats",purmats)
        WRITE(tapeout,*)" purmats = ",purmats

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "avec ou sans physique"
        physic=.true. ! default value
        call getin("physic",physic)
        WRITE(tapeout,*)" physic = ",physic

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "periode de la physique (en pas)"
        iphysiq=20 ! default value
        call getin("iphysiq",iphysiq)
        WRITE(tapeout,*)" iphysiq = ",iphysiq

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "choix d'une grille reguliere"
        grireg=.true.
        call getin("grireg",grireg)
        WRITE(tapeout,*)" grireg = ",grireg

ccc   .... P.Le Van, ajout le 03/01/96 pour l'ecriture phys ...
c
        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "frequence (en pas) de l'ecriture ",
     & "du fichier diagfi.nc"
        ecritphy=240
        call getin("ecritphy",ecritphy)
        WRITE(tapeout,*)" ecritphy = ",ecritphy

ccc  ....   P. Le Van , ajout  le 7/03/95 .pour le zoom ...
c     .........   (  modif  le 17/04/96 )   .........
c
        if (.not.etatinit ) then 

           clonn=63.
           call getin("clon",clonn)
           
           IF( ABS(clon - clonn).GE. 0.001 )  THEN
             PRINT *,' La valeur de clon passee par run.def est '
     *       ,'differente de celle lue sur le fichier start '
             STOP
           ENDIF
c
c
           clatt=0.
           call getin("clat",clatt)
  
           IF( ABS(clat - clatt).GE. 0.001 )  THEN
             PRINT *,' La valeur de clat passee par run.def est '
     *       ,'differente de celle lue sur le fichier start '
             STOP
           ENDIF

           grossismxx=1.0
           call getin("grossismx",grossismxx)

           if(grossismxx.eq.0) then  
             write(*,*)
             write(*,*)'ERREUR : dans run.def, grossismx =0'
             write(*,*)'Attention a ne pas utiliser une version de'
             write(*,*)'run.def avant le nouveau zoom LMDZ2.3 (06/2000)'
             write(*,*)'(Il faut ajouter grossismx,dzoomx,etc... a la'
             write(*,*)'place de alphax, alphay. cf. dyn3d). '
             write(*,*)
             stop
           end if

           IF( ABS(grossismx - grossismxx).GE. 0.001 )  THEN
             PRINT *,' La valeur de grossismx passee par run.def est '
     *       ,'differente de celle lue sur le fichier  start =',
     *        grossismx
             if (grossismx.eq.0) then
                  write(*,*) 'OK,Normal : c est un vieux start'
     *             , 'd avant le nouveau zoom LMDZ2.3 (06/2000)'
                 grossismx=grossismxx
             else
                   STOP
             endif
           ENDIF

           grossismyy=1.0
           call getin("grossismy",grossismyy)

           IF( ABS(grossismy - grossismyy).GE. 0.001 )  THEN
             PRINT *,' La valeur de grossismy passee par run.def est '
     *       ,'differente de celle lue sur le fichier  start =',
     *        grossismy
             if (grossismy.eq.0) then
                  write(*,*) 'OK,Normal : c est un vieux start'
     *             , 'd avant le nouveau zoom LMDZ2.3 (06/2000)'
                 grossismy=grossismyy
             else
                   STOP
             endif
           ENDIF


           IF( grossismx.LT.1. )  THEN
             PRINT *,' ***  ATTENTION !! grossismx < 1 .   *** '
             STOP
           ELSE
             alphax = 1. - 1./ grossismx
           ENDIF

           IF( grossismy.LT.1. )  THEN
             PRINT *,' ***  ATTENTION !! grossismy < 1 .   *** '
             STOP
           ELSE
             alphay = 1. - 1./ grossismy
           ENDIF

           PRINT *,' '
           PRINT *,' --> In defrun: alphax alphay  ',alphax,alphay
           PRINT *,' '
c
           fxyhypbb=.false.
           call getin("fxyhypbb",fxyhypbb)
  
           IF( .NOT.fxyhypb )  THEN
             IF( fxyhypbb )     THEN
                PRINT *,' ********  PBS DANS  DEFRUN  ******** '
                PRINT *,' *** fxyhypb lu sur le fichier start est F ',
     *          'alors  qu il est  T  sur  run.def  ***'
                STOP
             ENDIF
           ELSE
             IF( .NOT.fxyhypbb )   THEN
                PRINT *,' ********  PBS DANS  DEFRUN  ******** '
                PRINT *,' ***  fxyhypb lu sur le fichier start est T ',
     *         'alors  qu il est  F  sur  run.def  ****  '
                STOP
             ENDIF
           ENDIF
           dzoomxx=0.0
           call getin("dzoomx",dzoomxx)

           IF( fxyhypb )  THEN
              IF( ABS(dzoomx - dzoomxx).GE. 0.001 )  THEN
                PRINT *,' La valeur de dzoomx passee par run.def est '
     *          ,'differente de celle lue sur le fichier  start '
                STOP
              ENDIF
           ENDIF

           dzoomyy=0.0
           call getin("dzoomy",dzoomyy)

           IF( fxyhypb )  THEN
              IF( ABS(dzoomy - dzoomyy).GE. 0.001 )  THEN
                PRINT *,' La valeur de dzoomy passee par run.def est '
     *          ,'differente de celle lue sur le fichier  start '
                STOP
              ENDIF
           ENDIF

           tauxx=2.0
           call getin("taux",tauxx)

           tauyy=2.0
           call getin("tauy",tauyy)

           IF( fxyhypb )  THEN
              IF( ABS(taux - tauxx).GE. 0.001 )  THEN
                WRITE(6,*)' La valeur de taux passee par run.def est', 
     *             'differente de celle lue sur le fichier  start '
                CALL ABORT
              ENDIF

              IF( ABS(tauy - tauyy).GE. 0.001 )  THEN
                WRITE(6,*)' La valeur de tauy passee par run.def est',
     *          'differente de celle lue sur le fichier  start '
                CALL ABORT
              ENDIF
           ENDIF
  
        ELSE    ! Below, case when etainit=.true.

           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "longitude en degres du centre du zoom"
           clon=63. ! default value
           call getin("clon",clon)
           WRITE(tapeout,*)" clon = ",clon
           
c
           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "latitude en degres du centre du zoom "
           clat=0. ! default value
           call getin("clat",clat)
           WRITE(tapeout,*)" clat = ",clat

           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "facteur de grossissement du zoom,",
     & " selon longitude"
           grossismx=1.0 ! default value
           call getin("grossismx",grossismx)
           WRITE(tapeout,*)" grossismx = ",grossismx

           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "facteur de grossissement du zoom ,",
     & " selon latitude"
           grossismy=1.0 ! default value
           call getin("grossismy",grossismy)
           WRITE(tapeout,*)" grossismy = ",grossismy

           IF( grossismx.LT.1. )  THEN
            PRINT *,' ***  ATTENTION !! grossismx < 1 .   *** '
            STOP
           ELSE
             alphax = 1. - 1./ grossismx
           ENDIF

           IF( grossismy.LT.1. )  THEN
             PRINT *,' ***  ATTENTION !! grossismy < 1 .   *** '
             STOP
           ELSE
             alphay = 1. - 1./ grossismy
           ENDIF

           PRINT *,' Defrun  alphax alphay  ',alphax,alphay
c
           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "Fonction  f(y)  hyperbolique  si = .true.",
     &  ", sinon  sinusoidale"
           fxyhypb=.false. ! default value
           call getin("fxyhypb",fxyhypb)
           WRITE(tapeout,*)" fxyhypb = ",fxyhypb

           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "extension en longitude de la zone du zoom", 
     & " (fraction de la zone totale)"
           dzoomx=0. ! default value
           call getin("dzoomx",dzoomx)
           WRITE(tapeout,*)" dzoomx = ",dzoomx

           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "extension en latitude de la zone du zoom", 
     & " (fraction de la zone totale)"
           dzoomy=0. ! default value
           call getin("dzoomy",dzoomy)
           WRITE(tapeout,*)" dzoomy = ",dzoomy

           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "raideur du zoom en  X"
           taux=2. ! default value
           call getin("taux",taux)
           WRITE(tapeout,*)" taux = ",taux

           WRITE(tapeout,*) ""
           WRITE(tapeout,*) "raideur du zoom en  Y"
           tauy=2.0 ! default value
           call getin("tauy",tauy)
           WRITE(tapeout,*)" tauy = ",tauy

        END IF ! of if (.not.etatinit )

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "Use a sponge layer?"
        callsponge=.true. ! default value
        call getin("callsponge",callsponge)
        WRITE(tapeout,*)" callsponge = ",callsponge

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "Sponge: number of layers over which",
     &                    " sponge extends"
        nsponge=3 ! default value
        call getin("nsponge",nsponge)
        WRITE(tapeout,*)" nsponge = ",nsponge

        WRITE(tapeout,*)""
        WRITE(tapeout,*)"Sponge mode: (forcing is towards ..."
        WRITE(tapeout,*)"  over upper nsponge layers)"
        WRITE(tapeout,*)"  0: (h=hmean,u=v=0)"
        WRITE(tapeout,*)"  1: (h=hmean,u=umean,v=0)"
        WRITE(tapeout,*)"  2: (h=hmean,u=umean,v=vmean)"
        mode_sponge=2 ! default value
        call getin("mode_sponge",mode_sponge)
        WRITE(tapeout,*)" mode_sponge = ",mode_sponge

        WRITE(tapeout,*) ""
        WRITE(tapeout,*) "Sponge: characteristic time scale tetasponge"
        WRITE(tapeout,*) "(seconds) at topmost layer (time scale then "
        WRITE(tapeout,*) " doubles with decreasing layer index)."
        tetasponge=50000.0
        call getin("tetasponge",tetasponge)
        WRITE(tapeout,*)" tetasponge = ",tetasponge


      WRITE(tapeout,*) '-----------------------------------------------'
      WRITE(tapeout,*) ' '
      WRITE(tapeout,*) ' '
c

c       Unlike on Earth (cf LMDZ2.2) , always a regular grid on Mars :
        ysinus = .false. !Mars Mettre a jour


      WRITE(tapeout,*) '-----------------------------------------------'
      WRITE(tapeout,*) ' '
      WRITE(tapeout,*) ' '
cc
      ELSE
        write(tapeerr,*) ' WHERE IS run.def ? WE NEED IT !!!!!!!!!!!!!!'
        stop
      ENDIF ! of IF(ierr.eq.0)

c     Test sur le zoom

      if((grossismx.eq.1).and.(grossismy.eq.1)) then  
c        Pas de zoom :
         write(tapeout,*) 'No zoom ? -> fxyhypb set to False'
     &   ,'           (It runs better that way)'
         fxyhypb = .false.
      else     
c        Avec Zoom
         if (.not.fxyhypb) stop 'With zoom, fxyhypb should be set to T 
     &in run.def for this version... -> STOP ! '     
      end if

      RETURN
      END
