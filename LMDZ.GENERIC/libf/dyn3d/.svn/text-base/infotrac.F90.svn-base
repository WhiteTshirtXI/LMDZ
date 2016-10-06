MODULE infotrac

IMPLICIT NONE
! nqtot : total number of tracers and higher order of moment, water vapor and liquid included
  INTEGER, SAVE :: nqtot
  INTEGER,allocatable :: iadv(:)   ! tracer advection scheme number
  CHARACTER(len=20),allocatable ::  tname(:) ! tracer name

CONTAINS

      subroutine iniadvtrac(nq,numvanle)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine which initializes tracer names and advection schemes
! reads these infos from file 'traceur.def' but uses default values
! if that file is not found.
! Ehouarn Millour. Oct. 2008  (made this LMDZ4-like) for future compatibility 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

!#include "dimensions.h"
!#include "advtrac.h"
!#include "control.h"

! routine arguments:
      INTEGER,INTENT(out) :: nq ! number of tracers
      INTEGER,INTENT(out) :: numvanle

! local variables:
      LOGICAL :: first
      INTEGER :: iq
      INTEGER :: ierr
      CHARACTER(len=3) :: qname

! Look for file traceur.def
      OPEN(90,file='traceur.def',form='formatted',status='old', &
              iostat=ierr)
      IF (ierr.eq.0) THEN
        write(*,*) "iniadvtrac: Reading file traceur.def"
        ! read number of tracers:
        read(90,*,iostat=ierr) nq
        if (ierr.ne.0) then
          write(*,*) "iniadvtrac: error reading number of tracers"
          write(*,*) "   (first line of traceur.def) "
          stop
        endif
        
        ! allocate arrays:
        allocate(iadv(nq))
        allocate(tname(nq))
        
        ! initialize advection schemes to Van-Leer for all tracers
        do iq=1,nq
          iadv(iq)=3 ! Van-Leer
        enddo
        
        do iq=1,nq
        ! minimal version, just read in the tracer names, 1 per line
          read(90,*,iostat=ierr) tname(iq)
          if (ierr.ne.0) then
            write(*,*) 'iniadvtrac: error reading tracer names...'
            stop
          endif
        enddo !of do iq=1,nq
        close(90) ! done reading tracer names, close file
      ENDIF ! of IF (ierr.eq.0)

!  ....  Choix  des shemas d'advection pour l'eau et les traceurs  ...
!  ...................................................................
!
!     iadv = 1    shema  transport type "humidite specifique LMD"  
!     iadv = 2    shema   amont
!     iadv = 3    shema  Van-leer
!     iadv = 4    schema  Van-leer + humidite specifique
!                        Modif F.Codron
! 
!
      DO  iq = 1, nq-1
       IF( iadv(iq).EQ.1 ) PRINT *,' Choix du shema humidite specifique'&
       ,' pour le traceur no ', iq
       IF( iadv(iq).EQ.2 ) PRINT *,' Choix du shema  amont',' pour le'  &
       ,' traceur no ', iq
       IF( iadv(iq).EQ.3 ) PRINT *,' Choix du shema  Van-Leer ',' pour' &
       ,'le traceur no ', iq

       IF( iadv(iq).EQ.4 )  THEN
         PRINT *,' Le shema  Van-Leer + humidite specifique ',          &
       ' est  uniquement pour la vapeur d eau .'
         PRINT *,' Corriger iadv( ',iq, ')  et repasser ! '
         CALL ABORT
       ENDIF

       IF( iadv(iq).LE.0.OR.iadv(iq).GT.4 )   THEN
        PRINT *,' Erreur dans le choix de iadv (nqtot).Corriger et '    &
       ,' repasser car  iadv(iq) = ', iadv(iq)
         CALL ABORT
       ENDIF
      ENDDO

!       IF( iadv(nq).EQ.1 ) PRINT *,' Choix du shema humidite '          &
!       ,'specifique pour la vapeur d''eau'
!       IF( iadv(nq).EQ.2 ) PRINT *,' Choix du shema  amont',' pour la'  &
!       ,' vapeur d''eau '
!       IF( iadv(nq).EQ.3 ) PRINT *,' Choix du shema  Van-Leer '         &
!       ,' pour la vapeur d''eau'
!       IF( iadv(nq).EQ.4 ) PRINT *,' Choix du shema  Van-Leer + '       &
!       ,' humidite specifique pour la vapeur d''eau'
!
!       IF( (iadv(nq).LE.0).OR.(iadv(nq).GT.4) )   THEN
!        PRINT *,' Erreur dans le choix de iadv (nqtot).Corriger et '    &
!       ,' repasser car  iadv(nqtot) = ', iadv(nqtot)
!         CALL ABORT
!       ENDIF

      first = .TRUE.
      numvanle = nq + 1
      DO  iq = 1, nq
        IF(((iadv(iq).EQ.3).OR.(iadv(iq).EQ.4)).AND.first ) THEN
          numvanle = iq
          first    = .FALSE. 
        ENDIF
      ENDDO
!
      DO  iq = 1, nq

      IF( (iadv(iq).NE.3.AND.iadv(iq).NE.4).AND.iq.GT.numvanle )  THEN
          PRINT *,' Il y a discontinuite dans le choix du shema de ',   &
          'Van-leer pour les traceurs . Corriger et repasser . '
           CALL ABORT
      ENDIF

      ENDDO
!
      end subroutine iniadvtrac

END MODULE infotrac
