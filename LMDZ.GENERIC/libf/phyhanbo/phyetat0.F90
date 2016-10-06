subroutine phyetat0 (ngrid,nlayer,fichnom,tab0,Lmodif,nsoil,nq, &
                     day_ini,time,tsurf,tsoil, &
                     emis,q2,qsurf,cloudfrac,totcloudfrac,hice, &
                     rnat,pctsrf_sic,tslab,tsea_ice,sea_ice)


  USE infotrac, ONLY: tname
  USE surfdat_h, only: phisfi, albedodat, zmea, zstd, zsig, zgam, zthe
  use iostart, only: nid_start, open_startphy, close_startphy, &
                     get_field, get_var, inquire_field, &
                     inquire_dimension, inquire_dimension_length
  use slab_ice_h, only: noceanmx

  implicit none

!======================================================================
! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
!  Adaptation à Mars : Yann Wanherdrick 
! Objet: Lecture de l etat initial pour la physique
!======================================================================
#include "netcdf.inc"
!#include "dimensions.h"
!#include "dimphys.h"
!#include "planete.h"
#include "comcstfi.h"

!======================================================================
!  INTEGER nbsrf !Mars nbsrf a 1 au lieu de 4
!  PARAMETER (nbsrf=1) ! nombre de sous-fractions pour une maille
!======================================================================
!  Arguments:
!  ---------
!  inputs:
  integer,intent(in) :: ngrid
  integer,intent(in) :: nlayer
  character*(*),intent(in) :: fichnom ! "startfi.nc" file
  integer,intent(in) :: tab0
  integer,intent(in) :: Lmodif
  integer,intent(in) :: nsoil ! # of soil layers
  integer,intent(in) :: nq
  integer,intent(in) :: day_ini
  real,intent(in) :: time

!  outputs:
  real,intent(out) :: tsurf(ngrid) ! surface temperature
  real,intent(out) :: tsoil(ngrid,nsoil) ! soil temperature
  real,intent(out) :: emis(ngrid) ! surface emissivity
  real,intent(out) :: q2(ngrid,nlayer+1) ! 
  real,intent(out) :: qsurf(ngrid,nq) ! tracers on surface
! real co2ice(ngrid) ! co2 ice cover
  real,intent(out) :: cloudfrac(ngrid,nlayer)
  real,intent(out) :: hice(ngrid), totcloudfrac(ngrid)
  real,intent(out) :: pctsrf_sic(ngrid),tslab(ngrid,noceanmx)  
  real,intent(out) :: tsea_ice(ngrid),sea_ice(ngrid)
  real,intent(out) :: rnat(ngrid)

!======================================================================
!  Local variables:

!      INTEGER radpas
!      REAL co2_ppm
!      REAL solaire

      real xmin,xmax ! to display min and max of a field
!
      INTEGER ig,iq,lmax
      INTEGER nid, nvarid
      INTEGER ierr, i, nsrf
!      integer isoil 
!      INTEGER length
!      PARAMETER (length=100)
      CHARACTER*7 str7
      CHARACTER*2 str2
      CHARACTER*1 yes
!
      REAL p_rad,p_omeg,p_g,p_cpp,p_mugaz,p_daysec
      INTEGER nqold

! flag which identifies if 'startfi.nc' file is using old names (qsurf01,...)
!      logical :: oldtracernames=.false.
      integer :: count
      character(len=30) :: txt ! to store some text
      
      INTEGER :: indextime=1 ! index of selected time, default value=1
      logical :: found

!
! ALLOCATE ARRAYS IN surfdat_h
!
IF (.not. ALLOCATED(albedodat)) ALLOCATE(albedodat(ngrid))
IF (.not. ALLOCATED(phisfi)) ALLOCATE(phisfi(ngrid))
IF (.not. ALLOCATED(zmea)) ALLOCATE(zmea(ngrid))
IF (.not. ALLOCATED(zstd)) ALLOCATE(zstd(ngrid))
IF (.not. ALLOCATED(zsig)) ALLOCATE(zsig(ngrid))
IF (.not. ALLOCATED(zgam)) ALLOCATE(zgam(ngrid))
IF (.not. ALLOCATED(zthe)) ALLOCATE(zthe(ngrid))


! open physics initial state file:
call open_startphy(fichnom)


! possibility to modify tab_cntrl in tabfi
write(*,*)
write(*,*) 'TABFI in phyeta0: Lmodif=',Lmodif," tab0=",tab0
call tabfi (ngrid,nid_start,Lmodif,tab0,day_ini,lmax,p_rad, &
                   p_omeg,p_g,p_cpp,p_mugaz,p_daysec,time)

!c
!c Lecture des latitudes (coordonnees):
!c
!      ierr = NF_INQ_VARID (nid, "latitude", nvarid)
!      IF (ierr.NE.NF_NOERR) THEN
!         PRINT*, 'phyetat0: Le champ <latitude> est absent'
!         CALL abort
!      ENDIF
!#ifdef NC_DOUBLE
!      ierr = NF_GET_VARA_DOUBLE(nid,nvarid,sta,ngrid,lati)
!#else
!      ierr = NF_GET_VARA_REAL(nid,nvarid,sta,ngrid,lati)
!#endif
!      IF (ierr.NE.NF_NOERR) THEN
!         PRINT*, 'phyetat0: Lecture echouee pour <latitude>'
!         CALL abort
!      ENDIF
!c
!c Lecture des longitudes (coordonnees):
!c
!      ierr = NF_INQ_VARID (nid, "longitude", nvarid)
!      IF (ierr.NE.NF_NOERR) THEN
!         PRINT*, 'phyetat0: Le champ <longitude> est absent'
!         CALL abort
!      ENDIF
!#ifdef NC_DOUBLE
!      ierr = NF_GET_VARA_DOUBLE(nid,nvarid,sta,ngrid,long)
!#else
!      ierr = NF_GET_VARA_REAL(nid,nvarid,sta,ngrid,long)
!#endif
!      IF (ierr.NE.NF_NOERR) THEN
!         PRINT*, 'phyetat0: Lecture echouee pour <longitude>'
!         CALL abort
!      ENDIF
!c
!c Lecture des aires des mailles:
!c
!      ierr = NF_INQ_VARID (nid, "area", nvarid)
!      IF (ierr.NE.NF_NOERR) THEN
!         PRINT*, 'phyetat0: Le champ <area> est absent'
!         CALL abort
!      ENDIF
!#ifdef NC_DOUBLE
!      ierr = NF_GET_VARA_DOUBLE(nid,nvarid,sta,ngrid,area)
!#else
!      ierr = NF_GET_VARA_REAL(nid,nvarid,sta,ngrid,area)
!#endif
!      IF (ierr.NE.NF_NOERR) THEN
!         PRINT*, 'phyetat0: Lecture echouee pour <area>'
!         CALL abort
!      ENDIF
!      xmin = 1.0E+20
!      xmax = -1.0E+20
!      xmin = MINVAL(area)
!      xmax = MAXVAL(area)
!      PRINT*,'Aires des mailles <area>:', xmin, xmax

! Load surface geopotential:
call get_field("phisfi",phisfi,found)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <phisfi>"
  call abort
else
  write(*,*) "phyetat0: surface geopotential <phisfi> range:", &
             minval(phisfi), maxval(phisfi)
endif

! Load bare ground albedo:
call get_field("albedodat",albedodat,found)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <albedodat>"
  call abort
else
  write(*,*) "phyetat0: Bare ground albedo <albedodat> range:", &
             minval(albedodat), maxval(albedodat)
endif

! ZMEA
call get_field("ZMEA",zmea,found)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <ZMEA>"
  call abort
else
  write(*,*) "phyetat0: <ZMEA> range:", &
             minval(zmea), maxval(zmea)
endif

! ZSTD
call get_field("ZSTD",zstd,found)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <ZSTD>"
  call abort
else
  write(*,*) "phyetat0: <ZSTD> range:", &
             minval(zstd), maxval(zstd)
endif

! ZSIG
call get_field("ZSIG",zsig,found)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <ZSIG>"
  call abort
else
  write(*,*) "phyetat0: <ZSIG> range:", &
             minval(zsig), maxval(zsig)
endif

! ZGAM
call get_field("ZGAM",zgam,found)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <ZGAM>"
  call abort
else
  write(*,*) "phyetat0: <ZGAM> range:", &
             minval(zgam), maxval(zgam)
endif

! ZTHE
call get_field("ZTHE",zthe,found)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <ZTHE>"
  call abort
else
  write(*,*) "phyetat0: <ZTHE> range:", &
             minval(zthe), maxval(zthe)
endif

! Surface temperature :
call get_field("tsurf",tsurf,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <tsurf>"
  call abort
else
  write(*,*) "phyetat0: Surface temperature <tsurf> range:", &
             minval(tsurf), maxval(tsurf)
endif

! Surface emissivity
call get_field("emis",emis,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <emis>"
  call abort
else
  write(*,*) "phyetat0: Surface emissivity <emis> range:", &
             minval(emis), maxval(emis)
endif

! Cloud fraction (added by BC 2010)
call get_field("cloudfrac",cloudfrac,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <cloudfrac>"
  call abort
else
  write(*,*) "phyetat0: Cloud fraction <cloudfrac> range:", &
             minval(cloudfrac), maxval(cloudfrac)
endif

! Total cloud fraction (added by BC 2010)
call get_field("totcloudfrac",totcloudfrac,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <totcloudfrac>"
  call abort
else
  write(*,*) "phyetat0: Total cloud fraction <totcloudfrac> range:", &
             minval(totcloudfrac), maxval(totcloudfrac)
endif

! Height of oceanic ice (added by BC 2010)
call get_field("hice",hice,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <hice>"
!  call abort
      do ig=1,ngrid
      hice(ig)=0.
      enddo
else
  write(*,*) "phyetat0: Height of oceanic ice <hice> range:", &
             minval(hice), maxval(hice)
endif

! SLAB OCEAN (added by BC 2014)
! nature of the surface
call get_field("rnat",rnat,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <rnat>"
      do ig=1,ngrid
        rnat(ig)=1.
      enddo
else
      do ig=1,ngrid
        if((nint(rnat(ig)).eq.2).or.(nint(rnat(ig)).eq.0))then
          rnat(ig)=0.
        else
          rnat(ig)=1.
        endif      
      enddo

  write(*,*) "phyetat0: Nature of surface <rnat> range:", &
             minval(rnat), maxval(rnat)
endif
! Pourcentage of sea ice cover
call get_field("pctsrf_sic",pctsrf_sic,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <pctsrf_sic>"
      do ig=1,ngrid
      pctsrf_sic(ig)=0.
      enddo
else
  write(*,*) "phyetat0: Pourcentage of sea ice cover <pctsrf_sic> range:", &
             minval(pctsrf_sic), maxval(pctsrf_sic)
endif
! Slab ocean temperature (2 layers)
call get_field("tslab",tslab,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <tslab>"
      do ig=1,ngrid
      do iq=1,noceanmx
      tslab(ig,iq)=tsurf(ig)
      enddo
      enddo
else
  write(*,*) "phyetat0: Slab ocean temperature <tslab> range:", &
             minval(tslab), maxval(tslab)
endif
! Oceanic ice temperature
call get_field("tsea_ice",tsea_ice,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <tsea_ice>"
      do ig=1,ngrid
      tsea_ice(ig)=273.15-1.8
      enddo
else
  write(*,*) "phyetat0: Oceanic ice temperature <tsea_ice> range:", &
             minval(tsea_ice), maxval(tsea_ice)
endif
!  Oceanic ice quantity (kg/m^2)
call get_field("sea_ice",sea_ice,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <sea_ice>"
      do ig=1,ngrid
      tsea_ice(ig)=0.
      enddo
else
  write(*,*) "phyetat0: Oceanic ice quantity <sea_ice> range:", &
             minval(sea_ice), maxval(sea_ice)
endif




! pbl wind variance
call get_field("q2",q2,found,indextime)
if (.not.found) then
  write(*,*) "phyetat0: Failed loading <q2>"
  call abort
else
  write(*,*) "phyetat0: PBL wind variance <q2> range:", &
             minval(q2), maxval(q2)
endif

! tracer on surface
if (nq.ge.1) then
  do iq=1,nq
    txt=tname(iq)
    if (txt.eq."h2o_vap") then
      ! There is no surface tracer for h2o_vap;
      ! "h2o_ice" should be loaded instead
      txt="h2o_ice"
      write(*,*) 'phyetat0: loading surface tracer', &
                           ' h2o_ice instead of h2o_vap'
    endif
    call get_field(txt,qsurf(:,iq),found,indextime)
    if (.not.found) then
      write(*,*) "phyetat0: Failed loading <",trim(txt),">"
      write(*,*) "         ",trim(txt)," is set to zero"
    else
      write(*,*) "phyetat0: Surface tracer <",trim(txt),"> range:", &
                 minval(qsurf(:,iq)), maxval(qsurf(:,iq))
    endif
  enddo
endif ! of if (nq.ge.1)


! Call to soil_settings, in order to read soil temperatures,
! as well as thermal inertia and volumetric heat capacity
call soil_settings(nid_start,ngrid,nsoil,tsurf,tsoil,indextime)
!
! close file:
!
call close_startphy

END SUBROUTINE phyetat0
