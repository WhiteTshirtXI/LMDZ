subroutine iniphysiq(ngrid,nlayer, punjours, pdayref,ptimestep,          &
                     plat,plon,parea,pcu,pcv,                            &
                     prad,pg,pr,pcpp,iflag_phys)

use dimphy, only : klev ! number of atmospheric levels
use mod_grid_phy_lmdz, only : klon_glo ! number of atmospheric columns
                                       ! (on full grid)
use mod_phys_lmdz_para, only : klon_omp, & ! number of columns (on local omp grid)
                               klon_omp_begin, & ! start index of local omp subgrid
                               klon_omp_end, & ! end index of local omp subgrid
                               klon_mpi_begin ! start indes of columns (on local mpi grid)
use comgeomphy, only : airephy, & ! physics grid area (m2)
                       cuphy, & ! cu coeff. (u_covariant = cu * u)
                       cvphy, & ! cv coeff. (v_covariant = cv * v)
                       rlond, & ! longitudes
                       rlatd ! latitudes
use infotrac, only : nqtot ! number of advected tracers
use planete_mod, only: ini_planete_mod

implicit none
include "dimensions.h"
include "comvert.h"

real,intent(in) :: prad ! radius of the planet (m)
real,intent(in) :: pg ! gravitational acceleration (m/s2)
real,intent(in) :: pr ! ! reduced gas constant R/mu
real,intent(in) :: pcpp ! specific heat Cp
real,intent(in) :: punjours ! length (in s) of a standard day
integer,intent(in) :: ngrid ! number of horizontal grid points in the physics (full grid)
integer,intent(in) :: nlayer ! number of atmospheric layers
real,intent(in) :: plat(ngrid) ! latitudes of the physics grid
real,intent(in) :: plon(ngrid) ! longitudes of the physics grid
real,intent(in) :: parea(klon_glo) ! area (m2)
real,intent(in) :: pcu(klon_glo) ! cu coeff. (u_covariant = cu * u)
real,intent(in) :: pcv(klon_glo) ! cv coeff. (v_covariant = cv * v)
integer,intent(in) :: pdayref ! reference day of for the simulation
real,intent(in) :: ptimestep !physics time step (s)
integer,intent(in) :: iflag_phys ! type of physics to be called

integer :: ibegin,iend,offset
character(len=20) :: modname='iniphysiq'
character(len=80) :: abort_message

IF (nlayer.NE.klev) THEN
  write(*,*) 'STOP in ',trim(modname)
  write(*,*) 'Problem with dimensions :'
  write(*,*) 'nlayer     = ',nlayer
  write(*,*) 'klev   = ',klev
  abort_message = ''
  CALL abort_gcm (modname,abort_message,1)
ENDIF

IF (ngrid.NE.klon_glo) THEN
  write(*,*) 'STOP in ',trim(modname)
  write(*,*) 'Problem with dimensions :'
  write(*,*) 'ngrid     = ',ngrid
  write(*,*) 'klon   = ',klon_glo
  abort_message = ''
  CALL abort_gcm (modname,abort_message,1)
ENDIF

!$OMP PARALLEL PRIVATE(ibegin,iend) & 
	!$OMP SHARED(parea,pcu,pcv,plon,plat)
      
offset=klon_mpi_begin-1
airephy(1:klon_omp)=parea(offset+klon_omp_begin:offset+klon_omp_end)
cuphy(1:klon_omp)=pcu(offset+klon_omp_begin:offset+klon_omp_end)
cvphy(1:klon_omp)=pcv(offset+klon_omp_begin:offset+klon_omp_end)
rlond(1:klon_omp)=plon(offset+klon_omp_begin:offset+klon_omp_end)
rlatd(1:klon_omp)=plat(offset+klon_omp_begin:offset+klon_omp_end)

! copy over preff , ap() and bp() 
call ini_planete_mod(nlayer,preff,ap,bp)

! copy some fundamental parameters to physics 
! and do some initializations 
call inifis(klon_omp,nlayer,nqtot,pdayref,punjours,ptimestep, &
            rlatd,rlond,airephy,prad,pg,pr,pcpp)

!$OMP END PARALLEL


end subroutine iniphysiq
