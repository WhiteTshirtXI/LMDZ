subroutine iniwritesoil(nid,ngrid)

! initialization routine for 'writediagoil'. Here we create/define
! dimensions (longitude, latitude, depth and time) and other fixed
! (time-independent) parameters.

use comsoil_h, only: mlayer, inertiedat, nsoilmx

implicit none

#include"dimensions.h"
#include"paramet.h"
#include"comcstfi.h"
#include"comgeom.h"
#include"netcdf.inc"

! Arguments:
integer,intent(in) :: ngrid
integer,intent(in) :: nid ! NetCDF output file ID

! Local variables:

! NetCDF stuff:
integer :: ierr ! NetCDF routines return code
integer :: idim_rlatu ! ID of the 'latitude' dimension
integer :: idim_rlonv ! ID of the 'longitude' dimension
integer :: idim_depth ! ID of the 'depth' dimension
integer :: idim_time  ! ID of the 'time' dimension
integer :: varid ! to store NetCDF ID of a variable
integer,dimension(3) :: dimids ! to store IDs of dimensions of a variable
character(len=60) :: text ! to store some text
real,dimension(iip1,jjp1,nsoilmx) :: data3 ! to store 3D data
integer :: i,j,l,ig0

! 1. Define the dimensions
! Switch to NetCDF define mode
ierr=NF_REDEF(nid)

! Define the dimensions
ierr=NF_DEF_DIM(nid,"longitude",iip1,idim_rlonv)
! iip1 known from paramet.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define longitude dimension"
endif
ierr=NF_DEF_DIM(nid,"latitude",jjp1,idim_rlatu)
! jjp1 known from paramet.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define latitude dimension"
endif
ierr=NF_DEF_DIM(nid,"depth",nsoilmx,idim_depth)
! nsoilmx known from dimphys.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define depth dimension"
endif
ierr=NF_DEF_DIM(nid,"time",NF_UNLIMITED,idim_time)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define time dimension"
endif

! Switch out of NetCDF define mode
ierr=NF_ENDDEF(nid)

! 2. Define (as variables) and write dimensions, as well as their attributes
! 2.1. Longitude
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"longitude",NF_DOUBLE,1,idim_rlonv,varid)
#else
ierr=NF_DEF_VAR(nid,"longitude",NF_FLOAT,1,idim_rlonv,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define longitude variable"
endif

! Longitude attributes
text="East longitude"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="degrees_east"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write longitude to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,rlonv*(180./pi))
#else
ierr=NF_PUT_VAR_REAL(nid,varid,rlonv*(180./pi))
#endif
! Note: rlonv is known from comgeom.h and pi from comcstfi.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write longitude variable"
endif

! 2.2. Latitude
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"latitude",NF_DOUBLE,1,idim_rlatu,varid)
#else
ierr=NF_DEF_VAR(nid,"latitude",NF_FLOAT,1,idim_rlatu,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define latitude variable"
endif

! Latitude attributes
text="North latitude"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="degrees_north"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write latitude to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,rlatu*(180./pi))
#else
ierr=NF_PUT_VAR_REAL(nid,varid,rlatu*(180./pi))
#endif
! Note: rlatu is known from comgeom.h and pi from comcstfi.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write longitude variable"
endif

! 2.3. Depth
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"depth",NF_DOUBLE,1,idim_depth,varid)
#else
ierr=NF_DEF_VAR(nid,"depth",NF_FLOAT,1,idim_depth,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define depth variable"
endif

! Depth attributes
text="Soil mid-layer depth"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="m"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)
text="down"
ierr=NF_PUT_ATT_TEXT(nid,varid,"positive",len_trim(text),text)

! Write depth to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,mlayer)
#else
ierr=NF_PUT_VAR_REAL(nid,varid,mlayer)
#endif
! Note mlayer(0:nsoilmx-1) known from comsoil_h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write depth variable"
endif

! 2.4. Time
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"time",NF_DOUBLE,1,idim_time,varid)
#else
ierr=NF_DEF_VAR(nid,"time",NF_FLOAT,1,idim_time,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define depth variable"
endif

! time attributes
text="Time"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="days since 0000-01-01 00:00:00"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Note no need to write time variable here; it is done in writediagsoil.

! 3. Other variables to be included

! 3.1 mesh area surrounding each horizontal point
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
dimids(1)=idim_rlonv ! ID of the 'longitude' dimension
dimids(2)=idim_rlatu ! ID of the 'latitude' dimension
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"area",NF_DOUBLE,2,dimids,varid)
#else
ierr=NF_DEF_VAR(nid,"area",NF_FLOAT,2,dimids,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define area variable"
endif

! Area attributes
text="Mesh area"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="m2"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write area to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,aire)
#else
ierr=NF_PUT_VAR_REAL(nid,varid,aire)
#endif
! Note: aire is known from comgeom.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write area variable"
endif

! 3.2 Thermal inertia
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
dimids(1)=idim_rlonv ! ID of the 'longitude' dimension
dimids(2)=idim_rlatu ! ID of the 'latitude' dimension
dimids(3)=idim_depth ! ID of the 'depth' dimension
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,"th_inertia",NF_DOUBLE,3,dimids,varid)
#else
ierr=NF_DEF_VAR(nid,"th_inertia",NF_FLOAT,3,dimids,varid)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define th_inertia variable"
endif

! Attributes
text="Thermal inertia"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="J.s-1/2.m-2.K-1"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Recast data along 'dynamics' grid
! Note: inertiedat is known from comsoil_h

do l=1,nsoilmx
  ! handle the poles
  do i=1,iip1
    data3(i,1,l)=inertiedat(1,l)
    data3(i,jjp1,l)=inertiedat(ngrid,l)
  enddo
  ! rest of the grid
  do j=2,jjm
    ig0=1+(j-2)*iim
    do i=1,iim
      data3(i,j,l)=inertiedat(ig0+i,l)
    enddo
    data3(iip1,j,l)=data3(1,j,l) ! extra (modulo) longitude
  enddo
enddo ! of do l=1,nsoilmx

! Write data2 to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
#ifdef NC_DOUBLE
ierr=NF_PUT_VAR_DOUBLE(nid,varid,data3)
#else
ierr=NF_PUT_VAR_REAL(nid,varid,data3)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write th_inertia variable"
endif

end subroutine iniwritesoil
