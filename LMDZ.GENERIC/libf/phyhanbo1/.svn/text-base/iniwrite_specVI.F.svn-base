      SUBROUTINE iniwrite_specVI(nid,idayref)

      use radinc_h, only: L_NSPECTV
      use radcommon_h, only: WNOV,DWNV
      use comsoil_h

      implicit none

c=======================================================================
c
c   Auteur:  L. Fairhead  ,  P. Le Van, Y. Wanherdrick, F. Forget
c   -------
c
c   Objet:
c   ------
c
c   'Initialize' the diagfi_spec.nc file: write down dimensions as well
c   as time-independent fields (e.g: geopotential, mesh area, ...)
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
!include "comconst.h"
#include "comcstfi.h"
#include "comvert.h"
#include "comgeom.h"
#include "temps.h"
#include "ener.h"
#include "logic.h"
#include "netcdf.inc"
#include "serre.h"
!#include"dimphys.h"

c   Arguments:
c   ----------

      integer nid        ! NetCDF file ID
      INTEGER*4 idayref  ! date (initial date for this run)

c   Local:
c   ------
      INTEGER length,l
      parameter (length = 100)
      REAL tab_cntrl(length) ! run parameters are stored in this array
      INTEGER ierr

      integer :: nvarid,idim_index,idim_rlonu,idim_rlonv
      integer :: idim_rlatu,idim_rlatv,idim_llmp1,idim_llm
 !     integer :: idim_nsoilmx ! "subsurface_layers" dimension ID #
      integer :: idim_bandsVI ! "VI Wavenumber" dimension ID #
      integer, dimension(2) :: id  

c-----------------------------------------------------------------------

      DO l=1,length
         tab_cntrl(l)=0.
      ENDDO
      tab_cntrl(1)  = FLOAT(iim)
      tab_cntrl(2)  = FLOAT(jjm)
      tab_cntrl(3)  = FLOAT(llm)
      tab_cntrl(4)  = FLOAT(idayref)
      tab_cntrl(5)  = rad
      tab_cntrl(6)  = omeg
      tab_cntrl(7)  = g
      tab_cntrl(8)  = mugaz
      tab_cntrl(9)  = rcp
      tab_cntrl(10) = daysec
      tab_cntrl(11) = dtphys
      tab_cntrl(12) = etot0
      tab_cntrl(13) = ptot0
      tab_cntrl(14) = ztot0
      tab_cntrl(15) = stot0
      tab_cntrl(16) = ang0
c
c    ..........    P.Le Van  ( ajout le 8/04/96 )    .........
c         .....        parametres  pour le zoom          ......   
      tab_cntrl(17)  = clon
      tab_cntrl(18)  = clat
      tab_cntrl(19)  = grossismx
      tab_cntrl(20)  = grossismy
c
c     .....   ajout  le 6/05/97 et le 15/10/97  .......
c
      IF ( fxyhypb )   THEN
        tab_cntrl(21) = 1.
        tab_cntrl(22) = dzoomx
        tab_cntrl(23) = dzoomy
      ELSE
        tab_cntrl(21) = 0.
        tab_cntrl(22) = dzoomx
        tab_cntrl(23) = dzoomy
        tab_cntrl(24) = 0.
        IF( ysinus )  tab_cntrl(24) = 1.
      ENDIF

c    .........................................................

! Define dimensions
    
      ierr = NF_REDEF (nid)

      ierr = NF_DEF_DIM (nid, "index", length, idim_index)
      ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_rlatu)
      ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_rlonv)
      ierr = NF_DEF_DIM (nid, "VI Wavenumber",L_NSPECTV,idim_bandsVI)

      ierr = NF_ENDDEF(nid)

c  Contol parameters for this run
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "controle", NF_DOUBLE, 1, 
     .       idim_index,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "controle", NF_FLOAT, 1, 
     .       idim_index,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 18,
     .                       "Control parameters")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif

c --------------------------
c  longitudes and latitudes

      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid, "latitude", NF_DOUBLE, 1, idim_rlatu,nvarid)
#else
      ierr =NF_DEF_VAR (nid, "latitude", NF_FLOAT, 1, idim_rlatu,nvarid)
#endif
      ierr =NF_PUT_ATT_TEXT(nid,nvarid,'units',13,"degrees_north")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 14,
     .      "North latitude")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu/pi*180)
#endif
c
c --------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid,"longitude", NF_DOUBLE, 1, idim_rlonv,nvarid)
#else
      ierr = NF_DEF_VAR(nid,"longitude", NF_FLOAT, 1, idim_rlonv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 14,
     .      "East longitude")
      ierr = NF_PUT_ATT_TEXT(nid,nvarid,'units',12,"degrees_east")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv/pi*180)
#endif
c

!-------------------------------
! Number of bands in the IR
!-------------------------------

      ierr=NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      ! define variable
#ifdef NC_DOUBLE
      ierr=NF_DEF_VAR(nid,"VI Wavenumber",NF_DOUBLE,1,
     .				idim_bandsVI,nvarid)
#else
      ierr=NF_DEF_VAR(nid,"VI Wavenumber",NF_FLOAT,1,
     .				idim_bandsVI,nvarid)
#endif
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 33,
     .                        "Band mid frequency in the visible")
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"units",5,"cm^-1")
      ierr=NF_ENDDEF(nid) ! Leave NetCDF define mode
      ! write variable
#ifdef NC_DOUBLE
      ierr=NF_PUT_VAR_DOUBLE (nid,nvarid,dble(WNOV))
#else
      ierr=NF_PUT_VAR_REAL (nid,nvarid,real(WNOV))
#endif

c
!-------------------------------
! Width of bands in the Visible
!-------------------------------

      ierr=NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      ! define variable
#ifdef NC_DOUBLE
      ierr=NF_DEF_VAR(nid,"VI Bandwidth",NF_DOUBLE,1,
     .				idim_bandsVI,nvarid)
#else
      ierr=NF_DEF_VAR(nid,"VI Bandwidth",NF_FLOAT,1,
     .				idim_bandsVI,nvarid)
#endif
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 24,
     .                        "Bandwidth in the visible")
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"units",5,"cm^-1")
      ierr=NF_ENDDEF(nid) ! Leave NetCDF define mode
      ! write variable
#ifdef NC_DOUBLE
      ierr=NF_PUT_VAR_DOUBLE (nid,nvarid,dble(DWNV))
#else
      ierr=NF_PUT_VAR_REAL (nid,nvarid,real(DWNV))
#endif
c
c --------------------------
c   Mesh area
c --------------------------

      id(1)=idim_rlonv
      id(2)=idim_rlatu


      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "aire", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "aire", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Mesh area")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aire)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aire)
#endif



      END
