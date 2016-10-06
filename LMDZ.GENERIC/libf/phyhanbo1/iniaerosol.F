      SUBROUTINE iniaerosol()


      use radinc_h, only: naerkind
      use aerosol_mod

      IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Initialization related to aerosols 
c   (CO2 aerosols, dust, water, chemical species, ice...)   
c
c   author: Laura Kerber, S. Guerlet
c   ------
c        
c=======================================================================


#include "callkeys.h"
      integer ia

      ia=0
      if (aeroco2) then
         ia=ia+1
         iaero_co2=ia
      endif
      write(*,*) '--- CO2 aerosol = ', iaero_co2
 
      if (aeroh2o) then
         ia=ia+1
         iaero_h2o=ia
      endif
      write(*,*) '--- H2O aerosol = ', iaero_h2o

      if (dusttau.gt.0) then
         ia=ia+1
         iaero_dust=ia
      endif
      write(*,*) '--- Dust aerosol = ', iaero_dust

      if (aeroh2so4) then
         ia=ia+1
         iaero_h2so4=ia
      endif
      write(*,*) '--- H2SO4 aerosol = ', iaero_h2so4
      
      if (aeroback2lay) then
         ia=ia+1
         iaero_back2lay=ia
      endif
      write(*,*) '--- Two-layer aerosol = ', iaero_back2lay

      write(*,*) '=== Number of aerosols= ', ia
      
! For the zero aerosol case, we currently make a dummy co2 aerosol which is zero everywhere.
! (See aeropacity.F90 for how this works). A better solution would be to turn off the 
! aerosol machinery in the no aerosol case, but this would be complicated. LK

      if (ia.eq.0) then  !For the zero aerosol case. 
         ia = 1
         noaero = .true.
         iaero_co2=ia
      endif

      if (ia.ne.naerkind) then
          print*, 'Aerosols counted not equal to naerkind'
          print*, 'Compile with tag -s',ia,'to run'
          print*, 'or change options in callphys.def'
          print*, 'Abort in iniaerosol.F'
          call abort
      endif

      end
