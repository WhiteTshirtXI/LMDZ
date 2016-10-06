      module radinc_h

      implicit none

#include "dimensions.h"
#include "bands.h"
#include "scatterers.h"

!======================================================================
!
!     RADINC.H 
!
!     Includes for the radiation code; RADIATION LAYERS, LEVELS,
!     number of spectral intervals. . .
! 
!======================================================================

!     RADIATION parameters

!     In radiation code, layer 1 corresponds to the stratosphere.  Level
!     1 is the top of the stratosphere.  The dummy layer is at the same
!     temperature as the (vertically isothermal) stratosphere, and
!     any time it is explicitly needed, the appropriate quantities will
!     be dealt with (aka "top". . .)

!     L_NLEVRAD corresponds to the surface - i.e., the GCM Level that
!     is at the surface.  PLEV(L_NLEVRAD) = P(J,I)+PTROP, 
!     PLEV(2) = PTROP, PLEV(1) = ptop

!     L_NLAYRAD is the number of radiation code layers
!     L_NLEVRAD is the number of radiation code levels.  Level N is the
!               top of layer N. 
!
!     L_NSPECTI is the number of IR spectral intervals
!     L_NSPECTV is the number of Visual(or Solar) spectral intervals
!     L_NGAUSS  is the number of Gauss points for K-coefficients
!               GAUSS POINT 17 (aka the last one) is the special case
!
!     L_NPREF   is the number of reference pressures that the 
!               k-coefficients are calculated on
!     L_PINT    is the number of Lagrange interpolated reference
!               pressures for the gas k-coefficients - now for a
!		smaller p-grid than before
!     L_NTREF   is the number of reference temperatures for the
!               k-coefficients
!     L_TAUMAX  is the largest optical depth - larger ones are set
!               to this value
!
!     L_REFVAR  The number of different mixing ratio values for
!               the k-coefficients. Variable component of the mixture
!		can in princple be anything: currently it's H2O.
!
!     NAERKIND  The number of radiatively active aerosol types
!
!     NSIZEMAX  The maximum number of aerosol particle sizes
!
!----------------------------------------------------------------------

      integer, parameter :: L_NLAYRAD  = llm
      integer, parameter :: L_LEVELS   = 2*(llm-1)+3
      integer, parameter :: L_NLEVRAD  = llm+1

      ! These are set in sugas_corrk
      ! [uses allocatable arrays] -- AS 12/2011
      integer :: L_NPREF, L_NTREF, L_REFVAR, L_PINT   !L_NPREF, L_NTREF, L_REFVAR, L_PINT read by master in sugas_corrk

      integer, parameter :: L_NGAUSS  = 17

      integer, parameter :: L_NSPECTI = NBinfrared
      integer, parameter :: L_NSPECTV = NBvisible

!      integer, parameter :: NAERKIND  = 2 ! set in scatterers.h
      real,    parameter :: L_TAUMAX  = 35

      ! For Planck function integration: 
      ! equivalent temperatures are 1/NTfac of these values
      integer, parameter :: NTstar = 500
      integer, parameter :: NTstop = 15000 ! new default for all non hot Jupiter runs
      real*8, parameter :: NTfac = 1.0D+1  
      !integer, parameter :: NTstar = 1000
      !integer, parameter :: NTstop = 25000
      !real*8,parameter :: NTfac = 5.0D+1    
      !integer, parameter :: NTstar = 2000
      !integer, parameter :: NTstop = 50000
      !real*8,parameter :: NTfac = 1.0D+2    

      ! Maximum number of grain size classes for aerosol convolution:
      ! This must correspond to size of largest dataset used for aerosol 
      ! optical properties in datagcm folder.
      integer, parameter :: nsizemax = 60

      character (len=100) :: corrkdir
      save corrkdir
!$OMP THREADPRIVATE(corrkdir)

      character (len=100) :: banddir
      save banddir
!$OMP THREADPRIVATE(banddir)

      end module radinc_h
