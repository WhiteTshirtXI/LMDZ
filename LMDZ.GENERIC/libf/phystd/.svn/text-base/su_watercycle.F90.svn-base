      subroutine su_watercycle

      use watercommon_h
      implicit none
#include "comcstfi.h"
#include "callkeys.h"


!==================================================================
!
!     Purpose
!     -------
!     Set up relevant constants and parameters for the water cycle, and water cloud properties
!
!     Authors
!     -------
!     Robin Wordsworth (2010)
!     Jeremy Leconte (2012)
!
!==================================================================

      epsi   = mH2O / mugaz
      RCPD   = cpp 

      !RV = 1000.*R/mH2O
      RV = 1000.*8.314/mH2O ! caution! R is R/mugaz already!

      RCPV   = 1.88e3 ! specific heat capacity of water vapor at 350K

      RVTMP2 = RCPV/RCPD-1. ! not currently used...

      end subroutine su_watercycle
