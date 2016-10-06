      subroutine orbite(pls,pdist_star,pdecli)
      
      use planete_mod, only: p_elips, e_elips, timeperi, obliquit
      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Distance from star and declination as a function of the stellar
!     longitude Ls
!     
!     Inputs
!     ------
!     pls          Ls
!
!     Outputs
!     -------
!     pdist_star    Distance Star-Planet in UA
!     pdecli        declinaison ( in radians )
!
!=======================================================================

c   Declarations:
c   -------------

!#include "planete.h"
#include "comcstfi.h"

c arguments:
c ----------

      REAL pday,pdist_star,pdecli,pls,i

c-----------------------------------------------------------------------

c Star-Planet Distance

      pdist_star = p_elips/(1.+e_elips*cos(pls+timeperi))

c Stellar declination

c ********************* version before 01/01/2000 *******

      pdecli = asin (sin(pls)*sin(obliquit*pi/180.))

c ********************* version after 01/01/2000 *******
c     i=obliquit*pi/180.
c     pdecli=asin(sin(pls)*sin(i)/sqrt(sin(pls)**2+
c    & cos(pls)**2*cos(i)**2))
c ******************************************************

      RETURN
      END
