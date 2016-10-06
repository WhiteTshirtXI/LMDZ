subroutine watersat(T,p,qsat)

  use watercommon_h, only: T_h2O_ice_liq, epsi
  implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the water mass mixing ratio at saturation (kg/kg)
!     for a given pressure (Pa) and temperature (K)
!     A replacement for the old watersat.F in the Martian GCM.
!     Based on FCTTRE.h in the LMDTERRE model.
!
!     Authors
!     -------
!     Robin Wordsworth (2010)
!
!==================================================================

!   input
  real T, p
  
!   output
  real qsat

! checked vs. NIST data 22/06/2010 by RW.
! / by p gives partial pressure
! x by epsi converts to mass mixing ratio

  if (T.lt.T_h2O_ice_liq) then ! solid / vapour
     qsat = 100.0 * 10**(2.07023 - 0.00320991             &
          * T - 2484.896 / T + 3.56654 * alog10(T))
  else                 ! liquid / vapour
     qsat = 100.0 * 10**(23.8319 - 2948.964 / T - 5.028  &
          * alog10(T) - 29810.16 * exp( -0.0699382 * T)  &
          + 25.21935 * exp(-2999.924/T))
  endif
!  qsat=epsi*qsat/p
  if(qsat.gt.p) then
     qsat=1.
  else
     qsat=epsi*qsat/(p-(1.-epsi)*qsat)
  endif
  return
end subroutine watersat
