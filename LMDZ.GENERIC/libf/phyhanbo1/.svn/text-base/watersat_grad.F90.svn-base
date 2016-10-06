subroutine watersat_grad(T,qsat,dqsat)

  use watercommon_h, only: T_h2O_ice_liq, RLVTT, RCPD,T_coup
  implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the L/cp*d (q_sat)/d T
!     for a given temperature (K)
!
!     Authors
!     -------
!     Robin Wordsworth (2010)
!
!==================================================================

!   input
  real T,qsat
  
!   output
  real dqsat

!  if (T.lt.T_coup) then ! solid / vapour !why use T_coup?????????? JL12
  if (T.lt.T_h2O_ice_liq) then ! solid / vapour
     dqsat = RLVTT/RCPD*qsat*(3.56654/T             &
          +2484.896*LOG(10.)/T**2                   &
          -0.00320991*LOG(10.))
  else                 ! liquid / vapour
     dqsat = RLVTT/RCPD*qsat*LOG(10.)*              &
          (2948.964/T**2-5.028/LOG(10.)/T           &
          +25.21935*2999.924/T**2*EXP(-2999.924/T)  &
          +29810.16*0.0699382*EXP(-0.0699382*T))
  end if

  return
end subroutine watersat_grad

