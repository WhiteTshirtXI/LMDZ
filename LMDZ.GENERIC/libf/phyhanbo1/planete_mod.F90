MODULE planete_mod
  IMPLICIT NONE
  
  REAL :: apoastr ! maximum star-planet distance (AU)
  REAL :: periastr ! minimum star-planet distance (AU)
  REAL :: year_day ! length of year (sols)
  REAL :: peri_day ! date of periastron (sols since N. spring)
  REAL :: obliquit ! Obliquity of the planet (deg)
  REAL :: nres ! tidal resonance ratio
  REAL :: z0 ! surface roughness (m)
  REAL :: lmixmin ! mixing length
  REAL :: emin_turb ! minimal energy
  REAL :: coefvis
  REAL :: coefir
  REAL :: timeperi
  REAL :: e_elips
  REAL :: p_elips
  
  REAL :: preff ! reference surface pressure (Pa)	!read by master
  REAL,ALLOCATABLE :: ap(:) ! hybrid coordinate at layer interface	!read by master
  REAL,ALLOCATABLE :: bp(:) ! hybrid coordinate at layer interface 	!read by master
  
  CONTAINS
  
  subroutine ini_planete_mod(nlayer,preff_dyn,ap_dyn,bp_dyn)
  
  implicit none
  integer,intent(in) :: nlayer ! number of atmospheric layers
  real,intent(in) :: preff_dyn ! reference surface pressure (Pa)
  real,intent(in) :: ap_dyn(nlayer+1) ! hybrid coordinate at interfaces
  real,intent(in) :: bp_dyn(nlayer+1) ! hybrid coordinate at interfaces
  
!$OMP MASTER
  allocate(ap(nlayer+1))
  allocate(bp(nlayer+1))
  
  preff=preff_dyn
  ap(:)=ap_dyn(:)
  bp(:)=bp_dyn(:)
!$OMP END MASTER
!$OMP BARRIER
  
  end subroutine ini_planete_mod
  
END MODULE planete_mod
