module comsoil_h

implicit none
! nsoilmx : number of subterranean layers
!integer, parameter :: nsoilmx = 18 ! for z1=0.0002 m, depth = 18 m => mars case 
!integer, parameter :: nsoilmx = 13 ! for z1=0.03 m, depth = 104.8 m => earth case
  integer, parameter :: nsoilmx = 13

  real,save,allocatable,dimension(:) :: layer      ! soil layer depths
  real,save,allocatable,dimension(:) :: mlayer     ! soil mid-layer depths
  real,save,allocatable,dimension(:,:) :: inertiedat ! soil thermal inertia
  real,save :: volcapa    ! soil volumetric heat capacity
       ! NB: volcapa is read fromn control(35) from physicq start file
       !     in physdem (or set via tabfi, or initialized in
       !                 soil_settings.F)
!$OMP THREADPRIVATE(layer,mlayer,inertiedat,volcapa)

contains

  subroutine ini_comsoil_h(ngrid)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  
    if (.not.allocated(layer)) allocate(layer(nsoilmx)) !soil layer depths
    if (.not.allocated(mlayer)) allocate(mlayer(0:nsoilmx-1)) ! soil mid-layer depths
    if (.not.allocated(inertiedat)) allocate(inertiedat(ngrid,nsoilmx)) ! soil thermal inertia
  
  end subroutine ini_comsoil_h

end module comsoil_h

