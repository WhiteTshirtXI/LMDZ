
       module tracer_h

       implicit none

       character*20, DIMENSION(:), ALLOCATABLE :: noms   ! name of the tracer
       real, DIMENSION(:), ALLOCATABLE :: mmol     ! mole mass of tracer (g/mol-1) 
       real, DIMENSION(:), ALLOCATABLE :: radius   ! dust and ice particle radius (m)
       real, DIMENSION(:), ALLOCATABLE :: rho_q    ! tracer densities (kg.m-3)
       real, DIMENSION(:), ALLOCATABLE :: qext     ! Single Scat. Extinction coeff at 0.67 um
       real, DIMENSION(:), ALLOCATABLE :: alpha_lift  ! saltation vertical flux/horiz flux ratio (m-1)
       real, DIMENSION(:), ALLOCATABLE :: alpha_devil ! lifting coeeficient by dust devil
       real, DIMENSION(:), ALLOCATABLE :: qextrhor ! Intermediate for computing opt. depth from q

      real varian      ! Characteristic variance of log-normal distribution
      real r3n_q     ! used to compute r0 from number and mass mixing ratio
      real rho_dust     ! Mars dust density (kg.m-3)
      real rho_ice     ! Water ice density (kg.m-3)
      real rho_co2     ! CO2 ice density (kg.m-3)
      real ref_r0        ! for computing reff=ref_r0*r0 (in log.n. distribution)
!$OMP THREADPRIVATE(noms,mmol,radius,rho_q,qext,alpha_lift,alpha_devil,qextrhor, &
	!$OMP varian,r3n_q,rho_dust,rho_ice,rho_co2,ref_r0)

! tracer indexes: these are initialized in initracer and should be 0 if the
!                 corresponding tracer does not exist
      ! dust
      integer, DIMENSION(:), ALLOCATABLE :: igcm_dustbin ! for dustbin 'dust' tracers
      ! dust, special doubleq case
      integer :: igcm_dust_mass   ! dust mass mixing ratio (for transported dust)
      integer :: igcm_dust_number ! dust number mixing ratio (transported dust)
      ! water
      integer :: igcm_h2o_vap ! water vapour
      integer :: igcm_h2o_ice ! water ice
      ! chemistry:
      integer :: igcm_co2
      integer :: igcm_co
      integer :: igcm_o
      integer :: igcm_o1d
      integer :: igcm_o2
      integer :: igcm_o3
      integer :: igcm_h
      integer :: igcm_h2
      integer :: igcm_oh
      integer :: igcm_ho2
      integer :: igcm_h2o2
      integer :: igcm_n2
      integer :: igcm_ar
      ! other tracers
      integer :: igcm_ar_n2 ! for simulations using co2 +neutral gaz
      integer :: igcm_co2_ice ! CO2 ice 
!$OMP THREADPRIVATE(igcm_dustbin,igcm_dust_mass,igcm_dust_number,igcm_h2o_vap,igcm_h2o_ice, &
	!$OMP igcm_co2,igcm_co,igcm_o,igcm_o1d,igcm_o2,igcm_o3,igcm_h,igcm_h2,igcm_oh,	    &
	!$OMP igcm_ho2,igcm_h2o2,igcm_n2,igcm_ar,igcm_ar_n2,igcm_co2_ice)

       end module tracer_h

