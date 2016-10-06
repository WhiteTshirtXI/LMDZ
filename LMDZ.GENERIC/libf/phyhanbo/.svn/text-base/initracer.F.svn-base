      SUBROUTINE initracer(ngrid,nq,nametrac)

      use surfdat_h
      USE comgeomfi_h
      USE tracer_h
      IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Initialization related to tracer 
c   (transported dust, water, chemical species, ice...)
c
c   Name of the tracer
c
c   Test of dimension :
c   Initialize COMMON tracer in tracer.h, using tracer names provided
c   by the argument nametrac
c
c   author: F.Forget
c   ------
c            Ehouarn Millour (oct. 2008) identify tracers by their names
c=======================================================================


#include "dimensions.h"
!#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"

      integer :: ngrid,nq

!      real qsurf(ngrid,nq)       ! tracer on surface (e.g.  kg.m-2)
!      real co2ice(ngrid)           ! co2 ice mass on surface (e.g.  kg.m-2)
      character(len=20) :: txt ! to store some text
      integer iq,ig,count
      real r0_lift , reff_lift
!      logical :: oldnames ! =.true. if old tracer naming convention (q01,...)

      character*20 nametrac(nq)   ! name of the tracer from dynamics


c-----------------------------------------------------------------------
c  radius(nq)      ! aerosol particle radius (m)
c  rho_q(nq)       ! tracer densities (kg.m-3)
c  qext(nq)        ! Single Scat. Extinction coeff at 0.67 um
c  alpha_lift(nq)  ! saltation vertical flux/horiz flux ratio (m-1)
c  alpha_devil(nq) ! lifting coeeficient by dust devil
c  rho_dust          ! Mars dust density
c  rho_ice           ! Water ice density
c  doubleq           ! if method with mass (iq=1) and number(iq=2) mixing ratio
c  varian            ! Characteristic variance of log-normal distribution
c-----------------------------------------------------------------------

       !! we allocate once for all arrays in common in tracer_h.F90
       !! (supposedly those are not used before call to initracer)
       IF (.NOT.ALLOCATED(noms)) ALLOCATE(noms(nq))
       ALLOCATE(mmol(nq))
       ALLOCATE(radius(nq))
       ALLOCATE(rho_q(nq))
       ALLOCATE(qext(nq))
       ALLOCATE(alpha_lift(nq))
       ALLOCATE(alpha_devil(nq))
       ALLOCATE(qextrhor(nq))
       ALLOCATE(igcm_dustbin(nq))
       !! initialization
       alpha_lift(:)=0.
       alpha_devil(:)=0.

! Initialization: get tracer names from the dynamics and check if we are
!                 using 'old' tracer convention ('q01',q02',...)
!                 or new convention (full tracer names)
      ! check if tracers have 'old' names

! copy tracer names from dynamics
        do iq=1,nq
          noms(iq)=nametrac(iq)
        enddo


! Identify tracers by their names: (and set corresponding values of mmol)
      ! 0. initialize tracer indexes to zero:
      ! NB: igcm_* indexes are commons in 'tracer.h'
      do iq=1,nq
        igcm_dustbin(iq)=0
      enddo
      igcm_dust_mass=0
      igcm_dust_number=0
      igcm_h2o_vap=0
      igcm_h2o_ice=0
      igcm_co2=0
      igcm_co=0
      igcm_o=0
      igcm_o1d=0
      igcm_o2=0
      igcm_o3=0
      igcm_h=0
      igcm_h2=0
      igcm_oh=0
      igcm_ho2=0
      igcm_h2o2=0
      igcm_n2=0
      igcm_ar=0
      igcm_ar_n2=0
      igcm_co2_ice=0

      write(*,*) 'initracer: noms() ', noms


      !print*,'Setting dustbin = 0 in initracer.F'
      !dustbin=0

      ! 1. find dust tracers
      count=0
!      if (dustbin.gt.0) then
!        do iq=1,nq
!          txt=" "
!          write(txt,'(a4,i2.2)')'dust',count+1   
!          if (noms(iq).eq.txt) then
!            count=count+1
!            igcm_dustbin(count)=iq
!            mmol(iq)=100.
!          endif
!        enddo !do iq=1,nq
!      endif ! of if (dustbin.gt.0)


!      if (doubleq) then
!        do iq=1,nq
!          if (noms(iq).eq."dust_mass") then
!            igcm_dust_mass=iq
!            count=count+1
!          endif
!          if (noms(iq).eq."dust_number") then
!            igcm_dust_number=iq
!            count=count+1
!          endif
!        enddo
!      endif ! of if (doubleq)
      ! 2. find chemistry and water tracers
      do iq=1,nq
        if (noms(iq).eq."co2") then
          igcm_co2=iq
          mmol(igcm_co2)=44.
          count=count+1
!          write(*,*) 'co2: count=',count
        endif
        if (noms(iq).eq."co2_ice") then
          igcm_co2_ice=iq
          mmol(igcm_co2_ice)=44.
          count=count+1
!          write(*,*) 'co2_ice: count=',count
        endif
        if (noms(iq).eq."h2o_vap") then
          igcm_h2o_vap=iq
          mmol(igcm_h2o_vap)=18.
          count=count+1
!          write(*,*) 'h2o_vap: count=',count
        endif
        if (noms(iq).eq."h2o_ice") then
          igcm_h2o_ice=iq
          mmol(igcm_h2o_ice)=18.
          count=count+1
!          write(*,*) 'h2o_ice: count=',count
        endif
      enddo ! of do iq=1,nq
      
      ! check that we identified all tracers:
      if (count.ne.nq) then
        write(*,*) "initracer: found only ",count," tracers"
        write(*,*) "               expected ",nq
        do iq=1,count
          write(*,*)'      ',iq,' ',trim(noms(iq))
        enddo
!        stop
      else
        write(*,*) "initracer: found all expected tracers, namely:"
        do iq=1,nq
          write(*,*)'      ',iq,' ',trim(noms(iq))
        enddo
      endif


c------------------------------------------------------------
c     Initialisation tracers ....
c------------------------------------------------------------
      call zerophys(nq,rho_q)

      rho_dust=2500.  ! Mars dust density (kg.m-3)
      rho_ice=920.    ! Water ice density (kg.m-3)
      rho_co2=1620.   ! CO2 ice density (kg.m-3)



c$$$      if (doubleq) then
c$$$c       "doubleq" technique 
c$$$c       -------------------
c$$$c      (transport of mass and number mixing ratio)
c$$$c       iq=1: Q mass mixing ratio, iq=2: N number mixing ratio
c$$$
c$$$        if( (nq.lt.2).or.(water.and.(nq.lt.3)) ) then
c$$$          write(*,*)'initracer: nq is too low : nq=', nq
c$$$          write(*,*)'water= ',water,' doubleq= ',doubleq   
c$$$        end if
c$$$
c$$$        varian=0.637    ! Characteristic variance    
c$$$        qext(igcm_dust_mass)=3.04   ! reference extinction at 0.67 um for ref dust
c$$$        qext(igcm_dust_number)=3.04 ! reference extinction at 0.67 um for ref dust
c$$$        rho_q(igcm_dust_mass)=rho_dust
c$$$        rho_q(igcm_dust_number)=rho_dust
c$$$
c$$$c       Intermediate calcul for computing geometric mean radius r0
c$$$c       as a function of mass and number mixing ratio Q and N
c$$$c       (r0 = (r3n_q * Q/ N)
c$$$        r3n_q = exp(-4.5*varian**2)*(3./4.)/(pi*rho_dust)
c$$$
c$$$c       Intermediate calcul for computing effective radius reff
c$$$c       from geometric mean radius r0
c$$$c       (reff = ref_r0 * r0)
c$$$        ref_r0 = exp(2.5*varian**2)
c$$$        
c$$$c       lifted dust :
c$$$c       '''''''''''
c$$$        reff_lift = 3.e-6      !  Effective radius of lifted dust (m)
c$$$        alpha_devil(igcm_dust_mass)=9.e-9   !  dust devil lift mass coeff
c$$$        alpha_lift(igcm_dust_mass)=3.0e-15  !  Lifted mass coeff
c$$$
c$$$        r0_lift = reff_lift/ref_r0
c$$$        alpha_devil(igcm_dust_number)=r3n_q*
c$$$     &                        alpha_devil(igcm_dust_mass)/r0_lift**3
c$$$        alpha_lift(igcm_dust_number)=r3n_q*
c$$$     &                        alpha_lift(igcm_dust_mass)/r0_lift**3
c$$$
c$$$c       Not used:
c$$$        radius(igcm_dust_mass) = 0.
c$$$        radius(igcm_dust_number) = 0.
c$$$
c$$$      else


c$$$       if (dustbin.gt.1) then
c$$$        print*,'ATTENTION:',
c$$$     $   ' properties of dust need input in initracer !!!'
c$$$        stop
c$$$
c$$$       else if (dustbin.eq.1) then
c$$$
c$$$c       This will be used for 1 dust particle size:
c$$$c       ------------------------------------------
c$$$        radius(igcm_dustbin(1))=3.e-6
c$$$        Qext(igcm_dustbin(1))=3.04
c$$$        alpha_lift(igcm_dustbin(1))=0.0e-6
c$$$        alpha_devil(igcm_dustbin(1))=7.65e-9
c$$$        qextrhor(igcm_dustbin(1))=(3./4.)*Qext(igcm_dustbin(1))
c$$$     &                         /(rho_dust*radius(igcm_dustbin(1)))
c$$$        rho_q(igcm_dustbin(1))=rho_dust
c$$$
c$$$       endif
c$$$!      end if    ! (doubleq)

c     Initialization for water vapor
c     ------------------------------
      if(water) then
         radius(igcm_h2o_vap)=0.
         Qext(igcm_h2o_vap)=0.
         alpha_lift(igcm_h2o_vap) =0.
         alpha_devil(igcm_h2o_vap)=0.
	 qextrhor(igcm_h2o_vap)= 0.

c       "Dryness coefficient" controlling the evaporation and
c        sublimation from the ground water ice (close to 1)
c        HERE, the goal is to correct for the fact
c        that the simulated permanent water ice polar caps
c        is larger than the actual cap and the atmospheric
c        opacity not always realistic.


!         if(ngrid.eq.1)


!     to be modified for BC+ version?

         !! this is defined in surfdat_h.F90
         IF (.not.ALLOCATED(dryness)) ALLOCATE(dryness(ngrid))
         IF (.not.ALLOCATED(watercaptag)) ALLOCATE(watercaptag(ngrid))

         do ig=1,ngrid
           if (ngrid.ne.1) watercaptag(ig)=.false.
           dryness(ig) = 1.
         enddo




!         IF (caps) THEN
c Perennial H20 north cap defined by watercaptag=true (allows surface to be
c hollowed by sublimation in vdifc).
!         do ig=1,ngrid
!           if (lati(ig)*180./pi.gt.84) then
!             if (ngrid.ne.1) watercaptag(ig)=.true.
!             dryness(ig) = 1.
c Use the following cap definition for high spatial resolution (latitudinal bin <= 5 deg)
c	      if (lati(ig)*180./pi.lt.85.and.long(ig).ge.0) then
c	        if (ngrid.ne.1) watercaptag(ig)=.true.
c               dryness(ig) = 1.
c	      endif
c             if (lati(ig)*180./pi.ge.85) then
c               if (ngrid.ne.1) watercaptag(ig)=.true.
c               dryness(ig) = 1.
c	      endif
!           endif  ! (lati>80 deg)
!         end do ! (ngrid)
!        ENDIF ! (caps)

!         if(iceparty.and.(nq.ge.2)) then

           radius(igcm_h2o_ice)=3.e-6
           rho_q(igcm_h2o_ice)=rho_ice
           Qext(igcm_h2o_ice)=0.
!           alpha_lift(igcm_h2o_ice) =0.
!           alpha_devil(igcm_h2o_ice)=0.
           qextrhor(igcm_h2o_ice)= (3./4.)*Qext(igcm_h2o_ice) 
     $       / (rho_ice*radius(igcm_h2o_ice))



!         elseif(iceparty.and.(nq.lt.2)) then
!            write(*,*) 'nq is too low : nq=', nq
!            write(*,*) 'water= ',water,' iceparty= ',iceparty
!         endif

      end if  ! (water)

c     Output for records:
c     ~~~~~~~~~~~~~~~~~~
      write(*,*)
      Write(*,*) '******** initracer : dust transport parameters :'
      write(*,*) 'alpha_lift = ', alpha_lift
      write(*,*) 'alpha_devil = ', alpha_devil
      write(*,*) 'radius  = ', radius
!      if(doubleq) then
!        write(*,*) 'reff_lift (um) =  ', reff_lift
!        write(*,*) 'size distribution variance  = ', varian
!        write(*,*) 'r3n_q , ref_r0 : ', r3n_q , ref_r0
!      end if
      write(*,*) 'Qext  = ', qext 
      write(*,*)

      end
