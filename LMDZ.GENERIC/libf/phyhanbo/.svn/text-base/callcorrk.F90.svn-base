      subroutine callcorrk(ngrid,nlayer,pq,nq,qsurf,           &
          albedo,emis,mu0,pplev,pplay,pt,                      & 
          tsurf,fract,dist_star,aerosol,muvar,                 &
          dtlw,dtsw,fluxsurf_lw,                               &
          fluxsurf_sw,fluxtop_lw,fluxabs_sw,fluxtop_dn,        &
          OLR_nu,OSR_nu,                                       &
          tau_col,cloudfrac,totcloudfrac,                      &
          clearsky,firstcall,lastcall)

      use radinc_h
      use radcommon_h
      use watercommon_h
      use datafile_mod, only: datadir
!      use ioipsl_getincom 
      use ioipsl_getincom_p
      use gases_h
      use radii_mod, only : su_aer_radii,co2_reffrad,h2o_reffrad,dust_reffrad,h2so4_reffrad,back2lay_reffrad
      use aerosol_mod, only : iaero_co2,iaero_h2o,iaero_dust,iaero_h2so4, iaero_back2lay
      USE tracer_h

      implicit none

!==================================================================
!
!     Purpose
!     -------
!     Solve the radiative transfer using the correlated-k method for
!     the gaseous absorption and the Toon et al. (1989) method for
!     scatttering due to aerosols.
!
!     Authors
!     ------- 
!     Emmanuel 01/2001, Forget 09/2001
!     Robin Wordsworth (2009)
!
!==================================================================

!#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"

!-----------------------------------------------------------------------
!     Declaration of the arguments (INPUT - OUTPUT) on the LMD GCM grid
!     Layer #1 is the layer near the ground. 
!     Layer #nlayer is the layer at the top.

      INTEGER,INTENT(IN) :: ngrid ! number of atmospheric columns
      INTEGER,INTENT(IN) :: nlayer ! number of atmospheric layers
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq) ! tracers (.../kg_of_air)
      integer,intent(in) :: nq ! number of tracers
      REAL,INTENT(IN) :: qsurf(ngrid,nq) ! tracer on surface (kg.m-2)
      REAL,INTENT(IN) :: albedo(ngrid)   ! SW albedo
      REAL,INTENT(IN) :: emis(ngrid)     ! LW emissivity
      real,intent(in) :: mu0(ngrid) ! cosine of sun incident angle
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1)  ! inter-layer pressure (Pa)
      REAL,INTENT(IN) :: pplay(ngrid,nlayer)    ! mid-layer pressure (Pa)
      REAL,INTENT(IN) :: pt(ngrid,nlayer)  ! air temperature (K)
      REAL,INTENT(IN) :: tsurf(ngrid)        ! surface temperature (K)
      REAL,INTENT(IN) :: fract(ngrid)        ! fraction of day
      REAL,INTENT(IN) :: dist_star           ! distance star-planet (AU)
      REAL,INTENT(OUT) :: aerosol(ngrid,nlayer,naerkind) ! aerosol tau (kg/kg)
      real,intent(in) :: muvar(ngrid,nlayer+1)
      REAL,INTENT(OUT) :: dtlw(ngrid,nlayer) ! heating rate (K/s) due to LW
      REAL,INTENT(OUT) :: dtsw(ngrid,nlayer) ! heating rate (K/s) due to SW
      REAL,INTENT(OUT) :: fluxsurf_lw(ngrid)   ! incident LW flux to surf (W/m2)
      REAL,INTENT(OUT) :: fluxsurf_sw(ngrid)   ! incident SW flux to surf (W/m2)
      REAL,INTENT(OUT) :: fluxtop_lw(ngrid)    ! outgoing LW flux to space (W/m2)
      REAL,INTENT(OUT) :: fluxabs_sw(ngrid)    ! SW flux absorbed by planet (W/m2)
      REAL,INTENT(OUT) :: fluxtop_dn(ngrid)    ! incident top of atmosphere SW flux (W/m2)
      REAL,INTENT(OUT) :: OLR_nu(ngrid,L_NSPECTI)! Outgoing LW radition in each band (Normalized to the band width (W/m2/cm-1)
      REAL,INTENT(OUT) :: OSR_nu(ngrid,L_NSPECTV)! Outgoing SW radition in each band (Normalized to the band width (W/m2/cm-1)
      REAL,INTENT(OUT) :: tau_col(ngrid) ! diagnostic from aeropacity
!     for H2O cloud fraction in aeropacity
      real,intent(in) :: cloudfrac(ngrid,nlayer)
      real,intent(out) :: totcloudfrac(ngrid)
      logical,intent(in) :: clearsky
      logical,intent(in) :: firstcall ! signals first call to physics
      logical,intent(in) :: lastcall ! signals last call to physics

!     Globally varying aerosol optical properties on GCM grid
!     Not needed everywhere so not in radcommon_h
      REAL :: QVISsQREF3d(ngrid,nlayer,L_NSPECTV,naerkind)
      REAL :: omegaVIS3d(ngrid,nlayer,L_NSPECTV,naerkind)
      REAL :: gVIS3d(ngrid,nlayer,L_NSPECTV,naerkind)

      REAL :: QIRsQREF3d(ngrid,nlayer,L_NSPECTI,naerkind)
      REAL :: omegaIR3d(ngrid,nlayer,L_NSPECTI,naerkind)
      REAL :: gIR3d(ngrid,nlayer,L_NSPECTI,naerkind)

!      REAL :: omegaREFvis3d(ngrid,nlayer,naerkind)
!      REAL :: omegaREFir3d(ngrid,nlayer,naerkind) ! not sure of the point of these...

      REAL,ALLOCATABLE,SAVE :: reffrad(:,:,:) ! aerosol effective radius (m)
      REAL,ALLOCATABLE,SAVE :: nueffrad(:,:,:) ! aerosol effective variance
!$OMP THREADPRIVATE(reffrad,nueffrad)

!-----------------------------------------------------------------------
!     Declaration of the variables required by correlated-k subroutines
!     Numbered from top to bottom unlike in the GCM!

      REAL*8 tmid(L_LEVELS),pmid(L_LEVELS)
      REAL*8 tlevrad(L_LEVELS),plevrad(L_LEVELS)

!     Optical values for the optci/cv subroutines
      REAL*8 stel(L_NSPECTV),stel_fract(L_NSPECTV)
      REAL*8 dtaui(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      REAL*8 dtauv(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 cosbv(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 cosbi(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      REAL*8 wbari(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      REAL*8 wbarv(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 tauv(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 taucumv(L_LEVELS,L_NSPECTV,L_NGAUSS)
      REAL*8 taucumi(L_LEVELS,L_NSPECTI,L_NGAUSS)

      REAL*8 tauaero(L_LEVELS+1,naerkind)
      REAL*8 nfluxtopv,nfluxtopi,nfluxtop,fluxtopvdn
      real*8 nfluxoutv_nu(L_NSPECTV) ! outgoing band-resolved VI flux at TOA (W/m2)
      real*8 nfluxtopi_nu(L_NSPECTI) ! net band-resolved IR flux at TOA (W/m2)
      real*8 fluxupi_nu(L_NLAYRAD,L_NSPECTI) ! for 1D diagnostic
      REAL*8 fmneti(L_NLAYRAD),fmnetv(L_NLAYRAD)
      REAL*8 fluxupv(L_NLAYRAD),fluxupi(L_NLAYRAD)
      REAL*8 fluxdnv(L_NLAYRAD),fluxdni(L_NLAYRAD)
      REAL*8 albi,albv,acosz

      INTEGER ig,l,k,nw,iaer,irad
      INTEGER icount

      real szangle
      logical global1d
      save szangle,global1d
!$OMP THREADPRIVATE(szangle,global1d)
      real*8 taugsurf(L_NSPECTV,L_NGAUSS-1)
      real*8 taugsurfi(L_NSPECTI,L_NGAUSS-1)

      real*8 qvar(L_LEVELS)          ! mixing ratio of variable component (mol/mol)

!     Local aerosol optical properties for each column on RADIATIVE grid
      real*8,save ::  QXVAER(L_LEVELS+1,L_NSPECTV,naerkind)
      real*8,save ::  QSVAER(L_LEVELS+1,L_NSPECTV,naerkind)
      real*8,save ::  GVAER(L_LEVELS+1,L_NSPECTV,naerkind)
      real*8,save ::  QXIAER(L_LEVELS+1,L_NSPECTI,naerkind)
      real*8,save ::  QSIAER(L_LEVELS+1,L_NSPECTI,naerkind)
      real*8,save ::  GIAER(L_LEVELS+1,L_NSPECTI,naerkind)

      !REAL :: QREFvis3d(ngrid,nlayer,naerkind)
      !REAL :: QREFir3d(ngrid,nlayer,naerkind)
      !save QREFvis3d, QREFir3d 
      real, dimension(:,:,:), save, allocatable :: QREFvis3d
      real, dimension(:,:,:), save, allocatable :: QREFir3d
!$OMP THREADPRIVATE(QXVAER,QSVAER,GVAER,QXIAER,QSIAER,GIAER,QREFvis3d,QREFir3d)


!     Misc.
      logical nantest
      real*8  tempv(L_NSPECTV)
      real*8  tempi(L_NSPECTI)
      real*8  temp,temp1,temp2,pweight
      character(len=10) :: tmp1
      character(len=10) :: tmp2

!     for fixed water vapour profiles
      integer i_var
      real RH
      real*8 pq_temp(nlayer)
      real ptemp, Ttemp, qsat

!      real(KIND=r8) :: pq_temp(nlayer) ! better F90 way.. DOESNT PORT TO F77!!!

      !real ptime, pday
      logical OLRz
      real*8 NFLUXGNDV_nu(L_NSPECTV)


      ! for weird cloud test
      real pqtest(ngrid,nlayer,nq)

      real maxrad, minrad
            
      real,external :: CBRT

!     included by RW for runaway greenhouse 1D study
      real vtmp(nlayer)
      REAL*8 muvarrad(L_LEVELS)

!===============================================================
!     Initialization on first call

      qxvaer(:,:,:)=0.0
      qsvaer(:,:,:)=0.0
      gvaer(:,:,:) =0.0

      qxiaer(:,:,:)=0.0
      qsiaer(:,:,:)=0.0
      giaer(:,:,:) =0.0

      if(firstcall) then

         !!! ALLOCATED instances are necessary because of CLFvarying
         !!! strategy to call callcorrk twice in physiq...
         IF(.not.ALLOCATED(QREFvis3d)) ALLOCATE(QREFvis3d(ngrid,nlayer,naerkind))
         IF(.not.ALLOCATED(QREFir3d)) ALLOCATE(QREFir3d(ngrid,nlayer,naerkind))
         ! Effective radius and variance of the aerosols
         IF(.not.ALLOCATED(reffrad)) allocate(reffrad(ngrid,nlayer,naerkind))
         IF(.not.ALLOCATED(nueffrad)) allocate(nueffrad(ngrid,nlayer,naerkind))

         call system('rm -f surf_vals_long.out')

         if(naerkind.gt.4)then
            print*,'Code not general enough to deal with naerkind > 4 yet.'
            call abort
         endif
         call su_aer_radii(ngrid,nlayer,reffrad,nueffrad)
	 
	 

!--------------------------------------------------
!     set up correlated k
         print*, "callcorrk: Correlated-k data base folder:",trim(datadir)
         call getin_p("corrkdir",corrkdir)
         print*, "corrkdir = ",corrkdir
         write( tmp1, '(i3)' ) L_NSPECTI
         write( tmp2, '(i3)' ) L_NSPECTV
         banddir=trim(adjustl(tmp1))//'x'//trim(adjustl(tmp2))
         banddir=trim(adjustl(corrkdir))//'/'//trim(adjustl(banddir))

         call setspi            ! basic infrared properties
         call setspv            ! basic visible properties
         call sugas_corrk       ! set up gaseous absorption properties
         call suaer_corrk       ! set up aerosol optical properties


         if((igcm_h2o_vap.eq.0) .and. varactive)then
            print*,'varactive in callcorrk but no h2o_vap tracer.'
            stop
         endif

         OLR_nu(:,:) = 0.
         OSR_nu(:,:) = 0.

         if (ngrid.eq.1) then
           PRINT*, 'Simulate global averaged conditions ?'
           global1d = .false. ! default value
           call getin_p("global1d",global1d)
           write(*,*) "global1d = ",global1d
           ! Test of incompatibility:
           ! if global1d is true, there should not be any diurnal cycle
           if (global1d.and.diurnal) then
            print*,'if global1d is true, diurnal must be set to false'
            stop
           endif

           if (global1d) then
             PRINT *,'Solar Zenith angle (deg.) ?'
             PRINT *,'(assumed for averaged solar flux S/4)'
             szangle=60.0  ! default value
             call getin_p("szangle",szangle)
             write(*,*) "szangle = ",szangle
           endif
         endif

      end if ! of if (firstcall)

!=======================================================================
!     Initialization on every call    

!--------------------------------------------------
!     Effective radius and variance of the aerosols
      do iaer=1,naerkind

         if ((iaer.eq.iaero_co2).and.tracer.and.(igcm_co2_ice.gt.0)) then ! treat condensed co2 particles.
	    call co2_reffrad(ngrid,nlayer,nq,pq,reffrad(1,1,iaero_co2))
            print*,'Max. CO2 ice particle size = ',maxval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
            print*,'Min. CO2 ice particle size = ',minval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
	 end if
         if ((iaer.eq.iaero_h2o).and.water) then ! treat condensed water particles. to be generalized for other aerosols
	    call h2o_reffrad(ngrid,nlayer,pq(1,1,igcm_h2o_ice),pt, &
                             reffrad(1,1,iaero_h2o),nueffrad(1,1,iaero_h2o))
            print*,'Max. H2O cloud particle size = ',maxval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
            print*,'Min. H2O cloud particle size = ',minval(reffrad(1:ngrid,1:nlayer,iaer))/1.e-6,' um'
         endif
         if(iaer.eq.iaero_dust)then
	    call dust_reffrad(ngrid,nlayer,reffrad(1,1,iaero_dust))
            print*,'Dust particle size = ',reffrad(1,1,iaer)/1.e-6,' um'
         endif
         if(iaer.eq.iaero_h2so4)then
	    call h2so4_reffrad(ngrid,nlayer,reffrad(1,1,iaero_h2so4))
            print*,'H2SO4 particle size =',reffrad(1,1,iaer)/1.e-6,' um'
         endif
          if(iaer.eq.iaero_back2lay)then
	    call back2lay_reffrad(ngrid,reffrad(1,1,iaero_back2lay),nlayer,pplev)
         endif
     end do !iaer=1,naerkind


!     how much light we get
      do nw=1,L_NSPECTV
         stel(nw)=stellarf(nw)/(dist_star**2)
      end do

      call aeroptproperties(ngrid,nlayer,reffrad,nueffrad,         &
           QVISsQREF3d,omegaVIS3d,gVIS3d,                          &
           QIRsQREF3d,omegaIR3d,gIR3d,                             &
           QREFvis3d,QREFir3d)                                     ! get 3D aerosol optical properties

      call aeropacity(ngrid,nlayer,nq,pplay,pplev,pq,aerosol,      &
           reffrad,QREFvis3d,QREFir3d,                             & 
           tau_col,cloudfrac,totcloudfrac,clearsky)                ! get aerosol optical depths
	  
!-----------------------------------------------------------------------
!     Starting Big Loop over every GCM column
      do ig=1,ngrid

!=======================================================================
!     Transformation of the GCM variables

!-----------------------------------------------------------------------
!     Aerosol optical properties Qext, Qscat and g
!     The transformation in the vertical is the same as for temperature
           
!     shortwave
            do iaer=1,naerkind
               DO nw=1,L_NSPECTV 
                  do l=1,nlayer

                     temp1=QVISsQREF3d(ig,nlayer+1-l,nw,iaer)         &
                         *QREFvis3d(ig,nlayer+1-l,iaer)

                     temp2=QVISsQREF3d(ig,max(nlayer-l,1),nw,iaer)    &
                         *QREFvis3d(ig,max(nlayer-l,1),iaer)

                     qxvaer(2*l,nw,iaer)  = temp1
                     qxvaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=temp1*omegavis3d(ig,nlayer+1-l,nw,iaer)
                     temp2=temp2*omegavis3d(ig,max(nlayer-l,1),nw,iaer)

                     qsvaer(2*l,nw,iaer)  = temp1
                     qsvaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=gvis3d(ig,nlayer+1-l,nw,iaer)
                     temp2=gvis3d(ig,max(nlayer-l,1),nw,iaer)

                     gvaer(2*l,nw,iaer)  = temp1
                     gvaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                  end do

                  qxvaer(1,nw,iaer)=qxvaer(2,nw,iaer)
                  qxvaer(2*nlayer+1,nw,iaer)=0.

                  qsvaer(1,nw,iaer)=qsvaer(2,nw,iaer)
                  qsvaer(2*nlayer+1,nw,iaer)=0.

                  gvaer(1,nw,iaer)=gvaer(2,nw,iaer)
                  gvaer(2*nlayer+1,nw,iaer)=0.

               end do

!     longwave
               DO nw=1,L_NSPECTI 
                  do l=1,nlayer

                     temp1=QIRsQREF3d(ig,nlayer+1-l,nw,iaer)         &
                          *QREFir3d(ig,nlayer+1-l,iaer)

                     temp2=QIRsQREF3d(ig,max(nlayer-l,1),nw,iaer)    &
                          *QREFir3d(ig,max(nlayer-l,1),iaer)

                     qxiaer(2*l,nw,iaer)  = temp1
                     qxiaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=temp1*omegair3d(ig,nlayer+1-l,nw,iaer)
                     temp2=temp2*omegair3d(ig,max(nlayer-l,1),nw,iaer)

                     qsiaer(2*l,nw,iaer)  = temp1
                     qsiaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                     temp1=gir3d(ig,nlayer+1-l,nw,iaer)
                     temp2=gir3d(ig,max(nlayer-l,1),nw,iaer)

                     giaer(2*l,nw,iaer)  = temp1
                     giaer(2*l+1,nw,iaer)=(temp1+temp2)/2

                  end do

                  qxiaer(1,nw,iaer)=qxiaer(2,nw,iaer)
                  qxiaer(2*nlayer+1,nw,iaer)=0.

                  qsiaer(1,nw,iaer)=qsiaer(2,nw,iaer)
                  qsiaer(2*nlayer+1,nw,iaer)=0.

                  giaer(1,nw,iaer)=giaer(2,nw,iaer)
                  giaer(2*nlayer+1,nw,iaer)=0.

               end do
            end do

            ! test / correct for freaky s. s. albedo values
            do iaer=1,naerkind
               do k=1,L_LEVELS+1

                  do nw=1,L_NSPECTV
                     if(qsvaer(k,nw,iaer).gt.1.05*qxvaer(k,nw,iaer))then
                        print*,'Serious problems with qsvaer values' 
                        print*,'in callcorrk'
                        call abort
                     endif
                     if(qsvaer(k,nw,iaer).gt.qxvaer(k,nw,iaer))then
                        qsvaer(k,nw,iaer)=qxvaer(k,nw,iaer)
                     endif
                  end do

                  do nw=1,L_NSPECTI 
                     if(qsiaer(k,nw,iaer).gt.1.05*qxiaer(k,nw,iaer))then
                        print*,'Serious problems with qsiaer values'
                        print*,'in callcorrk'
                        call abort
                     endif
                     if(qsiaer(k,nw,iaer).gt.qxiaer(k,nw,iaer))then
                        qsiaer(k,nw,iaer)=qxiaer(k,nw,iaer)
                     endif
                  end do

               end do
            end do

!-----------------------------------------------------------------------
!     Aerosol optical depths
            
         do iaer=1,naerkind     ! a bug was here           
            do k=0,nlayer-1
               
               pweight=(pplay(ig,L_NLAYRAD-k)-pplev(ig,L_NLAYRAD-k+1))/   &
                        (pplev(ig,L_NLAYRAD-k)-pplev(ig,L_NLAYRAD-k+1))

               temp=aerosol(ig,L_NLAYRAD-k,iaer)/QREFvis3d(ig,L_NLAYRAD-k,iaer)

               tauaero(2*k+2,iaer)=max(temp*pweight,0.d0)
               tauaero(2*k+3,iaer)=max(temp-tauaero(2*k+2,iaer),0.d0)
!
            end do
            ! boundary conditions
            tauaero(1,iaer)          = tauaero(2,iaer)
            tauaero(L_LEVELS+1,iaer) = tauaero(L_LEVELS,iaer)
            !tauaero(1,iaer)          = 0.
            !tauaero(L_LEVELS+1,iaer) = 0.
         end do

!     Albedo and emissivity
         albi=1-emis(ig)        ! longwave
         albv=albedo(ig)        ! shortwave 

      if(nosurf.and.(albv.gt.0.0))then
         print*,'For open lower boundary in callcorrk must'
         print*,'have surface albedo set to zero!'
         call abort
      endif

      if ((ngrid.eq.1).and.(global1d)) then       ! fixed zenith angle 'szangle' in 1D simulations w/ globally-averaged sunlight
         acosz = cos(pi*szangle/180.0)
         print*,'acosz=',acosz,', szangle=',szangle
      else
         acosz=mu0(ig)          ! cosine of sun incident angle : 3D simulations or local 1D simulations using latitude
      endif

!!! JL13: in the following, I changed some indices in the interpolations so that the model rsults are less dependent on the number of layers
!!!    the older verions are commented with the commetn !JL13index 


!-----------------------------------------------------------------------
!     Water vapour (to be generalised for other gases eventually)
      
      if(varactive)then

         i_var=igcm_h2o_vap
         do l=1,nlayer
            qvar(2*l)   = pq(ig,nlayer+1-l,i_var)
            qvar(2*l+1) = pq(ig,nlayer+1-l,i_var)    
!JL13index            qvar(2*l+1) = (pq(ig,nlayer+1-l,i_var)+pq(ig,max(nlayer-l,1),i_var))/2    
!JL13index            ! Average approximation as for temperature...
         end do
         qvar(1)=qvar(2)

      elseif(varfixed)then

         do l=1,nlayer        ! here we will assign fixed water vapour profiles globally
            RH = satval * ((pplay(ig,l)/pplev(ig,1) - 0.02) / 0.98)
            if(RH.lt.0.0) RH=0.0
            
            ptemp=pplay(ig,l)
            Ttemp=pt(ig,l)
            call watersat(Ttemp,ptemp,qsat)

            !pq_temp(l) = qsat      ! fully saturated everywhere
            pq_temp(l) = RH * qsat ! ~realistic profile (e.g. 80% saturation at ground)
         end do
         
         do l=1,nlayer
            qvar(2*l)   = pq_temp(nlayer+1-l)
            qvar(2*l+1) = (pq_temp(nlayer+1-l)+pq_temp(max(nlayer-l,1)))/2
         end do
         qvar(1)=qvar(2)

         ! Lowest layer of atmosphere
         RH = satval * (1 - 0.02) / 0.98
         if(RH.lt.0.0) RH=0.0

!         ptemp = pplev(ig,1)
!         Ttemp = tsurf(ig)
!         call watersat(Ttemp,ptemp,qsat)

         qvar(2*nlayer+1)= RH * qsat ! ~realistic profile (e.g. 80% saturation at ground)
 
      else
         do k=1,L_LEVELS
            qvar(k) = 1.0D-7
         end do
      end if

      if(.not.kastprof)then
      ! IMPORTANT: Now convert from kg/kg to mol/mol
         do k=1,L_LEVELS
            qvar(k) = qvar(k)/(epsi+qvar(k)*(1.-epsi))
         end do
      end if

!-----------------------------------------------------------------------
!     kcm mode only
      if(kastprof)then

         ! initial values equivalent to mugaz
         DO l=1,nlayer
            muvarrad(2*l)   = mugaz
            muvarrad(2*l+1) = mugaz
         END DO

         if(ngasmx.gt.1)then

            DO l=1,nlayer
               muvarrad(2*l)   = muvar(ig,nlayer+2-l)
               muvarrad(2*l+1) = (muvar(ig,nlayer+2-l) + &
                                muvar(ig,max(nlayer+1-l,1)))/2
            END DO
      
            muvarrad(1) = muvarrad(2)
            muvarrad(2*nlayer+1)=muvar(ig,1)

            print*,'Recalculating qvar with VARIABLE epsi for kastprof'
            print*,'Assumes that the variable gas is H2O!!!'
            print*,'Assumes that there is only one tracer'
            !i_var=igcm_h2o_vap
            i_var=1
            if(nq.gt.1)then
               print*,'Need 1 tracer only to run kcm1d.e' 
               stop
            endif
            do l=1,nlayer
               vtmp(l)=pq(ig,l,i_var)/(epsi+pq(ig,l,i_var)*(1.-epsi)) 
               !vtmp(l)=pq(ig,l,i_var)*muvar(ig,l+1)/mH2O !JL to be changed
            end do

            do l=1,nlayer
               qvar(2*l)   = vtmp(nlayer+1-l)
               qvar(2*l+1) = vtmp(nlayer+1-l)
!               qvar(2*l+1) = ( vtmp(nlayer+1-l) + vtmp(max(nlayer-l,1)) )/2
            end do
            qvar(1)=qvar(2)

            print*,'Warning: reducing qvar in callcorrk.'
            print*,'Temperature profile no longer consistent ', &
                            'with saturated H2O. qsat=',satval
            do k=1,L_LEVELS
               qvar(k) = qvar(k)*satval
            end do

         endif
      else ! if kastprof
         DO l=1,nlayer
            muvarrad(2*l)   = muvar(ig,nlayer+2-l)
            muvarrad(2*l+1) = (muvar(ig,nlayer+2-l)+muvar(ig,max(nlayer+1-l,1)))/2
         END DO
      
         muvarrad(1) = muvarrad(2)
         muvarrad(2*nlayer+1)=muvar(ig,1)         
      endif
      
      ! Keep values inside limits for which we have radiative transfer coefficients
      if(L_REFVAR.gt.1)then ! there was a bug here!
         do k=1,L_LEVELS
            if(qvar(k).lt.wrefvar(1))then
               qvar(k)=wrefvar(1)+1.0e-8
            elseif(qvar(k).gt.wrefvar(L_REFVAR))then
               qvar(k)=wrefvar(L_REFVAR)-1.0e-8
            endif
         end do
      endif

!-----------------------------------------------------------------------
!     Pressure and temperature

      DO l=1,nlayer
         plevrad(2*l)   = pplay(ig,nlayer+1-l)/scalep
         plevrad(2*l+1) = pplev(ig,nlayer+1-l)/scalep
         tlevrad(2*l)   = pt(ig,nlayer+1-l)
         tlevrad(2*l+1) = (pt(ig,nlayer+1-l)+pt(ig,max(nlayer-l,1)))/2
      END DO
      
      plevrad(1) = 0.
      plevrad(2) = max(pgasmin,0.0001*plevrad(3))

      tlevrad(1) = tlevrad(2)
      tlevrad(2*nlayer+1)=tsurf(ig)
      
      tmid(1) = tlevrad(2)
      tmid(2) = tlevrad(2)
      pmid(1) = plevrad(2)
      pmid(2) = plevrad(2)
      
      DO l=1,L_NLAYRAD-1
         tmid(2*l+1) = tlevrad(2*l)
         tmid(2*l+2) = tlevrad(2*l)
         pmid(2*l+1) = plevrad(2*l)
         pmid(2*l+2) = plevrad(2*l)
!JL13index         tmid(2*l+1) = tlevrad(2*l+1)
!JL13index         tmid(2*l+2) = tlevrad(2*l+1)
!JL13index         pmid(2*l+1) = plevrad(2*l+1)
!JL13index         pmid(2*l+2) = plevrad(2*l+1)
      END DO
      pmid(L_LEVELS) = plevrad(L_LEVELS-1)
      tmid(L_LEVELS) = tlevrad(L_LEVELS-1)
!JL13index      pmid(L_LEVELS) = plevrad(L_LEVELS)
!JL13index      tmid(L_LEVELS) = tlevrad(L_LEVELS)

      ! test for out-of-bounds pressure
      if(plevrad(3).lt.pgasmin)then
         print*,'Minimum pressure is outside the radiative'
         print*,'transfer kmatrix bounds, exiting.'
         call abort
      elseif(plevrad(L_LEVELS).gt.pgasmax)then
         print*,'Maximum pressure is outside the radiative'
         print*,'transfer kmatrix bounds, exiting.'
         call abort
      endif

      ! test for out-of-bounds temperature
      do k=1,L_LEVELS
         if(tlevrad(k).lt.tgasmin)then
            print*,'Minimum temperature is outside the radiative'
            print*,'transfer kmatrix bounds'
            print*,"k=",k," tlevrad(k)=",tlevrad(k)
            print*,"tgasmin=",tgasmin
            if (strictboundcorrk) then
              call abort
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tlevrad=tgasmin' 
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              !tlevrad(k)=tgasmin
            endif
         elseif(tlevrad(k).gt.tgasmax)then
            print*,'Maximum temperature is outside the radiative'
            print*,'transfer kmatrix bounds, exiting.'
            print*,"k=",k," tlevrad(k)=",tlevrad(k)
            print*,"tgasmax=",tgasmax
            if (strictboundcorrk) then
              call abort
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tlevrad=tgasmax'  
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              !tlevrad(k)=tgasmax
            endif
         endif
      enddo
      do k=1,L_NLAYRAD+1
         if(tmid(k).lt.tgasmin)then
            print*,'Minimum temperature is outside the radiative'
            print*,'transfer kmatrix bounds, exiting.'
            print*,"k=",k," tmid(k)=",tmid(k)
            print*,"tgasmin=",tgasmin
            if (strictboundcorrk) then
              call abort
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tmid=tgasmin'
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              tmid(k)=tgasmin
            endif
         elseif(tmid(k).gt.tgasmax)then
            print*,'Maximum temperature is outside the radiative'
            print*,'transfer kmatrix bounds, exiting.'
            print*,"k=",k," tmid(k)=",tmid(k)
            print*,"tgasmax=",tgasmax
            if (strictboundcorrk) then
              call abort
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tmid=tgasmin'
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              tmid(k)=tgasmax
            endif
         endif
      enddo

!=======================================================================
!     Calling the main radiative transfer subroutines


         Cmk= 0.01 * 1.0 / (glat(ig) * mugaz * 1.672621e-27) ! q_main=1.0 assumed
	 glat_ig=glat(ig)

!-----------------------------------------------------------------------
!     Shortwave

         if(fract(ig) .ge. 1.0e-4) then ! only during daylight!
            if((ngrid.eq.1).and.(global1d))then
               do nw=1,L_NSPECTV
                  stel_fract(nw)= stel(nw)* 0.25 / acosz
                                ! globally averaged = divide by 4
                                ! but we correct for solar zenith angle
               end do
            else
               do nw=1,L_NSPECTV
                  stel_fract(nw)= stel(nw) * fract(ig)
               end do
            endif 
            call optcv(dtauv,tauv,taucumv,plevrad,                 &
                 qxvaer,qsvaer,gvaer,wbarv,cosbv,tauray,tauaero,   &
                 tmid,pmid,taugsurf,qvar,muvarrad)

            call sfluxv(dtauv,tauv,taucumv,albv,dwnv,wbarv,cosbv,  &
                 acosz,stel_fract,gweight,                         &
                 nfluxtopv,fluxtopvdn,nfluxoutv_nu,nfluxgndv_nu,              &
                 fmnetv,fluxupv,fluxdnv,fzerov,taugsurf)

         else                          ! during the night, fluxes = 0
            nfluxtopv       = 0.0d0
	    fluxtopvdn      = 0.0d0
            nfluxoutv_nu(:) = 0.0d0
            nfluxgndv_nu(:) = 0.0d0
            do l=1,L_NLAYRAD
               fmnetv(l)=0.0d0
               fluxupv(l)=0.0d0
               fluxdnv(l)=0.0d0
            end do
         end if

!-----------------------------------------------------------------------
!     Longwave

         call optci(plevrad,tlevrad,dtaui,taucumi,                  &
              qxiaer,qsiaer,giaer,cosbi,wbari,tauaero,tmid,pmid,    &
              taugsurfi,qvar,muvarrad)

         call sfluxi(plevrad,tlevrad,dtaui,taucumi,ubari,albi,      &
              wnoi,dwni,cosbi,wbari,gweight,nfluxtopi,nfluxtopi_nu, & 
              fmneti,fluxupi,fluxdni,fluxupi_nu,fzeroi,taugsurfi)

!-----------------------------------------------------------------------
!     Transformation of the correlated-k code outputs
!     (into dtlw, dtsw, fluxsurf_lw, fluxsurf_sw, fluxtop_lw, fluxtop_sw)

!     Flux incident at the top of the atmosphere
         fluxtop_dn(ig)=fluxtopvdn 

         fluxtop_lw(ig)  = real(nfluxtopi)
         fluxabs_sw(ig)  = real(-nfluxtopv)
         fluxsurf_lw(ig) = real(fluxdni(L_NLAYRAD))
         fluxsurf_sw(ig) = real(fluxdnv(L_NLAYRAD))

         if(fluxtop_dn(ig).lt.0.0)then
            print*,'Achtung! fluxtop_dn has lost the plot!'
            print*,'fluxtop_dn=',fluxtop_dn(ig)
            print*,'acosz=',acosz
            print*,'aerosol=',aerosol(ig,:,:)
            print*,'temp=   ',pt(ig,:)
            print*,'pplay=  ',pplay(ig,:)
            call abort
         endif

!     Spectral output, for exoplanet observational comparison
         if(specOLR)then
            do nw=1,L_NSPECTI 
               OLR_nu(ig,nw)=nfluxtopi_nu(nw)/DWNI(nw) !JL Normalize to the bandwidth
            end do
            do nw=1,L_NSPECTV 
               !GSR_nu(ig,nw)=nfluxgndv_nu(nw)
               OSR_nu(ig,nw)=nfluxoutv_nu(nw)/DWNV(nw) !JL Normalize to the bandwidth
            end do
         endif

!     Finally, the heating rates

         DO l=2,L_NLAYRAD
            dtsw(ig,L_NLAYRAD+1-l)=(fmnetv(l)-fmnetv(l-1))  &
                *glat(ig)/(cpp*scalep*(plevrad(2*l+1)-plevrad(2*l-1)))
            dtlw(ig,L_NLAYRAD+1-l)=(fmneti(l)-fmneti(l-1))  &
                *glat(ig)/(cpp*scalep*(plevrad(2*l+1)-plevrad(2*l-1)))
         END DO      

!     These are values at top of atmosphere
         dtsw(ig,L_NLAYRAD)=(fmnetv(1)-nfluxtopv)           &
             *glat(ig)/(cpp*scalep*(plevrad(3)-plevrad(1)))
         dtlw(ig,L_NLAYRAD)=(fmneti(1)-nfluxtopi)           &
             *glat(ig)/(cpp*scalep*(plevrad(3)-plevrad(1)))

!     ---------------------------------------------------------------
      end do                    ! end of big loop over every GCM column (ig = 1:ngrid)


!-----------------------------------------------------------------------
!     Additional diagnostics

!     IR spectral output, for exoplanet observational comparison


      if(lastcall.and.(ngrid.eq.1))then  ! could disable the 1D output, they are in the diagfi and diagspec... JL12

           print*,'Saving scalar quantities in surf_vals.out...'
           print*,'psurf = ', pplev(1,1),' Pa'
           open(116,file='surf_vals.out')
           write(116,*) tsurf(1),pplev(1,1),fluxtop_dn(1),         &
                real(-nfluxtopv),real(nfluxtopi) 
           close(116)

!          I am useful, please don`t remove me!
!           if(specOLR)then
!               open(117,file='OLRnu.out')
!               do nw=1,L_NSPECTI
!                  write(117,*) OLR_nu(1,nw)
!               enddo
!               close(117)
!
!               open(127,file='OSRnu.out')
!               do nw=1,L_NSPECTV
!                  write(127,*) OSR_nu(1,nw)
!               enddo
!               close(127)
!           endif

!     OLR vs altitude: do it as a .txt file
           OLRz=.false.
           if(OLRz)then
              print*,'saving IR vertical flux for OLRz...'
              open(118,file='OLRz_plevs.out')
              open(119,file='OLRz.out')
              do l=1,L_NLAYRAD
                 write(118,*) plevrad(2*l)
                 do nw=1,L_NSPECTI
                     write(119,*) fluxupi_nu(l,nw) 
                  enddo
              enddo 
              close(118)
              close(119)
           endif

      endif

      ! see physiq.F for explanations about CLFvarying. This is temporary.
      if (lastcall .and. .not.CLFvarying) then
        IF( ALLOCATED( gasi ) ) DEALLOCATE( gasi )
        IF( ALLOCATED( gasv ) ) DEALLOCATE( gasv )
!$OMP BARRIER
!$OMP MASTER
        IF( ALLOCATED( pgasref ) ) DEALLOCATE( pgasref )
        IF( ALLOCATED( tgasref ) ) DEALLOCATE( tgasref )
        IF( ALLOCATED( wrefvar ) ) DEALLOCATE( wrefvar )
        IF( ALLOCATED( pfgasref ) ) DEALLOCATE( pfgasref )
!$OMP END MASTER
!$OMP BARRIER	
        IF ( ALLOCATED(reffrad)) DEALLOCATE(reffrad)
        IF ( ALLOCATED(nueffrad)) DEALLOCATE(nueffrad)
      endif


    end subroutine callcorrk
