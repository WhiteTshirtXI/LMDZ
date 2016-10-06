      Subroutine aeropacity(ngrid,nlayer,nq,pplay,pplev,pq, &
         aerosol,reffrad,QREFvis3d,QREFir3d,tau_col,cloudfrac,totcloudfrac,clearsky)

       use radinc_h, only : L_TAUMAX,naerkind
       use aerosol_mod
       USE comgeomfi_h
       USE tracer_h, only: noms,rho_co2,rho_ice
                  
       implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Compute aerosol optical depth in each gridbox.
!     
!     Authors
!     ------- 
!     F. Forget
!     F. Montmessin (water ice scheme) 
!     update J.-B. Madeleine (2008)
!     dust removal, simplification by Robin Wordsworth (2009)
!
!     Input
!     ----- 
!     ngrid             Number of horizontal gridpoints
!     nlayer            Number of layers
!     nq                Number of tracers
!     pplev             Pressure (Pa) at each layer boundary
!     pq                Aerosol mixing ratio
!     reffrad(ngrid,nlayer,naerkind)         Aerosol effective radius
!     QREFvis3d(ngrid,nlayer,naerkind) \ 3d extinction coefficients
!     QREFir3d(ngrid,nlayer,naerkind)  / at reference wavelengths
!
!     Output
!     ------
!     aerosol            Aerosol optical depth in layer l, grid point ig
!     tau_col            Total column optical depth at grid point ig
!
!=======================================================================

!#include "dimensions.h"
!#include "dimphys.h"
#include "callkeys.h"
#include "comcstfi.h"
!#include "comvert.h"

      INTEGER,INTENT(IN) :: ngrid  ! number of atmospheric columns
      INTEGER,INTENT(IN) :: nlayer ! number of atmospheric layers
      INTEGER,INTENT(IN) :: nq     ! number of tracers
      REAL,INTENT(IN) :: pplay(ngrid,nlayer) ! mid-layer pressure (Pa)
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq) ! tracers (.../kg_of_air)
      REAL,INTENT(OUT) :: aerosol(ngrid,nlayer,naerkind) ! aerosol optical depth
      REAL,INTENT(IN) :: reffrad(ngrid,nlayer,naerkind) ! aerosol effective radius
      REAL,INTENT(IN) :: QREFvis3d(ngrid,nlayer,naerkind) ! extinction coefficient in the visible
      REAL,INTENT(IN) :: QREFir3d(ngrid,nlayer,naerkind)
      REAL,INTENT(OUT):: tau_col(ngrid) !column integrated visible optical depth
      ! BENJAMIN MODIFS
      real,intent(in) :: cloudfrac(ngrid,nlayer) ! cloud fraction
      real,intent(out) :: totcloudfrac(ngrid) ! total cloud fraction
      logical,intent(in) :: clearsky

      real aerosol0

      INTEGER l,ig,iq,iaer

      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)
      REAL CBRT
      EXTERNAL CBRT

      INTEGER,SAVE :: i_co2ice=0      ! co2 ice
      INTEGER,SAVE :: i_h2oice=0      ! water ice
!$OMP THREADPRIVATE(i_co2ice,i_h2oice)
      CHARACTER(LEN=20) :: tracername ! to temporarily store text

      ! for fixed dust profiles
      real topdust, expfactor, zp
      REAL taudusttmp(ngrid) ! Temporary dust opacity used before scaling
      REAL tauh2so4tmp(ngrid) ! Temporary h2so4 opacity used before scaling

      real CLFtot

      ! identify tracers
      IF (firstcall) THEN

        write(*,*) "Tracers found in aeropacity:"
        do iq=1,nq
          tracername=noms(iq)
          if (tracername.eq."co2_ice") then
            i_co2ice=iq
          write(*,*) "i_co2ice=",i_co2ice

          endif
          if (tracername.eq."h2o_ice") then
            i_h2oice=iq
            write(*,*) "i_h2oice=",i_h2oice
          endif
        enddo

        if (noaero) then
          print*, "No active aerosols found in aeropacity"
        else
          print*, "If you would like to use aerosols, make sure any old"
          print*, "start files are updated in newstart using the option"
          print*, "q=0"
          write(*,*) "Active aerosols found in aeropacity:"
        endif

        if ((iaero_co2.ne.0).and.(.not.noaero)) then
          print*, 'iaero_co2=  ',iaero_co2
        endif
        if (iaero_h2o.ne.0) then
          print*,'iaero_h2o=  ',iaero_h2o    
        endif
        if (iaero_dust.ne.0) then
          print*,'iaero_dust= ',iaero_dust
        endif
        if (iaero_h2so4.ne.0) then
          print*,'iaero_h2so4= ',iaero_h2so4
        endif
        if (iaero_back2lay.ne.0) then
          print*,'iaero_back2lay= ',iaero_back2lay
        endif

        firstcall=.false.
      ENDIF ! of IF (firstcall)


!     ---------------------------------------------------------
!==================================================================
!    CO2 ice aerosols
!==================================================================

      if (iaero_co2.ne.0) then
           iaer=iaero_co2
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.0
!       2. Opacity calculation
            if (noaero) then ! aerosol set to zero
             aerosol(1:ngrid,1:nlayer,iaer)=0.0
            elseif (aerofixco2.or.(i_co2ice.eq.0)) then !  CO2 ice cloud prescribed
               aerosol(1:ngrid,1:nlayer,iaer)=1.e-9
               !aerosol(1:ngrid,12,iaer)=4.0 ! single cloud layer option
            else
               DO ig=1, ngrid
                  DO l=1,nlayer-1 ! to stop the rad tran bug

                     aerosol0 =                         &
                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
                          ( rho_co2 * reffrad(ig,l,iaer) )  ) *   &
                          ( pq(ig,l,i_co2ice) + 1.E-9 ) *         &
                          ( pplev(ig,l) - pplev(ig,l+1) ) / g
                     aerosol0           = max(aerosol0,1.e-9)
                     aerosol0           = min(aerosol0,L_TAUMAX)
                     aerosol(ig,l,iaer) = aerosol0
!                     aerosol(ig,l,iaer) = 0.0
!                     print*, aerosol(ig,l,iaer)
!        using cloud fraction
!                     aerosol(ig,l,iaer) = -log(1 - CLF + CLF*exp(-aerosol0/CLF))
!                     aerosol(ig,l,iaer) = min(aerosol(ig,l,iaer),L_TAUMAX)


                  ENDDO
               ENDDO
            end if ! if fixed or varying
      end if ! if CO2 aerosols   
!==================================================================
!     Water ice / liquid 
!==================================================================

      if (iaero_h2o.ne.0) then 
           iaer=iaero_h2o
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.0
!       2. Opacity calculation
            if (aerofixh2o.or.(i_h2oice.eq.0).or.clearsky) then
               aerosol(1:ngrid,1:nlayer,iaer) =1.e-9

               ! put cloud at cloudlvl
               if(kastprof.and.(cloudlvl.ne.0.0))then
                  ig=1
                  do l=1,nlayer
                     if(int(cloudlvl).eq.l)then
                     !if(cloudlvl.gt.(pplay(ig,l)/pplev(ig,1)))then
                        print*,'Inserting cloud at level ',l
                        !aerosol(ig,l,iaer)=10.0

                        rho_ice=920.0

                        ! the Kasting approximation
                        aerosol(ig,l,iaer) =                      &
                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
                          ( rho_ice * reffrad(ig,l,iaer) )  ) *   &
                          !( pq(ig,l,i_h2oice) + 1.E-9 ) *         &
                          ( 4.0e-4 + 1.E-9 ) *         &
                          ( pplev(ig,l) - pplev(ig,l+1) ) / g


                        open(115,file='clouds.out',form='formatted')
                        write(115,*) l,aerosol(ig,l,iaer)
                        close(115)

                        return
                     endif
                  end do

                  call abort
               endif

            else

               do ig=1, ngrid
                  do l=1,nlayer-1 ! to stop the rad tran bug

                     aerosol(ig,l,iaer) =                                    & !modification by BC
                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
                          ( rho_ice * reffrad(ig,l,iaer) )  ) *   &
                          !  pq(ig,l,i_h2oice) *                   & !JL I dropped the +1e-9 here to have the same
                          !( pplev(ig,l) - pplev(ig,l+1) ) / g       !   opacity in the clearsky=true and the 
                                                                     !   clear=false/pq=0 case
                          ( pq(ig,l,i_h2oice) + 1.E-9 ) *         & ! Doing this makes the code unstable, so I have restored it (RW)
                          ( pplev(ig,l) - pplev(ig,l+1) ) / g 

                  enddo
               enddo

               if(CLFvarying)then
                  call totalcloudfrac(ngrid,nlayer,nq,cloudfrac,totcloudfrac,pplev,pq,aerosol(1,1,iaer))
                  do ig=1, ngrid
                     do l=1,nlayer-1 ! to stop the rad tran bug
                        CLFtot  = max(totcloudfrac(ig),0.01)
                        aerosol(ig,l,iaer)=aerosol(ig,l,iaer)/CLFtot
                        aerosol(ig,l,iaer) = max(aerosol(ig,l,iaer),1.e-9)
                     enddo
                  enddo
               else
                  do ig=1, ngrid
                     do l=1,nlayer-1 ! to stop the rad tran bug
                        CLFtot  = CLFfixval
                        aerosol(ig,l,iaer)=aerosol(ig,l,iaer)/CLFtot
                        aerosol(ig,l,iaer) = max(aerosol(ig,l,iaer),1.e-9)
                     enddo
                  enddo
              end if!(CLFvarying)
            endif !(aerofixed.or.(i_h2oice.eq.0).or.clearsky)
	      
      end if ! End if h2o aerosol

!==================================================================
!             Dust 
!==================================================================
      if (iaero_dust.ne.0) then
          iaer=iaero_dust
!         1. Initialization 
          aerosol(1:ngrid,1:nlayer,iaer)=0.0
          
          topdust=30.0 ! km  (used to be 10.0 km) LK

!       2. Opacity calculation

!           expfactor=0.
           DO l=1,nlayer-1
             DO ig=1,ngrid
!             Typical mixing ratio profile

                 zp=(pplev(ig,1)/pplay(ig,l))**(70./topdust)
                 expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)

!             Vertical scaling function
              aerosol(ig,l,iaer)= (pplev(ig,l)-pplev(ig,l+1)) &
               *expfactor


             ENDDO
           ENDDO

!          Rescaling each layer to reproduce the choosen (or assimilated)
!          dust extinction opacity at visible reference wavelength, which
!          is scaled to the surface pressure pplev(ig,1)

            taudusttmp(1:ngrid)=0.
              DO l=1,nlayer
                DO ig=1,ngrid
                   taudusttmp(ig) = taudusttmp(ig) &
                          +  aerosol(ig,l,iaer)
                ENDDO
              ENDDO
            DO l=1,nlayer-1
               DO ig=1,ngrid
                  aerosol(ig,l,iaer) = max(1E-20, &
                          dusttau &
                       *  pplev(ig,1) / pplev(ig,1) &
                       *  aerosol(ig,l,iaer) &
                       /  taudusttmp(ig))

              ENDDO
            ENDDO
      end if ! If dust aerosol   

!==================================================================
!           H2SO4 
!==================================================================
! added by LK
      if (iaero_h2so4.ne.0) then
         iaer=iaero_h2so4

!       1. Initialization
         aerosol(1:ngrid,1:nlayer,iaer)=0.0


!       2. Opacity calculation

!           expfactor=0.
         DO l=1,nlayer-1
            DO ig=1,ngrid
!              Typical mixing ratio profile

               zp=(pplev(ig,1)/pplay(ig,l))**(70./30) !emulating topdust
               expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)

!             Vertical scaling function
               aerosol(ig,l,iaer)= (pplev(ig,l)-pplev(ig,l+1))*expfactor

            ENDDO
         ENDDO
         tauh2so4tmp(1:ngrid)=0.
         DO l=1,nlayer
            DO ig=1,ngrid
               tauh2so4tmp(ig) = tauh2so4tmp(ig) + aerosol(ig,l,iaer)
            ENDDO
         ENDDO
         DO l=1,nlayer-1
            DO ig=1,ngrid
               aerosol(ig,l,iaer) = max(1E-20, &
                          1 &
                       *  pplev(ig,1) / pplev(ig,1) &
                       *  aerosol(ig,l,iaer) &
                       /  tauh2so4tmp(ig))

            ENDDO
         ENDDO

! 1/700. is assuming a "sulfurtau" of 1
! Sulfur aerosol routine to be improved.
!                     aerosol0 =                         &
!                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
!                          ( rho_h2so4 * reffrad(ig,l,iaer) )  ) *   &
!                          ( pq(ig,l,i_h2so4) + 1.E-9 ) *         &
!                          ( pplev(ig,l) - pplev(ig,l+1) ) / g
!                     aerosol0           = max(aerosol0,1.e-9)
!                     aerosol0           = min(aerosol0,L_TAUMAX)
!                     aerosol(ig,l,iaer) = aerosol0

!                  ENDDO
!               ENDDO
      end if
 
           
!     ---------------------------------------------------------
!==================================================================
!    Two-layer aerosols (unknown composition)
!    S. Guerlet (2013)
!==================================================================

      if (iaero_back2lay .ne.0) then
           iaer=iaero_back2lay
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.0
!       2. Opacity calculation
          DO ig=1,ngrid
           DO l=1,nlayer-1
             aerosol(ig,l,iaer) = ( pplev(ig,l) - pplev(ig,l+1) )
             !! 1. below tropospheric layer: no aerosols
             IF (pplev(ig,l) .gt. pres_bottom_tropo) THEN
               aerosol(ig,l,iaer) = 0.*aerosol(ig,l,iaer)
             !! 2. tropo layer
             ELSEIF (pplev(ig,l) .le. pres_bottom_tropo .and. pplev(ig,l) .ge. pres_top_tropo) THEN
               aerosol(ig,l,iaer) = obs_tau_col_tropo*aerosol(ig,l,iaer)
             !! 3. linear transition
             ELSEIF (pplev(ig,l) .lt. pres_top_tropo .and. pplev(ig,l) .gt. pres_bottom_strato) THEN
               expfactor=log(obs_tau_col_strato/obs_tau_col_tropo)/log(pres_bottom_strato/pres_top_tropo)
               aerosol(ig,l,iaer)= obs_tau_col_tropo*((pplev(ig,l)/pres_top_tropo)**expfactor)*aerosol(ig,l,iaer)/1.5
             !! 4. strato layer
             ELSEIF (pplev(ig,l) .le. pres_bottom_strato .and. pplev(ig,l) .gt. pres_top_strato) THEN
               aerosol(ig,l,iaer)= obs_tau_col_strato*aerosol(ig,l,iaer)
             !! 5. above strato layer: no aerosols
             ELSEIF (pplev(ig,l) .lt. pres_top_strato) THEN
               aerosol(ig,l,iaer) = 0.*aerosol(ig,l,iaer)
             ENDIF
	   ENDDO
          ENDDO

 !       3. Re-normalize to observed total column
         tau_col(:)=0.0
         DO l=1,nlayer
          DO ig=1,ngrid
               tau_col(ig) = tau_col(ig) &
                     + aerosol(ig,l,iaer)/(obs_tau_col_tropo+obs_tau_col_strato)
            ENDDO
         ENDDO

         DO ig=1,ngrid
           DO l=1,nlayer-1
                aerosol(ig,l,iaer)=aerosol(ig,l,iaer)/tau_col(ig)
           ENDDO
         ENDDO


      end if ! if Two-layer aerosols  


! --------------------------------------------------------------------------
! Column integrated visible optical depth in each point (used for diagnostic)

      tau_col(:)=0.0
      do iaer = 1, naerkind
         do l=1,nlayer
            do ig=1,ngrid
               tau_col(ig) = tau_col(ig) + aerosol(ig,l,iaer)
            end do
         end do
      end do

      do ig=1,ngrid
         do l=1,nlayer
            do iaer = 1, naerkind
               if(aerosol(ig,l,iaer).gt.1.e3)then
                  print*,'WARNING: aerosol=',aerosol(ig,l,iaer)
                  print*,'at ig=',ig,',  l=',l,', iaer=',iaer
                  print*,'QREFvis3d=',QREFvis3d(ig,l,iaer)
                  print*,'reffrad=',reffrad(ig,l,iaer)
               endif
            end do
         end do
      end do

      do ig=1,ngrid
         if(tau_col(ig).gt.1.e3)then
            print*,'WARNING: tau_col=',tau_col(ig)
            print*,'at ig=',ig
            print*,'aerosol=',aerosol(ig,:,:)
            print*,'QREFvis3d=',QREFvis3d(ig,:,:)
            print*,'reffrad=',reffrad(ig,:,:)
         endif
      end do
      return
    end subroutine aeropacity
      
