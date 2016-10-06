      SUBROUTINE inifis(ngrid,nlayer,nq,
     $           day_ini,pdaysec,ptimestep,
     $           plat,plon,parea,
     $           prad,pg,pr,pcpp)

      use radinc_h, only : naerkind
      use datafile_mod, only: datadir
      use comdiurn_h, only: sinlat, coslat, sinlon, coslon
      use comgeomfi_h, only: long, lati, area, totarea, totarea_planet
      use comsoil_h, only: ini_comsoil_h
      use control_mod, only: ecritphy
      use planete_mod, only: nres
      use planetwide_mod, only: planetwide_sumval

!=======================================================================
!
!   purpose:
!   -------
!
!   Initialisation for the physical parametrisations of the LMD 
!   Generic Model.
!
!   author: Frederic Hourdin 15 / 10 /93
!   -------
!   modified: Sebastien Lebonnois 11/06/2003 (new callphys.def)
!             Ehouarn Millour (oct. 2008) tracers are now identified
!              by their names and may not be contiguously
!              stored in the q(:,:,:,:) array
!             E.M. (june 2009) use getin routine to load parameters
!
!
!   arguments:
!   ----------
!
!   input:
!   ------
!
!    ngrid                 Size of the horizontal grid.
!                          All internal loops are performed on that grid.
!    nlayer                Number of vertical layers.
!    pdayref               Day of reference for the simulation
!    pday                  Number of days counted from the North. Spring
!                          equinoxe.
!
!=======================================================================
!
!-----------------------------------------------------------------------
!   declarations:
!   -------------
      use datafile_mod, only: datadir
! to use  'getin'
!      USE ioipsl_getincom 
      USE ioipsl_getincom_p
      IMPLICIT NONE
!#include "dimensions.h"
!#include "dimphys.h"
!#include "planete.h"
#include "comcstfi.h"
#include "callkeys.h"



      REAL,INTENT(IN) :: prad,pg,pr,pcpp,pdaysec,ptimestep
 
      INTEGER,INTENT(IN) :: ngrid,nlayer,nq
      REAL,INTENT(IN) :: plat(ngrid),plon(ngrid),parea(ngrid)
      integer day_ini
      INTEGER ig,ierr
 
      EXTERNAL iniorbit,orbite
      EXTERNAL SSUM
      REAL SSUM
 
      CHARACTER ch1*12
      CHARACTER ch80*80

      logical chem, h2o
      logical :: parameter, doubleq=.false.

      real psurf,pN2 ! added by RW for Gliese 581d N2+CO2

!$OMP MASTER
      rad=prad
      daysec=pdaysec
      dtphys=ptimestep
      cpp=pcpp
      g=pg
      r=pr
      rcp=r/cpp

      avocado = 6.02214179e23   ! added by RW

!$OMP END MASTER
!$OMP BARRIER

      ! read in 'ecritphy' (frequency of calls to physics, in dynamical steps)
      ! (also done in dyn3d/defrun_new but not in LMDZ.COMMON)
      call getin_p("ecritphy",ecritphy)

! --------------------------------------------------------------
!  Reading the "callphys.def" file controlling some key options
! --------------------------------------------------------------
     
!$OMP MASTER     
      ! check that 'callphys.def' file is around
      OPEN(99,file='callphys.def',status='old',form='formatted'
     &     ,iostat=ierr)
      CLOSE(99)
      IF(ierr.EQ.0) iscallphys=.true. !iscallphys initialised as false in callkeys.h
!$OMP END MASTER
!$OMP BARRIER
      
!!!      IF(ierr.EQ.0) THEN
      IF(iscallphys) THEN
         PRINT*
         PRINT*
         PRINT*,'--------------------------------------------'
         PRINT*,' inifis: Parametres pour la physique (callphys.def)'
         PRINT*,'--------------------------------------------'

         write(*,*) "Directory where external input files are:"
         ! default 'datadir' is set in "datadir_mod"
         call getin_p("datadir",datadir) ! default path
         write(*,*) " datadir = ",trim(datadir)

         write(*,*) "Run with or without tracer transport ?"
         tracer=.false. ! default value
         call getin_p("tracer",tracer)
         write(*,*) " tracer = ",tracer

         write(*,*) "Run with or without atm mass update ",
     &      " due to tracer evaporation/condensation?"
         mass_redistrib=.false. ! default value
         call getin_p("mass_redistrib",mass_redistrib)
         write(*,*) " mass_redistrib = ",mass_redistrib

         write(*,*) "Diurnal cycle ?"
         write(*,*) "(if diurnal=false, diurnal averaged solar heating)"
         diurnal=.true. ! default value
         call getin_p("diurnal",diurnal)
         write(*,*) " diurnal = ",diurnal

         write(*,*) "Seasonal cycle ?"
         write(*,*) "(if season=false, Ls stays constant, to value ",
     &   "set in 'start'"
         season=.true. ! default value
         call getin_p("season",season)
         write(*,*) " season = ",season

         write(*,*) "Tidally resonant rotation ?"
         tlocked=.false. ! default value
         call getin_p("tlocked",tlocked)
         write(*,*) "tlocked = ",tlocked

         write(*,*) "Saturn ring shadowing ?"
         rings_shadow = .false.
         call getin_p("rings_shadow", rings_shadow)
         write(*,*) "rings_shadow = ", rings_shadow
         
         write(*,*) "Compute latitude-dependent gravity field?"
         oblate = .false.
         call getin_p("oblate", oblate)
         write(*,*) "oblate = ", oblate

         write(*,*) "Flattening of the planet (a-b)/a "
         flatten = 0.0
         call getin_p("flatten", flatten)
         write(*,*) "flatten = ", flatten
         

         write(*,*) "Needed if oblate=.true.: J2"
         J2 = 0.0
         call getin_p("J2", J2)
         write(*,*) "J2 = ", J2
         
         write(*,*) "Needed if oblate=.true.: Planet mass (*1e24 kg)"
         MassPlanet = 0.0
         call getin_p("MassPlanet", MassPlanet)
         write(*,*) "MassPlanet = ", MassPlanet         

         write(*,*) "Needed if oblate=.true.: Planet mean radius (m)"
         Rmean = 0.0
         call getin_p("Rmean", Rmean)
         write(*,*) "Rmean = ", Rmean
         
! Test of incompatibility:
! if tlocked, then diurnal should be false
         if (tlocked.and.diurnal) then
           print*,'If diurnal=true, we should turn off tlocked.'
           stop
         endif

         write(*,*) "Tidal resonance ratio ?"
         nres=0          ! default value
         call getin_p("nres",nres)
         write(*,*) "nres = ",nres

         write(*,*) "Write some extra output to the screen ?"
         lwrite=.false. ! default value
         call getin_p("lwrite",lwrite)
         write(*,*) " lwrite = ",lwrite

         write(*,*) "Save statistics in file stats.nc ?"
         callstats=.true. ! default value
         call getin_p("callstats",callstats)
         write(*,*) " callstats = ",callstats

         write(*,*) "Test energy conservation of model physics ?"
         enertest=.false. ! default value
         call getin_p("enertest",enertest)
         write(*,*) " enertest = ",enertest

         write(*,*) "Check to see if cpp values used match gases.def ?"
         check_cpp_match=.true. ! default value
         call getin_p("check_cpp_match",check_cpp_match)
         write(*,*) " check_cpp_match = ",check_cpp_match

         write(*,*) "call radiative transfer ?"
         callrad=.true. ! default value
         call getin_p("callrad",callrad)
         write(*,*) " callrad = ",callrad

         write(*,*) "call correlated-k radiative transfer ?"
         corrk=.true. ! default value
         call getin_p("corrk",corrk)
         write(*,*) " corrk = ",corrk

         write(*,*) "prohibit calculations outside corrk T grid?"
         strictboundcorrk=.true. ! default value
         call getin_p("strictboundcorrk",strictboundcorrk)
         write(*,*) "strictboundcorrk = ",strictboundcorrk

         write(*,*) "call gaseous absorption in the visible bands?",
     &              "(matters only if callrad=T)"
         callgasvis=.false. ! default value
         call getin_p("callgasvis",callgasvis)
         write(*,*) " callgasvis = ",callgasvis
        
         write(*,*) "call continuum opacities in radiative transfer ?",
     &              "(matters only if callrad=T)"
         continuum=.true. ! default value
         call getin_p("continuum",continuum)
         write(*,*) " continuum = ",continuum

         write(*,*) "use analytic function for H2O continuum ?"
         H2Ocont_simple=.false. ! default value
         call getin_p("H2Ocont_simple",H2Ocont_simple)
         write(*,*) " H2Ocont_simple = ",H2Ocont_simple
 
         write(*,*) "call turbulent vertical diffusion ?"
         calldifv=.true. ! default value
         call getin_p("calldifv",calldifv)
         write(*,*) " calldifv = ",calldifv

         write(*,*) "use turbdiff instead of vdifc ?"
         UseTurbDiff=.true. ! default value
         call getin_p("UseTurbDiff",UseTurbDiff)
         write(*,*) " UseTurbDiff = ",UseTurbDiff

         write(*,*) "call convective adjustment ?"
         calladj=.true. ! default value
         call getin_p("calladj",calladj)
         write(*,*) " calladj = ",calladj

         write(*,*) "call CO2 condensation ?"
         co2cond=.false. ! default value
         call getin_p("co2cond",co2cond)
         write(*,*) " co2cond = ",co2cond
! Test of incompatibility
         if (co2cond.and.(.not.tracer)) then
            print*,'We need a CO2 ice tracer to condense CO2'
            call abort
         endif 
 
         write(*,*) "CO2 supersaturation level ?"
         co2supsat=1.0 ! default value
         call getin_p("co2supsat",co2supsat)
         write(*,*) " co2supsat = ",co2supsat

         write(*,*) "Radiative timescale for Newtonian cooling ?"
         tau_relax=30. ! default value
         call getin_p("tau_relax",tau_relax)
         write(*,*) " tau_relax = ",tau_relax
         tau_relax=tau_relax*24*3600 ! convert Earth days --> seconds

         write(*,*)"call thermal conduction in the soil ?"
         callsoil=.true. ! default value
         call getin_p("callsoil",callsoil)
         write(*,*) " callsoil = ",callsoil
         
         write(*,*)"Rad transfer is computed every iradia",
     &             " physical timestep"
         iradia=1 ! default value
         call getin_p("iradia",iradia)
         write(*,*)" iradia = ",iradia
       
         write(*,*)"Rayleigh scattering ?"
         rayleigh=.false.
         call getin_p("rayleigh",rayleigh)
         write(*,*)" rayleigh = ",rayleigh

         write(*,*) "Use blackbody for stellar spectrum ?"
         stelbbody=.false. ! default value
         call getin_p("stelbbody",stelbbody)
         write(*,*) " stelbbody = ",stelbbody

         write(*,*) "Stellar blackbody temperature ?"
         stelTbb=5800.0 ! default value
         call getin_p("stelTbb",stelTbb)
         write(*,*) " stelTbb = ",stelTbb

         write(*,*)"Output mean OLR in 1D?"
         meanOLR=.false.
         call getin_p("meanOLR",meanOLR)
         write(*,*)" meanOLR = ",meanOLR

         write(*,*)"Output spectral OLR in 3D?"
         specOLR=.false.
         call getin_p("specOLR",specOLR)
         write(*,*)" specOLR = ",specOLR

         write(*,*)"Operate in kastprof mode?"
         kastprof=.false.
         call getin_p("kastprof",kastprof)
         write(*,*)" kastprof = ",kastprof

         write(*,*)"Uniform absorption in radiative transfer?"
         graybody=.false.
         call getin_p("graybody",graybody)
         write(*,*)" graybody = ",graybody

! Slab Ocean 
         write(*,*) "Use slab-ocean ?"
         ok_slab_ocean=.false.         ! default value
         call getin_p("ok_slab_ocean",ok_slab_ocean)
         write(*,*) "ok_slab_ocean = ",ok_slab_ocean

         write(*,*) "Use slab-sea-ice ?"
         ok_slab_sic=.true.         ! default value
         call getin_p("ok_slab_sic",ok_slab_sic)
         write(*,*) "ok_slab_sic = ",ok_slab_sic

         write(*,*) "Use heat transport for the ocean ?"
         ok_slab_heat_transp=.true.   ! default value
         call getin_p("ok_slab_heat_transp",ok_slab_heat_transp)
         write(*,*) "ok_slab_heat_transp = ",ok_slab_heat_transp



! Test of incompatibility:
! if kastprof used, we must be in 1D
         if (kastprof.and.(ngrid.gt.1)) then
           print*,'kastprof can only be used in 1D!'
           call abort
         endif

         write(*,*)"Stratospheric temperature for kastprof mode?"
         Tstrat=167.0
         call getin_p("Tstrat",Tstrat)
         write(*,*)" Tstrat = ",Tstrat

         write(*,*)"Remove lower boundary?"
         nosurf=.false.
         call getin_p("nosurf",nosurf)
         write(*,*)" nosurf = ",nosurf

! Tests of incompatibility:
         if (nosurf.and.callsoil) then
           print*,'nosurf not compatible with soil scheme!'
           print*,'... got to make a choice!'
           call abort
         endif

         write(*,*)"Add an internal heat flux?",
     .             "... matters only if callsoil=F"
         intheat=0.
         call getin_p("intheat",intheat)
         write(*,*)" intheat = ",intheat

         write(*,*)"Use Newtonian cooling for radiative transfer?"
         newtonian=.false.
         call getin_p("newtonian",newtonian)
         write(*,*)" newtonian = ",newtonian

! Tests of incompatibility:
         if (newtonian.and.corrk) then
           print*,'newtonian not compatible with correlated-k!'
           call abort
         endif
         if (newtonian.and.calladj) then
           print*,'newtonian not compatible with adjustment!'
           call abort
         endif
         if (newtonian.and.calldifv) then
           print*,'newtonian not compatible with a boundary layer!'
           call abort
         endif

         write(*,*)"Test physics timescale in 1D?"
         testradtimes=.false.
         call getin_p("testradtimes",testradtimes)
         write(*,*)" testradtimes = ",testradtimes

! Test of incompatibility:
! if testradtimes used, we must be in 1D
         if (testradtimes.and.(ngrid.gt.1)) then
           print*,'testradtimes can only be used in 1D!'
           call abort
         endif

         write(*,*)"Default planetary temperature?"
         tplanet=215.0
         call getin_p("tplanet",tplanet)
         write(*,*)" tplanet = ",tplanet

         write(*,*)"Which star?"
         startype=1 ! default value = Sol
         call getin_p("startype",startype)
         write(*,*)" startype = ",startype

         write(*,*)"Value of stellar flux at 1 AU?"
         Fat1AU=1356.0 ! default value = Sol today
         call getin_p("Fat1AU",Fat1AU)
         write(*,*)" Fat1AU = ",Fat1AU


! TRACERS:

         write(*,*)"Varying H2O cloud fraction?"
         CLFvarying=.false.     ! default value
         call getin_p("CLFvarying",CLFvarying)
         write(*,*)" CLFvarying = ",CLFvarying

         write(*,*)"Value of fixed H2O cloud fraction?"
         CLFfixval=1.0                ! default value
         call getin_p("CLFfixval",CLFfixval)
         write(*,*)" CLFfixval = ",CLFfixval

         write(*,*)"fixed radii for Cloud particles?"
         radfixed=.false. ! default value
         call getin_p("radfixed",radfixed)
         write(*,*)" radfixed = ",radfixed

         if(kastprof)then
            radfixed=.true.
         endif  

	 write(*,*)"Number mixing ratio of CO2 ice particles:"
         Nmix_co2=1.e6 ! default value
         call getin_p("Nmix_co2",Nmix_co2)
         write(*,*)" Nmix_co2 = ",Nmix_co2

!         write(*,*)"Number of radiatively active aerosols:"
!         naerkind=0. ! default value
!         call getin_p("naerkind",naerkind)
!         write(*,*)" naerkind = ",naerkind

         write(*,*)"Opacity of dust (if used):"
         dusttau=0. ! default value
         call getin_p("dusttau",dusttau)
         write(*,*)" dusttau = ",dusttau

         write(*,*)"Radiatively active CO2 aerosols?"
         aeroco2=.false.     ! default value
         call getin_p("aeroco2",aeroco2)
         write(*,*)" aeroco2 = ",aeroco2

         write(*,*)"Fixed CO2 aerosol distribution?"
         aerofixco2=.false.     ! default value
         call getin_p("aerofixco2",aerofixco2)
         write(*,*)" aerofixco2 = ",aerofixco2

         write(*,*)"Radiatively active water ice?"
         aeroh2o=.false.     ! default value
         call getin_p("aeroh2o",aeroh2o)
         write(*,*)" aeroh2o = ",aeroh2o

         write(*,*)"Fixed H2O aerosol distribution?"
         aerofixh2o=.false.     ! default value
         call getin_p("aerofixh2o",aerofixh2o)
         write(*,*)" aerofixh2o = ",aerofixh2o

         write(*,*)"Radiatively active sulfuric acid aersols?"
         aeroh2so4=.false.     ! default value
         call getin_p("aeroh2so4",aeroh2so4)
         write(*,*)" aeroh2so4 = ",aeroh2so4
	 
!=================================

	 write(*,*)"Radiatively active two-layer aersols?"
         aeroback2lay=.false.     ! default value
         call getin_p("aeroback2lay",aeroback2lay)
         write(*,*)" aeroback2lay = ",aeroback2lay

         write(*,*)"TWOLAY AEROSOL: total optical depth ",
     &              "in the tropospheric layer (visible)"
         obs_tau_col_tropo=8.D0
         call getin_p("obs_tau_col_tropo",obs_tau_col_tropo)
         write(*,*)" obs_tau_col_tropo = ",obs_tau_col_tropo

         write(*,*)"TWOLAY AEROSOL: total optical depth ",
     &              "in the stratospheric layer (visible)"
         obs_tau_col_strato=0.08D0
         call getin_p("obs_tau_col_strato",obs_tau_col_strato)
         write(*,*)" obs_tau_col_strato = ",obs_tau_col_strato

         write(*,*)"TWOLAY AEROSOL: pres_bottom_tropo? in pa"
         pres_bottom_tropo=66000.0
         call getin_p("pres_bottom_tropo",pres_bottom_tropo)
         write(*,*)" pres_bottom_tropo = ",pres_bottom_tropo

         write(*,*)"TWOLAY AEROSOL: pres_top_tropo? in pa"
         pres_top_tropo=18000.0
         call getin_p("pres_top_tropo",pres_top_tropo)
         write(*,*)" pres_top_tropo = ",pres_top_tropo

         write(*,*)"TWOLAY AEROSOL: pres_bottom_strato? in pa"
         pres_bottom_strato=2000.0
         call getin_p("pres_bottom_strato",pres_bottom_strato)
         write(*,*)" pres_bottom_strato = ",pres_bottom_strato

         write(*,*)"TWOLAY AEROSOL: pres_top_strato? in pa"
         pres_top_strato=100.0
         call getin_p("pres_top_strato",pres_top_strato)
         write(*,*)" pres_top_strato = ",pres_top_strato

         write(*,*)"TWOLAY AEROSOL: particle size in the ",
     &              "tropospheric layer, in meters"
         size_tropo=2.e-6
         call getin_p("size_tropo",size_tropo)
         write(*,*)" size_tropo = ",size_tropo

         write(*,*)"TWOLAY AEROSOL: particle size in the ",
     &              "stratospheric layer, in meters"
         size_strato=1.e-7
         call getin_p("size_strato",size_strato)
         write(*,*)" size_strato = ",size_strato

!=================================

         write(*,*)"Cloud pressure level (with kastprof only):"
         cloudlvl=0. ! default value
         call getin_p("cloudlvl",cloudlvl)
         write(*,*)" cloudlvl = ",cloudlvl

         write(*,*)"Is the variable gas species radiatively active?"
         Tstrat=167.0
         varactive=.false.
         call getin_p("varactive",varactive)
         write(*,*)" varactive = ",varactive

         write(*,*)"Is the variable gas species distribution set?"
         varfixed=.false.
         call getin_p("varfixed",varfixed)
         write(*,*)" varfixed = ",varfixed

         write(*,*)"What is the saturation % of the variable species?"
         satval=0.8
         call getin_p("satval",satval)
         write(*,*)" satval = ",satval


! Test of incompatibility:
! if varactive, then varfixed should be false
         if (varactive.and.varfixed) then
           print*,'if varactive, varfixed must be OFF!'
           stop
         endif

         write(*,*) "Gravitationnal sedimentation ?"
         sedimentation=.false. ! default value
         call getin_p("sedimentation",sedimentation)
         write(*,*) " sedimentation = ",sedimentation

         write(*,*) "Compute water cycle ?"
         water=.false. ! default value
         call getin_p("water",water)
         write(*,*) " water = ",water
         
! Test of incompatibility:
! if water is true, there should be at least a tracer
         if (water.and.(.not.tracer)) then
           print*,'if water is ON, tracer must be ON too!'
           stop
         endif

         write(*,*) "Include water condensation ?"
         watercond=.false. ! default value
         call getin_p("watercond",watercond)
         write(*,*) " watercond = ",watercond

! Test of incompatibility:
! if watercond is used, then water should be used too
         if (watercond.and.(.not.water)) then
           print*,'if watercond is used, water should be used too'
           stop
         endif

         write(*,*) "Include water precipitation ?"
         waterrain=.false. ! default value
         call getin_p("waterrain",waterrain)
         write(*,*) " waterrain = ",waterrain

         write(*,*) "Include surface hydrology ?"
         hydrology=.false. ! default value
         call getin_p("hydrology",hydrology)
         write(*,*) " hydrology = ",hydrology

         write(*,*) "Evolve surface water sources ?"
         sourceevol=.false. ! default value
         call getin_p("sourceevol",sourceevol)
         write(*,*) " sourceevol = ",sourceevol

         write(*,*) "Ice evolution timestep ?"
         icetstep=100.0 ! default value
         call getin_p("icetstep",icetstep)
         write(*,*) " icetstep = ",icetstep

         write(*,*) "Snow albedo ?"
         albedosnow=0.5         ! default value
         call getin_p("albedosnow",albedosnow)
         write(*,*) " albedosnow = ",albedosnow

         write(*,*) "Maximum ice thickness ?"
         maxicethick=2.0         ! default value
         call getin_p("maxicethick",maxicethick)
         write(*,*) " maxicethick = ",maxicethick

         write(*,*) "Freezing point of seawater ?"
         Tsaldiff=-1.8          ! default value
         call getin_p("Tsaldiff",Tsaldiff)
         write(*,*) " Tsaldiff = ",Tsaldiff

         write(*,*) "Does user want to force cpp and mugaz?"
         force_cpp=.false. ! default value
         call getin_p("force_cpp",force_cpp)
         write(*,*) " force_cpp = ",force_cpp

         if (force_cpp) then
           mugaz = -99999.
           PRINT *,'MEAN MOLECULAR MASS in g mol-1 ?'
           call getin_p("mugaz",mugaz)
           IF (mugaz.eq.-99999.) THEN
               PRINT *, "mugaz must be set if force_cpp = T"
               STOP
           ELSE
               write(*,*) "mugaz=",mugaz
           ENDIF
           cpp = -99999.
           PRINT *,'SPECIFIC HEAT CAPACITY in J K-1 kg-1 ?'
           call getin_p("cpp",cpp)
           IF (cpp.eq.-99999.) THEN
               PRINT *, "cpp must be set if force_cpp = T"
               STOP
           ELSE
               write(*,*) "cpp=",cpp
           ENDIF
!         else
!           mugaz=8.314*1000./pr
         endif
         call su_gases
         call calc_cpp_mugaz

         PRINT*,'--------------------------------------------'
         PRINT*
         PRINT*
      ELSE
         write(*,*)
         write(*,*) 'Cannot read file callphys.def. Is it here ?'
         stop
      ENDIF

8000  FORMAT(t5,a12,l8)
8001  FORMAT(t5,a12,i8)

      PRINT*
      PRINT*,'inifis: daysec',daysec
      PRINT*
      PRINT*,'inifis: The radiative transfer is computed:'
      PRINT*,'           each ',iradia,' physical time-step'
      PRINT*,'        or each ',iradia*dtphys,' seconds'
      PRINT*


!-----------------------------------------------------------------------
!     Some more initialization:
!     ------------------------

      ! ALLOCATE ARRAYS IN comgeomfi_h
      IF (.not. ALLOCATED(lati)) ALLOCATE(lati(ngrid))
      IF (.not. ALLOCATED(long)) ALLOCATE(long(ngrid))
      IF (.not. ALLOCATED(area)) ALLOCATE(area(ngrid))

      CALL SCOPY(ngrid,plon,1,long,1)
      CALL SCOPY(ngrid,plat,1,lati,1)
      CALL SCOPY(ngrid,parea,1,area,1)
      totarea=SSUM(ngrid,area,1)
      call planetwide_sumval(area,totarea_planet)

      !! those are defined in comdiurn_h.F90
      IF (.not.ALLOCATED(sinlat)) ALLOCATE(sinlat(ngrid))
      IF (.not.ALLOCATED(coslat)) ALLOCATE(coslat(ngrid))
      IF (.not.ALLOCATED(sinlon)) ALLOCATE(sinlon(ngrid))
      IF (.not.ALLOCATED(coslon)) ALLOCATE(coslon(ngrid))

      DO ig=1,ngrid
         sinlat(ig)=sin(plat(ig))
         coslat(ig)=cos(plat(ig))
         sinlon(ig)=sin(plon(ig))
         coslon(ig)=cos(plon(ig))
      ENDDO

!$OMP MASTER
      pi=2.*asin(1.) ! NB: pi is a common in comcstfi.h
!$OMP END MASTER
!$OMP BARRIER

      ! allocate "comsoil_h" arrays
      call ini_comsoil_h(ngrid)
      
      END
