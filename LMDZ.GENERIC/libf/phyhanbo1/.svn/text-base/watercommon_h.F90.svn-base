module watercommon_h

      implicit none

      real, parameter :: T_coup = 234.0
      real, parameter :: T_h2O_ice_liq = 273.16
      real, parameter :: T_h2O_ice_clouds = T_h2O_ice_liq-15. 
      real, parameter :: mH2O = 18.01528   

      ! benjamin additions
      real, parameter :: RLVTT = 2.257E+6 ! Latent heat of vaporization (J kg-1) 
      real, parameter :: RLSTT = 2.257E+6 ! 2.591E+6 in reality ! Latent heat of sublimation (J kg-1)

      real, parameter :: RLFTT = 3.334E+5 ! Latent heat of fusion (J kg-1) ! entails an energy sink but better description of albedo
      real, parameter :: rhowater = 1.0E+3 ! mass of water (kg/m^3)
      real, parameter :: rhowaterice = 9.2E+2 ! mass of water (kg/m^3)
      real, parameter :: capcal_h2o_liq = 4181.3 ! specific heat capacity of liquid water J/kg/K
      real, parameter :: mx_eau_sol = 150 ! mass of water (kg/m^2)

      real, save :: epsi, RCPD, RCPV, RV, RVTMP2
!$OMP THREADPRIVATE(epsi,RCPD,RCPV,RV,RVTMP2)
      
      contains

      
!==================================================================
      subroutine Psat_water(T,p,psat,qsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the saturation vapor pressure and mass mixing ratio at saturation (kg/kg)
!     for a given pressure (Pa) and temperature (K)
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         real, intent(in) :: T, p
  
!        output
         real psat,qsat

! JL12 variables for tetens formula
         real,parameter :: Pref_solid_liquid=611.14
         real,parameter :: Trefvaporization=35.86
         real,parameter :: Trefsublimation=7.66
         real,parameter :: Tmin=8.
         real,parameter :: r3vaporization=17.269
         real,parameter :: r3sublimation=21.875

! checked vs. old watersat data 14/05/2012 by JL.

         if (T.gt.T_h2O_ice_liq) then 
            psat = Pref_solid_liquid*Exp(r3vaporization*(T-T_h2O_ice_liq)/(T-Trefvaporization)) ! liquid / vapour
         else if (T.lt.Tmin) then
	    print*, "careful, T<Tmin in psat water"
          !  psat = Pref_solid_liquid*Exp(r3sublimation*(Tmin-T_h2O_ice_liq)/(Tmin-Trefsublimation)) ! min psat  
         ! Ehouarn: gfortran says: Error: Result of EXP underflows its kind,
         !          so set psat to the smallest possible value instead
            psat=tiny(psat)
         else                 
            psat = Pref_solid_liquid*Exp(r3sublimation*(T-T_h2O_ice_liq)/(T-Trefsublimation)) ! solid / vapour
         endif
         if(psat.gt.p) then
            qsat=1.
         else
            qsat=epsi*psat/(p-(1.-epsi)*psat)
         endif
         return
      end subroutine Psat_water




!==================================================================
      subroutine Lcpdqsat_water(T,p,psat,qsat,dqsat,dlnpsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute dqsat=L/cp*d (q_sat)/d T and dlnpsat=L/cp d(ln Psat)/d T
!     for a given temperature (K)! 
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         real T, p, psat, qsat
  
!        output
         real dqsat,dlnpsat

! JL12 variables for tetens formula
         real,parameter :: Pref_solid_liquid=611.14
         real,parameter :: Trefvaporization=35.86
         real,parameter :: Tmin=8.
         real,parameter :: Trefsublimation=7.66
         real,parameter :: r3vaporization=17.269
         real,parameter :: r3sublimation=21.875

         real :: dummy

         if (psat.gt.p) then
	    dqsat=0.
	    return
	 endif

         if (T.gt.T_h2O_ice_liq) then
            dummy = r3vaporization*(T_h2O_ice_liq-Trefvaporization)/(T-Trefvaporization)**2  ! liquid / vapour
         else if (T.lt.Tmin) then
	    print*, "careful, T<Tmin in Lcp psat water"
            dummy = r3sublimation*(T_h2O_ice_liq-Trefsublimation)/(Tmin-Trefsublimation)**2  ! solid / vapour
         else               
            dummy = r3sublimation*(T_h2O_ice_liq-Trefsublimation)/(T-Trefsublimation)**2  ! solid / vapour
         endif

         dqsat=RLVTT/RCPD*qsat*(p/(p-(1.-epsi)*psat))*dummy
	 dlnpsat=RLVTT/RCPD*dummy
         return
      end subroutine Lcpdqsat_water




!==================================================================
      subroutine Tsat_water(p,Tsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the saturation temperature
!     for a given pressure (Pa) 
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         real p
  
!        output
         real Tsat

! JL12 variables for tetens formula
         real,parameter :: Pref_solid_liquid=611.14
         real,parameter :: Trefvaporization=35.86
         real,parameter :: Trefsublimation=7.66
         real,parameter :: r3vaporization=17.269
         real,parameter :: r3sublimation=21.875

         if (p.lt.Pref_solid_liquid) then ! solid / vapour
            Tsat =(T_h2O_ice_liq*r3sublimation- Trefsublimation*Log(p/Pref_solid_liquid))/(r3sublimation-Log(p/Pref_solid_liquid))
         else                 ! liquid / vapour
            Tsat =(T_h2O_ice_liq*r3vaporization- Trefvaporization*Log(p/Pref_solid_liquid))/(r3vaporization-Log(p/Pref_solid_liquid))
         endif

         return
      end subroutine Tsat_water

!==================================================================
      subroutine Psat_waterDP(T,p,psat,qsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the saturation vapor pressure and mass mixing ratio at saturation (kg/kg)
!     for a given pressure (Pa) and temperature (K)
!     DOUBLE PRECISION
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         double precision, intent(in) :: T, p
  
!        output
         double precision psat,qsat

! JL12 variables for tetens formula
         double precision,parameter :: Pref_solid_liquid=611.14d0
         double precision,parameter :: Trefvaporization=35.86D0
         double precision,parameter :: Trefsublimation=7.66d0
         double precision,parameter :: Tmin=8.d0
         double precision,parameter :: r3vaporization=17.269d0
         double precision,parameter :: r3sublimation=21.875d0

! checked vs. old watersat data 14/05/2012 by JL.

         if (T.gt.T_h2O_ice_liq) then 
            psat = Pref_solid_liquid*Exp(r3vaporization*(T-T_h2O_ice_liq)/(T-Trefvaporization)) ! liquid / vapour
         else if (T.lt.Tmin) then
	    print*, "careful, T<Tmin in psat water"
         !   psat = Pref_solid_liquid*Exp(r3sublimation*(Tmin-T_h2O_ice_liq)/(Tmin-Trefsublimation)) ! min psat  
         ! Ehouarn: gfortran says: Error: Result of EXP underflows its kind,
         !          so set psat to the smallest possible value instead
            psat=tiny(psat)
         else                 
            psat = Pref_solid_liquid*Exp(r3sublimation*(T-T_h2O_ice_liq)/(T-Trefsublimation)) ! solid / vapour
         endif
         if(psat.gt.p) then
            qsat=1.d0
         else
            qsat=epsi*psat/(p-(1.d0-epsi)*psat)
         endif
         return
      end subroutine Psat_waterDP




!==================================================================
      subroutine Lcpdqsat_waterDP(T,p,psat,qsat,dqsat,dlnpsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute dqsat=L/cp*d (q_sat)/d T and dlnpsat=L/cp d(ln Psat)/d T
!     for a given temperature (K)! 
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         double precision T, p, psat, qsat
  
!        output
         double precision dqsat,dlnpsat

! JL12 variables for tetens formula
         double precision,parameter :: Pref_solid_liquid=611.14d0
         double precision,parameter :: Trefvaporization=35.86d0
         double precision,parameter :: Tmin=8.d0
         double precision,parameter :: Trefsublimation=7.66d0
         double precision,parameter :: r3vaporization=17.269d0
         double precision,parameter :: r3sublimation=21.875d0

         double precision :: dummy

         if (psat.gt.p) then
	    dqsat=0.d0
	    return
	 endif

         if (T.gt.T_h2O_ice_liq) then
            dummy = r3vaporization*(T_h2O_ice_liq-Trefvaporization)/(T-Trefvaporization)**2  ! liquid / vapour
         else if (T.lt.Tmin) then
	    print*, "careful, T<Tmin in Lcp psat water"
            dummy = r3sublimation*(T_h2O_ice_liq-Trefsublimation)/(Tmin-Trefsublimation)**2  ! solid / vapour
         else               
            dummy = r3sublimation*(T_h2O_ice_liq-Trefsublimation)/(T-Trefsublimation)**2  ! solid / vapour
         endif

         dqsat=RLVTT/RCPD*qsat*(p/(p-(1.d0-epsi)*psat))*dummy
	 dlnpsat=RLVTT/RCPD*dummy
         return
      end subroutine Lcpdqsat_waterDP


end module watercommon_h
