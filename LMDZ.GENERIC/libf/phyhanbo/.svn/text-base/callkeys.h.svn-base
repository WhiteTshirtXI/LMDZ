!
! For Fortran 77/Fortran 90 compliance always use line continuation
! symbols '&' in columns 73 and 6
!
! Group commons according to their type for minimal performance impact

      COMMON/callkeys_l/callrad,corrk,calldifv,UseTurbDiff,calladj      &
     &   , co2cond,callsoil                                             &
     &   , season,diurnal,tlocked,rings_shadow,lwrite                   &
     &   , callstats,calleofdump                                        &
     &   , enertest                                                     &
     &   , callgasvis,continuum,H2Ocont_simple,graybody                 &
     &   , radfixed                                                     &
     &   , meanOLR, specOLR                                             &
     &   , kastprof                                                     &
     &   , nosurf, oblate                                               &     
     &   , newtonian, testradtimes                                      &
     &   , check_cpp_match, force_cpp                                   &
     &   , rayleigh                                                     &
     &   , stelbbody                                                    &
     &   , nearco2cond                                                  &
     &   , tracer, mass_redistrib, varactive, varfixed                  &
     &   , sedimentation,water,watercond,waterrain                      &
     &   , aeroco2,aeroh2o,aeroh2so4,aeroback2lay                       &
     &   , aerofixco2,aerofixh2o                                        &
     &   , hydrology, sourceevol                                        &
     &   , CLFvarying                                                   &
     &   , strictboundcorrk                                             &                                        
     &   , ok_slab_ocean                                                &
     &   , ok_slab_sic                                                  &
     &   , ok_slab_heat_transp                                          


      COMMON/callkeys_i/iaervar,iddist,iradia,startype
      
      COMMON/callkeys_r/topdustref,Nmix_co2,dusttau,Fat1AU,stelTbb,     &
     &                  Tstrat,tplanet,obs_tau_col_tropo,               &
     &                  obs_tau_col_strato,pres_bottom_tropo,           &
     &                  pres_top_tropo,pres_bottom_strato,              &
     &                  pres_top_strato,size_tropo,size_strato,satval,  &
     &                  CLFfixval,n2mixratio,co2supsat,pceil,albedosnow,&
     &                  maxicethick,Tsaldiff,tau_relax,cloudlvl,        &
     &                  icetstep,intheat,flatten,Rmean,J2,MassPlanet
      
      logical callrad,corrk,calldifv,UseTurbDiff                        &
     &   , calladj,co2cond,callsoil                                     &
     &   , season,diurnal,tlocked,rings_shadow,lwrite                   &
     &   , callstats,calleofdump                                        &
     &   , callgasvis,continuum,H2Ocont_simple,graybody                 &
     &   , strictboundcorrk                                             

      logical enertest
      logical nonideal
      logical meanOLR
      logical specOLR
      logical kastprof
      logical newtonian
      logical check_cpp_match
      logical force_cpp
      logical testradtimes
      logical rayleigh
      logical stelbbody
      logical ozone
      logical nearco2cond
      logical tracer
      logical mass_redistrib
      logical varactive
      logical varfixed
      logical radfixed
      logical sedimentation
      logical water,watercond,waterrain
      logical aeroco2,aeroh2o,aeroh2so4,aeroback2lay
      logical aerofixco2,aerofixh2o
      logical hydrology
      logical sourceevol
      logical CLFvarying
      logical nosurf
      logical oblate
      logical ok_slab_ocean
      logical ok_slab_sic
      logical ok_slab_heat_transp

      integer iddist
      integer iaervar
      integer iradia
      integer startype

      real topdustref
      real Nmix_co2
      real dusttau
      real Fat1AU
      real stelTbb
      real Tstrat
      real tplanet
      real obs_tau_col_tropo
      real obs_tau_col_strato
      real pres_bottom_tropo
      real pres_top_tropo
      real pres_bottom_strato
      real pres_top_strato
      real size_tropo
      real size_strato
      real satval
      real CLFfixval
      real n2mixratio
      real co2supsat
      real pceil
      real albedosnow
      real maxicethick
      real Tsaldiff
      real tau_relax
      real cloudlvl
      real icetstep
      real intheat
      real flatten
      real Rmean
      real J2
      real MassPlanet
      
      logical :: iscallphys=.false.!existence of callphys.def
