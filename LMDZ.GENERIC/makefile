#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Definitions de Macros pour Make
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# create_make_gcm arguments: dyn3d grid bibio phyhanbo

# Repertoires :

GCM     = /home/hanbo/GCM_program/LMDZ.GENERIC
MACHINE = sun4
LIBF    = $(GCM)/libf
LIBO    = $(LIBOGCM)/$(MACHINE)
LOCAL_DIR=/cita/d/raid-cita/hanbo/GCM/Thermaltides/10b_n2gray02
BIBIO    = $(LIBF)/bibio
FILTRE   = filtre
PHYS  = 
DYN  = $(DYNTYPE)
LIBPHY = $(LIBO)/libphy$(PHYS).a
DIRMAIN= $(DYN)$(DIM)d
RM=rm

OPLINK = 

# Les differentes librairies pour l"edition des liens:

dyn3d      = $(LIBO)/libdyn3d.a $(LIBO)/lib$(FILTRE).a
dyn2d      = $(LIBO)/libdyn2d.a
dyn1d      = $(LIBO)/libdyn1d.a
L_DYN      = -l$(DYN)$(DIM)d
L_FILTRE   = -l$(FILTRE)
L_PHY = -lphy$(PHYS) 
L_BIBIO    = -lbibio
L_ADJNT    =

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Option de compilation FORTRAN
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMPILE = $(F77) $(OPTIM) $(INCLUDE) -c
COMPILE90 = $(F90) $(OPTIM90) $(INCLUDE) -c
COMPTRU90 = $(F90) $(OPTIMTRU90) $(INCLUDE) -c
LINK    = 
AR      = 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creation des differents executables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Executables:
# ------------

PROG = code

main : $(DYN) bibio phys $(OPTION_DEP) 
	cd $(LIBO) ; $(RANLIB) lib*.a ; cd $(GCM) ;\
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/$(DIRMAIN)/$(SOURCE) -o $(PROG).o ; \
	$(LINK) $(PROG).o -L$(LIBO) $(L_DYN) $(L_ADJNT) $(L_FILTRE) $(L_PHY) $(L_BIBIO) $(L_DYN) $(OPLINK) $(OPTION_LINK) -o $(LOCAL_DIR)/$(PROG).e ; $(RM) $(PROG).o 

dyn : $(LIBO)/libdyn$(DIM)d.a $(FILTRE)$(DIM)d

olddyn : $(LIBO)/libolddyn$(DIM)d.a $(FILTRE)$(DIM)d

phys : $(LIBPHY)

bibio : $(LIBO)/libbibio.a

adjnt : $(LIBO)/libadjnt.a

$(FILTRE)3d : $(LIBO)/lib$(FILTRE).a

$(FILTRE)2d :

$(FILTRE)1d :

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Contenu des differentes bibliotheques
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory bibio
#=======================================================================

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cbrt.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cvmgp.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cvmgt.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ilsum.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(isamax.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(isitor.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ismax.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ismin.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxm.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxva.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxv.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scopy.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scopyi.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(sdot.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ssum.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(wheneq.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ioipsl_errioipsl.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ioipsl_getincom.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ioipsl_getincom_p.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ioipsl_stringop.o)

.PRECIOUS	: $(LIBO)/libbibio.a


# Compilation des membres de la bibliotheque libbibio.a

$(LIBO)/libbibio.a(cbrt.o) : $(LIBF)/bibio/cbrt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a cbrt.o ; \
	$(COMPILE) $(LIBF)/bibio/cbrt.F ; \
	$(AR) r $(LIBO)/libbibio.a cbrt.o ; $(RM) cbrt.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(cvmgp.o) : $(LIBF)/bibio/cvmgp.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a cvmgp.o ; \
	$(COMPILE) $(LIBF)/bibio/cvmgp.F ; \
	$(AR) r $(LIBO)/libbibio.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(cvmgt.o) : $(LIBF)/bibio/cvmgt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a cvmgt.o ; \
	$(COMPILE) $(LIBF)/bibio/cvmgt.F ; \
	$(AR) r $(LIBO)/libbibio.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ilsum.o) : $(LIBF)/bibio/ilsum.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ilsum.o ; \
	$(COMPILE) $(LIBF)/bibio/ilsum.F ; \
	$(AR) r $(LIBO)/libbibio.a ilsum.o ; $(RM) ilsum.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(isamax.o) : $(LIBF)/bibio/isamax.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a isamax.o ; \
	$(COMPILE) $(LIBF)/bibio/isamax.F ; \
	$(AR) r $(LIBO)/libbibio.a isamax.o ; $(RM) isamax.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(isitor.o) : $(LIBF)/bibio/isitor.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a isitor.o ; \
	$(COMPILE) $(LIBF)/bibio/isitor.F ; \
	$(AR) r $(LIBO)/libbibio.a isitor.o ; $(RM) isitor.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ismax.o) : $(LIBF)/bibio/ismax.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ismax.o ; \
	$(COMPILE) $(LIBF)/bibio/ismax.F ; \
	$(AR) r $(LIBO)/libbibio.a ismax.o ; $(RM) ismax.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ismin.o) : $(LIBF)/bibio/ismin.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ismin.o ; \
	$(COMPILE) $(LIBF)/bibio/ismin.F ; \
	$(AR) r $(LIBO)/libbibio.a ismin.o ; $(RM) ismin.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxm.o) : $(LIBF)/bibio/mxm.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a mxm.o ; \
	$(COMPILE) $(LIBF)/bibio/mxm.F ; \
	$(AR) r $(LIBO)/libbibio.a mxm.o ; $(RM) mxm.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxva.o) : $(LIBF)/bibio/mxva.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a mxva.o ; \
	$(COMPILE) $(LIBF)/bibio/mxva.F ; \
	$(AR) r $(LIBO)/libbibio.a mxva.o ; $(RM) mxva.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxv.o) : $(LIBF)/bibio/mxv.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a mxv.o ; \
	$(COMPILE) $(LIBF)/bibio/mxv.F ; \
	$(AR) r $(LIBO)/libbibio.a mxv.o ; $(RM) mxv.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scopy.o) : $(LIBF)/bibio/scopy.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a scopy.o ; \
	$(COMPILE) $(LIBF)/bibio/scopy.F ; \
	$(AR) r $(LIBO)/libbibio.a scopy.o ; $(RM) scopy.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scopyi.o) : $(LIBF)/bibio/scopyi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a scopyi.o ; \
	$(COMPILE) $(LIBF)/bibio/scopyi.F ; \
	$(AR) r $(LIBO)/libbibio.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(sdot.o) : $(LIBF)/bibio/sdot.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a sdot.o ; \
	$(COMPILE) $(LIBF)/bibio/sdot.F ; \
	$(AR) r $(LIBO)/libbibio.a sdot.o ; $(RM) sdot.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ssum.o) : $(LIBF)/bibio/ssum.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ssum.o ; \
	$(COMPILE) $(LIBF)/bibio/ssum.F ; \
	$(AR) r $(LIBO)/libbibio.a ssum.o ; $(RM) ssum.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(wheneq.o) : $(LIBF)/bibio/wheneq.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a wheneq.o ; \
	$(COMPILE) $(LIBF)/bibio/wheneq.F ; \
	$(AR) r $(LIBO)/libbibio.a wheneq.o ; $(RM) wheneq.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ioipsl_errioipsl.o) : $(LIBF)/bibio/ioipsl_errioipsl.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ioipsl_errioipsl.o ; \
	$(COMPTRU90) $(LIBF)/bibio/ioipsl_errioipsl.F90 ; \
	mv $(MOD_LOC_DIR)/ioipsl_errioipsl.$(MOD_SUFFIX) $(LIBO)/ioipsl_errioipsl.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libbibio.a ioipsl_errioipsl.o ; $(RM) ioipsl_errioipsl.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ioipsl_getincom.o) : $(LIBF)/bibio/ioipsl_getincom.F90 \
$(LIBO)/libbibio.a(ioipsl_errioipsl.o) \
$(LIBO)/libbibio.a(ioipsl_stringop.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ioipsl_getincom.o ; \
	$(COMPTRU90) $(LIBF)/bibio/ioipsl_getincom.F90 ; \
	mv $(MOD_LOC_DIR)/ioipsl_getincom.$(MOD_SUFFIX) $(LIBO)/ioipsl_getincom.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libbibio.a ioipsl_getincom.o ; $(RM) ioipsl_getincom.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ioipsl_getincom_p.o) : $(LIBF)/bibio/ioipsl_getincom_p.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ioipsl_getincom_p.o ; \
	$(COMPTRU90) $(LIBF)/bibio/ioipsl_getincom_p.F90 ; \
	mv $(MOD_LOC_DIR)/ioipsl_getincom_p.$(MOD_SUFFIX) $(LIBO)/ioipsl_getincom_p.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libbibio.a ioipsl_getincom_p.o ; $(RM) ioipsl_getincom_p.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ioipsl_stringop.o) : $(LIBF)/bibio/ioipsl_stringop.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ioipsl_stringop.o ; \
	$(COMPTRU90) $(LIBF)/bibio/ioipsl_stringop.F90 ; \
	mv $(MOD_LOC_DIR)/ioipsl_stringop.$(MOD_SUFFIX) $(LIBO)/ioipsl_stringop.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libbibio.a ioipsl_stringop.o ; $(RM) ioipsl_stringop.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory dyn3d
#=======================================================================

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(abort_gcm.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(addfi.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(advect.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(adv_h2o.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(amont.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(amont_qsat.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(bernoui.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(caldyn0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(caldyn.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(calfis.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(coefpoly.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(convflu.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(convmas.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(coordij.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(covcont.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(covnat.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(defrun_new.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(diagadv.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dissip.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(disvert.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(diverg.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divergf.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(diverg_gam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divergst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divgrad22.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divgrad2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divgrad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dteta1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dudv1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dudv2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dump2d.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dynetat0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dynredem.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(enercin.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(exner.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(exner_hyb.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(flumass.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fluxstoke.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxhyp.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxy.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxyhyper.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxysinus.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fyhyp.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(geopot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gradiv2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gradiv.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_dyn_fi.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_fi_dyn.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grid_atob.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grid_noro1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grid_noro.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_int_dyn.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(groupe.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(groupeun.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_u_scal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_v_scal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(heavyside.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniav.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniconst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(inidissip.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(inigeom.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniinterp_h.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(initial0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(integrd.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(interp_horiz.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(interp_vert.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien_gam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien_rot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien_rotgam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(lectflux.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(massbar.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(massbarxy.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(massdair.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(minmax.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(multipl.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(multscal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgrad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgrad_gam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgradst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgraro2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgrarot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(pression.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(psextbar.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(ran1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotat.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotatf.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotat_nfil.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotatst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(scal_wind.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(sortvarc0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(sortvarc.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(spline.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(splint.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(test_period.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(tourpot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(traceur.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(tracvl.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vanleer.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vitvert.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vlsplt.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vlspltqs.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(wind_scal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(control_mod.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(infotrac.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(mod_const_mpi.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(sponge_mod.o)

.PRECIOUS	: $(LIBO)/libdyn3d.a


# Compilation des membres de la bibliotheque libdyn3d.a

$(LIBO)/libdyn3d.a(abort_gcm.o) : $(LIBF)/dyn3d/abort_gcm.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a abort_gcm.o ; \
	$(COMPILE) $(LIBF)/dyn3d/abort_gcm.F ; \
	$(AR) r $(LIBO)/libdyn3d.a abort_gcm.o ; $(RM) abort_gcm.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(addfi.o) : $(LIBF)/dyn3d/addfi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a addfi.o ; \
	$(COMPILE) $(LIBF)/dyn3d/addfi.F ; \
	$(AR) r $(LIBO)/libdyn3d.a addfi.o ; $(RM) addfi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(advect.o) : $(LIBF)/dyn3d/advect.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/ener.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a advect.o ; \
	$(COMPILE) $(LIBF)/dyn3d/advect.F ; \
	$(AR) r $(LIBO)/libdyn3d.a advect.o ; $(RM) advect.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(adv_h2o.o) : $(LIBF)/dyn3d/adv_h2o.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a adv_h2o.o ; \
	$(COMPILE) $(LIBF)/dyn3d/adv_h2o.F ; \
	$(AR) r $(LIBO)/libdyn3d.a adv_h2o.o ; $(RM) adv_h2o.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(amont.o) : $(LIBF)/dyn3d/amont.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a amont.o ; \
	$(COMPILE) $(LIBF)/dyn3d/amont.F ; \
	$(AR) r $(LIBO)/libdyn3d.a amont.o ; $(RM) amont.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(amont_qsat.o) : $(LIBF)/dyn3d/amont_qsat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a amont_qsat.o ; \
	$(COMPILE) $(LIBF)/dyn3d/amont_qsat.F ; \
	$(AR) r $(LIBO)/libdyn3d.a amont_qsat.o ; $(RM) amont_qsat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(bernoui.o) : $(LIBF)/dyn3d/bernoui.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a bernoui.o ; \
	$(COMPILE) $(LIBF)/dyn3d/bernoui.F ; \
	$(AR) r $(LIBO)/libdyn3d.a bernoui.o ; $(RM) bernoui.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(caldyn0.o) : $(LIBF)/dyn3d/caldyn0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a caldyn0.o ; \
	$(COMPILE) $(LIBF)/dyn3d/caldyn0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a caldyn0.o ; $(RM) caldyn0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(caldyn.o) : $(LIBF)/dyn3d/caldyn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a caldyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d/caldyn.F ; \
	$(AR) r $(LIBO)/libdyn3d.a caldyn.o ; $(RM) caldyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(calfis.o) : $(LIBF)/dyn3d/calfis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/advtrac.h \
$(LIBO)/libdyn3d.a(infotrac.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a calfis.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/calfis.F ; \
	$(AR) r $(LIBO)/libdyn3d.a calfis.o ; $(RM) calfis.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(coefpoly.o) : $(LIBF)/dyn3d/coefpoly.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a coefpoly.o ; \
	$(COMPILE) $(LIBF)/dyn3d/coefpoly.F ; \
	$(AR) r $(LIBO)/libdyn3d.a coefpoly.o ; $(RM) coefpoly.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(convflu.o) : $(LIBF)/dyn3d/convflu.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a convflu.o ; \
	$(COMPILE) $(LIBF)/dyn3d/convflu.F ; \
	$(AR) r $(LIBO)/libdyn3d.a convflu.o ; $(RM) convflu.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(convmas.o) : $(LIBF)/dyn3d/convmas.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a convmas.o ; \
	$(COMPILE) $(LIBF)/dyn3d/convmas.F ; \
	$(AR) r $(LIBO)/libdyn3d.a convmas.o ; $(RM) convmas.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(coordij.o) : $(LIBF)/dyn3d/coordij.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a coordij.o ; \
	$(COMPILE) $(LIBF)/dyn3d/coordij.F ; \
	$(AR) r $(LIBO)/libdyn3d.a coordij.o ; $(RM) coordij.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(covcont.o) : $(LIBF)/dyn3d/covcont.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a covcont.o ; \
	$(COMPILE) $(LIBF)/dyn3d/covcont.F ; \
	$(AR) r $(LIBO)/libdyn3d.a covcont.o ; $(RM) covcont.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(covnat.o) : $(LIBF)/dyn3d/covnat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a covnat.o ; \
	$(COMPILE) $(LIBF)/dyn3d/covnat.F ; \
	$(AR) r $(LIBO)/libdyn3d.a covnat.o ; $(RM) covnat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(defrun_new.o) : $(LIBF)/dyn3d/defrun_new.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/comdissnew.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3d.a(sponge_mod.o) \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a defrun_new.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/defrun_new.F ; \
	$(AR) r $(LIBO)/libdyn3d.a defrun_new.o ; $(RM) defrun_new.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(diagadv.o) : $(LIBF)/dyn3d/diagadv.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a diagadv.o ; \
	$(COMPILE) $(LIBF)/dyn3d/diagadv.F ; \
	$(AR) r $(LIBO)/libdyn3d.a diagadv.o ; $(RM) diagadv.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dissip.o) : $(LIBF)/dyn3d/dissip.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comdissnew.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dissip.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dissip.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dissip.o ; $(RM) dissip.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(disvert.o) : $(LIBF)/dyn3d/disvert.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a disvert.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/disvert.F ; \
	$(AR) r $(LIBO)/libdyn3d.a disvert.o ; $(RM) disvert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(diverg.o) : $(LIBF)/dyn3d/diverg.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a diverg.o ; \
	$(COMPILE) $(LIBF)/dyn3d/diverg.F ; \
	$(AR) r $(LIBO)/libdyn3d.a diverg.o ; $(RM) diverg.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divergf.o) : $(LIBF)/dyn3d/divergf.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divergf.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divergf.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divergf.o ; $(RM) divergf.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(diverg_gam.o) : $(LIBF)/dyn3d/diverg_gam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a diverg_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/diverg_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a diverg_gam.o ; $(RM) diverg_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divergst.o) : $(LIBF)/dyn3d/divergst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divergst.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divergst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divergst.o ; $(RM) divergst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divgrad22.o) : $(LIBF)/dyn3d/divgrad22.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divgrad22.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divgrad22.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divgrad22.o ; $(RM) divgrad22.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divgrad2.o) : $(LIBF)/dyn3d/divgrad2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divgrad2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divgrad2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divgrad2.o ; $(RM) divgrad2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divgrad.o) : $(LIBF)/dyn3d/divgrad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divgrad.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divgrad.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divgrad.o ; $(RM) divgrad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dteta1.o) : $(LIBF)/dyn3d/dteta1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dteta1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dteta1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dteta1.o ; $(RM) dteta1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dudv1.o) : $(LIBF)/dyn3d/dudv1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dudv1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dudv1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dudv1.o ; $(RM) dudv1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dudv2.o) : $(LIBF)/dyn3d/dudv2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dudv2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dudv2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dudv2.o ; $(RM) dudv2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dump2d.o) : $(LIBF)/dyn3d/dump2d.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dump2d.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dump2d.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dump2d.o ; $(RM) dump2d.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dynetat0.o) : $(LIBF)/dyn3d/dynetat0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/logic.h \
$(LIBO)/libdyn3d.a(infotrac.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dynetat0.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/dynetat0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dynetat0.o ; $(RM) dynetat0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dynredem.o) : $(LIBF)/dyn3d/dynredem.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/advtrac.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libdyn3d.a(infotrac.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dynredem.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/dynredem.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dynredem.o ; $(RM) dynredem.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(enercin.o) : $(LIBF)/dyn3d/enercin.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a enercin.o ; \
	$(COMPILE) $(LIBF)/dyn3d/enercin.F ; \
	$(AR) r $(LIBO)/libdyn3d.a enercin.o ; $(RM) enercin.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(exner.o) : $(LIBF)/dyn3d/exner.F \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a exner.o ; \
	$(COMPILE) $(LIBF)/dyn3d/exner.F ; \
	$(AR) r $(LIBO)/libdyn3d.a exner.o ; $(RM) exner.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(exner_hyb.o) : $(LIBF)/dyn3d/exner_hyb.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a exner_hyb.o ; \
	$(COMPILE) $(LIBF)/dyn3d/exner_hyb.F ; \
	$(AR) r $(LIBO)/libdyn3d.a exner_hyb.o ; $(RM) exner_hyb.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(flumass.o) : $(LIBF)/dyn3d/flumass.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a flumass.o ; \
	$(COMPILE) $(LIBF)/dyn3d/flumass.F ; \
	$(AR) r $(LIBO)/libdyn3d.a flumass.o ; $(RM) flumass.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fluxstoke.o) : $(LIBF)/dyn3d/fluxstoke.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/tracstoke.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fluxstoke.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fluxstoke.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fluxstoke.o ; $(RM) fluxstoke.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxhyp.o) : $(LIBF)/dyn3d/fxhyp.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxhyp.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxhyp.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxhyp.o ; $(RM) fxhyp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxy.o) : $(LIBF)/dyn3d/fxy.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxy_new.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxy.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxy.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxy.o ; $(RM) fxy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxyhyper.o) : $(LIBF)/dyn3d/fxyhyper.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxyhyper.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxyhyper.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxyhyper.o ; $(RM) fxyhyper.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxysinus.o) : $(LIBF)/dyn3d/fxysinus.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxy_sin.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxysinus.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxysinus.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxysinus.o ; $(RM) fxysinus.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fyhyp.o) : $(LIBF)/dyn3d/fyhyp.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fyhyp.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fyhyp.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fyhyp.o ; $(RM) fyhyp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(geopot.o) : $(LIBF)/dyn3d/geopot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a geopot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/geopot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a geopot.o ; $(RM) geopot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grad.o) : $(LIBF)/dyn3d/grad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grad.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grad.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grad.o ; $(RM) grad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gradiv2.o) : $(LIBF)/dyn3d/gradiv2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gradiv2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gradiv2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gradiv2.o ; $(RM) gradiv2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gradiv.o) : $(LIBF)/dyn3d/gradiv.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gradiv.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gradiv.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gradiv.o ; $(RM) gradiv.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_dyn_fi.o) : $(LIBF)/dyn3d/gr_dyn_fi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_dyn_fi.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_dyn_fi.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_dyn_fi.o ; $(RM) gr_dyn_fi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_fi_dyn.o) : $(LIBF)/dyn3d/gr_fi_dyn.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_fi_dyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_fi_dyn.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_fi_dyn.o ; $(RM) gr_fi_dyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grid_atob.o) : $(LIBF)/dyn3d/grid_atob.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grid_atob.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grid_atob.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grid_atob.o ; $(RM) grid_atob.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grid_noro1.o) : $(LIBF)/dyn3d/grid_noro1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grid_noro1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grid_noro1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grid_noro1.o ; $(RM) grid_noro1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grid_noro.o) : $(LIBF)/dyn3d/grid_noro.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grid_noro.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grid_noro.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grid_noro.o ; $(RM) grid_noro.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_int_dyn.o) : $(LIBF)/dyn3d/gr_int_dyn.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_int_dyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_int_dyn.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_int_dyn.o ; $(RM) gr_int_dyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(groupe.o) : $(LIBF)/dyn3d/groupe.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a groupe.o ; \
	$(COMPILE) $(LIBF)/dyn3d/groupe.F ; \
	$(AR) r $(LIBO)/libdyn3d.a groupe.o ; $(RM) groupe.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(groupeun.o) : $(LIBF)/dyn3d/groupeun.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a groupeun.o ; \
	$(COMPILE) $(LIBF)/dyn3d/groupeun.F ; \
	$(AR) r $(LIBO)/libdyn3d.a groupeun.o ; $(RM) groupeun.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_u_scal.o) : $(LIBF)/dyn3d/gr_u_scal.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_u_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_u_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_u_scal.o ; $(RM) gr_u_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_v_scal.o) : $(LIBF)/dyn3d/gr_v_scal.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_v_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_v_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_v_scal.o ; $(RM) gr_v_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(heavyside.o) : $(LIBF)/dyn3d/heavyside.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a heavyside.o ; \
	$(COMPILE) $(LIBF)/dyn3d/heavyside.F ; \
	$(AR) r $(LIBO)/libdyn3d.a heavyside.o ; $(RM) heavyside.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniav.o) : $(LIBF)/dyn3d/iniav.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comav.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a iniav.o ; \
	$(COMPILE) $(LIBF)/dyn3d/iniav.F ; \
	$(AR) r $(LIBO)/libdyn3d.a iniav.o ; $(RM) iniav.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniconst.o) : $(LIBF)/dyn3d/iniconst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a iniconst.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/iniconst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a iniconst.o ; $(RM) iniconst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(inidissip.o) : $(LIBF)/dyn3d/inidissip.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a inidissip.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/inidissip.F ; \
	$(AR) r $(LIBO)/libdyn3d.a inidissip.o ; $(RM) inidissip.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(inigeom.o) : $(LIBF)/dyn3d/inigeom.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comdissnew.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a inigeom.o ; \
	$(COMPILE) $(LIBF)/dyn3d/inigeom.F ; \
	$(AR) r $(LIBO)/libdyn3d.a inigeom.o ; $(RM) inigeom.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniinterp_h.o) : $(LIBF)/dyn3d/iniinterp_h.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a iniinterp_h.o ; \
	$(COMPILE) $(LIBF)/dyn3d/iniinterp_h.F ; \
	$(AR) r $(LIBO)/libdyn3d.a iniinterp_h.o ; $(RM) iniinterp_h.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(initial0.o) : $(LIBF)/dyn3d/initial0.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a initial0.o ; \
	$(COMPILE) $(LIBF)/dyn3d/initial0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a initial0.o ; $(RM) initial0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(integrd.o) : $(LIBF)/dyn3d/integrd.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a integrd.o ; \
	$(COMPILE) $(LIBF)/dyn3d/integrd.F ; \
	$(AR) r $(LIBO)/libdyn3d.a integrd.o ; $(RM) integrd.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(interp_horiz.o) : $(LIBF)/dyn3d/interp_horiz.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a interp_horiz.o ; \
	$(COMPILE) $(LIBF)/dyn3d/interp_horiz.F ; \
	$(AR) r $(LIBO)/libdyn3d.a interp_horiz.o ; $(RM) interp_horiz.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(interp_vert.o) : $(LIBF)/dyn3d/interp_vert.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a interp_vert.o ; \
	$(COMPILE) $(LIBF)/dyn3d/interp_vert.F ; \
	$(AR) r $(LIBO)/libdyn3d.a interp_vert.o ; $(RM) interp_vert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien.o) : $(LIBF)/dyn3d/laplacien.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien.o ; $(RM) laplacien.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien_gam.o) : $(LIBF)/dyn3d/laplacien_gam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien_gam.o ; $(RM) laplacien_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien_rot.o) : $(LIBF)/dyn3d/laplacien_rot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien_rot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien_rot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien_rot.o ; $(RM) laplacien_rot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien_rotgam.o) : $(LIBF)/dyn3d/laplacien_rotgam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien_rotgam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien_rotgam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien_rotgam.o ; $(RM) laplacien_rotgam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(lectflux.o) : $(LIBF)/dyn3d/lectflux.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/tracstoke.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a lectflux.o ; \
	$(COMPILE) $(LIBF)/dyn3d/lectflux.F ; \
	$(AR) r $(LIBO)/libdyn3d.a lectflux.o ; $(RM) lectflux.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(massbar.o) : $(LIBF)/dyn3d/massbar.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a massbar.o ; \
	$(COMPILE) $(LIBF)/dyn3d/massbar.F ; \
	$(AR) r $(LIBO)/libdyn3d.a massbar.o ; $(RM) massbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(massbarxy.o) : $(LIBF)/dyn3d/massbarxy.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a massbarxy.o ; \
	$(COMPILE) $(LIBF)/dyn3d/massbarxy.F ; \
	$(AR) r $(LIBO)/libdyn3d.a massbarxy.o ; $(RM) massbarxy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(massdair.o) : $(LIBF)/dyn3d/massdair.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a massdair.o ; \
	$(COMPILE) $(LIBF)/dyn3d/massdair.F ; \
	$(AR) r $(LIBO)/libdyn3d.a massdair.o ; $(RM) massdair.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(minmax.o) : $(LIBF)/dyn3d/minmax.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a minmax.o ; \
	$(COMPILE) $(LIBF)/dyn3d/minmax.F ; \
	$(AR) r $(LIBO)/libdyn3d.a minmax.o ; $(RM) minmax.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(multipl.o) : $(LIBF)/dyn3d/multipl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a multipl.o ; \
	$(COMPILE) $(LIBF)/dyn3d/multipl.F ; \
	$(AR) r $(LIBO)/libdyn3d.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(multscal.o) : $(LIBF)/dyn3d/multscal.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a multscal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/multscal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a multscal.o ; $(RM) multscal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgrad.o) : $(LIBF)/dyn3d/nxgrad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgrad.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgrad.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgrad.o ; $(RM) nxgrad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgrad_gam.o) : $(LIBF)/dyn3d/nxgrad_gam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgrad_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgrad_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgrad_gam.o ; $(RM) nxgrad_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgradst.o) : $(LIBF)/dyn3d/nxgradst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgradst.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgradst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgradst.o ; $(RM) nxgradst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgraro2.o) : $(LIBF)/dyn3d/nxgraro2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgraro2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgraro2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgraro2.o ; $(RM) nxgraro2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgrarot.o) : $(LIBF)/dyn3d/nxgrarot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgrarot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgrarot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgrarot.o ; $(RM) nxgrarot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(pression.o) : $(LIBF)/dyn3d/pression.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a pression.o ; \
	$(COMPILE) $(LIBF)/dyn3d/pression.F ; \
	$(AR) r $(LIBO)/libdyn3d.a pression.o ; $(RM) pression.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(psextbar.o) : $(LIBF)/dyn3d/psextbar.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a psextbar.o ; \
	$(COMPILE) $(LIBF)/dyn3d/psextbar.F ; \
	$(AR) r $(LIBO)/libdyn3d.a psextbar.o ; $(RM) psextbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(ran1.o) : $(LIBF)/dyn3d/ran1.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a ran1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/ran1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a ran1.o ; $(RM) ran1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotat.o) : $(LIBF)/dyn3d/rotat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotat.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotat.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotat.o ; $(RM) rotat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotatf.o) : $(LIBF)/dyn3d/rotatf.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotatf.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotatf.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotatf.o ; $(RM) rotatf.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotat_nfil.o) : $(LIBF)/dyn3d/rotat_nfil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotat_nfil.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotat_nfil.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotat_nfil.o ; $(RM) rotat_nfil.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotatst.o) : $(LIBF)/dyn3d/rotatst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotatst.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotatst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotatst.o ; $(RM) rotatst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(scal_wind.o) : $(LIBF)/dyn3d/scal_wind.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a scal_wind.o ; \
	$(COMPILE) $(LIBF)/dyn3d/scal_wind.F ; \
	$(AR) r $(LIBO)/libdyn3d.a scal_wind.o ; $(RM) scal_wind.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(sortvarc0.o) : $(LIBF)/dyn3d/sortvarc0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a sortvarc0.o ; \
	$(COMPILE) $(LIBF)/dyn3d/sortvarc0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a sortvarc0.o ; $(RM) sortvarc0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(sortvarc.o) : $(LIBF)/dyn3d/sortvarc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a sortvarc.o ; \
	$(COMPILE) $(LIBF)/dyn3d/sortvarc.F ; \
	$(AR) r $(LIBO)/libdyn3d.a sortvarc.o ; $(RM) sortvarc.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(spline.o) : $(LIBF)/dyn3d/spline.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a spline.o ; \
	$(COMPILE) $(LIBF)/dyn3d/spline.F ; \
	$(AR) r $(LIBO)/libdyn3d.a spline.o ; $(RM) spline.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(splint.o) : $(LIBF)/dyn3d/splint.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a splint.o ; \
	$(COMPILE) $(LIBF)/dyn3d/splint.F ; \
	$(AR) r $(LIBO)/libdyn3d.a splint.o ; $(RM) splint.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(test_period.o) : $(LIBF)/dyn3d/test_period.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libdyn3d.a(infotrac.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a test_period.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/test_period.F ; \
	$(AR) r $(LIBO)/libdyn3d.a test_period.o ; $(RM) test_period.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(tourpot.o) : $(LIBF)/dyn3d/tourpot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a tourpot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/tourpot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a tourpot.o ; $(RM) tourpot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(traceur.o) : $(LIBF)/dyn3d/traceur.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a traceur.o ; \
	$(COMPILE) $(LIBF)/dyn3d/traceur.F ; \
	$(AR) r $(LIBO)/libdyn3d.a traceur.o ; $(RM) traceur.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(tracvl.o) : $(LIBF)/dyn3d/tracvl.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a tracvl.o ; \
	$(COMPILE) $(LIBF)/dyn3d/tracvl.F ; \
	$(AR) r $(LIBO)/libdyn3d.a tracvl.o ; $(RM) tracvl.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vanleer.o) : $(LIBF)/dyn3d/vanleer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vanleer.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vanleer.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vanleer.o ; $(RM) vanleer.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vitvert.o) : $(LIBF)/dyn3d/vitvert.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vitvert.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vitvert.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vitvert.o ; $(RM) vitvert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vlsplt.o) : $(LIBF)/dyn3d/vlsplt.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vlsplt.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vlsplt.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vlsplt.o ; $(RM) vlsplt.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vlspltqs.o) : $(LIBF)/dyn3d/vlspltqs.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vlspltqs.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vlspltqs.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vlspltqs.o ; $(RM) vlspltqs.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(wind_scal.o) : $(LIBF)/dyn3d/wind_scal.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a wind_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/wind_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a wind_scal.o ; $(RM) wind_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(control_mod.o) : $(LIBF)/dyn3d/control_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a control_mod.o ; \
	$(COMPTRU90) $(LIBF)/dyn3d/control_mod.F90 ; \
	mv $(MOD_LOC_DIR)/control_mod.$(MOD_SUFFIX) $(LIBO)/control_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libdyn3d.a control_mod.o ; $(RM) control_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(infotrac.o) : $(LIBF)/dyn3d/infotrac.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/advtrac.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a infotrac.o ; \
	$(COMPTRU90) $(LIBF)/dyn3d/infotrac.F90 ; \
	mv $(MOD_LOC_DIR)/infotrac.$(MOD_SUFFIX) $(LIBO)/infotrac.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libdyn3d.a infotrac.o ; $(RM) infotrac.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(mod_const_mpi.o) : $(LIBF)/dyn3d/mod_const_mpi.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a mod_const_mpi.o ; \
	$(COMPTRU90) $(LIBF)/dyn3d/mod_const_mpi.F90 ; \
	mv $(MOD_LOC_DIR)/mod_const_mpi.$(MOD_SUFFIX) $(LIBO)/mod_const_mpi.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libdyn3d.a mod_const_mpi.o ; $(RM) mod_const_mpi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(sponge_mod.o) : $(LIBF)/dyn3d/sponge_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissip.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a sponge_mod.o ; \
	$(COMPTRU90) $(LIBF)/dyn3d/sponge_mod.F90 ; \
	mv $(MOD_LOC_DIR)/sponge_mod.$(MOD_SUFFIX) $(LIBO)/sponge_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libdyn3d.a sponge_mod.o ; $(RM) sponge_mod.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory filtrez
#=======================================================================

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(acc.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(eigen.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(eigen_sort.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(filtreg.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(inifgn.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(inifilr.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(jacobi.o)

.PRECIOUS	: $(LIBO)/libfiltrez.a


# Compilation des membres de la bibliotheque libfiltrez.a

$(LIBO)/libfiltrez.a(acc.o) : $(LIBF)/filtrez/acc.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a acc.o ; \
	$(COMPILE) $(LIBF)/filtrez/acc.F ; \
	$(AR) r $(LIBO)/libfiltrez.a acc.o ; $(RM) acc.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(eigen.o) : $(LIBF)/filtrez/eigen.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a eigen.o ; \
	$(COMPILE) $(LIBF)/filtrez/eigen.F ; \
	$(AR) r $(LIBO)/libfiltrez.a eigen.o ; $(RM) eigen.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(eigen_sort.o) : $(LIBF)/filtrez/eigen_sort.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a eigen_sort.o ; \
	$(COMPILE) $(LIBF)/filtrez/eigen_sort.F ; \
	$(AR) r $(LIBO)/libfiltrez.a eigen_sort.o ; $(RM) eigen_sort.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(filtreg.o) : $(LIBF)/filtrez/filtreg.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/filtrez/parafilt.h \
$(LIBF)/filtrez/coefils.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a filtreg.o ; \
	$(COMPILE) $(LIBF)/filtrez/filtreg.F ; \
	$(AR) r $(LIBO)/libfiltrez.a filtreg.o ; $(RM) filtreg.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(inifgn.o) : $(LIBF)/filtrez/inifgn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/filtrez/coefils.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a inifgn.o ; \
	$(COMPILE) $(LIBF)/filtrez/inifgn.F ; \
	$(AR) r $(LIBO)/libfiltrez.a inifgn.o ; $(RM) inifgn.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(inifilr.o) : $(LIBF)/filtrez/inifilr.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/filtrez/parafilt.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/filtrez/coefils.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a inifilr.o ; \
	$(COMPILE) $(LIBF)/filtrez/inifilr.F ; \
	$(AR) r $(LIBO)/libfiltrez.a inifilr.o ; $(RM) inifilr.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(jacobi.o) : $(LIBF)/filtrez/jacobi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a jacobi.o ; \
	$(COMPILE) $(LIBF)/filtrez/jacobi.F ; \
	$(AR) r $(LIBO)/libfiltrez.a jacobi.o ; $(RM) jacobi.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phyhanbo
#=======================================================================

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(aerave.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(aerave_new.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(blackl.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(callsedim.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(convadj.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(cpdet.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(cvmgp.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(cvmgt.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(datareadnc.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(dsolver.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(dtridgl.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(evap.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(forceWCfn.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(gfluxi.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(gfluxv.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iniaerosol.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(ini_archive.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(inifis.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iniorbit.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(inistats.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(initracer.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iniwrite.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iniwrite_specIR.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iniwrite_specVI.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(interp_line.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(lagrange.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(lect_start_archive.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mucorr.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mufract.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(multipl.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(newsedim.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(orbite.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(profile.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(scopyi.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(sfluxi.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(sfluxv.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(soil.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(soil_settings.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(stelang.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(stellarlong.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(surface_nature.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(surfini.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(tabfi.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(tpindex.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(vdif_cd.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(vdifc.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(vdif_kc.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(vlz_fi.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(write_archive.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(writediagfi.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(writediagspecIR.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(writediagspecVI.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(writeg1d.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(zerophys.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(aeropacity.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(aeroptproperties.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(aerosol_mod.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(ave_stelspec.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(bilinearbig.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(bilinear.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(calc_cpp3d.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(calc_cpp_mugaz.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(calcenergy_kcm.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(calc_rayleigh.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(callcorrk.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(comdiurn_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(comgeomfi_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(comgeomphy.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(comsaison_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(comsoil_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(condense_cloud.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(cp_neutral.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(datafile_mod.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(def_var.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(dimphy.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(gases_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(gradients_kcm.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(hydrol.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iniphysiq.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(init_phys_lmdz.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iniwritesoil.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(interpolateH2H2.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(interpolateH2He.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(interpolateH2Ocont_CKD.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(interpolateH2Ocont_PPC.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(interpolateN2H2.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(interpolateN2N2.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(iostart.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(kcmprof_fn.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(largescale.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mass_redistribution.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mkstat.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_transfert.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_transfert.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(mod_phys_lmdz_transfert_para.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(moistadj.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(newtrelax.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(ocean_slab_mod.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(optci.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(optcv.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(params_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(phasechange_kcm.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(phyetat0.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(phyredem.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(physiq.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(planete_mod.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(planetwide_mod.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(radcommon_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(radii_mod.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(radinc_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(rain.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(rings.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(setspi.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(setspv.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(slab_ice_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(stokes.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(suaer_corrk.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(sugas_corrk.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(su_gases.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(surfdat_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(surf_heat_transp_mod.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(su_watercycle.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(totalcloudfrac.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(tracer_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(turbdiff.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(watercommon_h.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(watersat.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(watersat_grad.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(writediagsoil.o)

$(LIBO)/libphyhanbo.a :  $(LIBO)/libphyhanbo.a(wstats.o)

.PRECIOUS	: $(LIBO)/libphyhanbo.a


# Compilation des membres de la bibliotheque libphyhanbo.a

$(LIBO)/libphyhanbo.a(aerave.o) : $(LIBF)/phyhanbo/aerave.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a aerave.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/aerave.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a aerave.o ; $(RM) aerave.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(aerave_new.o) : $(LIBF)/phyhanbo/aerave_new.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a aerave_new.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/aerave_new.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a aerave_new.o ; $(RM) aerave_new.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(blackl.o) : $(LIBF)/phyhanbo/blackl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a blackl.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/blackl.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a blackl.o ; $(RM) blackl.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(callsedim.o) : $(LIBF)/phyhanbo/callsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radii_mod.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a callsedim.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/callsedim.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a callsedim.o ; $(RM) callsedim.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(convadj.o) : $(LIBF)/phyhanbo/convadj.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a convadj.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/convadj.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a convadj.o ; $(RM) convadj.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(cpdet.o) : $(LIBF)/phyhanbo/cpdet.F \
$(LIBF)/phyhanbo/cpdet.h \
$(LIBF)/phyhanbo/cpdet.h \
$(LIBF)/phyhanbo/cpdet.h \
$(LIBF)/phyhanbo/cpdet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a cpdet.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/cpdet.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a cpdet.o ; $(RM) cpdet.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(cvmgp.o) : $(LIBF)/phyhanbo/cvmgp.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a cvmgp.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/cvmgp.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(cvmgt.o) : $(LIBF)/phyhanbo/cvmgt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a cvmgt.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/cvmgt.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(datareadnc.o) : $(LIBF)/phyhanbo/datareadnc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxyprim.h \
$(LIBO)/libphyhanbo.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a datareadnc.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/datareadnc.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a datareadnc.o ; $(RM) datareadnc.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(dsolver.o) : $(LIBF)/phyhanbo/dsolver.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a dsolver.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/dsolver.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a dsolver.o ; $(RM) dsolver.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(dtridgl.o) : $(LIBF)/phyhanbo/dtridgl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a dtridgl.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/dtridgl.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a dtridgl.o ; $(RM) dtridgl.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(evap.o) : $(LIBF)/phyhanbo/evap.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a evap.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/evap.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a evap.o ; $(RM) evap.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(forceWCfn.o) : $(LIBF)/phyhanbo/forceWCfn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a forceWCfn.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/forceWCfn.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a forceWCfn.o ; $(RM) forceWCfn.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(gfluxi.o) : $(LIBF)/phyhanbo/gfluxi.F \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a gfluxi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/gfluxi.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a gfluxi.o ; $(RM) gfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(gfluxv.o) : $(LIBF)/phyhanbo/gfluxv.F \
$(LIBO)/libphyhanbo.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a gfluxv.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/gfluxv.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a gfluxv.o ; $(RM) gfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iniaerosol.o) : $(LIBF)/phyhanbo/iniaerosol.F \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iniaerosol.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/iniaerosol.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a iniaerosol.o ; $(RM) iniaerosol.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(ini_archive.o) : $(LIBF)/phyhanbo/ini_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a ini_archive.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/ini_archive.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a ini_archive.o ; $(RM) ini_archive.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(inifis.o) : $(LIBF)/phyhanbo/inifis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o) \
$(LIBO)/libphyhanbo.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo.a(planete_mod.o) \
$(LIBO)/libphyhanbo.a(planetwide_mod.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a inifis.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/inifis.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a inifis.o ; $(RM) inifis.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iniorbit.o) : $(LIBF)/phyhanbo/iniorbit.F \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iniorbit.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/iniorbit.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a iniorbit.o ; $(RM) iniorbit.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(inistats.o) : $(LIBF)/phyhanbo/inistats.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phyhanbo/statto.h \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a inistats.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/inistats.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a inistats.o ; $(RM) inistats.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(initracer.o) : $(LIBF)/phyhanbo/initracer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a initracer.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/initracer.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a initracer.o ; $(RM) initracer.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iniwrite.o) : $(LIBF)/phyhanbo/iniwrite.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iniwrite.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/iniwrite.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a iniwrite.o ; $(RM) iniwrite.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iniwrite_specIR.o) : $(LIBF)/phyhanbo/iniwrite_specIR.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iniwrite_specIR.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/iniwrite_specIR.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a iniwrite_specIR.o ; $(RM) iniwrite_specIR.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iniwrite_specVI.o) : $(LIBF)/phyhanbo/iniwrite_specVI.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iniwrite_specVI.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/iniwrite_specVI.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a iniwrite_specVI.o ; $(RM) iniwrite_specVI.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(interp_line.o) : $(LIBF)/phyhanbo/interp_line.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a interp_line.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/interp_line.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a interp_line.o ; $(RM) interp_line.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(lagrange.o) : $(LIBF)/phyhanbo/lagrange.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a lagrange.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/lagrange.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a lagrange.o ; $(RM) lagrange.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(lect_start_archive.o) : $(LIBF)/phyhanbo/lect_start_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a lect_start_archive.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/lect_start_archive.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a lect_start_archive.o ; $(RM) lect_start_archive.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mucorr.o) : $(LIBF)/phyhanbo/mucorr.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mucorr.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/mucorr.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a mucorr.o ; $(RM) mucorr.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mufract.o) : $(LIBF)/phyhanbo/mufract.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mufract.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/mufract.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a mufract.o ; $(RM) mufract.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(multipl.o) : $(LIBF)/phyhanbo/multipl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a multipl.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/multipl.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(newsedim.o) : $(LIBF)/phyhanbo/newsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a newsedim.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/newsedim.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a newsedim.o ; $(RM) newsedim.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(orbite.o) : $(LIBF)/phyhanbo/orbite.F \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a orbite.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/orbite.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a orbite.o ; $(RM) orbite.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(profile.o) : $(LIBF)/phyhanbo/profile.F \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a profile.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/profile.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a profile.o ; $(RM) profile.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(scopyi.o) : $(LIBF)/phyhanbo/scopyi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a scopyi.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/scopyi.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(sfluxi.o) : $(LIBF)/phyhanbo/sfluxi.F \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a sfluxi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/sfluxi.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a sfluxi.o ; $(RM) sfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(sfluxv.o) : $(LIBF)/phyhanbo/sfluxv.F \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a sfluxv.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/sfluxv.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a sfluxv.o ; $(RM) sfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(soil.o) : $(LIBF)/phyhanbo/soil.F \
$(LIBO)/libphyhanbo.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a soil.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/soil.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a soil.o ; $(RM) soil.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(soil_settings.o) : $(LIBF)/phyhanbo/soil_settings.F \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(iostart.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a soil_settings.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/soil_settings.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a soil_settings.o ; $(RM) soil_settings.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(stelang.o) : $(LIBF)/phyhanbo/stelang.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a stelang.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/stelang.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a stelang.o ; $(RM) stelang.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(stellarlong.o) : $(LIBF)/phyhanbo/stellarlong.F \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a stellarlong.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/stellarlong.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a stellarlong.o ; $(RM) stellarlong.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(surface_nature.o) : $(LIBF)/phyhanbo/surface_nature.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a surface_nature.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/surface_nature.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a surface_nature.o ; $(RM) surface_nature.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(surfini.o) : $(LIBF)/phyhanbo/surfini.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(planetwide_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a surfini.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/surfini.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a surfini.o ; $(RM) surfini.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(tabfi.o) : $(LIBF)/phyhanbo/tabfi.F \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(iostart.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a tabfi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/tabfi.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a tabfi.o ; $(RM) tabfi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(tpindex.o) : $(LIBF)/phyhanbo/tpindex.F \
$(LIBO)/libphyhanbo.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a tpindex.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/tpindex.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a tpindex.o ; $(RM) tpindex.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(vdif_cd.o) : $(LIBF)/phyhanbo/vdif_cd.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a vdif_cd.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/vdif_cd.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a vdif_cd.o ; $(RM) vdif_cd.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(vdifc.o) : $(LIBF)/phyhanbo/vdifc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a vdifc.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/vdifc.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a vdifc.o ; $(RM) vdifc.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(vdif_kc.o) : $(LIBF)/phyhanbo/vdif_kc.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a vdif_kc.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/vdif_kc.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a vdif_kc.o ; $(RM) vdif_kc.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(vlz_fi.o) : $(LIBF)/phyhanbo/vlz_fi.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a vlz_fi.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/vlz_fi.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a vlz_fi.o ; $(RM) vlz_fi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(write_archive.o) : $(LIBF)/phyhanbo/write_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBO)/libphyhanbo.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a write_archive.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/write_archive.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a write_archive.o ; $(RM) write_archive.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(writediagfi.o) : $(LIBF)/phyhanbo/writediagfi.F \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a writediagfi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/writediagfi.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a writediagfi.o ; $(RM) writediagfi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(writediagspecIR.o) : $(LIBF)/phyhanbo/writediagspecIR.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a writediagspecIR.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/writediagspecIR.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a writediagspecIR.o ; $(RM) writediagspecIR.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(writediagspecVI.o) : $(LIBF)/phyhanbo/writediagspecVI.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a writediagspecVI.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo/writediagspecVI.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a writediagspecVI.o ; $(RM) writediagspecVI.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(writeg1d.o) : $(LIBF)/phyhanbo/writeg1d.F \
$(LIBF)/phyhanbo/comg1d.h \
$(LIBF)/phyhanbo/comg1d.h \
$(LIBF)/phyhanbo/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a writeg1d.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/writeg1d.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a writeg1d.o ; $(RM) writeg1d.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(zerophys.o) : $(LIBF)/phyhanbo/zerophys.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a zerophys.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/zerophys.F ; \
	$(AR) r $(LIBO)/libphyhanbo.a zerophys.o ; $(RM) zerophys.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(aeropacity.o) : $(LIBF)/phyhanbo/aeropacity.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a aeropacity.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/aeropacity.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a aeropacity.o ; $(RM) aeropacity.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(aeroptproperties.o) : $(LIBF)/phyhanbo/aeroptproperties.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a aeroptproperties.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/aeroptproperties.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a aeroptproperties.o ; $(RM) aeroptproperties.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(aerosol_mod.o) : $(LIBF)/phyhanbo/aerosol_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a aerosol_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/aerosol_mod.F90 ; \
	mv $(MOD_LOC_DIR)/aerosol_mod.$(MOD_SUFFIX) $(LIBO)/aerosol_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a aerosol_mod.o ; $(RM) aerosol_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(ave_stelspec.o) : $(LIBF)/phyhanbo/ave_stelspec.F90 \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a ave_stelspec.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/ave_stelspec.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a ave_stelspec.o ; $(RM) ave_stelspec.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(bilinearbig.o) : $(LIBF)/phyhanbo/bilinearbig.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a bilinearbig.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/bilinearbig.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a bilinearbig.o ; $(RM) bilinearbig.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(bilinear.o) : $(LIBF)/phyhanbo/bilinear.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a bilinear.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/bilinear.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a bilinear.o ; $(RM) bilinear.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(calc_cpp3d.o) : $(LIBF)/phyhanbo/calc_cpp3d.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/phyhanbo/cpdet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a calc_cpp3d.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/calc_cpp3d.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a calc_cpp3d.o ; $(RM) calc_cpp3d.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(calc_cpp_mugaz.o) : $(LIBF)/phyhanbo/calc_cpp_mugaz.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a calc_cpp_mugaz.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/calc_cpp_mugaz.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a calc_cpp_mugaz.o ; $(RM) calc_cpp_mugaz.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(calcenergy_kcm.o) : $(LIBF)/phyhanbo/calcenergy_kcm.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(params_h.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a calcenergy_kcm.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/calcenergy_kcm.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a calcenergy_kcm.o ; $(RM) calcenergy_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(calc_rayleigh.o) : $(LIBF)/phyhanbo/calc_rayleigh.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a calc_rayleigh.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/calc_rayleigh.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a calc_rayleigh.o ; $(RM) calc_rayleigh.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(callcorrk.o) : $(LIBF)/phyhanbo/callcorrk.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo.a(gases_h.o) \
$(LIBO)/libphyhanbo.a(radii_mod.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a callcorrk.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/callcorrk.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a callcorrk.o ; $(RM) callcorrk.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(comdiurn_h.o) : $(LIBF)/phyhanbo/comdiurn_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a comdiurn_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/comdiurn_h.F90 ; \
	mv $(MOD_LOC_DIR)/comdiurn_h.$(MOD_SUFFIX) $(LIBO)/comdiurn_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a comdiurn_h.o ; $(RM) comdiurn_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(comgeomfi_h.o) : $(LIBF)/phyhanbo/comgeomfi_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a comgeomfi_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/comgeomfi_h.F90 ; \
	mv $(MOD_LOC_DIR)/comgeomfi_h.$(MOD_SUFFIX) $(LIBO)/comgeomfi_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a comgeomfi_h.o ; $(RM) comgeomfi_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(comgeomphy.o) : $(LIBF)/phyhanbo/comgeomphy.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a comgeomphy.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/comgeomphy.F90 ; \
	mv $(MOD_LOC_DIR)/comgeomphy.$(MOD_SUFFIX) $(LIBO)/comgeomphy.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a comgeomphy.o ; $(RM) comgeomphy.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(comsaison_h.o) : $(LIBF)/phyhanbo/comsaison_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a comsaison_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/comsaison_h.F90 ; \
	mv $(MOD_LOC_DIR)/comsaison_h.$(MOD_SUFFIX) $(LIBO)/comsaison_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a comsaison_h.o ; $(RM) comsaison_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(comsoil_h.o) : $(LIBF)/phyhanbo/comsoil_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a comsoil_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/comsoil_h.F90 ; \
	mv $(MOD_LOC_DIR)/comsoil_h.$(MOD_SUFFIX) $(LIBO)/comsoil_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a comsoil_h.o ; $(RM) comsoil_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(condense_cloud.o) : $(LIBF)/phyhanbo/condense_cloud.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(gases_h.o) \
$(LIBO)/libphyhanbo.a(radii_mod.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a condense_cloud.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/condense_cloud.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a condense_cloud.o ; $(RM) condense_cloud.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(cp_neutral.o) : $(LIBF)/phyhanbo/cp_neutral.F90 \
$(LIBO)/libphyhanbo.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a cp_neutral.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/cp_neutral.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a cp_neutral.o ; $(RM) cp_neutral.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(datafile_mod.o) : $(LIBF)/phyhanbo/datafile_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a datafile_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/datafile_mod.F90 ; \
	mv $(MOD_LOC_DIR)/datafile_mod.$(MOD_SUFFIX) $(LIBO)/datafile_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a datafile_mod.o ; $(RM) datafile_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(def_var.o) : $(LIBF)/phyhanbo/def_var.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a def_var.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/def_var.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a def_var.o ; $(RM) def_var.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(dimphy.o) : $(LIBF)/phyhanbo/dimphy.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a dimphy.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/dimphy.F90 ; \
	mv $(MOD_LOC_DIR)/dimphy.$(MOD_SUFFIX) $(LIBO)/dimphy.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a dimphy.o ; $(RM) dimphy.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(gases_h.o) : $(LIBF)/phyhanbo/gases_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a gases_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/gases_h.F90 ; \
	mv $(MOD_LOC_DIR)/gases_h.$(MOD_SUFFIX) $(LIBO)/gases_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a gases_h.o ; $(RM) gases_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(gradients_kcm.o) : $(LIBF)/phyhanbo/gradients_kcm.F90 \
$(LIBO)/libphyhanbo.a(params_h.o) \
$(LIBO)/libphyhanbo.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a gradients_kcm.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/gradients_kcm.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a gradients_kcm.o ; $(RM) gradients_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(hydrol.o) : $(LIBF)/phyhanbo/hydrol.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a hydrol.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/hydrol.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a hydrol.o ; $(RM) hydrol.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iniphysiq.o) : $(LIBF)/phyhanbo/iniphysiq.F90 \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(comgeomphy.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iniphysiq.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/iniphysiq.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a iniphysiq.o ; $(RM) iniphysiq.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(init_phys_lmdz.o) : $(LIBF)/phyhanbo/init_phys_lmdz.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a init_phys_lmdz.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/init_phys_lmdz.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a init_phys_lmdz.o ; $(RM) init_phys_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iniwritesoil.o) : $(LIBF)/phyhanbo/iniwritesoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libphyhanbo.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iniwritesoil.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/iniwritesoil.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a iniwritesoil.o ; $(RM) iniwritesoil.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(interpolateH2H2.o) : $(LIBF)/phyhanbo/interpolateH2H2.F90 \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a interpolateH2H2.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/interpolateH2H2.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a interpolateH2H2.o ; $(RM) interpolateH2H2.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(interpolateH2He.o) : $(LIBF)/phyhanbo/interpolateH2He.F90 \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a interpolateH2He.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/interpolateH2He.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a interpolateH2He.o ; $(RM) interpolateH2He.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(interpolateH2Ocont_CKD.o) : $(LIBF)/phyhanbo/interpolateH2Ocont_CKD.F90 \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a interpolateH2Ocont_CKD.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/interpolateH2Ocont_CKD.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a interpolateH2Ocont_CKD.o ; $(RM) interpolateH2Ocont_CKD.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(interpolateH2Ocont_PPC.o) : $(LIBF)/phyhanbo/interpolateH2Ocont_PPC.F90 \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a interpolateH2Ocont_PPC.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/interpolateH2Ocont_PPC.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a interpolateH2Ocont_PPC.o ; $(RM) interpolateH2Ocont_PPC.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(interpolateN2H2.o) : $(LIBF)/phyhanbo/interpolateN2H2.F90 \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a interpolateN2H2.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/interpolateN2H2.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a interpolateN2H2.o ; $(RM) interpolateN2H2.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(interpolateN2N2.o) : $(LIBF)/phyhanbo/interpolateN2N2.F90 \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a interpolateN2N2.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/interpolateN2N2.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a interpolateN2N2.o ; $(RM) interpolateN2N2.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(iostart.o) : $(LIBF)/phyhanbo/iostart.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a iostart.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/iostart.F90 ; \
	mv $(MOD_LOC_DIR)/iostart.$(MOD_SUFFIX) $(LIBO)/iostart.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a iostart.o ; $(RM) iostart.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(kcmprof_fn.o) : $(LIBF)/phyhanbo/kcmprof_fn.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(params_h.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a kcmprof_fn.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/kcmprof_fn.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a kcmprof_fn.o ; $(RM) kcmprof_fn.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(largescale.o) : $(LIBF)/phyhanbo/largescale.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a largescale.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/largescale.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a largescale.o ; $(RM) largescale.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mass_redistribution.o) : $(LIBF)/phyhanbo/mass_redistribution.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o) \
$(LIBO)/libphyhanbo.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mass_redistribution.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mass_redistribution.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a mass_redistribution.o ; $(RM) mass_redistribution.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mkstat.o) : $(LIBF)/phyhanbo/mkstat.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/statto.h \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mkstat.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mkstat.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a mkstat.o ; $(RM) mkstat.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) : $(LIBF)/phyhanbo/mod_grid_phy_lmdz.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mod_grid_phy_lmdz.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mod_grid_phy_lmdz.F90 ; \
	mv $(MOD_LOC_DIR)/mod_grid_phy_lmdz.$(MOD_SUFFIX) $(LIBO)/mod_grid_phy_lmdz.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a mod_grid_phy_lmdz.o ; $(RM) mod_grid_phy_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) : $(LIBF)/phyhanbo/mod_phys_lmdz_mpi_data.F90 \
$(LIBO)/libdyn3d.a(mod_const_mpi.o) \
$(LIBO)/libdyn3d.a(mod_const_mpi.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mod_phys_lmdz_mpi_data.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mod_phys_lmdz_mpi_data.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_mpi_data.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_mpi_data.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a mod_phys_lmdz_mpi_data.o ; $(RM) mod_phys_lmdz_mpi_data.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_transfert.o) : $(LIBF)/phyhanbo/mod_phys_lmdz_mpi_transfert.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mod_phys_lmdz_mpi_transfert.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mod_phys_lmdz_mpi_transfert.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_mpi_transfert.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_mpi_transfert.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a mod_phys_lmdz_mpi_transfert.o ; $(RM) mod_phys_lmdz_mpi_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) : $(LIBF)/phyhanbo/mod_phys_lmdz_omp_data.F90 \
$(LIBO)/libphyhanbo.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mod_phys_lmdz_omp_data.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mod_phys_lmdz_omp_data.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_omp_data.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_omp_data.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a mod_phys_lmdz_omp_data.o ; $(RM) mod_phys_lmdz_omp_data.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_transfert.o) : $(LIBF)/phyhanbo/mod_phys_lmdz_omp_transfert.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mod_phys_lmdz_omp_transfert.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mod_phys_lmdz_omp_transfert.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_omp_transfert.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_omp_transfert.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a mod_phys_lmdz_omp_transfert.o ; $(RM) mod_phys_lmdz_omp_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) : $(LIBF)/phyhanbo/mod_phys_lmdz_para.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_transfert_para.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mod_phys_lmdz_para.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mod_phys_lmdz_para.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_para.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_para.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a mod_phys_lmdz_para.o ; $(RM) mod_phys_lmdz_para.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(mod_phys_lmdz_transfert_para.o) : $(LIBF)/phyhanbo/mod_phys_lmdz_transfert_para.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_transfert.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_omp_transfert.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a mod_phys_lmdz_transfert_para.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/mod_phys_lmdz_transfert_para.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_transfert_para.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_transfert_para.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a mod_phys_lmdz_transfert_para.o ; $(RM) mod_phys_lmdz_transfert_para.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(moistadj.o) : $(LIBF)/phyhanbo/moistadj.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a moistadj.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/moistadj.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a moistadj.o ; $(RM) moistadj.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(newtrelax.o) : $(LIBF)/phyhanbo/newtrelax.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a newtrelax.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/newtrelax.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a newtrelax.o ; $(RM) newtrelax.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(ocean_slab_mod.o) : $(LIBF)/phyhanbo/ocean_slab_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(surf_heat_transp_mod.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a ocean_slab_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/ocean_slab_mod.F90 ; \
	mv $(MOD_LOC_DIR)/ocean_slab_mod.$(MOD_SUFFIX) $(LIBO)/ocean_slab_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a ocean_slab_mod.o ; $(RM) ocean_slab_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(optci.o) : $(LIBF)/phyhanbo/optci.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a optci.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/optci.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a optci.o ; $(RM) optci.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(optcv.o) : $(LIBF)/phyhanbo/optcv.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a optcv.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/optcv.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a optcv.o ; $(RM) optcv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(params_h.o) : $(LIBF)/phyhanbo/params_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a params_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/params_h.F90 ; \
	mv $(MOD_LOC_DIR)/params_h.$(MOD_SUFFIX) $(LIBO)/params_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a params_h.o ; $(RM) params_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(phasechange_kcm.o) : $(LIBF)/phyhanbo/phasechange_kcm.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a phasechange_kcm.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/phasechange_kcm.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a phasechange_kcm.o ; $(RM) phasechange_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(phyetat0.o) : $(LIBF)/phyhanbo/phyetat0.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(iostart.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a phyetat0.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/phyetat0.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a phyetat0.o ; $(RM) phyetat0.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(phyredem.o) : $(LIBF)/phyhanbo/phyredem.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(iostart.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(planete_mod.o) \
$(LIBO)/libphyhanbo.a(iostart.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a phyredem.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/phyredem.F90 ; \
	mv $(MOD_LOC_DIR)/phyredem.$(MOD_SUFFIX) $(LIBO)/phyredem.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a phyredem.o ; $(RM) phyredem.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(physiq.o) : $(LIBF)/phyhanbo/physiq.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(gases_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radii_mod.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo.a(comsaison_h.o) \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo.a(phyredem.o) \
$(LIBO)/libphyhanbo.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo.a(ocean_slab_mod.o) \
$(LIBO)/libphyhanbo.a(surf_heat_transp_mod.o) \
$(LIBO)/libphyhanbo.a(planetwide_mod.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a physiq.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/physiq.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a physiq.o ; $(RM) physiq.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(planete_mod.o) : $(LIBF)/phyhanbo/planete_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a planete_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/planete_mod.F90 ; \
	mv $(MOD_LOC_DIR)/planete_mod.$(MOD_SUFFIX) $(LIBO)/planete_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a planete_mod.o ; $(RM) planete_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(planetwide_mod.o) : $(LIBF)/phyhanbo/planetwide_mod.F90 \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo.a(dimphy.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a planetwide_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/planetwide_mod.F90 ; \
	mv $(MOD_LOC_DIR)/planetwide_mod.$(MOD_SUFFIX) $(LIBO)/planetwide_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a planetwide_mod.o ; $(RM) planetwide_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(radcommon_h.o) : $(LIBF)/phyhanbo/radcommon_h.F90 \
$(LIBO)/libphyhanbo.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a radcommon_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/radcommon_h.F90 ; \
	mv $(MOD_LOC_DIR)/radcommon_h.$(MOD_SUFFIX) $(LIBO)/radcommon_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a radcommon_h.o ; $(RM) radcommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(radii_mod.o) : $(LIBF)/phyhanbo/radii_mod.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a radii_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/radii_mod.F90 ; \
	mv $(MOD_LOC_DIR)/radii_mod.$(MOD_SUFFIX) $(LIBO)/radii_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a radii_mod.o ; $(RM) radii_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(radinc_h.o) : $(LIBF)/phyhanbo/radinc_h.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/bands.h \
$(LIBF)/phyhanbo/scatterers.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a radinc_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/radinc_h.F90 ; \
	mv $(MOD_LOC_DIR)/radinc_h.$(MOD_SUFFIX) $(LIBO)/radinc_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a radinc_h.o ; $(RM) radinc_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(rain.o) : $(LIBF)/phyhanbo/rain.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(radii_mod.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a rain.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/rain.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a rain.o ; $(RM) rain.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(rings.o) : $(LIBF)/phyhanbo/rings.F90 \
$(LIBO)/libphyhanbo.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a rings.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/rings.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a rings.o ; $(RM) rings.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(setspi.o) : $(LIBF)/phyhanbo/setspi.F90 \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a setspi.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/setspi.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a setspi.o ; $(RM) setspi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(setspv.o) : $(LIBF)/phyhanbo/setspv.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a setspv.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/setspv.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a setspv.o ; $(RM) setspv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(slab_ice_h.o) : $(LIBF)/phyhanbo/slab_ice_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a slab_ice_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/slab_ice_h.F90 ; \
	mv $(MOD_LOC_DIR)/slab_ice_h.$(MOD_SUFFIX) $(LIBO)/slab_ice_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a slab_ice_h.o ; $(RM) slab_ice_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(stokes.o) : $(LIBF)/phyhanbo/stokes.F90 \
$(LIBF)/phyhanbo/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a stokes.o ; \
	$(COMPILE) $(LIBF)/phyhanbo/stokes.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a stokes.o ; $(RM) stokes.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(suaer_corrk.o) : $(LIBF)/phyhanbo/suaer_corrk.F90 \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a suaer_corrk.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/suaer_corrk.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a suaer_corrk.o ; $(RM) suaer_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(sugas_corrk.o) : $(LIBF)/phyhanbo/sugas_corrk.F90 \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBO)/libphyhanbo.a(radinc_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(datafile_mod.o) \
$(LIBO)/libphyhanbo.a(gases_h.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a sugas_corrk.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/sugas_corrk.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a sugas_corrk.o ; $(RM) sugas_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(su_gases.o) : $(LIBF)/phyhanbo/su_gases.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a su_gases.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/su_gases.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a su_gases.o ; $(RM) su_gases.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(surfdat_h.o) : $(LIBF)/phyhanbo/surfdat_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a surfdat_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/surfdat_h.F90 ; \
	mv $(MOD_LOC_DIR)/surfdat_h.$(MOD_SUFFIX) $(LIBO)/surfdat_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a surfdat_h.o ; $(RM) surfdat_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(surf_heat_transp_mod.o) : $(LIBF)/phyhanbo/surf_heat_transp_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phyhanbo/comhdiff.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phyhanbo/comhdiff.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phyhanbo/comhdiff.h \
$(LIBO)/libphyhanbo.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a surf_heat_transp_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/surf_heat_transp_mod.F90 ; \
	mv $(MOD_LOC_DIR)/surf_heat_transp_mod.$(MOD_SUFFIX) $(LIBO)/surf_heat_transp_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a surf_heat_transp_mod.o ; $(RM) surf_heat_transp_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(su_watercycle.o) : $(LIBF)/phyhanbo/su_watercycle.F90 \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a su_watercycle.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/su_watercycle.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a su_watercycle.o ; $(RM) su_watercycle.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(totalcloudfrac.o) : $(LIBF)/phyhanbo/totalcloudfrac.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/comcstfi.h \
$(LIBF)/phyhanbo/callkeys.h \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a totalcloudfrac.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/totalcloudfrac.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a totalcloudfrac.o ; $(RM) totalcloudfrac.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(tracer_h.o) : $(LIBF)/phyhanbo/tracer_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a tracer_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/tracer_h.F90 ; \
	mv $(MOD_LOC_DIR)/tracer_h.$(MOD_SUFFIX) $(LIBO)/tracer_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a tracer_h.o ; $(RM) tracer_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(turbdiff.o) : $(LIBF)/phyhanbo/turbdiff.F90 \
$(LIBO)/libphyhanbo.a(watercommon_h.o) \
$(LIBO)/libphyhanbo.a(radcommon_h.o) \
$(LIBO)/libphyhanbo.a(surfdat_h.o) \
$(LIBO)/libphyhanbo.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a turbdiff.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/turbdiff.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a turbdiff.o ; $(RM) turbdiff.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(watercommon_h.o) : $(LIBF)/phyhanbo/watercommon_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a watercommon_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/watercommon_h.F90 ; \
	mv $(MOD_LOC_DIR)/watercommon_h.$(MOD_SUFFIX) $(LIBO)/watercommon_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo.a watercommon_h.o ; $(RM) watercommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(watersat.o) : $(LIBF)/phyhanbo/watersat.F90 \
$(LIBO)/libphyhanbo.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a watersat.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/watersat.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a watersat.o ; $(RM) watersat.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(watersat_grad.o) : $(LIBF)/phyhanbo/watersat_grad.F90 \
$(LIBO)/libphyhanbo.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a watersat_grad.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/watersat_grad.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a watersat_grad.o ; $(RM) watersat_grad.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(writediagsoil.o) : $(LIBF)/phyhanbo/writediagsoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libphyhanbo.a(comsoil_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a writediagsoil.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/writediagsoil.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a writediagsoil.o ; $(RM) writediagsoil.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo.a(wstats.o) : $(LIBF)/phyhanbo/wstats.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phyhanbo/statto.h \
$(LIBO)/libphyhanbo.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo.a wstats.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo/wstats.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo.a wstats.o ; $(RM) wstats.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phyhanbo1
#=======================================================================

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(aerave.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(aerave_new.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(blackl.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(callsedim.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(convadj.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(cpdet.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(cvmgp.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(cvmgt.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(datareadnc.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(dsolver.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(dtridgl.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(evap.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(forceWCfn.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(gfluxi.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(gfluxv.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iniaerosol.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(ini_archive.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(inifis.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iniorbit.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(inistats.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(initracer.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iniwrite.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iniwrite_specIR.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iniwrite_specVI.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(interp_line.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(lagrange.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(lect_start_archive.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mucorr.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mufract.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(multipl.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(newsedim.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(orbite.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(profile.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(scopyi.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(sfluxi.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(sfluxv.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(soil.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(soil_settings.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(stelang.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(stellarlong.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(surface_nature.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(surfini.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(tabfi.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(tpindex.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(vdif_cd.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(vdifc.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(vdif_kc.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(vlz_fi.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(write_archive.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(writediagfi.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(writediagspecIR.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(writediagspecVI.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(writeg1d.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(zerophys.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(aeropacity.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(aeroptproperties.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(aerosol_mod.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(ave_stelspec.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(bilinearbig.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(bilinear.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(calc_cpp3d.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(calc_cpp_mugaz.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(calcenergy_kcm.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(calc_rayleigh.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(callcorrk.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(comdiurn_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(comgeomfi_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(comgeomphy.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(comsaison_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(comsoil_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(condense_cloud.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(cp_neutral.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(datafile_mod.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(def_var.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(dimphy.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(gases_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(gradients_kcm.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(hydrol.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iniphysiq.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(init_phys_lmdz.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iniwritesoil.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(interpolateH2H2.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(interpolateH2He.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(interpolateH2Ocont_CKD.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(interpolateH2Ocont_PPC.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(interpolateN2H2.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(interpolateN2N2.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(iostart.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(kcmprof_fn.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(largescale.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mass_redistribution.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mkstat.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_transfert.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_transfert.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(mod_phys_lmdz_transfert_para.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(moistadj.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(newtrelax.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(ocean_slab_mod.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(optci.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(optcv.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(params_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(phasechange_kcm.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(phyetat0.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(phyredem.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(physiq.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(planete_mod.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(planetwide_mod.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(radcommon_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(radii_mod.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(radinc_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(rain.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(rings.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(setspi.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(setspv.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(slab_ice_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(stokes.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(suaer_corrk.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(sugas_corrk.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(su_gases.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(surfdat_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(surf_heat_transp_mod.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(su_watercycle.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(totalcloudfrac.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(tracer_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(turbdiff.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(watercommon_h.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(watersat.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(watersat_grad.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(writediagsoil.o)

$(LIBO)/libphyhanbo1.a :  $(LIBO)/libphyhanbo1.a(wstats.o)

.PRECIOUS	: $(LIBO)/libphyhanbo1.a


# Compilation des membres de la bibliotheque libphyhanbo1.a

$(LIBO)/libphyhanbo1.a(aerave.o) : $(LIBF)/phyhanbo1/aerave.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a aerave.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/aerave.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a aerave.o ; $(RM) aerave.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(aerave_new.o) : $(LIBF)/phyhanbo1/aerave_new.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a aerave_new.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/aerave_new.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a aerave_new.o ; $(RM) aerave_new.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(blackl.o) : $(LIBF)/phyhanbo1/blackl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a blackl.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/blackl.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a blackl.o ; $(RM) blackl.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(callsedim.o) : $(LIBF)/phyhanbo1/callsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radii_mod.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a callsedim.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/callsedim.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a callsedim.o ; $(RM) callsedim.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(convadj.o) : $(LIBF)/phyhanbo1/convadj.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a convadj.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/convadj.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a convadj.o ; $(RM) convadj.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(cpdet.o) : $(LIBF)/phyhanbo1/cpdet.F \
$(LIBF)/phyhanbo1/cpdet.h \
$(LIBF)/phyhanbo1/cpdet.h \
$(LIBF)/phyhanbo1/cpdet.h \
$(LIBF)/phyhanbo1/cpdet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a cpdet.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/cpdet.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a cpdet.o ; $(RM) cpdet.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(cvmgp.o) : $(LIBF)/phyhanbo1/cvmgp.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a cvmgp.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/cvmgp.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(cvmgt.o) : $(LIBF)/phyhanbo1/cvmgt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a cvmgt.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/cvmgt.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(datareadnc.o) : $(LIBF)/phyhanbo1/datareadnc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxyprim.h \
$(LIBO)/libphyhanbo1.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a datareadnc.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/datareadnc.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a datareadnc.o ; $(RM) datareadnc.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(dsolver.o) : $(LIBF)/phyhanbo1/dsolver.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a dsolver.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/dsolver.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a dsolver.o ; $(RM) dsolver.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(dtridgl.o) : $(LIBF)/phyhanbo1/dtridgl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a dtridgl.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/dtridgl.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a dtridgl.o ; $(RM) dtridgl.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(evap.o) : $(LIBF)/phyhanbo1/evap.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a evap.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/evap.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a evap.o ; $(RM) evap.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(forceWCfn.o) : $(LIBF)/phyhanbo1/forceWCfn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a forceWCfn.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/forceWCfn.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a forceWCfn.o ; $(RM) forceWCfn.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(gfluxi.o) : $(LIBF)/phyhanbo1/gfluxi.F \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a gfluxi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/gfluxi.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a gfluxi.o ; $(RM) gfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(gfluxv.o) : $(LIBF)/phyhanbo1/gfluxv.F \
$(LIBO)/libphyhanbo1.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a gfluxv.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/gfluxv.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a gfluxv.o ; $(RM) gfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iniaerosol.o) : $(LIBF)/phyhanbo1/iniaerosol.F \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iniaerosol.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/iniaerosol.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iniaerosol.o ; $(RM) iniaerosol.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(ini_archive.o) : $(LIBF)/phyhanbo1/ini_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo1.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a ini_archive.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/ini_archive.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a ini_archive.o ; $(RM) ini_archive.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(inifis.o) : $(LIBF)/phyhanbo1/inifis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o) \
$(LIBO)/libphyhanbo1.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo1.a(planete_mod.o) \
$(LIBO)/libphyhanbo1.a(planetwide_mod.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a inifis.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/inifis.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a inifis.o ; $(RM) inifis.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iniorbit.o) : $(LIBF)/phyhanbo1/iniorbit.F \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iniorbit.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/iniorbit.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iniorbit.o ; $(RM) iniorbit.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(inistats.o) : $(LIBF)/phyhanbo1/inistats.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phyhanbo1/statto.h \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a inistats.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/inistats.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a inistats.o ; $(RM) inistats.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(initracer.o) : $(LIBF)/phyhanbo1/initracer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a initracer.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/initracer.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a initracer.o ; $(RM) initracer.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iniwrite.o) : $(LIBF)/phyhanbo1/iniwrite.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo1.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iniwrite.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/iniwrite.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iniwrite.o ; $(RM) iniwrite.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iniwrite_specIR.o) : $(LIBF)/phyhanbo1/iniwrite_specIR.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iniwrite_specIR.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/iniwrite_specIR.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iniwrite_specIR.o ; $(RM) iniwrite_specIR.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iniwrite_specVI.o) : $(LIBF)/phyhanbo1/iniwrite_specVI.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iniwrite_specVI.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/iniwrite_specVI.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iniwrite_specVI.o ; $(RM) iniwrite_specVI.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(interp_line.o) : $(LIBF)/phyhanbo1/interp_line.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a interp_line.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/interp_line.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a interp_line.o ; $(RM) interp_line.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(lagrange.o) : $(LIBF)/phyhanbo1/lagrange.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a lagrange.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/lagrange.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a lagrange.o ; $(RM) lagrange.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(lect_start_archive.o) : $(LIBF)/phyhanbo1/lect_start_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a lect_start_archive.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/lect_start_archive.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a lect_start_archive.o ; $(RM) lect_start_archive.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mucorr.o) : $(LIBF)/phyhanbo1/mucorr.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mucorr.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/mucorr.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mucorr.o ; $(RM) mucorr.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mufract.o) : $(LIBF)/phyhanbo1/mufract.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mufract.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/mufract.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mufract.o ; $(RM) mufract.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(multipl.o) : $(LIBF)/phyhanbo1/multipl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a multipl.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/multipl.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(newsedim.o) : $(LIBF)/phyhanbo1/newsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a newsedim.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/newsedim.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a newsedim.o ; $(RM) newsedim.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(orbite.o) : $(LIBF)/phyhanbo1/orbite.F \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a orbite.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/orbite.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a orbite.o ; $(RM) orbite.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(profile.o) : $(LIBF)/phyhanbo1/profile.F \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a profile.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/profile.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a profile.o ; $(RM) profile.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(scopyi.o) : $(LIBF)/phyhanbo1/scopyi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a scopyi.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/scopyi.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(sfluxi.o) : $(LIBF)/phyhanbo1/sfluxi.F \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a sfluxi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/sfluxi.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a sfluxi.o ; $(RM) sfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(sfluxv.o) : $(LIBF)/phyhanbo1/sfluxv.F \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a sfluxv.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/sfluxv.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a sfluxv.o ; $(RM) sfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(soil.o) : $(LIBF)/phyhanbo1/soil.F \
$(LIBO)/libphyhanbo1.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a soil.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/soil.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a soil.o ; $(RM) soil.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(soil_settings.o) : $(LIBF)/phyhanbo1/soil_settings.F \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(iostart.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a soil_settings.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/soil_settings.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a soil_settings.o ; $(RM) soil_settings.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(stelang.o) : $(LIBF)/phyhanbo1/stelang.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a stelang.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/stelang.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a stelang.o ; $(RM) stelang.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(stellarlong.o) : $(LIBF)/phyhanbo1/stellarlong.F \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a stellarlong.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/stellarlong.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a stellarlong.o ; $(RM) stellarlong.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(surface_nature.o) : $(LIBF)/phyhanbo1/surface_nature.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a surface_nature.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/surface_nature.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a surface_nature.o ; $(RM) surface_nature.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(surfini.o) : $(LIBF)/phyhanbo1/surfini.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(planetwide_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a surfini.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/surfini.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a surfini.o ; $(RM) surfini.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(tabfi.o) : $(LIBF)/phyhanbo1/tabfi.F \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(iostart.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a tabfi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/tabfi.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a tabfi.o ; $(RM) tabfi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(tpindex.o) : $(LIBF)/phyhanbo1/tpindex.F \
$(LIBO)/libphyhanbo1.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a tpindex.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/tpindex.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a tpindex.o ; $(RM) tpindex.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(vdif_cd.o) : $(LIBF)/phyhanbo1/vdif_cd.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a vdif_cd.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/vdif_cd.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a vdif_cd.o ; $(RM) vdif_cd.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(vdifc.o) : $(LIBF)/phyhanbo1/vdifc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a vdifc.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/vdifc.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a vdifc.o ; $(RM) vdifc.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(vdif_kc.o) : $(LIBF)/phyhanbo1/vdif_kc.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a vdif_kc.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/vdif_kc.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a vdif_kc.o ; $(RM) vdif_kc.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(vlz_fi.o) : $(LIBF)/phyhanbo1/vlz_fi.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a vlz_fi.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/vlz_fi.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a vlz_fi.o ; $(RM) vlz_fi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(write_archive.o) : $(LIBF)/phyhanbo1/write_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBO)/libphyhanbo1.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a write_archive.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/write_archive.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a write_archive.o ; $(RM) write_archive.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(writediagfi.o) : $(LIBF)/phyhanbo1/writediagfi.F \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a writediagfi.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/writediagfi.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a writediagfi.o ; $(RM) writediagfi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(writediagspecIR.o) : $(LIBF)/phyhanbo1/writediagspecIR.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a writediagspecIR.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/writediagspecIR.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a writediagspecIR.o ; $(RM) writediagspecIR.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(writediagspecVI.o) : $(LIBF)/phyhanbo1/writediagspecVI.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a writediagspecVI.o ; \
	$(COMPILE90) $(LIBF)/phyhanbo1/writediagspecVI.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a writediagspecVI.o ; $(RM) writediagspecVI.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(writeg1d.o) : $(LIBF)/phyhanbo1/writeg1d.F \
$(LIBF)/phyhanbo1/comg1d.h \
$(LIBF)/phyhanbo1/comg1d.h \
$(LIBF)/phyhanbo1/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a writeg1d.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/writeg1d.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a writeg1d.o ; $(RM) writeg1d.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(zerophys.o) : $(LIBF)/phyhanbo1/zerophys.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a zerophys.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/zerophys.F ; \
	$(AR) r $(LIBO)/libphyhanbo1.a zerophys.o ; $(RM) zerophys.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(aeropacity.o) : $(LIBF)/phyhanbo1/aeropacity.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a aeropacity.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/aeropacity.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a aeropacity.o ; $(RM) aeropacity.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(aeroptproperties.o) : $(LIBF)/phyhanbo1/aeroptproperties.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a aeroptproperties.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/aeroptproperties.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a aeroptproperties.o ; $(RM) aeroptproperties.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(aerosol_mod.o) : $(LIBF)/phyhanbo1/aerosol_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a aerosol_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/aerosol_mod.F90 ; \
	mv $(MOD_LOC_DIR)/aerosol_mod.$(MOD_SUFFIX) $(LIBO)/aerosol_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a aerosol_mod.o ; $(RM) aerosol_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(ave_stelspec.o) : $(LIBF)/phyhanbo1/ave_stelspec.F90 \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a ave_stelspec.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/ave_stelspec.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a ave_stelspec.o ; $(RM) ave_stelspec.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(bilinearbig.o) : $(LIBF)/phyhanbo1/bilinearbig.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a bilinearbig.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/bilinearbig.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a bilinearbig.o ; $(RM) bilinearbig.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(bilinear.o) : $(LIBF)/phyhanbo1/bilinear.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a bilinear.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/bilinear.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a bilinear.o ; $(RM) bilinear.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(calc_cpp3d.o) : $(LIBF)/phyhanbo1/calc_cpp3d.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/phyhanbo1/cpdet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a calc_cpp3d.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/calc_cpp3d.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a calc_cpp3d.o ; $(RM) calc_cpp3d.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(calc_cpp_mugaz.o) : $(LIBF)/phyhanbo1/calc_cpp_mugaz.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a calc_cpp_mugaz.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/calc_cpp_mugaz.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a calc_cpp_mugaz.o ; $(RM) calc_cpp_mugaz.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(calcenergy_kcm.o) : $(LIBF)/phyhanbo1/calcenergy_kcm.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(params_h.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a calcenergy_kcm.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/calcenergy_kcm.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a calcenergy_kcm.o ; $(RM) calcenergy_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(calc_rayleigh.o) : $(LIBF)/phyhanbo1/calc_rayleigh.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a calc_rayleigh.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/calc_rayleigh.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a calc_rayleigh.o ; $(RM) calc_rayleigh.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(callcorrk.o) : $(LIBF)/phyhanbo1/callcorrk.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o) \
$(LIBO)/libphyhanbo1.a(radii_mod.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a callcorrk.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/callcorrk.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a callcorrk.o ; $(RM) callcorrk.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(comdiurn_h.o) : $(LIBF)/phyhanbo1/comdiurn_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a comdiurn_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/comdiurn_h.F90 ; \
	mv $(MOD_LOC_DIR)/comdiurn_h.$(MOD_SUFFIX) $(LIBO)/comdiurn_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a comdiurn_h.o ; $(RM) comdiurn_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) : $(LIBF)/phyhanbo1/comgeomfi_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a comgeomfi_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/comgeomfi_h.F90 ; \
	mv $(MOD_LOC_DIR)/comgeomfi_h.$(MOD_SUFFIX) $(LIBO)/comgeomfi_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a comgeomfi_h.o ; $(RM) comgeomfi_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(comgeomphy.o) : $(LIBF)/phyhanbo1/comgeomphy.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a comgeomphy.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/comgeomphy.F90 ; \
	mv $(MOD_LOC_DIR)/comgeomphy.$(MOD_SUFFIX) $(LIBO)/comgeomphy.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a comgeomphy.o ; $(RM) comgeomphy.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(comsaison_h.o) : $(LIBF)/phyhanbo1/comsaison_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a comsaison_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/comsaison_h.F90 ; \
	mv $(MOD_LOC_DIR)/comsaison_h.$(MOD_SUFFIX) $(LIBO)/comsaison_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a comsaison_h.o ; $(RM) comsaison_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(comsoil_h.o) : $(LIBF)/phyhanbo1/comsoil_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a comsoil_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/comsoil_h.F90 ; \
	mv $(MOD_LOC_DIR)/comsoil_h.$(MOD_SUFFIX) $(LIBO)/comsoil_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a comsoil_h.o ; $(RM) comsoil_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(condense_cloud.o) : $(LIBF)/phyhanbo1/condense_cloud.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o) \
$(LIBO)/libphyhanbo1.a(radii_mod.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a condense_cloud.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/condense_cloud.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a condense_cloud.o ; $(RM) condense_cloud.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(cp_neutral.o) : $(LIBF)/phyhanbo1/cp_neutral.F90 \
$(LIBO)/libphyhanbo1.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a cp_neutral.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/cp_neutral.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a cp_neutral.o ; $(RM) cp_neutral.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(datafile_mod.o) : $(LIBF)/phyhanbo1/datafile_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a datafile_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/datafile_mod.F90 ; \
	mv $(MOD_LOC_DIR)/datafile_mod.$(MOD_SUFFIX) $(LIBO)/datafile_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a datafile_mod.o ; $(RM) datafile_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(def_var.o) : $(LIBF)/phyhanbo1/def_var.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a def_var.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/def_var.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a def_var.o ; $(RM) def_var.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(dimphy.o) : $(LIBF)/phyhanbo1/dimphy.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a dimphy.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/dimphy.F90 ; \
	mv $(MOD_LOC_DIR)/dimphy.$(MOD_SUFFIX) $(LIBO)/dimphy.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a dimphy.o ; $(RM) dimphy.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(gases_h.o) : $(LIBF)/phyhanbo1/gases_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a gases_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/gases_h.F90 ; \
	mv $(MOD_LOC_DIR)/gases_h.$(MOD_SUFFIX) $(LIBO)/gases_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a gases_h.o ; $(RM) gases_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(gradients_kcm.o) : $(LIBF)/phyhanbo1/gradients_kcm.F90 \
$(LIBO)/libphyhanbo1.a(params_h.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a gradients_kcm.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/gradients_kcm.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a gradients_kcm.o ; $(RM) gradients_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(hydrol.o) : $(LIBF)/phyhanbo1/hydrol.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a hydrol.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/hydrol.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a hydrol.o ; $(RM) hydrol.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iniphysiq.o) : $(LIBF)/phyhanbo1/iniphysiq.F90 \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(comgeomphy.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo1.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iniphysiq.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/iniphysiq.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iniphysiq.o ; $(RM) iniphysiq.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(init_phys_lmdz.o) : $(LIBF)/phyhanbo1/init_phys_lmdz.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a init_phys_lmdz.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/init_phys_lmdz.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a init_phys_lmdz.o ; $(RM) init_phys_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iniwritesoil.o) : $(LIBF)/phyhanbo1/iniwritesoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libphyhanbo1.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iniwritesoil.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/iniwritesoil.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iniwritesoil.o ; $(RM) iniwritesoil.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(interpolateH2H2.o) : $(LIBF)/phyhanbo1/interpolateH2H2.F90 \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a interpolateH2H2.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/interpolateH2H2.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a interpolateH2H2.o ; $(RM) interpolateH2H2.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(interpolateH2He.o) : $(LIBF)/phyhanbo1/interpolateH2He.F90 \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a interpolateH2He.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/interpolateH2He.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a interpolateH2He.o ; $(RM) interpolateH2He.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(interpolateH2Ocont_CKD.o) : $(LIBF)/phyhanbo1/interpolateH2Ocont_CKD.F90 \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a interpolateH2Ocont_CKD.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/interpolateH2Ocont_CKD.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a interpolateH2Ocont_CKD.o ; $(RM) interpolateH2Ocont_CKD.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(interpolateH2Ocont_PPC.o) : $(LIBF)/phyhanbo1/interpolateH2Ocont_PPC.F90 \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a interpolateH2Ocont_PPC.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/interpolateH2Ocont_PPC.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a interpolateH2Ocont_PPC.o ; $(RM) interpolateH2Ocont_PPC.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(interpolateN2H2.o) : $(LIBF)/phyhanbo1/interpolateN2H2.F90 \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a interpolateN2H2.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/interpolateN2H2.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a interpolateN2H2.o ; $(RM) interpolateN2H2.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(interpolateN2N2.o) : $(LIBF)/phyhanbo1/interpolateN2N2.F90 \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a interpolateN2N2.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/interpolateN2N2.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a interpolateN2N2.o ; $(RM) interpolateN2N2.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(iostart.o) : $(LIBF)/phyhanbo1/iostart.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a iostart.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/iostart.F90 ; \
	mv $(MOD_LOC_DIR)/iostart.$(MOD_SUFFIX) $(LIBO)/iostart.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a iostart.o ; $(RM) iostart.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(kcmprof_fn.o) : $(LIBF)/phyhanbo1/kcmprof_fn.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(params_h.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a kcmprof_fn.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/kcmprof_fn.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a kcmprof_fn.o ; $(RM) kcmprof_fn.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(largescale.o) : $(LIBF)/phyhanbo1/largescale.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a largescale.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/largescale.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a largescale.o ; $(RM) largescale.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mass_redistribution.o) : $(LIBF)/phyhanbo1/mass_redistribution.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o) \
$(LIBO)/libphyhanbo1.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mass_redistribution.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mass_redistribution.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mass_redistribution.o ; $(RM) mass_redistribution.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mkstat.o) : $(LIBF)/phyhanbo1/mkstat.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/statto.h \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mkstat.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mkstat.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mkstat.o ; $(RM) mkstat.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) : $(LIBF)/phyhanbo1/mod_grid_phy_lmdz.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mod_grid_phy_lmdz.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mod_grid_phy_lmdz.F90 ; \
	mv $(MOD_LOC_DIR)/mod_grid_phy_lmdz.$(MOD_SUFFIX) $(LIBO)/mod_grid_phy_lmdz.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mod_grid_phy_lmdz.o ; $(RM) mod_grid_phy_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) : $(LIBF)/phyhanbo1/mod_phys_lmdz_mpi_data.F90 \
$(LIBO)/libdyn3d.a(mod_const_mpi.o) \
$(LIBO)/libdyn3d.a(mod_const_mpi.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mod_phys_lmdz_mpi_data.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mod_phys_lmdz_mpi_data.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_mpi_data.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_mpi_data.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mod_phys_lmdz_mpi_data.o ; $(RM) mod_phys_lmdz_mpi_data.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_transfert.o) : $(LIBF)/phyhanbo1/mod_phys_lmdz_mpi_transfert.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mod_phys_lmdz_mpi_transfert.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mod_phys_lmdz_mpi_transfert.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_mpi_transfert.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_mpi_transfert.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mod_phys_lmdz_mpi_transfert.o ; $(RM) mod_phys_lmdz_mpi_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) : $(LIBF)/phyhanbo1/mod_phys_lmdz_omp_data.F90 \
$(LIBO)/libphyhanbo1.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mod_phys_lmdz_omp_data.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mod_phys_lmdz_omp_data.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_omp_data.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_omp_data.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mod_phys_lmdz_omp_data.o ; $(RM) mod_phys_lmdz_omp_data.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_transfert.o) : $(LIBF)/phyhanbo1/mod_phys_lmdz_omp_transfert.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mod_phys_lmdz_omp_transfert.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mod_phys_lmdz_omp_transfert.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_omp_transfert.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_omp_transfert.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mod_phys_lmdz_omp_transfert.o ; $(RM) mod_phys_lmdz_omp_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) : $(LIBF)/phyhanbo1/mod_phys_lmdz_para.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_transfert_para.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mod_phys_lmdz_para.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mod_phys_lmdz_para.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_para.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_para.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mod_phys_lmdz_para.o ; $(RM) mod_phys_lmdz_para.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_transfert_para.o) : $(LIBF)/phyhanbo1/mod_phys_lmdz_transfert_para.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_transfert.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_omp_transfert.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a mod_phys_lmdz_transfert_para.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/mod_phys_lmdz_transfert_para.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_transfert_para.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_transfert_para.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a mod_phys_lmdz_transfert_para.o ; $(RM) mod_phys_lmdz_transfert_para.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(moistadj.o) : $(LIBF)/phyhanbo1/moistadj.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a moistadj.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/moistadj.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a moistadj.o ; $(RM) moistadj.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(newtrelax.o) : $(LIBF)/phyhanbo1/newtrelax.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a newtrelax.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/newtrelax.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a newtrelax.o ; $(RM) newtrelax.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(ocean_slab_mod.o) : $(LIBF)/phyhanbo1/ocean_slab_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(surf_heat_transp_mod.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a ocean_slab_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/ocean_slab_mod.F90 ; \
	mv $(MOD_LOC_DIR)/ocean_slab_mod.$(MOD_SUFFIX) $(LIBO)/ocean_slab_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a ocean_slab_mod.o ; $(RM) ocean_slab_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(optci.o) : $(LIBF)/phyhanbo1/optci.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a optci.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/optci.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a optci.o ; $(RM) optci.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(optcv.o) : $(LIBF)/phyhanbo1/optcv.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a optcv.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/optcv.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a optcv.o ; $(RM) optcv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(params_h.o) : $(LIBF)/phyhanbo1/params_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a params_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/params_h.F90 ; \
	mv $(MOD_LOC_DIR)/params_h.$(MOD_SUFFIX) $(LIBO)/params_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a params_h.o ; $(RM) params_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(phasechange_kcm.o) : $(LIBF)/phyhanbo1/phasechange_kcm.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a phasechange_kcm.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/phasechange_kcm.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a phasechange_kcm.o ; $(RM) phasechange_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(phyetat0.o) : $(LIBF)/phyhanbo1/phyetat0.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(iostart.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a phyetat0.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/phyetat0.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a phyetat0.o ; $(RM) phyetat0.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(phyredem.o) : $(LIBF)/phyhanbo1/phyredem.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(iostart.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(planete_mod.o) \
$(LIBO)/libphyhanbo1.a(iostart.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a phyredem.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/phyredem.F90 ; \
	mv $(MOD_LOC_DIR)/phyredem.$(MOD_SUFFIX) $(LIBO)/phyredem.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a phyredem.o ; $(RM) phyredem.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(physiq.o) : $(LIBF)/phyhanbo1/physiq.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radii_mod.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo1.a(comsaison_h.o) \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo1.a(phyredem.o) \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o) \
$(LIBO)/libphyhanbo1.a(ocean_slab_mod.o) \
$(LIBO)/libphyhanbo1.a(surf_heat_transp_mod.o) \
$(LIBO)/libphyhanbo1.a(planetwide_mod.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a physiq.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/physiq.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a physiq.o ; $(RM) physiq.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(planete_mod.o) : $(LIBF)/phyhanbo1/planete_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a planete_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/planete_mod.F90 ; \
	mv $(MOD_LOC_DIR)/planete_mod.$(MOD_SUFFIX) $(LIBO)/planete_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a planete_mod.o ; $(RM) planete_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(planetwide_mod.o) : $(LIBF)/phyhanbo1/planetwide_mod.F90 \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphyhanbo1.a(dimphy.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a planetwide_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/planetwide_mod.F90 ; \
	mv $(MOD_LOC_DIR)/planetwide_mod.$(MOD_SUFFIX) $(LIBO)/planetwide_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a planetwide_mod.o ; $(RM) planetwide_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(radcommon_h.o) : $(LIBF)/phyhanbo1/radcommon_h.F90 \
$(LIBO)/libphyhanbo1.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a radcommon_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/radcommon_h.F90 ; \
	mv $(MOD_LOC_DIR)/radcommon_h.$(MOD_SUFFIX) $(LIBO)/radcommon_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a radcommon_h.o ; $(RM) radcommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(radii_mod.o) : $(LIBF)/phyhanbo1/radii_mod.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a radii_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/radii_mod.F90 ; \
	mv $(MOD_LOC_DIR)/radii_mod.$(MOD_SUFFIX) $(LIBO)/radii_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a radii_mod.o ; $(RM) radii_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(radinc_h.o) : $(LIBF)/phyhanbo1/radinc_h.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo/bands.h \
$(LIBF)/phyhanbo/scatterers.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a radinc_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/radinc_h.F90 ; \
	mv $(MOD_LOC_DIR)/radinc_h.$(MOD_SUFFIX) $(LIBO)/radinc_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a radinc_h.o ; $(RM) radinc_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(rain.o) : $(LIBF)/phyhanbo1/rain.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(radii_mod.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a rain.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/rain.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a rain.o ; $(RM) rain.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(rings.o) : $(LIBF)/phyhanbo1/rings.F90 \
$(LIBO)/libphyhanbo1.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a rings.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/rings.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a rings.o ; $(RM) rings.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(setspi.o) : $(LIBF)/phyhanbo1/setspi.F90 \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a setspi.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/setspi.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a setspi.o ; $(RM) setspi.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(setspv.o) : $(LIBF)/phyhanbo1/setspv.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a setspv.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/setspv.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a setspv.o ; $(RM) setspv.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(slab_ice_h.o) : $(LIBF)/phyhanbo1/slab_ice_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a slab_ice_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/slab_ice_h.F90 ; \
	mv $(MOD_LOC_DIR)/slab_ice_h.$(MOD_SUFFIX) $(LIBO)/slab_ice_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a slab_ice_h.o ; $(RM) slab_ice_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(stokes.o) : $(LIBF)/phyhanbo1/stokes.F90 \
$(LIBF)/phyhanbo1/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a stokes.o ; \
	$(COMPILE) $(LIBF)/phyhanbo1/stokes.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a stokes.o ; $(RM) stokes.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(suaer_corrk.o) : $(LIBF)/phyhanbo1/suaer_corrk.F90 \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a suaer_corrk.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/suaer_corrk.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a suaer_corrk.o ; $(RM) suaer_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(sugas_corrk.o) : $(LIBF)/phyhanbo1/sugas_corrk.F90 \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBO)/libphyhanbo1.a(radinc_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(datafile_mod.o) \
$(LIBO)/libphyhanbo1.a(gases_h.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a sugas_corrk.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/sugas_corrk.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a sugas_corrk.o ; $(RM) sugas_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(su_gases.o) : $(LIBF)/phyhanbo1/su_gases.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a su_gases.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/su_gases.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a su_gases.o ; $(RM) su_gases.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(surfdat_h.o) : $(LIBF)/phyhanbo1/surfdat_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a surfdat_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/surfdat_h.F90 ; \
	mv $(MOD_LOC_DIR)/surfdat_h.$(MOD_SUFFIX) $(LIBO)/surfdat_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a surfdat_h.o ; $(RM) surfdat_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(surf_heat_transp_mod.o) : $(LIBF)/phyhanbo1/surf_heat_transp_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phyhanbo1/comhdiff.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phyhanbo1/comhdiff.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phyhanbo1/comhdiff.h \
$(LIBO)/libphyhanbo1.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a surf_heat_transp_mod.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/surf_heat_transp_mod.F90 ; \
	mv $(MOD_LOC_DIR)/surf_heat_transp_mod.$(MOD_SUFFIX) $(LIBO)/surf_heat_transp_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a surf_heat_transp_mod.o ; $(RM) surf_heat_transp_mod.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(su_watercycle.o) : $(LIBF)/phyhanbo1/su_watercycle.F90 \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a su_watercycle.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/su_watercycle.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a su_watercycle.o ; $(RM) su_watercycle.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(totalcloudfrac.o) : $(LIBF)/phyhanbo1/totalcloudfrac.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phyhanbo1/comcstfi.h \
$(LIBF)/phyhanbo1/callkeys.h \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(comdiurn_h.o) \
$(LIBO)/libphyhanbo1.a(comgeomfi_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a totalcloudfrac.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/totalcloudfrac.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a totalcloudfrac.o ; $(RM) totalcloudfrac.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(tracer_h.o) : $(LIBF)/phyhanbo1/tracer_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a tracer_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/tracer_h.F90 ; \
	mv $(MOD_LOC_DIR)/tracer_h.$(MOD_SUFFIX) $(LIBO)/tracer_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a tracer_h.o ; $(RM) tracer_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(turbdiff.o) : $(LIBF)/phyhanbo1/turbdiff.F90 \
$(LIBO)/libphyhanbo1.a(watercommon_h.o) \
$(LIBO)/libphyhanbo1.a(radcommon_h.o) \
$(LIBO)/libphyhanbo1.a(surfdat_h.o) \
$(LIBO)/libphyhanbo1.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a turbdiff.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/turbdiff.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a turbdiff.o ; $(RM) turbdiff.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(watercommon_h.o) : $(LIBF)/phyhanbo1/watercommon_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a watercommon_h.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/watercommon_h.F90 ; \
	mv $(MOD_LOC_DIR)/watercommon_h.$(MOD_SUFFIX) $(LIBO)/watercommon_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphyhanbo1.a watercommon_h.o ; $(RM) watercommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(watersat.o) : $(LIBF)/phyhanbo1/watersat.F90 \
$(LIBO)/libphyhanbo1.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a watersat.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/watersat.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a watersat.o ; $(RM) watersat.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(watersat_grad.o) : $(LIBF)/phyhanbo1/watersat_grad.F90 \
$(LIBO)/libphyhanbo1.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a watersat_grad.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/watersat_grad.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a watersat_grad.o ; $(RM) watersat_grad.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(writediagsoil.o) : $(LIBF)/phyhanbo1/writediagsoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libphyhanbo1.a(comsoil_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a writediagsoil.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/writediagsoil.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a writediagsoil.o ; $(RM) writediagsoil.o ; \
	cd $(GCM)

$(LIBO)/libphyhanbo1.a(wstats.o) : $(LIBF)/phyhanbo1/wstats.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phyhanbo1/statto.h \
$(LIBO)/libphyhanbo1.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphyhanbo1.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphyhanbo1.a wstats.o ; \
	$(COMPTRU90) $(LIBF)/phyhanbo1/wstats.F90 ; \
	$(AR) r $(LIBO)/libphyhanbo1.a wstats.o ; $(RM) wstats.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phystd
#=======================================================================

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aerave.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aerave_new.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(blackl.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(callsedim.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(convadj.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cpdet.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cvmgp.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cvmgt.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(datareadnc.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(dsolver.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(dtridgl.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(evap.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(forceWCfn.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(gfluxi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(gfluxv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniaerosol.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(ini_archive.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(inifis.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniorbit.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(inistats.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(initracer.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwrite.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwrite_specIR.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwrite_specVI.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interp_line.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(lagrange.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(lect_start_archive.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mucorr.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mufract.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(multipl.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(newsedim.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(orbite.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(profile.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(scopyi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(sfluxi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(sfluxv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(soil.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(soil_settings.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(stelang.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(stellarlong.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surface_nature.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surfini.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(tabfi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(tpindex.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vdif_cd.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vdifc.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vdif_kc.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vlz_fi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(write_archive.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagfi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagspecIR.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagspecVI.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writeg1d.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(zerophys.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aeropacity.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aeroptproperties.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aerosol_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(ave_stelspec.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(bilinearbig.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(bilinear.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(calc_cpp3d.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(calc_cpp_mugaz.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(calcenergy_kcm.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(calc_rayleigh.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(callcorrk.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comdiurn_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comgeomfi_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comgeomphy.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comsaison_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comsoil_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(condense_cloud.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cp_neutral.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(datafile_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(def_var.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(dimphy.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(gases_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(gradients_kcm.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(hydrol.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniphysiq.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(init_phys_lmdz.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwritesoil.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2H2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2He.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2Ocont_CKD.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2Ocont_PPC.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateN2H2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateN2N2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iostart.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(kcmprof_fn.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(largescale.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mass_redistribution.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mkstat.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mod_grid_phy_lmdz.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mod_phys_lmdz_mpi_transfert.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mod_phys_lmdz_omp_transfert.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mod_phys_lmdz_para.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mod_phys_lmdz_transfert_para.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(moistadj.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(newtrelax.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(ocean_slab_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(optci.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(optcv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(params_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(phasechange_kcm.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(phyetat0.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(phyredem.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(physiq.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(planete_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(planetwide_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(radcommon_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(radii_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(radinc_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(rain.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(rings.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(setspi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(setspv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(slab_ice_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(stokes.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(suaer_corrk.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(sugas_corrk.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(su_gases.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surfdat_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surf_heat_transp_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(su_watercycle.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(totalcloudfrac.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(tracer_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(turbdiff.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(watercommon_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(watersat.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(watersat_grad.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagsoil.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(wstats.o)

.PRECIOUS	: $(LIBO)/libphystd.a


# Compilation des membres de la bibliotheque libphystd.a

$(LIBO)/libphystd.a(aerave.o) : $(LIBF)/phystd/aerave.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a aerave.o ; \
	$(COMPILE) $(LIBF)/phystd/aerave.F ; \
	$(AR) r $(LIBO)/libphystd.a aerave.o ; $(RM) aerave.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aerave_new.o) : $(LIBF)/phystd/aerave_new.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a aerave_new.o ; \
	$(COMPILE) $(LIBF)/phystd/aerave_new.F ; \
	$(AR) r $(LIBO)/libphystd.a aerave_new.o ; $(RM) aerave_new.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(blackl.o) : $(LIBF)/phystd/blackl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a blackl.o ; \
	$(COMPILE) $(LIBF)/phystd/blackl.F ; \
	$(AR) r $(LIBO)/libphystd.a blackl.o ; $(RM) blackl.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(callsedim.o) : $(LIBF)/phystd/callsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a callsedim.o ; \
	$(COMPILE90) $(LIBF)/phystd/callsedim.F ; \
	$(AR) r $(LIBO)/libphystd.a callsedim.o ; $(RM) callsedim.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(convadj.o) : $(LIBF)/phystd/convadj.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a convadj.o ; \
	$(COMPILE90) $(LIBF)/phystd/convadj.F ; \
	$(AR) r $(LIBO)/libphystd.a convadj.o ; $(RM) convadj.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cpdet.o) : $(LIBF)/phystd/cpdet.F \
$(LIBF)/phystd/cpdet.h \
$(LIBF)/phystd/cpdet.h \
$(LIBF)/phystd/cpdet.h \
$(LIBF)/phystd/cpdet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a cpdet.o ; \
	$(COMPILE) $(LIBF)/phystd/cpdet.F ; \
	$(AR) r $(LIBO)/libphystd.a cpdet.o ; $(RM) cpdet.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cvmgp.o) : $(LIBF)/phystd/cvmgp.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a cvmgp.o ; \
	$(COMPILE) $(LIBF)/phystd/cvmgp.F ; \
	$(AR) r $(LIBO)/libphystd.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cvmgt.o) : $(LIBF)/phystd/cvmgt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a cvmgt.o ; \
	$(COMPILE) $(LIBF)/phystd/cvmgt.F ; \
	$(AR) r $(LIBO)/libphystd.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(datareadnc.o) : $(LIBF)/phystd/datareadnc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxyprim.h \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a datareadnc.o ; \
	$(COMPILE90) $(LIBF)/phystd/datareadnc.F ; \
	$(AR) r $(LIBO)/libphystd.a datareadnc.o ; $(RM) datareadnc.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(dsolver.o) : $(LIBF)/phystd/dsolver.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a dsolver.o ; \
	$(COMPILE) $(LIBF)/phystd/dsolver.F ; \
	$(AR) r $(LIBO)/libphystd.a dsolver.o ; $(RM) dsolver.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(dtridgl.o) : $(LIBF)/phystd/dtridgl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a dtridgl.o ; \
	$(COMPILE) $(LIBF)/phystd/dtridgl.F ; \
	$(AR) r $(LIBO)/libphystd.a dtridgl.o ; $(RM) dtridgl.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(evap.o) : $(LIBF)/phystd/evap.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a evap.o ; \
	$(COMPILE90) $(LIBF)/phystd/evap.F ; \
	$(AR) r $(LIBO)/libphystd.a evap.o ; $(RM) evap.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(forceWCfn.o) : $(LIBF)/phystd/forceWCfn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a forceWCfn.o ; \
	$(COMPILE90) $(LIBF)/phystd/forceWCfn.F ; \
	$(AR) r $(LIBO)/libphystd.a forceWCfn.o ; $(RM) forceWCfn.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(gfluxi.o) : $(LIBF)/phystd/gfluxi.F \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a gfluxi.o ; \
	$(COMPILE90) $(LIBF)/phystd/gfluxi.F ; \
	$(AR) r $(LIBO)/libphystd.a gfluxi.o ; $(RM) gfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(gfluxv.o) : $(LIBF)/phystd/gfluxv.F \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a gfluxv.o ; \
	$(COMPILE90) $(LIBF)/phystd/gfluxv.F ; \
	$(AR) r $(LIBO)/libphystd.a gfluxv.o ; $(RM) gfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniaerosol.o) : $(LIBF)/phystd/iniaerosol.F \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iniaerosol.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniaerosol.F ; \
	$(AR) r $(LIBO)/libphystd.a iniaerosol.o ; $(RM) iniaerosol.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(ini_archive.o) : $(LIBF)/phystd/ini_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphystd.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a ini_archive.o ; \
	$(COMPILE90) $(LIBF)/phystd/ini_archive.F ; \
	$(AR) r $(LIBO)/libphystd.a ini_archive.o ; $(RM) ini_archive.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(inifis.o) : $(LIBF)/phystd/inifis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(planetwide_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a inifis.o ; \
	$(COMPILE90) $(LIBF)/phystd/inifis.F ; \
	$(AR) r $(LIBO)/libphystd.a inifis.o ; $(RM) inifis.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniorbit.o) : $(LIBF)/phystd/iniorbit.F \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iniorbit.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniorbit.F ; \
	$(AR) r $(LIBO)/libphystd.a iniorbit.o ; $(RM) iniorbit.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(inistats.o) : $(LIBF)/phystd/inistats.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phystd/statto.h \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a inistats.o ; \
	$(COMPILE90) $(LIBF)/phystd/inistats.F ; \
	$(AR) r $(LIBO)/libphystd.a inistats.o ; $(RM) inistats.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(initracer.o) : $(LIBF)/phystd/initracer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a initracer.o ; \
	$(COMPILE90) $(LIBF)/phystd/initracer.F ; \
	$(AR) r $(LIBO)/libphystd.a initracer.o ; $(RM) initracer.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwrite.o) : $(LIBF)/phystd/iniwrite.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphystd.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iniwrite.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniwrite.F ; \
	$(AR) r $(LIBO)/libphystd.a iniwrite.o ; $(RM) iniwrite.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwrite_specIR.o) : $(LIBF)/phystd/iniwrite_specIR.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iniwrite_specIR.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniwrite_specIR.F ; \
	$(AR) r $(LIBO)/libphystd.a iniwrite_specIR.o ; $(RM) iniwrite_specIR.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwrite_specVI.o) : $(LIBF)/phystd/iniwrite_specVI.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iniwrite_specVI.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniwrite_specVI.F ; \
	$(AR) r $(LIBO)/libphystd.a iniwrite_specVI.o ; $(RM) iniwrite_specVI.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interp_line.o) : $(LIBF)/phystd/interp_line.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a interp_line.o ; \
	$(COMPILE) $(LIBF)/phystd/interp_line.F ; \
	$(AR) r $(LIBO)/libphystd.a interp_line.o ; $(RM) interp_line.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(lagrange.o) : $(LIBF)/phystd/lagrange.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a lagrange.o ; \
	$(COMPILE) $(LIBF)/phystd/lagrange.F ; \
	$(AR) r $(LIBO)/libphystd.a lagrange.o ; $(RM) lagrange.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(lect_start_archive.o) : $(LIBF)/phystd/lect_start_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a lect_start_archive.o ; \
	$(COMPILE90) $(LIBF)/phystd/lect_start_archive.F ; \
	$(AR) r $(LIBO)/libphystd.a lect_start_archive.o ; $(RM) lect_start_archive.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mucorr.o) : $(LIBF)/phystd/mucorr.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mucorr.o ; \
	$(COMPILE) $(LIBF)/phystd/mucorr.F ; \
	$(AR) r $(LIBO)/libphystd.a mucorr.o ; $(RM) mucorr.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mufract.o) : $(LIBF)/phystd/mufract.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mufract.o ; \
	$(COMPILE) $(LIBF)/phystd/mufract.F ; \
	$(AR) r $(LIBO)/libphystd.a mufract.o ; $(RM) mufract.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(multipl.o) : $(LIBF)/phystd/multipl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a multipl.o ; \
	$(COMPILE) $(LIBF)/phystd/multipl.F ; \
	$(AR) r $(LIBO)/libphystd.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(newsedim.o) : $(LIBF)/phystd/newsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a newsedim.o ; \
	$(COMPILE) $(LIBF)/phystd/newsedim.F ; \
	$(AR) r $(LIBO)/libphystd.a newsedim.o ; $(RM) newsedim.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(orbite.o) : $(LIBF)/phystd/orbite.F \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a orbite.o ; \
	$(COMPILE90) $(LIBF)/phystd/orbite.F ; \
	$(AR) r $(LIBO)/libphystd.a orbite.o ; $(RM) orbite.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(profile.o) : $(LIBF)/phystd/profile.F \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a profile.o ; \
	$(COMPILE90) $(LIBF)/phystd/profile.F ; \
	$(AR) r $(LIBO)/libphystd.a profile.o ; $(RM) profile.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(scopyi.o) : $(LIBF)/phystd/scopyi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a scopyi.o ; \
	$(COMPILE) $(LIBF)/phystd/scopyi.F ; \
	$(AR) r $(LIBO)/libphystd.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(sfluxi.o) : $(LIBF)/phystd/sfluxi.F \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a sfluxi.o ; \
	$(COMPILE90) $(LIBF)/phystd/sfluxi.F ; \
	$(AR) r $(LIBO)/libphystd.a sfluxi.o ; $(RM) sfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(sfluxv.o) : $(LIBF)/phystd/sfluxv.F \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a sfluxv.o ; \
	$(COMPILE90) $(LIBF)/phystd/sfluxv.F ; \
	$(AR) r $(LIBO)/libphystd.a sfluxv.o ; $(RM) sfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(soil.o) : $(LIBF)/phystd/soil.F \
$(LIBO)/libphystd.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a soil.o ; \
	$(COMPILE90) $(LIBF)/phystd/soil.F ; \
	$(AR) r $(LIBO)/libphystd.a soil.o ; $(RM) soil.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(soil_settings.o) : $(LIBF)/phystd/soil_settings.F \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(iostart.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a soil_settings.o ; \
	$(COMPILE90) $(LIBF)/phystd/soil_settings.F ; \
	$(AR) r $(LIBO)/libphystd.a soil_settings.o ; $(RM) soil_settings.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(stelang.o) : $(LIBF)/phystd/stelang.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a stelang.o ; \
	$(COMPILE) $(LIBF)/phystd/stelang.F ; \
	$(AR) r $(LIBO)/libphystd.a stelang.o ; $(RM) stelang.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(stellarlong.o) : $(LIBF)/phystd/stellarlong.F \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a stellarlong.o ; \
	$(COMPILE90) $(LIBF)/phystd/stellarlong.F ; \
	$(AR) r $(LIBO)/libphystd.a stellarlong.o ; $(RM) stellarlong.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surface_nature.o) : $(LIBF)/phystd/surface_nature.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a surface_nature.o ; \
	$(COMPILE90) $(LIBF)/phystd/surface_nature.F ; \
	$(AR) r $(LIBO)/libphystd.a surface_nature.o ; $(RM) surface_nature.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surfini.o) : $(LIBF)/phystd/surfini.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(planetwide_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a surfini.o ; \
	$(COMPILE90) $(LIBF)/phystd/surfini.F ; \
	$(AR) r $(LIBO)/libphystd.a surfini.o ; $(RM) surfini.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(tabfi.o) : $(LIBF)/phystd/tabfi.F \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(iostart.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a tabfi.o ; \
	$(COMPILE90) $(LIBF)/phystd/tabfi.F ; \
	$(AR) r $(LIBO)/libphystd.a tabfi.o ; $(RM) tabfi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(tpindex.o) : $(LIBF)/phystd/tpindex.F \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a tpindex.o ; \
	$(COMPILE90) $(LIBF)/phystd/tpindex.F ; \
	$(AR) r $(LIBO)/libphystd.a tpindex.o ; $(RM) tpindex.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vdif_cd.o) : $(LIBF)/phystd/vdif_cd.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a vdif_cd.o ; \
	$(COMPILE) $(LIBF)/phystd/vdif_cd.F ; \
	$(AR) r $(LIBO)/libphystd.a vdif_cd.o ; $(RM) vdif_cd.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vdifc.o) : $(LIBF)/phystd/vdifc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a vdifc.o ; \
	$(COMPILE90) $(LIBF)/phystd/vdifc.F ; \
	$(AR) r $(LIBO)/libphystd.a vdifc.o ; $(RM) vdifc.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vdif_kc.o) : $(LIBF)/phystd/vdif_kc.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a vdif_kc.o ; \
	$(COMPILE) $(LIBF)/phystd/vdif_kc.F ; \
	$(AR) r $(LIBO)/libphystd.a vdif_kc.o ; $(RM) vdif_kc.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vlz_fi.o) : $(LIBF)/phystd/vlz_fi.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a vlz_fi.o ; \
	$(COMPILE) $(LIBF)/phystd/vlz_fi.F ; \
	$(AR) r $(LIBO)/libphystd.a vlz_fi.o ; $(RM) vlz_fi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(write_archive.o) : $(LIBF)/phystd/write_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBO)/libphystd.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a write_archive.o ; \
	$(COMPILE90) $(LIBF)/phystd/write_archive.F ; \
	$(AR) r $(LIBO)/libphystd.a write_archive.o ; $(RM) write_archive.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagfi.o) : $(LIBF)/phystd/writediagfi.F \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a writediagfi.o ; \
	$(COMPILE90) $(LIBF)/phystd/writediagfi.F ; \
	$(AR) r $(LIBO)/libphystd.a writediagfi.o ; $(RM) writediagfi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagspecIR.o) : $(LIBF)/phystd/writediagspecIR.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a writediagspecIR.o ; \
	$(COMPILE90) $(LIBF)/phystd/writediagspecIR.F ; \
	$(AR) r $(LIBO)/libphystd.a writediagspecIR.o ; $(RM) writediagspecIR.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagspecVI.o) : $(LIBF)/phystd/writediagspecVI.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libdyn3d.a(control_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a writediagspecVI.o ; \
	$(COMPILE90) $(LIBF)/phystd/writediagspecVI.F ; \
	$(AR) r $(LIBO)/libphystd.a writediagspecVI.o ; $(RM) writediagspecVI.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writeg1d.o) : $(LIBF)/phystd/writeg1d.F \
$(LIBF)/phystd/comg1d.h \
$(LIBF)/phystd/comg1d.h \
$(LIBF)/phystd/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a writeg1d.o ; \
	$(COMPILE) $(LIBF)/phystd/writeg1d.F ; \
	$(AR) r $(LIBO)/libphystd.a writeg1d.o ; $(RM) writeg1d.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(zerophys.o) : $(LIBF)/phystd/zerophys.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a zerophys.o ; \
	$(COMPILE) $(LIBF)/phystd/zerophys.F ; \
	$(AR) r $(LIBO)/libphystd.a zerophys.o ; $(RM) zerophys.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aeropacity.o) : $(LIBF)/phystd/aeropacity.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a aeropacity.o ; \
	$(COMPTRU90) $(LIBF)/phystd/aeropacity.F90 ; \
	$(AR) r $(LIBO)/libphystd.a aeropacity.o ; $(RM) aeropacity.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aeroptproperties.o) : $(LIBF)/phystd/aeroptproperties.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a aeroptproperties.o ; \
	$(COMPTRU90) $(LIBF)/phystd/aeroptproperties.F90 ; \
	$(AR) r $(LIBO)/libphystd.a aeroptproperties.o ; $(RM) aeroptproperties.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aerosol_mod.o) : $(LIBF)/phystd/aerosol_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a aerosol_mod.o ; \
	$(COMPTRU90) $(LIBF)/phystd/aerosol_mod.F90 ; \
	mv $(MOD_LOC_DIR)/aerosol_mod.$(MOD_SUFFIX) $(LIBO)/aerosol_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a aerosol_mod.o ; $(RM) aerosol_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(ave_stelspec.o) : $(LIBF)/phystd/ave_stelspec.F90 \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a ave_stelspec.o ; \
	$(COMPTRU90) $(LIBF)/phystd/ave_stelspec.F90 ; \
	$(AR) r $(LIBO)/libphystd.a ave_stelspec.o ; $(RM) ave_stelspec.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(bilinearbig.o) : $(LIBF)/phystd/bilinearbig.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a bilinearbig.o ; \
	$(COMPILE) $(LIBF)/phystd/bilinearbig.F90 ; \
	$(AR) r $(LIBO)/libphystd.a bilinearbig.o ; $(RM) bilinearbig.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(bilinear.o) : $(LIBF)/phystd/bilinear.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a bilinear.o ; \
	$(COMPILE) $(LIBF)/phystd/bilinear.F90 ; \
	$(AR) r $(LIBO)/libphystd.a bilinear.o ; $(RM) bilinear.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(calc_cpp3d.o) : $(LIBF)/phystd/calc_cpp3d.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/phystd/cpdet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a calc_cpp3d.o ; \
	$(COMPILE) $(LIBF)/phystd/calc_cpp3d.F90 ; \
	$(AR) r $(LIBO)/libphystd.a calc_cpp3d.o ; $(RM) calc_cpp3d.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(calc_cpp_mugaz.o) : $(LIBF)/phystd/calc_cpp_mugaz.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a calc_cpp_mugaz.o ; \
	$(COMPTRU90) $(LIBF)/phystd/calc_cpp_mugaz.F90 ; \
	$(AR) r $(LIBO)/libphystd.a calc_cpp_mugaz.o ; $(RM) calc_cpp_mugaz.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(calcenergy_kcm.o) : $(LIBF)/phystd/calcenergy_kcm.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(params_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a calcenergy_kcm.o ; \
	$(COMPTRU90) $(LIBF)/phystd/calcenergy_kcm.F90 ; \
	$(AR) r $(LIBO)/libphystd.a calcenergy_kcm.o ; $(RM) calcenergy_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(calc_rayleigh.o) : $(LIBF)/phystd/calc_rayleigh.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a calc_rayleigh.o ; \
	$(COMPTRU90) $(LIBF)/phystd/calc_rayleigh.F90 ; \
	$(AR) r $(LIBO)/libphystd.a calc_rayleigh.o ; $(RM) calc_rayleigh.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(callcorrk.o) : $(LIBF)/phystd/callcorrk.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a callcorrk.o ; \
	$(COMPTRU90) $(LIBF)/phystd/callcorrk.F90 ; \
	$(AR) r $(LIBO)/libphystd.a callcorrk.o ; $(RM) callcorrk.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comdiurn_h.o) : $(LIBF)/phystd/comdiurn_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a comdiurn_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/comdiurn_h.F90 ; \
	mv $(MOD_LOC_DIR)/comdiurn_h.$(MOD_SUFFIX) $(LIBO)/comdiurn_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a comdiurn_h.o ; $(RM) comdiurn_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comgeomfi_h.o) : $(LIBF)/phystd/comgeomfi_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a comgeomfi_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/comgeomfi_h.F90 ; \
	mv $(MOD_LOC_DIR)/comgeomfi_h.$(MOD_SUFFIX) $(LIBO)/comgeomfi_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a comgeomfi_h.o ; $(RM) comgeomfi_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comgeomphy.o) : $(LIBF)/phystd/comgeomphy.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a comgeomphy.o ; \
	$(COMPTRU90) $(LIBF)/phystd/comgeomphy.F90 ; \
	mv $(MOD_LOC_DIR)/comgeomphy.$(MOD_SUFFIX) $(LIBO)/comgeomphy.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a comgeomphy.o ; $(RM) comgeomphy.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comsaison_h.o) : $(LIBF)/phystd/comsaison_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a comsaison_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/comsaison_h.F90 ; \
	mv $(MOD_LOC_DIR)/comsaison_h.$(MOD_SUFFIX) $(LIBO)/comsaison_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a comsaison_h.o ; $(RM) comsaison_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comsoil_h.o) : $(LIBF)/phystd/comsoil_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a comsoil_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/comsoil_h.F90 ; \
	mv $(MOD_LOC_DIR)/comsoil_h.$(MOD_SUFFIX) $(LIBO)/comsoil_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a comsoil_h.o ; $(RM) comsoil_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(condense_cloud.o) : $(LIBF)/phystd/condense_cloud.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a condense_cloud.o ; \
	$(COMPTRU90) $(LIBF)/phystd/condense_cloud.F90 ; \
	$(AR) r $(LIBO)/libphystd.a condense_cloud.o ; $(RM) condense_cloud.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cp_neutral.o) : $(LIBF)/phystd/cp_neutral.F90 \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a cp_neutral.o ; \
	$(COMPTRU90) $(LIBF)/phystd/cp_neutral.F90 ; \
	$(AR) r $(LIBO)/libphystd.a cp_neutral.o ; $(RM) cp_neutral.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(datafile_mod.o) : $(LIBF)/phystd/datafile_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a datafile_mod.o ; \
	$(COMPTRU90) $(LIBF)/phystd/datafile_mod.F90 ; \
	mv $(MOD_LOC_DIR)/datafile_mod.$(MOD_SUFFIX) $(LIBO)/datafile_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a datafile_mod.o ; $(RM) datafile_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(def_var.o) : $(LIBF)/phystd/def_var.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a def_var.o ; \
	$(COMPILE) $(LIBF)/phystd/def_var.F90 ; \
	$(AR) r $(LIBO)/libphystd.a def_var.o ; $(RM) def_var.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(dimphy.o) : $(LIBF)/phystd/dimphy.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a dimphy.o ; \
	$(COMPTRU90) $(LIBF)/phystd/dimphy.F90 ; \
	mv $(MOD_LOC_DIR)/dimphy.$(MOD_SUFFIX) $(LIBO)/dimphy.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a dimphy.o ; $(RM) dimphy.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(gases_h.o) : $(LIBF)/phystd/gases_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a gases_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/gases_h.F90 ; \
	mv $(MOD_LOC_DIR)/gases_h.$(MOD_SUFFIX) $(LIBO)/gases_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a gases_h.o ; $(RM) gases_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(gradients_kcm.o) : $(LIBF)/phystd/gradients_kcm.F90 \
$(LIBO)/libphystd.a(params_h.o) \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a gradients_kcm.o ; \
	$(COMPTRU90) $(LIBF)/phystd/gradients_kcm.F90 ; \
	$(AR) r $(LIBO)/libphystd.a gradients_kcm.o ; $(RM) gradients_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(hydrol.o) : $(LIBF)/phystd/hydrol.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a hydrol.o ; \
	$(COMPTRU90) $(LIBF)/phystd/hydrol.F90 ; \
	$(AR) r $(LIBO)/libphystd.a hydrol.o ; $(RM) hydrol.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniphysiq.o) : $(LIBF)/phystd/iniphysiq.F90 \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(comgeomphy.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iniphysiq.o ; \
	$(COMPTRU90) $(LIBF)/phystd/iniphysiq.F90 ; \
	$(AR) r $(LIBO)/libphystd.a iniphysiq.o ; $(RM) iniphysiq.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(init_phys_lmdz.o) : $(LIBF)/phystd/init_phys_lmdz.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a init_phys_lmdz.o ; \
	$(COMPTRU90) $(LIBF)/phystd/init_phys_lmdz.F90 ; \
	$(AR) r $(LIBO)/libphystd.a init_phys_lmdz.o ; $(RM) init_phys_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwritesoil.o) : $(LIBF)/phystd/iniwritesoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libphystd.a(comsoil_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iniwritesoil.o ; \
	$(COMPTRU90) $(LIBF)/phystd/iniwritesoil.F90 ; \
	$(AR) r $(LIBO)/libphystd.a iniwritesoil.o ; $(RM) iniwritesoil.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2H2.o) : $(LIBF)/phystd/interpolateH2H2.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2H2.o ; \
	$(COMPTRU90) $(LIBF)/phystd/interpolateH2H2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2H2.o ; $(RM) interpolateH2H2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2He.o) : $(LIBF)/phystd/interpolateH2He.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2He.o ; \
	$(COMPTRU90) $(LIBF)/phystd/interpolateH2He.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2He.o ; $(RM) interpolateH2He.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2Ocont_CKD.o) : $(LIBF)/phystd/interpolateH2Ocont_CKD.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2Ocont_CKD.o ; \
	$(COMPTRU90) $(LIBF)/phystd/interpolateH2Ocont_CKD.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2Ocont_CKD.o ; $(RM) interpolateH2Ocont_CKD.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2Ocont_PPC.o) : $(LIBF)/phystd/interpolateH2Ocont_PPC.F90 \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2Ocont_PPC.o ; \
	$(COMPTRU90) $(LIBF)/phystd/interpolateH2Ocont_PPC.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2Ocont_PPC.o ; $(RM) interpolateH2Ocont_PPC.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateN2H2.o) : $(LIBF)/phystd/interpolateN2H2.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a interpolateN2H2.o ; \
	$(COMPTRU90) $(LIBF)/phystd/interpolateN2H2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateN2H2.o ; $(RM) interpolateN2H2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateN2N2.o) : $(LIBF)/phystd/interpolateN2N2.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a interpolateN2N2.o ; \
	$(COMPTRU90) $(LIBF)/phystd/interpolateN2N2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateN2N2.o ; $(RM) interpolateN2N2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iostart.o) : $(LIBF)/phystd/iostart.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a iostart.o ; \
	$(COMPTRU90) $(LIBF)/phystd/iostart.F90 ; \
	mv $(MOD_LOC_DIR)/iostart.$(MOD_SUFFIX) $(LIBO)/iostart.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a iostart.o ; $(RM) iostart.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(kcmprof_fn.o) : $(LIBF)/phystd/kcmprof_fn.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(params_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a kcmprof_fn.o ; \
	$(COMPTRU90) $(LIBF)/phystd/kcmprof_fn.F90 ; \
	$(AR) r $(LIBO)/libphystd.a kcmprof_fn.o ; $(RM) kcmprof_fn.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(largescale.o) : $(LIBF)/phystd/largescale.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a largescale.o ; \
	$(COMPTRU90) $(LIBF)/phystd/largescale.F90 ; \
	$(AR) r $(LIBO)/libphystd.a largescale.o ; $(RM) largescale.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mass_redistribution.o) : $(LIBF)/phystd/mass_redistribution.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mass_redistribution.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mass_redistribution.F90 ; \
	$(AR) r $(LIBO)/libphystd.a mass_redistribution.o ; $(RM) mass_redistribution.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mkstat.o) : $(LIBF)/phystd/mkstat.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/statto.h \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mkstat.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mkstat.F90 ; \
	$(AR) r $(LIBO)/libphystd.a mkstat.o ; $(RM) mkstat.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) : $(LIBF)/phystd/mod_grid_phy_lmdz.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mod_grid_phy_lmdz.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mod_grid_phy_lmdz.F90 ; \
	mv $(MOD_LOC_DIR)/mod_grid_phy_lmdz.$(MOD_SUFFIX) $(LIBO)/mod_grid_phy_lmdz.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a mod_grid_phy_lmdz.o ; $(RM) mod_grid_phy_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) : $(LIBF)/phystd/mod_phys_lmdz_mpi_data.F90 \
$(LIBO)/libdyn3d.a(mod_const_mpi.o) \
$(LIBO)/libdyn3d.a(mod_const_mpi.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mod_phys_lmdz_mpi_data.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mod_phys_lmdz_mpi_data.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_mpi_data.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_mpi_data.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a mod_phys_lmdz_mpi_data.o ; $(RM) mod_phys_lmdz_mpi_data.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_transfert.o) : $(LIBF)/phystd/mod_phys_lmdz_mpi_transfert.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mod_phys_lmdz_mpi_transfert.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mod_phys_lmdz_mpi_transfert.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_mpi_transfert.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_mpi_transfert.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a mod_phys_lmdz_mpi_transfert.o ; $(RM) mod_phys_lmdz_mpi_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) : $(LIBF)/phystd/mod_phys_lmdz_omp_data.F90 \
$(LIBO)/libphystd.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mod_phys_lmdz_omp_data.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mod_phys_lmdz_omp_data.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_omp_data.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_omp_data.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a mod_phys_lmdz_omp_data.o ; $(RM) mod_phys_lmdz_omp_data.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mod_phys_lmdz_omp_transfert.o) : $(LIBF)/phystd/mod_phys_lmdz_omp_transfert.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mod_phys_lmdz_omp_transfert.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mod_phys_lmdz_omp_transfert.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_omp_transfert.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_omp_transfert.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a mod_phys_lmdz_omp_transfert.o ; $(RM) mod_phys_lmdz_omp_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) : $(LIBF)/phystd/mod_phys_lmdz_para.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_transfert_para.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mod_phys_lmdz_para.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mod_phys_lmdz_para.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_para.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_para.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a mod_phys_lmdz_para.o ; $(RM) mod_phys_lmdz_para.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mod_phys_lmdz_transfert_para.o) : $(LIBF)/phystd/mod_phys_lmdz_transfert_para.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_transfert.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_omp_transfert.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a mod_phys_lmdz_transfert_para.o ; \
	$(COMPTRU90) $(LIBF)/phystd/mod_phys_lmdz_transfert_para.F90 ; \
	mv $(MOD_LOC_DIR)/mod_phys_lmdz_transfert_para.$(MOD_SUFFIX) $(LIBO)/mod_phys_lmdz_transfert_para.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a mod_phys_lmdz_transfert_para.o ; $(RM) mod_phys_lmdz_transfert_para.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(moistadj.o) : $(LIBF)/phystd/moistadj.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a moistadj.o ; \
	$(COMPTRU90) $(LIBF)/phystd/moistadj.F90 ; \
	$(AR) r $(LIBO)/libphystd.a moistadj.o ; $(RM) moistadj.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(newtrelax.o) : $(LIBF)/phystd/newtrelax.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a newtrelax.o ; \
	$(COMPILE) $(LIBF)/phystd/newtrelax.F90 ; \
	$(AR) r $(LIBO)/libphystd.a newtrelax.o ; $(RM) newtrelax.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(ocean_slab_mod.o) : $(LIBF)/phystd/ocean_slab_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(surf_heat_transp_mod.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a ocean_slab_mod.o ; \
	$(COMPTRU90) $(LIBF)/phystd/ocean_slab_mod.F90 ; \
	mv $(MOD_LOC_DIR)/ocean_slab_mod.$(MOD_SUFFIX) $(LIBO)/ocean_slab_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a ocean_slab_mod.o ; $(RM) ocean_slab_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(optci.o) : $(LIBF)/phystd/optci.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a optci.o ; \
	$(COMPTRU90) $(LIBF)/phystd/optci.F90 ; \
	$(AR) r $(LIBO)/libphystd.a optci.o ; $(RM) optci.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(optcv.o) : $(LIBF)/phystd/optcv.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a optcv.o ; \
	$(COMPTRU90) $(LIBF)/phystd/optcv.F90 ; \
	$(AR) r $(LIBO)/libphystd.a optcv.o ; $(RM) optcv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(params_h.o) : $(LIBF)/phystd/params_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a params_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/params_h.F90 ; \
	mv $(MOD_LOC_DIR)/params_h.$(MOD_SUFFIX) $(LIBO)/params_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a params_h.o ; $(RM) params_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(phasechange_kcm.o) : $(LIBF)/phystd/phasechange_kcm.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a phasechange_kcm.o ; \
	$(COMPILE) $(LIBF)/phystd/phasechange_kcm.F90 ; \
	$(AR) r $(LIBO)/libphystd.a phasechange_kcm.o ; $(RM) phasechange_kcm.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(phyetat0.o) : $(LIBF)/phystd/phyetat0.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(iostart.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a phyetat0.o ; \
	$(COMPTRU90) $(LIBF)/phystd/phyetat0.F90 ; \
	$(AR) r $(LIBO)/libphystd.a phyetat0.o ; $(RM) phyetat0.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(phyredem.o) : $(LIBF)/phystd/phyredem.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(iostart.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(iostart.o) \
$(LIBO)/libdyn3d.a(infotrac.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a phyredem.o ; \
	$(COMPTRU90) $(LIBF)/phystd/phyredem.F90 ; \
	mv $(MOD_LOC_DIR)/phyredem.$(MOD_SUFFIX) $(LIBO)/phyredem.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a phyredem.o ; $(RM) phyredem.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(physiq.o) : $(LIBF)/phystd/physiq.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comsaison_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphystd.a(phyredem.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(ocean_slab_mod.o) \
$(LIBO)/libphystd.a(surf_heat_transp_mod.o) \
$(LIBO)/libphystd.a(planetwide_mod.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a physiq.o ; \
	$(COMPTRU90) $(LIBF)/phystd/physiq.F90 ; \
	$(AR) r $(LIBO)/libphystd.a physiq.o ; $(RM) physiq.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(planete_mod.o) : $(LIBF)/phystd/planete_mod.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a planete_mod.o ; \
	$(COMPTRU90) $(LIBF)/phystd/planete_mod.F90 ; \
	mv $(MOD_LOC_DIR)/planete_mod.$(MOD_SUFFIX) $(LIBO)/planete_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a planete_mod.o ; $(RM) planete_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(planetwide_mod.o) : $(LIBF)/phystd/planetwide_mod.F90 \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a planetwide_mod.o ; \
	$(COMPTRU90) $(LIBF)/phystd/planetwide_mod.F90 ; \
	mv $(MOD_LOC_DIR)/planetwide_mod.$(MOD_SUFFIX) $(LIBO)/planetwide_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a planetwide_mod.o ; $(RM) planetwide_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(radcommon_h.o) : $(LIBF)/phystd/radcommon_h.F90 \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a radcommon_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/radcommon_h.F90 ; \
	mv $(MOD_LOC_DIR)/radcommon_h.$(MOD_SUFFIX) $(LIBO)/radcommon_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a radcommon_h.o ; $(RM) radcommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(radii_mod.o) : $(LIBF)/phystd/radii_mod.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a radii_mod.o ; \
	$(COMPTRU90) $(LIBF)/phystd/radii_mod.F90 ; \
	mv $(MOD_LOC_DIR)/radii_mod.$(MOD_SUFFIX) $(LIBO)/radii_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a radii_mod.o ; $(RM) radii_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(radinc_h.o) : $(LIBF)/phystd/radinc_h.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/bands.h \
$(LIBF)/phystd/scatterers.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a radinc_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/radinc_h.F90 ; \
	mv $(MOD_LOC_DIR)/radinc_h.$(MOD_SUFFIX) $(LIBO)/radinc_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a radinc_h.o ; $(RM) radinc_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(rain.o) : $(LIBF)/phystd/rain.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a rain.o ; \
	$(COMPTRU90) $(LIBF)/phystd/rain.F90 ; \
	$(AR) r $(LIBO)/libphystd.a rain.o ; $(RM) rain.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(rings.o) : $(LIBF)/phystd/rings.F90 \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a rings.o ; \
	$(COMPTRU90) $(LIBF)/phystd/rings.F90 ; \
	$(AR) r $(LIBO)/libphystd.a rings.o ; $(RM) rings.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(setspi.o) : $(LIBF)/phystd/setspi.F90 \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a setspi.o ; \
	$(COMPTRU90) $(LIBF)/phystd/setspi.F90 ; \
	$(AR) r $(LIBO)/libphystd.a setspi.o ; $(RM) setspi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(setspv.o) : $(LIBF)/phystd/setspv.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a setspv.o ; \
	$(COMPTRU90) $(LIBF)/phystd/setspv.F90 ; \
	$(AR) r $(LIBO)/libphystd.a setspv.o ; $(RM) setspv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(slab_ice_h.o) : $(LIBF)/phystd/slab_ice_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a slab_ice_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/slab_ice_h.F90 ; \
	mv $(MOD_LOC_DIR)/slab_ice_h.$(MOD_SUFFIX) $(LIBO)/slab_ice_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a slab_ice_h.o ; $(RM) slab_ice_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(stokes.o) : $(LIBF)/phystd/stokes.F90 \
$(LIBF)/phystd/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a stokes.o ; \
	$(COMPILE) $(LIBF)/phystd/stokes.F90 ; \
	$(AR) r $(LIBO)/libphystd.a stokes.o ; $(RM) stokes.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(suaer_corrk.o) : $(LIBF)/phystd/suaer_corrk.F90 \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(aerosol_mod.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a suaer_corrk.o ; \
	$(COMPTRU90) $(LIBF)/phystd/suaer_corrk.F90 ; \
	$(AR) r $(LIBO)/libphystd.a suaer_corrk.o ; $(RM) suaer_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(sugas_corrk.o) : $(LIBF)/phystd/sugas_corrk.F90 \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libbibio.a(ioipsl_getincom_p.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a sugas_corrk.o ; \
	$(COMPTRU90) $(LIBF)/phystd/sugas_corrk.F90 ; \
	$(AR) r $(LIBO)/libphystd.a sugas_corrk.o ; $(RM) sugas_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(su_gases.o) : $(LIBF)/phystd/su_gases.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a su_gases.o ; \
	$(COMPTRU90) $(LIBF)/phystd/su_gases.F90 ; \
	$(AR) r $(LIBO)/libphystd.a su_gases.o ; $(RM) su_gases.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surfdat_h.o) : $(LIBF)/phystd/surfdat_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a surfdat_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/surfdat_h.F90 ; \
	mv $(MOD_LOC_DIR)/surfdat_h.$(MOD_SUFFIX) $(LIBO)/surfdat_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a surfdat_h.o ; $(RM) surfdat_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surf_heat_transp_mod.o) : $(LIBF)/phystd/surf_heat_transp_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phystd/comhdiff.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phystd/comhdiff.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phystd/comhdiff.h \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a surf_heat_transp_mod.o ; \
	$(COMPTRU90) $(LIBF)/phystd/surf_heat_transp_mod.F90 ; \
	mv $(MOD_LOC_DIR)/surf_heat_transp_mod.$(MOD_SUFFIX) $(LIBO)/surf_heat_transp_mod.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a surf_heat_transp_mod.o ; $(RM) surf_heat_transp_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(su_watercycle.o) : $(LIBF)/phystd/su_watercycle.F90 \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a su_watercycle.o ; \
	$(COMPTRU90) $(LIBF)/phystd/su_watercycle.F90 ; \
	$(AR) r $(LIBO)/libphystd.a su_watercycle.o ; $(RM) su_watercycle.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(totalcloudfrac.o) : $(LIBF)/phystd/totalcloudfrac.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phystd/comcstfi.h \
$(LIBF)/phystd/callkeys.h \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a totalcloudfrac.o ; \
	$(COMPTRU90) $(LIBF)/phystd/totalcloudfrac.F90 ; \
	$(AR) r $(LIBO)/libphystd.a totalcloudfrac.o ; $(RM) totalcloudfrac.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(tracer_h.o) : $(LIBF)/phystd/tracer_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a tracer_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/tracer_h.F90 ; \
	mv $(MOD_LOC_DIR)/tracer_h.$(MOD_SUFFIX) $(LIBO)/tracer_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a tracer_h.o ; $(RM) tracer_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(turbdiff.o) : $(LIBF)/phystd/turbdiff.F90 \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a turbdiff.o ; \
	$(COMPTRU90) $(LIBF)/phystd/turbdiff.F90 ; \
	$(AR) r $(LIBO)/libphystd.a turbdiff.o ; $(RM) turbdiff.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(watercommon_h.o) : $(LIBF)/phystd/watercommon_h.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a watercommon_h.o ; \
	$(COMPTRU90) $(LIBF)/phystd/watercommon_h.F90 ; \
	mv $(MOD_LOC_DIR)/watercommon_h.$(MOD_SUFFIX) $(LIBO)/watercommon_h.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphystd.a watercommon_h.o ; $(RM) watercommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(watersat.o) : $(LIBF)/phystd/watersat.F90 \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a watersat.o ; \
	$(COMPTRU90) $(LIBF)/phystd/watersat.F90 ; \
	$(AR) r $(LIBO)/libphystd.a watersat.o ; $(RM) watersat.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(watersat_grad.o) : $(LIBF)/phystd/watersat_grad.F90 \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a watersat_grad.o ; \
	$(COMPTRU90) $(LIBF)/phystd/watersat_grad.F90 ; \
	$(AR) r $(LIBO)/libphystd.a watersat_grad.o ; $(RM) watersat_grad.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagsoil.o) : $(LIBF)/phystd/writediagsoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libdyn3d.a(control_mod.o) \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a writediagsoil.o ; \
	$(COMPTRU90) $(LIBF)/phystd/writediagsoil.F90 ; \
	$(AR) r $(LIBO)/libphystd.a writediagsoil.o ; $(RM) writediagsoil.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(wstats.o) : $(LIBF)/phystd/wstats.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phystd/statto.h \
$(LIBO)/libphystd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphystd.a wstats.o ; \
	$(COMPTRU90) $(LIBF)/phystd/wstats.F90 ; \
	$(AR) r $(LIBO)/libphystd.a wstats.o ; $(RM) wstats.o ; \
	cd $(GCM)


