!
! $Id$
!
SUBROUTINE init_phys_lmdz(iim,jjp1,llm,nb_proc,distrib)
  USE mod_phys_lmdz_omp_data, only: klon_omp
  USE mod_grid_phy_lmdz, only: nbp_lev
  USE mod_phys_lmdz_para, only: init_phys_lmdz_para
  USE mod_grid_phy_lmdz, only: init_grid_phy_lmdz
  USE dimphy, ONLY : init_dimphy

  IMPLICIT NONE
  
    INTEGER,INTENT(in) :: iim
    INTEGER,INTENT(in) :: jjp1
    INTEGER,INTENT(in) :: llm
    INTEGER,INTENT(in) :: nb_proc
    INTEGER,INTENT(in) :: distrib(0:nb_proc-1)


    CALL init_grid_phy_lmdz(iim,jjp1,llm)
    CALL init_phys_lmdz_para(iim,jjp1,nb_proc,distrib)
!$OMP PARALLEL
    CALL init_dimphy(klon_omp,nbp_lev)

!$OMP END PARALLEL
 
END SUBROUTINE init_phys_lmdz  
