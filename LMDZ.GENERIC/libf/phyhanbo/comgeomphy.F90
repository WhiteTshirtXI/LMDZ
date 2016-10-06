module comgeomphy
   real,save,allocatable :: airephy(:)
   real,save,allocatable :: cuphy(:)
   real,save,allocatable :: cvphy(:)
   real,save,allocatable :: rlatd(:)
   real,save,allocatable :: rlond(:)
!$OMP THREADPRIVATE(airephy,cuphy,cvphy,rlatd,rlond)
contains
  
  subroutine initcomgeomphy
  USE mod_phys_lmdz_para, only: klon_omp
  implicit none
    
 
    allocate(airephy(klon_omp))
    allocate(cuphy(klon_omp))
    allocate(cvphy(klon_omp))
    allocate(rlatd(klon_omp))
    allocate(rlond(klon_omp))

  end subroutine initcomgeomphy
  
end module comgeomphy
