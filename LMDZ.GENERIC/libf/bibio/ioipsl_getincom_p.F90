!
! $Id$
!
! Module/Routines extracted from IOIPSL v2_1_8
!
MODULE ioipsl_getincom_p
!-
!$Id: getincom.f90 536 2009-01-30 11:46:27Z bellier $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
USE ioipsl_getincom ! to use  'getin'
USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
USE mod_phys_lmdz_omp_data, ONLY :  is_omp_root
use mod_phys_lmdz_para, only : bcast
!-
IMPLICIT NONE
!-
PRIVATE
PUBLIC :: getin_p
!-
INTERFACE getin_p

  MODULE PROCEDURE getinrs_p, getinr1d_p, getinr2d_p, &
 &                 getinis_p, getini1d_p, getini2d_p, &
 &                 getincs_p, 		              &
 &                 getinls_p, getinl1d_p, getinl2d_p
END INTERFACE
!-
CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Definition des getin -> bcast      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Les chaines de caracteres -- !!
  
  SUBROUTINE getincs_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    CHARACTER(LEN=*),INTENT(INOUT) :: VarOut    

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getincs_p

!! -- Les entiers -- !!
  
  SUBROUTINE getinis_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut    

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getinis_p

  SUBROUTINE getini1d_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut(:)

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getini1d_p

  SUBROUTINE getini2d_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut(:,:)

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getini2d_p

!! -- Les flottants -- !!
  
  SUBROUTINE getinrs_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getinrs_p

  SUBROUTINE getinr1d_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut(:)

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getinr1d_p

  SUBROUTINE getinr2d_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut(:,:)

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getinr2d_p

!! -- Les Booleens -- !!
  
  SUBROUTINE getinls_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getinls_p

  SUBROUTINE getinl1d_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut(:)

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getinl1d_p

  SUBROUTINE getinl2d_p(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut(:,:)

!$OMP BARRIER
    IF (is_mpi_root .AND. is_omp_root) THEN
    	CALL getin(VarIn,VarOut)
    ENDIF
    CALL bcast(VarOut)
  END SUBROUTINE getinl2d_p
!-
!-----------------------------
!-----------------------------
!-----------------------------

END MODULE ioipsl_getincom_p
