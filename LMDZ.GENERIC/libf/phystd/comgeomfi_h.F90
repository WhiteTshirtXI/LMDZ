
       module comgeomfi_h

       implicit none

       REAL,ALLOCATABLE,DIMENSION(:) :: long,lati,area
       REAL :: totarea, totarea_planet
!$OMP THREADPRIVATE(long,lati,area,totarea)

       end module comgeomfi_h

