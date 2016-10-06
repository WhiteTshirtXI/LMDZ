      subroutine writediagfi(ngrid,nom,titre,unite,dim,px)

!  Ecriture de variables diagnostiques au choix dans la physique 
!  dans un fichier NetCDF nomme  'diagfi'. Ces variables peuvent etre
!  3d (ex : temperature), 2d (ex : temperature de surface), ou
!  0d (pour un scalaire qui ne depend que du temps : ex : la longitude
!  solaire)
!  (ou encore 1d, dans le cas de testphys1d, pour sortir une colonne)
!  La periode d'ecriture est donnee par 
!  "ecritphy " regle dans le fichier de controle de run :  run.def
!
!    writediagfi peut etre appele de n'importe quelle subroutine
!    de la physique, plusieurs fois. L'initialisation et la creation du
!    fichier se fait au tout premier appel.
!
! WARNING : les variables dynamique (u,v,t,q,ps)
!  sauvees par writediagfi avec une
! date donnee sont legerement differentes que dans le fichier histoire car 
! on ne leur a pas encore ajoute de la dissipation et de la physique !!!
! IL est  RECOMMANDE d'ajouter les tendance physique a ces variables
! avant l'ecriture dans diagfi (cf. physiq.F)
!  
! Modifs: Aug.2010 Ehouarn: enforce outputs to be real*4
!         Oct 2011 Francois: enable having a 'diagfi.def' file to select
!                            at runtime, which variables to put in file
!
!  parametres (input) :
!  ----------
!      ngrid : nombres de point ou est calcule la physique
!                (ngrid = 2+(jjm-1)*iim - 1/jjm)
!                 (= nlon ou klon dans la physique terrestre)
!      
!      unit : unite logique du fichier de sortie (toujours la meme)
!      nom  : nom de la variable a sortir (chaine de caracteres)
!      titre: titre de la variable (chaine de caracteres)
!      unite : unite de la variable (chaine de caracteres)
!      px : variable a sortir (real 0, 1, 2, ou 3d)
!      dim : dimension de px : 0, 1, 2, ou 3 dimensions
!
!=================================================================
      use surfdat_h, only: phisfi
      use control_mod, only: ecritphy, day_step, iphysiq
      USE mod_phys_lmdz_para, only : is_parallel, is_mpi_root,
     &                               is_master, gather
      USE mod_grid_phy_lmdz, only : klon_glo, Grid1Dto2D_glo
      implicit none

! Commons
      include "dimensions.h"
      include "paramet.h"
      include "comvert.h"
      include "comgeom.h"
      include "netcdf.inc"
      include "temps.h"

! Arguments on input:
      integer,intent(in) :: ngrid
      character (len=*),intent(in) :: nom,titre,unite
      integer,intent(in) :: dim
      real,intent(in) :: px(ngrid,llm)

! Local variables:

      real*4 dx3(iip1,jjp1,llm) ! to store a 3D data set
      real*4 dx2(iip1,jjp1)     ! to store a 2D (surface) data set
      real*4 dx1(llm)           ! to store a 1D (column) data set
      real*4 dx0

      real*4,save :: date
!$OMP THREADPRIVATE(date)

      REAL phis(ip1jmp1)

      integer irythme
      integer ierr,ierr2
      integer iq
      integer i,j,l,zmax , ig0

      integer,save :: zitau=0
      character(len=20),save :: firstnom='1234567890'
!$OMP THREADPRIVATE(zitau,firstnom)

! Ajouts
      integer, save :: ntime=0
!$OMP THREADPRIVATE(ntime)
      integer :: idim,varid
      integer :: nid
      character(len=*),parameter :: fichnom="diagfi.nc"
      integer, dimension(4) :: id
      integer, dimension(4) :: edges,corner

! Added to use diagfi.def to select output variable
      logical,save :: diagfi_def
      logical :: getout
      integer,save :: n_nom_def
      integer :: n
      integer,parameter :: n_nom_def_max=199
      character(len=120),save :: nom_def(n_nom_def_max)
      logical,save :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall) 	!diagfi_def,n_nom_def,nom_def read in diagfi.def
      
#ifndef MESOSCALE

#ifdef CPP_PARA
! Added to work in parallel mode
      real dx3_glop(klon_glo,llm)
      real dx3_glo(iim,jjp1,llm) ! to store a global 3D data set
      real dx2_glop(klon_glo)
      real dx2_glo(iim,jjp1)     ! to store a global 2D (surface) data set
      real px2(ngrid)
!      real dx1_glo(llm)          ! to store a 1D (column) data set
!      real dx0_glo
      real phisfi_glo(klon_glo) ! surface geopotential on global physics grid
#else
      real phisfi_glo(ngrid) ! surface geopotential on global physics grid
#endif

!***************************************************************
!Sortie des variables au rythme voulu

      irythme = ecritphy ! sortie au rythme de ecritphy
!     irythme = iconser  ! sortie au rythme des variables de controle
!     irythme = iphysiq  ! sortie a tous les pas physique
!     irythme = iecri*day_step ! sortie au rythme des fichiers histoires
!     irythme = periodav*day_step ! sortie au rythme des fichiers histmoy

!***************************************************************

! At very first call, check if there is a "diagfi.def" to use and read it
! -----------------------------------------------------------------------
      IF (firstcall) THEN
         firstcall=.false.

!$OMP MASTER
  !      Open diagfi.def definition file if there is one:
         open(99,file="diagfi.def",status='old',form='formatted',
     s   iostat=ierr2)

         if(ierr2.eq.0) then
            diagfi_def=.true.
            write(*,*) "******************"
            write(*,*) "Reading diagfi.def"
            write(*,*) "******************"
            do n=1,n_nom_def_max
              read(99,fmt='(a)',end=88) nom_def(n)
              write(*,*) 'Output in diagfi: ', trim(nom_def(n))
            end do 
 88         continue
            if (n.ge.n_nom_def_max) then
               write(*,*)"n_nom_def_max too small in writediagfi.F:",n
               stop
            end if 
            n_nom_def=n-1
            close(99)
         else
            diagfi_def=.false.
         endif
!$OMP END MASTER
!$OMP BARRIER
      END IF ! of IF (firstcall)

! Get out of write_diagfi if there is diagfi.def AND variable not listed
!  ---------------------------------------------------------------------
      if (diagfi_def) then
          getout=.true.
          do n=1,n_nom_def
             if(trim(nom_def(n)).eq.nom) getout=.false.
          end do
          if (getout) return
      end if

! Initialisation of 'firstnom' and create/open the "diagfi.nc" NetCDF file
! ------------------------------------------------------------------------
! (at very first call to the subroutine, in accordance with diagfi.def)

      if (firstnom.eq.'1234567890') then ! .true. for the very first valid
      !   call to this subroutine; now set 'firstnom'
         firstnom = nom
         ! just to be sure, check that firstnom is large enough to hold nom
         if (len_trim(firstnom).lt.len_trim(nom)) then
           write(*,*) "writediagfi: Error !!!"
           write(*,*) "   firstnom string not long enough!!"
           write(*,*) "   increase its size to at least ",len_trim(nom)
           stop
         endif
         
         zitau = -1 ! initialize zitau

#ifdef CPP_PARA
          ! Gather phisfi() geopotential on physics grid
          call Gather(phisfi,phisfi_glo)
#else
         phisfi_glo(:)=phisfi(:)
#endif

         !! parallel: we cannot use the usual writediagfi method
!!         call iophys_ini
         if (is_master) then
         ! only the master is required to do this

         ! Create the NetCDF file
         ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
         ! Define the 'Time' dimension
         ierr = nf_def_dim(nid,"Time",NF_UNLIMITED,idim)
         ! Define the 'Time' variable
!#ifdef NC_DOUBLE
!         ierr = NF_DEF_VAR (nid, "Time", NF_DOUBLE, 1, idim,varid)
!#else
         ierr = NF_DEF_VAR (nid, "Time", NF_FLOAT, 1, idim,varid)
!#endif
         ! Add a long_name attribute
         ierr = NF_PUT_ATT_TEXT (nid, varid, "long_name",
     .          4,"Time")
         ! Add a units attribute
         ierr = NF_PUT_ATT_TEXT(nid, varid,'units',29,
     .          "days since 0000-00-0 00:00:00")
         ! Switch out of NetCDF Define mode
         ierr = NF_ENDDEF(nid)

         ! write "header" of file (longitudes, latitudes, geopotential, ...)
         call gr_fi_dyn(1,size(phisfi_glo),iip1,jjp1,phisfi_glo,phis)
         call iniwrite(nid,day_ini,phis)

         endif ! of if (is_master)

      else

         if (is_master) then
           ! only the master is required to do this

           ! Open the NetCDF file
           ierr = NF_OPEN(fichnom,NF_WRITE,nid)
         endif ! of if (is_master)

      endif ! if (firstnom.eq.'1234567890')

      if (ngrid.eq.1) then
        ! in testphys1d, for the 1d version of the GCM, iphysiq and irythme
        ! are undefined; so set them to 1
        iphysiq=1 
        irythme=1
        ! NB: 
      endif

! Increment time index 'zitau' if it is the "fist call" (at given time level)
! to writediagfi
!------------------------------------------------------------------------
      if (nom.eq.firstnom) then
          zitau = zitau + iphysiq
      end if

!--------------------------------------------------------
! Write the variables to output file if it's time to do so
!--------------------------------------------------------

      if ( MOD(zitau+1,irythme) .eq.0.) then

! Compute/write/extend 'Time' coordinate (date given in days)
! (done every "first call" (at given time level) to writediagfi)
! Note: date is incremented as 1 step ahead of physics time
!--------------------------------------------------------

        if (is_master) then
           ! only the master is required to do this
        if (nom.eq.firstnom) then
        ! We have identified a "first call" (at given date)
           ntime=ntime+1 ! increment # of stored time steps
           ! compute corresponding date (in days and fractions thereof)
           date= float (zitau +1)/float (day_step)
           ! Get NetCDF ID of 'Time' variable
           ierr= NF_INQ_VARID(nid,"Time",varid)
           ! Write (append) the new date to the 'Time' array
!#ifdef NC_DOUBLE
!           ierr= NF_PUT_VARA_DOUBLE(nid,varid,ntime,1,date)
!#else
           ierr= NF_PUT_VARA_REAL(nid,varid,ntime,1,date)
!#endif
           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi_nc"
              write(*,*) "***** with time"
              write(*,*) 'ierr=', ierr   
c             call abort
           endif

           write(6,*)'WRITEDIAGFI: date= ', date
        end if ! of if (nom.eq.firstnom)

        endif ! of if (is_master)

!Case of a 3D variable
!---------------------
        if (dim.eq.3) then

#ifdef CPP_PARA
          ! Gather field on a "global" (without redundant longitude) array
          call Gather(px,dx3_glop)
!$OMP MASTER
          if (is_mpi_root) then
            call Grid1Dto2D_glo(dx3_glop,dx3_glo)
            ! copy dx3_glo() to dx3(:) and add redundant longitude
            dx3(1:iim,:,:)=dx3_glo(1:iim,:,:)
            dx3(iip1,:,:)=dx3(1,:,:)
          endif
!$OMP END MASTER
!$OMP BARRIER
#else
!         Passage variable physique -->  variable dynamique
!         recast (copy) variable from physics grid to dynamics grid
           DO l=1,llm
             DO i=1,iip1
                dx3(i,1,l)=px(1,l)
                dx3(i,jjp1,l)=px(ngrid,l)
             ENDDO
             DO j=2,jjm
                ig0= 1+(j-2)*iim
                DO i=1,iim
                   dx3(i,j,l)=px(ig0+i,l)
                ENDDO
                dx3(iip1,j,l)=dx3(1,j,l)
             ENDDO
           ENDDO
#endif
!         Ecriture du champs

          if (is_master) then
           ! only the master writes to output
! name of the variable
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"altitude",id(3))
              ierr= NF_INQ_DIMID(nid,"Time",id(4))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",nom
              call def_var(nid,nom,titre,unite,4,id,varid,ierr)

           endif 

           corner(1)=1
           corner(2)=1
           corner(3)=1
           corner(4)=ntime

           edges(1)=iip1
           edges(2)=jjp1
           edges(3)=llm
           edges(4)=1
!#ifdef NC_DOUBLE
!           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx3)
!#else
!           write(*,*)"test:  nid=",nid," varid=",varid
!           write(*,*)"       corner()=",corner
!           write(*,*)"       edges()=",edges
!           write(*,*)"       dx3()=",dx3
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx3)
!#endif

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR problem in writediagfi"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
c             call abort
           endif 

          endif !of if (is_master)

!Case of a 2D variable
!---------------------

        else if (dim.eq.2) then

#ifdef CPP_PARA
          ! Gather field on a "global" (without redundant longitude) array
          px2(:)=px(:,1)
          call Gather(px2,dx2_glop)
!$OMP MASTER
          if (is_mpi_root) then
            call Grid1Dto2D_glo(dx2_glop,dx2_glo)
            ! copy dx2_glo() to dx2(:) and add redundant longitude
            dx2(1:iim,:)=dx2_glo(1:iim,:)
            dx2(iip1,:)=dx2(1,:)
          endif
!$OMP END MASTER
!$OMP BARRIER
#else

!         Passage variable physique -->  physique dynamique
!         recast (copy) variable from physics grid to dynamics grid

             DO i=1,iip1
                dx2(i,1)=px(1,1)
                dx2(i,jjp1)=px(ngrid,1)
             ENDDO
             DO j=2,jjm
                ig0= 1+(j-2)*iim
                DO i=1,iim
                   dx2(i,j)=px(ig0+i,1)
                ENDDO
                dx2(iip1,j)=dx2(1,j)
             ENDDO
#endif

          if (is_master) then
           ! only the master writes to output
!         write (*,*) 'In  writediagfi, on sauve:  ' , nom
!         write (*,*) 'In  writediagfi. Estimated date = ' ,date
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"Time",id(3))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",nom

              call def_var(nid,nom,titre,unite,3,id,varid,ierr)

           endif

           corner(1)=1
           corner(2)=1
           corner(3)=ntime
           edges(1)=iip1
           edges(2)=jjp1
           edges(3)=1


!#ifdef NC_DOUBLE
!           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,dx2) 
!#else         
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx2)
!#endif     

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
c             call abort
           endif 

          endif !of if (is_master)

!Case of a 1D variable (ie: a column)
!---------------------------------------------------

       else if (dim.eq.1) then
         if (is_parallel) then
           write(*,*) "writediagfi error: dim=1 not implemented ",
     &                 "in parallel mode"
           stop
         endif
!         Passage variable physique -->  physique dynamique
!         recast (copy) variable from physics grid to dynamics grid
          do l=1,llm
            dx1(l)=px(1,l)
          enddo
          
          ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"altitude",id(1))
              ierr= NF_INQ_DIMID(nid,"Time",id(2))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",nom

              call def_var(nid,nom,titre,unite,2,id,varid,ierr)
              
           endif
           
           corner(1)=1
           corner(2)=ntime
           
           edges(1)=llm
           edges(2)=1
!#ifdef NC_DOUBLE
!           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx1)
!#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx1)
!#endif

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR problem in writediagfi"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
c             call abort
           endif 

!Case of a 0D variable (ie: a time-dependent scalar)
!---------------------------------------------------

        else if (dim.eq.0) then

           dx0 = px (1,1)

          if (is_master) then
           ! only the master writes to output
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"Time",id(1))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",nom

              call def_var(nid,nom,titre,unite,1,id,varid,ierr)

           endif

           corner(1)=ntime
           edges(1)=1

!#ifdef NC_DOUBLE
!           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,dx0)  
!#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx0)
!#endif
           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
c             call abort
           endif 

          endif !of if (is_master)

        endif ! of if (dim.eq.3) elseif(dim.eq.2)...

      endif ! of if ( MOD(zitau+1,irythme) .eq.0.)

      if (is_master) then
        ierr= NF_CLOSE(nid)
      endif

#endif
      end
