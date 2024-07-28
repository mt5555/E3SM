#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_restart_mod 
   !------------------
   use kinds , only : real_kind, int_kind
   !------------------
   use dimensions_mod, only : nelemd,np,nlev,nlevp,qsize
   !------------------
   use parallel_mod, only : parallel_t, MPIreal_t, abortmp
   !------------------
   use time_mod, only : timelevel_t
   !------------------
   use element_state, only : elem_state_t
   !------------------
   use element_mod,      only: element_t
   !------------------
   use element_ops,      only: get_R_star
   !------------------
   use restart_io_mod, only : nwordsRestartBuffer_t, RestartBuffer,  File_elem_t, &
        StateDesc_t, createstatedescriptor, AddStateField, constructelementfile, &
        collective_io_read, collective_io_write, printstatedescriptor
   !------------------
   use interpolate_driver_mod, only : pio_read_var
   !------------------
   use physical_constants, only : Rgas
   implicit none

private 
   ! =========================================
   !  Some variables used by all MPI routines
   ! =========================================
   integer                         :: errorcode,errorlen,ierr
   character(len=80)               :: errorstring
   ! ====================================================
   !  Routines for Restart files
   ! ====================================================
   public :: initRestartFile,readstate_uniquepts

contains 
! =========================================================
! initRestartFile:
!
!  Initalizes MPI I-O  by seting up some MPI datastructures
! =========================================================
   subroutine initRestartFile(state, par,File)
    type (elem_state_t), intent(in) :: state
    type (parallel_t),intent(in)    :: par
    type (File_elem_t),intent(out) :: File

    integer                      :: ie,ig,ierr

    integer                      :: count
    integer,allocatable          :: blklen(:),disp(:),oldtype(:)
    integer                      :: len
    type (StateDesc_t)           :: RestDesc
    integer                      :: NumComponents
    integer                      :: type
    integer                      :: i
    logical,parameter            :: Debug = .FALSE.

    !=========================================
    !  Copy over the parallel_t datastructure
    !=========================================
    File%par = par

    !=========================================
    !  Allocate restart buffer
    !=========================================
    if(Debug) print *,'initRestartFile: nelemd: ',nelemd

    collective_io_read = .true.
    collective_io_write = .true.

    allocate(RestartBuffer(nelemd))

    !================================================================
    !   Construct the descriptor of the state variable for MPI I/O
    !================================================================
    NumComponents = 9   ! THIS NUMBER MUST MATCH NUMBER OF TIMES AddStateField is called below
    RestDesc = CreateStateDescriptor(NumComponents)

    type = MPIReal_t       !  All the types are Real
    !=========================================
    ! for PRESTART, must add *all* the fields in the State variable
    !=========================================
    len = SIZE(state%v)
    call AddStateField(RestDesc,len,type)
#ifdef MODEL_THETA_L
    len = SIZE(state%w_i)
#else
    len = SIZE(state%w)
#endif
    call AddStateField(RestDesc,len,type)
#ifdef MODEL_THETA_L
    len = SIZE(state%vtheta_dp)
#else
    len = SIZE(state%theta_dp_cp)
#endif
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%ps_v)
    call AddStateField(RestDesc,len,type)

#ifdef MODEL_THETA_L
    len = SIZE(state%phinh_i)
#else
    len = SIZE(state%phinh)
#endif
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%phis)
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%Q)
    call AddStateField(RestDesc,len,type)
    len = SIZE(state%Qdp)
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%dp3d)
    call AddStateField(RestDesc,len,type)

#if defined(_MPI) && defined(_PRESTART)
    if(Debug) call PrintStateDescriptor(RestDesc)
    call ConstructElementFile(RestDesc,File,ierr)
#endif
    nwordsRestartBuffer_t=RestDesc%nwords

    end subroutine initRestartFile



!
!   read state from a regular uniquepts output file
!
    subroutine readstate_uniquepts(elem, par,tl,infilenames_index)
    implicit none    
    type (element_t),   intent(inout)     :: elem(:)
    type (TimeLevel_t), intent(in)        :: tl     ! time level struct
    type(parallel_t),   intent(in)        :: par
    integer,            intent(in)        :: infilenames_index

    ! dont do this - will allocate for every thread:
    !real(kind=real_kind) :: temp3d(np,np,1,nelemd)
    real(kind=real_kind),allocatable :: temp3d(:,:,:,:)
    real(kind=real_kind) :: Rstar(np,np,nlev)
    integer :: ie,qindex
    character(len=2) :: vname

!$OMP BARRIER
!$OMP MASTER
    ! only coded for single threaded
    ! read Q
    ! read (into tl%n0): u,v,w_i,vtheta, phinh_i, dp3d, ps_v
    ! convert vtheta to vtheta_dp
    ! phis set during init.  (dont read to restart with different phis)
    ! driver code will populate Qdp from Q

    allocate(temp3d(np,np,1,nelemd))
    call pio_read_var(temp3d,elem,par,'ps',1,infilenames_index)
    do ie=1,nelemd
       elem(ie)%state%ps_v(:,:,tl%n0) = temp3d(:,:,1,ie)
    enddo
    deallocate(temp3d)

    allocate(temp3d(np,np,nlevp,nelemd))
    call pio_read_var(temp3d,elem,par,'geo_i',nlevp,infilenames_index)
    do ie=1,nelemd
       elem(ie)%state%phinh_i(:,:,:,tl%n0) = temp3d(:,:,:,ie)
    enddo
    call pio_read_var(temp3d,elem,par,'w_i',nlevp,infilenames_index)
    do ie=1,nelemd
       elem(ie)%state%w_i(:,:,:,tl%n0) = temp3d(:,:,:,ie)
    enddo
    deallocate(temp3d)

    allocate(temp3d(np,np,nlev,nelemd))
    call pio_read_var(temp3d,elem,par,'u',nlev,infilenames_index)
    do ie=1,nelemd
       elem(ie)%state%v(:,:,1,:,tl%n0) = temp3d(:,:,:,ie)
    enddo
    call pio_read_var(temp3d,elem,par,'v',nlev,infilenames_index)
    do ie=1,nelemd
       elem(ie)%state%v(:,:,2,:,tl%n0) = temp3d(:,:,:,ie)
    enddo
    call pio_read_var(temp3d,elem,par,'dp3d',nlev,infilenames_index)
    do ie=1,nelemd
       elem(ie)%state%dp3d(:,:,:,tl%n0) = temp3d(:,:,:,ie)
    enddo
    do qindex=1,min(qsize,4)  ! prim_movie_mod() can only output Q1..Q4
       write(vname,'(a1,i1)') 'Q',qindex
       call pio_read_var(temp3d,elem,par,vname,nlev,infilenames_index)
       do ie=1,nelemd
          elem(ie)%state%Q(:,:,:,qindex) = temp3d(:,:,:,ie)
       enddo
    enddo
    call pio_read_var(temp3d,elem,par,'Th',nlev,infilenames_index)
    do ie=1,nelemd
       call get_R_star(Rstar,elem(ie)%state%Q(:,:,:,1))
       elem(ie)%state%vtheta_dp(:,:,:,tl%n0) = &
            Rstar(:,:,:)*elem(ie)%state%dp3d(:,:,:,tl%n0)*temp3d(:,:,:,ie)/Rgas
    enddo
    deallocate(temp3d)

!$OMP END MASTER
!$OMP BARRIER
    end subroutine readstate_uniquepts




end module prim_restart_mod
