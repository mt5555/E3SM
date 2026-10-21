module elem_ops_interface

  use iso_c_binding,  only: c_int, c_bool, c_ptr, c_f_pointer
  use dimensions_mod, only: nlev, nlevp, np
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use parallel_mod,   only: abortmp

  implicit none

  type(hvcoord_t) :: hvcoord
  public :: init_f90
  public :: compute_r_star_f90
  public :: compute_r_star_wl_f90
contains

  subroutine init_f90 (hyai, ps0) bind(c)
    !
    ! Inputs
    !
    real (kind=real_kind), intent(in) :: hyai(nlevp)
    real (kind=real_kind), intent(in) :: ps0

    hvcoord%hyai = hyai
    hvcoord%ps0 = ps0

  end subroutine init_f90

  subroutine compute_r_star_f90(num_elems, moist, Q_ptr, R_ptr) bind(c)
    use element_ops, only: get_R_star
    use control_mod, only: use_moisture
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: num_elems
    logical (kind=c_bool), intent(in) :: moist
    type (c_ptr), intent(in) :: Q_ptr, R_ptr

    !
    ! Locals
    !
    real (kind=real_kind), dimension(:,:,:,:),  pointer :: Q, R
    integer :: ie

    use_moisture = moist

    call c_f_pointer(Q_ptr, Q, [np,np,nlev,num_elems])
    call c_f_pointer(R_ptr, R, [np,np,nlev,num_elems])

    do ie=1,num_elems
      call get_R_star(R(:,:,:,ie),Q(:,:,:,ie))
    enddo
  end subroutine compute_R_star_f90

  ! R* from the dry mass, with 3 tracers: vapor (1) and condensates (3,2)
  subroutine compute_r_star_wl_f90(num_elems, moist, wl, dp_ptr, q1_ptr, q2_ptr, q3_ptr, R_ptr) bind(c)
    use element_ops, only: get_dry_mass, get_R_star
    use control_mod, only: use_moisture, water_loading, num_condensates, condensate_idx
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: num_elems
    logical (kind=c_bool), intent(in) :: moist, wl
    type (c_ptr), intent(in) :: dp_ptr, q1_ptr, q2_ptr, q3_ptr, R_ptr

    !
    ! Locals
    !
    real (kind=real_kind), dimension(:,:,:,:),  pointer :: dp, q1, q2, q3, R
    real (kind=real_kind) :: Q(np,np,nlev,3), dp_dry(np,np,nlev)
    integer :: ie

    use_moisture = moist
    water_loading = wl
    num_condensates = 2
    condensate_idx = 0
    condensate_idx(1:2) = [3,2]

    call c_f_pointer(dp_ptr, dp, [np,np,nlev,num_elems])
    call c_f_pointer(q1_ptr, q1, [np,np,nlev,num_elems])
    call c_f_pointer(q2_ptr, q2, [np,np,nlev,num_elems])
    call c_f_pointer(q3_ptr, q3, [np,np,nlev,num_elems])
    call c_f_pointer(R_ptr,  R,  [np,np,nlev,num_elems])

    do ie=1,num_elems
      Q(:,:,:,1) = q1(:,:,:,ie)
      Q(:,:,:,2) = q2(:,:,:,ie)
      Q(:,:,:,3) = q3(:,:,:,ie)
      call get_dry_mass(dp_dry,dp(:,:,:,ie),Q)
      call get_R_star(R(:,:,:,ie),dp_dry,dp(:,:,:,ie),Q(:,:,:,1))
    enddo

    water_loading = .false.
    num_condensates = 0
    condensate_idx = 0
  end subroutine compute_r_star_wl_f90

end module elem_ops_interface
