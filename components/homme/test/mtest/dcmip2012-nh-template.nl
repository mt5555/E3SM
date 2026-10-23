!
! namelist for dcmip2012 test2-1: nonhydro mountain waves without shear
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2012_test2_1"       ! test identifier
  dcmip_rayleigh_depth = @DCMIP_RAYLEIGH_DEPTH@  ! Rayleigh friction sponge depth (m) below model top, <0: test default, 0: off
  theta_hydrostatic_mode = .false.
  rsplit            = 1
  ne                = 20                        ! number of elements per cube face
  qsize             = 0                         ! num tracer fields
  nmax              = 72000                     ! 7200s / 0.1s per step = 72000 steps
  statefreq         = 360                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.1                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5                         ! 1 => default method
  nu                = 3.2e7                     ! reduced planet hyperviscosity hv/500^3
  nu_s              = 3.2e7
  nu_p              = 3.2e7                     ! like 4e15
  nu_top            = @NU_TOP@                   ! 2.5e5 / 500 = 500    8x stronger=4000
  tom_sponge_start  = @SPONGE_START@            ! TOM is about 30mb
  hypervis_subcycle_tom = 1
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 12752.0                   ! reduced planet radius rearth = a/500.0
  omega             = 0.0                       ! earth angular speed = 0.0
  dcmip2_x_ueq      = 20.0                      ! windspeed at equator  (m/s)
  dcmip2_x_h0       = 250.0                     ! mountain height       (m)
  dcmip2_x_d        = 5000.0                    ! mountain half-width   (m)
  dcmip2_x_xi       = 4000.0                    ! mountain wavelength   (m)
/
&vert_nl
  vanalytic         = @VANALYTIC@                ! 1=analytic vcoords (default), 0=read from vfile_mid/vfile_int
  vtop              = 3.2818e-2                 ! vertical coordinate at top of atm (z=30km), ignored when vanalytic=0
  vfile_mid = '@VFILE_MID@'                     ! ignored when vanalytic=1
  vfile_int = '@VFILE_INT@'                     ! ignored when vanalytic=1
/
&analysis_nl
  output_dir        = "./@MOVIES_DIR@/"         ! destination dir for netcdf file
  output_timeunits  = 3,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 720,                      ! 720 seconds (10+1 outputs)
  output_varnames1  ='T','ps','u','v','w','geo_i','pnh'   ! variables to write to file
  interp_type       = 1                         ! 0=native grid, 1=bilinear
  interp_nlat       = 7                         ! just a few lon bands needed
  interp_nlon       = 512                       ! good for NE20, NE30, bump up to 1024 for NE60
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
