! M1 tests as in KSP2015, runs with 60 levels.
!_______________________________________________________________________
&ctl_nl
  theta_hydrostatic_mode = .false.
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "mtest1"       ! test identifier: Quasi Two-Dimensional Mountain Ridge
  dcmip_rayleigh_depth = @DCMIP_RAYLEIGH_DEPTH@  ! Rayleigh friction sponge depth (m) below model top, <0: test default, 0: off
  ne                = 30                        ! number of elements per cube face
  qsize             = 0                         ! num tracer fields
  ndays             = 0                         ! num simulation days: 0 => use nmax steps
  nmax              = 4000                      ! 7200s / 1.8s per step = 4000 steps
  statefreq         = 360                      ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 1.8                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 9                         ! 1 => default method
  nu                = 2.2e8                       ! reduced planet hyperviscosity hv/166.7^3
  nu_s              = 2.2e8
  nu_p              = 2.2e8
  nu_top            = @NU_TOP@                  ! top-layer Laplacian sponge coefficient
                                                ! 2.5e5 / 166.7 = 1500    8x stronger=12000
  tom_sponge_start  = @SPONGE_START@             ! reference pressure (hPa) where top sponge starts
  hypervis_subcycle_tom = 1
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 38248.0                   ! reduced planet radius rearth = a/166.7
  omega             = 0.0                       ! earth angular speed = 0.0
/
&vert_nl
  vanalytic         = @VANALYTIC@                ! 1=analytic vcoords (default), 0=read from vfile_mid/vfile_int
  vfile_mid = '@VFILE_MID@'                     ! ignored when vanalytic=1
  vfile_int = '@VFILE_INT@'                     ! ignored when vanalytic=1
/
&analysis_nl
  output_dir        = "./@MOVIES_DIR@/"         ! destination dir for netcdf file
  output_timeunits  = 3,                        ! 1=days, 2=hours, 0=timesteps
  output_frequency  = 360,                     ! 100s /0.1s = 1000 steps between outputs
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
