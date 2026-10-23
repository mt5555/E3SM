! M1 tests as in KSP2015, runs with 60 levels.
!_______________________________________________________________________
&ctl_nl
  theta_hydrostatic_mode = .false.
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "mtest1"       ! test identifier
  dcmip_rayleigh_depth = -1                  ! Rayleigh friction sponge depth (m) below model top, <0: test default, 0: off
  ne                = 30                        ! number of elements per cube face
  qsize             = 0                         ! num tracer fields
  ndays             = 0                         ! num simulation days: 0 => use nmax steps
  nmax              = 4000                      ! 7200s / 1.8s per step = 4000 steps
  statefreq         = 360                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 1.8                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 9                         ! 1 => default method
  nu                = 2.2e8                       ! reduced planet hyperviscosity hv/166.7^3
  nu_s              = 2.2e8
  nu_p              = 2.2e8
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 38248.0                   ! reduced planet radius rearth = a/166.7
  omega             = 0.0                       ! earth angular speed = 0.0
/
&vert_nl
  vanalytic         = 1                         ! set vcoords in initialization routine
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 3,                        ! 1=days, 2=hours, 3=seconds, 0=timesteps
  output_frequency  = 360,                      ! seconds between outputs
  output_varnames1  ='T','ps','u','v','w','geo_i','pnh'   ! variables to write to file
  interp_type       = 1                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16         
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
