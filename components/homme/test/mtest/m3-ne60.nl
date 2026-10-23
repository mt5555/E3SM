! M3 test as in KSP2015, runs with 60 levels. 
!_______________________________________________________________________
&ctl_nl
  theta_hydrostatic_mode = .false.
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "mtest3"       ! test identifier: Circular Mountain, with shear
  dcmip_rayleigh_depth = -1                  ! Rayleigh friction sponge depth (m) below model top, <0: test default, 0: off
  ne                = 60                        ! number of elements per cube face
  qsize             = 0                         ! num tracer fields
  ndays             = 0                         ! num simulation days: 0 => use nmax steps
  nmax              = 8000                      ! 7200s / 0.9s per step = 8000 steps
  statefreq         = 400                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.9                       ! largest timestep in seconds
  dt_remap_factor   = 1                         ! turn on vertical remap
  dt_tracer_factor  = 1
  vert_remap_q_alg  = 10
  vert_remap_u_alg  = 11
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 9                         ! 1 => default method
  nu                = 2.2e7                       ! reduced planet hyperviscosity hv/166.7^3
  nu_s              = 2.2e7
  nu_p              = 2.2e7
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
