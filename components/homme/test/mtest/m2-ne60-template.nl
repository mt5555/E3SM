! M2 test as in KSP2015, using 60 levels by default.
!_______________________________________________________________________
&ctl_nl
  theta_hydrostatic_mode = .false.
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "mtest2"                  ! test identifier: Circular Mountain, no shear
  dcmip_rayleigh_depth = @DCMIP_RAYLEIGH_DEPTH@  ! Rayleigh friction sponge depth (m) below model top, <0: test default, 0: off
  ne                = 60                        ! number of elements per cube face
  qsize             = 0                         ! num tracer fields
  ndays             = 0                         ! num simulation days: 0 => use nmax steps
  nmax              = 8000                      ! 7200s / 0.9s per step = 8000 steps
  statefreq         = 400                       ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.9                       ! largest timestep in seconds  0.9 works for 60L case
                                                ! L144 with sponger layer

  dt_remap_factor   = 1                         ! turn on vertical remap
  dt_tracer_factor  = 1
  vert_remap_q_alg    = 10
  vert_remap_u_alg    = 11   
                                                                
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 9                         ! 1 => default method
  nu                = 2.2e7                     ! reduced planet hyperviscosity hv/166.7^3
  nu_s              = 2.2e7
  nu_p              = 2.2e7
  nu_top            = @NU_TOP@                  ! top-layer Laplacian sponge coefficient
  tom_sponge_start  = @SPONGE_START@             ! reference pressure (hPa) where top sponge starts
                                                ! For 30km reference case, 20km = 100hPa.  10km=320hPa
                                                !  layer 10 = 50 hPa
  hypervis_subcycle_tom = 6
                                                ! 2.5e5 / 166.7 = 1500    8x stronger=12000
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
  output_frequency  = 360,                      ! 100s /0.1s = 1000 steps between outputs
  output_varnames1  ='T','ps','u','v','w','geo_i','pnh','omega'   ! variables to write to file
  interp_type       = 1                         ! 0=native grid, 1=bilinear
  interp_nlat       = 7                         ! just a few lon bands needed
  interp_nlon       = 512                       ! good for NE20, NE30, bump up to 1024 for NE60
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file  = .true.
/
