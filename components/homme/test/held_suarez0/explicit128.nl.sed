&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
sub_case = 2
ne            = NE
ndays         = 400
statefreq     = SFREQ
!theta_hydrostatic_mode = .true.
!tstep_type    = 5
theta_hydrostatic_mode = .false.
tstep_type    = 9
qsize         = 1
theta_advect_form = 1
pgrad_correction=1
hv_ref_profiles=6
hv_theta_correction=0
limiter_option = 9
dt_remap_factor = 1
dt_tracer_factor = 1
vert_remap_q_alg = 10
!vert_remap_u_alg = 11
restartfreq   =  1
restartfile   = "restart/R0001"
restartdir    = "./restart/"
runtype       = RUNTYPE
tstep         = TSTEP
integration   = "explicit"
nu            = NU1
!nu_s        = 4.5e-8
!nu_p        = 4.5e-8
nu_top = 0e5  ! default 2.5e5    HSV1 1.5ok.  2.0 bad
                ! timesplit version.  5e5 works.  10e5 crashes.  
hypervis_scaling = 3  ! 0 for constant coeff HV
hypervis_order = 2
hypervis_subcycle = 1
hypervis_subcycle_tom = 1
/
&vert_nl
vfile_mid     = "../vcoord/scream-128m.ascii"
vfile_int     = "../vcoord/scream-128i.ascii"
/
&analysis_nl
infilenames=''
output_timeunits=2,0,2    ! 1=days, 2=hours, 3=seconds
output_frequency=0,0,0    ! 0 to disable
output_start_time=2520,0,2424        ! 200h*24 = 4800
output_end_time=3240,999999999,2424
!output_varnames1='u','v','T','zeta','div','ps','geos','omega'
output_varnames1='u','v','T','omega','ps'
! debug output
output_varnames2='u','v','T','zeta','div','ps','geo','dp3d','geos','Th'
! output3: hourly data for 20 days  
output_varnames3='omega','ps'
io_stride = 32
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/
