&ctl_nl
NThreads      = -1
partmethod    = 4
topology      = "cube"
test_case     = "held_suarez0"
ne            = NE
ndays         = 400
statefreq     = SFREQ
!theta_hydrostatic_mode = .true.
tstep_type    = 9
qsize         = 11
theta_advect_form = 1
pgrad_correction=1
hv_ref_profiles=6
hv_theta_correction=0
limiter_option = 9
restartfreq   =  1
restartfile   = "restart/R0001"
restartdir    = "./restart/"
runtype       = RUNTYPE
tstep         = TSTEP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EUL transport
!dt_remap_factor = 1
!dt_tracer_factor = 1
! SL transport:
transport_alg=12
semi_lagrange_cdr_alg=3 
semi_lagrange_nearest_point_lev = 256
semi_lagrange_hv_q = 1
dt_tracer_factor  = 6
hypervis_subcycle_q = 6
dt_remap_factor = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vert_remap_q_alg = 10
integration   = "explicit"
nu            = NU1
nu_top = 2.5e5
hypervis_scaling = 3
hypervis_order = 2
hypervis_subcycle = 1
hypervis_subcycle_tom = 1
se_ftype=0
/
&vert_nl
vfile_mid     = "../vcoord/camm-58.ascii"
vfile_int     = "../vcoord/cami-58.ascii"
/
&analysis_nl
infilenames=''
output_timeunits=1,0,2    ! 1=days, 2=hours, 3=seconds
output_frequency=0,0,0    ! 0 to disable
output_start_time=200,0,0
output_end_time=500,999999999,0
!output_varnames1='u','v','T','zeta','div','ps','geos','omega'
output_varnames1='u','v','T','omega','ps'
! debug output
output_varnames2='u','v','T','zeta','div','ps','geo','dp3d','geos','Th'
! output3: hourly data for 20 days  
output_varnames3='geos','omega','zeta','ps','div'
io_stride = 32
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/
