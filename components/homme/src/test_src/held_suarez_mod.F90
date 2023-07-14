#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module held_suarez_mod


  use coordinate_systems_mod, only: spherical_polar_t
  use dimensions_mod,         only: nlev,np,qsize,nlevp
  use element_mod,            only: element_t
  use element_state,          only: timelevels
  use element_ops,            only: set_thermostate, get_temperature, get_phi
  use hybrid_mod,             only: hybrid_t
  use hybvcoord_mod,          only: hvcoord_t
  use kinds,                  only: real_kind, iulog
  use physical_constants,     only: p0, kappa,g, dd_pi, Rgas, TREF, g
  use control_mod,            only: sub_case
  use physics_mod,            only: prim_condense
  use time_mod,               only: secpday
#ifndef HOMME_WITHOUT_PIOLIBRARY
  use common_io_mod,          only: infilenames
#endif

implicit none
private

  real (kind=real_kind), public, parameter :: sigma_b = 0.70D0
  real (kind=real_kind), public, parameter :: k_a     = 1.0D0/(40.0D0*secpday)
  real (kind=real_kind), public, parameter :: k_f     = 1.0D0/(1.0D0*secpday)
  real (kind=real_kind), public, parameter :: k_s     = 1.0D0/(4.0D0*secpday)
  real (kind=real_kind), public, parameter :: dT_y    = 60.0D0
  real (kind=real_kind), public, parameter :: dtheta_z= 10.0D0

  public :: hs_v_forcing
  public :: hs_T_forcing
  public :: hs0_init_state
  public :: hs_forcing

contains

  subroutine hs_forcing(elemin,hvcoord,nm1,nm1_Q,dt)

    type (element_t)                 :: elemin
    type (hvcoord_t)                 :: hvcoord
    integer                          :: nm1,nm1_Q  ! timelevel to use
    real (kind=real_kind)                   :: dt

    ! local
    real (kind=real_kind)                   :: pmid,r0,r1,dtf_q,dp,rdp,FQ
    real (kind=real_kind), dimension(np,np) :: psfrc 
    real (kind=real_kind)                   :: temperature(np,np,nlev)
    real (kind=real_kind)                   :: v(np,np,3,nlev)
    real (kind=real_kind)                   :: fv(np,np,3,nlev)
    integer                                 :: i,j,k,q

    ! tms
    real (kind=real_kind) :: pi(np,np,nlev)
    real (kind=real_kind) :: zm(np,np,nlev)
    real (kind=real_kind) :: phi(np,np,nlev)
    real (kind=real_kind) :: phi_i(np,np,nlevp)
    real (kind=real_kind) :: ksrf(np,np)
    real (kind=real_kind) :: r_exner(np,np,nlev)
    real (kind=real_kind) :: ucomp(np,np,nlev)
    real (kind=real_kind) :: vcomp(np,np,nlev)
    

    dtf_q = dt
    call get_temperature(elemin,temperature,hvcoord,nm1)
        
    do j=1,np
       do i=1,np
          if (sub_case==4) then
             psfrc(i,j) = hvcoord%ps0 * exp ( -elemin%state%phis(i,j)/(Rgas*TREF)) 
          else
             psfrc(i,j) = (elemin%state%ps_v(i,j,nm1))
          endif
       end do
    end do

    elemin%derived%FT(:,:,:) = elemin%derived%FT(:,:,:) + &
         hs_T_forcing(hvcoord,psfrc(1,1),               &
         temperature,elemin%spherep,np, nlev)


    v(:,:,1:2,:) = elemin%state%v(:,:,1:2,:,nm1)
#if ( defined MODEL_THETA_L ) 
    v(:,:,3,:) = elemin%state%w_i(:,:,1:nlev,nm1)  ! dont apply at surface
#else
    v(:,:,3,:) = 0
#endif

    fv = hs_v_forcing(hvcoord,psfrc(1,1),v,np,nlev)

#if ( defined MODEL_THETA_L ) 
    elemin%derived%FM(:,:,1:3,:) = elemin%derived%FM(:,:,1:3,:) + fv(:,:,1:3,:)
#else
    elemin%derived%FM(:,:,1:2,:) = elemin%derived%FM(:,:,1:2,:) + fv(:,:,1:2,:)
#endif


#if 0
Notes on TMS:
compute_tms() output is ksrf, in units of kg/s/m2
ksrf is of the form rho*tau*|V|
in CAMs vertial diffusion, it is applied implicitly, but if it was applied explictly,
it would be via:
           u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol)*taux(:ncol)
           v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol)*tauy(:ncol)
           taux = ksrf*u
           tauy = ksrf*v
           tmp1 = ztod*g/pdel
and this is consistent with units of kg/s/m2 for tau
              
clubb scales the output of compute_tms by: ksrf/rho_ds_zm(1)
  dz_g(k) = state1%zi(i,k)-state1%zi(i,k+1)  
  rho(i,k+1)           = invrs_gravit*state1%pdel(i,pver-k+1)/dz_g(pver-k+1)   
  rho_ds_zt(k+1)       = real(rho(i,k+1), kind = core_rknd)                    
  rho_ds_zt(1)       = rho_ds_zt(2)                                           
  rho_ds_zm       = zt2zm_api(rho_ds_zt)
                  = zt2zm(azt) = spine based interpolation to midpoints from interfaces

  CLUBB scalings looks to be   dz*g/pdel

sgh30 ranges from 0 up to about 400, but is mostly around 100 over rough topo
order of magnitude esimates:
! Determine z0m for orography              ! sgh=1          sgh=100    
z0oro = min( z0fac * horo, z0max )           .1          10
! Calculate neutral drag coefficient
cd = ( karman / log( ( zm(i,pver) + z0oro ) / z0oro) )**2
take zm=10                                   .0075        .33
multiply by rho*g/dp3d  = 1/dz =            .00075        .033        
Damping time in seconds:                  1333s            30s       

so very unstable for sgh=100
#endif


#undef USE_TMS
#ifdef USE_TMS
    call get_phi(elemin,phi,phi_i,hvcoord,nm1)
    do k=1,nlev
       pi(:,:,k)    = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*psfrc(:,:)
       ! cam's exner as used by clubb when calling tms
       r_exner(:,:,k) = (psfrc(:,:)/pi(:,:,k))**kappa
       zm(:,:,k) = (phi(:,:,k)-elemin%state%phis(:,:))/g
    end do
    ucomp(:,:,:)=v(:,:,1,:)
    vcomp(:,:,:)=v(:,:,2,:)
    call compute_tms(np*np,nlev,np*np,ucomp,vcomp,temperature,pi,r_exner,zm,elemin%derived%sgh30,ksrf)

    ! scaling used by CAM's vertical_diffusion
    ! if ksrf is in unts of kg/s/m2 (as claimed), then this scales version
    ! will be in units of 1/s, as desired:
    ksrf(:,:)=ksrf(:,:)*g/elemin%state%dp3d(:,:,nlev,nm1)

#if 0
    if (maxval(ksrf)>0) then
       if ( maxval(ksrf) > 0.1d0*dt)  then
          ! stronger then 1/3h
          print *,"TMS limiter ON:  damping time seconds = ", 1/maxval(ksrf(:,:))
       else
          print *,"TMS limiter OFF: damping time seconds = ", 1/maxval(ksrf(:,:))
       endif
    endif
#endif
    ! CFL limiter
    where ( dt*ksrf(:,:) > 0.1d0 ) 
       ksrf(:,:)=0.1d0/dt
    end where
    ! 3h limiter
    where ( ksrf(:,:) > 1d0/(3*3600) ) 
       ksrf(:,:)=1d0/(3*3600)
    end where

    elemin%derived%FM(:,:,1,nlev) = elemin%derived%FM(:,:,1,nlev) - ksrf(:,:)*ucomp(:,:,nlev)
    elemin%derived%FM(:,:,2,nlev) = elemin%derived%FM(:,:,2,nlev) - ksrf(:,:)*vcomp(:,:,nlev)

#endif


    if (qsize>=1) then
       ! HS with tracer  (Galewsky type forcing, with flux of  2.3e-5 kg/m^2/s
       ! MASS in kg/m^2   = < Q dp_in_Pa / g >   
       ! flux in kg/m^2/s = < FQ dp_in_Pa / g >   
       ! We want < FQ dp_in_Pa / g > = 2.3e-5  so:  FQ = 2.3e-5*g/dp_in_Pa 

       ! lowest layer thickness, in Pa
       dp = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(nlev) ) + &
               ( hvcoord%hybi(nlev+1) - hvcoord%hybi(nlev) )*1000*100
       rdp = 1./ dp
       q=1
       do j=1,np
          do i=1,np
             FQ = rdp * g * 2.3E-5 * COS(elemin%spherep(i,j)%lat)**2
             elemin%derived%FQ(i,j,nlev,q) =elemin%derived%FQ(i,j,nlev,q)+FQ
          enddo
       enddo

       do j=1,np
          do i=1,np
             do k=1,nlev
                pmid = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*(elemin%state%ps_v(i,j,nm1))
                r0=elemin%state%Q(i,j,k,q)
                r1=r0
                call Prim_Condense(r1,temperature(i,j,k),pmid)
                elemin%derived%FQ(i,j,k,q) = elemin%derived%FQ(i,j,k,q) + &
                     (r1-r0)/(dtf_q)
             enddo
          enddo
       enddo
    endif

  end subroutine hs_forcing

  function hs_v_forcing(hvcoord,ps,v,npts,nlevels) result(hs_v_frc)

    integer, intent(in)               :: npts
    integer, intent(in)               :: nlevels
    type (hvcoord_t), intent(in)       :: hvcoord
    real (kind=real_kind), intent(in) :: ps(npts,npts)

    real (kind=real_kind), intent(in) :: v(npts,npts,3,nlevels)
    real (kind=real_kind)             :: hs_v_frc(npts,npts,3,nlevels)

    ! Local variables

    integer i,j,k
    real (kind=real_kind) :: k_v
    real (kind=real_kind) :: p,etam

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p    = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             etam      = hvcoord%hyam(k) + hvcoord%hybm(k)
             k_v = k_f*MAX(0.0_real_kind,(etam - sigma_b )/(1.0_real_kind - sigma_b))
             if (sub_case==5) &
                  k_v = 8*k_f*MAX(0.0_real_kind,(etam - sigma_b )/(1.0_real_kind - sigma_b))

             ! surface drag: 8x stronger 
             if (sub_case==6 .and. k==nlevels)  k_v = 8*k_v

             hs_v_frc(i,j,1,k) = -k_v*v(i,j,1,k)
             hs_v_frc(i,j,2,k) = -k_v*v(i,j,2,k)

             etam      = hvcoord%hyai(k) + hvcoord%hybi(k)
             k_v = k_f*MAX(0.0_real_kind,(etam - sigma_b )/(1.0_real_kind - sigma_b))
             hs_v_frc(i,j,3,k) = -k_v*v(i,j,3,k)
          end do
       end do
    end do

  end function hs_v_forcing

  function hs_T_forcing(hvcoord,ps,T,sphere,npts,nlevels) result(hs_T_frc)

    integer, intent(in) :: npts
    integer, intent(in) :: nlevels

    type (hvcoord_t), intent(in)          :: hvcoord
    real (kind=real_kind), intent(in)    :: ps(npts,npts)
    real (kind=real_kind), intent(in)    :: T(npts,npts,nlevels)
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)

    real (kind=real_kind)             :: hs_T_frc(npts,npts,nlevels)

    ! Local variables

    real (kind=real_kind) :: p,logprat,pratk,Teq
    real (kind=real_kind) :: logps0,etam
    real (kind=real_kind) :: lat,snlat

    real (kind=real_kind) :: k_t(npts,npts)
    real (kind=real_kind) :: snlatsq(npts,npts)
    real (kind=real_kind) :: cslatsq(npts,npts)

    real (kind=real_kind) :: rec_one_minus_sigma_b

    integer i,j,k

    logps0    = LOG(hvcoord%ps0 )

    do j=1,npts
       do i=1,npts
         snlat        = SIN(sphere(i,j)%lat)
         snlatsq(i,j) = snlat*snlat
         cslatsq(i,j) = 1.0D0 - snlatsq(i,j)
       end do
    end do

    rec_one_minus_sigma_b = 1.0D0/(1.0D0 - sigma_b)

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p         = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             logprat   = LOG(p)-logps0
             pratk     = EXP(kappa*(logprat))
             etam      = hvcoord%hyam(k) + hvcoord%hybm(k)

             k_t(i,j)= k_a + (k_s-k_a)*cslatsq(i,j)*cslatsq(i,j)* &
                       MAX(0.0D0,(etam - sigma_b)/(1.0D0 - sigma_b))
             Teq     = MAX(200.0D0,(315.0D0 - dT_y*snlatsq(i,j) - dtheta_z*logprat*cslatsq(i,j))*pratk)

#if 0
             ! ======================================
             ! This is a smooth forcing 
             ! for debugging purposes only...
             ! ======================================

             k_t(i,j)= k_a 
             pratk   = EXP(0.081*(logprat))
             Teq     = (315.0D0 - dT_y*snlatsq(i,j))*pratk
#endif
             hs_T_frc(i,j,k)= -k_t(i,j)*(T(i,j,k)-Teq)
             if (sub_case==2 .or. sub_case==5 .or. sub_case==6) &
                  hs_T_frc(i,j,k)= -k_a*(T(i,j,k)-Teq)
             if (sub_case==3) then
                if (etam > sigma_b) hs_T_frc(i,j,k)=0
             endif
          end do
       end do
    end do
      
  end function hs_T_forcing

  subroutine hs0_init_state(elem, hybrid, hvcoord,nets,nete,Tinit)

    type(element_t),        intent(inout) :: elem(:)
    type(hybrid_t),         intent(in)    :: hybrid                   ! hybrid parallel structure
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: nets
    integer,                intent(in)    :: nete
    real (kind=real_kind),  intent(in)    :: Tinit

    ! Local variables
    
    integer ie,i,j,k,q,tl
    integer :: nm1 
    integer :: n0 
    integer :: np1
    real (kind=real_kind) :: lat_mtn,lon_mtn,r_mtn,h_mtn,rsq,lat,lon
    real (kind=real_kind) :: temperature(np,np,nlev),p(np,np),exner(np,np),ps(np,np)

    if (hybrid%masterthread) write(iulog,*) 'initializing Held-Suarez primitive equations test'

    nm1= 1
    n0 = 2
    np1= 3

    do ie=nets,nete

       elem(ie)%state%ps_v(:,:,n0) =hvcoord%ps0
       elem(ie)%state%ps_v(:,:,nm1)=hvcoord%ps0
       elem(ie)%state%ps_v(:,:,np1)=hvcoord%ps0

       elem(ie)%state%v(:,:,:,:,n0) =0.0D0
       elem(ie)%state%v(:,:,:,:,nm1)=elem(ie)%state%v(:,:,:,:,n0)
       elem(ie)%state%v(:,:,:,:,np1)=elem(ie)%state%v(:,:,:,:,n0)

#ifdef MODEL_THETA_L
       elem(ie)%state%w_i = 0.0
#endif

       temperature(:,:,:)=Tinit

       ! if topo file was given in the namelist, PHIS was initilized in prim_main
       ! otherwise assume 0
#ifndef HOMME_WITHOUT_PIOLIBRARY
       if (infilenames(1)=='') then
          elem(ie)%state%phis(:,:)=0.0D0
       endif
#endif

#undef HS_TOPO1
#ifdef HS_TOPO1
       lat_mtn = dd_pi/6
       lon_mtn = 3*dd_pi/2
       r_mtn = dd_pi/9
       h_mtn = 4000
       do i=1,np
          do j=1,np
             lat = elem(ie)%spherev(i,j)%lat
             lon = elem(ie)%spherev(i,j)%lon
             rsq=MIN((lat-lat_mtn)**2 + (lon-lon_mtn)**2,R_mtn**2)
             elem(ie)%state%phis(i,j)=g*h_mtn*(1.0D0 - SQRT(rsq)/R_mtn)
          enddo
       enddo
#endif


       ! initialize surface pressure to be consistent with topo
       elem(ie)%state%ps_v(:,:,n0) = elem(ie)%state%ps_v(:,:,n0)*&
            exp(-elem(ie)%state%phis(:,:) / (Rgas*Tinit))
       elem(ie)%state%ps_v(:,:,nm1)=elem(ie)%state%ps_v(:,:,n0)
       elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,n0)

#if 0
       do k=1,nlev
          p(:,:)=hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,n0)
          exner(:,:) = (p(:,:)/hvcoord%ps0)**kappa
          temperature(:,:,k)=(Tinit-150)+150*exner(:,:)
       enddo
#endif


       if (qsize>=1) then
       q=1
       elem(ie)%state%Q(:,:,:,q) =0  ! moist HS tracer IC=0
       do q=2,qsize
          elem(ie)%state%Q(:,:,:,q) =temperature(:,:,:)/400
       enddo
       endif
       ps=elem(ie)%state%ps_v(:,:,n0)
       call set_thermostate(elem(ie),ps,temperature,hvcoord)

    end do

  end subroutine hs0_init_state


  ! pcols=ncol=np*np
  ! pver=nlev
  ! zm  = (phi()-phis())/g
  subroutine compute_tms( pcols    , pver    , ncol    ,                     &
                          u        , v       , t       , pmid    , exner   , &
                          zm       , sgh     , ksrf    )
  integer,  parameter :: r8 = selected_real_kind(12) ! 8 byte real

  real(r8), parameter :: horomin= 1._r8       ! Minimum value of subgrid orographic height for mountain stress [ m ]
  real(r8), parameter :: z0max  = 100._r8     ! Maximum value of z_0 for orography [ m ]
  real(r8), parameter :: dv2min = 0.01_r8     ! Minimum shear squared [ m2/s2 ]
  real(r8)            :: orocnst              ! Converts from standard deviation to height [ no unit ]
  real(r8)            :: z0fac                ! Factor determining z_0 from orographic standard deviation [ no unit ] 
  real(r8)            :: karman               ! von Karman constant
  real(r8)            :: gravit               ! Acceleration due to gravity
  real(r8)            :: rair                 ! Gas constant for dry air

    !------------------------------------------------------------------------------ !
    ! Turbulent mountain stress parameterization                                    !  
    !                                                                               !
    ! Returns surface drag coefficient and stress associated with subgrid mountains !
    ! For points where the orographic variance is small ( including ocean ),        !
    ! the returned surface drag coefficient and stress is zero.                     !
    !                                                                               !
    ! Lastly arranged : Sungsu Park. Jan. 2010.                                     !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments ! 
    ! ---------------------- !


    integer,  intent(in)  :: pcols                 ! Number of columns (dimension)
    integer,  intent(in)  :: ncol                 ! Number of columns
    integer,  intent(in)  :: pver                  ! Number of model layers

    real(r8), intent(in)  :: u(pcols,pver)         ! Layer mid-point zonal wind [ m/s ]
    real(r8), intent(in)  :: v(pcols,pver)         ! Layer mid-point meridional wind [ m/s ]
    real(r8), intent(in)  :: t(pcols,pver)         ! Layer mid-point temperature [ K ]
    real(r8), intent(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: exner(pcols,pver)     ! Layer mid-point exner function [ no unit ]
    real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
    real(r8), intent(in)  :: sgh(pcols)            ! Standard deviation of orography [ m ]
    
    real(r8), intent(out) :: ksrf(pcols)           ! Surface drag coefficient [ kg/s/m2 ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer  :: i                                  ! Loop index
    integer  :: kb, kt                             ! Bottom and top of source region
    
    real(r8) :: horo                               ! Orographic height [ m ]
    real(r8) :: z0oro                              ! Orographic z0 for momentum [ m ]
    real(r8) :: dv2                                ! (delta v)**2 [ m2/s2 ]
    real(r8) :: ri                                 ! Richardson number [ no unit ]
    real(r8) :: stabfri                            ! Instability function of Richardson number [ no unit ]
    real(r8) :: rho                                ! Density [ kg/m3 ]
    real(r8) :: cd                                 ! Drag coefficient [ no unit ]
    real(r8) :: vmag                               ! Velocity magnitude [ m /s ]

    orocnst  = 1            ! namelist tms_orocnst 
    z0fac    = 0.1d0        ! namelist tms_z0fac  0.1, 0.075
    karman   = 0.4d0        ! karman_in    SHR_CONST_KARMAN
    gravit   = g            !gravit_in
    rair     = Rgas         !rair_in
    

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
    do i=1,ncol

     ! determine subgrid orgraphic height ( mean to peak )
     ! SGH30 ranges from 0..120m in NE30 x6t
       horo = orocnst * sgh(i)

     ! No mountain stress if horo is too small

       if( horo < horomin ) then

           ksrf(i) = 0._r8

       else

         ! Determine z0m for orography

           z0oro = min( z0fac * horo, z0max )   ! ranges from [1,112]*.1 = [.1,11.2]

         ! Calculate neutral drag coefficient

           cd = ( karman / log( ( zm(i,pver) + z0oro ) / z0oro) )**2

         ! Calculate the Richardson number over the lowest 2 layers

           kt  = pver - 1
           kb  = pver
           dv2 = max( ( u(i,kt) - u(i,kb) )**2 + ( v(i,kt) - v(i,kb) )**2, dv2min )

         ! Modification : Below computation of Ri is wrong. Note that 'Exner' function here is
         !                inverse exner function. Here, exner function is not multiplied in
         !                the denominator. Also, we should use moist Ri not dry Ri.
         !                Also, this approach using the two lowest model layers can be potentially
         !                sensitive to the vertical resolution.  
         ! OK. I only modified the part associated with exner function.

           ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
                                / ( ( t(i,kt) * exner(i,kt) + t(i,kb) * exner(i,kb) ) * dv2 )

         ! ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
         !                      / ( ( t(i,kt) + t(i,kb) ) * dv2 )

         ! Calculate the instability function and modify the neutral drag cofficient.
         ! We should probably follow more elegant approach like Louis et al (1982) or Bretherton and Park (2009) 
         ! but for now we use very crude approach : just 1 for ri < 0, 0 for ri > 1, and linear ramping.

           stabfri = max( 0._r8, min( 1._r8, 1._r8 - ri ) )
           cd      = cd * stabfri

         ! Compute density, velocity magnitude and stress using bottom level properties

           ! MT: remove rho scalings
           rho     = pmid(i,pver) / ( rair * t(i,pver) ) 
           vmag    = sqrt( u(i,pver)**2 + v(i,pver)**2 )
           ksrf(i) = rho * cd * vmag 
           !taux(i) = -ksrf(i) * u(i,pver)
           !tauy(i) = -ksrf(i) * v(i,pver)

       end if

    end do
    return
  end subroutine 




end module held_suarez_mod

