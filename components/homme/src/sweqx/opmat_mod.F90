#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module opmat_mod
  ! ===========================================================================
  !
  ! Extract the matrices of two linear operators of sweqx, for eigenvalue
  ! analysis (see test/sw_conservative/python/opmat_eigs.py):
  !
  !   L : the swtc1 transport operator, dp/dt = L p, exactly as computed by
  !       advance_mod:compute_and_apply_rhs for the prescribed velocity
  !       (analytic or DSS(curl psi), see shallow_water_mod:tc1_init_velocity)
  !       in the namelist theta_advect_form, with pmean=0.  Includes the DSS
  !       and the inverse mass matrix.  L is the tendency (L p = -div(u p) for
  !       the conservation form).
  !
  !   B : the hyperviscosity operator as applied in advance_hypervis with
  !       hypervis_order=2:  p -> p - dt*nu_s*B p,  B p = rspheremp*DSS(biharmonic_wk(p)).
  !       B approximates del^4 (positive semidefinite); it excludes dt and nu_s.
  !
  ! Column j is computed by applying the operator to the continuous delta
  ! function at unique GLL node j (p=1 at every element copy of node j, 0
  ! elsewhere).  The nlev levels are independent, so nlev probes are applied
  ! with each operator application.  By default many columns are probed at once
  ! (coloring): the columns of the same local GLL node of elements of one color,
  ! with elements of the same color on the same cube face and >= 4 (L) or 6 (B)
  ! elements apart, so that the supports of their columns (the nodes within 1 (L)
  ! or 2 (B) elements of the node) are disjoint.  The support of each probed
  ! column is labeled with DSS operations, which also check that supports do
  ! not overlap.  This needs 6*16*np^2/nlev (L) and 6*36*np^2/nlev (B)
  ! operator applications independent of ne, instead of n/nlev.
  !
  ! Enabled with environment variables (so that normal runs and the shared
  ! namelist code are unaffected):
  !   SWEQX_OPMAT=<prefix>        extract, write <prefix>.*, then stop.
  !   SWEQX_OPMAT_OPS=L|B|LB      operators to extract (default LB)
  !   SWEQX_OPMAT_COLOR=0         probe one column per level at a time (no coloring)
  ! Requires integration="runge_kutta", 1 thread, kmass=-1, and (for L)
  ! test_case="swtc1".  Any number of MPI ranks.
  !
  ! Output.  dof = 0-based contiguous global index of the unique GLL nodes,
  ! 0 <= dof < n = 6 ne^2 (np-1)^2 + 2, in HOMME's unique-point ordering
  ! (elem%idxP%UniquePtOffset, dof_mod:CreateMetaData), independent of the
  ! MPI decomposition.  Binary files are stream (no record markers), native
  ! (little) endian, and each MPI rank r writes the rows it owns:
  !   <prefix>.meta              text, "key = value"
  !   <prefix>.dofs.r<r>.bin     (int32 dof, int64 gdofP, real64 lat, lon, mass)
  !                              mass = 1/rspheremp, the assembled GLL mass (diagonal of
  !                              the mass matrix M; elem%spheremp is the element's contribution)
  !   <prefix>.L.r<r>.bin        (int32 row, int32 col, real64 value), exact zeros skipped
  !   <prefix>.B.r<r>.bin        same
  !   <prefix>.Ltest.r<r>.bin    (int32 dof, real64 x(nlev), real64 y(nlev)), y = L x computed
  !   <prefix>.Btest.r<r>.bin    directly, for test fields x, to validate the extracted matrix
  !
  ! ===========================================================================
  use kinds, only : real_kind, int_kind, long_kind
  use dimensions_mod, only : np, nlev, nelemd, ne
  use element_mod, only : element_t
  use edgetype_mod, only : EdgeBuffer_t
  use derivative_mod, only : derivative_t
  use hybrid_mod, only : hybrid_t
  use time_mod, only : timelevel_t
  use parallel_mod, only : abortmp
#ifdef _MPI
  use parallel_mod, only : MPI_MAX, MPI_SUM, MPI_IN_PLACE, MPIinteger_t, MPIreal_t
#endif
  implicit none
  private

  public :: opmat_requested
  public :: opmat_extract

  type :: triplet_t
     sequence
     integer(kind=4) :: row, col
     real(kind=8)    :: val
  end type triplet_t

  integer, parameter :: bufsize = 65536

contains

  ! true if the environment variable SWEQX_OPMAT is set (nonempty)
  logical function opmat_requested()
    character(len=1024) :: prefix
    integer :: stat
    call get_environment_variable("SWEQX_OPMAT",prefix,status=stat)
    opmat_requested = (stat==0 .and. len_trim(prefix)>0)
  end function opmat_requested


  subroutine opmat_extract(elem,edge3,deriv,hybrid,tl,nets,nete)
    use dof_mod, only : CreateMetaData
    use control_mod, only : test_case, sub_case, theta_advect_form, kmass, hypervis_order, &
         hypervis_scaling, topology
    use physical_constants, only : dd_pi
    use shallow_water_mod, only : tc1_parameters
    use advance_mod, only : compute_and_apply_rhs
    use viscosity_mod, only : biharmonic_wk
    use edge_mod, only : edgeVpack, edgeVunpack
    use bndry_mod, only : bndry_exchangeV

    type (element_t)     , intent(inout), target :: elem(:)
    type (EdgeBuffer_t)  , intent(inout) :: edge3
    type (derivative_t)  , intent(in) :: deriv
    type (hybrid_t)      , intent(in) :: hybrid
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets,nete

    character(len=1024) :: prefix
    character(len=8)    :: ops
    character(len=16)   :: rankstr
    integer(kind=int_kind) :: fdof(np,np,nelemd)   ! 1-based contiguous global id of each GLL node
    real (kind=real_kind) :: ptens(np,np,nlev,nets:nete), vtens(np,np,2,nlev,nets:nete)
    real (kind=real_kind) :: y(np,np,nlev,nets:nete)
    real (kind=real_kind), allocatable :: xg(:,:)
    real (kind=real_kind) :: u0,alpha,rearth,zero,t0,t1
    logical :: curlpsi
    integer :: n,nloc,stat,ie,i,j,ii,k,iop,ierr,unit
    character(len=1) :: op
    integer :: ncall,ic,col,row
    integer(kind=8) :: nnz
    type (triplet_t), allocatable :: buf(:)
    integer :: nbuf
    logical :: color, owned(np,np,nelemd)
    character(len=8) :: colorstr
    integer :: eix(nelemd),eiy(nelemd),ecolor(nelemd),kc,nring,nset,iset,c,a,b
    real (kind=real_kind) :: lab(np,np,nlev,nets:nete), dxe

    call get_environment_variable("SWEQX_OPMAT",prefix,status=stat)
    call get_environment_variable("SWEQX_OPMAT_OPS",ops,status=stat)
    if (stat/=0 .or. len_trim(ops)==0) ops = "LB"
    if (hybrid%hthreads /= 1) call abortmp("SWEQX_OPMAT: requires 1 thread")
    if (nets/=1 .or. nete/=nelemd) call abortmp("SWEQX_OPMAT: requires nets=1,nete=nelemd")
    if (kmass /= -1) call abortmp("SWEQX_OPMAT: requires kmass=-1 (independent levels)")
    if (index(ops,"L")>0 .and. test_case/="swtc1") &
         call abortmp("SWEQX_OPMAT: L requires test_case=swtc1 (prescribed steady velocity)")
    if (index(ops,"B")>0 .and. hypervis_order/=2) &
         call abortmp("SWEQX_OPMAT: B requires hypervis_order=2")

    call get_environment_variable("SWEQX_OPMAT_COLOR",colorstr,status=stat)
    color = .not. (stat==0 .and. trim(colorstr)=="0")

    write(rankstr,'(a,i4.4,a)') '.r',hybrid%par%rank,'.bin'

    ! contiguous global ids of the unique points, on every element copy
    call CreateMetaData(hybrid%par,elem,fdofp=fdof)
    nloc = maxval(fdof)
    n = nloc
#ifdef _MPI
    call MPI_Allreduce(nloc,n,1,MPIinteger_t,MPI_MAX,hybrid%par%comm,ierr)
#endif
    call tc1_parameters(u0,alpha,rearth,curlpsi)

    owned = .false.
    do ie=nets,nete
       do ii=1,elem(ie)%idxP%NumUniquePts
          owned(elem(ie)%idxP%ia(ii),elem(ie)%idxP%ja(ii),ie) = .true.
       enddo
    enddo
    if (color) then
       ! element (ix,iy) index on its cube face, from the equiangular corner coordinates
       if (topology/="cube") call abortmp("SWEQX_OPMAT: coloring requires topology=cube, set SWEQX_OPMAT_COLOR=0")
       dxe = (dd_pi/2)/ne
       do ie=nets,nete
          eix(ie) = floor( (sum(elem(ie)%corners(:)%x)/4 + dd_pi/4)/dxe )
          eiy(ie) = floor( (sum(elem(ie)%corners(:)%y)/4 + dd_pi/4)/dxe )
          if (min(eix(ie),eiy(ie))<0 .or. max(eix(ie),eiy(ie))>=ne .or. &
               elem(ie)%FaceNum<1 .or. elem(ie)%FaceNum>6) &
               call abortmp("SWEQX_OPMAT: unexpected element coordinates, set SWEQX_OPMAT_COLOR=0")
       enddo
    endif

    if (hybrid%par%masterproc) then
       print *,'SWEQX_OPMAT: extracting operators ',trim(ops),' to ',trim(prefix),'.*'
       print *,'SWEQX_OPMAT: n=',n,' expected 6*ne^2*(np-1)^2+2=',6*ne*ne*(np-1)**2+2
       open(newunit=unit,file=trim(prefix)//'.meta',status='replace',action='write')
       write(unit,'(a,i0)') 'ne = ',ne
       write(unit,'(a,i0)') 'np = ',np
       write(unit,'(a,i0)') 'nlev = ',nlev
       write(unit,'(a,i0)') 'n = ',n
       write(unit,'(a,i0)') 'nranks = ',hybrid%par%nprocs
       write(unit,'(a,a)')  'ops = ',trim(ops)
       write(unit,'(a,a)')  'test_case = ',trim(test_case)
       write(unit,'(a,i0)') 'sub_case = ',sub_case
       write(unit,'(a,es25.17)') 'alpha = ',alpha
       write(unit,'(a,es25.17)') 'u0 = ',u0
       write(unit,'(a,es25.17)') 'rearth = ',rearth
       write(unit,'(a,i0)') 'theta_advect_form = ',theta_advect_form
       if (curlpsi) then
          write(unit,'(a)') 'velocity = curlpsi'
       else
          write(unit,'(a)') 'velocity = analytic'
       endif
       write(unit,'(a,es25.17)') 'hypervis_scaling = ',hypervis_scaling
       write(unit,'(a,l1)') 'color = ',color
       close(unit)
    endif

    ! coordinates and assembled mass of the unique points owned by this rank
    open(newunit=unit,file=trim(prefix)//'.dofs'//trim(rankstr),access='stream', &
         form='unformatted',status='replace',action='write')
    do ie=nets,nete
       do ii=1,elem(ie)%idxP%NumUniquePts
          i=elem(ie)%idxP%ia(ii); j=elem(ie)%idxP%ja(ii)
          write(unit) int(fdof(i,j,ie)-1,4), int(elem(ie)%gdofP(i,j),8), &
               real(elem(ie)%spherep(i,j)%lat,8), real(elem(ie)%spherep(i,j)%lon,8), &
               real(1/elem(ie)%rspheremp(i,j),8)
       enddo
    enddo
    close(unit)

    ! test fields (continuous: one value per global dof), for validation
    allocate(xg(n,nlev))
    call test_fields(xg)

    allocate(buf(bufsize))
    zero = 0

    do iop=1,2
       op = 'LB'(iop:iop)
       if (index(ops,op)==0) cycle
       call cpu_time(t0)

       ! validation: y = op(x) for the test fields
       do ie=nets,nete
          do k=1,nlev
             do j=1,np
                do i=1,np
                   elem(ie)%state%p(i,j,k,tl%n0) = xg(fdof(i,j,ie),k)
                enddo
             enddo
          enddo
       enddo
       call apply_op(op)
       open(newunit=unit,file=trim(prefix)//'.'//op//'test'//trim(rankstr),access='stream', &
            form='unformatted',status='replace',action='write')
       do ie=nets,nete
          do ii=1,elem(ie)%idxP%NumUniquePts
             i=elem(ie)%idxP%ia(ii); j=elem(ie)%idxP%ja(ii)
             write(unit) int(fdof(i,j,ie)-1,4), real(xg(fdof(i,j,ie),:),8), real(y(i,j,:,ie),8)
          enddo
       enddo
       close(unit)

       ! the matrix.  Probe sets: without coloring, set s = column s+1.  With coloring,
       ! set s = (element color c, local GLL index (a,b)): the columns of the nodes (a,b)
       ! owned by the elements of color c, whose supports are disjoint.  nlev sets per
       ! operator application.  lab = column whose support contains each node (0: none)
       if (color) then
          kc = merge(4,6,op=='L')      ! same-color elements are >= kc elements apart
          nring = merge(1,2,op=='L')   ! support of a column: nodes within nring elements
          do ie=nets,nete
             ecolor(ie) = (elem(ie)%FaceNum-1)*kc*kc + mod(eix(ie),kc)*kc + mod(eiy(ie),kc)
          enddo
          nset = 6*kc*kc*np*np
       else
          nset = n
       endif
       ncall = (nset+nlev-1)/nlev
       open(newunit=unit,file=trim(prefix)//'.'//op//trim(rankstr),access='stream', &
            form='unformatted',status='replace',action='write')
       nbuf=0; nnz=0
       do ic=0,ncall-1
          if (color) then
             lab = 0
             do k=1,nlev
                iset = ic*nlev+k-1
                if (iset>=nset) cycle
                c = iset/(np*np); a = mod(iset,np)+1; b = mod(iset/np,np)+1
                do ie=nets,nete
                   if (ecolor(ie)==c .and. owned(a,b,ie)) lab(a,b,k,ie) = fdof(a,b,ie)
                enddo
             enddo
             call dss_labels(.false.)       ! label every element copy of the probed nodes
             do ie=nets,nete
                elem(ie)%state%p(:,:,:,tl%n0) = merge(1d0,0d0,lab(:,:,:,ie)>0)
             enddo
             do i=1,nring
                call dss_labels(.true.)     ! extend labels to the support
             enddo
          else
             do ie=nets,nete
                do k=1,nlev
                   col = ic*nlev+k
                   elem(ie)%state%p(:,:,k,tl%n0) = merge(1d0,0d0,fdof(:,:,ie)==col)
                   lab(:,:,k,ie) = merge(col,0,col<=n)
                enddo
             enddo
          endif
          call apply_op(op)
          do ie=nets,nete
             do ii=1,elem(ie)%idxP%NumUniquePts
                i=elem(ie)%idxP%ia(ii); j=elem(ie)%idxP%ja(ii)
                row = fdof(i,j,ie)
                do k=1,nlev
                   if (y(i,j,k,ie)/=0) then
                      col = nint(lab(i,j,k,ie))
                      if (col<1) call abortmp("SWEQX_OPMAT: nonzero outside the predicted support")
                      nbuf=nbuf+1
                      buf(nbuf)%row = row-1
                      buf(nbuf)%col = col-1
                      buf(nbuf)%val = y(i,j,k,ie)
                      if (nbuf==bufsize) then
                         write(unit) buf(1:nbuf)
                         nnz=nnz+nbuf; nbuf=0
                      endif
                   endif
                enddo
             enddo
          enddo
          if (hybrid%par%masterproc .and. ncall>=2000 .and. mod(ic+1,max(1,ncall/10))==0) then
             call cpu_time(t1)
             print '(a,a,a,i7,a,i7,a,f9.1,a)','SWEQX_OPMAT: ',op,' call ',ic+1,' of ',ncall, &
                  ' (',t1-t0,' s)'
          endif
       enddo
       if (nbuf>0) write(unit) buf(1:nbuf)
       nnz=nnz+nbuf
       close(unit)
       call cpu_time(t1)
       print '(a,a,a,i6,a,i10,a,i7,a,f9.1,a)','SWEQX_OPMAT: ',op,' rank ',hybrid%par%rank, &
            ': nnz=',nnz,'  operator applications=',ncall,'  time=',t1-t0,' s'
    enddo

  contains

    ! y = op(p(n0)).  Overwrites the state velocity and, for B, ps.
    subroutine apply_op(op)
      character(len=1), intent(in) :: op
      integer :: ie
      if (op=='L') then
         ! dt2=0: return the DSS'd tendency in ptens; pmean=0, real_time=0
         call compute_and_apply_rhs(tl%np1,tl%nm1,tl%n0,zero,zero,edge3,elem,zero,hybrid, &
              deriv,vtens,ptens,nets,nete)
         y = ptens
      else
         ! as in advance_hypervis (hypervis_order=2): weak biharmonic, DSS, inverse mass.
         ! biharmonic_wk acts on p+ps and on the lat-lon velocity
         do ie=nets,nete
            elem(ie)%state%ps = 0
            elem(ie)%state%v(:,:,:,:,tl%n0) = 0
         enddo
         call biharmonic_wk(elem,ptens,vtens,deriv,edge3,hybrid,tl%n0,nets,nete)
         do ie=nets,nete
            call edgeVpack(edge3,ptens(1,1,1,ie),nlev,0,ie)
         enddo
         call bndry_exchangeV(hybrid,edge3)
         do ie=nets,nete
            call edgeVunpack(edge3,ptens(1,1,1,ie),nlev,0,ie)
            do k=1,nlev
               y(:,:,k,ie) = elem(ie)%rspheremp(:,:)*ptens(:,:,k,ie)
            enddo
         enddo
      endif
    end subroutine apply_op

    ! DSS of the labels lab: every element copy of a labeled node gets the label.
    ! With extend, first label all nodes of each element containing a labeled node.
    ! Aborts if two different labels meet (supports of probed columns overlap).
    subroutine dss_labels(extend)
      logical, intent(in) :: extend
      real (kind=real_kind) :: w(np,np,3*nlev,nets:nete), lmax, lmin
      integer :: ie,k
      do ie=nets,nete
         do k=1,nlev
            if (extend .and. any(lab(:,:,k,ie)>0)) then
               lmax = maxval(lab(:,:,k,ie))
               lmin = minval(lab(:,:,k,ie),mask=lab(:,:,k,ie)>0)
               if (lmax/=lmin) call abortmp("SWEQX_OPMAT: coloring conflict in an element")
               lab(:,:,k,ie) = lmax
            endif
            w(:,:,k,ie) = lab(:,:,k,ie)
            w(:,:,nlev+k,ie) = merge(1d0,0d0,lab(:,:,k,ie)>0)
            w(:,:,2*nlev+k,ie) = lab(:,:,k,ie)**2
         enddo
         call edgeVpack(edge3,w(1,1,1,ie),3*nlev,0,ie)
      enddo
      call bndry_exchangeV(hybrid,edge3)
      do ie=nets,nete
         call edgeVunpack(edge3,w(1,1,1,ie),3*nlev,0,ie)
         do k=1,nlev
            where (w(:,:,nlev+k,ie)>0)
               lab(:,:,k,ie) = w(:,:,k,ie)/w(:,:,nlev+k,ie)
            elsewhere
               lab(:,:,k,ie) = 0
            end where
            ! all copies carry the same label iff sum(l^2) = count*mean^2
            if (any(w(:,:,2*nlev+k,ie) /= w(:,:,nlev+k,ie)*lab(:,:,k,ie)**2)) &
                 call abortmp("SWEQX_OPMAT: coloring conflict at a shared node")
         enddo
      enddo
    end subroutine dss_labels

    ! continuous test fields: smooth (low-degree polynomials in x,y,z on the
    ! sphere, cosine bell) and pseudo-random (hash of the dof id)
    subroutine test_fields(xg)
      real (kind=real_kind), intent(out) :: xg(:,:)
      real (kind=real_kind) :: lat,lon,xx,yy,zz,r
      integer :: ie,ii,i,j,d,k
      xg = 0
      do ie=nets,nete
         do ii=1,elem(ie)%idxP%NumUniquePts
            i=elem(ie)%idxP%ia(ii); j=elem(ie)%idxP%ja(ii)
            d = fdof(i,j,ie)
            lat = elem(ie)%spherep(i,j)%lat; lon = elem(ie)%spherep(i,j)%lon
            xx = cos(lat)*cos(lon); yy = cos(lat)*sin(lon); zz = sin(lat)
            do k=1,nlev
               select case (mod(k-1,4))
               case (0)
                  xg(d,k) = 1 + xx + 0.5d0*yy*zz - 0.3d0*xx*xx*yy + 0.2d0*zz**4
               case (1)
                  ! cosine bell centered at (lat,lon)=(0.3,1.0), radius 0.6 rad
                  r = acos(max(-1d0,min(1d0,sin(0.3d0)*zz + cos(0.3d0)*cos(lat)*cos(lon-1d0))))
                  xg(d,k) = merge(0.5d0*(1+cos(acos(-1d0)*r/0.6d0)),0d0,r<0.6d0)
               case (2)
                  xg(d,k) = exp(-4*((xx-0.5d0)**2 + (yy+0.2d0)**2 + (zz-0.4d0)**2))
               case (3)
                  xg(d,k) = sin(real(d,real_kind)*12.9898d0 + k)*43758.5453d0
                  xg(d,k) = xg(d,k) - floor(xg(d,k)) - 0.5d0
               end select
            enddo
         enddo
      enddo
#ifdef _MPI
      ! each dof is owned by exactly one rank: sum = broadcast
      call MPI_Allreduce(MPI_IN_PLACE,xg,size(xg),MPIreal_t,MPI_SUM,hybrid%par%comm,ierr)
#endif
    end subroutine test_fields

  end subroutine opmat_extract

end module opmat_mod
