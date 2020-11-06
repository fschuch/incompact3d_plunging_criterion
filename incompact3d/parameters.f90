subroutine parameter()

  USE param
  USE variables
  USE complex_geometry
  USE decomp_2d

  implicit none

  integer :: is

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter start'
#endif

  if (nrank==0) then
     print *,'==========================================================='
     print *,'======================Incompact3d=========================='
     print *,'===Copyright (c) 2018 Eric Lamballais and Sylvain Laizet==='
     print *,'===Modified by Felipe Schuch and Ricardo Frantz============'
     print *,'==========================================================='
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif

  allocate(nsc(nphi),uset(nphi),cp(nphi),ri(nphi),group(nphi))

  angle = zero
  beta = zero
  cont_phi = 0
  cp = one
  !datapath = './data/'
  dt = zpone
  fpi2 = four
  ifirst = 1
  iin = 0
  ilag = 0
  ilast = 0
  ilit = 0
  imodulo = 100
  imodulo2 = 100
  iprocessing = 100
  irotation = 0
  isave = 10000
  istret = 0
  itest = 1
  itrip = 0
  ivirt = 0
  izap = 1
  jLES = 0
  nclx1 = 0
  nclxn = 0
  nclxS1 = 0
  nclxSn = 0
  ncly1 = 0
  nclyn = 0
  nclyn = 0
  nclyS1 = 0
  nclySn = 0
  nclz1 = 0
  nclzn = 0
  nclzS1 = 0
  nclzSn = 0
  nobjmax = 1
  noise = zero
  noise1 = zero
  npif = 2
  nraf = 10
  nsc = one
  nscheme = 2
  re = one
  ri = zero
  ro = 99999999._mytype
  save_dmap = 0
  save_dphidx = 0
  save_dphidy = 0
  save_dphidz = 0
  save_dudx = 0
  save_dudy = 0
  save_dudz = 0
  save_dvdx = 0
  save_dvdy = 0
  save_dvdz = 0
  save_dwdx = 0
  save_dwdy = 0
  save_dwdz = 0
  save_ibm = 0
  save_pc = 0
  save_phi = 0
  save_phim = 0
  save_pre = 0
  save_prem = 0
  save_qc = 0
  save_utmap = 0
  save_ux = 0
  save_uxm = 0
  save_uy = 0
  save_uym = 0
  save_uz = 0
  save_uzm = 0
  save_V = 0
  save_w = 0
  save_w1 = 0
  save_w2 = 0
  save_w3 = 0
  u1 = 2
  u2 = one
  uset = zero
  wrotation = zero
  xlx = one
  yly = one
  zlz = one

  call ft_parameter(1)

  if (nclx1.eq.0.and.nclxn.eq.0) then
     nclx=.true.
     nxm=nx
  else
     nclx=.false.
     nxm=nx-1
  endif
  if (ncly1.eq.0.and.nclyn.eq.0) then
     ncly=.true.
     nym=ny
  else
     ncly=.false.
     nym=ny-1
  endif
  if (nclz1.eq.0.and.nclzn.eq.0) then
     nclz=.true.
     nzm=nz
  else
     nclz=.false.
     nzm=nz-1
  endif

  dx=xlx/real(nxm,mytype)
  dy=yly/real(nym,mytype)
  dz=zlz/real(nzm,mytype)

  if (nrank==0) then
     print *,''
     print *,'(lx,ly,lz)=',xlx,yly,zlz
     print *,'(nx,ny,nz)=',nx,ny,nz
     print *,'(dx,dy,dz)=',dx,dy,dz
     print *,'(nx*ny*nz)=',nx*ny*nz
     print *,'(p_row,p_col)=',p_row,p_col
     print *,''
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
     print *,'Numerical precision: Double, saving in single'
#else
     print *,'Numerical precision: Double'
#endif
#else
     print *,'Numerical precision: Single'
#endif
     write(*,"(' Boundary condition : (nclx1 ,nclxn )=(',I1,',',I1,')')") nclx1,nclxn
     write(*,"('                      (ncly1 ,nclyn )=(',I1,',',I1,')')") ncly1,nclyn
     write(*,"('                      (nclz1 ,nclzn )=(',I1,',',I1,')')") nclz1,nclzn
     write(*,"(' High and low speed : u1=',F6.2,' and u2=',F6.2)") u1,u2
     write(*,"(' Reynolds number Re : ',F15.8)") re
     write(*,"(' Time step dt       : ',F15.8)") dt
     write (*,"(' Spatial scheme     : ',F15.8)") fpi2
     if (jLES.eq.0) print *,'                   : DNS'
     if (jLES.eq.1) print *,'                   : iLES'
     if (jLES.eq.2) print *,'                   : Explicit Simple Smagorinsky'
     if (jLES.eq.3) print *,'                   : Explicit Wall-Adaptive LES'
     if (jLES.eq.4) print *,'                   : Explicit Dynamic Smagorinsky LES'
     if (nscheme.eq.1) print *,'Temporal scheme    : Adams-bashforth 2'
     if (nscheme.eq.2) print *,'Temporal scheme    : Adams-bashforth 3'
     if (nscheme.eq.3) print *,'Temporal scheme    : Runge-Kutta 3'
     if (iscalar.eq.0) print *,'Scalar             : off'
     if (iscalar.eq.1) then
        print *,'Scalar             : on'
        write(*,"(' Boundary condition : (nclxS1,nclxSn)=(',I1,',',I1,')')") nclxS1,nclxSn
        write(*,"('                      (nclyS1,nclySn)=(',I1,',',I1,')')") nclyS1,nclySn
        write(*,"('                      (nclzS1,nclzSn)=(',I1,',',I1,')')") nclzS1,nclzSn
        do is=1, nphi
           write (*,"(' Particle fraction  : #',I1)") is
           write (*,"(' Concentration      : ',F15.8)") cp(is)
           write (*,"(' Richardson number  : ',F15.8)") ri(is)
           write (*,"(' Settling velocity  : ',F15.8)") uset(is)
           write (*,"(' Schmidt number     : ',F15.8)") nsc(is)
        end do
     endif
     if (ivirt.eq.0) print *,'Immersed boundary  : off'
     if (ivirt.eq.1) print *,'Immersed boundary  : old school'
     if (ivirt.eq.2) print *,'Immersed boundary  : on with Lagrangian Poly'
     if (angle.ne.0.) write(*,"(' Solid rotation     : ',F6.2)") angle
     print *,''
  endif

  xnu=one/re

#ifdef DOUBLE_PREC
  anglex = dsin(pi*angle/180._mytype)
  angley = dcos(pi*angle/180._mytype)
#else
  anglex = sin(pi*angle/180._mytype)
  angley = cos(pi*angle/180._mytype)
#endif

  dx2=dx*dx
  dy2=dy*dy
  dz2=dz*dz

  call system('mkdir '//trim(datapath)//' 2> /dev/null')
  call system('mkdir '//trim(datapath)//'/3d-snapshot 2> /dev/null')
  call system('mkdir '//trim(datapath)//'/xy-plane 2> /dev/null')
  call system('mkdir '//trim(datapath)//'/xz-plane 2> /dev/null')
  call system('mkdir '//trim(datapath)//'/yz-plane 2> /dev/null')
  !call system('mkdir '//trim(datapath)//'/mean 2> /dev/null')

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter done'
#endif

  return
end subroutine parameter
