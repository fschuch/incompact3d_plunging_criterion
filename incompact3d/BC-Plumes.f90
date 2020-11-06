module flow_type
  use decomp_2d
  !Geometry
  real(mytype) :: x0ramp, y0ramp, layer, xlx_pi, xlx_pf, declramp
  real(mytype),dimension(:),allocatable :: ramp
  !Postprocessing
  integer :: nx_pi, nx_pf
  integer :: initstats1
  !Sponge Zone
  type(decomp_info) :: sponge_tyzplane
  real(mytype),dimension(:,:,:),allocatable :: sponge_ux1
  real(mytype),dimension(:,:,:,:),allocatable :: sponge_phi1
  !
  integer :: sponge_type, sponge_type2, sponge_init
  integer :: sponge_nt, nx_xref
  real(mytype) :: sponge_xloc, sponge_xref
  real(mytype),dimension(:),allocatable :: abs_coef_x, abs_coef2_x
  real(mytype),dimension(:),allocatable :: sponge_ux
  real(mytype),dimension(:,:),allocatable :: sponge_phi
  !Other
  integer,dimension(:),allocatable :: phi_type
  real(mytype),dimension(:),allocatable :: phi_clip
  !
  real(mytype),dimension(:,:),allocatable :: area_yz1, area_yzn
  !
  integer :: icrfile, file1, filen, ssfile
  integer :: comp_visu, comp_post
  logical :: flag_probes = .False.
  !
end module flow_type

#include "cla_kinds.f90"
#include "cla.f90"

subroutine ft_parameter(arg)

  use kinds
  use cla !http://web.cecs.pdx.edu/~zaron/pub/CLA-F90.html
  USE param
  USE variables
  USE flow_type
  USE complex_geometry
  USE decomp_2d

  implicit none

  integer,intent(in) :: arg
  integer :: i,j,jj,k,is
  real(mytype) :: x,y,ym,tmp
  character(len=100) :: buffer, label
  integer :: pos, ios, line
  logical :: flag

  !Recomended:
  !arg = 0 for the first call, i.e., inicialization and define nx, ny, nz, nphi, p_row and p_col
  !arg = 1 to read the all the other parameters
  !arg = 2 to re-read parameters on the go

  if (arg.eq.0) then
    iscalar = 1
    call cla_init
    !call cla_register('-p','--plot','Plot with a python file on the go', cla_char,'visu.py') !not working
    call cla_register('-i','--input','Input parameter file', cla_char, 'BC-Plumes.prm')
    call cla_register('-o','--output','Output diretory', cla_char, './data/')
    !call cla_register('-r','--restart','Automatically restart from last checkpoint', cla_flag,'f') !not working
    call cla_register('-0','--inicialization','Set ilast to zero ', cla_flag,'f')
    call cla_register('-w','--2decomp-prow','2D pencil decomposition by 2DECOMP&FFT', cla_int, '0')
    call cla_register('-l','--2decomp-pcol','2D pencil decomposition by 2DECOMP&FFT', cla_int, '0')
    call cla_validate("incompact3d")
    !
    call cla_get('-i',filepath)
    call cla_get('-o',datapath)
    call cla_get('-w',p_row)
    call cla_get('-l',p_col)
  endif

  if (arg.eq.1) then !Defining default values
    declramp=zero
    sponge_type = 0
    sponge_type2 = 0
    allocate(phi_type(nphi))
    phi_type = 1
    initstats1 = 10000000
    sponge_init = 100000000 !subroutine front2d
    sponge_nt = 1
    sponge_xloc = three
    sponge_xref = three*ten*ten

    !velocity
    nclx1 = 2 !Boundary condition in x=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
    nclxn = 2 !Boundary condition in x=Lx (0: Periodic, 1:Free-slip, 2: Dirichlet)
    ncly1 = 2 !Boundary condition in y=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
    nclyn = 1 !Boundary condition in y=Ly (0: Periodic, 1:Free-slip, 2: Dirichlet)
    nclz1 = 1 !Boundary condition in z=0  (0: Periodic, 1:Free-slip, 2: Dirichlet)
    nclzn = 1 !Boundary condition in z=Lz (0: Periodic, 1:Free-slip, 2: Dirichlet)
    !scalar
    nclxS1 = 2 !Boundary condition in x=0  (0: Periodic, 1:No-flux, 2: Dirichlet)
    nclxSn = 2 !Boundary condition in x=Lx (0: Periodic, 1:No-flux, 2: Dirichlet)
    nclyS1 = 2 !Boundary condition in y=0  (0: Periodic, 1:No-flux, 2: Dirichlet)
    nclySn = 1 !Boundary condition in y=Ly (0: Periodic, 1:No-flux, 2: Dirichlet)
    nclzS1 = 1 !Boundary condition in z=0  (0: Periodic, 1:No-flux, 2: Dirichlet)
    nclzSn = 1 !Boundary condition in z=Lz (0: Periodic, 1:No-flux, 2: Dirichlet)
  endif

  open(10,file=filepath,status='old',form='formatted')
  ios=0; line=0
  do while (ios == 0)
    read(10, '(A)', iostat=ios) buffer
    if (ios == 0) then
      line = line + 1

      ! Find the first instance of whitespace.  Split label and data.
      pos = scan(buffer, '    ')
      label = buffer(1:pos)
      buffer = buffer(pos+1:)
!sponge_init1, sponge_init2
      select case (label)
      case ('declramp')
         if (arg.eq.1) read(buffer, *, iostat=ios) declramp
      case ('initstats1')
         if (arg.eq.1) read(buffer, *, iostat=ios) initstats1
      case ('layer')
         if (arg.eq.1) read(buffer, *, iostat=ios) layer
      case ('phi_type')
         ! (1: particle-laden, 2: ambient fluid, 3: Numerical dye)
         if (arg.eq.1) read(buffer, *, iostat=ios) phi_type
      case ('sponge_init')
         if (arg.eq.1) read(buffer, *, iostat=ios) sponge_init
      case ('sponge_nt')
         if (arg.eq.1) read(buffer, *, iostat=ios) sponge_nt
      case ('sponge_type')
         if (arg.eq.1) read(buffer, *, iostat=ios) sponge_type
      case ('sponge_type2')
         if (arg.eq.1) read(buffer, *, iostat=ios) sponge_type2
      case ('sponge_xloc')
         if (arg.eq.1) read(buffer, *, iostat=ios) sponge_xloc
      case ('sponge_xref')
         if (arg.eq.1) read(buffer, *, iostat=ios) sponge_xref
      case ('x0ramp')
         if (arg.eq.1) read(buffer, *, iostat=ios) x0ramp
      case ('xlx_pf')
         if (arg.eq.1) read(buffer, *, iostat=ios) xlx_pf
      case ('xlx_pi')
         if (arg.eq.1) read(buffer, *, iostat=ios) xlx_pi
#include "parameters_input.f90"
      end select
    endif
  enddo
  close(10)

  if (arg.eq.1) then

    y0ramp = yly - one

    call cla_get('-0',flag)
    if (flag) ilast = 0

    !!$  call cla_get('-r',flag)
    !!$  if (flag) then
    !!$        ilit = 1
    !!$        call system('ls -t sauve* | head -1 | tail -c 8 > .temp.out')
    !!$        open(10,file='.temp.out',status='unknown',form='formatted')
    !!$        read(10,*) ifirst
    !!$        ifirst = ifirst + 1
    !!$        close(10,status='delete')
    !!$  endif

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

    !
    !Computatation of the flow rate Inflow/Outflow
    !
    allocate(area_yz1(xsize(2),xsize(3)), area_yzn(xsize(2),xsize(3)))
    area_yz1=dy*dz
    if (.not.ncly) then
      if (xstart(2) .eq. 1 ) area_yz1(1         ,1:xsize(3))=area_yz1(1         ,1:xsize(3))*half
      if (xend(2)   .eq. ny) area_yz1(xsize(2)  ,1:xsize(3))=area_yz1(xsize(2)  ,1:xsize(3))*half
    endif
    if (.not.nclz) then
      if (xstart(3) .eq. 1 ) area_yz1(1:xsize(2),1         )=area_yz1(1:xsize(2),1         )*half
      if (xend(3)   .eq. nz) area_yz1(1:xsize(2),xsize(3)  )=area_yz1(1:xsize(2),xsize(3)  )*half
    endif
    area_yzn = area_yz1
    !Adjusting variables at the bottom before integrate
    jj = int(y0ramp / dy) + 2
    do j = xstart(2), xend(2)
      if (j.lt.jj) then
        area_yz1(j-xstart(2)+1,:) = zero
      elseif (j.eq.jj) then
        tmp = real(j-1,mytype) - y0ramp/dy + half
        area_yz1(j-xstart(2)+1,:) = area_yz1(j-xstart(2)+1,:) * tmp
      endif
    enddo
    jj = int(layer / dy) + 2
    do j = xstart(2), xend(2)
      if (j.lt.jj) then
        area_yzn(j-xstart(2)+1,:) = zero
      elseif (j.eq.jj) then
        tmp = real(j-1,mytype) - layer/dy + half
        area_yzn(j-xstart(2)+1,:) = area_yzn(j-xstart(2)+1,:) * tmp
      endif
    enddo
    !
    !
    !
    nx_pi = int(xlx_pi / dx)+1
    nx_pf = int(xlx_pf / dx)+1
    nx_xref = int(sponge_xref / dx)+1

    allocate(sponge_ux(xsize(2))) !reading velocity profile
    sponge_ux = zero

    allocate(sponge_phi(xsize(2),nphi)) !reading phi profile
    sponge_phi = zero

    if ( sponge_type .eq. 2 ) then !Read from file
      open(10,file=trim(datapath)//'sponge_zone.csv',status='old',form='formatted')
      !write(*,*) 'nrank', nrank, xstart(2)
      do j=1,xstart(2)-1
        read (10,*) sponge_ux(1), sponge_phi(1,:)
        !write(*,*) 'nrank', nrank, j, sponge_ux(1)
      enddo
      do j=1,xsize(2)
        read (10,*) sponge_ux(j), sponge_phi(j,:)
        !write(*,*) 'nrank', nrank, j, sponge_ux(j)
      enddo
      close(10)
      do is=1, nphi
        if ( phi_type(is) .eq. 2) stop 'not prepared for sponge_type = 2'
      enddo
    elseif ( sponge_type .eq. 3 ) then !Schuch et al. (2018) style
      !
    end if

    allocate(ramp(nx))
    ramp = zero
    do i=1,nx
      x = real(i-1,mytype) * dx
      ramp(i) = y0ramp - declramp * (x-x0ramp)
      if (ramp(i) .gt. y0ramp) ramp(i) = y0ramp
      if (ramp(i) .lt. layer)  ramp(i) = layer
    enddo

    allocate(abs_coef_x(1:xsize(1)))
    abs_coef_x = zero
    do i=1,xsize(1)
      x=real( i - 1 ,mytype)*dx
      abs_coef_x(i) = half * (one + tanh(x - xlx + sponge_xloc))
    enddo

    allocate(abs_coef2_x(1:xsize(1)))
    abs_coef2_x = zero
    do i=1,xsize(1)
      x=real( i - 1 ,mytype)*dx
      abs_coef2_x(i) = half * (one + tanh(x - xlx + sponge_xloc - two))
    enddo

    call decomp_info_init(sponge_nt, ny, nz, sponge_tyzplane)

    if ( sponge_type .eq. 1 ) then
      call alloc_x(sponge_ux1, sponge_tyzplane)
      sponge_ux1 = zero
    endif

    if ( sponge_type2 .eq. 1 ) then
      allocate(sponge_phi1(sponge_tyzplane%xsz(1),sponge_tyzplane%xsz(2),sponge_tyzplane%xsz(3),nphi))
      sponge_phi1 = zero
    endif

    allocate(phi_clip(nphi))
    phi_clip(:) = cp(:)

    if (nrank==0) then
      print *,'========================Plumes============================='
      write(*,"(' x0, y0, layer      : (',F8.4,',',F8.4,',',F8.4,')')") x0ramp,&
         y0ramp, layer
      write(*,"(' xlx_pi, xlx_pf     : (',F8.4,',',F8.4,')')") xlx_pi, xlx_pf
      write(*,"(' declramp           : (',F8.4,')')") declramp
      write(*,*) 'Input File         : ', trim(filepath)
      write(*,*) 'Output diretory    : ', trim(datapath)
      write(*,"(' initstats1         : ',I15)") initstats1
      if (sponge_type.eq.0) print *,'Sponge Zone Vel.   : Off'
      if (sponge_type.eq.1) print *,'Sponge Zone Vel.   : Intrinsic value'
      if (sponge_type.eq.2) print *,'Sponge Zone Vel.   : Read from file'
      if (sponge_type.eq.3) print *,'Sponge Zone Vel.   : Schuch et al. (2018) style'
      if (sponge_type2.eq.0) print *,'Sponge Zone Scalar : Off'
      if (sponge_type2.eq.1) print *,'Sponge Zone Scalar : Intrinsic value'
      if (sponge_type2.eq.2) print *,'Sponge Zone Scalar : Read from file'
      if (sponge_type2.eq.3) print *,'Sponge Zone Scalar : Schuch et al. (2018) style'
      if (flag_probes) print *,'Write probes       : On'
      print *,'==========================================================='
    endif
  endif
  !sponge_type, sponge_init, sponge_position
  return
end subroutine ft_parameter
!********************************************************************
subroutine geomcomplex(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
  use flow_type, only : declramp,x0ramp,y0ramp,layer
  use decomp_2d, only : mytype
  use param, only : one,two
  implicit none
  !
  real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
  real(mytype),dimension(ny) :: yp
  integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
  real(mytype)               :: dx,dz
  real(mytype)               :: remp
  integer                    :: i,j,k
  real(mytype)               :: x,y,fy
  real(mytype)               :: zeromach

  zeromach=one
  do while ((one + zeromach / two) .gt. one)
    zeromach = zeromach/two
  end do
  zeromach = 1.0e1*zeromach

  do k=nzi,nzf
    do j=nyi,nyf
      y=yp(j)
      do i=nxi,nxf
        x=real(i - 1,mytype)*dx-x0ramp
        fy = y0ramp-declramp*x
        if (fy .gt. y0ramp) fy=y0ramp
        if (fy .lt. layer)  fy=layer
        if (y .le. fy+zeromach) epsi(i,j,k)=remp
      enddo
    enddo
  enddo
  !
  return
end subroutine geomcomplex
!********************************************************************
subroutine boundary_conditions (ux,uy,uz,phi1,ep1)

  USE param
  USE variables
  USE decomp_2d
  USE flow_type
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2
  !
  real(mytype),dimension(xsize(2),xsize(3)) :: avg_tmp
  !
  real(mytype),dimension(ny) :: sponge_tmp,sponge_tmp1
  !
  real(mytype) :: cx, ct, ct2, uinflow, y, udx, trap, temp, tmp
  integer :: i,ii,j,jj,k,is, code

  if (iin.ne.2) then !if not LST
    !time coefficient for inflow, it provides a smooth start
    ct=one!half+half*tanh(t-three) !se mudar, mudar tambem no initt
    !time coefficient for sponge zone, it provides a smooth start
    !after the front leaves the domain (check front subroutine)
    ct2=half+half*tanh(t-real(sponge_init, mytype)*dt)
  else !if LST
    ct = one
    ct2 = one
  endif

  if (nrank.eq.0) print *, 'ct2', ct2, sponge_init

  !velocity
  !inflow
  if (istret.ne.0) stop 'not prepared for istret.ne.0'
  !Flow rate control, it should be equal to one per width unity
  trap = zero
  do i=1, nx_pi
    call flow_rate(ux(i,:,:), temp, area_yz1)
    trap = trap + temp
  enddo
  !trap = zlz / (trap / real(nx_pi, mytype))
  trap = real(nx_pi, mytype) * zlz / trap
  if (nrank.eq.0) write(*,*) ' trap ux1 =',trap
  do i=1, nx_pi
    ux(i,:,:) = ux(i,:,:) * trap
  enddo
  !Recycling
  ii = nx_pi*2/3-1
  do i=1,5
    ux(i,:,:) = ux(ii+i,:,:)
    uy(i,:,:) = uy(ii+i,:,:)
    uz(i,:,:) = uz(ii+i,:,:)
  enddo
  bxx1(:,:)=ux(1,:,:)
  bxy1(:,:)=uy(1,:,:)
  bxz1(:,:)=uz(1,:,:)

  !Outflow

  udx = one / dx
  if (u1==0) stop 'not prepared for u1 = 0'
  if (u1==1) stop 'not prepared for u1 = 1'
  if (u1==2) cx=u2*gdt(itr)*udx !works better

  if ( sponge_type .eq. 1) then !intrinsic value
    i = mod(itime,sponge_nt)+1
    !
    sponge_ux1(i,:,:) = ux(nx_xref,:,:)
    sponge_phi1(i,:,:,:) = phi1(nx_xref,:,:,:)
    !
    if (ct2 .gt. zero) then
      call mean_plane_x(sponge_ux1, sponge_nt, xsize(2), xsize(3), avg_tmp)
      sponge_tmp = zero
      sponge_tmp1 = zero
      do j=xstart(2),xend(2)
        do k=1,xsize(3)
          sponge_tmp(j) = sponge_tmp(j) + avg_tmp(j-xstart(2)+1,k)
        enddo
      enddo
      call MPI_ALLREDUCE(sponge_tmp,sponge_tmp1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      sponge_tmp1 = sponge_tmp1/real(nz, mytype)
      !
      jj = int(layer / dy) + 2                 !first point above the vertical
      tmp = real(jj-1,mytype) - layer/dy + half!layer. Trapezoidal rule for
      trap = sponge_tmp1(jj)*tmp*dy            !integration is adjusted properly
      do j=jj+1,ny-1
        trap = trap + sponge_tmp1(j)*dy
      enddo
      trap = trap + sponge_tmp1(ny)*half*dy
      if (nrank.eq.0) print *,'trap ux sz = ', trap
      !
      do j=1,xsize(2)
        sponge_ux(j) = sponge_tmp1(j+xstart(2)-1)
      enddo
      trap = one / trap !flow rate per width unit at outlet should be one
      sponge_ux = sponge_ux * trap
    endif
  endif
  !
  if ( sponge_type .eq. 1 .or. sponge_type .eq. 2) then !read from a file
    !Sponge Zone
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=nx_pf,xsize(1)
          ux(i,j,k) = ux(i,j,k)+(sponge_ux(j)-ux(i,j,k))*ct2*abs_coef_x(i)
        enddo
      enddo
    enddo
  elseif ( sponge_type .eq. 3) then !Schuch et al. (2018) style
    !Sponge Zone
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=nx_pf,xsize(1)
          if (ux(i,j,k).lt.0.) then
            ux(i,j,k)= ux(i,j,k)*(one-abs_coef_x(i))
          endif
        enddo
      enddo
    enddo
  end if
  !Convective outflow
  do k=1,xsize(3)
    do j=1,xsize(2)
      bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
      bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
      bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
    enddo
  enddo
  !
  do k=1,xsize(3)
    do j=1,xsize(2)
      bxxn(j,k)=bxxn(j,k)*(one-ep1(nx,j,k))
      bxyn(j,k)=bxyn(j,k)*(one-ep1(nx,j,k))
      bxzn(j,k)=bxzn(j,k)*(one-ep1(nx,j,k))
    enddo
  enddo

  !Scalar
  if (iscalar.eq.1) then
    !Inflow
    do is=1, nphi
      if ( phi_type(is) == 1 ) then !particle-laden
        avg_tmp = zero
        do j=1,xsize(2)
          y=(j+xstart(2)-1-1)*dy
          if (y.ge.yly-one) then
            !Inflow profile
            uinflow = half*(one - tanh(two*re*nsc(is)*uset(is)*(y-yly)))
            do k=1,xsize(3)
              avg_tmp(j,k) = uinflow
            enddo
          endif
        enddo
        call flow_rate(avg_tmp, trap, area_yz1)
        trap = one - trap / zlz
        avg_tmp = avg_tmp + trap
        !
        avg_tmp=cp(is)*avg_tmp
        !for clipping, see subroutine scalar, at convdiff.f90
        !Max value for phi is where y = yly - one
        y = yly - one
        uinflow = half*(one - tanh(two*re*nsc(is)*uset(is)*(y-yly)))
        phi_clip(is) = cp(is)*(uinflow + trap)
#ifdef DEBG
        call flow_rate(avg_tmp, trap, area_yz1)
        if (nrank.eq.0) write(*,*) 'phi trap, is = ', trap, is
#endif
      elseif ( phi_type(is) == 2 ) then !ambient fluid
        avg_tmp=zero
      elseif ( phi_type(is) == 3 ) then !Numerical dye
        avg_tmp=cp(is)
      endif
      ! Recycling technique
      do i=1, nx_pi
        phi1(i,:,:,is) = avg_tmp
      enddo
    enddo

    !BC - Top for scalar
    do is=1, nphi
      call transpose_x_to_y(phi1(:,:,:,is),phi2)
      if (uset(is) .gt. zero) then
        temp = one / (six*re*nsc(is)*dy*uset(is) + eleven)
        do k=1,ysize(3)
          do i=1,ysize(1)
            !Robin on top BC
            phi2(i,ysize(2),k)= temp*(two*phi2(i,ysize(2)-3,k)&
                                      -nine*phi2(i,ysize(2)-2,k)&
                                      +nine*two*phi2(i,ysize(2)-1,k))
          enddo
        enddo
      endif
      call lagpolphiy(phi2)
      call transpose_y_to_x(phi2,phi1(:,:,:,is))
    enddo
    !
    if ( sponge_type2 .eq. 1) then !intrinsic value
      i = mod(itime,sponge_nt)+1
      !
      sponge_phi1(i,:,:,:) = phi1(nx_xref,:,:,:)
      !
      if (ct2 .gt. zero) then
        !
        do is=1,nphi
          call mean_plane_x(sponge_phi1(:,:,:,is), sponge_nt, xsize(2), xsize(3), avg_tmp)
          sponge_tmp = zero
          sponge_tmp1 = zero
          do j=xstart(2),xend(2)
            do k=1,xsize(3)
              sponge_tmp(j) = sponge_tmp(j) + avg_tmp(j-xstart(2)+1,k)
            enddo
          enddo
          call MPI_ALLREDUCE(sponge_tmp,sponge_tmp1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          sponge_tmp1 = sponge_tmp1/real(nz, mytype)
          do j=1,xsize(2)
            sponge_phi(j,is) = sponge_tmp1(j+xstart(2)-1)
          enddo
          !
          if ( phi_type(is) .eq. 1 ) then
            do j=1,xsize(2)
              if (sponge_phi(j,is).lt.0.001) sponge_phi(j,is) = zero
            enddo
          endif
        enddo
      endif
    endif
    !
    if ( sponge_type2 .eq. 1 .or. sponge_type2 .eq. 2) then !read from a file
      !Sponge zone
      do is=1, nphi
        if ( phi_type(is) == 1 ) then !particle-laden
          do k=1,xsize(3)
            do j=1,xsize(2)
              do i=nx_pf,xsize(1)
                phi1(i,j,k,is) = phi1(i,j,k,is)+(sponge_phi(j,is)-phi1(i,j,k,is))*ct2*abs_coef_x(i)
              enddo
            enddo
          enddo
        elseif ( phi_type(is) == 2 ) then !ambient fluid
          do k=1,xsize(3)
            do j=1,xsize(2)
              do i=nx_pf,xsize(1)
                phi1(i,j,k,is) = phi1(i,j,k,is)+(sponge_phi(j,is)-phi1(i,j,k,is))*ct2*abs_coef_x(i)
                !phi1(i,j,k,is) = phi1(i,j,k,1)+(one - phi1(i,j,k,1))*abs_coef_x(i)
              enddo
            enddo
          enddo
        elseif ( phi_type(is) == 3 ) then !Numerical dye
          do k=1,xsize(3)
            do j=1,xsize(2)
              do i=nx_pf,xsize(1)
                phi1(i,j,k,is) = zero
              enddo
            enddo
          enddo
        endif
      enddo
    elseif ( sponge_type2 .eq. 3) then !Schuch et al. (2018) style
      !Sponge zone
      do is=1, nphi
        if ( phi_type(is) == 1 ) then !particle-laden
          do k=1,xsize(3)
            do j=1,xsize(2)
              do i=nx_pf,xsize(1)
                phi1(i,j,k,is) = phi1(i,j,k,is)*(one-ct2*abs_coef2_x(i))
              enddo
            enddo
          enddo
        elseif ( phi_type(is) == 2 ) then !ambient fluid
          do k=1,xsize(3)
            do j=1,xsize(2)
              do i=nx_pf,xsize(1)
                phi1(i,j,k,is) = phi1(i,j,k,is)+(one-phi1(i,j,k,is))*ct2*abs_coef2_x(i)
              enddo
            enddo
          enddo
        elseif ( phi_type(is) == 3 ) then !Numerical dye
          do k=1,xsize(3)
            do j=1,xsize(2)
              do i=nx_pf,xsize(1)
                phi1(i,j,k,is) = zero
              enddo
            enddo
          enddo
        endif
      enddo
    endif
  endif
  !Outflow
  do k=1,xsize(3)
    do j=1,xsize(2)
      phi1(nx,j,k,:)=phi1(nx,j,k,:)-cx*(phi1(nx,j,k,:)-phi1(nx-1,j,k,:))
    enddo
  enddo

  return
end subroutine boundary_conditions
!********************************************************************
subroutine flow_rate(plane1, q1, area)

  USE decomp_2d
  USE MPI
  USE param, only: zero

  implicit none

  real(mytype), intent(in), dimension(xsize(2),xsize(3)) :: plane1, area
  real(mytype), intent(out) :: q1
  real(mytype) :: q
  !
  integer :: j,k, code

  q = zero; q1 = zero
  do k=1,xsize(3)
    do j=1,xsize(2)
      q=q+plane1(j,k)*area(j,k)
    enddo
  enddo

  call MPI_ALLREDUCE(q,q1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

end subroutine flow_rate
!********************************************************************
subroutine init (ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI
  USE flow_type

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2

  real(mytype) :: y,ym,r,um,r1,r2,r3,x,z,h,ct,uinflow,trap
  integer :: k,j,i,ii,ijk,fh,ierror,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp
  !
  character(120) :: filename
  !
  integer, dimension (:), allocatable :: seed

  ct=one!half+half*tanh(t-three) !Acrescentado por Felipe Schuch
  phi1 = zero
  ux1=zero; uy1=zero; uz1=zero

  if (iin.eq.0) then !empty domain
    if (nrank==0) write(*,*) "Empty initial domain!"
  endif

  if (iin.eq.1) then !generation of a random noise
    ux1=zero; uy1=zero; uz1=zero
    !call ecoule(ux1,uy1,uz1)
    if (nrank==0) write(*,*) "Filled initial domain!"

    call system_clock(count=code)
    !if (iin.eq.2) code=0
    call random_seed(size = ii)
    call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)
    !
    ux1 = (two*ux1-one)*noise1*ct
    uy1 = (two*uy1-one)*noise1*ct
    uz1 = (two*uz1-one)*noise1*ct
    !
    do i=nx_pi+1,xsize(1)
      ux1(i,:,:) = zero
      uy1(i,:,:) = zero
      uz1(i,:,:) = zero
    enddo
    !Mudança por Felipe Schuch
    bxx1 = zero
    do j=1,xsize(2)
      y=(j+xstart(2)-1-1)*dy
      if (y.ge.yly-one) then
        !Inflow profile
        uinflow=tanh(sqrt(pi)/zpone*(y-(yly-one)))
        do k=1,xsize(3)
          bxx1(j,k)=uinflow
        enddo
      endif
    enddo
    !
    call flow_rate(bxx1, trap, area_yz1) !flow rate per width unit at inlet should be one
    trap = ct * zlz / trap
    if (nrank.eq.0) write(*,*) ' trap ux1 =',trap,'init'
    bxx1(:,:) = bxx1(:,:) * trap
    do i=1,nx_pi
      ux1(i,:,:) = ux1(i,:,:) + bxx1(:,:)
    enddo
  endif

  if (iin.eq.2) then
    !
    if (nrank==0) write(*,*) "Reading initial domain!"
    !
    do is=1, nphi
       write(filename,"('phi',I1.1,I4.4)") is, 0
       filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
       call decomp_2d_read_one(1,phi1(:,:,:,is),filename)
    enddo
    !
    write(filename,"('ux',I4.4)") 0
    filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
    call decomp_2d_read_one(1,ux1,filename)
    !
    write(filename,"('uy',I4.4)") 0
    filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
    call decomp_2d_read_one(1,uy1,filename)
    !
    write(filename,"('uz',I4.4)") 0
    filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
    call decomp_2d_read_one(1,uz1,filename)
    !
    ! Sponge zone
    !
    do i=1, sponge_nt
      sponge_ux1(i,:,:) = ux1(nx_xref,:,:)
    enddo
    do i=1, sponge_nt
      sponge_phi1(i,:,:,:) = phi1(nx_xref,:,:,:)
    enddo
  else
    if (iscalar==1) then
      do is=1,nphi
        if ( phi_type(is) == 2 ) then !ambient
          phi1(:,:,:,is) = cp(is)
        end if
        call transpose_x_to_y(phi1(:,:,:,is),phi2)
        call lagpolphiy(phi2)
        call transpose_y_to_x(phi2,phi1(:,:,:,is))
      enddo
    endif
  endif

  do is=1,nphi
    do ijk=1,xsize(1)*xsize(2)*xsize(3)
      phis1(ijk,1,1,is)=phi1(ijk,1,1,is)
      phiss1(ijk,1,1,is)=phis1(ijk,1,1,is)
    enddo
  enddo

  do k=1,xsize(3)
    do j=1,xsize(2)
      do i=1,xsize(1)
        gx1(i,j,k)=ux1(i,j,k)
        gy1(i,j,k)=uy1(i,j,k)
        gz1(i,j,k)=uz1(i,j,k)
        hx1(i,j,k)=gx1(i,j,k)
        hy1(i,j,k)=gy1(i,j,k)
        hz1(i,j,k)=gz1(i,j,k)
      enddo
    enddo
  enddo

#ifdef DEBG
  if (nrank .eq. 0) print *,'# init end ok'
#endif

  return
end subroutine init
!*******************************************************************
subroutine pre_correc(ux,uy,uz)

  USE decomp_2d
  USE variables
  USE param
  USE MPI
  USE flow_type, only: area_yz1, area_yzn, layer, y0ramp

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  integer :: i,j,k

  real(mytype) :: qx1,qxn,qmiss

  !********NCLZ==2*************************************
  if (nclz1==2) then
    if (xstart(3)==1) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
          dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
        enddo
      enddo
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux(i,j,1)=bzx1(i,j)+dpdxz1(i,j)
          uy(i,j,1)=bzy1(i,j)+dpdyz1(i,j)
          uz(i,j,1)=bzz1(i,j)
        enddo
      enddo
    endif
  endif

  if (nclzn==2) then
    if (xend(3)==nz) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
          dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
        enddo
      enddo
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux(i,j,xsize(3))=bzxn(i,j)+dpdxzn(i,j)
          uy(i,j,xsize(3))=bzyn(i,j)+dpdyzn(i,j)
          uz(i,j,xsize(3))=bzzn(i,j)
        enddo
      enddo
    endif
  endif
  !********NCLZ==1************************************* !just to reforce free-slip condition
  if (nclz1==1) then
    if (xstart(3)==1) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          uz(i,j,1)=zero
        enddo
      enddo
    endif
  endif

  if (nclzn==1) then
    if (xend(3)==nz) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          uz(i,j,xsize(3))=zero
        enddo
      enddo
    endif
  endif

  !********NCLX==2*************************************
  !we are in X pencils:
  if (nclx1==2.and.nclxn==2) then
    call flow_rate(bxx1, qx1, area_yz1)
    call flow_rate(bxxn, qxn, area_yzn)

    if (qxn.eq.zero) then
      qmiss = (qx1 - qxn)/((yly - layer)*zlz)
      bxxn(:,:) = bxxn(:,:) + qmiss
      qmiss = one + qmiss
    else
      qmiss = qx1/qxn
      bxxn(:,:) = bxxn(:,:) * qmiss
    endif

    if (nrank==0) print *,'Q/Lz   x1, xn, miss',real(qx1/zlz,4),&
                                                real(qxn/zlz,4),&
                                                real(qmiss,4)
  endif

  if (nclx1==2) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
        dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
      enddo
    enddo
    do k=1,xsize(3)
      do j=1,xsize(2)
        ux(1 ,j,k)=bxx1(j,k)
        uy(1 ,j,k)=bxy1(j,k)+dpdyx1(j,k)
        uz(1 ,j,k)=bxz1(j,k)+dpdzx1(j,k)
      enddo
    enddo
  endif
  if (nclxn==2) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
        dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
      enddo
    enddo
    do k=1,xsize(3)
      do j=1,xsize(2)
        ux(nx,j,k)=bxxn(j,k)
        uy(nx,j,k)=bxyn(j,k)+dpdyxn(j,k)
        uz(nx,j,k)=bxzn(j,k)+dpdzxn(j,k)
      enddo
    enddo
  endif

  !********NCLX==1*************************************
  if (nclx1==1) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        ux(1 ,j,k)=zero
      enddo
    enddo
  endif
  if (nclxn==1) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        ux(nx,j,k)=zero
      enddo
    enddo
  endif

  !********NCLY==2*************************************
  if (ncly1==2) then
    if (xstart(2)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
          dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
        enddo
      enddo
      do k=1,xsize(3)
        do i=1,xsize(1)
          ux(i,1,k)=byx1(i,k)+dpdxy1(i,k)
          uy(i,1,k)=byy1(i,k)
          uz(i,1,k)=byz1(i,k)+dpdzy1(i,k)
        enddo
      enddo
    endif
  endif

  if (nclyn==2) then
    if (xend(2)==ny) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
          dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
        enddo
      enddo
      do k=1,xsize(3)
        do i=1,xsize(1)
          ux(i,xsize(2),k)=byxn(i,k)+dpdxyn(i,k)
          uy(i,xsize(2),k)=byyn(i,k)
          uz(i,xsize(2),k)=byzn(i,k)+dpdzyn(i,k)
        enddo
      enddo
    endif
  endif

  !********NCLY==1*************************************
  if (ncly1==1) then
    if (xstart(2)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          uy(i,1,k)=zero
        enddo
      enddo
    endif
  endif

  if (nclyn==1) then
    if (xend(2)==ny) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          uy(i,xsize(2),k)=zero
        enddo
      enddo
    endif
  endif

  return
end subroutine pre_correc
