#ifdef POST
module post_processing

  USE decomp_2d
  USE variables
  USE param
  USE flow_type

  implicit none
  !
  type(decomp_info) :: xyplane
  type(decomp_info) :: xzplane
  !
  integer, save :: nxi, nxf
  real(mytype), save, allocatable, dimension(:,:,:) :: vol1
  real(mytype), save, allocatable, dimension(:,:) :: area2
  !
  integer :: FS
  character(len=160) :: fileformat, filename
  character(len=1),parameter :: NL=char(10) !new line character
  !
  !probes !só vai funcionar se a impressão for em relação ao lapis X!
  integer :: nprobes
  integer, save, allocatable, dimension(:) :: rankprobes
  integer, save, allocatable, dimension(:) :: nxprobes, nyprobes, nzprobes
  !
  real(mytype),save,allocatable,dimension(:,:,:) :: usum,vsum,wsum,uusum,uvsum
  real(mytype),save,allocatable,dimension(:,:,:) :: uwsum,vvsum,vwsum,wwsum
  real(mytype),save,allocatable,dimension(:,:,:,:) :: nphisum
  real(mytype),save,allocatable,dimension(:,:,:) :: phisum,uphisum,vphisum
  real(mytype),save,allocatable,dimension(:,:,:) :: wphisum,phiphisum
  !
  real(mytype),save,allocatable,dimension(:,:,:) :: presum, prepresum
  real(mytype),save,allocatable,dimension(:,:,:) :: preusum, prevsum, prewsum
contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI
    USE decomp_2d
    USE decomp_2d_io

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: dxdydz, dxdz, x, xprobes, yprobes, zprobes
    integer :: i,is,j,k,code
    character :: a

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post start'
#endif

    call decomp_info_init(nx, ny, 1, xyplane, .True.)
    call decomp_info_init(nx, 1, nz, xzplane, .True.)

    call alloc_z(usum, xyplane)
    call alloc_z(vsum, xyplane)
    call alloc_z(wsum, xyplane)
    call alloc_z(uusum, xyplane)
    call alloc_z(uvsum, xyplane)
    call alloc_z(uwsum, xyplane)
    call alloc_z(vvsum, xyplane)
    call alloc_z(vwsum, xyplane)
    call alloc_z(wwsum, xyplane)
    !
    allocate(nphisum(zsize(1),zsize(2),1,nphi))
    !
    call alloc_z(phisum, xyplane)
    call alloc_z(uphisum, xyplane)
    call alloc_z(vphisum, xyplane)
    call alloc_z(wphisum, xyplane)
    call alloc_z(phiphisum, xyplane)
    !
    call alloc_z(presum, xyplane)
    call alloc_z(prepresum, xyplane)
    call alloc_z(preusum, xyplane)
    call alloc_z(prevsum, xyplane)
    call alloc_z(prewsum, xyplane)
    !
    usum=zero;vsum=zero;wsum=zero
    uusum=zero;uvsum=zero;uwsum=zero
    vvsum=zero;vwsum=zero;wwsum=zero
    nphisum = zero
    phisum=zero;uphisum=zero;vphisum=zero;wphisum=zero;phiphisum=zero
    !
    presum = zero; prepresum = zero; preusum = zero
    prevsum = zero; prewsum = zero
    !
    if (nclx) then
      nxi=int(xlx_pi/dx) !Physical domain start point
      nxf=int(xlx_pf/dx) !Physical domain final point
    else
      nxi=int(xlx_pi/dx+1) !Physical domain start point
      nxf=int(xlx_pf/dx+1) !Physical domain final point
    end if
    if (nxi .lt. 0)  nxi = 0
    if (nxf .gt. nx) nxf = nx

    call alloc_x(vol1, opt_global=.true.)
    vol1 = zero

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico (método de Simpson)
    dxdydz=dx*dy*dz
    do k=xstart(3),xend(3)
      do j=xstart(2),xend(2)
        do i=nxi,nxf
          vol1(i,j,k)=dxdydz
          if (i .eq. nxi .or. i .eq. nxf) vol1(i,j,k) = vol1(i,j,k) * five/twelve
          if (i .eq. nxi+1 .or. i .eq. nxf-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
          if (j .eq. 1   .or. j .eq. ny)  vol1(i,j,k) = vol1(i,j,k) * five/twelve
          if (j .eq. 2   .or. j .eq. ny-1)  vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
          if (.not.nclz) then
            if (k .eq. 1   .or. k .eq. nz)  vol1(i,j,k) = vol1(i,j,k) * five/twelve
            if (k .eq. 2   .or. k .eq. nz-1)  vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
          endif
        end do
      end do
    end do

#ifdef IBM
    vol1 = vol1 * (one - ep1)
#endif

    !Y PENCILS
    allocate(area2(ystart(1):yend(1),ystart(3):yend(3)))
    dxdz=dx*dz
    area2=zero
    do k=ystart(3),yend(3)
      do i=ystart(1),yend(1)
        if (i .ge. nxi .and. i .le. nxf) then
          area2(i,k)=dxdz
          if (i .eq. nxi .or. i .eq. nxf) area2(i,k) = area2(i,k)/two
          if (k .eq. 1   .or. k .eq. nz)  area2(i,k) = area2(i,k)/two
          !if (ramp(i) .lt. y0ramp) area2(i,k) = area2(i,k) * ( - cos(declramp1))
          !if (ramp(i) .gt. layer)  area2(i,k) = area2(i,k) * ( - cos(declramp1))
        end if
      end do
    end do

    ! Writing CSV header
    open(67,file=trim(datapath)//'statistics.csv',status='unknown',&
    form='formatted',access='direct',recl=14)

    i=1
    write(67,'(13A,A)',rec=i) '            t', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '         xp_x', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '         xp_y', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '         xp_z', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '       xp2d_x', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '       xp2d_y', ','
    i=i+1
    do is=1,nphi
      write(67,'(A12,I1.1,A)',rec=i) '          mp', is, ','
      i=i+1
    enddo
    write(67,'(13A,A)',rec=i) '           vl', ','
    i=i+1
    do is=1,nphi
      write(67,'(A12,I1.1,A)',rec=i) '         dms', is, ','
      i=i+1
    enddo
    write(67,'(13A,A)',rec=i) '         xf_x', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '         xf_y', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '         xf_z', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '       xf2d_x', ','
    i=i+1
    write(67,'(13A,A)',rec=i) '       xf2d_y', NL

    close(67)

    if (flag_probes) then
    !probes
    !WORK X-PENCILS
      open(10,file='probes.csv',status='unknown',form='formatted')
      read (10,*) nprobes
      if (nrank.eq.0) write(*,"(' nprobes            : ',I15)") nprobes
      if (nprobes .gt. 0) then
         allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
         rankprobes(:)=0
         do i=1, nprobes
            read (10,*) nxprobes(i), nyprobes(i), nzprobes(i)
            !probes come from python, so we should add one
            nxprobes(i) = nxprobes(i) + 1
            nyprobes(i) = nyprobes(i) + 1
            nzprobes(i) = nzprobes(i) + 1
            if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
               if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
                  if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
                    rankprobes(i) = 1
                    !
                    write(filename,"('probe',I4.4,'.csv')") i
                    ! Writing CSV header
                    open(67,file=trim(datapath)//trim(filename),&
                         status='unknown',form='formatted',&
                         access='direct',recl=15)
                    k = 1
                    write(67,'(14A,A)',rec=k) '            t', ','
                    k = k + 1
                    write(67,'(14A,A)',rec=k) '           ux', ','
                    k = k + 1
                    write(67,'(14A,A)',rec=k) '           uy', ','
                    k = k + 1
                    write(67,'(14A,A)',rec=k) '           uz', ','
                    k = k + 1
                    do is=1,nphi
                      if (is.eq.nphi) then
                        write(67,'(A13,I1.1,A)',rec=k) '         phi', is, NL
                      else
                        write(67,'(A13,I1.1,A)',rec=k) '         phi', is, ','
                      endif
                      k = k + 1
                    enddo
                    close(67)
                    !
                  endif
               endif
            endif
         enddo
      endif
      close(10)
    endif

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine restart_post(irestart)

    USE decomp_2d_io
    USE variables
    USE param
    USE MPI
#ifdef IBM
    USE flow_type, only : sponge_ux1, sponge_phi1, sponge_init, sponge_type
#endif

    implicit none

    integer :: irestart,fh,ierror,code,is
    integer :: ierror_o=0 !error to open sauve file during restart
    character(len=30) :: filename, filestart
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: vec

#ifdef IBM
    if ( sponge_type .eq. 1 ) then
      vec = 0
      !
      write(filename, "('sauve-post',I7.7)") itime
      filename = trim(datapath)//trim(filename)
      write(filestart,"('sauve-post',I7.7)") ifirst-1
      filestart = trim(datapath)//trim(filestart)
      !
      if (irestart==1) then !write
        call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
        filesize = 0_MPI_OFFSET_KIND
        call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
        disp = 0_MPI_OFFSET_KIND
        vec(1) = sponge_init
        call decomp_2d_write_scalar(fh,disp,3,vec)
        if ( sponge_type .eq. 1 ) then
          call decomp_2d_write_var(fh,disp,1,sponge_ux1,sponge_tyzplane)
        endif
        if ( sponge_type2 .eq. 1 ) then
          do is=1, nphi
            call decomp_2d_write_var(fh,disp,1,sponge_phi1(:,:,:,is),sponge_tyzplane)
          enddo
        endif
        call MPI_FILE_CLOSE(fh,ierror)
      else !read
        call MPI_FILE_OPEN(MPI_COMM_WORLD, filestart, &
            MPI_MODE_RDONLY, MPI_INFO_NULL, &
            fh, ierror_o)
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_read_scalar(fh,disp,3,vec)
        sponge_init = vec(1)
        if ( sponge_type .eq. 1 ) then
          call decomp_2d_read_var(fh,disp,1,sponge_ux1,sponge_tyzplane)
        endif
        if ( sponge_type2 .eq. 1 ) then
          do is=1, nphi
            call decomp_2d_read_var(fh,disp,1,sponge_phi1(:,:,:,is),sponge_tyzplane)
          enddo
        endif
        call MPI_FILE_CLOSE(fh,ierror_o)
      endif

      if (nrank.eq.0) then
         if (ierror_o .ne. 0) then !Included by Felipe Schuch
            print *,'==========================================================='
            print *,'Error: Impossible to read '//trim(filestart)
            print *,'==========================================================='
            call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
         endif
      endif
    endif
#endif

  end subroutine restart_post
  !############################################################################
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    USE decomp_2d_io
    USE var, only : pre1

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    !
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2, uy2, uz2
    !
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3, uy3, uz3
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: phisum3, pre3
    !
    real(mytype),dimension(zsize(1),zsize(2),1) :: tmp_plane3
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1, di1, temp11, utau1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: temp2, di2, temp22, utau2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: temp3
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phisum1
    real(mytype),dimension(zsize(1),zsize(2)) :: phimsum3
    !
    real(mytype),dimension(ystart(1):yend(1),1,ystart(3):yend(3),nphi) :: dep2
    real(mytype),dimension(ysize(1),1,ysize(3)) :: utaumap2
    !
    integer :: i,j,k,is, j0, j1, tt
    character(len=120) :: filename
    real(mytype) :: thetas, thetac, dj
    real(mytype) :: xp(1:2,1:3),xp2d(1:2,1:2),mp(nphi),vl,dms(nphi),xf(1:2,1:3),xf2d(1:2,1:2)
    !
    logical :: flag

    if (mod(itime,iprocessing).eq.0) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Saving averaged planes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DOUBLE_PREC
      thetas = dsin(-atan(declramp))
      thetac = dcos(-atan(declramp))
#else
      thetas = sin(-atan(declramp))
      thetac = cos(-atan(declramp))
#endif
      call postprocessing_aux(ux1,ux2,ux3,'ux')
      call postprocessing_aux(uy1,uy2,uy3,'uy')
      call postprocessing_aux(uz1,uz2,uz3,'uz')

      phisum1=zero
      phimsum3=zero
      do is=1,nphi
        call postprocessing_aux(phi1(:,:,:,is),temp2,temp3,'phi'//char(is+48))
        if ( phi_type(is).eq.1) then
          phisum1(:,:,:)=phisum1(:,:,:)+phi1(:,:,:,is)
          phimsum3(:,:) = phimsum3(:,:) + temp3(:,:,1)
        endif
      enddo

      !utau - https://math.stackexchange.com/questions/1179818/gradient-in-a-rotated-reference-frame
      call derx (utau1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
      utau1 = thetas * thetac * utau1
      call derx (temp1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
      temp1 = thetas * thetas * temp1
      utau1 = utau1 - temp1
      call transpose_x_to_y(utau1,utau2)
      !
      call dery (temp2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
      temp2 = thetac * thetac * temp2
      utau2 = utau2 + temp2
      call dery (temp2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
      temp2 = thetac * thetas * temp2
      utau2 = utau2 - temp2
      !
      utau2 = utau2 * utau2
      !
      call derx (temp1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
      temp1 = thetas * temp1
      call transpose_x_to_y(temp1,temp2)
      !
      call dery (temp22,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
      !
      temp22 = thetac * temp22
      temp2 = temp2 + temp22
      temp2 = temp2 * temp2
      !
      utau2 = utau2 + temp2
      utau2 = sqrt(utau2)
      utau2 = xnu*utau2
      utau2 = sqrt(utau2)
      !
#ifdef IBM
      do k=1,ysize(3)
        do i=1,ysize(1)
          j0=int(ramp(i+ystart(1)-1)/dy)+1
          j1=j0+1
          dj = ramp(i+ystart(1)-1) - real(j0 - 1, mytype)*dy
          utaumap2(i,1,k) = utau2(i,j0,k) + dj*(utau2(i,j1,k)&
          &-utau2(i,j0,k))/dy!interpolação entre j0 e j1
        end do
      end do
#else
      utaumap2(:,1,:) = utau2(:,1,:)
#endif
      write(filename,"('utau',I4.4)") itime/iprocessing
      filename = trim(datapath)//'/xz-plane/'//trim(filename)
      call decomp_2d_write_one(2,utaumap2,filename,xzplane)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Computing classic post-processing
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      xp=zero; xp2d=zero; mp=zero; vl=zero; dms=zero; xf=zero; xf2d=zero

#ifdef IBM
      call dep_ibm(phi1,dep2)
      call plunge_point(phisum1,xp)
      call plunge_point2d(phimsum3,xp2d)
#else
      call dep_no_ibm(phi1,dep2)
      call budget(ux1,uy1,uz1,phi1,vol1)
#endif
      call suspended(phi1,vol1,mp,vl)
      call depositrate (dep2,dms)
      call front(phisum1,xf)
      call front2d(phimsum3,xf2d)

      if (nrank .eq. 0) then
        FS = 1+3+2+nphi+1+nphi+3+2 !Number of columns
        write(fileformat, '( "(",I4,"(E13.6,:,"",""),E13.6,A)" )' ) FS-1
        FS = FS*14 !Line width
        open(67,file=trim(datapath)//'statistics.csv',status='unknown',&
        form='formatted',access='direct',recl=FS)
        write(67,fileformat,rec=itime/iprocessing+2) t,& !1
        xp(1,:),&                                        !3
        xp2d(1,:),&                                      !2
        mp,&                                             !nphi
        vl,&                                             !1
        dms,&                                            !nphi
        xf(1,:),&                                        !3
        xf2d(1,:),NL                                     !2
        close(67)
      end if
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Computing classic stats
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (itime.ge.initstats1) then
      !
      flag = .False.
      if (iin.eq.2) then !LST simulation writes variables more often
        if (mod(itime,iprocessing).eq.0) flag = .True.
        tt = itime/iprocessing
      else
        if (mod(itime,imodulo).eq.0) flag = .True.
        tt = itime/imodulo
      endif
      !
      if (iin.eq.2) then !LST simulation does not sum variables
        usum = zero; vsum = zero; wsum = zero; phisum = zero
        uusum = zero; uvsum = zero; uwsum = zero
        vvsum = zero; vwsum = zero; wwsum = zero
        nphisum = zero
        uphisum = zero; vphisum = zero
        wphisum = zero; phiphisum = zero
      endif
      !
      temp1 = ux1*(one-ep1)
      call transpose_x_to_y(temp1,temp2)
      call transpose_y_to_z(temp2,ux3)
      temp1 = uy1*(one-ep1)
      call transpose_x_to_y(temp1,temp2)
      call transpose_y_to_z(temp2,uy3)
      temp1 = uz1*(one-ep1)
      call transpose_x_to_y(temp1,temp2)
      call transpose_y_to_z(temp2,uz3)
      !
      phisum3 = zero
      do is=1, nphi
        temp1 = phi1(:,:,:,is)*(one-ep1)
        call transpose_x_to_y(temp1,temp2)
        call transpose_y_to_z(temp2,temp3)
        call avg_and_sum(nphisum(:,:,:,is),temp3,flag,'phi'//char(is+48)//'sum',tt)
        if (phi_type(is).eq.1) then
          phisum3 = phisum3 + temp3
        endif
      enddo
      !
      temp1 = pre1*(one-ep1)
      call transpose_x_to_y(temp1,temp2)
      call transpose_y_to_z(temp2,pre3)
      !
      call avg_and_sum(usum,ux3,flag,'usum',tt)
      call avg_and_sum(vsum,uy3,flag,'vsum',tt)
      call avg_and_sum(wsum,uz3,flag,'wsum',tt)
      call avg_and_sum(phisum,phisum3,flag,'phisum',tt)
      !
      call avg_and_sum(uusum,ux3*ux3,flag,'uusum',tt)
      call avg_and_sum(uvsum,ux3*uy3,flag,'uvsum',tt)
      call avg_and_sum(uwsum,ux3*uz3,flag,'uwsum',tt)
      call avg_and_sum(vvsum,uy3*uy3,flag,'vvsum',tt)
      call avg_and_sum(vwsum,uy3*uz3,flag,'vwsum',tt)
      call avg_and_sum(wwsum,uz3*uz3,flag,'wwsum',tt)
      !
      call avg_and_sum(uphisum,ux3*phisum3,flag,'uphisum',tt)
      call avg_and_sum(vphisum,uy3*phisum3,flag,'vphisum',tt)
      call avg_and_sum(wphisum,uz3*phisum3,flag,'wphisum',tt)
      call avg_and_sum(phiphisum,phisum3*phisum3,flag,'phiphisum',tt)
      !
      call avg_and_sum(presum,pre3,flag,'psum',tt)
      call avg_and_sum(prepresum,pre3*pre3,flag,'ppsum',tt)
      call avg_and_sum(preusum,pre3*ux3,flag,'pusum',tt)
      call avg_and_sum(prevsum,pre3*uy3,flag,'pvsum',tt)
      call avg_and_sum(prewsum,pre3*uz3,flag,'pwsum',tt)
      !
      if (flag) then !write results
!         !!
!         ! Derivatives
!         !!
!         temp3 = zero
!         !ux
!         temp3(:,:,1) = usum(:,:,1)/real(ntimes1,mytype)
!         call transpose_z_to_y(temp3, temp2)
!         call transpose_y_to_x(temp2, temp1)
!         call derx (temp11,temp1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
!         call dery (temp22,temp2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!         write(filename,"('dudx',I4.4)") tt
!         filename = trim(datapath)//'/xy-plane/'//trim(filename)
!         call decomp_2d_write_plane(1,temp11,3,1,filename)
!         write(filename,"('dudy',I4.4)") tt
!         filename = trim(datapath)//'/xy-plane/'//trim(filename)
!         call decomp_2d_write_plane(2,temp22,3,1,filename)
!         !uy
!         temp3(:,:,1) = vsum(:,:,1)/real(ntimes1,mytype)
!         call transpose_z_to_y(temp3, temp2)
!         call transpose_y_to_x(temp2, temp1)
!         call derx (temp11,temp1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!         call dery (temp22,temp2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
!         write(filename,"('dvdx',I4.4)") tt
!         filename = trim(datapath)//'/xy-plane/'//trim(filename)
!         call decomp_2d_write_plane(1,temp11,3,1,filename)
!         write(filename,"('dvdy',I4.4)") tt
!         filename = trim(datapath)//'/xy-plane/'//trim(filename)
!         call decomp_2d_write_plane(2,temp22,3,1,filename)
!         !phi
!         if (ivirt==2) ilag=0
!         temp3(:,:,1) = phisum(:,:,1)/real(ntimes1,mytype)
!         call transpose_z_to_y(temp3, temp2)
! #ifdef IBM
!         call lagpolphiy(temp2)
! #endif
!         call transpose_y_to_x(temp2, temp1)
!         call derxS (temp11,temp1,di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1)
!         call deryS (temp22,temp2,di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1)
!         write(filename,"('dphidx',I4.4)") tt
!         filename = trim(datapath)//'/xy-plane/'//trim(filename)
!         call decomp_2d_write_plane(1,temp11,3,1,filename)
!         write(filename,"('dphidy',I4.4)") tt
!         filename = trim(datapath)//'/xy-plane/'//trim(filename)
!         call decomp_2d_write_plane(2,temp22,3,1,filename)
!         if (ivirt==2) ilag=1
        ! Reset accumulators after writing
        usum = zero; vsum = zero; wsum = zero; phisum = zero
        uusum = zero; uvsum = zero; uwsum = zero
        vvsum = zero; vwsum = zero; wwsum = zero
        nphisum = zero
        uphisum = zero; vphisum = zero
        wphisum = zero; phiphisum = zero
      endif
    endif
    return
  end subroutine postprocessing
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

    integer :: i, code

    if (itime.ge.initstats1) then
      FS = 1+3+nphi !Number of columns
      write(fileformat, '( "(",I4,"(E14.7,:,"",""),E14.7,A)" )' ) FS-1
      FS = FS*15  !Line width
      do i=1, nprobes
         if (rankprobes(i) .eq. 1) then
            write(filename,"('probe',I4.4,'.csv')") i
            open(67,file=trim(datapath)//trim(filename),status='unknown',&
                form='formatted',access='direct',recl=FS)
            write(67,fileformat,rec=itime-initstats1+2) t,&            !1
                  ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
                  uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
                  uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
                  phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),NL        !nphi
            close(67)
         endif
      enddo
    endif

  end subroutine write_probes
  !############################################################################
  subroutine suspended(phi1,vol1,mp1,vl1)

    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: vol1
    real(mytype),intent(out) :: mp1(1:nphi),vl1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
    real(mytype) :: mp(1:nphi),vl
    integer :: is,code

    mp=zero; mp1=zero
    vl=zero; vl1=zero

    do is=1, nphi
      temp1 = phi1(:,:,:,is)*vol1(:,:,:)
      mp(is)= sum(temp1)
    end do
    vl=sum(vol1)

    call MPI_REDUCE(mp,mp1,nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(vl,vl1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    return
  end subroutine suspended
  !############################################################################
  subroutine plunge_point ( phi1, xp)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: phi1

    real(mytype),intent(out) :: xp(1:2,1:3)
    real(mytype) :: xp1(1:2)

    integer :: i, j ,k,code

    xp(2,:) = real(nrank,mytype)
    xp(1,1:3)=zero
    xp1=zero

    if (xend(2) .eq. ny) then
      j=ny
      kloop: do k =xstart(3),xend(3)
        iloop: do i=xstart(1),xend(1)
          if ( phi1(i,j,k) .lt. 0.01_mytype ) then
            xp(1,1:3) = (/ real (i-1,4)*dx, real (j-1,4)*dy, real(k-1,4)*dz /)
            exit kloop
          end if
        end do iloop
      end do kloop
    end if

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 3,real2_type, int(xp1(2)), MPI_COMM_WORLD,code)

    return
  end subroutine plunge_point
  !############################################################################
  subroutine plunge_point2d ( phimsum3, xp )

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(zstart(1):zend(1),zstart(2):zend(2)) :: phimsum3
    real(mytype),intent(out) :: xp(1:2,1:2)

    real(mytype) :: xp1(1:2)
    integer :: i, j ,code

    xp(2,:) = real(nrank,mytype)
    xp(1,1:2)=zero
    xp1=zero

    if (zend(2) .eq. ny .and. phimsum3(zstart(1),zend(2)) .ge. 0.01_mytype) then
      j=ny
      iloop: do i=zstart(1),zend(1)
        if ( phimsum3(i,j) .lt. 0.01_mytype ) then
          xp(1,1:2) = (/ real (i-1,4)*dx, real (j-1,4)*dy /)
          exit iloop
        end if
      end do iloop
    end if

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 2,real2_type, int(xp1(2)), MPI_COMM_WORLD,code)

    return
  end subroutine plunge_point2d
  !############################################################################
  subroutine front ( phisum1, xp )

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: phisum1
    real(mytype),intent(out) :: xp(1:2,1:3)

    real(mytype) :: xp1(1:2),y
    integer :: is, i, j ,k, code

    xp(2,:) = real(nrank,mytype)
    xp(1,1:3)=zero
    xp1=zero
    kloop: do k=xstart(3),xend(3)
      jloop: do j=xstart(2),xend(2)
        y = ( j - 1 )*dy
        iloop: do i=nxf, nxi, -1
          if ( phisum1(i,j,k) .ge. 0.01_mytype .and. y .ge. ramp(i)) then
            xp(1,1:3) = (/ real (i-1,4)*dx, real (j-1,4)*dy, real(k-1,4)*dz /)
            exit kloop
          end if
        end do iloop
      end do jloop
    end do kloop

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 3,real2_type, int(xp1(2)), MPI_COMM_WORLD,code)


  end subroutine front
  !############################################################################
  subroutine front2d ( phimsum3, xp )

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(zstart(1):zend(1),zstart(2):zend(2)) :: phimsum3
    real(mytype),intent(out) :: xp(1:2,1:2)
    real(mytype) :: xp1(1:2),y
    integer :: is, i, j, code

    xp(2,:) = real(nrank,mytype)
    xp(1,1:2)=zero
    xp1=zero
    jloop: do j=zstart(2),zend(2)
      y = ( j - 1 )*dy
      iloop: do i=zend(1), zstart(1), -1
        if ( phimsum3(i,j) .ge. 0.01_mytype .and. y .ge. ramp(i)) then
          xp(1,1:2) = (/ real (i-1,4)*dx, real (j-1,4)*dy /)
          exit jloop
        end if
      end do iloop
    end do jloop

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 2,real2_type, int(xp1(2)), MPI_COMM_WORLD,code)

#ifdef IBM
    if ( xp(1,1) .eq. xlx .and. sponge_init .eq. 100000000) then
      sponge_init = itime + sponge_nt * 2
    endif
#endif

  end subroutine front2d
  !############################################################################
  subroutine dep_ibm(phi1,dep2)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    real(mytype),intent(out),dimension(ystart(1):yend(1),1,ystart(3):yend(3),nphi) :: dep2
    real(mytype),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: phi2
    !
    real(mytype) :: dj, x
    integer :: i, j, k, j0, j1, is, code
    character(len=30) :: filename

    dep2 = zero

    do is=1, nphi
      if (uset(is) .eq. 0.) cycle
      call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:))

      !Y PENCILS
      do k=ystart(3),yend(3)
        do i=ystart(1),yend(1)
          j0=int(ramp(i)/dy)+1
          j1=j0+1
          dj = ramp(i) - real(j0 - 1, mytype)*dy
          dep2(i,1,k,is) = phi2(i,j0,k) + dj*(phi2(i,j1,k)&
          &-phi2(i,j0,k))/dy!interpolação entre j0 e j1
          dep2(i,1,k,is) = dep2(i,1,k,is)*uset(is)
        end do
      end do

      write(filename,"('dep',I1.1,I4.4)") is, itime/iprocessing
      filename = trim(datapath)//'/xz-plane/'//trim(filename)
      call decomp_2d_write_one(2,dep2(:,:,:,is),filename,xzplane)
    enddo

  end subroutine dep_ibm
  !############################################################################
  subroutine dep_no_ibm(phi1,dep2)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nphi) :: phi1
    real(mytype),intent(out),dimension(ystart(1):yend(1),1,ystart(3):yend(3),nphi) :: dep2
    !
    real(mytype),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),nphi) :: phi2

    integer :: i, j, k, is
    character(len=30) :: filename

    dep2 = 0.
    do is=1, nphi
      if (uset(is) .eq. 0.) cycle
      call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

      do k=ystart(3),yend(3)
        do i=ystart(1),yend(1)
          dep2(i,1,k,is) = phi2(i,1,k,is)*uset(is)
        end do
      end do

      write(filename,"('dep',I1.1,I4.4)") is, itime/iprocessing
      filename = trim(datapath)//'/xz-plane/'//trim(filename)
      call decomp_2d_write_one(2,dep2(:,:,:,is),filename,xzplane)
    enddo

  end subroutine dep_no_ibm
  !############################################################################
  subroutine depositrate ( dep2, dms1)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(ystart(1):yend(1),1,ystart(3):yend(3),nphi) :: dep2
    real(mytype),intent(out) :: dms1(nphi)

    real(mytype) :: dms(nphi)
    integer :: i,k,is,code

    dms=zero; dms1=zero
    do is=1, nphi
      if (uset(is) .eq. 0.) cycle
      do k=ystart(3),yend(3)
        do i=ystart(1),yend(1)
          dms(is)=dms(is)+dep2(i,1,k,is)*area2(i,k)
        end do
      end do
    enddo

    call MPI_REDUCE(dms,dms1,nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)


  end subroutine depositrate
  !############################################################################
  subroutine budget(ux1,uy1,uz1,phi1,vol1)

    USE decomp_2d
    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,vol1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1,dphiy1, dphixx1, dphiyy1, dphizz1, ddphi1, di1, temp1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2, dphiy2, dphixx2, dphiyy2, dphizz2, ddphi2, di2, vol2, temp2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: phi3, dphiy3, dphixx3, dphiyy3, dphizz3, ddphi3, di3, temp3

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3
    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A

    real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

    real(8) :: ek,ek1,dek,dek1,ep,ep1,dep,dep1,xvol
    integer :: ijk,i,j,k,l,m,is,code
    character(len=30) :: filename

    ek=zero;ek1=zero;dek=zero;dek1=zero;ep=zero;ep1=zero;dep=zero;dep1=zero;diss1=zero

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
    call transpose_x_to_y(ux1,td2)
    call transpose_x_to_y(uy1,te2)
    call transpose_x_to_y(uz1,tf2)
    call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    !!z-derivatives
    call transpose_y_to_z(td2,td3)
    call transpose_y_to_z(te2,te3)
    call transpose_y_to_z(tf2,tf3)
    call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)

    A(:,:,:,:,:)=zero
    A(1,1,:,:,:)=ta1(:,:,:)
    A(2,1,:,:,:)=tb1(:,:,:)
    A(3,1,:,:,:)=tc1(:,:,:)
    A(1,2,:,:,:)=td1(:,:,:)
    A(2,2,:,:,:)=te1(:,:,:)
    A(3,2,:,:,:)=tf1(:,:,:)
    A(1,3,:,:,:)=tg1(:,:,:)
    A(2,3,:,:,:)=th1(:,:,:)
    A(3,3,:,:,:)=ti1(:,:,:)

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             do m=1,3
                do l=1,3
                   diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**two
                enddo
             enddo
          enddo
       enddo
    enddo

    do ijk=1,xsize(1)*xsize(2)*xsize(3)
       xvol=real(vol1(ijk,1,1),8)
       ek = ek + half * xvol * (ux1(ijk,1,1)**two+uy1(ijk,1,1)**two+uz1(ijk,1,1)**two)
       dek = dek + xvol * diss1(ijk,1,1)
    enddo

    call transpose_x_to_y(vol1,vol2)

    if (ivirt==2) then
       ilag=0
    endif
    do is=1, nphi

       if (ri(is) .eq. 0.) cycle

       if (uset(is).eq.0..and.nclySn.ne.1) then
          nclySn=1
          if (nrank.eq.0) write(*,*) 'nclySn 2 -> 1   is =', is
          call schemes(.False.)
       endif
       if (uset(is).ne.0..and.nclySn.ne.2) then
          nclySn=2
          if (nrank.eq.0) write(*,*) 'nclySn 1 -> 2   is =', is
          call schemes(.False.)
       endif

       call derxxS (dphixx1,phi1(:,:,:,is),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1)

       call transpose_x_to_y(dphixx1,dphixx2)

       call transpose_x_to_y(phi1(:,:,:,is),phi2)

       call deryS (dphiy2,phi2,di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),1)

       call deryyS (dphiyy2,phi2,di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1)

       call transpose_y_to_z(phi2,phi3)

       call derzzS (dphizz3,phi3,di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1)

       call transpose_z_to_y(dphizz3,dphizz2)

       do ijk=1,ysize(1)*ysize(2)*ysize(3)
          ddphi2(ijk,1,1)=dphixx2(ijk,1,1)+dphiyy2(ijk,1,1)+dphizz2(ijk,1,1)
       enddo

       do k=1,ysize(3)
          do j=1,ysize(2)
             do i=1,ysize(1)
                xvol=real(vol2(i,j,k),8)
                ep=ep + xvol * (phi2(i,j,k)*(j-1)*dy)
                dep=dep - xvol * (ddphi2(i,j,k)*xnu/nsc(is)+uset(is)*dphiy2(i,j,k))*(j-1)*dy
             enddo
          enddo
       enddo
    enddo
    if (ivirt==2) then
       ilag=1
    endif

    call MPI_REDUCE(ek,ek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dek,dek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(ep,ep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dep,dep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)

    if (nrank .eq. 0) then

      if (itime .eq. 0) then
        ! Writing CSV header
        open(67,file=trim(datapath)//'budget.csv',status='unknown',&
        form='formatted',access='direct',recl=14)
        i=1
        write(67,'(13A,A)',rec=i) '            t', ','
        i=i+1
        write(67,'(13A,A)',rec=i) '           ek', ','
        i=i+1
        write(67,'(13A,A)',rec=i) '          dek', ','
        i=i+1
        write(67,'(13A,A)',rec=i) '           ep', ','
        i=i+1
        write(67,'(13A,A)',rec=i) '          dep', NL
        close(67)
      endif

      FS = 5 !Number of columns
      write(fileformat, '( "(",I4,"(E13.6,:,"",""),E13.6,A)" )' ) FS-1
      FS = FS*14 !Line width
      open(67,file=trim(datapath)//'budget.csv',status='unknown',&
      form='formatted',access='direct',recl=FS)
      write(67,fileformat,rec=itime/iprocessing+2) t,& !1
      ek1,&                                            !1
      dek1,&                                           !1
      ep1,&                                            !1
      dep1,NL                                          !1
      close(67)
    end if

  end subroutine budget
end module post_processing

subroutine postprocessing_aux(u1,u2,u3,name)

  USE decomp_2d
  USE decomp_2d_io
  USE param, only : datapath, one, itime, iprocessing
  USE var, only : ep1
  USE variables, only: ny, nz
  USE post_processing, only: filename

  implicit none

  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: u1
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u2
  real(mytype),intent(out),dimension(zsize(1),zsize(2),zsize(3)) :: u3
  character(len=*), intent(in) :: name
  !
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tmp1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tmp2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: tmp3
  character(len=4) :: num
  !
  write(num,"(I4.4)") itime/iprocessing
  !
  tmp1 = u1*(one-ep1)
  call transpose_x_to_y (tmp1,u2)
  call transpose_y_to_z (u2,u3)
  !depth averaged
  call mean_plane_y(u2,ysize(1),ysize(2),ysize(3),tmp2(:,1,:))
  filename = trim(datapath)//'/xz-plane/'//trim(name)//trim(num)
  call decomp_2d_write_plane(2,tmp2,2,1,filename)
  !spanwise averaged
  call mean_plane_z(u3,zsize(1),zsize(2),zsize(3),tmp3(:,:,1))
  filename = trim(datapath)//'/xy-plane/'//trim(name)//trim(num)
  call decomp_2d_write_plane(3,tmp3,3,1,filename)
  !center plane
  filename = trim(datapath)//'/xy-plane/'//trim(name)//'c'//trim(num)
  call decomp_2d_write_plane(3,u3,3,nz/2,filename)
  !top plane
  filename = trim(datapath)//'/xz-plane/'//trim(name)//'t'//trim(num)
  call decomp_2d_write_plane(2,u2,2,ny,filename)


end subroutine postprocessing_aux
!###############################################################################
subroutine avg_and_sum(plane,field,flag,prefix,num)

  USE decomp_2d_io
  USE decomp_2d, only : mytype, zsize
  USE post_processing, only : xyplane, filename
  USE param, only : datapath

  implicit none

  real(mytype),intent(inout),dimension(zsize(1),zsize(2),1) :: plane
  real(mytype),intent(in),dimension(zsize(1),zsize(2),zsize(3)) :: field
  logical, intent(in) :: flag
  character(len=*), intent(in) :: prefix
  integer, intent(in) :: num
  !
  real(mytype),dimension(zsize(1),zsize(2),1) :: tmp


  call mean_plane_z(field,zsize(1),zsize(2),zsize(3),tmp(:,:,1))
  plane = plane + tmp

  if (flag) then
    write(filename,"(I4.4)") num
    filename = trim(datapath)//'/xy-plane/'//prefix//trim(filename)
    call decomp_2d_write_one(3,plane,filename,xyplane)
  endif

end subroutine avg_and_sum
#endif
