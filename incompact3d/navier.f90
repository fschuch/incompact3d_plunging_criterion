!******************************************************************
subroutine  intt (ux,uy,uz,gx,gy,gz,hx,hy,hz,ta1,tb1,tc1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  integer :: ijk,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx,gy,gz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx,hy,hz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1

#ifdef DEBG
  if (nrank .eq. 0) print *,'# intt start'
#endif

  nxyz=xsize(1)*xsize(2)*xsize(3)

  if ((nscheme.eq.1).or.(nscheme.eq.3)) then !AB2 or RK3
    if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
    (nscheme.eq.3.and.itr.eq.1)) then
    do ijk=1,nxyz
      ux(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+ux(ijk,1,1)
      uy(ijk,1,1)=gdt(itr)*tb1(ijk,1,1)+uy(ijk,1,1)
      uz(ijk,1,1)=gdt(itr)*tc1(ijk,1,1)+uz(ijk,1,1)
      gx(ijk,1,1)=ta1(ijk,1,1)
      gy(ijk,1,1)=tb1(ijk,1,1)
      gz(ijk,1,1)=tc1(ijk,1,1)
    enddo
  else
    if (nz.gt.1) then
      do ijk=1,nxyz
        ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
        uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)
        uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+uz(ijk,1,1)
        gx(ijk,1,1)=ta1(ijk,1,1)
        gy(ijk,1,1)=tb1(ijk,1,1)
        gz(ijk,1,1)=tc1(ijk,1,1)
      enddo
    else
      do ijk=1,nxyz
        ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
        uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)
        gx(ijk,1,1)=ta1(ijk,1,1)
        gy(ijk,1,1)=tb1(ijk,1,1)
      enddo
    endif
  endif
endif

if (nscheme.eq.2) then !AB3
  if ((itime.eq.1).and.(ilit.eq.0)) then
    if (nrank==0) print *,'start with Euler',itime
    do ijk=1,nxyz !start with Euler
      ux(ijk,1,1)=dt*ta1(ijk,1,1)+ux(ijk,1,1)
      uy(ijk,1,1)=dt*tb1(ijk,1,1)+uy(ijk,1,1)
      uz(ijk,1,1)=dt*tc1(ijk,1,1)+uz(ijk,1,1)
      gx(ijk,1,1)=ta1(ijk,1,1)
      gy(ijk,1,1)=tb1(ijk,1,1)
      gz(ijk,1,1)=tc1(ijk,1,1)
    enddo
  else
    if  ((itime.eq.2).and.(ilit.eq.0)) then
      if (nrank==0) print *,'then with AB2',itime
        do ijk=1,nxyz
          ux(ijk,1,1)=onepfive*dt*ta1(ijk,1,1)-half*dt*gx(ijk,1,1)+ux(ijk,1,1)
          uy(ijk,1,1)=onepfive*dt*tb1(ijk,1,1)-half*dt*gy(ijk,1,1)+uy(ijk,1,1)
          uz(ijk,1,1)=onepfive*dt*tc1(ijk,1,1)-half*dt*gz(ijk,1,1)+uz(ijk,1,1)
          hx(ijk,1,1)=gx(ijk,1,1)
          hy(ijk,1,1)=gy(ijk,1,1)
          hz(ijk,1,1)=gz(ijk,1,1)
          gx(ijk,1,1)=ta1(ijk,1,1)
          gy(ijk,1,1)=tb1(ijk,1,1)
          gz(ijk,1,1)=tc1(ijk,1,1)
        enddo
      else
        do ijk=1,nxyz
          ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+&
          cdt(itr)*hx(ijk,1,1)+ux(ijk,1,1)
          uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+&
          cdt(itr)*hy(ijk,1,1)+uy(ijk,1,1)
          uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+&
          cdt(itr)*hz(ijk,1,1)+uz(ijk,1,1)
          hx(ijk,1,1)=gx(ijk,1,1)
          hy(ijk,1,1)=gy(ijk,1,1)
          hz(ijk,1,1)=gz(ijk,1,1)
          gx(ijk,1,1)=ta1(ijk,1,1)
          gy(ijk,1,1)=tb1(ijk,1,1)
          gz(ijk,1,1)=tc1(ijk,1,1)
        enddo
      endif
    endif
  endif

  return

#ifdef DEBG
  if (nrank .eq. 0) print *,'# intt done'
#endif

end subroutine intt
!********************************************************************
subroutine corgp (ux,gx,uy,uz,px,py,pz)

  USE decomp_2d
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  integer :: ijk,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx

  do ijk=1,xsize(1)*xsize(2)*xsize(3)
    ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
    uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1)
    uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1)
  enddo

  return
end subroutine corgp
!********************************************************************
subroutine divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
  td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
  nxmsize,nymsize,nzmsize,ph1,ph3,ph4,nlock)

  USE param
  USE decomp_2d
  USE variables
  USE MPI

  implicit none

  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

  !X PENCILS NX NY NZ  -->NXM NY NZ
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1,ux1,uy1,uz1,ep1
  real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1
  !Y PENCILS NXM NY NZ  -->NXM NYM NZ
  real(mytype),dimension(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)) :: td2,te2,tf2,di2
  real(mytype),dimension(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)) :: ta2,tb2,tc2
  !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
  real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)) :: ta3,tb3,tc3,di3
  real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: td3,te3,tf3,pp3

  integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nlock
  integer :: nxmsize,nymsize,nzmsize,code
  real(mytype) :: tmax,tmoy,tmax1,tmoy1


  nvect1=xsize(1)*xsize(2)*xsize(3)
  nvect2=ysize(1)*ysize(2)*ysize(3)
  nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

  if (ivirt.eq.0.and.ilag.eq.0) ep1(:,:,:)=zero
  do ijk=1,nvect1
    ta1(ijk,1,1)=(one-ep1(ijk,1,1))*ux1(ijk,1,1)
    tb1(ijk,1,1)=(one-ep1(ijk,1,1))*uy1(ijk,1,1)
    tc1(ijk,1,1)=(one-ep1(ijk,1,1))*uz1(ijk,1,1)
  enddo

  !WORK X-PENCILS
  call decx6(td1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
  call inter6(te1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
  call inter6(tf1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

  call transpose_x_to_y(td1,td2,ph4)!->NXM NY NZ
  call transpose_x_to_y(te1,te2,ph4)
  call transpose_x_to_y(tf1,tf2,ph4)


  !WORK Y-PENCILS
  call intery6(ta2,td2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
  call decy6(tb2,te2,di2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)
  call intery6(tc2,tf2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

  call transpose_y_to_z(ta2,ta3,ph3)!->NXM NYM NZ
  call transpose_y_to_z(tb2,tb3,ph3)
  call transpose_y_to_z(tc2,tc3,ph3)


  !WORK Z-PENCILS
  call interz6(td3,ta3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
  (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
  call interz6(te3,tb3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
  (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
  call decz6(tf3,tc3,di3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
  (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

  do k=1,nzmsize
    do j=ph1%zst(2),ph1%zen(2)
      do i=ph1%zst(1),ph1%zen(1)
        pp3(i,j,k)=td3(i,j,k)+te3(i,j,k)+tf3(i,j,k)
      enddo
    enddo
  enddo

  if (nlock==2) then
    pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
  endif

  tmax=-1609._mytype
  tmoy=0._mytype
  do k=1,nzmsize
    do j=ph1%zst(2),ph1%zen(2)
      do i=ph1%zst(1),ph1%zen(1)
        if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
        tmoy=tmoy+abs(pp3(i,j,k))
      enddo
    enddo
  enddo
  tmoy=tmoy/nvect3

  call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

  if (nrank==0) then
    if (nlock==2) then
      print *,'DIV U  max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
    else
      print *,'DIV U* max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
    endif
  endif

  return
end subroutine divergence
!*******************************************************************
subroutine gradp(ta1,tb1,tc1,di1,td2,tf2,ta2,tb2,tc2,di2,&
  ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

  USE param
  USE decomp_2d
  USE variables
  USE MPI

  implicit none

  TYPE(DECOMP_INFO) :: ph2,ph3
  integer :: i,j,k,ijk,nxmsize,nymsize,nzmsize,code
  integer, dimension(2) :: dims, dummy_coords
  logical, dimension(2) :: dummy_periods

  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
  !Z PENCILS NXM NYM NZM-->NXM NYM NZ
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,tc3,di3
  !Y PENCILS NXM NYM NZ -->NXM NY NZ
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2,tc2
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,td2,tf2,di2
  !X PENCILS NXM NY NZ  -->NX NY NZ
  real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1

  !WORK Z-PENCILS
  call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
  (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
  call deciz6(tc3,pp3,di3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
  (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

  !WORK Y-PENCILS
  call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
  call transpose_z_to_y(tc3,tc2,ph3)

  call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
  (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  call deciy6(td2,ta2,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
  (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  call interiy6(tf2,tc2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
  (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

  !WORK X-PENCILS

  call transpose_y_to_x(tb2,td1,ph2) !nxm ny nz
  call transpose_y_to_x(td2,te1,ph2)
  call transpose_y_to_x(tf2,tf1,ph2)

  call deci6(ta1,td1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)
  call interi6(tb1,te1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)
  call interi6(tc1,tf1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
  nxmsize,xsize(1),xsize(2),xsize(3),1)

  !we are in X pencils:
  if (nclx1.eq.2) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        dpdyx1(j,k)=tb1(1,j,k)/gdt(itr)
        dpdzx1(j,k)=tc1(1,j,k)/gdt(itr)
      enddo
    enddo
  endif
  if (nclxn.eq.2) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        dpdyxn(j,k)=tb1(nx,j,k)/gdt(itr)
        dpdzxn(j,k)=tc1(nx,j,k)/gdt(itr)
      enddo
    enddo
  endif

  if (ncly1.eq.2) then
    if (xsize(2)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxy1(i,k)=ta1(i,1,k)/gdt(itr)
          dpdzy1(i,k)=tc1(i,1,k)/gdt(itr)
        enddo
      enddo
    endif
  endif
  if (nclyn.eq.2) then
    if (xsize(2)==ny) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxyn(i,k)=ta1(i,ny,k)/gdt(itr)
          dpdzyn(i,k)=tc1(i,ny,k)/gdt(itr)
        enddo
      enddo
    endif
  endif

  if (nclz1.eq.2) then
    if (xstart(3)==1) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxz1(i,j)=tb1(i,j,1)/gdt(itr)
          dpdyz1(i,j)=tc1(i,j,1)/gdt(itr)
        enddo
      enddo
    endif
  endif
  if (nclzn.eq.2) then
    if (xend(3)==nz) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxzn(i,j)=tb1(i,j,xsize(3))/gdt(itr)
          dpdyzn(i,j)=tc1(i,j,xsize(3))/gdt(itr)
        enddo
      enddo
    endif
  endif

  return
end subroutine gradp
!*******************************************************************
subroutine corgp_IBM (ux,uy,uz,px,py,pz,nlock)
  USE param
  USE decomp_2d
  USE variables
  implicit none
  integer :: ijk,nlock,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz
  nxyz=xsize(1)*xsize(2)*xsize(3)
  if (nlock.eq.1) then
    if (nz.gt.1) then
      do ijk=1,nxyz
        uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1)
        uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1)
        ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
    else
      do ijk=1,nxyz
        uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1)
        ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
    endif
  endif
  if (nlock.eq.2) then
    if (nz.gt.1) then
      do ijk=1,nxyz
        uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1)
        uz(ijk,1,1)=pz(ijk,1,1)+uz(ijk,1,1)
        ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
    else
      do ijk=1,nxyz
        uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1)
        ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
    endif
  endif

  return
end subroutine corgp_IBM
!*******************************************************************
#ifdef IBM
subroutine body(ux1,uy1,uz1,ep1,arg)
  USE param
  USE decomp_2d
  USE decomp_2d_io
  USE variables
  implicit none
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  integer :: arg,ijk,nvect1

#ifdef DEBG
  if (nrank .eq. 0) print *,'# body start'
#endif


  if (arg==0) then !First execution, initt epsi
     ep1(:,:,:)=zero
     call geomcomplex(ep1,xstart(1),xend(1),ny,xstart(2),xend(2),xstart(3),xend(3),dx,yp,dz,one)
  elseif (arg==1) then  !Any other iteration
     nvect1=xsize(1)*xsize(2)*xsize(3)
     do ijk=1,nvect1
        ux1(ijk,1,1)=(one-ep1(ijk,1,1))*ux1(ijk,1,1)
        uy1(ijk,1,1)=(one-ep1(ijk,1,1))*uy1(ijk,1,1)
        uz1(ijk,1,1)=(one-ep1(ijk,1,1))*uz1(ijk,1,1)
     enddo
  endif

  !X PENCILS

#ifdef DEBG
  if (nrank .eq. 0) print *,'# body done'
#endif

  return
end subroutine body
#endif
!*******************************************************************
subroutine square(ycenter,zcenter,xthick,xlenght,ux,uy,uz,esp)
  USE param
  USE decomp_2d
  USE variables
  implicit none
  real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,esp
  real(mytype) :: ycenter,zcenter,xthick,xlenght,slenght,xz1,xz2,xy1,xy2
  integer :: j1,j2,z1,z2,iep,i,j,k,k1,k2,ilen

  iep=int(xthick*ny/yly)
  ilen=int(xlenght*ny/yly)

  j1=int(ycenter*ny/yly)
  k1=int(zcenter*nz/zlz)+ilen
  k2=int(zcenter*nz/zlz)-ilen

  do k=xstart(3),xend(3)
    do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
        if ((k.ge.(k2-iep)).and.(k.le.(k1+iep))) then
          if ((j.ge.(j1+ilen-iep)).and.(j.le.(j1+ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
        if ((k.ge.(k2-iep)).and.(k.le.(k1+iep))) then
          if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1-ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
        if ((k.ge.(k1-iep)).and.(k.le.(k1+iep))) then
          if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1+ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
        if ((k.ge.(k2-iep)).and.(k.le.(k2+iep))) then
          if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1+ilen+iep))) then
            if ((i.ge.75).and.(i.le.83)) then
              esp(i,j,k)=one
              ux(i,j,k)=zero
              uy(i,j,k)=zero
              uz(i,j,k)=zero
            endif
          endif
        endif
      enddo
    enddo
  enddo


  return
end subroutine square
!*******************************************************************
subroutine forcage_square(ux,uy,uz,esp)
  USE param
  USE decomp_2d
  USE variables
  implicit none
  real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,esp
  real(mytype),dimension(nz) :: xx1
  real(mytype),dimension(ny) :: yy1
  integer :: j, i, k, np,i1 ,ii
  real(mytype) :: ep0,ep1,ep2,ep3
  real(mytype) :: l0,l1,l2,l3,l4
  real(mytype), dimension(4) :: l5,l6
  real(mytype) :: y,z
  esp(:,:,:)=zero

  ep3=two/five !0.4  Tr=8.5 grid tmin=2.1, t2=4.2, t1=8.6, tmax=17.5
  ep2=four/five !0.8
  ep1=nine/five!
  ep0=seventeen/five!

  l0=yly/two
  l1=yly/four
  l2=yly/eight
  l3=yly/sixteen
  l4=yly/thirytwo

  !####fisrt fractal iteration##########################
  call square(l0,l0,ep0,l1,ux,uy,uz,esp)
  !####Second fractal iteration#########################
  call square(l0+l1,l0+l1,ep1,l2,ux,uy,uz,esp)
  call square(l0+l1,l1,ep1,l2,ux,uy,uz,esp)
  call square(l1,l0+l1,ep1,l2,ux,uy,uz,esp)
  call square(l1,l1,ep1,l2,ux,uy,uz,esp)
  !####third fractal iteration########################
  call square(l0-l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0-l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l1+l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0+l2,l0-l2,ep2,l3,ux,uy,uz,esp)
  call square(l0-l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l0-l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l0+l2,ep2,l3,ux,uy,uz,esp)
  call square(l2,l0+l1+l2,ep2,l3,ux,uy,uz,esp)
  !###fourth fractal iteration
  l5(1)=l3;l5(2)=-l3;l5(3)=l3;l5(4)=-l3
  l6(1)=l3;l6(2)=l3;l6(3)=-l3;l6(4)=-l3
  do ii=1,4
    call square(l0-l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0-l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l1+l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0+l2+l5(ii),l0-l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0-l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l0-l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l0+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
    call square(l2+l5(ii),l0+l1+l2+l6(ii),ep3,l4,ux,uy,uz,esp)
  enddo

  return
end subroutine forcage_square
