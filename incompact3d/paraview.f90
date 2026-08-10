#ifdef VISU

subroutine subparaview(nx, ny, nz, io, num, filename, label)

  use variables, only : prec

  implicit none

  integer, intent(in) :: nx, ny, nz, io
  character(len=*), intent(in) :: num
  character(len=*), intent(in) :: filename
  character(len=*), intent(IN), optional :: label
  character(60) :: filelabel

  if (present(label)) then
    filelabel = label
  else
    filelabel = filename
  endif

  ! write(*,*) nx
  ! write(*,*) ny
  ! write(*,*) nz
  ! write(*,*) io
  ! write(*,*) num
  ! write(*,*) filename
  ! write(*,*) filelabel
  ! write(*,*) '======= yeah ======'

  write(io,*)'            <Attribute Name="'//trim(filelabel)//'" Center="Node">'
  write(io,*)'               <DataItem Format="Binary"'
  write(io,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  write(io,*)'                Dimensions="',nz,ny,nx,'">'
  write(io,*)'                  '//filename//num
  write(io,*)'               </DataItem>'
  write(io,*)'            </Attribute>'

end subroutine subparaview

SUBROUTINE paraview()

  USE variables
  USE param
  USE decomp_2d
  USE flow_type

  IMPLICIT NONE

  INTEGER :: i,k,ifile,is,nfil=48
  real(mytype) :: xp(nx),zp(nz)
  CHARACTER(4) :: chits

  9999 FORMAT(I4.4)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview start'
#endif

  PRINT *,'Writing XDMF files'

  do i=1,nx
    xp(i) = real(i-1,mytype)*dx
  enddo
  do k=1,nz
    zp(k) = real(k-1,mytype)*dz
  enddo

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview snapshots.xdmf start'
#endif

  OPEN(nfil,file=trim(datapath)//'1_snapshots.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  if (istret.ne.0) then
    write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
    write(nfil,*)'        Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'    </Topology>'
    write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
    write(nfil,*)'        <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',xp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',yp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="',nz,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',zp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'    </Geometry>'
  else
    WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
    WRITE(nfil,*)'        Dimensions="',nz/nvisu,ny/nvisu,nx/nvisu,'">'
    WRITE(nfil,*)'    </Topology>'
    WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
    WRITE(nfil,*)'        <!-- Origin -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        0.0 0.0 0.0'
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'        <!-- DxDyDz -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        ',dz*nvisu,dy*nvisu,dx*nvisu
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'    </Geometry>'
  endif
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'            <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',REAL(imodulo*dt,4)
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  DO ifile = 0, INT(ilast/imodulo)
    WRITE(chits, 9999) ifile

    WRITE(nfil,'(/)')
    WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
    WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    IF (save_ibm.ne.0) THEN
      IF (save_ibm.eq.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,"0000","./3d-snapshot/ibm","ibm")
      IF (save_ibm.eq.2) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/ibm","ibm")
    END IF
    IF (save_ux.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/ux","ux")
    IF (save_uy.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/uy","uy")
    IF (save_uz.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/uz","uz")
    IF (iscalar.EQ.1) THEN
      DO is=1, nphi
        IF (save_phi.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/phi"//CHAR(48+is), 'phi'//CHAR(48+is))
        IF (save_dphidx.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dphi"//CHAR(48+is)//"dx", 'dphi'//CHAR(48+is)//"dx")
        IF (save_dphidy.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dphi"//CHAR(48+is)//"dy", 'dphi'//CHAR(48+is)//"dy")
        IF (save_dphidz.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dphi"//CHAR(48+is)//"dz", 'dphi'//CHAR(48+is)//"dz")
      END DO
    ENDIF
    IF (save_qc.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/qc","q-criterion")
    IF (save_pc.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/pc","p-criterion")
    IF (save_V.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/V","V")
    IF (save_pre.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/pre","pre")
    IF (save_w.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/w","w")
    IF (save_w1.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/w1","w1")
    IF (save_w2.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/w2","w2")
    IF (save_w3.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/w3","w3")
    IF (save_dudx.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dudx","du/dx")
    IF (save_dudy.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dudy","du/dy")
    IF (save_dudz.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dudz","du/dz")
    IF (save_dvdx.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dvdx","dv/dx")
    IF (save_dvdy.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dvdy","dv/dy")
    IF (save_dvdz.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dvdz","dv/dz")
    IF (save_dwdx.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dwdx","dw/dx")
    IF (save_dwdy.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dwdy","dw/dy")
    IF (save_dwdz.EQ.1) call subparaview(nx/nvisu,ny/nvisu,nz/nvisu,nfil,chits,"./3d-snapshot/dwdz","dw/dz")
    WRITE(nfil,*)'        </Grid>'
  ENDDO
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview snapshots.xdmf done'
  IF (nrank .EQ. 0) PRINT *,'# paraview planes.xdmf start'
#endif

  !###############################################################################
  OPEN(nfil,file=trim(datapath)//'2_mean-planes-z.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  if (istret.ne.0) then
    write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
    write(nfil,*)'        Dimensions="',1,ny,nx,'">'
    write(nfil,*)'    </Topology>'
    write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
    write(nfil,*)'        <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',xp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',yp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        0.0'
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'    </Geometry>'
  else
    WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
    WRITE(nfil,*)'        Dimensions="',1,ny,nx,'">'
    WRITE(nfil,*)'    </Topology>'
    WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
    WRITE(nfil,*)'        <!-- Origin -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        0.0 0.0 0.0'
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'        <!-- DxDyDz -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        ',0.0,dy,dx
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'    </Geometry>'
  endif
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'           <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',REAL(imodulo*dt,4)
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  DO ifile = 0, INT(ilast/imodulo)
    WRITE(chits, 9999) ifile

    WRITE(nfil,'(/)')
    WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
    WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

    ! IF (iscalar.EQ.1.AND.save_phim.EQ.1) THEN
    !   DO is=1, nphi
    !     call subparaview(nx,ny,1,nfil,chits,"../xy-plane/phim"//CHAR(48+is), 'phim'//CHAR(48+is))
    !   END DO
    ! ENDIF
    ! IF (save_uxm.EQ.1) call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uxm", 'uxm')
    ! IF (save_uym.EQ.1) call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uym", 'uym')
    ! IF (save_uzm.EQ.1) call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uzm", 'uzm')
    ! IF (save_prem.EQ.1) call subparaview(nx,ny,1,nfil,chits,"./xy-plane/prem", 'prem')

    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/usum", 'usum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vsum", 'vsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/wsum", 'wsum')
    do is=1, nphi
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phi"//char(is+48)//"sum",&
      "phi"//char(is+48)//"sum")
    enddo
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phisum", 'phisum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uusum", 'uusum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uvsum", 'uvsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uwsum", 'uwsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vvsum", 'vvsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vwsum", 'vwsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/wwsum", 'wwsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uphisum", 'uphisum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vphisum", 'vphisum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/wphisum", 'wphisum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phiphisum", 'phiphisum')
    !
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/psum", 'psum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/ppsum", 'ppsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/pusum", 'pusum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/pvsum", 'pvsum')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/pwsum", 'pwsum')
    !
    ! call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dudx", 'dudx')
    ! call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dudy", 'dudy')
    ! call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dvdx", 'dvdx')
    ! call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dvdy", 'dvdy')
    ! call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dphidx", 'dphidx')
    ! call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dphidy", 'dphidy')

    WRITE(nfil,*)'        </Grid>'
  ENDDO

  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview planes.xdmf done'
  IF (nrank .EQ. 0) PRINT *,'# paraview stats_sum.xdmf start'
#endif

  !###############################################################################
  IF(save_dmap.EQ.1.OR.save_utmap.EQ.1) THEN

    OPEN(nfil,file=trim(datapath)//'3_maps.xdmf')
    WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
    WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    WRITE(nfil,*)'<Domain>'
    WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
    WRITE(nfil,*)'        Dimensions="',nz,1,nx,'">'
    WRITE(nfil,*)'    </Topology>'
    WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
    WRITE(nfil,*)'        <!-- Origin -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        0.0 0.0 0.0'
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'        <!-- DxDyDz -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        ',dz,0.0,dx
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'    </Geometry>'
    WRITE(nfil,'(/)')
    WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
    WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
    WRITE(nfil,*)'           <!--Start, Stride, Count-->'
    WRITE(nfil,*)'            0.0',REAL(imodulo*dt,4)
    WRITE(nfil,*)'            </DataItem>'
    WRITE(nfil,*)'        </Time>'

    DO ifile = 0, INT(ilast/imodulo)
      WRITE(chits, 9999) ifile

      WRITE(nfil,'(/)')
      WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
      WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

      IF (iscalar.EQ.1) THEN
        DO is=1, nphi
          call subparaview(nx,1,nz,nfil,chits,"./xz-plane/phif"//CHAR(48+is), 'phif'//CHAR(48+is))
          IF (save_dmap.EQ.1) call subparaview(nx,1,nz,nfil,chits,"./xz-plane/dep"//CHAR(48+is), 'dep'//CHAR(48+is))
        ENDDO
      ENDIF

      IF (save_utmap.EQ.1) call subparaview(nx,1,nz,nfil,chits,"./xz-plane/utmap", 'utmap')

      WRITE(nfil,*)'        </Grid>'
    ENDDO
    WRITE(nfil,'(/)')
    WRITE(nfil,*)'    </Grid>'
    WRITE(nfil,*)'</Domain>'
    WRITE(nfil,'(A7)')'</Xdmf>'
    CLOSE(nfil)

  ENDIF

  !###############################################################################
#ifdef STATS

  OPEN(nfil,file=trim(datapath)//'4_stats_sum.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
  WRITE(nfil,*)'        Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'    </Topology>'
  WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
  WRITE(nfil,*)'        <!-- Origin -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        0.0 0.0 0.0'
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'        <!-- DxDyDz -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        ',dz*nstat,dy*nstat,dx*nstat
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'    </Geometry>'
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'           <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',ilast
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  WRITE(nfil,'(/)')
  WRITE(nfil,*)'        <Grid Name="" GridType="Uniform">'
  WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
  WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

  IF (save_pre.EQ.1) THEN
    WRITE(nfil,*)'            <Attribute Name="pre.sum" Center="Node">'
    WRITE(nfil,*)'               <DataItem Format="Binary" '
    WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
    WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
    WRITE(nfil,*)'                  ./3d-snapshot/pre.sum'
    WRITE(nfil,*)'               </DataItem>'
    WRITE(nfil,*)'            </Attribute>'
  ENDIF

  WRITE(nfil,*)'            <Attribute Name="c1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/c1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="c1c1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/c1c1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="uc1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/uc.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="vc1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/vc1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="wc1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/wc1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="u1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/u1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/v1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w1.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/w1.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="u2.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/u2.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v2.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/v2.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w2.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/w2.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="uv.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/uv.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="uw.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/uv.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="vw.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/vw.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="k.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/k.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  IF (save_dudx.EQ.0) THEN
    WRITE(nfil,*)'            <Attribute Name="dudx.sum" Center="Node">'
    WRITE(nfil,*)'               <DataItem Format="Binary" '
    WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
    WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
    WRITE(nfil,*)'                  ./3d-snapshot/dudx.sum'
    WRITE(nfil,*)'               </DataItem>'
    WRITE(nfil,*)'            </Attribute>'
  ENDIF

  WRITE(nfil,*)'            <Attribute Name="u3.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/u3.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v3.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/v3.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w3.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/w3.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="u4.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/u4.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="v4.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/v4.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'            <Attribute Name="w4.sum" Center="Node">'
  WRITE(nfil,*)'               <DataItem Format="Binary" '
  WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
  WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'                  ./3d-snapshot/w4.sum'
  WRITE(nfil,*)'               </DataItem>'
  WRITE(nfil,*)'            </Attribute>'

  WRITE(nfil,*)'        </Grid>'
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview stats_sum.xdmf done'
  IF (nrank .EQ. 0) PRINT *,'# paraview stats_mean.xdmf start'
#endif

  !###############################################################################

  OPEN(nfil,file=trim(datapath)//'5_stats_mean_convergence.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
  WRITE(nfil,*)'        Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
  WRITE(nfil,*)'    </Topology>'
  WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
  WRITE(nfil,*)'        <!-- Origin -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        0.0 0.0 0.0'
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'        <!-- DxDyDz -->'
  WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
  WRITE(nfil,*)'        ',dz*nstat,dy*nstat,dx*nstat
  WRITE(nfil,*)'        </DataItem>'
  WRITE(nfil,*)'    </Geometry>'
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'           <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',REAL(ntik*dt,4)
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  DO ifile = 0, INT(ilast/imodulo)
    WRITE(chits, 9999) ifile

    WRITE(nfil,'(/)')
    WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
    WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

    WRITE(nfil,*)'            <Attribute Name="ux" Center="Node">'
    WRITE(nfil,*)'               <DataItem Format="Binary" '
    WRITE(nfil,*)'                DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
    WRITE(nfil,*)'                Dimensions="',nz/nstat,ny/nstat,nx/nstat,'">'
    WRITE(nfil,*)'                  ./mean/uxmean'//chits
    WRITE(nfil,*)'               </DataItem>'
    WRITE(nfil,*)'            </Attribute>'

    WRITE(nfil,*)'        </Grid>'
  ENDDO
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

#ifdef DEBG
  IF (nrank .EQ. 0) PRINT *,'# paraview stats_mean.xdmf end'
  IF (nrank .EQ. 0) PRINT *,'# paraview end'
#endif
#endif

#ifdef POST
  !#######################################################################################
  OPEN(nfil,file=trim(datapath)//'xy-planes.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  if (istret.ne.0) then
    write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
    write(nfil,*)'        Dimensions="',1,ny,nx,'">'
    write(nfil,*)'    </Topology>'
    write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
    write(nfil,*)'        <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',xp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',yp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        0.0'
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'    </Geometry>'
  else
    WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
    WRITE(nfil,*)'        Dimensions="',1,ny,nx,'">'
    WRITE(nfil,*)'    </Topology>'
    WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
    WRITE(nfil,*)'        <!-- Origin -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        0.0 0.0 0.0'
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'        <!-- DxDyDz -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        ',0.0,dy,dx
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'    </Geometry>'
  endif
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'           <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',REAL(iprocessing*dt,4)
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  DO ifile = 0, INT(ilast/iprocessing)
    WRITE(chits, 9999) ifile

    WRITE(nfil,'(/)')
    WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
    WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

    DO is=1, nphi
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phi"//CHAR(48+is), 'phi'//CHAR(48+is)//"m")
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phi"//CHAR(48+is)//'c', 'phi'//CHAR(48+is)//"c")
    END DO
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/ux", 'uxm')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uy", 'uym')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uz", 'uzm')
    !
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uxc", 'uxc')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uyc", 'uyc')
    call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uzc", 'uzc')

    if (iin.eq.2) then
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/usum", 'usum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vsum", 'vsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/wsum", 'wsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phisum", 'phisum')
      do is=1, nphi
        call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phi"//char(is+48)//"sum",&
        "phi"//char(is+48)//"sum")
      enddo
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uusum", 'uusum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uvsum", 'uvsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uwsum", 'uwsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vvsum", 'vvsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vwsum", 'vwsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/wwsum", 'wwsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/uphisum", 'uphisum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/vphisum", 'vphisum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/wphisum", 'wphisum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/phiphisum", 'phiphisum')
      !
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dudx", 'dudx')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dudy", 'dudy')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dvdx", 'dvdx')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dvdy", 'dvdy')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dphidx", 'dphidx')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/dphidy", 'dphidy')
      !
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/psum", 'psum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/ppsum", 'ppsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/pusum", 'pusum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/pvsum", 'pvsum')
      call subparaview(nx,ny,1,nfil,chits,"./xy-plane/pwsum", 'pwsum')
    endif

    WRITE(nfil,*)'        </Grid>'
  ENDDO

  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)

  !#######################################################################################
  OPEN(nfil,file=trim(datapath)//'xz-planes.xdmf')
  WRITE(nfil,'(A22)')'<?xml version="1.0" ?>'
  WRITE(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  WRITE(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  WRITE(nfil,*)'<Domain>'
  if (istret.ne.0) then
    write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
    write(nfil,*)'        Dimensions="',nz,1,nx,'">'
    write(nfil,*)'    </Topology>'
    write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
    write(nfil,*)'        <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',xp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        0.0'
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'        <DataItem Dimensions="',nz,'" NumberType="Float" Precision="4" Format="XML">'
    write(nfil,*)'        ',zp(:)
    write(nfil,*)'        </DataItem>'
    write(nfil,*)'    </Geometry>'
  else
    WRITE(nfil,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
    WRITE(nfil,*)'        Dimensions="',nz,1,nx,'">'
    WRITE(nfil,*)'    </Topology>'
    WRITE(nfil,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
    WRITE(nfil,*)'        <!-- Origin -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        0.0 0.0 0.0'
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'        <!-- DxDyDz -->'
    WRITE(nfil,*)'        <DataItem Format="XML" Dimensions="3">'
    WRITE(nfil,*)'        ',dz,0.0,dx
    WRITE(nfil,*)'        </DataItem>'
    WRITE(nfil,*)'    </Geometry>'
  endif
  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  WRITE(nfil,*)'        <Time TimeType="HyperSlab">'
  WRITE(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  WRITE(nfil,*)'           <!--Start, Stride, Count-->'
  WRITE(nfil,*)'            0.0',REAL(iprocessing*dt,4)
  WRITE(nfil,*)'            </DataItem>'
  WRITE(nfil,*)'        </Time>'

  DO ifile = 0, INT(ilast/iprocessing)
    WRITE(chits, 9999) ifile

    WRITE(nfil,'(/)')
    WRITE(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
    WRITE(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    WRITE(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

    call subparaview(nx,1,nz,nfil,chits,"./xz-plane/ux", 'uxm')
    call subparaview(nx,1,nz,nfil,chits,"./xz-plane/uy", 'uym')
    call subparaview(nx,1,nz,nfil,chits,"./xz-plane/uz", 'uzm')
    !
    call subparaview(nx,1,nz,nfil,chits,"./xz-plane/uxt", 'uxt')
    call subparaview(nx,1,nz,nfil,chits,"./xz-plane/uyt", 'uyt')
    call subparaview(nx,1,nz,nfil,chits,"./xz-plane/uzt", 'uzt')
    !
    call subparaview(nx,1,nz,nfil,chits,"./xz-plane/utau", 'utau')
    !
    DO is=1, nphi
      call subparaview(nx,1,nz,nfil,chits,"./xz-plane/phi"//CHAR(48+is)//"t", 'phi'//CHAR(48+is)//"t")
      call subparaview(nx,1,nz,nfil,chits,"./xz-plane/phi"//CHAR(48+is), 'phi'//CHAR(48+is)//"m")
      if (uset(is) .ne. 0.) then
        call subparaview(nx,1,nz,nfil,chits,"./xz-plane/dep"//CHAR(48+is), 'dep'//CHAR(48+is))
      endif
    END DO

    WRITE(nfil,*)'        </Grid>'
  ENDDO

  WRITE(nfil,'(/)')
  WRITE(nfil,*)'    </Grid>'
  WRITE(nfil,*)'</Domain>'
  WRITE(nfil,'(A7)')'</Xdmf>'
  CLOSE(nfil)
#endif

  return
END SUBROUTINE paraview
#endif
