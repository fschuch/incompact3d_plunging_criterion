PROGRAM post

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  USE post_processing

  implicit none

  integer :: code,i,j,k,is,ii

  integer :: ifile, nt, iwrn, num, ttsize
  integer :: read_phi, read_u, read_ibm     !! 3D fields
  real(8) :: tstart,t1,trank,tranksum,ttotal,tremaining,telapsed,trstart, trend
  character(120) :: filename

  call ft_parameter(0)

  CALL MPI_INIT(code)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.) !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.) !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  call parameter()
  call init_variables
  call schemes(.True.)

  ux1=zero; uxm1=zero
  uy1=zero; uym1=zero
  uz1=zero;
  phi1=zero; phim1=zero
  diss1=zero; dissm1=zero
  pre1=zero; prem1=zero

  read_phi=0; read_u=0; read_ibm=0

  nt = (filen-file1)/icrfile+1

  if(comp_post .eq. 1) then
     read_phi=1; read_ibm=1; read_u=1
  endif

  if(comp_visu .eq. 1) then
     call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
     read_phi=1; read_u=1; read_ibm=1
  endif
  if ( iscalar.eq.0) read_phi=0
  if ( ivirt  .eq.0) read_ibm=0
  if ( read_ibm .eq. 1 ) then
     write(filename,"('ibm',I4.4)") 0 !itime/imodulo
     filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
     call decomp_2d_read_one(1,ep1,filename)
  endif

  call init_post(ep1)

  ttsize=(read_phi*nphi+read_u*3)*nx*ny*nz
  tstart=0.;t1=0.;trank=0.;tranksum=0.;ttotal=0.
  call cpu_time(tstart)

  do ii=1, nt
     call cpu_time(t1)
     ifile = (ii-1)*icrfile+file1
     write(filename,"(I4.4)") ifile
     t=dt*real(imodulo*ifile,mytype)
     itime=imodulo*ifile
     if (nrank==0) then
        print *,'--------------------------------------------'
        print *,'Snapshot',ifile, t
     endif

     call ft_parameter(2)

     !READ DATA
     call cpu_time(trstart)
     if ( read_phi .eq. 1 ) then
        do is=1, nphi
           write(filename,"('phi',I1.1,I4.4)") is, ifile
           filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
           call decomp_2d_read_one(1,phi1(:,:,:,is),filename)
        enddo
        call test_scalar_min_max(phi1)
     endif
     if ( read_u .eq. 1 ) then
        write(filename,"('ux',I4.4)") itime/imodulo
        filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
        call decomp_2d_read_one(1,ux1,filename)
        !
        write(filename,"('uy',I4.4)") itime/imodulo
        filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
        call decomp_2d_read_one(1,uy1,filename)
        !
        write(filename,"('uz',I4.4)") itime/imodulo
        filename = trim(datapath)//'/3d-snapshot/'//trim(filename)
        call decomp_2d_read_one(1,uz1,filename)
        !
        call test_speed_min_max(ux1,uy1,uz1)
     endif
     call cpu_time(trend)

#ifdef VISU
     if (comp_visu .eq. 1) then
        call VISU_INSTA(ux1,uy1,uz1,phi1,ep1,.True.)
     endif
#endif
#ifdef POST
     if (comp_post .eq. 1) then
        call postprocessing(ux1,uy1,uz1,phi1,ep1)
     endif
#endif
     call cpu_time(trank)

     telapsed = (trank-tstart)/3600.
     tremaining  =  telapsed*(nt-ii)/(ii)

     if (nrank==0) then
        print *,'Time per this snapshot (s):',real(trank-t1)
        print *,'Reading speed (MB/s)',real((ttsize*1e-6)/(trend-trstart),4)
        write(*,"(' Remaining time:',I8,' h ',I2,' min')"), int(tremaining), int((tremaining-int(tremaining))*60.)
        write(*,"(' Elapsed time:',I8,' h ',I2,' min')"), int(telapsed), int((telapsed-int(telapsed))*60.)
     endif
  enddo

  call cpu_time(trank)
  ttotal=trank-tstart

  if (nrank==0) then
     print *,'==========================================================='
     print *,''
     print *,'Post-processing finished successfully!'
     print *,''
     print *,'2DECOMP with p_row*p_col=',p_row,p_col
     print *,''
     print *,'nx*ny*nz=',nx*ny*nz
     print *,'nx,ny,nz=',nx,ny,nz
     print *,'dx,dy,dz=',dx,dy,dz
     print *,''
     print *,'Averaged time per snapshot (s):',real(ttotal/nt,4)
     print *,'Total wallclock (s):',real(ttotal,4)
     print *,'Total wallclock (m):',real(ttotal/60.,4)
     print *,'Total wallclock (h):',real(ttotal/3600.,4)
     print *,'Total wallclock (d):',real(ttotal*1.1574e-5,4)
     print *,''
  endif

  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)
end PROGRAM post
