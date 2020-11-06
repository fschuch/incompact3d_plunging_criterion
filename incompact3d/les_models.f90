#ifdef ELES

subroutine init_explicit_les

  USE param
  USE variables
  USE decomp_2d

  implicit none

  integer :: j

  !you can specify other metric e.g. max(dx,dy,dz) or min(dx,dy,dz)

  do j = 1, ny - 1
     del(j) = (dx * (yp(j + 1) - yp(j)) * dz)**(one / three)
  enddo
  del(ny) = del(ny - 1)

  if(nrank==0) then
     write(*, *) ' '
     write(*, *) '++++++++++++++++++++++++++++++++'
     write(*, *) 'Initializing explicit LES Filter'
     if(jLES==2) then
        write(*, *) ' Classic Smagorinsky is used ... '
        write(*, *) ' Smagorinsky constant = ', smagcst
        write(*, *) ' Filter Size / Grid Size = ', FSGS
     else if (jLES==3) then
        write(*, *) ' Wall-adaptive LES (WALES) is used ... '
        write(*, *) ' WALE constant = ', walecst
     else if (jLES==4) then
        write(*, *) ' Dynamic Smagorinsky is used ... '
        write(*, *) ' Max value for the dynamic constant field = ', maxdsmagcst
     endif
     write(*, *) '++++++++++++++++++++++++++++++++'
     if (nrank==0) print *, "Del y min max= ", minval(del), maxval(del)
     write(*, *) '++++++++++++++++++++++++++++++++'
     write(*, *) ' '
  endif

end subroutine init_explicit_les
!************************************************************
subroutine smag(ux1, uy1, uz1, gxx1, gyx1, gzx1, gxy1, gyy1, gzy1, gxz1, gyz1, gzz1, &
     sxx1, syy1, szz1, sxy1, sxz1, syz1, srt_smag, nut1, ta2, ta3, di1, di2, di3)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2, uy2, uz2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ux3, uy3, uz3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: di1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ta2, di2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ta3, di3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sxx1, syy1, szz1, sxy1, sxz1, syz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: gxx1, gyx1, gzx1, gxy1, gyy1, gzy1, gxz1, gyz1, gzz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: srt_smag, nut1

  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: syy2, szz2, sxy2, syz2, srt_smag2, nut2
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: gxy2, gyy2, gzy2, gxz2, gyz2, gzz2

  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: szz3, syz3
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: gxz3, gyz3, gzz3

  integer :: i, j, k
  character(len = 30) :: filename

  call derx (gxx1, ux1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0)
  call derx (gyx1, uy1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)
  call derx (gzx1, uz1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)

  sxx1 = gxx1

  !WORK Y-PENCILS
  call transpose_x_to_y(ux1, ux2)
  call transpose_x_to_y(uy1, uy2)
  call transpose_x_to_y(uz1, uz2)
  call transpose_x_to_y(gyx1, ta2)

  call dery (gxy2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)
  call dery (gyy2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)
  call dery (gzy2, uz2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)

  sxy2 = half * (gxy2 + ta2)
  syy2 = gyy2

  !WORK Z-PENCILS
  call transpose_y_to_z(ux2, ux3)
  call transpose_y_to_z(uy2, uy3)
  call transpose_y_to_z(uz2, uz3)
  call transpose_y_to_z(gzy2, ta3)

  call derz(gxz3, ux3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
  call derz(gyz3, uy3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
  call derz(gzz3, uz3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0)

  szz3 = gzz3
  syz3 = half * (gyz3 + ta3)

  !WORK Y-PENCILS
  call transpose_z_to_y(syz3, syz2)
  call transpose_z_to_y(szz3, szz2)

  call transpose_z_to_y(gxz3, gxz2)
  call transpose_z_to_y(gyz3, gyz2)
  call transpose_z_to_y(gzz3, gzz2)

  !WORK X-PENCILS
  call transpose_y_to_x(sxy2, sxy1)
  call transpose_y_to_x(syy2, syy1)
  call transpose_y_to_x(syz2, syz1)
  call transpose_y_to_x(szz2, szz1)

  call transpose_y_to_x(gxy2, gxy1)
  call transpose_y_to_x(gyy2, gyy1)
  call transpose_y_to_x(gzy2, gzy1)
  call transpose_y_to_x(gxz2, gxz1)
  call transpose_y_to_x(gyz2, gyz1)
  call transpose_y_to_x(gzz2, gzz1)

  sxz1 = half * (gzx1 + gxz1)

  srt_smag = zero
  srt_smag = sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1

  nut1 = zero; nut2 = zero
  call transpose_x_to_y(srt_smag, srt_smag2)
  do k = 1, ysize(3)
     do j = 1, ysize(2)
        do i = 1, ysize(1)
           nut2(i, j, k) = ((smagcst * FSGS * del(j))**two) * sqrt(two * srt_smag2(i, j, k))
        enddo
     enddo
  enddo
  call transpose_y_to_x(nut2, nut1)

  if (nrank==0) print *, "smag srt_smag min max= ", minval(srt_smag), maxval(srt_smag)
  if (nrank==0) print *, "smag nut1     min max= ", minval(nut1), maxval(nut1)

  if (mod(itime, imodulo).eq.0) then

     write(filename, "('./data/srt_smag',I4.4)") itime / imodulo
     call decomp_2d_write_one(1, srt_smag, filename, 2)

     write(filename, "('./data/nut_smag',I4.4)") itime / imodulo
     call decomp_2d_write_one(1, nut1, filename, 2)

  endif

end subroutine smag
!************************************************************
subroutine wale(gxx1, gyx1, gzx1, gxy1, gyy1, gzy1, gxz1, gyz1, gzz1, srt_smag, nut1)


  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: gxx1, gyx1, gzx1, &
       gxy1, gyy1, gzy1, gxz1, gyz1, gzz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: srt_wale, srt_smag
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: srt_wale2, srt_smag2, nut2
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sxxd1, syyd1, szzd1, sxyd1, sxzd1, syzd1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: nut1

  character(len = 30) :: filename
  integer :: i, j, k

  sxxd1 = gxx1 * gxx1 + gxy1 * gyx1 + gxz1 * gzx1 - (one / three) * (gxx1 * gxx1 + gyy1 * gyy1 + gzz1 * gzz1 + two * gxy1 * gyx1 + two * gxz1 * gzx1 + two * gzy1 * gyz1)
  sxyd1 = half * (gxx1 * gxy1 + gxy1 * gyy1 + gxz1 * gzy1 + gyx1 * gxx1 + gyy1 * gyx1 + gyz1 * gzx1)
  sxzd1 = half * (gxx1 * gxz1 + gxy1 * gyz1 + gxz1 * gzz1 + gzx1 * gxx1 + gzy1 * gyx1 + gzz1 * gzx1)
  syyd1 = gyx1 * gxy1 + gyy1 * gyy1 + gyz1 * gzy1 - (one / three) * (gxx1 * gxx1 + gyy1 * gyy1 + gzz1 * gzz1 + two * gxy1 * gyx1 + two * gxz1 * gzx1 + two * gzy1 * gyz1)
  syzd1 = half * (gyx1 * gxz1 + gyy1 * gyz1 + gyz1 * gzz1 + gzx1 * gxy1 + gzy1 * gyy1 + gzz1 * gzy1)
  szzd1 = gzx1 * gxz1 + gzy1 * gyz1 + gzz1 * gzz1 - (one / three) * (gxx1 * gxx1 + gyy1 * gyy1 + gzz1 * gzz1 + two * gxy1 * gyx1 + two * gxz1 * gzx1 + two * gzy1 * gyz1)

  srt_wale = zero
  srt_wale = sxxd1 * sxxd1 + syyd1 * syyd1 + szzd1 * szzd1 + two * sxyd1 * sxyd1 + two * sxzd1 * sxzd1 + two * syzd1 * syzd1

  nut1 = zero; nut2 = zero
  call transpose_x_to_y(srt_smag, srt_smag2)
  call transpose_x_to_y(srt_wale, srt_wale2)
  do k = 1, ysize(3)
     do j = 1, ysize(2)
        do i = 1, ysize(1)
           nut2(i, j, k) = ((walecst * del(j))**two) * ((srt_wale2(i, j, k)**(three / two)) / ((srt_smag2(i, j, k)**(five / two)) + (srt_wale2(i, j, k)**(five / four))))
        enddo
     enddo
  enddo
  call transpose_y_to_x(nut2, nut1)

  if (nrank==0) print *, "wale srt_wale min max= ", minval(srt_wale), maxval(srt_wale)
  if (nrank==0) print *, "wale nut1     min max= ", minval(nut1), maxval(nut1)

  if (mod(itime, imodulo).eq.0) then

     write(filename, "('./data/srt_wale',I4.4)") itime / imodulo
     call decomp_2d_write_one(1, srt_wale, filename, 2)

     write(filename, "('./data/nut_wale',I4.4)") itime / imodulo
     call decomp_2d_write_one(1, nut1, filename, 2)

  endif

end subroutine wale
!************************************************************
subroutine dynsmag(ux1,uy1,uz1,ep1,sxx1,syy1,szz1,sxy1,sxz1,syz1,&
     srt_smag,nut1,di1,ta1,tb1,tc1,td1,ta2,tb2,tc2,td2,te2,tf2,&
     tg2,th2,ti2,di2,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE MPI

  implicit none

  TYPE(DECOMP_INFO) :: ph3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1, srt_smag
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: nut1

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1f, uy1f, uz1f
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2f, uy2f, uz2f, ep2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ux3f, uy3f, uz3f

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: uxx1, uyy1, uzz1, uxy1, uxz1, uyz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: uxx1f, uyy1f, uzz1f, uxy1f, uxz1f, uyz1f
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: uxx2f, uyy2f, uzz2f, uxy2f, uxz2f, uyz2f
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: uxx3f, uyy3f, uzz3f, uxy3f, uxz3f, uyz3f

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ta1, tb1, tc1, td1, di1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, di2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ta3, tb3, tc3, td3, te3, tf3, tg3, th3, ti3, di3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sxx1, syy1, szz1, sxy1, sxz1, syz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sxx1f, syy1f, szz1f, sxy1f, sxz1f, syz1f
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: syy2f, szz2f, sxy2f, syz2f
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: szz3f, syz3f

  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: nut2, srt_smag2

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: gxx1f, gyx1f, gzx1f, gxz1f
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: gxy2f, gyy2f, gzy2f, gxz2f
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: gxz3f, gyz3f, gzz3f

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: axx1, ayy1, azz1, axy1, axz1, ayz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: axx1f, ayy1f, azz1f, axy1f, axz1f, ayz1f
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: axx2f, ayy2f, azz2f, axy2f, axz2f, ayz2f
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: axx3f, ayy3f, azz3f, axy3f, axz3f, ayz3f

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: bbxx1, bbyy1, bbzz1, bbxy1, bbxz1, bbyz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: lxx1, lyy1, lzz1, lxy1, lxz1, lyz1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mxx1, myy1, mzz1, mxy1, mxz1, myz1

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: smagC1, smagC1f, dsmagcst1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: smagC2, smagC2f, dsmagcst2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: smagC3, smagC3f, dsmagcst3

  character(len = 30) :: filename

  integer :: i, j, k, ijk, code, nvect1
  nvect1 = xsize(1) * xsize(2) * xsize(3)

#ifdef IBM
  ta1 = ux1 * (one - ep1)
  tb1 = uy1 * (one - ep1)
  tc1 = uz1 * (one - ep1)
#else
  ta1 = ux1
  tb1 = uy1
  tc1 = uz1
#endif

  !    !!!average TA1
  !    call transpose_x_to_y(ta1,ta2)
  !    call transpose_y_to_z(ta2,ta3)
  !    dsmagcst3 = zero
  !    call mean_plane_z(ta3,zsize(1),zsize(2),zsize(3),dsmagcst3(:,:,1))
  !    do k=2,zsize(3)
  !      dsmagcst3(:,:,k) = dsmagcst3(:,:,1)
  !    enddo
  !    call transpose_z_to_y(dsmagcst3,ta2)
  !    ta1=zero
  !    call transpose_y_to_x(ta2,ta1)

  !    !!!average TB1
  !    call transpose_x_to_y(tb1,tb2)
  !    call transpose_y_to_z(tb2,tb3)
  !    dsmagcst3 = zero
  !    call mean_plane_z(tb3,zsize(1),zsize(2),zsize(3),dsmagcst3(:,:,1))
  !    do k=2,zsize(3)
  !      dsmagcst3(:,:,k) = dsmagcst3(:,:,1)
  !    enddo
  !    call transpose_z_to_y(dsmagcst3,tb2)
  !    tb1=zero
  !    call transpose_y_to_x(tb2,tb1)

  !    !!!average TC1
  !    call transpose_x_to_y(tc1,tc2)
  !    call transpose_y_to_z(tc2,tc3)
  !    dsmagcst3 = zero
  !    call mean_plane_z(tc3,zsize(1),zsize(2),zsize(3),dsmagcst3(:,:,1))
  !    do k=2,zsize(3)
  !      dsmagcst3(:,:,k) = dsmagcst3(:,:,1)
  !    enddo
  !    call transpose_z_to_y(dsmagcst3,tc2)
  !    tc1=zero
  !    call transpose_y_to_x(tc2,tc1)

  uxx1 = ta1 * ta1
  uyy1 = tb1 * tb1
  uzz1 = tc1 * tc1
  uxy1 = ta1 * tb1
  uxz1 = ta1 * tc1
  uyz1 = tb1 * tc1

  if (nclx1.eq.0.and.nclxn.eq.0) then

     call filx(ux1f, ta1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0) !ux1
     call filx(uy1f, tb1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uy1
     call filx(uz1f, tc1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uz1

     call filx(uxx1f, uxx1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !ux1*ux1
     call filx(uyy1f, uyy1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uy1*uy1
     call filx(uzz1f, uzz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uz1*uz1

     call filx(uxy1f, uxy1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0) !ux1*uy1
     call filx(uxz1f, uxz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0) !ux1*uz1
     call filx(uyz1f, uyz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uy1*uz1

  else

     call filx(ux1f, ta1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0) !ux1
     call filx(uy1f, tb1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uy1
     call filx(uz1f, tc1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uz1

     call filx(uxx1f, uxx1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !ux1*ux1
     call filx(uyy1f, uyy1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uy1*uy1
     call filx(uzz1f, uzz1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uz1*uz1

     call filx(uxy1f, uxy1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0) !ux1*uy1
     call filx(uxz1f, uxz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0) !ux1*uz1
     call filx(uyz1f, uyz1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1) !uy1*uz1

  endif

  if (nrank==0) print *, "filx ux= ", maxval(ta1), maxval(ux1f), maxval(ta1) - maxval(ux1f)

#ifdef IBM
  ux1f = ux1f * (one - ep1)
  uy1f = uy1f * (one - ep1)
  uz1f = uz1f * (one - ep1)
  uxx1f = uxx1f * (one - ep1)
  uyy1f = uyy1f * (one - ep1)
  uzz1f = uzz1f * (one - ep1)
  uxy1f = uxy1f * (one - ep1)
  uxz1f = uxz1f * (one - ep1)
  uyz1f = uyz1f * (one - ep1)
  call transpose_x_to_y(ep1, ep2)
#endif

  call transpose_x_to_y(ux1f, ta2)
  call transpose_x_to_y(uy1f, tb2)
  call transpose_x_to_y(uz1f, tc2)
  call transpose_x_to_y(uxx1f, td2)
  call transpose_x_to_y(uyy1f, te2)
  call transpose_x_to_y(uzz1f, tf2)
  call transpose_x_to_y(uxy1f, tg2)
  call transpose_x_to_y(uxz1f, th2)
  call transpose_x_to_y(uyz1f, ti2)

  if (ncly1.eq.0.and.nclyn.eq.0) then

     call fily(ux2f, ta2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !ux2
     call fily(uy2f, tb2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0) !uy2
     call fily(uz2f, tc2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !uz2

     call fily(uxx2f, td2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !ux2*ux2
     call fily(uyy2f, te2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !uy2*uy2
     call fily(uzz2f, tf2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !uz2*uz2

     call fily(uxy2f, tg2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0) !ux2*uy2
     call fily(uxz2f, th2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !ux2*uz2
     call fily(uyz2f, ti2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0) !uy2*uz2

  else

     call fily(ux2f, ta2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !ux2
     call fily(uy2f, tb2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0) !uy2
     call fily(uz2f, tc2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !uz2

     call fily(uxx2f, td2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !ux2*ux2
     call fily(uyy2f, te2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !uy2*uy2
     call fily(uzz2f, tf2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !uz2*uz2

     call fily(uxy2f, tg2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0) !ux2*uy2
     call fily(uxz2f, th2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1) !ux2*uz2
     call fily(uyz2f, ti2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0) !uy2*uz2

  endif

  if (nrank==0) print *, "fily ux= ", maxval(ta2), maxval(ux2f), maxval(ta2) - maxval(ux2f)

#ifdef IBM
  ux2f = ux2f * (one - ep2)
  uy2f = uy2f * (one - ep2)
  uz2f = uz2f * (one - ep2)
  uxx2f = uxx2f * (one - ep2)
  uyy2f = uyy2f * (one - ep2)
  uzz2f = uzz2f * (one - ep2)
  uxy2f = uxy2f * (one - ep2)
  uxz2f = uxz2f * (one - ep2)
  uyz2f = uyz2f * (one - ep2)
#endif

  ta2 = zero; tb2 = zero; tc2 = zero
  td2 = zero; te2 = zero; tf2 = zero
  tg2 = zero; th2 = zero; ti2 = zero

  call transpose_y_to_z(ux2f, ta3)
  call transpose_y_to_z(uy2f, tb3)
  call transpose_y_to_z(uz2f, tc3)
  call transpose_y_to_z(uxx2f, td3)
  call transpose_y_to_z(uyy2f, te3)
  call transpose_y_to_z(uzz2f, tf3)
  call transpose_y_to_z(uxy2f, tg3)
  call transpose_y_to_z(uxz2f, th3)
  call transpose_y_to_z(uyz2f, ti3)

  if (nclz1.eq.0.and.nclzn.eq.0) then

     call filz(ux3f, ta3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !ux3
     call filz(uy3f, tb3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !uy3
     call filz(uz3f, tc3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0) !uz3

     call filz(uxx3f, td3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !ux3*ux3
     call filz(uyy3f, te3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !uy3*uy3
     call filz(uzz3f, tf3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !uz3*uz3

     call filz(uxy3f, tg3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !ux3*uy3
     call filz(uxz3f, th3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0) !ux3*uz3
     call filz(uyz3f, ti3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0) !uy3*uz3

  else

     call filz(ux3f, ta3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !ux3
     call filz(uy3f, tb3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !uy3
     call filz(uz3f, tc3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0) !uz3

     call filz(uxx3f, td3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !ux3*ux3
     call filz(uyy3f, te3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !uy3*uy3
     call filz(uzz3f, tf3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !uz3*uz3

     call filz(uxy3f, tg3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1) !ux3*uy3
     call filz(uxz3f, th3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0) !ux3*uz3
     call filz(uyz3f, ti3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0) !uy3*uz3

  endif

  if (nrank==0) print *, "filz ux= ", maxval(ta3), maxval(ux3f), maxval(ta3) - maxval(ux3f)

  ta3 = zero; tb3 = zero; tc3 = zero
  td3 = zero; te3 = zero; tf3 = zero
  tg3 = zero; th3 = zero; ti3 = zero

  ux2f = zero;uy2f = zero;uz2f = zero

  call transpose_z_to_y(ux3f, ux2f)
  call transpose_z_to_y(uy3f, uy2f)
  call transpose_z_to_y(uz3f, uz2f)
  call transpose_z_to_y(uxx3f, uxx2f)
  call transpose_z_to_y(uyy3f, uyy2f)
  call transpose_z_to_y(uzz3f, uzz2f)
  call transpose_z_to_y(uxy3f, uxy2f)
  call transpose_z_to_y(uxz3f, uxz2f)
  call transpose_z_to_y(uyz3f, uyz2f)

  ux1f = zero;uy1f = zero;uz1f = zero

  call transpose_y_to_x(ux2f, ux1f)
  call transpose_y_to_x(uy2f, uy1f)
  call transpose_y_to_x(uz2f, uz1f)
  call transpose_y_to_x(uxx2f, uxx1f)
  call transpose_y_to_x(uyy2f, uyy1f)
  call transpose_y_to_x(uzz2f, uzz1f)
  call transpose_y_to_x(uxy2f, uxy1f)
  call transpose_y_to_x(uxz2f, uxz1f)
  call transpose_y_to_x(uyz2f, uyz1f)

#ifdef IBM
  ux1f = ux1f * (one - ep1)
  uy1f = uy1f * (one - ep1)
  uz1f = uz1f * (one - ep1)
  uxx1f = uxx1f * (one - ep1)
  uyy1f = uyy1f * (one - ep1)
  uzz1f = uzz1f * (one - ep1)
  uxy1f = uxy1f * (one - ep1)
  uxz1f = uxz1f * (one - ep1)
  uyz1f = uyz1f * (one - ep1)
#endif

  !Lij tensor OK
  lxx1 = uxx1f - ux1f * ux1f
  lyy1 = uyy1f - uy1f * uy1f
  lzz1 = uzz1f - uz1f * uz1f
  lxy1 = uxy1f - ux1f * uy1f
  lxz1 = uxz1f - ux1f * uz1f
  lyz1 = uyz1f - uy1f * uz1f

  call derx (gxx1f, ux1f, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0)
  call derx (gyx1f, uy1f, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)
  call derx (gzx1f, uz1f, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)

  sxx1f = gxx1f

  !WORK Y-PENCILS
  call transpose_x_to_y(ux1f, ux2f)
  call transpose_x_to_y(uy1f, uy2f)
  call transpose_x_to_y(uz1f, uz2f)
  call transpose_x_to_y(gyx1f, ta2)

  call dery (gxy2f, ux2f, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)
  call dery (gyy2f, uy2f, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)
  call dery (gzy2f, uz2f, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)

  sxy2f = half * (gxy2f + ta2)
  syy2f = gyy2f

  !WORK Z-PENCILS
  call transpose_y_to_z(ux2f, ux3f)
  call transpose_y_to_z(uy2f, uy3f)
  call transpose_y_to_z(uz2f, uz3f)
  call transpose_y_to_z(gzy2f, ta3)

  call derz(gxz3f, ux3f, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
  call derz(gyz3f, uy3f, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
  call derz(gzz3f, uz3f, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0)

  szz3f = gzz3f
  syz3f = half * (gyz3f + ta3)

  !WORK Y-PENCILS
  call transpose_z_to_y(syz3f, syz2f)
  call transpose_z_to_y(szz3f, szz2f)
  call transpose_z_to_y(gxz3f, gxz2f)

  !WORK X-PENCILS
  call transpose_y_to_x(sxy2f, sxy1f)
  call transpose_y_to_x(syy2f, syy1f)
  call transpose_y_to_x(syz2f, syz1f)
  call transpose_y_to_x(szz2f, szz1f)
  call transpose_y_to_x(gxz2f, gxz1f)

  sxz1f = half * (gzx1f + gxz1f)

  !Bij tensor with u test filtered OK
  bbxx1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * sxx1f
  bbyy1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * syy1f
  bbzz1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * szz1f

  bbxy1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * sxy1f
  bbxz1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * sxz1f
  bbyz1 = -two * sqrt(two * (sxx1f * sxx1f + syy1f * syy1f + szz1f * szz1f + two * sxy1f * sxy1f + two * sxz1f * sxz1f + two * syz1f * syz1f)) * syz1f

  !Aij tensor with u
  axx1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * sxx1
  ayy1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * syy1
  azz1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * szz1

  axy1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * sxy1
  axz1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * sxz1
  ayz1 = -two * sqrt(two * (sxx1 * sxx1 + syy1 * syy1 + szz1 * szz1 + two * sxy1 * sxy1 + two * sxz1 * sxz1 + two * syz1 * syz1)) * syz1

#ifdef IBM
  bbxx1 = bbxx1 * (one - ep1)
  bbyy1 = bbyy1 * (one - ep1)
  bbzz1 = bbzz1 * (one - ep1)
  bbxy1 = bbxy1 * (one - ep1)
  bbxz1 = bbxz1 * (one - ep1)
  bbyz1 = bbyz1 * (one - ep1)
  axx1 = axx1 * (one - ep1)
  ayy1 = ayy1 * (one - ep1)
  azz1 = azz1 * (one - ep1)
  axy1 = axy1 * (one - ep1)
  axz1 = axz1 * (one - ep1)
  ayz1 = ayz1 * (one - ep1)
#endif

  call transpose_x_to_y(axx1, ta2)
  call transpose_x_to_y(ayy1, tb2)
  call transpose_x_to_y(azz1, tc2)
  call transpose_x_to_y(axy1, td2)
  call transpose_x_to_y(axz1, te2)
  call transpose_x_to_y(ayz1, tf2)

  do j = 1, ysize(2)
     ta2(:, j, :) = ta2(:, j, :) * (del(j))**two
     tb2(:, j, :) = tb2(:, j, :) * (del(j))**two
     tc2(:, j, :) = tc2(:, j, :) * (del(j))**two
     td2(:, j, :) = td2(:, j, :) * (del(j))**two
     te2(:, j, :) = te2(:, j, :) * (del(j))**two
     tf2(:, j, :) = tf2(:, j, :) * (del(j))**two
  enddo

  call transpose_y_to_x(ta2, axx1)
  call transpose_y_to_x(tb2, ayy1)
  call transpose_y_to_x(tc2, azz1)
  call transpose_y_to_x(td2, axy1)
  call transpose_y_to_x(te2, axz1)
  call transpose_y_to_x(tf2, ayz1)

  call transpose_x_to_y(bbxx1, ta2)
  call transpose_x_to_y(bbyy1, tb2)
  call transpose_x_to_y(bbzz1, tc2)
  call transpose_x_to_y(bbxy1, td2)
  call transpose_x_to_y(bbxz1, te2)
  call transpose_x_to_y(bbyz1, tf2)

  do j = 1, ysize(2)
     ta2(:, j, :) = ta2(:, j, :) * (two * del(j))**two
     tb2(:, j, :) = tb2(:, j, :) * (two * del(j))**two
     tc2(:, j, :) = tc2(:, j, :) * (two * del(j))**two
     td2(:, j, :) = td2(:, j, :) * (two * del(j))**two
     te2(:, j, :) = te2(:, j, :) * (two * del(j))**two
     tf2(:, j, :) = tf2(:, j, :) * (two * del(j))**two
  enddo

  call transpose_y_to_x(ta2, bbxx1)
  call transpose_y_to_x(tb2, bbyy1)
  call transpose_y_to_x(tc2, bbzz1)
  call transpose_y_to_x(td2, bbxy1)
  call transpose_y_to_x(te2, bbxz1)
  call transpose_y_to_x(tf2, bbyz1)

  !Need to filter Aij components

  if (nclx1.eq.0.and.nclxn.eq.0) then

     call filx(axx1f, axx1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)
     call filx(ayy1f, ayy1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)
     call filx(azz1f, azz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)

     call filx(axy1f, axy1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0)
     call filx(axz1f, axz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0)
     call filx(ayz1f, ayz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)

  else

     call filx(axx1f, axx1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)
     call filx(ayy1f, ayy1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)
     call filx(azz1f, azz1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)

     call filx(axy1f, axy1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0)
     call filx(axz1f, axz1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0)
     call filx(ayz1f, ayz1, di1, sx, vx, fiffxp, fifxp, ficxp, fibxp, fibbxp, filaxp, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 1)

  endif

  if (nrank==0) print *, "filx axx1= ", maxval(axx1), maxval(axx1f), maxval(axx1) - maxval(axx1f)

#ifdef IBM
  axx1f = axx1f * (one - ep1)
  ayy1f = ayy1f * (one - ep1)
  azz1f = azz1f * (one - ep1)
  axy1f = axy1f * (one - ep1)
  axz1f = axz1f * (one - ep1)
  ayz1f = ayz1f * (one - ep1)
#endif

  call transpose_x_to_y(axx1f, ta2)
  call transpose_x_to_y(ayy1f, tb2)
  call transpose_x_to_y(azz1f, tc2)
  call transpose_x_to_y(axy1f, td2)
  call transpose_x_to_y(axz1f, te2)
  call transpose_x_to_y(ayz1f, tf2)

  if (ncly1.eq.0.and.nclyn.eq.0) then

     call fily(axx2f, ta2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)
     call fily(ayy2f, tb2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)
     call fily(azz2f, tc2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)

     call fily(axy2f, td2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0)
     call fily(axz2f, te2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)
     call fily(ayz2f, tf2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0)

  else

     call fily(axx2f, ta2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)
     call fily(ayy2f, tb2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)
     call fily(azz2f, tc2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)

     call fily(axy2f, td2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0)
     call fily(axz2f, te2, di2, sy, vy, fiffyp, fifyp, ficyp, fibyp, fibbyp, filayp, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 1)
     call fily(ayz2f, tf2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0)

  endif

  if (nrank==0) print *, "fily axx2= ", maxval(ta2), maxval(axx2f), maxval(ta2) - maxval(axx2f)

#ifdef IBM
  axx2f = axx2f * (one - ep2)
  ayy2f = ayy2f * (one - ep2)
  azz2f = azz2f * (one - ep2)
  axy2f = axy2f * (one - ep2)
  axz2f = axz2f * (one - ep2)
  ayz2f = ayz2f * (one - ep2)
#endif

  ta2 = zero; tb2 = zero; tc2 = zero
  td2 = zero; te2 = zero; tf2 = zero

  call transpose_y_to_z(axx2f, ta3)
  call transpose_y_to_z(ayy2f, tb3)
  call transpose_y_to_z(azz2f, tc3)
  call transpose_y_to_z(axy2f, td3)
  call transpose_y_to_z(axz2f, te3)
  call transpose_y_to_z(ayz2f, tf3)

  if (nclz1.eq.0.and.nclzn.eq.0) then

     call filz(axx3f, ta3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)
     call filz(ayy3f, tb3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)
     call filz(azz3f, tc3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)

     call filz(axy3f, td3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)
     call filz(axz3f, te3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0)
     call filz(ayz3f, tf3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0)

  else

     call filz(axx3f, ta3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)
     call filz(ayy3f, tb3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)
     call filz(azz3f, tc3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)

     call filz(axy3f, td3, di3, sz, vz, fiffzp, fifzp, ficzp, fibzp, fibbzp, filazp, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 1)
     call filz(axz3f, te3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0)
     call filz(ayz3f, tf3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0)

  endif

  if (nrank==0) print *, "filz axx3= ", maxval(ta3), maxval(axx3f), maxval(ta3) - maxval(axx3f)

  ta3 = zero; tb3 = zero; tc3 = zero
  td3 = zero; te3 = zero; tf3 = zero

  call transpose_z_to_y(axx3f, axx2f)
  call transpose_z_to_y(ayy3f, ayy2f)
  call transpose_z_to_y(azz3f, azz2f)
  call transpose_z_to_y(axy3f, axy2f)
  call transpose_z_to_y(axz3f, axz2f)
  call transpose_z_to_y(ayz3f, ayz2f)

  call transpose_y_to_x(axx2f, axx1f)
  call transpose_y_to_x(ayy2f, ayy1f)
  call transpose_y_to_x(azz2f, azz1f)
  call transpose_y_to_x(axy2f, axy1f)
  call transpose_y_to_x(axz2f, axz1f)
  call transpose_y_to_x(ayz2f, ayz1f)

#ifdef IBM
  axx1f = axx1f * (one - ep1)
  ayy1f = ayy1f * (one - ep1)
  azz1f = azz1f * (one - ep1)
  axy1f = axy1f * (one - ep1)
  axz1f = axz1f * (one - ep1)
  ayz1f = ayz1f * (one - ep1)
#endif

  !Mij tensor OK !!!M_{ij} = B_{ij} - Aij
  mxx1 = bbxx1 - axx1f
  myy1 = bbyy1 - ayy1f
  mzz1 = bbzz1 - azz1f
  mxy1 = bbxy1 - axy1f
  mxz1 = bbxz1 - axz1f
  myz1 = bbyz1 - ayz1f

  !Lij deviator
  lxx1 = lxx1 - (lxx1 + lyy1 + lzz1) / three
  lyy1 = lyy1 - (lxx1 + lyy1 + lzz1) / three
  lzz1 = lzz1 - (lxx1 + lyy1 + lzz1) / three

#ifdef IBM
  do ijk = 1, nvect1
     if (ep1(ijk, 1, 1) .eq. one) then
        ta1(ijk, 1, 1) = zero
        tb1(ijk, 1, 1) = one
     endif
  enddo
#endif

  !MODEL OF LILLY (1992)
  smagC1 = (lxx1 * mxx1 + lyy1 * myy1 + lzz1 * mzz1 + two * (lxy1 * mxy1 + lxz1 * mxz1 + lyz1 * myz1)) / &
       (mxx1 * mxx1 + myy1 * myy1 + mzz1 * mzz1 + two * (mxy1 * mxy1 + mxz1 * mxz1 + myz1 * myz1)) !l/M

  do ijk = 1, nvect1 !ERIC LIMITEUR SI BESOIN
     if (smagC1(ijk, 1, 1).gt. maxdsmagcst) smagC1(ijk, 1, 1) = zero
     if (smagC1(ijk, 1, 1).lt. 0.0) smagC1(ijk, 1, 1) = zero
  enddo

  !FILTERING THE NON-CONSTANT CONSTANT
  call filx(smagC1f, smagC1, di1, sx, vx, fiffx, fifx, ficx, fibx, fibbx, filax, fiz1x, fiz2x, xsize(1), xsize(2), xsize(3), 0)
  if (nrank==0) print *, "filx smagC1= ", maxval(smagC1), maxval(smagC1f), maxval(smagC1) - maxval(smagC1f)

  call transpose_x_to_y(smagC1f, ta2)
  call fily(smagC2f, ta2, di2, sy, vy, fiffy, fify, ficy, fiby, fibby, filay, fiz1y, fiz2y, ysize(1), ysize(2), ysize(3), 0)
  if (nrank==0) print *, "fily smagC1= ", maxval(ta2), maxval(smagC2f), maxval(ta2) - maxval(smagC2f)

  call transpose_y_to_z(smagC2f, ta3)
  call filz(smagC3f, ta3, di3, sz, vz, fiffz, fifz, ficz, fibz, fibbz, filaz, fiz1z, fiz2z, zsize(1), zsize(2), zsize(3), 0)
  if (nrank==0) print *, "filz smagC1= ", maxval(ta3), maxval(smagC3f), maxval(ta3) - maxval(smagC3f)

  dsmagcst3 = zero
  call mean_plane_z(smagC3f, zsize(1), zsize(2), zsize(3), dsmagcst3(:, :, 1))

  do k = 2, zsize(3)
     dsmagcst3(:, :, k) = dsmagcst3(:, :, 1)
  enddo

  call transpose_z_to_y(dsmagcst3, dsmagcst2)
  !call transpose_z_to_y(smagC3f,dsmagcst2)

  call transpose_y_to_x(dsmagcst2, dsmagcst1)

  do ijk = 1, nvect1 !ERIC LIMITEUR SI BESOIN
     if (dsmagcst1(ijk, 1, 1).gt. maxdsmagcst) dsmagcst1(ijk, 1, 1) = zero
     if (dsmagcst1(ijk, 1, 1).lt. 0.0) dsmagcst1(ijk, 1, 1) = zero
  enddo

  nut1 = zero; nut2 = zero
  call transpose_x_to_y(srt_smag, srt_smag2)
  call transpose_x_to_y(dsmagcst1, dsmagcst2)
  do k = 1, ysize(3)
     do j = 1, ysize(2)
        do i = 1, ysize(1)
           nut2(i, j, k) = dsmagcst2(i, j, k) * ((FSGS * del(j))**two) * sqrt(two * srt_smag2(i, j, k))
        enddo
     enddo
  enddo
  call transpose_y_to_x(nut2, nut1)

  if (nrank==0) print *, "dsmagc init   min max= ", minval(smagC1), maxval(smagC1)
  if (nrank==0) print *, "dsmagc final  min max= ", minval(dsmagcst1), maxval(dsmagcst1)
  if (nrank==0) print *, "dsmag nut1    min max= ", minval(nut1), maxval(nut1)

  if (mod(itime, imodulo).eq.0) then

     write(filename, "('./data/dsmagcst_initial',I4.4)") itime / imodulo
     call decomp_2d_write_one(1, smagC1, filename, 2)

     write(filename, "('./data/dsmagcst_final',I4.4)") itime / imodulo
     call decomp_2d_write_one(1, dsmagcst1, filename, 2)

     write(filename, "('./data/nut_dynsmag',I4.4)") itime / imodulo
     call decomp_2d_write_one(1, nut1, filename, 2)

  endif

end subroutine dynsmag
!************************************************************
subroutine lesdiff(ux1, uy1, uz1, ep1, sxx1, syy1, szz1, sxy1, sxz1, syz1, nut1, &
     sgsx1, sgsy1, sgsz1, ta1, td1, te1, tf1, di1, ta2, td2, te2, tf2, tj2, di2, &
     ta3, td3, te3, tf3, di3)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, nut1, ep1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2, uy2, uz2, nut2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ux3, uy3, uz3, nut3
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ta1, td1, te1, tf1, di1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ta2, td2, te2, tf2, di2, tj2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ta3, td3, te3, tf3, di3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sgsx1, sgsy1, sgsz1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: sgsx2, sgsy2, sgsz2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: sgsx3, sgsy3, sgsz3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: sxx1, syy1, szz1, sxy1, sxz1, syz1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: syy2, szz2, sxy2, sxz2, syz2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: szz3, sxz3, syz3

  integer :: i, j, k, ijk, nvect1

  ta1 = zero; ta2 = zero; ta3 = zero
  nvect1 = xsize(1) * xsize(2) * xsize(3)

  !WORK X-PENCILS

  call derx (ta1, nut1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0)

  call derxx (td1, ux1, di1, sx, sfx, ssx, swx, xsize(1), xsize(2), xsize(3), 0)
  call derxx (te1, uy1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)
  call derxx (tf1, uz1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)

  !nu_t constant
  !sgsx1 = td1*nut1 + gxx1*ta1
  !sgsy1 = te1*nut1 + gyx1*ta1
  !sgsz1 = tf1*nut1 + gzx1*ta1

  sgsx1 = td1 * nut1 + two * sxx1 * ta1
  sgsy1 = te1 * nut1 + two * sxy1 * ta1
  sgsz1 = tf1 * nut1 + two * sxz1 * ta1

  !WORK Y-PENCILS
  call transpose_x_to_y(sgsx1, sgsx2)
  call transpose_x_to_y(sgsy1, sgsy2)
  call transpose_x_to_y(sgsz1, sgsz2)
  call transpose_x_to_y(sxz1, sxz2)
  call transpose_x_to_y(sxy1, sxy2)
  call transpose_x_to_y(syy1, syy2)
  call transpose_x_to_y(syz1, syz2)
  call transpose_x_to_y(szz1, szz2)
  call transpose_x_to_y(nut1, nut2)
  call transpose_x_to_y(ux1, ux2)
  call transpose_x_to_y(uy1, uy2)
  call transpose_x_to_y(uz1, uz2)

  call dery (ta2, nut2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)

  !-->for ux
  td2 = zero
  if (istret.ne.0) then
     call deryy (td2, ux2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
     call dery (te2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)
     do k = 1, ysize(3)
        do j = 1, ysize(2)
           do i = 1, ysize(1)
              td2(i, j, k) = td2(i, j, k) * pp2y(j) - pp4y(j) * te2(i, j, k)
           enddo
        enddo
     enddo
  else
     call deryy (td2, ux2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
  endif

  !-->for uy
  te2 = zero
  if (istret.ne.0) then
     call deryy (te2, uy2, di2, sy, sfy, ssy, swy, ysize(1), ysize(2), ysize(3), 0)
     call dery (tf2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)
     do k = 1, ysize(3)
        do j = 1, ysize(2)
           do i = 1, ysize(1)
              te2(i, j, k) = te2(i, j, k) * pp2y(j) - pp4y(j) * tf2(i, j, k)
           enddo
        enddo
     enddo
  else
     call deryy (te2, uy2, di2, sy, sfy, ssy, swy, ysize(1), ysize(2), ysize(3), 0)
  endif

  !-->for uz
  tf2 = zero
  if (istret.ne.0) then
     call deryy (tf2, uz2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
     call dery (tj2, uz2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)
     do k = 1, ysize(3)
        do j = 1, ysize(2)
           do i = 1, ysize(1)
              tf2(i, j, k) = tf2(i, j, k) * pp2y(j) - pp4y(j) * tj2(i, j, k)
           enddo
        enddo
     enddo
  else
     call deryy (tf2, uz2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
  endif

  !nu_t constant
  !sgsx2=sgsx2+nut2*td2+gxy2*ta2
  !sgsy2=sgsy2+nut2*te2+gyy2*ta2
  !sgsz2=sgsz2+nut2*tf2+gzy2*ta2

  sgsx2 = sgsx2 + nut2 * td2 + two * sxy2 * ta2
  sgsy2 = sgsy2 + nut2 * te2 + two * syy2 * ta2
  sgsz2 = sgsz2 + nut2 * tf2 + two * syz2 * ta2

  !WORK Z-PENCILS
  call transpose_y_to_z(sgsx2, sgsx3)
  call transpose_y_to_z(sgsy2, sgsy3)
  call transpose_y_to_z(sgsz2, sgsz3)
  call transpose_y_to_z(sxz2, sxz3)
  call transpose_y_to_z(syz2, syz3)
  call transpose_y_to_z(szz2, szz3)
  call transpose_y_to_z(nut2, nut3)
  call transpose_y_to_z(ux2, ux3)
  call transpose_y_to_z(uy2, uy3)
  call transpose_y_to_z(uz2, uz3)

  call derz (ta3, nut3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0)

  call derzz (td3, ux3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)
  call derzz (te3, uy3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)
  call derzz (tf3, uz3, di3, sz, sfz, ssz, swz, zsize(1), zsize(2), zsize(3), 0)

  !nu_t constant
  !sgsx3=sgsx3+nut3*td3+gxz3*ta3
  !sgsy3=sgsy3+nut3*te3+gyz3*ta3
  !sgsz3=sgsz3+nut3*tf3+gzz3*ta3

  sgsx3 = sgsx3 + nut3 * td3 + two * sxz3 * ta3
  sgsy3 = sgsy3 + nut3 * te3 + two * syz3 * ta3
  sgsz3 = sgsz3 + nut3 * tf3 + two * szz3 * ta3

  call transpose_z_to_y(sgsx3, sgsx2)
  call transpose_z_to_y(sgsy3, sgsy2)
  call transpose_z_to_y(sgsz3, sgsz2)

  call transpose_y_to_x(sgsx2, sgsx1)
  call transpose_y_to_x(sgsy2, sgsy1)
  call transpose_y_to_x(sgsz2, sgsz1)

#ifdef IBM
  do ijk = 1, nvect1
     if(ep1(ijk, 1, 1).eq.1) then
        sgsx1(ijk, 1, 1) = zero
        sgsy1(ijk, 1, 1) = zero
        sgsz1(ijk, 1, 1) = zero
     endif
  enddo
#endif

  !call transpose_x_to_y(sgsx1,sgsx2)
  !call transpose_x_to_y(sgsy1,sgsy2)
  !call transpose_x_to_y(sgsz1,sgsz2)

  !sgsx2(:,1,:) = zero; sgsx2(:,ysize(2),:) = zero
  !sgsy2(:,1,:) = zero; sgsy2(:,ysize(2),:) = zero
  !sgsz2(:,1,:) = zero; sgsz2(:,ysize(2),:) = zero

  !call transpose_y_to_x(sgsx2,sgsx1)
  !call transpose_y_to_x(sgsy2,sgsy1)
  !call transpose_y_to_x(sgsz2,sgsz1)

end subroutine lesdiff

!************************************************************
subroutine lesdiff_scalar(phi1, di1, di2, di3, kappat1, sgsphi1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  real(mytype), dimension(xsize(1), xsize(2), xsize(3), nphi) :: phi1, sgsphi1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: phi2, sgsphi2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: phi3, sgsphi3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: di1, tb1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: di2, tb2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: di3, tb3

  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: kappat1, dkappat1
  real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: kappat2, dkappat2
  real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: kappat3, dkappat3

  integer :: is

  sgsphi1 = zero; sgsphi2 = zero; sgsphi3 = zero

  call derxS (dkappat1, kappat1, di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1)
  call transpose_x_to_y(kappat1, kappat2)
  call deryS (dkappat2, kappat2, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1)
  call transpose_y_to_z(kappat2, kappat3)
  call derzS (dkappat3, kappat3, di3, sz, ffzpS, fszpS, fwzpS, zsize(1), zsize(2), zsize(3), 1)

  do is = 1, nphi
     call derxS (tb1, phi1(:, :, :, is), di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1)
     sgsphi1(:, :, :, is) = tb1 * dkappat1 !d(phi)/dx * d(kappa_t)/dx

     call transpose_x_to_y(phi1(:, :, :, is), phi2)
     call transpose_x_to_y(sgsphi1(:, :, :, is),sgsphi2)

     call deryS (tb2, phi2, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1)
     sgsphi2 = sgsphi2 + tb2 * dkappat2 !d(phi)/dy * d(kappa_t)/dy

     call transpose_y_to_z(phi2, phi3)
     call transpose_y_to_z(sgsphi2, sgsphi3)

     call derzS (tb3, phi3, di3, sz, ffzpS, fszpS, fwzpS, zsize(1), zsize(2), zsize(3), 1)
     sgsphi3 = sgsphi3 + tb3 * dkappat3 !d(phi)/dz * d(kappa_t)/dz

     call transpose_z_to_y(sgsphi3, sgsphi2)
     call transpose_y_to_x(sgsphi2, sgsphi1(:, :, :, is))
  end do

end subroutine lesdiff_scalar

#endif
