  case ('#')
     !It is just a comment
  case ('angle')
     if (arg.eq.1) read(buffer, *, iostat=ios) angle
  case ('beta')
     if (arg.eq.1) read(buffer, *, iostat=ios) beta
  case ('cont_phi')
     if (arg.ge.1) read(buffer, *, iostat=ios) cont_phi
  case ('cp')
     if (arg.eq.1) read(buffer, *, iostat=ios) cp
  case ('datapath')
     if (arg.eq.1) read(buffer, *, iostat=ios) datapath
  case ('dt')
     if (arg.eq.1) read(buffer, *, iostat=ios) dt
  case ('fpi2')
     if (arg.eq.1) read(buffer, *, iostat=ios) fpi2
  case ('ifirst')
     if (arg.eq.1) read(buffer, *, iostat=ios) ifirst
  case ('iin')
     if (arg.eq.1) read(buffer, *, iostat=ios) iin
  case ('ilast')
     if (arg.ge.1) read(buffer, *, iostat=ios) ilast
  case ('ilit')
     if (arg.eq.1) read(buffer, *, iostat=ios) ilit
  case ('imodulo')
     if (arg.eq.1) read(buffer, *, iostat=ios) imodulo
  case ('imodulo2')
     if (arg.eq.1) read(buffer, *, iostat=ios) imodulo2
  case ('iprocessing')
     if (arg.ge.1) read(buffer, *, iostat=ios) iprocessing
  case ('irotation')
     if (arg.eq.1) read(buffer, *, iostat=ios) irotation
  case ('isave')
     if (arg.ge.1) read(buffer, *, iostat=ios) isave
  case ('iscalar')
     if (arg.eq.0) read(buffer, *, iostat=ios) iscalar
  case ('istret')
     if (arg.eq.1) read(buffer, *, iostat=ios) istret
  case ('istat')
     if (arg.ge.1) read(buffer, *, iostat=ios) istat
  case ('itest')
     if (arg.ge.1) read(buffer, *, iostat=ios) itest
  case ('itrip')
     if (arg.eq.1) read(buffer, *, iostat=ios) itrip
  case ('jLES')
     if (arg.eq.1) read(buffer, *, iostat=ios) jLES
  case ('nclx1')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclx1
  case ('nclxn')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclxn
  case ('nclxS1')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclxS1
  case ('nclxSn')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclxSn
  case ('ncly1')
     if (arg.eq.1) read(buffer, *, iostat=ios) ncly1
  case ('nclyn')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclyn
  case ('nclyS1')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclyS1
  case ('nclySn')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclySn
  case ('nclz1')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclz1
  case ('nclzn')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclzn
  case ('nclzS1')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclzS1
  case ('nclzSn')
     if (arg.eq.1) read(buffer, *, iostat=ios) nclzSn
  case ('noise')
     if (arg.ge.1) read(buffer, *, iostat=ios) noise
  case ('noise1')
     if (arg.ge.1) read(buffer, *, iostat=ios) noise1
  case ('nphi')
     if (arg.eq.0) read(buffer, *, iostat=ios) nphi
  case ('npif')
     if (arg.eq.1) read(buffer, *, iostat=ios) npif
  case ('nsc')
     if (arg.eq.1) read(buffer, *, iostat=ios) nsc
  case ('nscheme')
     if (arg.eq.1) read(buffer, *, iostat=ios) nscheme
  case ('nx')
     if (arg.eq.0) read(buffer, *, iostat=ios) nx
  case ('ny')
     if (arg.eq.0) read(buffer, *, iostat=ios) ny
  case ('nz')
     if (arg.eq.0) read(buffer, *, iostat=ios) nz
  case ('p_col')
     if (arg.eq.0.and.p_col.eq.0) read(buffer, *, iostat=ios) p_col
  case ('p_row')
     if (arg.eq.0.and.p_col.eq.0) read(buffer, *, iostat=ios) p_row
  case ('re')
     if (arg.eq.1) read(buffer, *, iostat=ios) re
  case ('ri')
     if (arg.eq.1) read(buffer, *, iostat=ios) ri
  case ('ro')
     if (arg.eq.1) read(buffer, *, iostat=ios) ro
  case ('save_dmap')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dmap
  case ('save_dphidx')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dphidx
  case ('save_dphidy')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dphidy
  case ('save_dphidz')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dphidz
  case ('save_dudx')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dudx
  case ('save_dudy')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dudy
  case ('save_dudz')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dudz
  case ('save_dvdx')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dvdx
  case ('save_dvdy')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dvdy
  case ('save_dvdz')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dvdz
  case ('save_dwdx')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dwdx
  case ('save_dwdy')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dwdy
  case ('save_dwdz')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_dwdz
  case ('save_pc')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_pc
  case ('save_phi')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_phi
  case ('save_phim')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_phim
  case ('save_pre')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_pre
  case ('save_prem')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_prem
  case ('save_qc')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_qc
  case ('save_utmap')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_utmap
  case ('save_ux')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_ux
  case ('save_uxm')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_uxm
  case ('save_uy')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_uy
  case ('save_uym')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_uym
  case ('save_uz')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_uz
  case ('save_uzm')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_uzm
  case ('save_V')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_V
  case ('save_w')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_w
  case ('save_w1')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_w1
  case ('save_w2')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_w2
  case ('save_w3')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_w3
  case ('u1')
     if (arg.eq.1) read(buffer, *, iostat=ios) u1
  case ('u2')
     if (arg.eq.1) read(buffer, *, iostat=ios) u2
  case ('uset')
     if (arg.eq.1) read(buffer, *, iostat=ios) uset
  case ('wrotation')
     if (arg.eq.1) read(buffer, *, iostat=ios) wrotation
  case ('xlx')
     if (arg.eq.1) read(buffer, *, iostat=ios) xlx
  case ('yly')
     if (arg.eq.1) read(buffer, *, iostat=ios) yly
  case ('zlz')
     if (arg.eq.1) read(buffer, *, iostat=ios) zlz
#ifdef IBM
  case ('ilag')
     if (arg.eq.1) read(buffer, *, iostat=ios) ilag
  case ('ivirt')
     if (arg.eq.1) read(buffer, *, iostat=ios) ivirt
  case ('izap')
     if (arg.eq.1) read(buffer, *, iostat=ios) izap
  case ('nobjmax')
     if (arg.eq.1) read(buffer, *, iostat=ios) nobjmax
  case ('nraf')
     if (arg.eq.1) read(buffer, *, iostat=ios) nraf
  case ('save_ibm')
     if (arg.ge.1) read(buffer, *, iostat=ios) save_ibm
#endif
#ifdef FORCES
  case ('xld')
     if (arg.eq.1) read(buffer, *, iostat=ios) xld
  case ('xrd')
     if (arg.eq.1) read(buffer, *, iostat=ios) xrd
  case ('yld')
     if (arg.eq.1) read(buffer, *, iostat=ios) yld
  case ('yud')
     if (arg.eq.1) read(buffer, *, iostat=ios) yud
#endif
#ifdef POST
  case ('icrfile')
    if (arg.eq.1) read(buffer, *, iostat=ios) icrfile
  case ('file1')
    if (arg.eq.1) read(buffer, *, iostat=ios) file1
  case ('filen')
    if (arg.eq.1) read(buffer, *, iostat=ios) filen
  case ('ssfile')
    if (arg.eq.1) read(buffer, *, iostat=ios) ssfile
  case ('comp_visu')
    if (arg.eq.1) read(buffer, *, iostat=ios) comp_visu
  case ('comp_post')
    if (arg.eq.1) read(buffer, *, iostat=ios) comp_post
  case ('probes')
    if (arg.eq.0) flag_probes = .True.
#endif
  case ('stop')
     if (arg.ge.1) then
       ilast = itime
       if (itime.ne.0) isave = itime
     endif
  case default
     if (nrank.eq.0) then
        print *,'Invalid label      : ', trim(label)
        print *,'    file           : ', trim(filepath)
        print *,'    line           : ', line
     endif
