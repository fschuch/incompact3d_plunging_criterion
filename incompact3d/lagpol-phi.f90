#ifdef IBM
subroutine lagpolphiy(u)
  !
  USE flow_type, only : ramp
  USE param
  USE complex_geometry
  USE decomp_2d
  USE variables
  !
  implicit none
  !
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
  integer                                            :: i,j,k
  real(mytype)                                       :: x,y,z
  integer                                            :: jy              != position du point "zappé"
  integer                                            :: jpif,jpol,nypif
  integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
  real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
  real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en
  integer                                            :: ia,na           !|double précision
  !
  do k=1,ysize(3)
    do i=1,ysize(1)
      if(nobjy(i,k).ne.0)then
        ia=0
        do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
          !1ère frontière
          ia=ia+1
          xa(ia)=zero
          ya(ia)=zero
          jpoli=1

          !2ème frontière
          nypif=npif
          !ia=ia+1
          !xa(ia)=yf(j,i,k)
          !ya(ia)=zero
          if(yf(j,i,k).lt.yly)then!objet immergé
            jy=1!jy=(yf(j,i,k)+dy)/dy+1
            do while(yp(jy).lt.yf(j,i,k))  !there was a bug here yi<-->yf
              jy=jy+1
            enddo
            jpolf=jy-1
            if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
            do jpif=1,nypif
              ia=ia+1
              if(izap.eq.1)then!zapping
                xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
                ya(ia)=u(i,jy+jpif,k)
                !
                ia=ia+1
                xa(ia)=two*ramp(i+ystart(1)-1)-yp(jy+jpif)
                ya(ia)=u(i,jy+jpif,k)
              else             !no zapping
                xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
                ya(ia)=u(i,jy+jpif-1,k)
                !
                ia=ia+1
                xa(ia)=two*ramp(i+ystart(1)-1)-yp(jy+jpif-1)
                ya(ia)=u(i,jy+jpif-1,k)
              endif
            enddo
          else                   !objet semi-immergé
            jpolf=ny
          endif
          !calcul du polynôme
          na=ia
          do jpol=jpoli,jpolf
            xpol=yp(jpol)!dy*(jpol-1)
            call polint(xa,ya,na,xpol,ypol,dypol)
            u(i,jpol,k)=ypol
          enddo
        enddo
      endif
    enddo
  enddo
  !
  return
end subroutine lagpolphiy
#endif
