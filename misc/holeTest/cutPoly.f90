!
! perform a 2-d hole cut for an arbitrary polygon
! and find weights and indices of all points inside
! (weights are correct only if polygon is triangle or
! quad)
!
subroutine cutPoly(xlo,dx,xp,jkloc,fjk,jmax,kmax,nvert,nsect)
implicit none
integer, intent(in):: nvert,jmax,kmax !number of vertices of polygon and dims of cart grid
integer, intent(out) :: nsect         !number of points that are inside
real*8, intent(in):: xlo(2),dx        !corner point and dx of cart grid 
real*8, intent(in) :: xp(2,nvert)     !coords of the polygon
real*8, intent(out) :: fjk(2,jmax*kmax) ! weights for all points inside
integer, intent(out) :: jkloc(2,jmax*kmax) ! indices of all the point that are inside
integer, allocatable :: ipointer(:)
real*8, allocatable :: dd(:,:)
!
! local variables
!
integer :: i,j,k,id,ip,js,je,k1,k2,ipp,itr,jj,ii,js1,je1
real*8 :: x1(2),x2(2),tmp(2),xd(nvert+1),f1(2),f2(2),f3(2),xx(2)
real*8 :: a1(2),a2(2)
real*8 :: dyy,dxx,ds1,ds2,slope,xj,eps2,det
!
! begin
!
allocate(ipointer(jmax))
allocate(dd(2,(nvert+1)*jmax))
!
id=0
ipointer=0
eps2=1e-10
nsect=0
js1=jmax
je1=1
!
do i=1,nvert
   ip=mod(i,nvert)+1
   x1=xp(:,i)
   x2=xp(:,ip)

   if (x1(1) > x2(1)) then
      tmp=x1
      x1=x2
      x2=tmp
   endif
!
   dyy=x2(2)-x1(2)
   dxx=x2(1)-x1(1)
!   
   if (dxx > 0) then

      ds1=(x1(1)-xlo(1))/dx
      ds2=(x2(1)-xlo(1))/dx
      
      js=int(ds1)+1+ceiling(ds1-int(ds1))
      je=int(ds2)+1+ceiling(ds2-int(ds2))
      js=max(1,js)
      js=min(js,jmax)
      je=max(1,je)
      je=min(je,jmax)

      slope=dyy/dxx
   
      do j=js,je
         xj=xlo(1)+dx*(j-1)
         if  ((xj-x1(1))*(xj-x2(1)) <= 0.) then
            id=id+1
            dd(1,id)=x1(2)+slope*(xj-x1(1))
            dd(2,id)=ipointer(j)
            ipointer(j)=id
            if (js1.gt.j) js1=j
            if (je1.lt.j) je1=j
         endif
      enddo
   endif
!
enddo
!
if (js1.gt.je1) return
!
! do some preprocessing on the polygon
!
if (nvert==3) then
   a1=xp(:,2)-xp(:,1)
   a2=xp(:,3)-xp(:,1)
   det=1./(a1(1)*a2(2)-a1(2)*a2(1))
   f1(1)=a2(2)*det
   f1(2)=-a1(2)*det
   f2(1)=-a2(1)*det
   f2(2)=a1(1)*det
else
   f1=xp(:,2)-xp(:,1)
   f2=xp(:,4)-xp(:,1)
   f3=xp(:,1)-xp(:,2)+xp(:,3)-xp(:,4)
endif
!
do j=js1,je1
   ipp=ipointer(j)
   itr=0
   do while(ipp > 0) 
      itr=itr+1
      xd(itr)=dd(1,ipp)
      ipp=dd(2,ipp)
   enddo
   ii=2
   cullLoop: do while(ii.le.itr)
      if (abs(xd(ii)-xd(ii-1)).lt.eps2) then
         do jj=ii+1,itr
            xd(jj-1)=xd(jj)
         enddo
         itr=itr-1
         ii=ii-1
      endif
      ii=ii+1
   enddo cullLoop

   sortLoop: do ii=1,itr
      do jj=ii+1,itr
         if (xd(jj).lt.xd(ii)) then
            tmp(1)=xd(ii)
            xd(ii)=xd(jj)
            xd(jj)=tmp(1)
         endif
      enddo
   enddo sortLoop
   ii=1
   do while(ii.lt.itr)
      k1=max(ceiling((xd(ii)-xlo(2))/dx)+1,1)
      k2=min(floor((xd(ii+1)-xlo(2))/dx)+1,kmax)
      do k=k1,k2
         nsect=nsect+1
         jkloc(1,nsect)=j
         jkloc(2,nsect)=k
         xx(1)=(j-1)*dx+xlo(1)-xp(1,1)
         xx(2)=(k-1)*dx+xlo(2)-xp(2,1)
         call solveLoc(fjk(1,nsect),fjk(2,nsect),f1,f2,f3,&
              xx,nvert)
      enddo
      ii=ii+2
   enddo
enddo
!
deallocate(ipointer,dd)
!
return
end subroutine cutPoly
!
! solve parametric coordinates for a 
! contained point
!
!   
subroutine solveLoc(u,v,f1,f2,f3,xx,nvert) 
implicit none
integer :: nvert
real*8, intent(out) :: u,v
real*8, intent(in) :: f1(2),f2(2),f3(2),xx(2)
!
! local variables
!
real*8 f(2),fd(2,2),fdi(2,2),det,du,dv,tol
integer :: iter,itmax
!
if (nvert==3) then
   u=f1(1)*xx(1)+f2(1)*xx(2)
   v=f1(2)*xx(1)+f2(2)*xx(2)
else
   f=1e10
   itmax=50
   tol=1E-12
   do iter=1,itmax

      f=f1*u+f2*v+f3*u*v-xx
      
      if (sqrt(f(1)**2+f(2)**2) <=tol) exit
      
      fd(:,1)=f1+f3*v
      fd(:,2)=f2+f3*u
      det=1./(fd(1,1)*fd(2,2)-fd(1,2)*fd(2,1))
      fdi(1,1)=fd(2,2)*det
      fdi(1,2)=-fd(1,2)*det
      fdi(2,1)=-fd(2,1)*det
      fdi(2,2)=fd(1,1)*det
      
      du=fdi(1,1)*f(1)+fdi(1,2)*f(2)
      dv=fdi(2,1)*f(1)+fdi(2,2)*f(2)

      u=u-du
      v=v-dv
   enddo
endif
!
return
end subroutine solveLoc

