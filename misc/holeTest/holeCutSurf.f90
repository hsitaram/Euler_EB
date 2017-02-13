!====================================================================
!
!   Routine to cut holes of a uniform isotropic cartesian mesh
!   using a unstructured triangle surface as the cutter
!   
!   Author Jay Sitaraman
!   Created :: 08/21/07 
!   Modified for pure iblanking on:: 12/03/09
!
!====================================================================
subroutine holeCutSurf(jmax,kmax,lmax,ibl,xl,ds,xsurf,ncf3,nfringe,nnodes,ntri,holeFlag,holeCount)
!
implicit none
!
! subroutine arguments
!
integer,intent(in) :: jmax,kmax,lmax,nfringe,nnodes,ntri !dimensions of cart grid, number of fringes
                                                         !number of nodes and triangles in the cutter mesh
integer, intent(inout) :: holeCount                      ! number of hole points
integer, intent(inout) :: ibl((jmax+2*nfringe)*(kmax+2*nfringe) & !iblank array that masks all the point inside
                              *(lmax+2*nfringe))                  ! with provided holeFlag
real*8, intent(in)    :: xl(3)                                    ! lower-left-bottom corner of cartesian grid
real*8, intent(in)    :: ds                                       ! size of cartesian grid cell 
real*8, intent(in)    :: xsurf(3,nnodes)                          ! surface grid nodes 
integer    :: ncf3(3,ntri)                                        ! surface grid connectivity 
integer    :: holeFlag                                            ! hole flag 
!
! local variables
!
real*8, allocatable :: iblank(:,:,:)
real*8, allocatable :: xdd(:)
integer, allocatable :: ipointer(:,:),idd(:)
integer, allocatable :: iltri(:)
real*8, allocatable :: cof(:,:),xd(:)
real*8, allocatable :: xrecvtmp(:,:)
integer, allocatable :: irecvtmp(:)
!
integer :: ibcount,buff,maxintersect,ipp,ixd
integer :: js,je,ks,ke,ls,le
integer :: jp1,jm1,kp1,km1,lp1,lm1
integer :: js1,je1,ks1,ke1,ls1,le1
integer :: ka,kb,la,lb
integer :: itri,maxtri
integer :: i,j,k,l,m,n,i0,j1,j2,i1,i2,ii,jj,itr,d0,d1,d2,inside,iflag
integer :: jt,kt,lt
real*8   :: l1,l2,l3
real*8   :: xmax(3),xmin(3),areamax
real*8   :: lenmax,xx,yy,zz,x0,x1,x2,y1,z1,tmp,fac,eps1,eps2,den,xlen,xcm
integer :: kk,ll,mm,ib1
integer :: jdim,jkdim,nf1,nf2
!
! begin
!
!write(6,*) 'tri1=',ncf3(:,1)
ibcount=0
holeCount=0
!
! initialize local iblanking
!
allocate(iblank(jmax,kmax,lmax))
iblank=1
!
!
do m=1,3
   xmin(m)=xsurf(m,ncf3(1,1))
   xmax(m)=xmin(m)
enddo
!
!     find bounds of the given triangular facet data
!
do i=1,ntri
   do n=1,3
      do m=1,3
         xmax(m)=max(xmax(m),xsurf(m,ncf3(n,i)))
         xmin(m)=min(xmin(m),xsurf(m,ncf3(n,i)))
      enddo
   enddo
enddo
!
!     increase the boundaries by 5 percent for preventing
!     round off errors
!
fac=1.05
lenmax=0.
d0=1
do m=1,3
   xlen=(xmax(m)-xmin(m))*0.5*fac
   xcm=(xmax(m)+xmin(m))*0.5
   xmin(m)=xcm-xlen
   xmax(m)=xcm+xlen
   if (xlen.gt.lenmax) then
      lenmax=xlen
      d0=m
   endif
enddo
!
!     now find bounds of the coordinates which are contained
!     inside the bounding box
!

js=jmax
je=1
ks=kmax
ke=1
ls=lmax
le=1

do l=1,lmax
   zz=ds*(l-1)+xl(3)
   do k=1,kmax
      yy=ds*(k-1)+xl(2)
      do j=1,jmax
         xx=ds*(j-1)+xl(1)

         if ((xx-xmin(1))*(xx-xmax(1)).le.0.and.&
             (yy-xmin(2))*(yy-xmax(2)).le.0.and.&
             (zz-xmin(3))*(zz-xmax(3)).le.0) then

            je=max(j,je)
            js=min(j,js)
            ke=max(k,ke)
            ks=min(k,ks)
            le=max(l,le)
            ls=min(l,ls)
            
         endif
         
      enddo
   enddo
enddo
!
! if no intersection return
!
if (js.gt.je.or.ks.gt.ke.or.ls.gt.le) then
 deallocate(iblank)
 return
endif
!
! find the triangles which fall inside the Y-Z projection of 
! current cartesian box
!
allocate(iltri(ntri))

itri=0
do i=1,ntri
   vertexLoop1: do n=1,3

      kt=int((xsurf(2,ncf3(n,i))-xl(2))/ds)+1
      lt=int((xsurf(3,ncf3(n,i))-xl(3))/ds)+1
      
      if ((kt-ke-1)*(kt-ks+1).le.0 .and. &
           (lt-le-1)*(lt-ls+1).le.0) then
         itri=itri+1
         iltri(itri)=i
         exit vertexLoop1
      endif
   enddo vertexLoop1
enddo      
!
if (itri.eq.0) then
 deallocate(iltri)
 return
endif
!
!     preprocess the projected areas in Y-Z plane
!     note :: have to change to 
!             generalize and pick the plane with maximum
!             exposure
!
eps1=1e-10
eps2=1e-10
allocate(cof(5,itri),xd(itri))
areamax=0.

do ii=1,itri
   i=iltri(ii)
   cof(1,ii)=xsurf(2,ncf3(1,i))-xsurf(2,ncf3(3,i))
   cof(2,ii)=xsurf(2,ncf3(2,i))-xsurf(2,ncf3(3,i))
   cof(3,ii)=xsurf(3,ncf3(1,i))-xsurf(3,ncf3(3,i))
   cof(4,ii)=xsurf(3,ncf3(2,i))-xsurf(3,ncf3(3,i))
   den=cof(1,ii)*cof(4,ii)-cof(2,ii)*cof(3,ii)
   if (den.gt.areamax) areamax=den
!
!     if zero projected area then don't check this triangle
!
   if (den.ne.0) then
      cof(5,ii)=1./den
   else
      cof(5,ii)=0.
   endif
enddo
!
! allocate pointer arrays
!
maxintersect=40*(ke-ks+1)*(le-ls+1)
allocate(idd(maxintersect))
allocate(xdd(maxintersect))
allocate(ipointer(kmax,lmax))
ipointer=0
ixd=0
!
! find intersections for triangles in iltri
!
do ii=1,itri
   if (cof(5,ii).ne.0) then
      i=iltri(ii)
      ka=ke
      kb=ks
      la=le
      lb=ls

      do n=1,3
         kt=int((xsurf(2,ncf3(n,i))-xl(2))/ds)+1
         lt=int((xsurf(3,ncf3(n,i))-xl(3))/ds)+1
         
         ka=min(kt,ka)
         kb=max(kt,kb)
         
         la=min(lt,la)
         lb=max(lt,lb)
      enddo
   
      ka=max(ka,ks)
      kb=min(kb,ke)
      la=max(la,ls)
      lb=min(lb,le)
      
      do l=la,lb
         z1=xl(3)+(l-1)*ds
         zz=(xsurf(3,ncf3(3,i))-z1)
         do k=ka,kb
            y1=xl(2)+(k-1)*ds
            yy=(xsurf(2,ncf3(3,i))-y1)
            !     area coordinates
            l1=(cof(2,ii)*zz-cof(4,ii)*yy)*cof(5,ii)
            l2=(cof(3,ii)*yy-cof(1,ii)*zz)*cof(5,ii)
            l3=1-l1-l2
            !     if convex set the point is inside
            if (l1.ge.-eps1.and.l2.ge.-eps1.and.l3.ge.-eps1) then
               ixd=ixd+1
               xdd(ixd)=l1*xsurf(1,ncf3(1,i))+&
                    l2*xsurf(1,ncf3(2,i))+&
                    l3*xsurf(1,ncf3(3,i))
               idd(ixd)=ipointer(k,l)
               ipointer(k,l)=ixd
            endif
         enddo
      enddo
   endif
enddo
!
js1=jmax
je1=1
ks1=kmax
ke1=1
ls1=lmax
le1=1
!
! perform the actual iblanking
!
lloop: do l=ls,le
   kloop: do k=ks,ke
      
      itr=0
      ipp=ipointer(k,l)
      do while(ipp > 0) 
         itr=itr+1
         xd(itr)=xdd(ipp)
         ipp=idd(ipp)
      enddo
      
      itrChk: if (itr.gt.1) then
         !
         !     sort the distance vector
         !
         sortLoop: do ii=1,itr
            do jj=ii+1,itr
               if (xd(jj).lt.xd(ii)) then
                  tmp=xd(ii)
                  xd(ii)=xd(jj)
                  xd(jj)=tmp
               endif
            enddo
         enddo sortLoop
         !
         !     cull out same distances
         !
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
         
         ii=1
         j=js
         inside=0

         js1=js
         je1=je
         ks1=ks
         ke1=ke
         ls1=ls
         le1=le
         !
         !     iblank if inside
         !
         do while(ii.lt.itr)
            j1=max(ceiling((xd(ii)-xl(1))/ds)+1,js)
            j2=min(floor((xd(ii+1)-xl(1))/ds)+1,je)
            if (j2.ge.j1) then
               iblank(j1:j2,k,l)=holeFlag
               ibcount=ibcount+(j2-j1+1)

               js1=min(js1,j1)
               je1=max(je1,j2)
               ks1=min(ks1,k)
               ke1=max(ke1,k)
               ls1=min(ls1,l)
               le1=max(le1,l)
            endif
            ii=ii+2
         enddo

      endif itrChk
   enddo kloop
enddo lloop
!
nf1=nfringe-1
nf2=2*nfringe
jdim=(jmax+nf2)
jkdim=(kmax+nf2)*jdim
!
do l=1,lmax
   do k=1,kmax
      do j=1,jmax
         ii=(l+nf1)*jkdim+(k+nf1)*jdim+j+nfringe
	 if (iblank(j,k,l).eq.holeFlag) then
            ibl(ii)=iblank(j,k,l)
            holeCount=holeCount+1
         endif
      enddo
   enddo
enddo
!
! deallocate local memory
!
deallocate(cof,xd)
deallocate(iltri)
deallocate(iblank)
deallocate(ipointer,xdd,idd)
!write(6,*) 'ibcount=',ibcount

return
end subroutine holeCutSurf


