!====================================================================
!
!   Routine to cut holes of a uniform isotropic cartesian mesh
!   using a unstructured triangle surface as the cutter
!   
!   Author: Jay Sitaraman
!   Created :: 08/21/07 
!   Modified for pure iblanking on:: 12/03/09
!   Adapted by Hari Sitaraman for BoxLib : 12/13/16
!
!====================================================================
subroutine holeCutSurf(glo,ghi,lo,hi,xl,ds,ng,xsurf,ncf3,nnodes,ntri,solidPoints,vfrac)
!
implicit none
include 'f_constants.h'
!
! subroutine arguments
!
integer, intent(in) :: glo(THREEDIM),ghi(THREEDIM) !global extents of domain
integer,intent(in) :: lo(THREEDIM),hi(THREEDIM),ng,nnodes,ntri                    ! dimensions of cart grid, number of fringes
double precision, intent(in)    :: xl(THREEDIM)                                   ! lower-left-bottom corner of cartesian grid
double precision, intent(in)    :: ds(THREEDIM)                                   ! size of cartesian grid cell 
                                                                                  ! number of nodes and triangles in the cutter mesh
integer, intent(inout) :: solidPoints                                               ! number of hole points
double precision, intent(inout)  :: vfrac(lo(DIRX)-ng:hi(DIRX)+ng,& 
				lo(DIRY)-ng:hi(DIRY)+ng,lo(DIRZ)-ng:hi(DIRZ)+ng)  ! iblank array that masks all the point inside
                                                                                  ! with provided holeFlag
double precision, intent(in)    :: xsurf(THREEDIM,nnodes)                         ! surface grid nodes 
integer    :: ncf3(THREEDIM,ntri)                                                 ! surface grid connectivity 
!
! local variables
!
double precision, allocatable :: xdd(:)
integer, allocatable :: idd(:)
integer :: ipointer(glo(DIRY)-ng:ghi(DIRY)+ng,glo(DIRZ)-ng:ghi(DIRZ)+ng)
integer, allocatable :: iltri(:)
double precision, allocatable :: cof(:,:),xd(:)
!
integer :: ibcount,maxintersect,ipp,ixd
integer :: js,je,ks,ke,ls,le
integer :: js1,je1,ks1,ke1,ls1,le1
integer :: ka,kb,la,lb
integer :: itri
integer :: i,j,k,l,m,n,i0,j1,j2,ii,jj,itr,d0,inside
integer :: jt,kt,lt
double precision   :: l1,l2,l3
double precision   :: xmax(THREEDIM),xmin(THREEDIM),areamax
double precision   :: lenmax,xx,yy,zz,y1,z1,tmp,fac,eps1,eps2,den,xlen,xcm
integer :: jdim,jkdim,nf1,nf2
double precision :: vertexsum
integer :: j1_local,j2_local
!
! begin
!
!write(6,*) 'tri1=',ncf3(:,1)
ibcount=0
solidPoints=0
vfrac=1.d0
!
! initialize local iblanking
!
do m=1,THREEDIM
   xmin(m)=xsurf(m,ncf3(1,1))
   xmax(m)=xmin(m)
enddo
!
!     find bounds of the given triangular facet data
!
do i=1,ntri
   do n=1,THREEDIM
      do m=1,THREEDIM
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
do m=1,THREEDIM
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

!note: the +1 is for number of points instead of cells
js=ghi(DIRX)+1
je=glo(DIRX)

ks=ghi(DIRY)+1
ke=glo(DIRY)

ls=ghi(DIRZ)+1
le=glo(DIRZ)

do l=glo(DIRZ),ghi(DIRZ)+1
   zz=ds(DIRZ)*dble(l)+xl(DIRZ)

   do k=glo(DIRY),ghi(DIRY)+1
      yy=ds(DIRY)*dble(k)+xl(DIRY)

      do j=glo(DIRX),ghi(DIRX)+1
         xx=ds(DIRX)*dble(j)+xl(DIRX)

         if ((xx-xmin(DIRX))*(xx-xmax(DIRX)).le.0.and.&
             (yy-xmin(DIRY))*(yy-xmax(DIRY)).le.0.and.&
             (zz-xmin(DIRZ))*(zz-xmax(DIRZ)).le.0) then

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
 return
endif
!
! find the triangles which fall inside the Y-Z projection of 
! current cartesian box
!
allocate(iltri(ntri))
!print *,"here",glo,ghi,lo,hi
itri=0
do i=1,ntri
   vertexLoop1: do n=1,THREEDIM

      kt=int((xsurf(DIRY,ncf3(n,i))-xl(DIRY))/ds(DIRY))
      lt=int((xsurf(DIRZ,ncf3(n,i))-xl(DIRZ))/ds(DIRZ))
      
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

   cof(1,ii)=xsurf(DIRY,ncf3(1,i)) - xsurf(DIRY,ncf3(3,i))
   cof(2,ii)=xsurf(DIRY,ncf3(2,i)) - xsurf(DIRY,ncf3(3,i))

   cof(3,ii)=xsurf(DIRZ,ncf3(1,i)) - xsurf(DIRZ,ncf3(3,i))
   cof(4,ii)=xsurf(DIRZ,ncf3(2,i)) - xsurf(DIRZ,ncf3(3,i))

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
ipointer=0
ixd=0
xdd=0
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

      do n=1,THREEDIM
         kt=int((xsurf(DIRY,ncf3(n,i))-xl(DIRY))/ds(DIRY))
         lt=int((xsurf(DIRZ,ncf3(n,i))-xl(DIRZ))/ds(DIRZ))
         
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
         z1=xl(DIRZ)+dble(l)*ds(DIRZ)
         zz=(xsurf(DIRZ,ncf3(3,i))-z1)

         do k=ka,kb
            y1=xl(DIRY)+dble(k)*ds(DIRY)
            yy=(xsurf(DIRY,ncf3(3,i))-y1)

            !     area coordinates
            l1=(cof(2,ii)*zz-cof(4,ii)*yy)*cof(5,ii)
            l2=(cof(3,ii)*yy-cof(1,ii)*zz)*cof(5,ii)
            l3=1-l1-l2
            !     if convex set the point is inside
            if (l1.ge.-eps1.and.l2.ge.-eps1.and.l3.ge.-eps1) then

               ixd=ixd+1
               xdd(ixd) = l1*xsurf(DIRX,ncf3(1,i))+&
                          l2*xsurf(DIRX,ncf3(2,i))+&
                          l3*xsurf(DIRX,ncf3(3,i))

               idd(ixd)=ipointer(k,l)
               ipointer(k,l)=ixd

            endif
         enddo
      enddo
   endif
enddo
!print *,"here",glo,ghi,lo,hi
!
js1=ghi(DIRX)+1
je1=glo(DIRX)

ks1=ghi(DIRY)+1
ke1=glo(DIRY)

ls1=ghi(DIRZ)+1
le1=glo(DIRZ)

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

            j1=max(ceiling((xd(ii)-xl(DIRX))/ds(DIRX)),js)
            j2=min(floor((xd(ii+1)-xl(DIRX))/ds(DIRX)),je)

            if (j2.ge.j1) then
               
               j1_local=max(j1,lo(DIRX))
               j2_local=min(j2,hi(DIRX))

               !print *,j1_local,j2_local,j1,j2

               if(j2_local .ge. j1_local) then

                    if( (k .ge. lo(DIRY)) .and. (k .le. hi(DIRY)) .and. (l .ge. lo(DIRZ)) .and. (l .le. hi(DIRZ)) ) then
                        !print *,"j,k,l:",j1_local,j2_local,k,l,lo,hi
                        vfrac(j1_local:j2_local,k,l)=0.d0
                    endif
               endif

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
! deallocate local memory
!
deallocate(cof,xd)
deallocate(iltri)
deallocate(xdd,idd)
!write(6,*) 'ibcount=',ibcount

return
end subroutine holeCutSurf
