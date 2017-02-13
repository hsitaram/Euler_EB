module paraviewfilewrite

	implicit none

	integer, parameter :: DOUBLESIZE=8
	integer, parameter :: INTSIZE=4
	integer, parameter :: NDIM=3
	integer, parameter :: NPTS_PER_TRI  = 3
	integer, parameter :: NPTS_PER_QUAD = 4
	integer, parameter :: NPTS_PER_HEX  = 8
	integer, parameter :: NPTS_PER_TET  = 4
	integer, parameter :: ZERO=0
	integer, parameter :: MAXSTRSIZE=500
	integer, parameter :: SMALLSTRSIZE=8

	!VTK codes for shapes
	integer, parameter :: CODE_TRI  = 5  !triangle
	integer, parameter :: CODE_QUAD = 9  !quadrilateral
	integer, parameter :: CODE_TET  = 10 !tetrahedron
	integer, parameter :: CODE_HEX  = 12 !hexahedron
	integer, parameter :: CODE_WEDG = 13 !tri prism
	integer, parameter :: CODE_PRMD = 14 !pyramid



	contains

!===============================================================================================
subroutine Writecartmeshfile(fname,solnvars,nscalars,scalarnames,&
				Nx,Ny,Nz,dx,dy,dz,offx,offy,offz)

character(LEN=*), intent(in) :: fname
integer, intent(in) :: Nx,Ny,Nz
integer, intent(in) :: nscalars
character(LEN=*), intent(in) :: scalarnames(nscalars)
real*8, intent(in) :: solnvars(Nx,Ny,Nz,nscalars)
real*8, intent(in) :: dx,dy,dz
real*8, intent(in) :: offx,offy,offz

character(LEN=MAXSTRSIZE) :: str
character(LEN=SMALLSTRSIZE) :: offset
character(LEN=SMALLSTRSIZE) :: xmin_str
character(LEN=SMALLSTRSIZE) :: ymin_str
character(LEN=SMALLSTRSIZE) :: zmin_str
character(LEN=SMALLSTRSIZE) :: xmax_str
character(LEN=SMALLSTRSIZE) :: ymax_str
character(LEN=SMALLSTRSIZE) :: zmax_str
character(LEN=MAXSTRSIZE) :: global_extent
character(LEN=MAXSTRSIZE) :: local_extent

character :: slashn
integer :: filenum
integer :: i,j,k
integer :: byte_offset
integer :: scalars_size
real*8,allocatable :: local_xc(:),local_yc(:),local_zc(:)

filenum = 9
slashn  = char(10)


allocate(local_xc(Nx))
allocate(local_yc(Ny))
allocate(local_zc(Nz))
	
do i=1,Nx
	local_xc(i) = offx + (i-1)*dx
enddo

do i=1,Ny
	local_yc(i) = offy + (i-1)*dy
enddo

do i=1,Nz
	local_zc(i) = offz + (i-1)*dz
enddo

scalars_size   = Nx*Ny*Nz * DOUBLESIZE

write(xmin_str(1:8),'(i8)') 0
write(ymin_str(1:8),'(i8)') 0
write(zmin_str(1:8),'(i8)') 0
write(xmax_str(1:8),'(i8)') Nx-1
write(ymax_str(1:8),'(i8)') Ny-1
write(zmax_str(1:8),'(i8)') Nz-1

global_extent = trim(xmin_str)//' '//trim(xmax_str)//' '//trim(ymin_str)//' '&
		//trim(ymax_str)//trim(zmin_str)//' '//trim(zmax_str)

write(xmin_str(1:8),'(i8)') 0
write(ymin_str(1:8),'(i8)') 0
write(zmin_str(1:8),'(i8)') 0
write(xmax_str(1:8),'(i8)') Nx-1
write(ymax_str(1:8),'(i8)') Ny-1
write(zmax_str(1:8),'(i8)') Nz-1

local_extent = trim(xmin_str)//' '//trim(xmax_str)//' '//trim(ymin_str)//' '&
		//trim(ymax_str)//trim(zmin_str)//' '//trim(zmax_str)

byte_offset = 0                             
open(unit=filenum,file=fname,form='unformatted',access='stream')

!=========================================================
str = '<?xml version="1.0"?>'//slashn                                                                                                   
write(filenum) trim(str)
str = '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'//slashn		                               
write(filenum) trim(str)
str = '  <RectilinearGrid WholeExtent="'//trim(global_extent)//'">'//slashn                                       
write(filenum) trim(str)
str = '    <Piece Extent="'//trim(local_extent)//'">'//slashn
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '      <PointData>  '//slashn
write(filenum) trim(str)
do i=1,nscalars

	write(offset(1:8),'(i8)') byte_offset
	str = '         <DataArray type="Float64" Name="'//trim(scalarnames(i))//&
		'" format="appended" offset="'//offset//'"           />'//slashn
	write(filenum) trim(str)
	byte_offset = byte_offset + intsize + scalars_size

enddo
str = '      </PointData>'//slashn            
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '      <CellData> '//slashn                                                       
write(filenum) trim(str)
str = '      </CellData>'//slashn            
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '      <Coordinates>'//slashn                                                                                   
write(filenum) trim(str)
write(offset(1:8),'(i8)') byte_offset

str = '<DataArray type="Float64" Name="X"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + Nx*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Y"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + Ny*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Z"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

str = '      </Coordinates>'//slashn       
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '    </Piece>'//slashn                                                                                                  
write(filenum) trim(str)
str = '  </RectilinearGrid>'//slashn                                                                                                   
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '  <AppendedData encoding="raw">'//slashn                                                                                         
write(filenum) trim(str)
str = '_'                                                                                                                           
write(filenum) trim(str)
!=========================================================

!=========================================================
do i=1,nscalars
	write(filenum) scalars_size  ,  solnvars(1:Nx,1:Ny,1:Nz,i)
enddo
!=========================================================

!=========================================================
write(filenum) Nx*DOUBLESIZE  , local_xc(:)
write(filenum) Ny*DOUBLESIZE  , local_yc(:)
write(filenum) Nz*DOUBLESIZE  , local_zc(:)
!=========================================================

!=========================================================
str = slashn//'  </AppendedData>'//slashn                                                                             
write(filenum) trim(str)
str = '</VTKFile>'//slashn                                                                                                             
write(filenum) trim(str)
!=========================================================

close(filenum)
deallocate(local_xc,local_yc,local_zc)

end subroutine Writecartmeshfile
!===============================================================================================
subroutine Writeumeshfile(fname,nodes,connectivity,nnodes,ncells)

implicit none

character(LEN=*),intent(in)  :: fname
integer, intent(in)   ::  nnodes,ncells
real*8,intent(in)     ::  nodes(NDIM,nnodes)
integer,intent(in)    ::  connectivity(NDIM,ncells)

character(LEN=MAXSTRSIZE) :: str
character(LEN=SMALLSTRSIZE) :: offset
character(LEN=SMALLSTRSIZE) :: str1,str2
character :: slashn
integer   :: filenum = 9
integer   :: npts_per_cell
integer   :: cell_code

integer :: i,j,k
integer :: byte_offset,intval

slashn = char(10)
byte_offset = 0                            
npts_per_cell = NPTS_PER_TRI
cell_code     = CODE_TRI

open(unit=filenum,file=fname,form='unformatted',access='stream')

str = '<?xml version="1.0"?>'//slashn                                                                                                   
write(filenum) trim(str)

str = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//slashn		                               
write(filenum) trim(str)

str = '  <UnstructuredGrid>'//slashn                                       
write(filenum) trim(str)

write(str1(1:8),'(i8)') nnodes
write(str2(1:8),'(i8)') ncells
str = '    <Piece NumberOfPoints="'//str1//'" NumberOfCells="'//str2//'">'//slashn
write(filenum) trim(str)

!IF YOU WANT TO ADD SCALARS AT POINTS========================================
!str = '      <PointData> '//slashn                                                       
!write(filenum) trim(str)

!write(offset(1:8),'(i8)') byte_offset
!str = '         <DataArray type="Float64" Name="scal1" format="appended" offset="'//offset//'"           />'//slashn
!write(filenum) trim(str)

!str = '      </PointData>'//slashn            
!write(filenum) trim(str)
!=============================================================================

!IF YOU WANT TO ADD SCALARS AT CELL CENTERS===================================
!str = '      <CellData>  </CellData>'//slashn
!write(filenum) trim(str)
!=============================================================================
!=============================================================================
!give appropriate offsets if there are scalars
!byte_offset = byte_offset + INTSIZE + scalars_size
!=============================================================================

write(offset(1:8),'(i8)') byte_offset
str = '      <Points>'//slashn                                                                                   
write(filenum) trim(str)
str = '<DataArray type="Float64" Name="Pcoord" NumberOfComponents="3" format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)
str = '      </Points>'//slashn       
write(filenum) trim(str)

str = '      <Cells>'//slashn                                                                                                 
write(filenum) trim(str)
byte_offset = byte_offset + INTSIZE + NDIM*nnodes*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '        <DataArray type="Int32" Name="connectivity" format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + INTSIZE + npts_per_cell*INTSIZE*ncells 
write(offset(1:8),'(i8)') byte_offset
str = '        <DataArray type="Int32" Name="offsets" format="appended" offset="'//offset//'" />'//slashn                  
write(filenum) trim(str)

byte_offset = byte_offset + INTSIZE + INTSIZE*ncells
write(offset(1:8),'(i8)') byte_offset
str = '        <DataArray type="Int32" Name="types" format="appended" offset="'//offset//'" />'//slashn                       
write(filenum) trim(str)

str = '      </Cells>'//slashn                         
write(filenum) trim(str)

str = '    </Piece>'//slashn                                                                                                  
write(filenum) trim(str)

str = '  </UnstructuredGrid>'//slashn                                                                                                   
write(filenum) trim(str)

str = '  <AppendedData encoding="raw">'//slashn                                                                                         
write(filenum) trim(str)

str = '_'                                                                                                                           
write(filenum) trim(str)

!ADD SCALAR DATA====================================================
!write(filenum) scalars_size  , scalars(:)
!===================================================================
write(filenum) NDIM*nnodes*DOUBLESIZE   , nodes(:,:)

!writing connectivity
write(filenum) npts_per_cell*INTSIZE*ncells
do i=1,ncells
	do j=1,npts_per_cell

		!node numbers start from 0, not 1
		intval = connectivity(j,i)-1
		write(filenum) intval

	enddo
enddo

!writing offsets (goes as 8,16,24..) 
!paraview will read connectivities as 1-8,8-16 and so on.
write(filenum) INTSIZE*ncells
do i=1,ncells
	intval=i*npts_per_cell
	write(filenum) intval
enddo


!writing types of cells
write(filenum) INTSIZE*ncells
intval=cell_code !cell type code for hex
do i=1,ncells
	write(filenum) intval
enddo

str = slashn//'  </AppendedData>'//slashn                                                                             
write(filenum) trim(str)

str = '</VTKFile>'//slashn                                                                                                             
write(filenum) trim(str)

close(filenum)

end subroutine Writeumeshfile
!===============================================================================================
subroutine Writepointfile(fname,npoints,coords,pointdata)

implicit none

character(LEN=*),intent(in)  :: fname
integer, intent(in)   ::  npoints
real*8,intent(in)     ::  coords(NDIM,npoints)
real*8,intent(in)     ::  pointdata(npoints)

character(LEN=MAXSTRSIZE) :: str
character(LEN=SMALLSTRSIZE) :: offset
character(LEN=SMALLSTRSIZE) :: str1,str2
character :: slashn
integer   :: filenum = 9
integer   :: npts_per_cell
integer   :: cell_code

integer :: i,j,k
integer :: byte_offset,intval

slashn = char(10)
byte_offset = 0                            

open(unit=filenum,file=fname,form='unformatted',access='stream')

!=========================================================
str = '<?xml version="1.0"?>'//slashn                                                                                                   
write(filenum) trim(str)
str = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//slashn		                               
write(filenum) trim(str)
str = '  <UnstructuredGrid>'//slashn                                       
write(filenum) trim(str)
write(str1(1:8),'(i8)') npoints
str = '    <Piece NumberOfPoints="'//str1//'" NumberOfCells="0">'//slashn
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '      <PointData> '//slashn                                                       
write(filenum) trim(str)
write(offset(1:8),'(i8)') byte_offset
str = '         <DataArray type="Float64" Name="iblank" format="appended" offset="'//offset//'"           />'//slashn
write(filenum) trim(str)
str = '      </PointData>'//slashn            
write(filenum) trim(str)
!=========================================================

!=========================================================
byte_offset = byte_offset + INTSIZE + npoints*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '      <Points>'//slashn                                                                                   
write(filenum) trim(str)
str = '<DataArray type="Float64" Name="Pcoord" NumberOfComponents="3" format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)
str = '      </Points>'//slashn       
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '      <Cells>'//slashn                                                                                                 
write(filenum) trim(str)
byte_offset = byte_offset + INTSIZE + NDIM*npoints*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '        <DataArray type="Int32" Name="connectivity" format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + INTSIZE
write(offset(1:8),'(i8)') byte_offset
str = '        <DataArray type="Int32" Name="offsets" format="appended" offset="'//offset//'" />'//slashn                  
write(filenum) trim(str)

byte_offset = byte_offset + INTSIZE
write(offset(1:8),'(i8)') byte_offset
str = '        <DataArray type="Int32" Name="types" format="appended" offset="'//offset//'" />'//slashn                       
write(filenum) trim(str)
str = '      </Cells>'//slashn                         
write(filenum) trim(str)
!=========================================================


!=========================================================
str = '    </Piece>'//slashn                                                                                                  
write(filenum) trim(str)
str = '  </UnstructuredGrid>'//slashn                                                                                                   
write(filenum) trim(str)
!=========================================================

!=========================================================
str = '  <AppendedData encoding="raw">'//slashn                                                                                         
write(filenum) trim(str)
str = '_'                                                                                                                           
write(filenum) trim(str)
!=========================================================

!=========================================================
write(filenum) npoints*DOUBLESIZE  , pointdata(:)
!=========================================================

!=========================================================
write(filenum) NDIM*npoints*DOUBLESIZE   , coords(:,:)
!=========================================================

!=========================================================
write(filenum) ZERO
write(filenum) ZERO
write(filenum) ZERO
!=========================================================

!=========================================================
str = slashn//'  </AppendedData>'//slashn                                                                             
write(filenum) trim(str)
str = '</VTKFile>'//slashn                                                                                                             
write(filenum) trim(str)
!=========================================================

close(filenum)

end subroutine Writepointfile
!===============================================================================================

end module paraviewfilewrite
