#include"ebmesher.h"
#include"eb_fortran_functions.h"

//==================================================================================
void ebmesher::read_tri_file()
{
	char hash;
	std::string temp;
	double x,y,z;

	std::ifstream infile("surface_triangulation.plt");
	
	infile>>hash>>num_of_nodes>>num_of_triangles;
	infile>>temp;   //title..
	infile>>temp;   //variables..
	infile>>temp;   //zone..
	std::cout<<"tri nums:"<<num_of_nodes<<"\t"<<num_of_triangles<<"\n";

	surface_coordinates = new Real[DIM3*num_of_nodes];
	connectivity        = new int[DIM3*num_of_triangles];

	read_triangulation_file(surface_coordinates,connectivity,&num_of_nodes,&num_of_triangles);

	infile.close();
}
//==================================================================================
void ebmesher::update_volfrac(Geometry &geom,int &nghost,int ncx,int ncy,int ncz)
{
	int glo[DIM3],ghi[DIM3];
	int solidpoints,ng_global;

	glo[XDIR]=0;      glo[YDIR]=0;      glo[ZDIR]=0;
	ghi[XDIR]=ncx-1;  ghi[YDIR]=ncy-1;  ghi[ZDIR]=ncz-1;

	std::cout<<"updating volfrac\n";

#ifndef EB_SPHERE
	read_tri_file();
#endif

	for(MFIter mfi(*eb_volfrac);mfi.isValid();++mfi)
	{
		const Box& bx = mfi.validbox();

#ifdef EB_SPHERE
		std::cout<<"using eb sphere\n";
		find_volfrac_in_box_forsphere((*eb_volfrac)[mfi].dataPtr(),
				bx.loVect(), bx.hiVect(), &nghost,
				geom.CellSize(), geom.ProbLo(), geom.ProbHi(),
				&rad,pos);
#else
		find_volfrac_in_box((*eb_volfrac)[mfi].dataPtr(),
				bx.loVect(), bx.hiVect(),&nghost,
				geom.CellSize(), geom.ProbLo(), geom.ProbHi(),
				surface_coordinates,connectivity,&num_of_nodes,
				&num_of_triangles);
#endif

	}
} 
//==================================================================================
