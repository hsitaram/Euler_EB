#include <fstream>
#include <iomanip>
#include <ParmParse.H>
#include <Geometry.H>
#include <VisMF.H>
#include "writePlotFile.h"
#include"c_constants.h"
class ebmesher
{
	private:
		//triangulation parameters
		Real *surface_coordinates;
		int *connectivity;
		int num_of_nodes,num_of_triangles;

		void read_tri_file();

	public:
		MultiFab *eb_volfrac;
		Real rad;
		Real pos[BL_SPACEDIM];
		void setparams(Real sphere_radius,Real sphere_position[BL_SPACEDIM])
		{
			rad=sphere_radius;
			
			pos[XDIR] = sphere_position[XDIR];
			pos[YDIR] = sphere_position[YDIR];
			pos[ZDIR] = sphere_position[ZDIR];
		}

		void update_volfrac(Geometry &geom,int &nghost,int ncellx,int ncelly,int ncellz);
};
