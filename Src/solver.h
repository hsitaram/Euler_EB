#include <fstream>
#include <iomanip>
#include <ParmParse.H>
#include <Geometry.H>
#include <VisMF.H>
#include <PhysBCFunct.H>
#include "writePlotFile.h"
#include "fortran_functions.h"
#include "ebmesher.h"
#include "boundarycondition.h" 
class solver: public ebmesher,public boundarycondition
{
	private:
		int n_cell_x,n_cell_y,n_cell_z,max_grid_size;
		int stdoutstepnum;
		int nsteps,nghost,ncomp;
		const Real *dx;
		Real vars_init[NCVARS];
		Geometry geom;
		BoxArray ba;

		MultiFab *consvars;
		MultiFab *consvars_nm1;
		MultiFab *dens;
		MultiFab *velx;
		MultiFab *vely;
		MultiFab *velz;
		MultiFab *pres;

		MultiFab *residual;

		Real lcorner[BL_SPACEDIM];
		Real ucorner[BL_SPACEDIM];

		Real tfinal;
		Real delt,cfl;
		bool cflcontrolled;
		Real dt_out;
		int temporal_order;

		std::vector<PhysBCFunct> physbcf; //note: should I try PArray?

		std::vector<std::string> varnames;

		void writeoutputfile(int n,double time);
		void updateResidual();
                void update_conserved_variables();
                void update_primitive_variables();
		void update_primvar_bcs(Real current_time);
		double find_timestep();
		
		int getvarindex(std::string vname);

		void advance_conserved_variables(Real current_time,Real &timestep,int nstep);

		bool sod_stube_flag;

	public:
		void timestepping();
		void readinputs();
		void makefabs();
		void intialize_vars();

};
