#include <fstream>
#include <iomanip>
#include <ParmParse.H>
#include <Geometry.H>
#include <VisMF.H>
#include <PhysBCFunct.H>
#include "c_constants.h"
#include "fortran_functions.h" 
class boundarycondition
{
	private:
		void assignbc(Array<Real> bcarray,int dir);
		int my_bl_bc(std::string fieldname,int solver_bc);
		void set_fortran_bc();
	public:
		Real bcparams[NB][NCVARS+1]; //1 extra for the kind of bc
		void readbc();
		void printbcparams();
		void assign_solver_bcs(Array<int> &lo_bc,Array<int> &hi_bc);
		void assign_bl_bcs(std::string fieldname,Array<int> &lo_bc,Array<int> &hi_bc);
		void updatebc(MultiFab *field,PhysBCFunct &physbcf,Geometry &geom,Real time);
};
