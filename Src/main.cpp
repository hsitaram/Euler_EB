#include "solver.h"
//***********************************************
int main (int argc, char* argv[])
{
    Real pos[BL_SPACEDIM];
    BoxLib::Initialize(argc,argv);
    std::cout<<"Hello\n";
	
    pos[0]=0.5; pos[1]=0.5; pos[2]=0.5;
    solver obj;
    obj.setparams(0.2,pos);

    obj.readinputs();
    obj.makefabs();
    obj.intialize_vars();
    std::cout<<"Initialized solution\n";

    std::cout<<"starting time stepper\n";
    obj.timestepping();

    BoxLib::Finalize();
    return 0;
}
