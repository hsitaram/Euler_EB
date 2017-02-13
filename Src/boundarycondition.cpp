#include"boundarycondition.h"

//===============================================================
void boundarycondition::assignbc(Array<Real> bcarray,int dir)
{
	for(int i=0;i<NCVARS+1;i++)
	{
		bcparams[dir][i]=bcarray[i];
	}
}
//===============================================================
void boundarycondition::readbc()
{
	ParmParse pp;
        Array<Real> bcarray(NCVARS+1);

	for(int i=0;i<NB;i++)
	{
		for(int j=0;j<NCVARS+1;j++)
		{
			bcparams[i][j]=0.0;
		}
	}	
		
	pp.queryarr("Left_bc"  ,bcarray,0,NCVARS+1);
	assignbc(bcarray,LFT);
	pp.queryarr("Right_bc" ,bcarray,0,NCVARS+1);
	assignbc(bcarray,RGT);

	pp.queryarr("Bottom_bc",bcarray,0,NCVARS+1);
	assignbc(bcarray,BOT);
	pp.queryarr("Top_bc"   ,bcarray,0,NCVARS+1);
	assignbc(bcarray,TOP);

	pp.queryarr("Back_bc"  ,bcarray,0,NCVARS+1);
	assignbc(bcarray,BCK);
	pp.queryarr("Front_bc" ,bcarray,0,NCVARS+1);
	assignbc(bcarray,FRT);

	std::cout<<"bcarray:"<<bcarray[0]<<"\n";

	set_fortran_bc();
}
//===============================================================
void boundarycondition::set_fortran_bc()
{
	//note: why separate variables? not sure how 2d arrays will
	//behave from C to fortran
	//I can guess!
	Real dens[NB];
	Real velx[NB];
	Real vely[NB];
	Real velz[NB];
	Real p[NB];

	for(int i=0;i<NB;i++)
	{
		dens[i] = bcparams[i][RHO_INDEX];
		velx[i] = bcparams[i][RHO_U_INDEX];
		vely[i] = bcparams[i][RHO_V_INDEX];
		velz[i] = bcparams[i][RHO_W_INDEX];
		p[i]    = bcparams[i][RHO_E_INDEX];
	}

	set_dirichlet_bcs(dens,velx,vely,velz,p);
}
//===============================================================
void boundarycondition::printbcparams()
{
	for(int i=0;i<NB;i++)
	{
		std::cout<<"side "<<i<<":\t";
		for(int j=0;j<NCVARS+1;j++)
		{
			std::cout<<bcparams[i][j]<<"\t";
		}
		std::cout<<"\n";
	}
}
//===============================================================
void boundarycondition::assign_solver_bcs(Array<int> &lo_bc,Array<int> &hi_bc)
{
	//these are solver specific bcs defined in constantvars
	//and not boxlib bcs like EXT_DIR, FOEXTRAP etc..

	lo_bc[XDIR]=int(bcparams[LFT][0]);
	hi_bc[XDIR]=int(bcparams[RGT][0]);

	lo_bc[YDIR]=int(bcparams[BOT][0]);
	hi_bc[YDIR]=int(bcparams[TOP][0]);

	lo_bc[ZDIR]=int(bcparams[BCK][0]);
	hi_bc[ZDIR]=int(bcparams[FRT][0]);
}
//===============================================================
void boundarycondition::updatebc(MultiFab *field,PhysBCFunct &physbcf,
Geometry &geom,Real time)
{
	field->FillBoundary(geom.periodicity());
	physbcf.FillBoundary(*field,time);
}
//===============================================================
void boundarycondition::assign_bl_bcs(std::string fieldname,Array<int>& lo_bc,
Array<int>& hi_bc)
{
	lo_bc[XDIR] = my_bl_bc(fieldname,int(bcparams[LFT][0]));
	hi_bc[XDIR] = my_bl_bc(fieldname,int(bcparams[RGT][0]));

	lo_bc[YDIR] = my_bl_bc(fieldname,int(bcparams[BOT][0]));
	hi_bc[YDIR] = my_bl_bc(fieldname,int(bcparams[TOP][0]));

	lo_bc[ZDIR] = my_bl_bc(fieldname,int(bcparams[BCK][0]));
	hi_bc[ZDIR] = my_bl_bc(fieldname,int(bcparams[FRT][0]));
}
//===============================================================
int boundarycondition::my_bl_bc(std::string fieldname,int solver_bc)
{
	int bctype;

	if(fieldname=="density")
	{
		switch(solver_bc)
		{
			case BC_SUPINLET:      {bctype = EXT_DIR;   break;}
			case BC_SUBINLET:      {bctype = FOEXTRAP;  break;}
			case BC_SUPOUTLET:     {bctype = FOEXTRAP;  break;}
			case BC_SUBOUTLET:     {bctype = FOEXTRAP;  break;}
			case BC_WALL:          {bctype = FOEXTRAP;  break;}
			case BC_PERIODIC:      {bctype = INT_DIR;   break;}
			default: {bctype = INT_DIR; break;}
		}
	}
	else if((fieldname=="vel_x") || (fieldname=="vel_y") || (fieldname=="vel_z"))
	{
		switch(solver_bc)
		{
			case BC_SUPINLET:      {bctype = EXT_DIR;  break;}
			case BC_SUBINLET:      {bctype = EXT_DIR;  break;}
			case BC_SUPOUTLET:     {bctype = FOEXTRAP; break;}
			case BC_SUBOUTLET:     {bctype = FOEXTRAP; break;}
			case BC_WALL:          {bctype = EXT_DIR;  break;}
			case BC_PERIODIC:      {bctype = INT_DIR;  break;}
			default: {bctype = INT_DIR; break;}
		}

	}	
	else if(fieldname=="pressure")
	{
		switch(solver_bc)
		{
			case BC_SUPINLET:  {bctype = EXT_DIR;  break;}
			case BC_SUBINLET:  {bctype = EXT_DIR;  break;}
			case BC_SUPOUTLET: {bctype = FOEXTRAP; break;}
			case BC_SUBOUTLET: {bctype = EXT_DIR;  break;}
			case BC_WALL:      {bctype = FOEXTRAP; break;}
			case BC_PERIODIC:  {bctype = INT_DIR;  break;}
			default: {bctype = INT_DIR; break;}
		}
	
	}
	else
	{
		std::cerr<<"Field name "<<fieldname<<" does not exist\n";
		exit(0);
	}

	return(bctype);
}
//===============================================================
