#include"solver.h"

//===============================================================
void solver::readinputs()
{
	int sodflag;
	ParmParse pp;
	Array<Real> lwcorner(BL_SPACEDIM),upcorner(BL_SPACEDIM);

	pp.get("n_cell_x",n_cell_x);
	pp.get("n_cell_y",n_cell_y);
	pp.get("n_cell_z",n_cell_z);
	pp.get("max_grid_size",max_grid_size);

	nsteps=0;
	pp.query("nsteps",nsteps);

	lwcorner[0]=0.0; lwcorner[1]=0.0;
	upcorner[0]=1.0; upcorner[1]=1.0;

#if (BL_SPACEDIM == 3)
	lwcorner[2]=0.0;
	upcorner[2]=1.0;
#endif

	pp.queryarr("lower_corner",lwcorner,0,BL_SPACEDIM);
	pp.queryarr("upper_corner",upcorner,0,BL_SPACEDIM);

	lcorner[0]=lwcorner[0]; lcorner[1]=lwcorner[1];
	ucorner[0]=upcorner[0]; ucorner[1]=upcorner[1];

#if (BL_SPACEDIM == 3)
	lcorner[2]=lwcorner[2];
	ucorner[2]=upcorner[2];
#endif

	nghost=1;
	pp.query("nghost",nghost);
	
	ncomp=1;
	pp.query("ncomp",ncomp);


	pp.query("density_initial",vars_init[0]);
	pp.query("velx_initial",vars_init[1]);
	pp.query("vely_initial",vars_init[2]);
	pp.query("velz_initial",vars_init[3]);
	pp.query("pressure_initial",vars_init[4]);

	pp.query("final_time",tfinal);
	pp.query("time_step",delt);
	pp.query("cfl_num",cfl);
	pp.query("write_time",dt_out);

	stdoutstepnum=1;
	pp.query("stdout_step_num",stdoutstepnum);

	temporal_order=1;
	pp.query("time_stepping_order",temporal_order);

	if(delt > 0.0)
	{
		cflcontrolled=false;
	}
	else
	{
		cflcontrolled=true;
	}

	sodflag=0;
	pp.query("sod_shock_tube",sodflag);

	if(sodflag > 0) 
	{
		sod_stube_flag=true;
	}
	else
	{
		sod_stube_flag=false;
	}

	readbc();
	printbcparams();

}
//===============================================================
int solver::getvarindex(std::string vname)
{
	int index;
	
	index=-1;
	for(unsigned int i=0;i<varnames.size();i++)
	{
		if(varnames[i].compare(vname) == 0) 
		{
			index=i;
			break;
		}
	}

	if(index < 0)
	{
		std::cerr<<"Variable name "<<vname<<" does not exist\n";
		exit(0);
	}
	
	return(index);
	
}
//===============================================================
void solver::makefabs()
{
	varnames.resize(0);

	//add variables
	varnames.push_back("density");
	varnames.push_back("vel_x");
	varnames.push_back("vel_y");
	varnames.push_back("vel_z");
	varnames.push_back("pressure");
	varnames.push_back("EB_volfrac");

	IntVect dom_lo(IntVect(D_DECL(0,0,0)));
	IntVect dom_hi(IntVect(D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)));
	
	Box domain(dom_lo,dom_hi);
	ba.define(domain);
	ba.maxSize(max_grid_size);

	RealBox real_box;
	for (int n=0;n<BL_SPACEDIM;n++)
	{
		real_box.setLo(n,lcorner[n]);
		real_box.setHi(n,ucorner[n]);
	}

	int coord=0;

	//Boundary conditions
	Array<int> lo_bc(BL_SPACEDIM),hi_bc(BL_SPACEDIM);
	assign_solver_bcs(lo_bc,hi_bc);

	int is_periodic[BL_SPACEDIM];
	for(int i=0;i<BL_SPACEDIM;i++)
	{
      	  is_periodic[i] = 0;
          if (lo_bc[i] == 0 && hi_bc[i] == 0) 
	  {
		is_periodic[i] = 1;
	  }
	}	

	geom.define(domain,&real_box,coord,is_periodic);
	dx = geom.CellSize();

	physbcf.resize(NCVARS);
	for(int i=0;i<NCVARS;i++)
	{
		assign_bl_bcs(varnames[i],lo_bc,hi_bc);
		BCRec bcr(&lo_bc[0],&hi_bc[0]);
  		physbcf[i].define(geom, bcr, BndryFunctBase(fillboundary));
	}

	for(int i=0;i<BL_SPACEDIM;i++)
	{
		std::cout<<lo_bc[i]<<"\t"<<hi_bc[i]<<"\n";
	}
 	
	dens   = new MultiFab(ba,1,nghost);
	dens->setVal(0.0);

	velx   = new MultiFab(ba,1,nghost);
	velx->setVal(0.0);

	vely   = new MultiFab(ba,1,nghost);
	vely->setVal(0.0);

	velz   = new MultiFab(ba,1,nghost);
	velz->setVal(0.0);

	pres   = new MultiFab(ba,1,nghost);
	pres->setVal(0.0);

	eb_volfrac   = new MultiFab(ba,1,nghost);
	eb_volfrac->setVal(1.0);

	consvars = new MultiFab(ba,NCVARS,nghost);
	consvars->setVal(0.0);

	consvars_nm1 = new MultiFab(ba,NCVARS,nghost);
	consvars_nm1->setVal(0.0);

	residual = new MultiFab(ba,NCVARS,nghost);
	residual->setVal(0.0);
}
//===============================================================
void solver::writeoutputfile(int n,double time)
{
	int nvars;

        nvars = varnames.size();
	MultiFab plotvars(ba,nvars,nghost);
	plotvars.setVal(0.0);

	MultiFab::Copy(plotvars,*dens,0,getvarindex("density"),1,nghost);
	MultiFab::Copy(plotvars,*velx,0,getvarindex("vel_x"),1,nghost);
	MultiFab::Copy(plotvars,*vely,0,getvarindex("vel_y"),1,nghost);
	MultiFab::Copy(plotvars,*velz,0,getvarindex("vel_z"),1,nghost);
	MultiFab::Copy(plotvars,*pres,0,getvarindex("pressure"),1,nghost);
	MultiFab::Copy(plotvars,*eb_volfrac,0,getvarindex("EB_volfrac"),1,nghost);

	const std::string& pltfile = BoxLib::Concatenate("plt",n,5);
	writePlotFile(pltfile, plotvars, varnames, geom, time);
}
//===============================================================
void solver::intialize_vars()
{
	int box_ind,nf,ns;
	int comp;

	if(!sod_stube_flag)
	{
		update_volfrac(geom,nghost,n_cell_x,n_cell_y,n_cell_z);
	}

	dens->setVal(vars_init[0]);
	velx->setVal(vars_init[1]);
	vely->setVal(vars_init[2]);
	velz->setVal(vars_init[3]);
	pres->setVal(vars_init[4]);

	MultiFab::Multiply(*velx,*eb_volfrac,0,0,1,nghost);
	MultiFab::Multiply(*vely,*eb_volfrac,0,0,1,nghost);
	MultiFab::Multiply(*velz,*eb_volfrac,0,0,1,nghost);

	if(sod_stube_flag)
	{
		for(MFIter mfi(*dens);mfi.isValid();++mfi)
		{
			const Box& bx = mfi.validbox();

			initialize_sod((*dens)[mfi].dataPtr(),(*velx)[mfi].dataPtr(),(*vely)[mfi].dataPtr(),
			(*velz)[mfi].dataPtr(), (*pres)[mfi].dataPtr(), bx.loVect(), bx.hiVect(),
			&nghost,geom.CellSize(),geom.ProbLo(),geom.ProbHi());
		}
	}

	update_primvar_bcs(0.0);

	update_conserved_variables();
	MultiFab::Copy(*consvars_nm1,*consvars,0,0,NCVARS,nghost);
	updateResidual();
	writeoutputfile(0,0.0);
}
//===============================================================
void solver::update_primvar_bcs(Real current_time)
{
	int comp;

	comp=1;
	update_currentcomponent(&comp);
	updatebc(dens,physbcf[comp-1],geom,current_time);

	comp=2;
	update_currentcomponent(&comp);
	updatebc(velx,physbcf[comp-1],geom,current_time);
	
	comp=3;
	update_currentcomponent(&comp);
	updatebc(vely,physbcf[comp-1],geom,current_time);
	
	comp=4;
	update_currentcomponent(&comp);
	updatebc(velz,physbcf[comp-1],geom,current_time);
	
	comp=5;
	update_currentcomponent(&comp);
	updatebc(pres,physbcf[comp-1],geom,current_time);
}
//===============================================================
void solver::updateResidual()
{
	for(MFIter mfi(*dens);mfi.isValid();++mfi)
	{
		const Box& bx = mfi.validbox();

		update_inviscidResidual((*consvars)[mfi].dataPtr(),
			(*residual)[mfi].dataPtr(),bx.loVect(), bx.hiVect(), 
			&nghost,geom.CellSize(), geom.ProbLo(), geom.ProbHi());

	}

}
//===============================================================
void solver::update_conserved_variables()
{
	for(MFIter mfi(*dens);mfi.isValid();++mfi)
	{
		const Box& bx = mfi.validbox();

   		update_conservative_vars((*dens)[mfi].dataPtr(),(*velx)[mfi].dataPtr(),
		(*vely)[mfi].dataPtr(),(*velz)[mfi].dataPtr(),(*pres)[mfi].dataPtr(),
		(*consvars)[mfi].dataPtr(),bx.loVect(),bx.hiVect(),&nghost,geom.CellSize(),
		geom.ProbLo(),geom.ProbHi());
	}

}
//===============================================================
void solver::update_primitive_variables()
{
	for(MFIter mfi(*dens);mfi.isValid();++mfi)
	{
		const Box& bx = mfi.validbox();

		update_primitive_vars((*dens)[mfi].dataPtr(),(*velx)[mfi].dataPtr(),
				(*vely)[mfi].dataPtr(),(*velz)[mfi].dataPtr(),(*pres)[mfi].dataPtr(),
				(*consvars)[mfi].dataPtr(),bx.loVect(),bx.hiVect(),&nghost,geom.CellSize(),
				geom.ProbLo(),geom.ProbHi());
	}

}
//===============================================================
double solver::find_timestep()
{
	Real tstep,tstep_global;
	int nanflag;
	Real tstep_min;


	if(cflcontrolled)
	{
		tstep_min= BIGNUM;
		nanflag  = 0;
		for(MFIter mfi(*dens);mfi.isValid();++mfi)
		{
			const Box& bx = mfi.validbox();

			findcfltimestep((*consvars)[mfi].dataPtr(),
					bx.loVect(),bx.hiVect(),&nghost,geom.CellSize(),
					geom.ProbLo(),geom.ProbHi(),&cfl,&tstep,&nanflag);
			if(nanflag)
			{
				break;
			}

			if(tstep < tstep_min)
			{
				tstep_min = tstep;
			}
		}
		if(nanflag)
		{
			std::cerr<<"\nNaN detected in finding time step\n";
			std::cerr<<"*********************************\n";
#ifdef BL_USE_MPI
			MPI_Abort(MPI_COMM_WORLD,911);
#else
			exit(0);
#endif
		}

#ifdef BL_USE_MPI
		MPI_Allreduce(&tstep_min,&tstep_global,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#else
		tstep_global = tstep_min;
#endif

	}
	else
	{
		tstep_global=delt;
	}

	return(tstep_global);
}
//===============================================================
void solver::advance_conserved_variables(Real current_time,Real &timestep,int nstep)
{
	Real timefactor,tstep;

	timefactor = 1.0;

	update_primvar_bcs(current_time);
	update_conserved_variables();
	MultiFab::Copy(*consvars_nm1,*consvars,0,0,NCVARS,nghost);
	tstep = find_timestep();
	timestep = tstep;

	if (ParallelDescriptor::IOProcessor() && (nstep%stdoutstepnum == 0) )	
	{
		std::cout<<"\n===================================================\n";
		std::cout<<"time step num:"<<nstep<<"\tcurrent time:"<<current_time<<"\n";
		std::cout<<"time step:"<<tstep<<"\n";
	}

	if(temporal_order == 1)
	{
		updateResidual();

		for(MFIter mfi(*dens);mfi.isValid();++mfi)
		{
			const Box& bx = mfi.validbox();

			advance_fwdeuler((*consvars_nm1)[mfi].dataPtr(),(*consvars)[mfi].dataPtr(),
					(*residual)[mfi].dataPtr(),(*eb_volfrac)[mfi].dataPtr(),
					bx.loVect(),bx.hiVect(),&nghost,geom.CellSize(),
					&tstep,&timefactor,geom.ProbLo(),geom.ProbHi());
		}

		update_primitive_variables();
	}
	else if(temporal_order == 4)
	{
		for(int st=0;st<NRKSTAGES;st++)
		{
			if (ParallelDescriptor::IOProcessor() && (nstep%stdoutstepnum == 0) )	
			{
				std::cout<<"RK stage:"<<st<<"\n";
			}

			update_primvar_bcs(current_time);
			update_conserved_variables();
			updateResidual();

			timefactor = RK4COEFFS[st];

			for(MFIter mfi(*dens);mfi.isValid();++mfi)
			{
				const Box& bx = mfi.validbox();

				advance_fwdeuler((*consvars_nm1)[mfi].dataPtr(),(*consvars)[mfi].dataPtr(),
						(*residual)[mfi].dataPtr(),(*eb_volfrac)[mfi].dataPtr(),
						bx.loVect(),bx.hiVect(),&nghost,geom.CellSize(),
						&tstep,&timefactor,geom.ProbLo(),geom.ProbHi());
			
			}

			update_primitive_variables();

		}	
	}
	else
	{
		if (ParallelDescriptor::IOProcessor())
		{
			std::cerr<<"temporal order "<<temporal_order<<" not implemented yet\n";

#ifdef BL_USE_MPI
			MPI_Abort(MPI_COMM_WORLD,911);
#else
			exit(0);
#endif
		}

	}
	if (ParallelDescriptor::IOProcessor() && (nstep%stdoutstepnum == 0) )
	{	
		std::cout<<"===================================================\n";
	}

}
//===============================================================
void solver::timestepping()
{
	int nstep,outputnum;
	Real current_time,timestep;
	Real output_time;

	current_time = 0.0;
	output_time = 0.0;
	outputnum = 1;
	nstep = 0;

	while(current_time < tfinal)
	{
		if(output_time >= dt_out)
		{
			writeoutputfile(outputnum++,current_time);
			output_time=0.0;
		}

		advance_conserved_variables(current_time,timestep,nstep);

		current_time += timestep;	
		output_time  += timestep;
		nstep++;
	}
	writeoutputfile(outputnum,current_time);
}
//===============================================================
