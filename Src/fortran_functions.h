#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <BLFort.H>

extern "C"
{
  void print_var(Real* data, const int* lo, const int* hi, const int* ng, 
		const Real* dx, const Real* prob_lo, const Real* prob_hi);

  void init_blanking_field(Real *data,const int *lo,const int *hi,const int *ng,
	const Real *dx,const Real *prob_lo,const Real *prob_hi,
	int *findices,int *sindices,int *nfindices,int *nsindices);

  void fillboundary(BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const Real* dx, const Real* glo, 
		  const Real* time, const int* bc);

   void update_currentcomponent(int *val);
   void set_dirichlet_bcs(Real *d,Real *vx,Real *vy,Real *vz,Real *p);

   void update_inviscidResidual(Real *consvar,Real *residual,const int *lo,const int *hi,const int *ng,
	const Real *dx,const Real *prob_lo,const Real *prob_hi);

   void update_conservative_vars(Real *dens,Real *velx,Real *vely,Real *velz,Real *pres,Real *consvar,
	const int *lo,const int *hi,const int *ng,const Real *dx,const Real *prob_lo,const Real *prob_hi);

   void update_primitive_vars(Real *dens,Real *velx,Real *vely,Real *velz,Real *pres,Real *consvar,
	const int *lo,const int *hi,const int *ng,const Real *dx,const Real *prob_lo,const Real *prob_hi);

   void findcfltimestep(Real *consvar,const int *lo,const int *hi,const int *ng,const Real *dx,
	const Real *prob_lo,const Real *prob_hi,Real *cfl,Real *tstep,int *nanflag);

   void advance_fwdeuler(Real *consvars_nm1,Real *consvars,Real *residual,Real *vfraction,const int *lo,const int *hi,
	const int *ng,const Real *dx,Real *tstep,Real *tfactor,const Real *problo,const Real *probhi);

  void initialize_sod(Real* dens,Real *velx,Real *vely,Real *velz,Real *pres,const int* lo, 
		const int* hi, const int* ng, const Real* dx, const Real* prob_lo, const Real* prob_hi);
}
#endif
