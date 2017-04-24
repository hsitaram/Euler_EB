#ifndef IB_FORTRAN_FUNC_H
#define IB_FORTRAN_FUNC_H

#include <BLFort.H>

extern "C"
{
  void spherefindindices(const int* lo, const int* hi,const Real* dx, const Real* prob_lo, 
		 const Real* prob_hi,double* rad,double* pos,int 
		*fldindices,int *sldindices,int* nindmax,int *nfldindices,int *nsldindices);

  void find_volfrac_in_box_forsphere(Real *vfrac,const int *lo,const int *hi,int *ng,const Real *dx,
	const Real *problo,const Real *probhi,Real *rad,Real *pos);

  void find_volfrac_in_box(Real *vfrac,const int *glo,const int *ghi,const int *lo,const int *hi,
	int *ng,const Real *dx,const Real *problo,const Real *probhi,Real *surfcoord,int *conn,int *nnodes,int *ntri);

   void read_triangulation_file(Real *surfcoord,int *conn,int *nnodes,int *ntri);
}

#endif
