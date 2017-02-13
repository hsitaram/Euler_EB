#ifndef C_CONSTANTS_H
#define C_CONSTANTS_H

#define XDIR 0
#define YDIR 1
#define ZDIR 2
#define DIM3 3
#define BC_PERIODIC 0 
#define BC_SUPINLET 1 
#define BC_SUBINLET 2 
#define BC_SUPOUTLET 3
#define BC_SUBOUTLET 4
#define BC_WALL 5
#define NB 6
#define NCVARS 5
#define LFT 0
#define RGT 1
#define BOT 2
#define TOP 3
#define BCK 4
#define FRT 5

#define RHO_INDEX   1
#define RHO_U_INDEX 2
#define RHO_V_INDEX 3
#define RHO_W_INDEX 4
#define RHO_E_INDEX 5

#define BIGNUM 3e8
const int NRKSTAGES=5;
const Real RK4COEFFS[]={0.0533,0.1263,0.2375,0.4414,1.0};

#endif
