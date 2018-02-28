#ifndef _GENERIC_U1
#define _GENERIC_U1

/************************ generic_ks.h **********************************
*									*
*  Macros and declarations for generic_u1 routines                      *
*  This header is for codes that call generic_u1 routines               *
*  MIMD version 7 							*
*									*
*/
/* gauge_stuff_u1.c */
void g_measure_nc_u1(void);

/* u1link.c */
Real *create_u1_A_field(void);
void destroy_u1_A_field(Real *A);
void u1phase_on(Real charge, Real *A);
void u1phase_off(void);
void construct_u1_links_from_A(Real charge, Real *A, complex *clink);

/* u1plaq.c */
complex u1ploop(void);
void u1plaq(Real *ssplq,Real *stplq);
void u1Fmunu(Real *ssFmunu, Real *stFmunu);

/* gauge_force_u1_cpu.c */
void gauge_force_u1_cpu( Real eps);
#ifdef USE_GF_GPU
#define gauge_force_u1 gauge_force_u1_cpu
#else
#define gauge_force_u1 gauge_force_u1_cpu
#endif

#endif /* _GENERIC_U1 */
