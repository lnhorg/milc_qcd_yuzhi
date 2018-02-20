/*************** update_u.c ************************************/
/* MIMD version 7 */

/* update the link matrices					*
 *  								*
 *  Go to eight order in the exponential of the momentum		*
 *  matrices, since higher order integrators use large step	*
 *  Evaluation is done as:					*
 *	exp(H) * U = ( 1 + H + H^2/2 + H^3/3 ...)*U		*
 *	= U + H*( U + (H/2)*( U + (H/3)*( ... )))		*
 *								*
 */

#include "ks_imp_includes.h"	/* definitions files and prototypes */

#ifdef USE_GF_GPU

#include "../include/generic_quda.h"

/* QUDA version */

void update_u(Real eps){

  int i,dir;
  site *s;
  int j;

#ifdef FN
  invalidate_fermion_links(fn_links);
#endif

  initialize_quda();
  Real *momentum = (Real*)malloc(sites_on_node*4*sizeof(anti_hermitmat));
  Real *gauge = (Real*)malloc(sites_on_node*4*sizeof(su3_matrix));
  // Populate gauge and momentum fields
  FORALLSITES(i,s){
    for(dir=XUP; dir<=TUP; ++dir){
      for(j=0; j<18; ++j){
        gauge[(4*i + dir)*18 + j] = ((Real*)(&(s->link[dir])))[j];
      }
      for(j=0; j<10; ++j){
        momentum[(4*i + dir)*10 + j] = ((Real*)(&(s->mom[dir])))[j];
      }
    } // dir
  }

  qudaUpdateU(PRECISION, eps, momentum, gauge);

  // Copy updated gauge field back to site structure
  FORALLSITES(i,s){
    for(dir=XUP; dir<=TUP; ++dir){
      for(j=0; j<18; ++j){
        ((Real*)(&(s->link[dir])))[j] = gauge[(4*i + dir)*18 + j];
      }
    }
  }
  free(momentum);
  free(gauge);
  return;
}

#else

/* CPU version */

void update_u( Real eps ){
#ifndef U1_ONLY
  register int i,dir;
  register site *s;
  su3_matrix *link,temp1,temp2,htemp;
  register Real t2,t3,t4,t5,t6,t7,t8;
  /**TEMP**
    Real gf_x,gf_av,gf_max;
    int gf_i,gf_j;
   **END TEMP **/

  /**double dtime,dtime2,dclock();**/
  /**dtime = -dclock();**/

  /* Take divisions out of site loop (can't be done by compiler) */
  t2 = eps/2.0;
  t3 = eps/3.0;
  t4 = eps/4.0;
  t5 = eps/5.0;
  t6 = eps/6.0;
  t7 = eps/7.0;
  t8 = eps/8.0;

  /** TEMP **
    gf_av=gf_max=0.0;
   **END TEMP**/
#ifdef FN
  invalidate_fermion_links(fn_links);
#endif

  FORALLSITES(i,s){
    for(dir=XUP; dir <=TUP; dir++){
      uncompress_anti_hermitian( &(s->mom[dir]) , &htemp );
      link = &(s->link[dir]);
      mult_su3_nn(&htemp,link,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,t8,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,t7,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,t6,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,t5,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,t4,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,t3,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,t2,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,eps    ,&temp2); 
      su3mat_copy(&temp2,link);
    }
  }

  /**dtime += dclock();
    node0_printf("LINK_UPDATE: time = %e  mflops = %e\n",
    dtime, (double)(5616.0*volume/(1.0e6*dtime*numnodes())) );**/
#endif /* U1_ONLY */
#ifdef HAVE_U1
  update_u_u1(eps);
#endif
} /* update_u */

#endif

/* Added by YL on 05/27/2016 */
void update_u_u1(Real eps) {
#ifdef HAVE_U1
  register int i, dir;
  register site *s;

  eps *= (-1.0 / 3.0 / pseudo_charges[0]);
  FORALLSITES(i, s) {
#ifdef SCHROED_FUN
    for (dir = XUP; dir <= TUP; dir++) if (dir == TUP || s->t > 0) {
#else
    for (dir = XUP; dir <= TUP; dir++) {
#endif
        u1_A[4 * i + dir] += eps * (s->mom_u1[dir]); /* dA/d\tau = mom_u1 */
#ifdef U1_DEBUG
        if (dir == XUP && i == 0) {
          node0_printf("update_u.c s->mom_u1[XUP] eps, s->mom_u1[dir], "
                       "u1_A[4*i+dir] %e %e %e\n",
                       eps, s->mom_u1[dir], u1_A[4 * i + dir]);
        }
#endif
      }
  }

#else /* HAVE_U1 */
  printf("function update_u_u1 should not be called in a QCD code\n");
  terminate(-1);
#endif

} /* update_u_u1 */
/* YL end */
