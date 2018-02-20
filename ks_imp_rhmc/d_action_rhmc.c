/*************** d_action_rhmc.c ****************************************/
/* MIMD version 7 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.  */

#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "../include/fermion_links.h"
static Real ahmat_mag_sq(anti_hermitmat *pt);

/*DEBUG*/
double old_g, old_h, old_f, old_a;
#ifdef HAVE_U1
double old_g_u1, old_h_u1, old_a_u1;
#endif
/*ENDDEBUG*/

double d_action_rhmc( su3_vector **multi_x, su3_vector *sumvec){
  double ssplaq,stplaq,g_action,h_action,f_action;
#ifdef HAVE_U1
  double g_action_u1, h_action_u1;
  Real s2Fmunu_u1, t2Fmunu_u1;
#endif

  d_plaquette(&ssplaq,&stplaq);
  ssplaq *= -1.0; stplaq *= -1.0;
  g_action = -beta*volume*(ssplaq+stplaq);
  node0_printf("PLAQUETTE ACTION: %e\n",g_action);

  rephase(OFF);
  g_action = (beta/3.0)*imp_gauge_action();
  rephase(ON);
  h_action = hmom_action();
  f_action = fermion_action(multi_x,sumvec);

  node0_printf("ACTION: g,h,f = %.14e  %.14e  %.14e  %.14e\n",
	       g_action, h_action, f_action, g_action+h_action+f_action );

#ifdef HAVE_U1
  u1Fmunu(&s2Fmunu_u1, &t2Fmunu_u1);
  /* The U1 gauge action is 1 / 2 \sum_{n, \mu < \nu} F_{\mu\nu}F^{\mu\nu}. */
  g_action_u1 = (beta_u1 / 2.0) * (s2Fmunu_u1 + t2Fmunu_u1);
  h_action_u1 = hmom_action_u1();
  node0_printf("ACTION U1: g_u1, h_u1, g_u1+h_u1, "
               "s2_u1, t2_u1 = %.14e  %.14e  %.14e %.14e %.14e\n",
               g_action_u1,
               h_action_u1,
               g_action_u1 + h_action_u1,
               s2Fmunu_u1,
               t2Fmunu_u1);
  node0_printf("ACTION SU3+U1: f = %.14e\n",
               g_action + h_action + g_action_u1 + h_action_u1 + f_action );
#endif


  /*DEBUG*/
  node0_printf("DG = %e, DH = %e, DF = %e, D = %e\n",
	       g_action-old_g, h_action-old_h, f_action-old_f,
	       g_action+h_action+f_action-old_a);
#ifdef HAVE_U1
  node0_printf("DG_U1 = %e, DH_U1 = %e, D_U1 = %e\n",
               g_action_u1 - old_g_u1,
               h_action_u1 - old_h_u1,
               g_action_u1 + h_action_u1 - old_a_u1);

  node0_printf("DG_SU3_U1 = %e, DH_SU3_U1 = %e, "
               "DF_SU3_U1 = %e D_SU3_U1 = %e\n",
               g_action + g_action_u1 - old_g - old_g_u1,
               h_action + h_action_u1 - old_h - old_h_u1,
               f_action - old_f,
               g_action + g_action_u1 + h_action + h_action_u1 + f_action - old_a - old_a_u1);
  old_g_u1 = g_action_u1;
  old_h_u1 = h_action_u1;
  old_a_u1 = g_action_u1 + h_action_u1;
#endif
  old_g=g_action; old_h=h_action; old_f=f_action;
  old_a=g_action+h_action+f_action;
  /*ENDDEBUG*/

#ifdef HAVE_U1
  return (g_action + h_action + g_action_u1 + h_action_u1 + f_action);
#else
  return(g_action+h_action+f_action);
#endif
}

/* fermion contribution to the action */
double fermion_action( su3_vector **multi_x, su3_vector *sumvec) {
  register int i;
  register site *s;
  Real final_rsq;
  double sum;
  int iphi, inaik, jphi, n;
  imp_ferm_links_t **fn;
  sum=0.0;
  iphi=0;
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  n = fermion_links_get_n_naiks(fn_links);
#else
  n = 1;
#endif
  for( inaik=0; inaik<n; inaik++ ) {
    for( jphi=0; jphi<n_pseudo_naik[inaik]; jphi++ ) {
#ifdef HAVE_U1
      /* can be improved by checking whether the charge is changed */
      current_charge_u1 = 1.0 * pseudo_charges[iphi];
      u1phase_on(current_charge_u1, u1_A);
      invalidate_fermion_links(fn_links);
#endif
      restore_fermion_links_from_site(fn_links, prec_fa[iphi]);
      fn = get_fm_links(fn_links);
      node0_printf("fnlinks %f\n", fn[inaik]->fat[0].e[0][0].real);
      ks_ratinv( F_OFFSET(phi[iphi]), multi_x, rparam[iphi].FA.pole,
		 rparam[iphi].FA.order, niter_fa[iphi], rsqmin_fa[iphi],
		 prec_fa[iphi], EVEN, &final_rsq, fn[inaik],
		 inaik, rparam[iphi].naik_term_epsilon );
      ks_rateval( sumvec, F_OFFSET(phi[iphi]), multi_x,
  		rparam[iphi].FA.res, rparam[iphi].FA.order, EVEN );
      FOREVENSITES(i,s){ /* phi is defined on even sites only */
        sum += magsq_su3vec( &(sumvec[i]) );
      }
#ifdef HAVE_U1
      /* Unapply the U(1) field phases */
      u1phase_off();
      invalidate_fermion_links(fn_links);
#endif
      iphi++;
    }
  }

  g_doublesum( &sum );
  return(sum);
}

/* gauge momentum contribution to the action */
double hmom_action(void) {
  register int i,dir;
  register site *s;
  double sum;

  sum=0.0;
  FORALLSITES(i,s){
    for(dir=XUP;dir<=TUP;dir++){
      sum += (double)ahmat_mag_sq( &(s->mom[dir]) ) - 4.0;
      /* subtract 1/2 per d.o.f. to help numerical acc. in sum */
    }
  }
  g_doublesum( &sum );
  return(sum);
}

/* magnitude squared of an antihermition matrix */
Real ahmat_mag_sq(anti_hermitmat *pt){
  register Real x,sum;
  x = pt->m00im; sum  = 0.5*x*x;
  x = pt->m11im; sum += 0.5*x*x;
  x = pt->m22im; sum += 0.5*x*x;
  x = pt->m01.real; sum += x*x;
  x = pt->m01.imag; sum += x*x;
  x = pt->m02.real; sum += x*x;
  x = pt->m02.imag; sum += x*x;
  x = pt->m12.real; sum += x*x;
  x = pt->m12.imag; sum += x*x;
  return(sum);
}

/* gauge momentum contribution to the U(1) action */
#ifdef HAVE_U1
double hmom_action_u1(void) {
  register int i, dir;
  register site *s;
  double sum;

  sum = 0.0;
  FORALLSITES(i, s) {
    for (dir = XUP; dir <= TUP; dir++) {
      sum += 0.5 * (double)( (s->mom_u1[dir]) * (s->mom_u1[dir]));
      /* mom_u1 initially generated in ../generic/ranmom.c */
    }
  }
  g_doublesum( &sum );
  return (sum);
}
#endif
