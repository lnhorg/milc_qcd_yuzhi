#ifndef _FLAVOR_OPS_H
#define _FLAVOR_OPS_H

#include "../include/gammatypes.h"

void mult_pion5_field( int r0[], const su3_vector *const src, su3_vector *dest );
void mult_pion05_field( int r0[], const su3_vector *const src, su3_vector *dest );
void mult_pioni5_field( int fdir, int r0[], const su3_vector *const src, su3_vector *dest,
			const su3_matrix *const links,int *refresh_links );
void mult_pionij_field( int fdir, int r0[], const su3_vector *const src, su3_vector *dest,
			const su3_matrix *const linkss,	int *refresh_links );
void mult_pioni_field( int fdir, int r0[], const su3_vector *const src, su3_vector *dest,
		       const su3_matrix *const linkss, int *refresh_links );
void mult_pioni0_field( int fdir, int r0[], const su3_vector *const src, su3_vector *dest,
			const su3_matrix *const linkss,	int *refresh_links );
void mult_pions_field(int r0[], const su3_vector *const src, su3_vector *dest,
		      const su3_matrix *const linkss, int *refresh_links );
void mult_pion0_field(int r0[], const su3_vector *const src, su3_vector *dest,
		      const su3_matrix *const linkss, int *refresh_links );
void mult_rhoi_field( int pdir,  int r0[], const su3_vector *const src, su3_vector *dest );
void mult_rhoi0_field( int pdir,  int r0[], const su3_vector *const src, su3_vector *dest );
void mult_rhos_field( int fdir,  int r0[], const su3_vector *const src, su3_vector *dest,
		      const su3_matrix *const linkss, int *refresh_links );
void mult_rho0_field( int fdir,  int r0[], const su3_vector *const src, su3_vector *dest,
		      const su3_matrix *const linkss, int *refresh_links );
void spin_taste_op(int index, int r0[], su3_vector *dest, const su3_vector *const src);

int spin_taste_index(char *label);
const char *spin_taste_label(int index);
int is_rhosfn_index(int index);
int is_rhosffn_index(int index);
int is_rhosbfn_index(int index);
int is_rhosape_index(int index);
int is_rhosfape_index(int index);
int is_rhosbape_index(int index);
int is_fn_index(int index);
int forward_index(int index);
int backward_index(int index);
int is_fn_index(int index);
  
#endif /* _FLAVOR_OPS_H */
