#ifndef _FN_LINKS_QOP_H
#define _FN_LINKS_QOP_H
/************************ fn_links_qop.h ***************************/
/* MILC version 7 */

#include "../include/complex.h"
#include "../include/precision.h"
#include "../include/su3.h"
#include "../include/link_phase_info.h"
#include "../include/ks_action_coeffs_qop.h"
#include "../include/imp_ferm_links.h"
#include <qop.h>

/* The fn_links_qop_t "class" */

typedef struct {
  int preserve;
  link_phase_info_t *phase;
  int al_F_allocated;
  QOP_F3_FermionLinksAsqtad *al_F;
  QOP_D3_FermionLinksAsqtad *al_D;
  /* Derived, read-only values extracted from al_F or al_D if needed: */
  su3_matrix *fat, *lng;
} fn_links_qop_t;

#if FERM_ACTION != HISQ  /* To avoid conflicts with hisq_links_qop_t */

void load_fn_links_qop(fn_links_qop_t *fn, QOP_asqtad_coeffs_t *ac, 
		       int precision, su3_matrix *links, int want_back);

#define destroy_fn_links destroy_fn_links_qop

#endif

QOP_F3_FermionLinksAsqtad *get_F_asqtad_links(fn_links_qop_t *fn);
QOP_D3_FermionLinksAsqtad *get_D_asqtad_links(fn_links_qop_t *fn);

void destroy_fn_links_qop(fn_links_qop_t *fn);
fn_links_qop_t *create_fn_links_qop(void);
void free_fn_links_qop(fn_links_qop_t *fn);

#endif /* _FN_LINKS_QOP_H */
