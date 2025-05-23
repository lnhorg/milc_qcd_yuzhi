/****************** ks_imp_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_ksprop.h"
#include "../include/io_lat.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/params_rhmc.h"
#include "../include/imp_ferm_links.h"
#include "quark_action.h"
#ifdef HAVE_QUDA
#include "../include/generic_quda.h"
#endif
#include "../include/dirs.h"
#include "../include/dirs.h"
#ifdef HAVE_QIO
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#endif

#ifdef FN
#define dslash_site dslash_fn_site
#define dslash_field dslash_fn_field
#endif
#ifdef EO
#define dslash_site dslash_eo_site
#define dslash_field dslash_eo_field
#endif

#ifdef PRTIME
#define STARTTIME dtime = -dclock();
#define ENDTIME(string) dtime += dclock(); node0_printf("Aggregate time to %s %e\n",(string),dtime);  fflush(stdout);
#else
#define STARTTIME
#define ENDTIME(string)
#endif


/* prototypes for functions in high level code */
int setup();
int readin(int prompt);

int ask_color_vector( int prompt, int *flag, char *filename );
int ask_color_matrix( int prompt, int *flag, char *filename );
void check_link_fattening( char *lngansfile, int lngansflag, char *fatansfile, int fatansflag );
void check_fermion_force( char srcfile[MAX_MASS][MAXFILENAME], int srcflag,
			  char *ansfile, int ansflag, int nmass, ks_param *ksp);
void check_ks_invert( char *srcfile, int srcflag, 
		      char ansfile[MAX_MASS][MAXFILENAME],
		      int ansflag[MAX_MASS],
		      int nmass, ks_param ksp[], 
		      quark_invert_control qic[]);
void check_invert2( su3_vector *src, su3_vector *dest, 
		    Real mass, Real tol, int parity,
		    imp_ferm_links_t *fn);
void check_reunitarization_derivative( char *ansfilein, int ansflagin,
				       char *ansfileout, int ansflagout );
char *create_QCDML();
void free_QCDML(char *qcdml);

