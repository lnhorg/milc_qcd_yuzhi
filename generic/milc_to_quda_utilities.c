/******************* milc_to_quda_utilities.c ************************/
/* For the QUDA/GPU interface */

//#include <cuda.h>  // DEBUG-JNS
//#include <cuda_runtime.h> // DEBUG-JNS
#include "generic_includes.h"
#include "../include/generic_quda.h"
#include <string.h>

static int is_quda_initialized = 0;


int initialize_quda(void){

  QudaInitArgs_t init_args;

#if defined(SET_QUDA_SILENT)
  init_args.verbosity = QUDA_SILENT;
#elif defined(SET_QUDA_VERBOSE)
  init_args.verbosity = QUDA_VERBOSE;
#elif defined(SET_QUDA_DEBUG_VERBOSE)
  init_args.verbosity = QUDA_DEBUG_VERBOSE;
#else
  init_args.verbosity = QUDA_SUMMARIZE; /* default */
#endif

  const int dim[4] = {nx, ny, nz, nt};
  int status = 0;

  if(is_quda_initialized)return status;

  init_args.layout.device = 0; 								// only valid for single-gpu build
  init_args.layout.latsize = dim;
  init_args.layout.machsize = get_logical_dimensions();

  /* Tell QUDA which communicator we are using, in case we have split it */
  qudaSetMPICommHandle(mycomm());
  qudaInit(init_args);

  //cudaDeviceSetLimit(cudaLimitPrintfFifoSize,128*1024*1024); // DEBUG-JNS


  if(status == 0)
    is_quda_initialized = 1;

  return status;

} /* milc_to_quda_utilities */

void finalize_quda(void){
#ifdef USE_CG_GPU
  qudaCleanUpDeflationSpace();
#ifdef MULTIGRID
  mat_invert_mg_cleanup();
#endif
#endif
  qudaFinalize();
}
