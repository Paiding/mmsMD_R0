#define __P_COUPLE__
/* -------------------------------------------------------- *\
 *   DEMms is a software package for multiscale simulation  *
 *   of discrete element method.                            *
 * -------------------------------------------------------- *

License: This file is part of DEMms.
Maintainer: Ji Xu

 * -------------------------------------------------------- */
 

 
#define IN_LINUX
// tmp ** #define IN_VC6_0
// tmp ** #define IN_VS2008

#define DEBUG_ME

// tmp ** #define USE_DOUBLE

// tmp ** #define SOCKET

// mpi
#define USE_MPI
#define USE_INTEL_MPI

// GPU
//#define USE_GPU
#define USE_STREAM
#define K80
// tmp ** #define C2050
#define MAX_GPUS 24 // maximum GPUs in one node

#ifndef USE_DOUBLE
#define USE_TEXTURE
#endif
 
// thermal
#define USE_THERMAL

// cfd-dem
#define TWO_PHASE_COUPLE
 
// shared memory
#ifdef TWO_PHASE_COUPLE
#define USE_SYSTEM_SHARED_MEMORY
#endif
 
// reaction
#define USE_REACTION
#define USE_REACTION_MTO // MTO reaction, 7 lumped / 2 circle
// tmp ** #define USE_REACTION_LINEAR // A -> B, A -> B air
// tmp ** #define USE_REACTION_MTO_COKE_BURN

// vapor
// tmp ** #define USE_VAPOR

// cfd-dem thermal
#define USE_TWO_PHASE_THERMO

// tmp ** #define RUNTIME_MODIFY

// tmp ** #define REMOTE_COKE

// tmp ** #define USE_TIMER

