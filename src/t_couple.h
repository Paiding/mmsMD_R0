#ifndef _T_COUPLE_H
#define _T_COUPLE_H

#include "includes.h"
#include "top.h"
#include "file.h"

typedef struct {
  int     N;							 /* T-Coupl groups */
  float   *nrdf;	         /* Nr of degrees of freedom in a group */
  float   *ref_t;	         /* Coupling temperature per group */
  float   *tau_t;	         /* Tau coupling time */
} t_tcpara; // it is only needed in CPU

typedef struct {
	int      N;              /* number of T-coupling groups */
  array<float> Th;            /* Dim: N. Temperature at half step */
  array<float> T;             /* Dim: N. Temperature at full step */
  array<float> ekinh_old;     /* Dim: DIM*DIM*N, Kinetic energy at old half step */
	/* This array is needed in GPU */ 
  array<float> ekinh;        /* Dim: DIM*DIM*N, Kinetic energy at half step.*/ 
  array<float> ekin;          /* Dim: DIM*DIM*N, Kinetic energy at full step */
	/* This array is needed in GPU */ 
  float    *lambda;        /* Dim: N, Berendsen coupling lambda.*/
	/* This array is needed in GPU */ 
  //float    *nh_xi;         /* for Nose-Hoover tcoupl (ngtc) */

	float    ekin_tot[3]; // dim: 9        cross not used , 3 dims only 
	//float    *T_tot;    // dim: 1

  float T_sum = 0;
  float nrdf_sum = 0;
// tmp **   double   *therm_integral; /* for N-H/V-rescale tcoupl (ngtc) */
} t_tcstat;

void ekin_sum(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics,
              t_tcstat & tcstat,t_tcpara & tcpara,t_group & tcgrp);
float Berendsen_lambda(t_tcpara & tcpara, t_tcstat & tcstat,
                        float dt, int index);

void gen_vel(t_atom_kinematics & atom_kinematics,
             t_atom_statics & atom_statics,
             array<t_molecule> & molelist,
             t_tcpara & tcpara
            );

void init_tc_group(t_group & tcgrp, int grp_num);
void init_tc_nogroup(t_group & tcgrp,int atom_size);












#endif