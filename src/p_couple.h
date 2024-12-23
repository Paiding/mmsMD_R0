#ifndef _P_COUPLE_H_
#define _P_COUPLE_H_

#include "includes.h"
#include "virial.h"
#include "t_couple.h"

/* 
 * This file implements temperature and pressure coupling algorithms:
 * For now only the Weak coupling and the modified weak coupling.
 * Furthermore computation of pressure and temperature is done here.
 */

typedef struct {
	int do_pc = 0;
	float tau_p; 
	float dt;
	float p_ref[9] = {1.0,0.0,0.0,
                    0.0,1.0,0.0,
                    0.0,0.0,1.0}; // dim: 9
	float compress[9] = {4.5e-5,0.0,0.0,
                       0.0,4.5e-5,0.0,
                       0.0,0.0,4.5e-5}; // dim: 9
} t_pcpara; // only cpu

//extern double pcpara_init(t_pcpara *d_pcpara, t_pcpara *h_pcpara, t_para *h_para);
//extern void pcpara_free(t_pcpara *d_pcpara, t_pcpara *h_pcpara);

typedef struct {
  float	pres_scalar = 1.0f;  /* dim: 1*/

  float	pres_prev[3];	/* dim: 9, the pressure used for pressure coupling*/
  float	pres[3];	/* dim: 9, the pressure used for pressure coupling*/
// tmp **   float	*ekin;		/* dim: 9, kinetic energy*/
// tmp **   float	*vir;		/* dim: 9, virial of the whole system */
  float	mu[9];		/* dim:9,  the scaling matrix of pressure couping */
} t_pcstat;
void calc_vir(t_virial & virial,
              t_force & force,
              t_atom_kinematics & atom_kinematics
             );
/* The pressure calculating need the kinetic energy and virial tensor
 * and also the long-range electrostatic correction if not with PME or Ewald.
 */
void calc_pres(t_tcstat & tcstat,
               t_pcstat & pcstat,
               t_virial & virial,
               float *box
              );

/* we get the mu matrix in this function */
void Berendsen_pc_mu(t_pcstat & pcstat,
                     t_pcpara & pcpara
                    );
/* use the mu matrix to scale the particles' position and box. Also copy the box to the last box */
void Berendsen_pc_rx_scale(t_pcstat & pcstat,
                           t_atom_kinematics & atom_kinematics
                          );
void Berendsen_pc_box_scale(t_pcstat & pcstat,
                           float * box
                           );

#endif