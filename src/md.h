#ifndef _MD_H
#define _MD_H
#include "includes.h"
#include "top.h"
#include "vdwcf.h"
#include "update.h"
#include "nblist.h"
#include "file.h"
#include "pretreat.h"
#include "bond.h"
#include "angle.h"
#include "settle.h"
#include "dihedral.h"
#include "lj14.h"
#include "lincs.h"
#include "p_couple.h"
#include "vcm.h"
#include "tools.h"

int init_all(t_atom_kinematics & atom_kinematics,
              t_atom_kinematics & old_k,
              t_atom_statics & atom_statics,
              t_tcstat & tcstat, t_tcpara & tcpara,
              t_pcstat & pcstat, t_pcpara & pcpara,
              t_force & force,
              t_virial & virial,
              t_nstep & step,
              t_radius & radius,
              array<t_molecule> & mole_array,
              t_mole_map & mole_map,
              t_dsf_para & dsf_para,
              t_vcm & vcm
              );

void one_step(int step, t_nstep & nstep, float dt,
              float *box, array<float> & lj_para,
              t_dsf_para & dsf_para,
              t_shift_switch_para & sspara,
              t_ub_para & ubpara,t_g96a_para & g96a_para, t_pdih_para & pdpara, t_ipdih_para & idpara, t_lj14_para & lj14_para,
              int atom_size,t_radius & radius,
              t_nblist & nblist,
              t_atom_kinematics & atom_kinematics,
              t_atom_kinematics & old_k,
              t_atom_statics & atom_statics,
              t_tcstat & tcstat, t_tcpara & tcpara,
              t_pcstat & pcstat, t_pcpara & pcpara,
              t_lincs_data & lincs_data,
              t_force & force,
              t_virial & virial,
              int tctype,
              t_excl_list & excl_list,
              t_mole_map & mole_map,
              array<t_molecule> & molelist,
              t_vcm & vcm,
              t_group & tcgrp
              );

void kinematics_copy(t_atom_kinematics & src,t_atom_kinematics & dest);
void force_log(t_force & force);
void fclear(t_force & force, t_virial & virial);
void csv_print_data_2(char* filename, t_atom_kinematics & atom_kinematics,
                      t_atom_statics & atom_statics,
                      t_atom_map & atom_map
                     );
void init_box(float *box);

#endif