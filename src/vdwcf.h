#pragma once//avoid PCH warning when using #ifndef code in the last .h file

//#ifndef _VDWCF_H
//#define _VDWCF_H

#include "includes.h"
#include "top.h"
#include "nblist.h"

#include "virial.h"

struct t_dsf_para
{
    float alpha;
    float Rc_force;
};

struct t_shift_switch_para
{
    float LJ_12_A = 0.0;
    float LJ_12_B = 0.0;
    float LJ_12_C = 0.0;
    float LJ_12_A_3 = 0.0;
    float LJ_12_B_4 = 0.0;

    float LJ_6_A = 0.0;
    float LJ_6_B = 0.0;
    float LJ_6_C = 0.0;
    float LJ_6_A_3 = 0.0;
    float LJ_6_B_4 = 0.0;
};

float dump_shift(float alpha, float r);

void init_shift_switch_para(t_radius & radius,
                            t_shift_switch_para & sspara
                            );

void vdwcf_dsf_ss_img(const array<float> & lj_para,
                      t_dsf_para & dsf_para,
                      t_shift_switch_para & sspara,
                      const t_atom_statics & atom_statics,
                      const t_atom_kinematics & atom_kinematics,
                      const t_nblist & nblist,
                      t_force & force, 
                      t_virial & virial,
                      float * box,
                      int atom_size
                      );

void vdwcf_dsf_img(const array<float> & lj_para,
               t_dsf_para & dsf_para,
               const t_atom_statics & atom_statics,
               const t_atom_kinematics & atom_kinematics,
               const t_nblist & nblist,
               t_force & force, 
               t_virial & virial,
               float * box
              );





//#endif
