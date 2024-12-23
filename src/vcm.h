#ifndef _VCM_H
#define _VCM_H
#include "top.h"

struct t_vcm
{
    int n = 0;
    float m = 0.0;
    //momentums
    float group_p_x = 0.0;
    float group_p_y = 0.0;
    float group_p_z = 0.0;

    float vcm_x = 0.0;
    float vcm_y = 0.0;
    float vcm_z = 0.0;

    float l_x = 0.0;
    float l_y = 0.0;
    float l_z = 0.0;

    float mr2[6] = {0};

    float mx = 0.0;
    float my = 0.0;
    float mz = 0.0;

    float w_x;
    float w_y;
    float w_z;

};

void calc_mass_grp(t_atom_statics & atom_statics,
                   t_vcm & vcm
                  );
void calc_vcm_grp(t_atom_statics & atom_statics,
                  t_atom_kinematics & atom_kinematics,
                  t_vcm & vcm
                 );
void vcm_v_calc(t_vcm & vcm);

void calc_vcm_ang_grp(t_atom_statics & atom_statics,
                  t_atom_kinematics & atom_kinematics,
                  t_vcm & vcm
                 );
void vcm_w_calc(t_vcm & vcm);

void do_stopcm_grp(t_atom_kinematics & atom_kinematics,
                   t_vcm & vcm,
                   int is_angular
                  );




#endif