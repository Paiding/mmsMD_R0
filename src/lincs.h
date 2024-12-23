#ifndef _LINCS_H
#define _LINCS_H
#include <includes.h>
#include <top.h>
#include <virial.h>
#define WANGLE 30
#define WFAC (0.5236) // = DEG2RAD * wangle

struct t_lincs_data
{
    int cons_num = 0;

    array<int> A_list;
    array<int> B_list;
    array<float> len;

    array<float> dir_x;
    array<float> dir_y;
    array<float> dir_z;

    array<int> nbr_cons;
    array<int> nbr_cons_loc;

    array<float> Sdiag;
    array<float> coef;

    array<float> blcc;
    array<float> rhs1;
    array<float> rhs2;
    array<float> sol;

    array<float> lambda;

};

void init_lincs_data(array<t_molecule> & molelist,
                    t_lincs_data & lincs_data,
                    t_atom_statics & atom_statics,
                    t_constr_para & cons_para
                    );

void calc_direct(t_atom_kinematics & old_k,
                 t_lincs_data & lincs_data
                );      

void b4expand_0(t_lincs_data & ld,
                t_atom_kinematics & ak
                );

void lincs_mat_exp_0(t_lincs_data & ld,int flag);

void lincs_mat_exp_1(t_lincs_data & ld);

void b4expand_1(t_lincs_data & ld,
                t_atom_kinematics & ak
                );

void lincs_mat_exp_2(t_lincs_data & ld);

void lincs_solve(t_lincs_data & ld,
                 t_atom_kinematics & ak,
                 t_atom_statics & as,
                 float invdt
                );

void lincs_virial(t_virial & virial,
                  t_lincs_data & ld,
                  float invdt
                 );

void do_lincs(t_lincs_data & lincs_data,
              t_atom_statics & atom_statics,
              t_atom_kinematics & atom_kinemetics,
              t_atom_kinematics & old_k,
              t_virial & virial,
              float invdt
            );

#endif