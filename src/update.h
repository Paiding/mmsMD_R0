#ifndef _UPDATE_H
#define _UPDATE_H
#include "top.h"
#include "t_couple.h"

void vv_update_1(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics, t_force & force, float dt);
void vv_update_2(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics, t_force & force,float* box, float dt);


void Berendsen_update(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics,
                      t_tcstat & tcstat, t_tcpara & tcpara,t_group & tcgrp,
                      t_force & force,float *box, float dt);

void pbc_with_water(t_atom_kinematics & atom_kinematics,
                    t_atom_statics & atom_statics,
                    float* box,
                    int size,
                    array<t_molecule> & molelist
                    );



#endif
