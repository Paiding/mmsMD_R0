#ifndef _ANGLE_H
#define _ANGLE_H
#include "includes.h"

#include "top.h"
#include "virial.h"

void do_angleforce(array<t_molecule> & mole,
                  t_atom_kinematics & atom_kinematics,
                  t_atom_statics & atom_statics,
                  t_force & force,
                  t_virial & virial);

void do_ub_angle(t_ub_para & ubpara, 
                 array<t_molecule> & mole,
                 t_atom_kinematics & atom_kinematics,
                 t_force & force);

void do_g96_angle(t_g96a_para & g96a_para, 
                 array<t_molecule> & mole,
                 t_atom_kinematics & atom_kinematics,
                 t_force & force);
#endif