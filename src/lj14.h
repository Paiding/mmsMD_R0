#ifndef _LJ14_H
#define _LJ14_H
#include "includes.h"
#include "top.h"

void do_lj14_force(array<t_molecule> & mole,
                   t_lj14_para & lj14_para,
                   t_force & force,
                   t_atom_kinematics & atom_kinematics
                  );












#endif