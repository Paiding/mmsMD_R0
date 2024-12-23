#ifndef _BOND_H
#define _BOND_H
#include "includes.h"

#include "top.h"
#include "virial.h"

void do_bondforce(array<t_molecule> & mole,
                  t_atom_kinematics & atom_kinematics,
                  t_atom_statics & atom_statics,
                  t_force & force,
                  t_virial & virial);

#endif