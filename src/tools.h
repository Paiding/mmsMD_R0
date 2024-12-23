#ifndef _TOOLS_H_
#define _TOOLS_H_

#include "includes.h"
#include "top.h"
void check_broken_molecule(float maxl,t_molecule & molecule,
                           t_atom_kinematics & atom_kinematics
                          );

void fix_broken_water(float maxl,t_molecule & molecule,
                      t_atom_kinematics & atom_kinematics,
                      float *box
                     );

#endif