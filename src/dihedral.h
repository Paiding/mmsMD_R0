#ifndef _DIHEDRAL_H
#define _DEHEDRAL_H
#include "includes.h"
#include "top.h"
void do_pdih_force(array<t_molecule> & mole,
                   t_pdih_para & pdpara,
                   t_force & force,
                   t_atom_kinematics & atom_kinematics
                  );

void do_idih_force(array<t_molecule> & mole,
                   t_ipdih_para & idpara,
                   t_force & force,
                   t_atom_kinematics & atom_kinematics
                  );
#endif
