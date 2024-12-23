#ifndef _SETTLE_H
#define _SETTLE_H
#include "includes.h"

#include "top.h"
#include "virial.h"
void settle(t_molecule & water,
            t_atom_kinematics & atom_kinematics,
            t_atom_kinematics & old_k,
            t_atom_statics & atom_statics,
			t_virial & virial,
			float invdt);

#endif