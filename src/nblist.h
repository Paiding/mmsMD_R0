#ifndef _NBLIST_H
#define _NBLIST_H

#include "includes.h"
#include "top.h"
#include "pretreat.h"
#define EXCL_MAX 32

struct t_nblist
{
    int ** va_list = NULL;
    int ** bin_list;
    int size = 0;

};

struct t_excl_list
{//list of exclusion atoms
    int** bonded_ranks = NULL;
    int sizei = 0;
};
void build_excl_list(array<t_molecule> & molelist,
                    t_mole_map & mole_map,
                     t_excl_list & excl_list,
                     int atom_size
                    );

void build_nblist_img(float *box, 
                        t_nblist & nblist_basic,
                        int atom_size,const t_radius & radius,
                        t_atom_kinematics & atom_kinematics,
                        t_excl_list & excl_list);

void build_nblist_bin(float *box, 
                      t_nblist & nblist,
                      int atom_size,const t_radius & radius,
                      t_atom_kinematics & atom_kinematics,
                      t_excl_list & excl_list);


#endif
