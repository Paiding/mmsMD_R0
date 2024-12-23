#include "lj14.h"
void do_lj14_force(array<t_molecule> & mole,
                   t_lj14_para & lj14_para,
                   t_force & force,
                   t_atom_kinematics & atom_kinematics
                  )
{
    int mole_size = mole.getsize();
    //printf("angle to be compute is %d",mole_size);
    for (int i = 0; i < mole_size; i++)
    {
        for (int j = 0; j < mole[i].lj14_num; j++)
        {//cycle through all angles in one molecule
            int atomi = mole[i].atom_ranks[mole[i].lj14i[j]];
            int atomj = mole[i].atom_ranks[mole[i].lj14j[j]];

            float dx = atom_kinematics.rx[atomi] - atom_kinematics.rx[atomj];
            float dy = atom_kinematics.ry[atomi] - atom_kinematics.ry[atomj];
            float dz = atom_kinematics.rz[atomi] - atom_kinematics.rz[atomj];

            int type = mole[i].lj14Type[j] - lj14_para.type_index[0];
            float lj6 = lj14_para.c6[type];
            float lj12 = lj14_para.c12[type];

            float rsq = dx*dx + dy*dy + dz*dz;
		    float r2inv = rsqrtf(rsq);//calculate distance

            float r6inv = r2inv * r2inv * r2inv;
            float fforce = r2inv * r6inv * (12.0f*lj12 * r6inv - 6.0f*lj6);
            //printf("%d,%d,%f\n",atomi,atomj,fforce);

            force.fx[atomi] += fforce * dx;
            force.fy[atomi] += fforce * dy;
            force.fz[atomi] += fforce * dz;
            force.fx[atomj] -= fforce * dx;
            force.fy[atomj] -= fforce * dy;
            force.fz[atomj] -= fforce * dz;
        }
    }
}