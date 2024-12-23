#include "bond.h"

void do_bondforce(array<t_molecule> & mole,
                  t_atom_kinematics & atom_kinematics,
                  t_atom_statics & atom_statics,
                  t_force & force,
                  t_virial & virial)
{
    int mole_size = mole.getsize();

#ifdef __P_COUPLE__
	// initialize the virial to 0
	float vxx = 0.0f;
// tmp ** 	float vxy = 0.0f;
// tmp ** 	float vxz = 0.0f;
// tmp ** 	float vyx = 0.0f;
	float vyy = 0.0f;
// tmp ** 	float vyz = 0.0f;
// tmp ** 	float vzx = 0.0f; 
// tmp ** 	float vzy = 0.0f;
	float vzz = 0.0f;
#endif

    for (int i = 0; i < mole_size; i++)
    {
        for (int j = 0; j < mole[i].bond_num; j++)
        {//cycle through all bonds in the molecule
            int atomi = mole[i].atom_ranks[mole[i].bondi[j]];//find the two atoms of the bond
            int atomj = mole[i].atom_ranks[mole[i].bondj[j]];

            float fx = 0.0f;
            float fy = 0.0f;
            float fz = 0.0f;

            float dx = atom_kinematics.rx[atomi] - atom_kinematics.rx[atomj];
            float dy = atom_kinematics.ry[atomi] - atom_kinematics.ry[atomj];
            float dz = atom_kinematics.rz[atomi] - atom_kinematics.rz[atomj];

            float rsq = dx*dx + dy*dy + dz*dz;
		    float rinv = rsqrtf(rsq);//calculate distance

            float r0 = mole[i].bondr0[j];//reference distance
            float rK = mole[i].bondrK[j];//bond energy parameter

            float fforce = rK*(r0*rinv-1.0f);//laplacian of potential energy
            
            force.esum_bond += fforce;
            

            float tfx = dx * fforce;//x portion of bond force
		    float tfy = dy * fforce;
		    float tfz = dz * fforce;
		    fx += tfx;
		    fy += tfy;
		    fz += tfz;
            
            force.fx[atomi] += fx;
            force.fy[atomi] += fy;
            force.fz[atomi] += fz;
            force.fx[atomj] -= fx;//!to be checked!
            force.fy[atomj] -= fy;
            force.fz[atomj] -= fz;

            #ifdef __P_COUPLE__
			// sum the virial
			//vxx += dx*tfx;

			//vyy += dy*tfy;

			//vzz += dz*tfz;
            //virial.vir_x[atomi] -= 0.5f * vxx;// (-dx)*(-dx) = dx * dx
            //virial.vir_y[atomi] -= 0.5f * vyy;
            //virial.vir_z[atomi] -= 0.5f * vzz;
            #endif
            
        }
    }
}