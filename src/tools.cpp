#include "tools.h"

void check_broken_molecule(float maxl,t_molecule & molecule,
                           t_atom_kinematics & atom_kinematics
                          )
{
    float sx = 0;
    float sy = 0;
    float sz = 0;
    float dx = 0;
    float dy = 0;
    float dz = 0;
    int size = molecule.atom_ranks.getsize();
    for (int i = 0; i < size; i++)
    {
        sx += atom_kinematics.rx[molecule.atom_ranks[i]];
        sy += atom_kinematics.ry[molecule.atom_ranks[i]];
        sz += atom_kinematics.rz[molecule.atom_ranks[i]];
    }
    //compute center of shape
    sx /= size;
    sy /= size;
    sz /= size;
    for (int i = 0; i < size; i++)
    {
        dx = atom_kinematics.rx[molecule.atom_ranks[i]] - sx;
        dy = atom_kinematics.ry[molecule.atom_ranks[i]] - sy;
        dz = atom_kinematics.rz[molecule.atom_ranks[i]] - sz;
        dx = fabs(dx);
        dy = fabs(dy);
        dz = fabs(dz);
        if (dx > maxl || dy > maxl || dz > maxl)
        {
            printf("%d\n",molecule.atom_ranks[i]);
        }
    }
}

void fix_broken_water(float maxl,t_molecule & molecule,
                      t_atom_kinematics & atom_kinematics,
                      float *box
                     )
{
    if (molecule.iswater != 0)
    {
        float ix[9];
        float ox = atom_kinematics.rx[molecule.atom_ranks[0]];
        float oy = atom_kinematics.ry[molecule.atom_ranks[0]];
        float oz = atom_kinematics.rz[molecule.atom_ranks[0]];

        float h1x = atom_kinematics.rx[molecule.atom_ranks[1]];
        float h1y = atom_kinematics.ry[molecule.atom_ranks[1]];
        float h1z = atom_kinematics.rz[molecule.atom_ranks[1]];

        float h2x = atom_kinematics.rx[molecule.atom_ranks[2]];
        float h2y = atom_kinematics.ry[molecule.atom_ranks[2]];
        float h2z = atom_kinematics.rz[molecule.atom_ranks[2]];

        if (h1x - ox > maxl)
        {
            atom_kinematics.rx[molecule.atom_ranks[1]] -= box[XX];
        }
        else if (ox - h1x > maxl)
        {
            atom_kinematics.rx[molecule.atom_ranks[1]] += box[XX];
        }

        if (h1y - oy > maxl)
        {
            atom_kinematics.ry[molecule.atom_ranks[1]] -= box[YY];
        }
        else if (oy - h1y > maxl)
        {
            atom_kinematics.ry[molecule.atom_ranks[1]] += box[YY];
        }

        if (h1z - oz > maxl)
        {
            atom_kinematics.rz[molecule.atom_ranks[1]] -= box[ZZ];
        }
        else if (oz - h1z > maxl)
        {
            atom_kinematics.rz[molecule.atom_ranks[1]] += box[ZZ];
        }

        if (h2x - ox > maxl)
        {
            atom_kinematics.rx[molecule.atom_ranks[2]] -= box[XX];
        }
        else if (ox - h2x > maxl)
        {
            atom_kinematics.rx[molecule.atom_ranks[2]] += box[XX];
        }

        if (h2y - oy > maxl)
        {
            atom_kinematics.ry[molecule.atom_ranks[2]] -= box[YY];
        }
        else if (oy - h2y > maxl)
        {
            atom_kinematics.ry[molecule.atom_ranks[2]] += box[YY];
        }

        if (h2z - oz > maxl)
        {
            atom_kinematics.rz[molecule.atom_ranks[2]] -= box[ZZ];
        }
        else if (oz - h2z > maxl)
        {
            atom_kinematics.rz[molecule.atom_ranks[2]] += box[ZZ];
        }
    }
}