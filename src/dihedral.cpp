#include "dihedral.h"

void do_pdih_force(array<t_molecule> & mole,
                   t_pdih_para & pdpara,
                   t_force & force,
                   t_atom_kinematics & atom_kinematics
                  )
{
    int mole_size = mole.getsize();
    //printf("angle to be compute is %d",mole_size);
    for (int i = 0; i < mole_size; i++)
    {
        for (int j = 0; j < mole[i].pdih_num; j++)
        {//cycle through all angles in one molecule
            int atomi = mole[i].atom_ranks[mole[i].pdihi[j]];
            int atomj = mole[i].atom_ranks[mole[i].pdihj[j]];
            int atomk = mole[i].atom_ranks[mole[i].pdihk[j]];
            int atoml = mole[i].atom_ranks[mole[i].pdihl[j]];
            //printf("i = %d,j = %d,k = %d,l = %d\n",atomi,atomj,atomk,atoml);


            float dxij = atom_kinematics.rx[atomi] - atom_kinematics.rx[atomj];//distance portion
            float dyij = atom_kinematics.ry[atomi] - atom_kinematics.ry[atomj];
            float dzij = atom_kinematics.rz[atomi] - atom_kinematics.rz[atomj];
            float dxkj = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomj];
            float dykj = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomj];
            float dzkj = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomj];
            float dxkl = atom_kinematics.rx[atomk] - atom_kinematics.rx[atoml];
            float dykl = atom_kinematics.ry[atomk] - atom_kinematics.ry[atoml];
            float dzkl = atom_kinematics.rz[atomk] - atom_kinematics.rz[atoml];

            float mx = dyij * dzkj - dzij * dykj;
            float my = dzij * dxkj - dxij * dzkj;
            float mz = dxij * dykj - dyij * dxkj;
            float nx = dykj * dzkl - dzkj * dykl;
            float ny = dzkj * dxkl - dxkj * dzkl;
            float nz = dxkj * dykl - dykj * dxkl;

            float mw = mx * mx + my * my + mz * mz;
            float nw = nx * nx + ny * ny + nz * nz;
            float ip = mx * nx + my * ny + mz * nz;

            float cos_phi = 0.0f;
            float ipmn = mw * nw;

            if(ipmn > 0)
                cos_phi = ip * rsqrtf(ipmn);
            else
                cos_phi = 1.0f;

            if (cos_phi > 1.0f)
                cos_phi = 1.0f;
            
            if (cos_phi < -1.0f)
                cos_phi = -1.0f;
            
            float phi = acosf(cos_phi);
            ip = dxij * nx + dyij * ny + dzij * nz;
            float sign = (ip<0.0f) ? -1.0f : 1.0f;
            phi *= sign;

            //compute scalar force
            int type = mole[i].pdihType[j];
            type = type - pdpara.type_index[0];
            float phi0 = pdpara.phi0[type];
            float k0 = pdpara.cp0[type];
            float mult = pdpara.mult[type];
            float ddphi = k0 * mult * sinf(phi0 - mult*phi);
            force.esum_pdih += k0 * (1.0f + cosf(mult*phi-phi0));

            float rwkj = dxkj * dxkj + dykj * dykj + dzkj * dzkj;
            float r_kj_1 = 0.0f;
            if (rwkj > 0)
            {   
                r_kj_1 = sqrtf(rwkj);
            }
            float inv_rwkj = 1 / rwkj;
            float inv_r_kj_1 = 1 / r_kj_1;
            
            //common
            float f1,f4;
            float p = dxij * dxkj + dyij * dykj + dzij * dzkj;
            float q = dxkl * dxkj + dykl * dykj + dzkl * dzkj;
            p /= rwkj;
            q /= rwkj;
            float a = -ddphi * r_kj_1 / mw;
            float b = ddphi * r_kj_1 / nw;
            
            f1 = a*mx;
            f4 = b*nx;
            force.fx[atomi] += f1;
            force.fx[atoml] += f4;
            force.fx[atomj] += p*f1 - q*f4 - f1;
            force.fx[atomk] += q*f4 - p*f1 - f4;

            f1 = a*my;
            f4 = b*ny;
            
            force.fy[atomi] += f1;
            force.fy[atoml] += f4;
            force.fy[atomj] += p*f1 - q*f4 - f1;
            force.fy[atomk] += q*f4 - p*f1 - f4;

            f1 = a*mz;
            f4 = b*nz;
            force.fz[atomi] += f1;
            force.fz[atoml] += f4;
            force.fz[atomj] += p*f1 - q*f4 - f1;
            force.fz[atomk] += q*f4 - p*f1 - f4;
        }
    }
}

void do_idih_force(array<t_molecule> & mole,
                   t_ipdih_para & idpara,
                   t_force & force,
                   t_atom_kinematics & atom_kinematics
                  )
{
    int mole_size = mole.getsize();
    //printf("angle to be compute is %d",mole_size);
    for (int i = 0; i < mole_size; i++)
    {
        for (int j = 0; j < mole[i].ipdih_num; j++)
        {
            //copied from proper dihedral
            int atomi = mole[i].atom_ranks[mole[i].ipdihi[j]];
            int atomj = mole[i].atom_ranks[mole[i].ipdihj[j]];
            int atomk = mole[i].atom_ranks[mole[i].ipdihk[j]];
            int atoml = mole[i].atom_ranks[mole[i].ipdihl[j]];
            //printf("i = %d,j = %d,k = %d,l = %d\n",atomi,atomj,atomk,atoml);
            float dxij = atom_kinematics.rx[atomi] - atom_kinematics.rx[atomj];//distance portion
            float dyij = atom_kinematics.ry[atomi] - atom_kinematics.ry[atomj];
            float dzij = atom_kinematics.rz[atomi] - atom_kinematics.rz[atomj];
            float dxkj = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomj];
            float dykj = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomj];
            float dzkj = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomj];
            float dxkl = atom_kinematics.rx[atomk] - atom_kinematics.rx[atoml];
            float dykl = atom_kinematics.ry[atomk] - atom_kinematics.ry[atoml];
            float dzkl = atom_kinematics.rz[atomk] - atom_kinematics.rz[atoml];

            float mx = dyij * dzkj - dzij * dykj;
            float my = dzij * dxkj - dxij * dzkj;
            float mz = dxij * dykj - dyij * dxkj;
            float nx = dykj * dzkl - dzkj * dykl;
            float ny = dzkj * dxkl - dxkj * dzkl;
            float nz = dxkj * dykl - dykj * dxkl;

            float mw = mx * mx + my * my + mz * mz;
            float nw = nx * nx + ny * ny + nz * nz;
            float ip = mx * nx + my * ny + mz * nz;

            float cos_phi = 0.0f;
            float ipmn = mw * nw;

            if(ipmn > 0)
                cos_phi = ip * rsqrtf(ipmn);
            else
                cos_phi = 1.0f;

            if (cos_phi > 1.0f)
                cos_phi = 1.0f;
            
            if (cos_phi < -1.0f)
                cos_phi = -1.0f;
            
            float phi = acosf(cos_phi);
            ip = dxij * nx + dyij * ny + dzij * nz;
            float sign = (ip<0.0f) ? -1.0f : 1.0f;
            phi *= sign;


            //ipdih only
            //scalar force
            int type = mole[i].ipdihType[j];
            type = type - idpara.type_index[0];
            float e0 = idpara.xi0[type];
            float k0 = idpara.cx0[type];
            float dp = phi - e0;
            //printf("e0:%f,k0:%f,dp:%f\n",e0,k0,dp);
            if (dp >= float(M_PI))
            {
                dp -= 2.0f*float(M_PI);
            }
            else if (dp < - float(M_PI))
            {
                dp += 2.0f*float(M_PI);
            }
            
            float ddphi = k0 * dp;
            force.esum_idih += 0.5f *k0 * dp * dp;

            float rwkj = dxkj * dxkj + dykj * dykj + dzkj * dzkj;
            float r_kj_1 = 0.0f;
            if (rwkj > 0)
            {   
                r_kj_1 = sqrtf(rwkj);
            }
            
            float f1,f4;
            float p = dxij * dxkj + dyij * dykj + dzij * dzkj;
            p /= rwkj;
            float q = dxkl * dxkj + dykl * dykj + dzkl * dzkj;
            q /= rwkj;
            float a = -ddphi * r_kj_1 / mw;
            float b = ddphi * r_kj_1 / nw;
            
            f1 = a*mx;
            f4 = b*nx;
            force.fx[atomi] += f1;
            force.fx[atoml] += f4;
            force.fx[atomj] += p*f1 - q*f4 - f1;
            force.fx[atomk] += q*f4 - p*f1 - f4;
            //printf("[%d,%d,%d,%d] %f,%f,%f,%f\n",atomi,atomj,atomk,atoml,f1,p*f1 - q*f4 - f1,q*f4 - p*f1 - f4,f4);
            //printf("%f,%f,%f,%f\n",1/r_kj_1,mw,nw,1/rwkj);
            //printf("%f,%f,%f,%f\n",a,b,p,q);
            
            f1 = a*my;
            f4 = b*ny;
            
            force.fy[atomi] += f1;
            force.fy[atoml] += f4;
            force.fy[atomj] += p*f1 - q*f4 - f1;
            force.fy[atomk] += q*f4 - p*f1 - f4;

            f1 = a*mz;
            f4 = b*nz;
            force.fz[atomi] += f1;
            force.fz[atoml] += f4;
            force.fz[atomj] += p*f1 - q*f4 - f1;
            force.fz[atomk] += q*f4 - p*f1 - f4;
        }
    }
}