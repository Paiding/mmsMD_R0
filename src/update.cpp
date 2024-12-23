#include "update.h"

void vv_update_1(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics,
                 t_force & force, float dt)
{
    int size = atom_kinematics.ax.getsize();
    
    for (int i = 0; i < size; i++)
    {
        atom_kinematics.ax[i] = force.fx[i] * atom_statics.m_inv[i];
        atom_kinematics.ay[i] = force.fy[i] * atom_statics.m_inv[i];
        atom_kinematics.az[i] = force.fz[i] * atom_statics.m_inv[i];
        atom_kinematics.vx[i] = atom_kinematics.vx[i] + 0.5 * atom_kinematics.ax[i] * dt;
        atom_kinematics.vy[i] = atom_kinematics.vy[i] + 0.5 * atom_kinematics.ay[i] * dt;
        atom_kinematics.vz[i] = atom_kinematics.vz[i] + 0.5 * atom_kinematics.az[i] * dt;
        atom_kinematics.rx[i] = atom_kinematics.rx[i] + atom_kinematics.vx[i] * dt;
        atom_kinematics.ry[i] = atom_kinematics.ry[i] + atom_kinematics.vy[i] * dt;
        atom_kinematics.rz[i] = atom_kinematics.rz[i] + atom_kinematics.vz[i] * dt;
    }
}

void vv_update_2(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics, t_force & force,float* box, float dt)
{
    int size = atom_kinematics.ax.getsize();
    float bx = box[XX];
	float by = box[YY];
	float bz = box[ZZ];
	float bxinv = 1.0f / bx;
	float byinv = 1.0f / by;
	float bzinv = 1.0f / bz;
    
    for (int i = 0; i < size; i++)
    {
        //atom_kinematics.ax[i] = force.fx[i] * atom_statics.m_inv[i];
        //atom_kinematics.ay[i] = force.fy[i] * atom_statics.m_inv[i];
        //atom_kinematics.az[i] = force.fz[i] * atom_statics.m_inv[i];
        atom_kinematics.vx[i] = atom_kinematics.vx[i] + 0.5 * atom_kinematics.ax[i] * dt;
        atom_kinematics.vy[i] = atom_kinematics.vy[i] + 0.5 * atom_kinematics.ay[i] * dt;
        atom_kinematics.vz[i] = atom_kinematics.vz[i] + 0.5 * atom_kinematics.az[i] * dt;

        //do pbc here
        int sx = rintf( (atom_kinematics.rx[i]-0.5f*bx)*bxinv );
        int sy = rintf( (atom_kinematics.ry[i]-0.5f*by)*byinv );
        int sz = rintf( (atom_kinematics.rz[i]-0.5f*bz)*bzinv );
  // PBC
        atom_kinematics.rx[i] -= sx * bx;
        atom_kinematics.ry[i] -= sy * by;
        atom_kinematics.rz[i] -= sz * bz;
    }
}

void Berendsen_update(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics,
                      t_tcstat & tcstat, t_tcpara & tcpara,t_group & tcgrp,
                      t_force & force,float *box, float dt)
{
    int size = atom_kinematics.ax.getsize();
    int tcgrp_id;
    
    //tmp
    float lambda;

    for (int i = 0; i < size; i++)
    {
        tcgrp_id = tcgrp.list[i];
        lambda = tcstat.lambda[tcgrp_id];
        lambda = 1.0;
        //printf("i=%d,rx=%f,vx=%f,ax=%f\n",i,atom_kinematics.rx[i],atom_kinematics.vx[i],atom_kinematics.ax[i]);
        atom_kinematics.ax[i] = force.fx[i] * atom_statics.m_inv[i];
        atom_kinematics.ay[i] = force.fy[i] * atom_statics.m_inv[i];
        atom_kinematics.az[i] = force.fz[i] * atom_statics.m_inv[i];
        atom_kinematics.vx[i] = lambda * (atom_kinematics.vx[i] + atom_kinematics.ax[i] * dt);
        atom_kinematics.vy[i] = lambda * (atom_kinematics.vy[i] + atom_kinematics.ay[i] * dt);
        atom_kinematics.vz[i] = lambda * (atom_kinematics.vz[i] + atom_kinematics.az[i] * dt);

        // printf("(%f,%f,%f)\n",dt,lambda,atom_statics.m_inv[i]);

        
        //printf("(%f,%f,%f)\n",atom_kinematics.vx[i],atom_kinematics.vy[i],atom_kinematics.vz[i]);
        atom_kinematics.rx[i] = atom_kinematics.rx[i] + atom_kinematics.vx[i] * dt;
        atom_kinematics.ry[i] = atom_kinematics.ry[i] + atom_kinematics.vy[i] * dt;
        atom_kinematics.rz[i] = atom_kinematics.rz[i] + atom_kinematics.vz[i] * dt;
        //printf("i=%d,rx=%f,vx=%f,ax=%f\n",i,atom_kinematics.rx[i],atom_kinematics.vx[i],atom_kinematics.ax[i]);

        
    }

}

void pbc_with_water(t_atom_kinematics & atom_kinematics,
                    t_atom_statics & atom_statics,
                    float* box,
                    int size,
                    array<t_molecule> & molelist
                    )
{
    size = atom_kinematics.rx.getsize();
    float bx = box[XX];
	float by = box[YY];
	float bz = box[ZZ];
	float bxinv = 1.0f / bx;
	float byinv = 1.0f / by;
	float bzinv = 1.0f / bz;
    for (int i = 0; i < size; i++)
    {
        if (!atom_statics.iswater[i])
        {
            //do pbc here
            float sx = rintf( (atom_kinematics.rx[i]-0.5f*bx)*bxinv );
            float sy = rintf( (atom_kinematics.ry[i]-0.5f*by)*byinv );
            float sz = rintf( (atom_kinematics.rz[i]-0.5f*bz)*bzinv );
        //printf("inner,%f,rx,%f\n",(atom_kinematics.rx[i]-0.5f*bx)*bxinv,atom_kinematics.rx[i]);
        //printf("sr = (%f,%f,%f)",sx,sy,sz);
  // PBC
            atom_kinematics.rx[i] -= sx * bx;
            atom_kinematics.ry[i] -= sy * by;
            atom_kinematics.rz[i] -= sz * bz;
        }
    }

    int mol_size = molelist.getsize();
    for (int i = 0; i < mol_size; i++)
    {
        if (molelist[i].iswater)
        {
            int iOW = molelist[i].atom_ranks[0];
            int iHW1 = molelist[i].atom_ranks[1];
            int iHW2 = molelist[i].atom_ranks[2];
            float ox = atom_kinematics.rx[iOW];
            float oy = atom_kinematics.ry[iOW];
            float oz = atom_kinematics.rz[iOW];

            if ((fabs(ox - 0.5*bx) > 0.5*bx)||(fabs(oy - 0.5*by) > 0.5*by)||(fabs(oz - 0.5*bz) > 0.5*bz))
            {
                float sx = rintf( (ox-0.5f*bx)*bxinv );
                float sy = rintf( (oy-0.5f*by)*byinv );
                float sz = rintf( (oz-0.5f*bz)*bzinv );
        
                atom_kinematics.rx[iOW] -= sx * bx;
                atom_kinematics.ry[iOW] -= sy * by;
                atom_kinematics.rz[iOW] -= sz * bz;
                atom_kinematics.rx[iHW1] -= sx * bx;
                atom_kinematics.ry[iHW1] -= sy * by;
                atom_kinematics.rz[iHW1] -= sz * bz;
                atom_kinematics.rx[iHW2] -= sx * bx;
                atom_kinematics.ry[iHW2] -= sy * by;
                atom_kinematics.rz[iHW2] -= sz * bz;
            }
        }
    }
}