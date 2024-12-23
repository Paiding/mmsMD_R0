#include "vcm.h"

void calc_mass_grp(t_atom_statics & atom_statics,
              t_vcm & vcm
             )
{
    vcm.m = 0.0f;
    for (int i = 0; i < vcm.n; i++)
    {
        vcm.m += atom_statics.m[i];
    }
}

void calc_vcm_grp(t_atom_statics & atom_statics,
                  t_atom_kinematics & atom_kinematics,
                  t_vcm & vcm
                 )
{
    //tmp cal all atoms 
    double sx;
    double sy;
    double sz;
    array<double> ipx;
    array<double> ipy;
    array<double> ipz;
    for (int i = 0; i < vcm.n; i++)
    {
        ipx.put(atom_statics.m[i] * atom_kinematics.vx[i]);
        ipy.put(atom_statics.m[i] * atom_kinematics.vy[i]);
        ipz.put(atom_statics.m[i] * atom_kinematics.vz[i]);
    }
    sx = ipx.reduce();
    sy = ipy.reduce();
    sz = ipz.reduce();
    vcm.group_p_x = (float)sx;
    vcm.group_p_y = (float)sy;
    vcm.group_p_z = (float)sz;
}

void vcm_v_calc(t_vcm & vcm)
{
    vcm.vcm_x = vcm.group_p_x / vcm.m;
    vcm.vcm_y = vcm.group_p_y / vcm.m;
    vcm.vcm_z = vcm.group_p_z / vcm.m;
}

void calc_vcm_ang_grp(t_atom_statics & atom_statics,
                  t_atom_kinematics & atom_kinematics,
                  t_vcm & vcm
                 )
{
    double sx;
    double sy;
    double sz;
    array<double> ilx;
    array<double> ily;
    array<double> ilz;
    array<double> mx;
    array<double> my;
    array<double> mz;
    array<double> mxx;
    array<double> mxy;
    array<double> mxz;
    array<double> myy;
    array<double> myz;
    array<double> mzz;
    double jx;
    double jy;
    double jz;
    for (int i = 0; i < vcm.n; i++)
    {
        jx = atom_kinematics.ry[i]*atom_kinematics.vz[i] - atom_kinematics.rz[i]*atom_kinematics.vy[i];
        jy = atom_kinematics.rz[i]*atom_kinematics.vx[i] - atom_kinematics.rx[i]*atom_kinematics.vz[i];
        jz = atom_kinematics.rx[i]*atom_kinematics.vy[i] - atom_kinematics.ry[i]*atom_kinematics.vx[i];
        ilx.put(atom_statics.m[i]*jx);
        ily.put(atom_statics.m[i]*jy);
        ilz.put(atom_statics.m[i]*jz);
        mx.put(atom_statics.m[i]*atom_kinematics.rx[i]);
        my.put(atom_statics.m[i]*atom_kinematics.ry[i]);
        mz.put(atom_statics.m[i]*atom_kinematics.rz[i]);
        mxx.put(atom_statics.m[i]*atom_kinematics.rx[i]*atom_kinematics.rx[i]);
        mxy.put(atom_statics.m[i]*atom_kinematics.rx[i]*atom_kinematics.ry[i]);
        mxz.put(atom_statics.m[i]*atom_kinematics.rx[i]*atom_kinematics.rz[i]);
        myy.put(atom_statics.m[i]*atom_kinematics.ry[i]*atom_kinematics.ry[i]);
        myz.put(atom_statics.m[i]*atom_kinematics.ry[i]*atom_kinematics.rz[i]);
        mzz.put(atom_statics.m[i]*atom_kinematics.rz[i]*atom_kinematics.rz[i]);
    }
    vcm.l_x = ilx.reduce();
    vcm.l_y = ily.reduce();
    vcm.l_z = ilz.reduce();
    vcm.mx = mx.reduce();
    vcm.my = my.reduce();
    vcm.mz = mz.reduce();
    vcm.mr2[0] = mxx.reduce();
    vcm.mr2[1] = mxy.reduce();
    vcm.mr2[2] = mxz.reduce();
    vcm.mr2[3] = myy.reduce();
    vcm.mr2[4] = myz.reduce();
    vcm.mr2[5] = mzz.reduce();
}

void vcm_w_calc(t_vcm & vcm)
{
    float i_mass = vcm.m;
    float inv_mass = 1.0f / i_mass;
    float group_rx = vcm.mx;
    float group_ry = vcm.my;
    float group_rz = vcm.mz;
    group_rx *= inv_mass;
    group_ry *= inv_mass;
    group_rz *= inv_mass;

    float group_vx = vcm.vcm_x;
    float group_vy = vcm.vcm_y;
    float group_vz = vcm.vcm_z;

    float cprodx = group_ry*group_vz - group_rz*group_vy;
    float cprody = group_rz*group_vx - group_rx*group_vz;
    float cprodz = group_rx*group_vy - group_ry*group_vx;

		// wrute back
		vcm.mx = group_rx;
		vcm.my = group_ry;
		vcm.mz = group_rz;

		vcm.l_x -= i_mass*cprodx;
		vcm.l_y -= i_mass*cprody;
		vcm.l_z -= i_mass*cprodz;
    
    float i_tmp[9];
    i_tmp[XXXX] = vcm.mr2[0];
    i_tmp[XXYY] = vcm.mr2[1];
    i_tmp[XXZZ] = vcm.mr2[2];
    i_tmp[YYXX] = vcm.mr2[1];
    i_tmp[YYYY] = vcm.mr2[3];
    i_tmp[YYZZ] = vcm.mr2[4];
    i_tmp[ZZXX] = vcm.mr2[2];
    i_tmp[ZZYY] = vcm.mr2[4];
    i_tmp[ZZZZ] = vcm.mr2[5];

    float ixx = i_tmp[XXXX] - group_rx*group_rx*i_mass;
    float ixy = i_tmp[XXYY] - group_rx*group_ry*i_mass;
    float ixz = i_tmp[XXZZ] - group_rx*group_rz*i_mass;
    float iyy = i_tmp[YYYY] - group_ry*group_ry*i_mass;
    float iyz = i_tmp[YYZZ] - group_ry*group_rz*i_mass;
    float izz = i_tmp[ZZZZ] - group_rz*group_rz*i_mass;
    float iyx = ixy;
    float izx = ixz;
    float izy = iyz;

    float tmp[DIM][DIM];
    tmp[0][0] = iyy + izz;
    tmp[1][0] = -ixy;
    tmp[2][0] = -ixz;
    tmp[0][1] = -ixy;
    tmp[1][1] = ixx + izz;
    tmp[2][1] = -iyz;
    tmp[0][2] = -ixz;
    tmp[1][2] = -iyz;
    tmp[2][2] = ixx + iyy;

    float fac = 3.0f / (tmp[0][0]+tmp[1][1]+tmp[2][2]);
    tmp[0][0] *= fac;
    tmp[1][0] *= fac;
    tmp[2][0] *= fac;
    tmp[0][1] *= fac;
    tmp[1][1] *= fac;
    tmp[2][1] *= fac;
    tmp[0][2] *= fac;
    tmp[1][2] *= fac;
    tmp[2][2] *= fac;

    float det_1 = 1.0f / ( tmp[0][0]*(tmp[1][1]*tmp[2][2]-tmp[2][1]*tmp[1][2])
				                  -tmp[1][0]*(tmp[0][1]*tmp[2][2]-tmp[2][1]*tmp[0][2])
                          +tmp[2][0]*(tmp[0][1]*tmp[1][2]-tmp[1][1]*tmp[0][2]) );

    ixx =  det_1*(tmp[1][1]*tmp[2][2] - tmp[2][1]*tmp[1][2]);
    ixy = -det_1*(tmp[0][1]*tmp[2][2] - tmp[2][1]*tmp[0][2]);
    ixz =  det_1*(tmp[0][1]*tmp[1][2] - tmp[1][1]*tmp[0][2]);
    iyx = -det_1*(tmp[1][0]*tmp[2][2] - tmp[2][0]*tmp[1][2]);
    iyy =  det_1*(tmp[0][0]*tmp[2][2] - tmp[2][0]*tmp[0][2]);
    iyz = -det_1*(tmp[0][0]*tmp[1][2] - tmp[1][0]*tmp[0][2]);
    izx =  det_1*(tmp[1][0]*tmp[2][1] - tmp[2][0]*tmp[1][1]);
    izy = -det_1*(tmp[0][0]*tmp[2][1] - tmp[2][0]*tmp[0][1]);
    izz =  det_1*(tmp[0][0]*tmp[1][1] - tmp[1][0]*tmp[0][1]);

    ixx *= fac;
    ixy *= fac;
    ixz *= fac;
    iyx *= fac;
    iyy *= fac;
    iyz *= fac;
    izx *= fac;
    izy *= fac;
    izz *= fac;

    float jx = vcm.l_x;
    float jy = vcm.l_y;
    float jz = vcm.l_z;

    vcm.w_x = ixx*jx + ixy*jy + ixz*jz;
    vcm.w_y = iyx*jx + iyy*jy + iyz*jz;
    vcm.w_z = izx*jx + izy*jy + izz*jz;
}

void do_stopcm_grp(t_atom_kinematics & atom_kinematics,
                   t_vcm & vcm,
                   int is_angular
                  )
{
    int size = atom_kinematics.vx.getsize();
    float dvx = vcm.vcm_x;
    float dvy = vcm.vcm_y;
    float dvz = vcm.vcm_z;
    for (int i = 0; i < size; i++)
    {
        atom_kinematics.vx[i] -= dvx;
        atom_kinematics.vy[i] -= dvy;
        atom_kinematics.vz[i] -= dvz;
    }

    if (is_angular == 1)
    {
        float grx = vcm.mx;
        float gry = vcm.my;
        float grz = vcm.mz;

        float gwx = vcm.w_x;
        float gwy = vcm.w_y;
        float gwz = vcm.w_z;

        for (int i = 0; i < size; i++)
        {
            float drx = atom_kinematics.rx[i] - grx;
            float dry = atom_kinematics.rx[i] - grx;
            float drz = atom_kinematics.rx[i] - grx;

            float dvx = gwy*drz - gwz*dry;
		    float dvy = gwz*drx - gwx*drz;
		    float dvz = gwx*dry - gwy*drx;

		    atom_kinematics.vx[i] -= dvx;
            atom_kinematics.vy[i] -= dvy;
            atom_kinematics.vz[i] -= dvz;
        }
        
    }
    
}
