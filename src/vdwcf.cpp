#include "vdwcf.h"

float dump_shift(float alpha, float r)
{
    float A = erfcf(alpha * r) / r * r;
    float B1 = 2 * alpha / sqrtf(M_PI);
    float B2 = expf(- alpha * alpha * r * r) / r;
    return A + (B1 * B2);

}

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))

void init_shift_switch_para(t_radius & radius,
                            t_shift_switch_para & sspara
                            )
{
    double p,rc,r1=0.0;
    double A=0,B=0,C=0,A_3=0,B_4=0;
    
    p = 6.0;
    rc = radius.r_vdw;
    //printf("rc=%f\n",rc);
    //r1 = 0.0;
	A = p * ((p+1)*r1-(p+4)*rc)/(pow(rc,p+2)*pow2(rc-r1));
	B = -p * ((p+1)*r1-(p+3)*rc)/(pow(rc,p+2)*pow3(rc-r1));
	C = 1.0/pow(rc,p)-A/3.0*pow3(rc-r1)-B/4.0*pow4(rc-r1);
	A=-A;
	B=-B;
	C=-C;
	A_3=A/3.0;
	B_4=B/4.0;

    sspara.LJ_6_A = A;
    sspara.LJ_6_B = B;
    sspara.LJ_6_C = C;
    sspara.LJ_6_A_3 = A_3;
    sspara.LJ_6_B_4 = B_4;
    //printf("(%f,%f,%f,%f,%f)\n",A,B,C,A_3,B_4);
    // LJ12
	p=12.0;
	rc = radius.r_vdw;
    //r1 = 0.0;
	A = p * ((p+1)*r1-(p+4)*rc)/(pow(rc,p+2)*pow2(rc-r1));
	B = -p * ((p+1)*r1-(p+3)*rc)/(pow(rc,p+2)*pow3(rc-r1));
	C = 1.0/pow(rc,p)-A/3.0*pow3(rc-r1)-B/4.0*pow4(rc-r1);
	A_3=A/3.0;
	B_4=B/4.0;

    sspara.LJ_12_A = A;
    sspara.LJ_12_B = B;
    sspara.LJ_12_C = C;
    sspara.LJ_12_A_3 = A_3;
    sspara.LJ_12_B_4 = B_4;
    //printf("(%f,%f,%f,%f,%f\n)",A,B,C,A_3,B_4);
}

void vdwcf_dsf_ss_img(const array<float> & lj_para,
                      t_dsf_para & dsf_para,
                      t_shift_switch_para & sspara,
                      const t_atom_statics & atom_statics,
                      const t_atom_kinematics & atom_kinematics,
                      const t_nblist & nblist,
                      t_force & force, 
                      t_virial & virial,
                      float * box,
                      int atom_size
                      )
{
    float bx = box[XX];
	float by = box[YY];
	float bz = box[ZZ];
	float bxinv = 1.0f / bx;
	float byinv = 1.0f / by;
	float bzinv = 1.0f / bz;
    //int size = nblist.size;
    int size = atom_size;
    
    int atnr = atom_statics.atnr;

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

    float LJ_12_A = sspara.LJ_12_A;
    float LJ_12_B = sspara.LJ_12_B;
    float LJ_12_C = sspara.LJ_12_C;
    float LJ_12_A_3 = sspara.LJ_12_A_3;
    float LJ_12_B_4 = sspara.LJ_12_B_4;
    float LJ_6_A = sspara.LJ_6_A;
    float LJ_6_B = sspara.LJ_6_B;
    float LJ_6_C = sspara.LJ_6_C;
    float LJ_6_A_3 = sspara.LJ_6_A_3;
    float LJ_6_B_4 = sspara.LJ_6_B_4;


    for (int i = 0; i < size; i++)
    {
        int typei = atom_statics.type[i];
        //float lj1 = lj_para[2 * type];
        //float lj2 = lj_para[2 * type + 1];
        
        float ix = atom_kinematics.rx[i];
        float iy = atom_kinematics.ry[i];
        float iz = atom_kinematics.rz[i];
        
        //printf("%d\n",inner_size);
        //printf("{x=%f,y=%f,z=%f}\n",ix,iy,iz);
        for (int k = 0; k < size; k++)
        {
            int t = nblist.va_list[i][k];
            
            if (t == -1)
                continue;
            
            int tx = t / 9;
            int ty = (t % 9) / 3;
            int tz = t % 3;
            tx -= 1;
            ty -= 1;
            tz -= 1;
            //printf("%d,%d,%d\n",tx,ty,tz);
            int typek = atom_statics.type[k];
            float c6 = lj_para[2 * (atnr * typei + typek)];
            float c12 = lj_para[2 * (atnr * typei + typek) + 1];
            int type = (atnr * typei + typek);
            
            float dx = ix - atom_kinematics.rx[k] - tx * bx;
            float dy = iy - atom_kinematics.ry[k] - ty * by;
            float dz = iz - atom_kinematics.rz[k] - tz * bz;
            //printf("i %d,j %d,k %d [%d to %d] [%f,%f,%f]\n",tx,ty,tz,i,k,dx,dy,dz);
            //dx -= bx * rintf(dx * bxinv);
			//dy -= by * rintf(dy * byinv);
			//dz -= bz * rintf(dz * bzinv);
            //printf("type %d [%d,%d] = (%f,%f,%f)\n",t,i,k,dx,dy,dz);
            
			float rsq = dx*dx + dy*dy + dz*dz;
            
			float rsq_1 = sqrtf(rsq);
            float r3 = rsq_1 * rsq;
            float r4 = rsq * rsq;
			float r1inv = rsqrtf(rsq);
			float r2inv = 1.0f / rsq;
            

            // vdw force
			float r6inv = r2inv*r2inv*r2inv;
            float r12inv = r6inv * r6inv;
			
            float fforce = 
					c12*(12.0f*r12inv*r1inv + LJ_12_A*rsq + LJ_12_B*r3)*r1inv // repulsion
					- c6*(6.0f*r6inv*r1inv - LJ_6_A*rsq - LJ_6_B*r3)*r1inv; // dispersion
            fforce = r2inv * r6inv *(12.0f*c12 * r6inv - 6.0f*c6);
            //fforce = 0.0;
            //printf("%f\n",fforce);

			force.esum_vdw += 
					c6*(-r6inv - LJ_6_A_3*r3 - LJ_6_B_4*r4 - LJ_6_C) // dispersion
					+ c12*(r12inv - LJ_12_A_3*r3 - LJ_12_B_4*r4 - LJ_12_C); // repulsion

            
            // dsf Coulumb force
            
            
            
            fforce += K_COULUMB * r2inv * r1inv * atom_statics.q[i] * atom_statics.q[k] * (dump_shift(dsf_para.alpha, rsq_1) - dsf_para.Rc_force);
            //fforce += K_COULUMB * r2inv * r1inv * atom_statics.q[i] * atom_statics.q[k];
            //printf("%f,%f\n",dump_shift(dsf_para.alpha, rsq_1),dsf_para.Rc_force);
            
            force.esum_qq += atom_statics.q[i] * atom_statics.q[k] * r1inv;

            float tfx = dx * fforce;
            float tfy = dy * fforce;
            float tfz = dz * fforce;
            //printf("f = %f\n",tfx);
            force.fx[i] += tfx; 
            force.fy[i] += tfy; 
            force.fz[i] += tfz; 

        }
    }
}
