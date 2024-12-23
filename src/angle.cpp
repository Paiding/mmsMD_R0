#include "angle.h"

void do_angleforce(array<t_molecule> & mole,
                  t_atom_kinematics & atom_kinematics,
                  t_atom_statics & atom_statics,
                  t_force & force,
                  t_virial & virial)
{
    int mole_size = mole.getsize();
    //printf("angle to be compute is %d",mole_size);
    for (int i = 0; i < mole_size; i++)
    {
        for (int j = 0; j < mole[i].angle_num; j++)
        {//cycle through all angles in one molecule
            if (mole[i].angleType[j] == 0)
            {
            
            int atomi = mole[i].atom_ranks[mole[i].anglei[j]];//find the three atoms of the angle
            int atomj = mole[i].atom_ranks[mole[i].anglej[j]];
            int atomk = mole[i].atom_ranks[mole[i].anglek[j]];
            //printf("i = %d,j = %d,k = %d\n",atomi,atomj,atomk);

            float dxij = atom_kinematics.rx[atomi] - atom_kinematics.rx[atomj];//distance portion
            float dyij = atom_kinematics.ry[atomi] - atom_kinematics.ry[atomj];
            float dzij = atom_kinematics.rz[atomi] - atom_kinematics.rz[atomj];
            float dxki = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomi];
            float dyki = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomi];
            float dzki = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomi];
            float dxkj = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomj];
            float dykj = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomj];
            float dzkj = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomj];

            float d2ij = dxij * dxij + dyij * dyij + dzij * dzij;//square distance between three atoms
            float d2ki = dxki * dxki + dyki * dyki + dzki * dzki;
            float d2jk = dxkj * dxkj + dykj * dykj + dzkj * dzkj;

            //printf("d2ij = %f,d2ki = %f,d2jk = %f",d2ij,d2ki,d2jk);

            float ipab = d2ij * d2jk;//square of the denominator of cos-vec
            float ipab_inv_1 = 0.0f;// 1 / ||a|| * ||b||
            float cos_theta = 0.0f;

            if(ipab>0) 
            {
			    ipab_inv_1 = rsqrtf(ipab);
			    cos_theta = (dxij*dxkj + dyij*dykj + dzij*dzkj) * ipab_inv_1;
		    } 
            else 
            {
			    cos_theta = 1.0f;
		    }

		    if(cos_theta>1.0f) 
            {
			    cos_theta = 1.0f;
		    } 
            else if(cos_theta<-1.0f) 
            {
			    cos_theta = -1.0f;
		    }
            //printf("costheta = %f\n",cos_theta);
            float theta = 180 * acosf(cos_theta) / M_PI;
            

            float a0 = mole[i].anglea0[j];//reference angle
            float aK = mole[i].angleaK[j];//standard angle force potential
            //printf("theta = %f,a0 = %f,aK = %f\n",theta,a0,aK);
            //float dtheta = theta - a0;//for virial

            //V_ijk = 0.5aK * (a0 - theta) ^ 2
            
            float fforce =aK * (a0 - theta);//differential of potential energy   dV_ijk
            
            force.esum_angle += fforce;
            
            //printf("angle fforce = %f\n",fforce);
            float st = 0.0f;// dV/sintheta
            float sin_theta2 = 1.0f - cos_theta * cos_theta;
            if(sin_theta2 > 0)
			st = fforce*rsqrtf(sin_theta2);
		    float sth = st * cos_theta;
            
            float cik = st * ipab_inv_1;
      		float cii = sth / d2ij;
		    float ckk = sth / d2jk;

            //atom i
            float fxi = cii * dxij - cik * dxkj;
            float fyi = cii * dyij - cik * dykj;
            float fzi = cii * dzij - cik * dzkj;

            //atom k
            float fxk = ckk * dxkj - cik * dxij;
            float fyk = ckk * dykj - cik * dyij;
            float fzk = ckk * dzkj - cik * dzij;

            //atom j
            float fxj = -fxi - fxk;
            float fyj = -fyi - fyk;
            float fzj = -fzi - fzk;
            
            //printf("angles:[%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",fxi,fyi,fzi,fxk,fyk,fzk,fxj,fyj,fzj);
            force.fx[atomi] += fxi;
            force.fy[atomi] += fyi;
            force.fz[atomi] += fzi;
            force.fx[atomj] += fxj;
            force.fy[atomj] += fyj;
            force.fz[atomj] += fzj;
            force.fx[atomk] += fxk;
            force.fy[atomk] += fyk;
            force.fz[atomk] += fzk;
            }

#ifdef __P_COUPLE__
	        //virial.vir_x[atomi] -= 0.5 * fxi * dxij;
            //virial.vir_y[atomi] -= 0.5 * fyi * dyij;
            //virial.vir_z[atomi] -= 0.5 * fzi * dzij;
            //second atom has no virial contribution
            //virial.vir_x[atomk] -= 0.5 * fxk * dxkj;
            //virial.vir_y[atomk] -= 0.5 * fyk * dykj;
            //virial.vir_z[atomk] -= 0.5 * fzk * dzkj;
#endif
            
        }   
    }
}

void do_ub_angle(t_ub_para & ubpara, 
                 array<t_molecule> & mole,
                 t_atom_kinematics & atom_kinematics,
                 t_force & force)
{
    int mole_size = mole.getsize();

    for (int i = 0; i < mole_size; i++)
    {
        for (int j = 0;  j < mole[i].angle_num; j++)
        {//cycle through all angles in one molecule
            if (mole[i].angleType[j] == 1)
            {
        
            int atomi = mole[i].atom_ranks[mole[i].anglei[j]];//find the three atoms of the angle
            int atomj = mole[i].atom_ranks[mole[i].anglej[j]];
            int atomk = mole[i].atom_ranks[mole[i].anglek[j]];
            //printf("i = %d,j = %d,k = %d\n",atomi,atomj,atomk);

        /*

       i|\               Do some triangle calculation:
        | \              
    a   |  \                               ax*bx + ay*by + az*bz
        |   \              cos-vec (a,b) = ---------------------
       j|____\ k                                ||a|| * ||b||
           b                        ->      ->
                                a = ji; b = jk;
        */

            float dxij = atom_kinematics.rx[atomi] - atom_kinematics.rx[atomj];//distance portion
            float dyij = atom_kinematics.ry[atomi] - atom_kinematics.ry[atomj];
            float dzij = atom_kinematics.rz[atomi] - atom_kinematics.rz[atomj];
            float dxki = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomi];
            float dyki = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomi];
            float dzki = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomi];
            float dxkj = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomj];
            float dykj = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomj];
            float dzkj = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomj];

            float d2ij = dxij * dxij + dyij * dyij + dzij * dzij;//square distance between three atoms
            float d2ki = dxki * dxki + dyki * dyki + dzki * dzki;
            float d2jk = dxkj * dxkj + dykj * dykj + dzkj * dzkj;

            //printf("d2ij = %f,d2ki = %f,d2jk = %f",d2ij,d2ki,d2jk);

            float ipab = d2ij * d2jk;//square of the denominator of cos-vec
            float ipab_inv_1 = 0.0f;// 1 / ||a|| * ||b||
            float cos_theta = 0.0f;

            if(ipab>0) 
            {
			    ipab_inv_1 = rsqrtf(ipab);
			    cos_theta = (dxij*dxkj + dyij*dykj + dzij*dzkj) * ipab_inv_1;
		    } 
            else 
            {
			    cos_theta = 1.0f;
		    }

		    if(cos_theta>1.0f) 
            {
			    cos_theta = 1.0f;
		    } 
            else if(cos_theta<-1.0f) 
            {
			    cos_theta = -1.0f;
		    }
            //printf("costheta = %f\n",cos_theta);
            float theta = acosf(cos_theta);
            float r13 = sqrtf(d2ki);
            int ftype = mole[i].funcType_a[j] - ubpara.type_index[0];
            //printf("%f\n",theta-ubpara.theta0[ftype]);
            float fforce = ubpara.ktheta[ftype]*(theta-ubpara.theta0[ftype]*DEG2RAD);
            float fforce_b =  ubpara.kUB[ftype]*(r13-ubpara.r130[ftype]);//differential of potential energy   dV_ijk
            //printf("[%d,%d,%d]theta=%f;va = %f;vb = %f\n",atomi,atomj,atomk,theta,fforce,fforce_b);
            force.esum_angle += fforce*(theta-ubpara.theta0[ftype]*DEG2RAD) + fforce_b*(r13-ubpara.r130[ftype]);
            
            //printf("angle fforce = %f\n",fforce);
            float st = 0.0f;// dV/sintheta
            float sin_theta2 = 1.0f - cos_theta * cos_theta;
            if(sin_theta2 > 0)
			st = fforce*rsqrtf(sin_theta2);
		    float sth = st * cos_theta;
            
            float cik = st * ipab_inv_1;
      		float cii = sth / d2ij;
		    float ckk = sth / d2jk;

            //atom i
            float fxi = cii * dxij - cik * dxkj;
            float fyi = cii * dyij - cik * dykj;
            float fzi = cii * dzij - cik * dzkj;

            //atom k
            float fxk = ckk * dxkj - cik * dxij;
            float fyk = ckk * dykj - cik * dyij;
            float fzk = ckk * dzkj - cik * dzij;

            //atom j
            float fxj = -fxi - fxk;
            float fyj = -fyi - fyk;
            float fzj = -fzi - fzk;
            
            //printf("angles:[%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",fxi,fyi,fzi,fxk,fyk,fzk,fxj,fyj,fzj);
            force.fx[atomi] += fxi;
            force.fy[atomi] += fyi;
            force.fz[atomi] += fzi;
            force.fx[atomj] += fxj;
            force.fy[atomj] += fyj;
            force.fz[atomj] += fzj;
            force.fx[atomk] += fxk;
            force.fy[atomk] += fyk;
            force.fz[atomk] += fzk;
            /*
            force.fx[atomi] -= fxi;
            force.fy[atomi] -= fyi;
            force.fz[atomi] -= fzi;
            force.fx[atomj] -= fxj;
            force.fy[atomj] -= fyj;
            force.fz[atomj] -= fzj;
            force.fx[atomk] -= fxk;
            force.fy[atomk] -= fyk;
            force.fz[atomk] -= fzk;
*/
            //bond phase
            fforce_b *= rsqrtf(d2ki);
            float fxik = fforce_b * dxki;
            float fyik = fforce_b * dyki;
            float fzik = fforce_b * dzki;
            //printf("bonds:[%f,%f,%f]\n",fxik,fyik,fzik);
            //printf("dr:[%f,%f,%f]\n",dxki,dyki,dzki);
            force.fx[atomi] += fxik;
            force.fy[atomi] += fyik;
            force.fz[atomi] += fzik;
            force.fx[atomk] -= fxik;
            force.fy[atomk] -= fyik;
            force.fz[atomk] -= fzik;



            }
        }
    }
    
}

void do_g96_angle(t_g96a_para & g96a_para, 
                 array<t_molecule> & mole,
                 t_atom_kinematics & atom_kinematics,
                 t_force & force)
{
    int mole_size = mole.getsize();

    for (int i = 0; i < mole_size; i++)
    {
        for (int j = 0;  j < mole[i].angle_num; j++)
        {//cycle through all angles in one molecule
            if (mole[i].angleType[j] == 2)
            {
        
                int atomi = mole[i].atom_ranks[mole[i].anglei[j]];//find the three atoms of the angle
                int atomj = mole[i].atom_ranks[mole[i].anglej[j]];
                int atomk = mole[i].atom_ranks[mole[i].anglek[j]];
                //printf("i = %d,j = %d,k = %d\n",atomi,atomj,atomk);

                float dxij = atom_kinematics.rx[atomi] - atom_kinematics.rx[atomj];//distance portion
                float dyij = atom_kinematics.ry[atomi] - atom_kinematics.ry[atomj];
                float dzij = atom_kinematics.rz[atomi] - atom_kinematics.rz[atomj];
                float dxki = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomi];
                float dyki = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomi];
                float dzki = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomi];
                float dxkj = atom_kinematics.rx[atomk] - atom_kinematics.rx[atomj];
                float dykj = atom_kinematics.ry[atomk] - atom_kinematics.ry[atomj];
                float dzkj = atom_kinematics.rz[atomk] - atom_kinematics.rz[atomj];

                float d2ij = dxij * dxij + dyij * dyij + dzij * dzij;//square distance between three atoms
                float d2ki = dxki * dxki + dyki * dyki + dzki * dzki;
                float d2jk = dxkj * dxkj + dykj * dykj + dzkj * dzkj;
                
                //g96 angle special
                /*
                float c11 = d2ij;
                float c12 = dxij * dxki + dyij * dyki + dzij * dzki;
                float c22 = d2jk;
                float cD = sqrtf(c11*c22);

                float cos_theta = 0.0f;

                if(cD>0) 
                {
			        cos_theta = c12/cD;
		        } 
                else 
                {
			        cos_theta = 1.0f;
		        }*/
                float ipab = d2ij * d2jk;//square of the denominator of cos-vec
                float ipab_inv_1 = 0.0f;// 1 / ||a|| * ||b||
                float cos_theta = 0.0f;

                if(ipab>0) 
                {
			        ipab_inv_1 = rsqrtf(ipab);
			        cos_theta = (dxij*dxkj + dyij*dykj + dzij*dzkj) * ipab_inv_1;
		        } 
                else 
                {
			        cos_theta = 1.0f;
		        }

		        if(cos_theta>1.0f) 
                {
			        cos_theta = 1.0f;
		        } 
                else if(cos_theta<-1.0f) 
                {
			        cos_theta = -1.0f;
		        }
            
                int ftype = mole[i].funcType_a[j] - g96a_para.type_index[0];
                float a0 = g96a_para.theta0[ftype];
                float aK = g96a_para.ct[ftype];
                float fforce = aK * (a0 - cos_theta);
                //printf("[%d,%d,%d]theta=%f;va = %f;vb = %f\n",atomi,atomj,atomk,theta,fforce,fforce_b);
                force.esum_angle += 0.5 * aK * (a0 - cos_theta) * (a0 - cos_theta);
            
                float rij_1 = 1 / sqrtf(d2ij);
		        float rkj_1 = 1 / sqrtf(d2jk);
		        float rij_2 = rij_1*rij_1;
		        float rkj_2 = rkj_1*rkj_1;
		        float rijrkj_1 = rij_1*rkj_1;

                //atom i
                float fxi = fforce*(dxkj*rijrkj_1 - dxij*rij_2*cos_theta);
                float fyi = fforce*(dykj*rijrkj_1 - dyij*rij_2*cos_theta);
                float fzi = fforce*(dzkj*rijrkj_1 - dzij*rij_2*cos_theta);

                //atom k
                float fxk = fforce*(dxij*rijrkj_1 - dxkj*rij_2*cos_theta);
                float fyk = fforce*(dyij*rijrkj_1 - dykj*rij_2*cos_theta);
                float fzk = fforce*(dzij*rijrkj_1 - dzkj*rij_2*cos_theta);

                //atom j
                float fxj = -fxi - fxk;
                float fyj = -fyi - fyk;
                float fzj = -fzi - fzk;
            
                //printf("angles:[%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",fxi,fyi,fzi,fxk,fyk,fzk,fxj,fyj,fzj);
                force.fx[atomi] += fxi;
                force.fy[atomi] += fyi;
                force.fz[atomi] += fzi;
                force.fx[atomj] += fxj;
                force.fy[atomj] += fyj;
                force.fz[atomj] += fzj;
                force.fx[atomk] += fxk;
                force.fy[atomk] += fyk;
                force.fz[atomk] += fzk;
            }
        }
    }
    
}