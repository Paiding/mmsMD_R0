#include "t_couple.h"

void ekin_sum(t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics,
              t_tcstat & tcstat,t_tcpara & tcpara,t_group & tcgrp)
{
    int atom_size = atom_kinematics.ax.getsize();

//prepare for group coupling
    int ngroup = tcstat.N;
    
    tcstat.ekinh_old.clear();
    for (int i = 0; i < ngroup; i++)
    {//write the old ekin before refreshing
        tcstat.ekinh_old.put(tcstat.ekinh[i]);
    }
    
    double* ekin = new double[3*ngroup];
    double* ekin_sum = new double[ngroup];
    for (int i = 0; i < ngroup; i++)
    {
        ekin[3*i] = 0.0;
        ekin[3*i+1] = 0.0;
        ekin[3*i+2] = 0.0;
        ekin_sum[i] = 0.0;
    }
    int grp_id;
    
    tcstat.ekinh.clear();
    for (int i = 0; i < atom_size; i++)
    {//in KJ/mol if /MILLION
        grp_id = tcgrp.list[i];
        
        ekin[3*grp_id] += atom_kinematics.vx[i] * atom_kinematics.vx[i] * 0.5f * atom_statics.m[i];
        //ekin += atom_kinematics.vx[i] * atom_kinematics.vy[i];
        //ekin += atom_kinematics.vx[i] * atom_kinematics.vz[i];
        ekin[3*grp_id+1] += atom_kinematics.vy[i] * atom_kinematics.vy[i] * 0.5f * atom_statics.m[i];
        //ekin += atom_kinematics.vy[i] * atom_kinematics.vz[i];
        ekin[3*grp_id+2] += atom_kinematics.vz[i] * atom_kinematics.vz[i] * 0.5f * atom_statics.m[i];
        
    }
    
    tcstat.ekin_tot[XX] = 0;
    tcstat.ekin_tot[YY] = 0;
    tcstat.ekin_tot[ZZ] = 0;
    for (int i = 0; i < ngroup; i++)
    {
        tcstat.ekin_tot[XX] += ekin[3*i];
        tcstat.ekin_tot[YY] += ekin[3*i+1];
        tcstat.ekin_tot[ZZ] += ekin[3*i+2];
        ekin_sum[i] += ekin[3*i];
        ekin_sum[i] += ekin[3*i+1];
        ekin_sum[i] += ekin[3*i+2];
        
    }
    
    //printf("ekin:%f,%f,%f\n",tcstat.ekin_tot[XX] ,tcstat.ekin_tot[YY] ,tcstat.ekin_tot[ZZ]);
    
    
    
    
    tcstat.ekin.clear();
    tcstat.T.clear();
    tcstat.Th.clear();
    tcstat.T_sum = 0;
    tcstat.nrdf_sum = 0;
    for (int i = 0; i < ngroup; i++)
    {
        tcstat.ekinh.put(ekin_sum[i]);
        ekin_sum[i] = 0.5 * (tcstat.ekinh[i] + tcstat.ekinh_old[i]);
        tcstat.ekin.put(ekin_sum[i]);
        tcstat.T.put(2.0f * ekin_sum[i] / (tcpara.nrdf[i] * float(BOLTZ)));
        tcstat.Th.put(2.0f * tcstat.ekinh[i] / (tcpara.nrdf[i] * float(BOLTZ)));
        
        tcstat.T_sum += tcstat.T[i];
        tcstat.nrdf_sum += tcpara.nrdf[i];
    }
}
//to be arranged for group coupling
float Berendsen_lambda(t_tcpara & tcpara, t_tcstat & tcstat,
                        float dt, int index)
{
    float ref_t = tcpara.ref_t[index];
    float tau_t = tcpara.tau_t[index];
    //origin lambda to be fixed in (0.8,1.25)
    float lambda_o = sqrtf(1.0f + ((dt/tau_t)*((ref_t/tcstat.T[index]) - 1.0f)));
    //printf("%d,%f,%f,%f,%f\n",index,ref_t,tau_t,tcstat.T[index],dt);
    //float lambda = ((lambda_o <=1.25f ? lambda_o : 1.25f) >= 0.8f ? (lambda_o <=1.25f ? lambda_o : 1.25f) : 0.8f);
    if (lambda_o > 1.25)
    {
        return 1.25;
    }
    else if(lambda_o < 0.8)
    {
        return 0.80;
    }
    else
    {
        return lambda_o;
    }
}

void gen_vel(t_atom_kinematics & atom_kinematics,
             t_atom_statics & atom_statics,
             array<t_molecule> & molelist,
             t_tcpara & tcpara
            )
{
    int mole_size = molelist.getsize();
    float Tg = tcpara.ref_t[0];
    Tg /= 3.0;
    float ekg = Tg * BOLTZ;/*g/mol * nm^2 / ns^2 */

    int seed = ini_get_int("GENV_SEED","in/settings.ini");
    if (seed == -1)
    {
        return;
    }
    
    srand(seed);
    for (int i = 0; i < mole_size; i++)
    {
        float m = 0.0f;
        int n = molelist[i].atom_ranks.getsize();
        for (int j = 0; j < n; j++)
        {
            int rank = molelist[i].atom_ranks[j];
            m += atom_statics.m[rank];
        }
        float vg = sqrtf(ekg / m);
        float vgx = 1 - (0.00002 * (rand() % 100000));
        float vgy = 1 - (0.00002 * (rand() % 100000));
        float vgz = 1 - (0.00002 * (rand() % 100000));
        for (int j = 0; j < n; j++)
        {
            int rank = molelist[i].atom_ranks[j];
            atom_kinematics.vx[rank] += vgx * vg;
            atom_kinematics.vy[rank] += vgy * vg;
            atom_kinematics.vz[rank] += vgz * vg;
        }
    }
}

void init_tc_group(t_group & tcgrp, int grp_num)
{
    char filename[64];
    char line[256];
    char dest1[32];
    char dest2[32];
    char dest3[32];
    int size;
    int tmp;

    char *tik;
    ini_get_str("dump_file","in/settings.ini",filename);
    tcgrp.index = (int*)malloc(sizeof(int)*grp_num);
    for (int i = 0; i < grp_num; i++)
    {
        tcgrp.index[i] = 0;
    }

    FILE *fp;
    if (fp = fopen(filename,"r"))
    {
        while (true)
	    {
		    fgets(line,128,fp);
            sscanf(line,"%s %d",dest1,&size);
            if (strcmp(dest1,"allocated") == 0)
            {
                break;
            }
        }
        tcgrp.list = (int*)malloc(sizeof(int)*size);
        for (int i = 0; i < size; i++)
        {
            fgets(line,128,fp);//name = xx
            sscanf(line,"%s %s %s %d",dest1,dest2,dest3,&tmp);
            tcgrp.list[i] = tmp;
            tcgrp.index[tmp] ++;
        }
    }
}

void init_tc_nogroup(t_group & tcgrp,int atom_size)
{
    tcgrp.index = (int*)malloc(sizeof(int));
    tcgrp.index[0] = 0;
    tcgrp.list = (int*)malloc(sizeof(int)*atom_size);
    for (int i = 0; i < atom_size; i++)
    {
        tcgrp.list[i] = 0;
    }
    

}