#include <lincs.h>
void init_lincs_data(array<t_molecule> & molelist,
                    t_lincs_data & lincs_data,
                    t_atom_statics & atom_statics,
                    t_constr_para & cons_para
                    )
{//catch not-water bonds in molelist
    int size = molelist.getsize();
    int count = 0;
    for (int i = 0; i < size; i++)
    {
        if (molelist[i].iswater == 0)
        {
            int num = molelist[i].cons_num;
            for (int j = 0; j < num; j++)
            {
                int atomi = molelist[i].atom_ranks[molelist[i].consi[j]];
                int atomj = molelist[i].atom_ranks[molelist[i].consj[j]];
                lincs_data.A_list.put(atomi);
                lincs_data.B_list.put(atomj);
                int type = molelist[i].consType[j] - cons_para.type_index[0];
                lincs_data.len.put(cons_para.d0[type]);
                lincs_data.Sdiag.put(rsqrtf(atom_statics.m_inv[atomi] + atom_statics.m_inv[atomj]));
            }
            count += num;
        }
    }
    lincs_data.cons_num = count;
    int loc = 0;
    lincs_data.nbr_cons_loc.put(loc);
    float tmp;
    //i,j are ranks of constrains; a,b are ranks of (bridge) atoms
    for (int i = 0; i < count; i++)
    {
        int a = lincs_data.A_list[i];
        int b = lincs_data.B_list[i];

        //init all
        lincs_data.dir_x.put(0);
        lincs_data.dir_y.put(0);
        lincs_data.dir_z.put(0);
        
        lincs_data.rhs1.put(0);
        lincs_data.rhs2.put(0);
        lincs_data.sol.put(0);
        lincs_data.lambda.put(0);

        
        for (int j = 0; j < count; j++)
        {
            if (j == i) continue;
            if (lincs_data.A_list[j] == a)
            {
                lincs_data.nbr_cons.put(j);
                tmp = - atom_statics.m_inv[a] * lincs_data.Sdiag[i] * lincs_data.Sdiag[j];
                lincs_data.coef.put(tmp);
                lincs_data.blcc.put(0);
                loc ++;
            }
            else if (lincs_data.A_list[j] == b)
            {
                lincs_data.nbr_cons.put(j);
                tmp = atom_statics.m_inv[b] * lincs_data.Sdiag[i] * lincs_data.Sdiag[j];
                lincs_data.coef.put(tmp);
                lincs_data.blcc.put(0);
                loc ++;
            }
            else if (lincs_data.B_list[j] == b)
            {
                lincs_data.nbr_cons.put(j);
                tmp = -atom_statics.m_inv[b] * lincs_data.Sdiag[i] * lincs_data.Sdiag[j];
                lincs_data.coef.put(tmp);
                lincs_data.blcc.put(0);
                loc ++;
            }
            else if (lincs_data.B_list[j] == a)
            {
                lincs_data.nbr_cons.put(j);
                tmp = atom_statics.m_inv[a] * lincs_data.Sdiag[i] * lincs_data.Sdiag[j];
                lincs_data.coef.put(tmp);
                lincs_data.blcc.put(0);
                loc ++;
            }
        }

        lincs_data.nbr_cons_loc.put(loc);
        
    }
    
}

void calc_direct(t_atom_kinematics & old_k,
                 t_lincs_data & lincs_data
                )
{//calc direct for every bonds
    int atomi,atomj;
    float rlen,dx,dy,dz;
    //printf("%d",lincs_data.cons_num);
    for (int i = 0; i < lincs_data.cons_num; i++)
    {
        atomi = lincs_data.A_list[i];
        atomj = lincs_data.B_list[i];
        dx = old_k.rx[atomi] - old_k.rx[atomj];
        dy = old_k.ry[atomi] - old_k.ry[atomj];
        dz = old_k.rz[atomi] - old_k.rz[atomj];
        rlen = rsqrtf(dx*dx + dy*dy + dz*dz);
        dx *= rlen;
        dy *= rlen;
        dz *= rlen;
        //printf("dir(%f,%f,%f)\n",dx,dy,dz);
        lincs_data.dir_x[i] = dx;
        lincs_data.dir_y[i] = dy;
        lincs_data.dir_z[i] = dz;
    }
    

}

void b4expand_0(t_lincs_data & ld,
                t_atom_kinematics & ak
                )
{
    for (int i = 0; i < ld.cons_num; i++)
    {
        int begin_blnr = ld.nbr_cons_loc[i];
        int end_blnr = ld.nbr_cons_loc[i+1];
        int atomi = ld.A_list[i];
        int atomj = ld.B_list[i];
        for (int j = begin_blnr; j < end_blnr; j++)
        {
            int k = ld.nbr_cons[j];
            ld.blcc[j] = (ld.coef[j]*(ld.dir_x[i]*ld.dir_x[k] 
            + ld.dir_y[i]*ld.dir_y[k] 
            + ld.dir_z[i]*ld.dir_z[k]));
            
        }
        float mvb = ld.Sdiag[i] * (ld.dir_x[i] * (ak.rx[atomi] - ak.rx[atomj]) +
                                   ld.dir_y[i] * (ak.ry[atomi] - ak.ry[atomj]) +
                                   ld.dir_z[i] * (ak.rz[atomi] - ak.rz[atomj]) - ld.len[i]);
        ld.rhs1[i] = mvb;
        ld.sol[i] = mvb;
    }
    
}

void lincs_mat_exp_0(t_lincs_data & ld, int flag)
{
    if (flag == 1)
    {
       for (int i = 0; i < ld.cons_num; i++)
        {
            int begin_blnr = ld.nbr_cons_loc[i];
            int end_blnr = ld.nbr_cons_loc[i+1];

            float mvb = 0.0f;
            for (int j = begin_blnr; j < end_blnr; j++)
            {
                int k = ld.nbr_cons[j];
                mvb += ld.blcc[j] * ld.rhs1[k];
                //printf("%f,",ld.blcc[i]);
            }
            ld.rhs2[i] = mvb;//kill these puts
            ld.sol[i] += mvb;
            //printf("%d:%f,%d,%d\n",i,mvb,begin_blnr,end_blnr);
        
        }
    }
    else
    {
        for (int i = 0; i < ld.cons_num; i++)
        {
            int begin_blnr = ld.nbr_cons_loc[i];
            int end_blnr = ld.nbr_cons_loc[i+1];

            float mvb = 0.0f;
            for (int j = begin_blnr; j < end_blnr; j++)
            {
                int k = ld.nbr_cons[j];
                mvb += ld.blcc[j] * ld.rhs2[k];
                //printf("%f,",ld.blcc[i]);
            }
            ld.rhs1[i] = mvb;//kill these puts
            ld.sol[i] += mvb;
            //printf("%d:%f,%d,%d\n",i,mvb,begin_blnr,end_blnr);
        
        }
    }
}

void lincs_mat_exp_1(t_lincs_data & ld)
{
    for (int i = 0; i < ld.cons_num; i++)
    {
        int begin_blnr = ld.nbr_cons_loc[i];
        int end_blnr = ld.nbr_cons_loc[i+1];

        float mvb = 0.0f;
        for (int j = begin_blnr; j < end_blnr; j++)
        {
            int k = ld.nbr_cons[j];
            mvb += ld.blcc[j] * ld.rhs1[k];
        }
        ld.rhs2[i] = mvb;
        ld.sol[i] += mvb;
        ld.lambda[i] = -ld.sol[i] * ld.Sdiag[i];
    }
    
}
// correction for centripetal effects
void b4expand_1(t_lincs_data & ld,
                t_atom_kinematics & ak
                )
{
    int atomi,atomj;
    double dx,dy,dz,len,len2,dlen2;
    float mvb;
    for (int i = 0; i < ld.cons_num; i++)
    {
        atomi = ld.A_list[i];
        atomj = ld.B_list[i];
        dx = ak.rx[atomi] - ak.rx[atomj];
        dy = ak.ry[atomi] - ak.ry[atomj];
        dz = ak.rz[atomi] - ak.rz[atomj];

        len = ld.len[i];
        len2 = len*len;
        dlen2 = 2 * len2 - (dx*dx + dy*dy + dz*dz);
    
        mvb = 0.0f;
        if (dlen2 > 0)
        {
            mvb = ld.Sdiag[i] * (len - dlen2*rsqrtf(dlen2));
            //mvb = ld.Sdiag[i] * (len - sqrtf(dlen2));
            //printf(">%d,%f,%f,%f,%f,%f\n",i,mvb,ld.len[i],dlen2,rsqrtf(dlen2),len - dlen2*rsqrtf(dlen2));
            //printf("dx=%f,dy=%f,dz=%f\n",dx,dy,dz);
        }
        else{
            mvb = ld.Sdiag[i] * len;
        }

        ld.rhs1[i] = mvb;
        ld.sol[i] = mvb;
        
    }
    
}

void lincs_mat_exp_2(t_lincs_data & ld)
{
    for (int i = 0; i < ld.cons_num; i++)
    {
        int begin_blnr = ld.nbr_cons_loc[i];
        int end_blnr = ld.nbr_cons_loc[i+1];

        float mvb = 0.0f;
        for (int j = begin_blnr; j < end_blnr; j++)
        {
            int k = ld.nbr_cons[j];
            mvb += ld.blcc[j] * ld.rhs1[k];
        }
        ld.sol[i] += mvb;
        ld.lambda[i] += -ld.sol[i] * ld.Sdiag[i];
    }
}

void lincs_solve(t_lincs_data & ld,
                 t_atom_kinematics & ak,
                 t_atom_statics & as,
                 float invdt
                )
{
    int atomi,atomj;
    
    float mvb;
    for (int i = 0; i < ld.cons_num; i++)
    {
        atomi = ld.A_list[i];
        atomj = ld.B_list[i];
        //printf("$%d[%f,%f,%f],%f,%f,  %f\n",i,ld.dir_x[i],ld.dir_y[i],ld.dir_z[i],ld.Sdiag[i],ld.sol[i],ld.lambda[i]);
        //printf(">%f(%d,%d) %f,%f,%f\n",ld.lambda[i],atomi,atomj,ld.dir_x[i],ld.dir_y[i],ld.dir_z[i]);
        //printf(">%f(%d,%d) %f,%f,%f\n",ld.lambda[i],atomi,atomj,ld.dir_x[i] * ld.lambda[i] * invdt,ld.dir_y[i] * ld.lambda[i] * invdt,ld.dir_z[i] * ld.lambda[i] * invdt);
        ak.rx[atomi] -= as.m_inv[atomi] * ld.dir_x[i] * ld.Sdiag[i] * ld.sol[i];
        ak.ry[atomi] -= as.m_inv[atomi] * ld.dir_y[i] * ld.Sdiag[i] * ld.sol[i];
        ak.rz[atomi] -= as.m_inv[atomi] * ld.dir_z[i] * ld.Sdiag[i] * ld.sol[i];
        ak.rx[atomj] += as.m_inv[atomj] * ld.dir_x[i] * ld.Sdiag[i] * ld.sol[i];
        ak.ry[atomj] += as.m_inv[atomj] * ld.dir_y[i] * ld.Sdiag[i] * ld.sol[i];
        ak.rz[atomj] += as.m_inv[atomj] * ld.dir_z[i] * ld.Sdiag[i] * ld.sol[i];

        if (true)
        {
            ak.vx[atomi] -= as.m_inv[atomi] * ld.dir_x[i] * ld.lambda[i] * invdt;
            ak.vy[atomi] -= as.m_inv[atomi] * ld.dir_y[i] * ld.lambda[i] * invdt;
            ak.vz[atomi] -= as.m_inv[atomi] * ld.dir_z[i] * ld.lambda[i] * invdt;
            ak.vx[atomj] += as.m_inv[atomj] * ld.dir_x[i] * ld.lambda[i] * invdt;
            ak.vy[atomj] += as.m_inv[atomj] * ld.dir_y[i] * ld.lambda[i] * invdt;
            ak.vz[atomj] += as.m_inv[atomj] * ld.dir_z[i] * ld.lambda[i] * invdt;  
        }
        
    }
}

void lincs_virial(t_virial & virial,
                  t_lincs_data & ld,
                  float invdt
                 )
{
    for (int i = 0; i < ld.cons_num; i++)
    {
        float tmp0 = ld.len[i] * ld.lambda[i];
        float tmp1 = 0.0f;

        tmp1 = tmp0 * ld.dir_x[i];
        float vxx = tmp1 * ld.dir_x[i];

        tmp1 = tmp0 * ld.dir_y[i];
        float vyy = tmp1 * ld.dir_y[i];

        tmp1 = tmp0 * ld.dir_z[i];
        float vzz = tmp1 * ld.dir_z[i];

        float vir_factor = 0.5f * invdt * invdt;
        virial.vir_x[i] -= vir_factor * vxx;
        virial.vir_y[i] -= vir_factor * vyy;
        virial.vir_z[i] -= vir_factor * vzz;
    }
    
}

void do_lincs(t_lincs_data & lincs_data,
              t_atom_statics & atom_statics,
              t_atom_kinematics & atom_kinemetics,
              t_atom_kinematics & old_k,
              t_virial & virial,
              float invdt
            )
{
    int nOrder = 4;
    int nIter = 1;
    calc_direct(old_k, lincs_data);

    b4expand_0(lincs_data,atom_kinemetics);

    int flag = 1;
    for (int i = 0; i < (nOrder-1); i++)
    {
        lincs_mat_exp_0(lincs_data,flag);
        flag *= -1;
    }
    /*for (int i = 0; i < lincs_data.cons_num; i++)
    {
        printf("%d:%f\n",i,lincs_data.blcc[i]);
    }*/
    
    
    
    lincs_mat_exp_1(lincs_data);

    

    lincs_solve(lincs_data,atom_kinemetics,atom_statics,invdt);

    //centripetal

    for (int i = 0; i < nIter; i++)
    {
        b4expand_1(lincs_data,atom_kinemetics);
        flag = 1;
        for (int j = 0; j < nOrder; j++)
        {
            lincs_mat_exp_0(lincs_data,flag);
            flag *= -1;
        }

        lincs_mat_exp_2(lincs_data);
        lincs_solve(lincs_data,atom_kinemetics,atom_statics,invdt);
    }

    //last iteration
    b4expand_1(lincs_data,atom_kinemetics);
    flag = 1;
    for (int j = 0; j < nOrder; j++)
    {
        lincs_mat_exp_0(lincs_data,flag);
        flag *= -1;
    }
    
    lincs_mat_exp_2(lincs_data);
    lincs_solve(lincs_data,atom_kinemetics,atom_statics,invdt);

    lincs_virial(virial,lincs_data,invdt);

    
    
}