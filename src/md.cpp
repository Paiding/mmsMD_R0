#include "md.h"

int init_all(t_atom_kinematics & atom_kinematics,
              t_atom_kinematics & old_k,
              t_atom_statics & atom_statics,
              t_tcstat & tcstat, t_tcpara & tcpara,
              t_pcstat & pcstat, t_pcpara & pcpara,
              t_force & force,
              t_virial & virial,
              t_nstep & step,
              t_radius & radius,
              array<t_molecule> & mole_array,
              t_mole_map & mole_map,
              t_dsf_para & dsf_para,
              t_vcm & vcm
              )
{

    //init atom data 
    int atom_size = 0;
    char strline[256];
    char filename[64];
    char tmp[256];
    bool fileNotEnd = true;
    ini_get_str("dump_file","in/settings.ini",filename);
    FILE *fp = NULL;

    if(fp = fopen(filename,"r"))
    {
        while (fileNotEnd)
        {
            if (fgets(strline, MAX_STRLEN, fp)!=NULL)
            {//reaching the end/
                if (strncmp("x", strline, 1) == 0)
                {//leftvalue = key?
                    strcpy(tmp,strline);
                    char *tik = strtok(strline, " ");
                    strTrim(tik);
                    if (strcmp(tik,"x") == 0)
                    {//case: "keynote" to "key" cannot be read as data
                        tik = strtok(tmp,"(");
                        tik = strtok(NULL,"x");
                        atom_size = atoi(tik);
                        int timer = 0;//count the serial of atom
                        int moltimer = 0;//count the serial of molecule
                        int moloffset;
                        int angle_loc = 0;
                        int pdih_loc = 0;
                        int ipdih_loc = 0;
                        int lj14_loc = 0;
                        int cons_loc = 0;
                        for (int i = 0; i < mole_map.moletypes; i++)
                        {//rank is one type of molecule
                            int loc = mole_map.atom_loc[i];
                            
                            bool iswater = false;
                            if (strcmp("SOL",mole_map.mole_type[i]) == 0)
                            {   
                                iswater = true;
                            }
                            moloffset = moltimer;
                            for (int j = moloffset; j < moloffset + mole_map.mole_num[i]; j++)
                            {//rank is serial of molecule array
                                mole_array[j].mole_map_rank = i;
                                mole_array[j].iswater = iswater;
                                strcpy(mole_array[j].mole_type,mole_map.mole_type[i]);
                                
                                for (int k = 0; k < mole_map.atom_num[i]; k++)
                                {//rank is atom serial in one molecule
                                    fgets(strline, MAX_STRLEN, fp);
                                    tik = strtok(strline,"{");
	                                tik = strtok(NULL,",");
	                                float xx = atof(tik);
	                                
	                                tik = strtok(NULL,",");
	                                float yy = atof(tik);
	                                
	                                tik = strtok(NULL,"}");
	                                float zz = atof(tik);
	                                
                                    atom_kinematics.rx.put(xx);
                                    atom_kinematics.ry.put(yy);
                                    atom_kinematics.rz.put(zz);
                                    mole_array[j].atom_ranks.put(timer);
                                    
                                    float mass = mole_map.atom_mass[loc+k];

                                    atom_statics.m.put(mass);
                                    atom_statics.m_inv.put(1/mass);
                                    atom_statics.type.put(mole_map.atom_type[loc+k]);
                                    atom_statics.q.put(mole_map.atom_q[loc+k]);
                                    atom_statics.mole_serial.put(j);
                                    atom_statics.iswater.put(iswater);
                                    timer ++;
                                }
                                mole_array[j].angle_num = mole_map.angle_num[i];
                                for (int ii = 0; ii < mole_map.angle_num[i]; ii++)
                                {
                                    mole_array[j].funcType_a.put(mole_map.funcType_a[ii+angle_loc]);
                                    mole_array[j].angleType.put(mole_map.angleType[ii+angle_loc]);
                                    mole_array[j].anglei.put(mole_map.anglei[ii+angle_loc]);
                                    mole_array[j].anglej.put(mole_map.anglej[ii+angle_loc]);
                                    mole_array[j].anglek.put(mole_map.anglek[ii+angle_loc]);
                                }
                                mole_array[j].pdih_num = mole_map.pdih_num[i];
                                for (int ii = 0; ii < mole_map.pdih_num[i]; ii++)
                                {
                                    mole_array[j].pdihType.put(mole_map.funcType_pdih[ii+pdih_loc]);
                                    
                                    mole_array[j].pdihi.put(mole_map.pdihi[ii+pdih_loc]);
                                    mole_array[j].pdihj.put(mole_map.pdihj[ii+pdih_loc]);
                                    mole_array[j].pdihk.put(mole_map.pdihk[ii+pdih_loc]);
                                    mole_array[j].pdihl.put(mole_map.pdihl[ii+pdih_loc]);
                                }
                                mole_array[j].ipdih_num = mole_map.ipdih_num[i];
                                for (int ii = 0; ii < mole_map.ipdih_num[i]; ii++)
                                {
                                    mole_array[j].ipdihType.put(mole_map.funcType_ipdih[ii+ipdih_loc]);
                                    
                                    mole_array[j].ipdihi.put(mole_map.ipdihi[ii+ipdih_loc]);
                                    mole_array[j].ipdihj.put(mole_map.ipdihj[ii+ipdih_loc]);
                                    mole_array[j].ipdihk.put(mole_map.ipdihk[ii+ipdih_loc]);
                                    mole_array[j].ipdihl.put(mole_map.ipdihl[ii+ipdih_loc]);
                                }
                                mole_array[j].lj14_num = mole_map.lj14_num[i];
                                for (int ii = 0; ii < mole_map.lj14_num[i]; ii++)
                                {
                                    mole_array[j].lj14Type.put(mole_map.funcType_lj14[ii+lj14_loc]);
                                    
                                    mole_array[j].lj14i.put(mole_map.lj14i[ii+lj14_loc]);
                                    mole_array[j].lj14j.put(mole_map.lj14j[ii+lj14_loc]);
                                }
                                mole_array[j].cons_num = mole_map.cons_num[i];
                                for (int ii = 0; ii < mole_map.cons_num[i]; ii++)
                                {
                                    mole_array[j].consType.put(mole_map.funcType_constr[ii+lj14_loc]);
                                    
                                    mole_array[j].consi.put(mole_map.constri[ii+cons_loc]);
                                    mole_array[j].consj.put(mole_map.constrj[ii+cons_loc]);
                                }
                                moltimer ++;
                            }
                            angle_loc += mole_map.angle_num[i];
                            pdih_loc += mole_map.pdih_num[i];
                            ipdih_loc += mole_map.ipdih_num[i];
                            lj14_loc += mole_map.lj14_num[i];
                            cons_loc += mole_map.cons_num[i];
                            
                        }
                        
                        fileNotEnd = false;
                    }
                }
            }
            else
            {
                fileNotEnd = false;
            }
        }

        /*
        if (false)
        {
            while (fileNotEnd)
            {
                if (fgets(strline, MAX_STRLEN, fp)!=NULL)
                {//reaching the end/
                    if (strncmp("v", strline, 1) == 0)
                    {//leftvalue = key?
                        strcpy(tmp,strline);
                        char *tik = strtok(strline, " ");
                        strTrim(tik);
                        if (strcmp(tik,"v") == 0)
                        {//case: "keynote" to "key" cannot be read as data
                            tik = strtok(tmp,"(");
                            tik = strtok(NULL,"v");
                            atom_size = atoi(tik);
                            for (int i = 0; i < atom_size; i++)
                            {
                                fgets(strline, MAX_STRLEN, fp);
                                tik = strtok(strline,"{");
	                            tik = strtok(NULL,",");
	                            float vx = atof(tik);
	                                
	                            tik = strtok(NULL,",");
	                            float vy = atof(tik);
	                                
	                            tik = strtok(NULL,"}");
	                            float vz = atof(tik);
                                atom_kinematics.vx[i] = vx;
                                atom_kinematics.vy[i] = vy;
                                atom_kinematics.vz[i] = vz;
                            }
                        }
                    }
                }
            }
                            

        }*/
        
        
    }
    


    

    //simulation basic parameter
    step.full_step = ini_get_int("step","in/settings.ini");
    step.n_log = ini_get_int("step_log","in/settings.ini");
    step.n_nblist = ini_get_int("step_nblist","in/settings.ini");
    radius.r_vdw = ini_get_float("rvdw","in/settings.ini");
    radius.r_dsf = ini_get_float("rdsf","in/settings.ini");
    dsf_para.alpha = ini_get_float("dsf_alpha","in/settings.ini");
    dsf_para.Rc_force = dump_shift(dsf_para.alpha, radius.r_dsf);

    //temperature coupling init
    tcpara.N = ini_get_int("ngroup","in/settings.ini");
    tcstat.N = tcpara.N;
    tcpara.nrdf = (float*)malloc(sizeof(float)*tcpara.N);
    //init tc group

    tcpara.nrdf[0] = 2 * atom_size - 3;
    tcpara.ref_t = (float*)malloc(sizeof(float)*tcpara.N);
    
    tcpara.tau_t = (float*)malloc(sizeof(float)*tcpara.N);
    
    tcstat.lambda = (float*)malloc(sizeof(float)*tcpara.N);

    find_float_list(filename,"nrdf:",tcpara.N,tcpara.nrdf);
    find_float_list(filename,"ref-t:",tcpara.N,tcpara.ref_t);
    find_float_list(filename,"tau-t:",tcpara.N,tcpara.tau_t);
    for (int i = 0; i < tcpara.N; i++)
    {
        tcstat.lambda[i] = 1.0;
        tcstat.ekin.put(0.0);
        tcstat.ekinh.put(0.0);
        tcstat.ekinh_old.put(0.0);
        tcstat.T.put(1.0);
        tcstat.Th.put(1.0);
    }
    
    

    //pressure coupling init
    pcpara.do_pc = ini_get_int("do_pc","in/settings.ini");
    pcpara.dt = ini_get_float("dt","in/settings.ini");
    pcpara.tau_p = ini_get_float("tau_p","in/settings.ini");

    vcm.n = atom_size;
    vcm.m = 0.0;
    vcm.group_p_x = 0.0;
    vcm.group_p_y = 0.0;
    vcm.group_p_z = 0.0;
    
    //create empty array
    for (int i = 0; i < atom_size; i++)
    {
        force.fx.put(0);
        force.fy.put(0);
        force.fz.put(0);
        virial.vir_x.put(0.0);
        virial.vir_y.put(0.0);
        virial.vir_z.put(0.0);
        atom_kinematics.vx.put(0.0);
        atom_kinematics.vy.put(0.0);
        atom_kinematics.vz.put(0.0);
        atom_kinematics.ax.put(0.0);
        atom_kinematics.ay.put(0.0);
        atom_kinematics.az.put(0.0);

    }
    //kinematics_copy(atom_kinematics,old_k);
    return atom_size;
}

void one_step(int step, t_nstep & nstep, float dt,
              float *box, array<float> & lj_para,
              t_dsf_para & dsf_para,
              t_shift_switch_para & sspara,
              t_ub_para & ubpara,t_g96a_para & g96a_para, t_pdih_para & pdpara, t_ipdih_para & idpara,t_lj14_para & lj14_para,
              int atom_size,t_radius & radius,
              t_nblist & nblist,
              t_atom_kinematics & atom_kinematics,
              t_atom_kinematics & old_k,
              t_atom_statics & atom_statics,
              t_tcstat & tcstat, t_tcpara & tcpara,
              t_pcstat & pcstat, t_pcpara & pcpara,
              t_lincs_data & lincs_data,
              t_force & force,
              t_virial & virial,
              int tctype,
              t_excl_list & excl_list,
              t_mole_map & mole_map,
              array<t_molecule> & molelist,
              t_vcm & vcm,
              t_group & tcgrp
              )
{
/* 1 . neighbor list */
	if (step % nstep.n_nblist == 1)
    {
        //build_excl_list(molelist,mole_map,excl_list,atom_size);
       
        //build_nblist_img(box, nblist, atom_size,radius, atom_kinematics,excl_list);
        build_nblist_bin(box, nblist, atom_size,radius, atom_kinematics,excl_list);
        
        /*
        for (int i = 0; i < nblist.size; i++)
	{
		int sizej = nblist.neighbor_ranks[i].getsize();
		for (int j = 0; j < sizej; j++)
		{
			printf("[%d,%d] ",nblist.neighbor_ranks[i][j],nblist.image_type[i][j]);
		}
		printf("\n");
	}*/
    }
    

    
    
    fclear(force,virial);
/* 2 . do force */
    
    
    vdwcf_dsf_ss_img(lj_para, dsf_para, sspara, atom_statics, atom_kinematics, nblist, force, virial, box, atom_size);
    
    
    //do_bondforce(molelist,atom_kinematics,atom_statics,force,virial);


    //do_angleforce(molelist,atom_kinematics,atom_statics,force,virial);
    do_lj14_force(molelist,lj14_para,force,atom_kinematics);
    
    do_ub_angle(ubpara,molelist,atom_kinematics,force);
    do_g96_angle(g96a_para,molelist,atom_kinematics,force);

    do_pdih_force(molelist,pdpara,force,atom_kinematics);
    do_idih_force(molelist,idpara,force,atom_kinematics);
    
    //force_log(force);
    
            
/* 3 . update */
    if (tctype == 1)
    {
        Berendsen_update(atom_kinematics,atom_statics,
                            tcstat,tcpara,tcgrp,force,box,dt);
    }
    pbc_with_water(atom_kinematics,atom_statics,box,atom_size,molelist);
    float invdt = 1.0 / dt;
    
    //printf("%.23f,%.23f,%.23f\n",atom_kinematics.rx[1280],atom_kinematics.ry[1280],atom_kinematics.rz[1280]);
    //printf("%.23f,%.23f,%.23f\n",atom_kinematics.rx[1281],atom_kinematics.ry[1281],atom_kinematics.rz[1281]);
    //printf("%.23f,%.23f,%.23f\n",atom_kinematics.rx[1282],atom_kinematics.ry[1282],atom_kinematics.rz[1282]);
    //printf("%.23f,%.23f,%.23f\n",old_k.rx[1280],old_k.ry[1280],old_k.rz[1280]);
    //printf("%.23f,%.23f,%.23f\n",old_k.rx[1281],old_k.ry[1281],old_k.rz[1281]);
    //printf("%.23f,%.23f,%.23f\n",old_k.rx[1282],old_k.ry[1282],old_k.rz[1282]);

    //log_print_data(step,"out/beforesettle.log",atom_kinematics);
    for (int i = 0; i < molelist.getsize(); i++)
    {
        if (molelist[i].iswater)
        {
            settle(molelist[i],atom_kinematics,old_k,atom_statics,virial,invdt);
        }
    }
    //printf("lincs\n");
    //do_lincs(lincs_data,atom_statics,atom_kinematics,old_k,virial,invdt);

    if ((pcpara.do_pc != 0) && step > 1000)
    {
        Berendsen_pc_rx_scale(pcstat,atom_kinematics);
        Berendsen_pc_box_scale(pcstat,box);
    }
            
    kinematics_copy(atom_kinematics,old_k);
    //tmp
    //log_print_data(step,"out/beforeekin.log",atom_kinematics);
    ekin_sum(atom_kinematics, atom_statics, tcstat,tcpara,tcgrp);
    for (int i = 0; i < tcpara.N; i++)
    {
        tcstat.lambda[i] = Berendsen_lambda(tcpara, tcstat, dt, i);
    }
    //printf("ekin\n");
            
    //calc_vcm_grp(atom_statics,atom_kinematics,vcm);
    if (pcpara.do_pc != 0)
    {
        calc_vir(virial,force,atom_kinematics);
        calc_pres(tcstat,pcstat,virial,box);
        Berendsen_pc_mu(pcstat,pcpara);
    }
    /*
    vcm_v_calc(vcm);
    do_stopcm_grp(atom_kinematics,vcm,0);
    }*/

    //if (step % nstep.n_log == 0 || step > 1100)
    if (step % nstep.n_log == 0 )
    {
        printf("step %d\n",step);
        
        FILE *fp = NULL;
        if (fp = fopen("out/lambda.txt","a"))
        {
            fprintf(fp,"step %d\n",step);
            for (int i = 0; i < tcpara.N; i++)
            {
                fprintf(fp,"group %d\n",i);
                fprintf(fp,"lbd = %f\n",tcstat.lambda[i]);
                fprintf(fp,"Tc = %f\n",tcstat.T[i]);
                fprintf(fp,"Tch = %f\n",tcstat.Th[i]);
                fprintf(fp,"Ekin = %f\n",tcstat.ekin[i]);
            }
            
            
            fprintf(fp,"Evdw = %f\n",force.esum_vdw);
            fprintf(fp,"Eqq = %f\n",force.esum_qq);
            fprintf(fp,"Ebond = %f\n",force.esum_bond);
            fprintf(fp,"Eangle = %f\n",force.esum_angle);
            fprintf(fp,"Epdih = %f\n",force.esum_pdih);
            fprintf(fp,"Eidih = %f\n",force.esum_idih);
            fclose(fp);
        }
        if (fp = fopen("out/press.txt","a"))
        {
            fprintf(fp,"step %d\n",step);
            fprintf(fp,"p = %f\n",pcstat.pres_scalar);
            fprintf(fp,"Ekin = %f,%f,%f\n",tcstat.ekin_tot[0],tcstat.ekin_tot[1],tcstat.ekin_tot[2]);
            fprintf(fp,"Vir = %f,%f,%f\n",virial.sum[0],virial.sum[1],virial.sum[2]);
            fprintf(fp,"mu = %f,%f,%f\n",pcstat.mu[XXXX],pcstat.mu[YYYY],pcstat.mu[YYYY]);
            fprintf(fp,"box = %f,%f,%f\n",box[0],box[1],box[2]);
            fclose(fp);
        }
        if (fp = fopen("out/p.txt","a"))
        {
            fprintf(fp,"%f\n",pcstat.pres_scalar);
            fclose(fp);
        }
        if (true)
        {
            //printf("vcm:%f,%f,%f\n",vcm.vcm_x,vcm.vcm_y,vcm.vcm_z);
            //log_print_data(step,"out/md.log",atom_kinematics);
            log_print_data_2(step,"out/md_2.log",atom_kinematics,atom_statics);

            //force_log(force);
        }
        
        
    }
}

void kinematics_copy(t_atom_kinematics & src,t_atom_kinematics & dest)
{
    int size = src.rx.getsize();
    dest.ax.clear();
    dest.ay.clear();
    dest.az.clear();
    dest.vx.clear();
    dest.vy.clear();
    dest.vz.clear();
    dest.rx.clear();
    dest.ry.clear();
    dest.rz.clear();
    for (int i = 0; i < size; i++)
    {
        dest.ax.put(src.ax[i]);
        dest.ay.put(src.ay[i]);
        dest.az.put(src.az[i]);
        dest.vx.put(src.vx[i]);
        dest.vy.put(src.vy[i]);
        dest.vz.put(src.vz[i]);
        dest.rx.put(src.rx[i]);
        dest.ry.put(src.ry[i]);
        dest.rz.put(src.rz[i]);
    }
}

void csv_print_data_2(char* filename, t_atom_kinematics & atom_kinematics,
                      t_atom_statics & atom_statics,
                      t_atom_map & atom_map
                     )
{
    FILE *fp = NULL;
    int l = atom_kinematics.ax.getsize();
    if (fp = fopen(filename,"wt"))
    {
        for (int i = 0; i < l; i++)
        {
            int k = atom_statics.type[i];
            int s = atom_statics.mole_serial[i];
            fprintf(fp,"%d,%d,%s,%s,%f,%f,%f\n",
            s+1,k+1,atom_map.mole_name[k],atom_map.atom_name[k],
            atom_kinematics.rx[i],atom_kinematics.ry[i],atom_kinematics.rz[i]);
        }
        fclose(fp);
    }
}

void force_log(t_force & force)
{
    int n = force.fx.getsize();
    FILE *fp = NULL;
    fp = fopen("out/force.log","a");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp,"atom %d,fx = %f ,fy = %f ,fz = %f\n",i,force.fx[i],force.fy[i],force.fz[i]);
    }
    fclose(fp);
}

void fclear(t_force & force, t_virial & virial)
{
    int count = force.fx.getsize();
    for (int i = 0; i < count; i++)
    {
        force.fx[i] = 0.0f;
        force.fy[i] = 0.0f;
        force.fz[i] = 0.0f;
        virial.vir_x[i] = 0.0f;
        virial.vir_y[i] = 0.0f;
        virial.vir_z[i] = 0.0f;
    }
    force.esum_qq = 0.0f;
    force.esum_vdw = 0.0f;
    force.esum_bond = 0.0f;
    force.esum_angle = 0.0f;
    force.esum_pdih = 0.0f;
    force.esum_idih = 0.0f;
    
}

void init_box(float *box)
{
    int atom_size = 0;
    char strline[256];
    char filename[64];
    char tmp[256];
    bool fileNotEnd = true;
    char *tik;
    ini_get_str("dump_file","in/settings.ini",filename);
    FILE *fp = NULL;

    if(fp = fopen(filename,"r"))
    {
        while (fileNotEnd)
        {
            if (fgets(strline, MAX_STRLEN, fp)!=NULL)
            {//reaching the end/
                if (strncmp("box", strline, 3) == 0)
                {//leftvalue = key?
                    strcpy(tmp,strline);
                    tik = strtok(strline, " ");
                    strTrim(tik);
                    if (strcmp(tik,"box") == 0)
                    {
                        fgets(strline, MAX_STRLEN, fp);
                        tik = strtok(strline,"{");
                        tik = strtok(NULL,",");
                        box[XX] = atof(tik);
                        fgets(strline, MAX_STRLEN, fp);
                        tik = strtok(strline,",");
                        tik = strtok(NULL,",");
                        box[YY] = atof(tik);
                        fgets(strline, MAX_STRLEN, fp);
                        tik = strtok(strline,",");
                        tik = strtok(NULL,",");
                        tik = strtok(NULL,"}");
                        box[ZZ] = atof(tik);
                    }
                }
            }
            else
                fileNotEnd = false;
        }
        fclose(fp);
    }
}

