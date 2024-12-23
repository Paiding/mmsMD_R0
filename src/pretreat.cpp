#include "pretreat.h"



void t_mole_map::build_mole_map()
{
    char filename[64];
    char line[256];
    char *tik;
    ini_get_str("dump_file","in/settings.ini",filename);
    FILE *fp;
    if (fp = fopen(filename,"r"))
    {
        /* code */
        int k = 1;
        while (k != 0)
	    {
		    fgets(line,128,fp);
            if (strcmp(line,"topology:\n") == 0)
            {
                break;
            }
        }
        fgets(line,128,fp);//name = xx
        fgets(line,128,fp);//#atoms
        fgets(line,128,fp);//#molblock
        tik = strtok(line,"=");
        tik = strtok(NULL,"=");
        moletypes = atoi(tik);
        for (int i = 0; i < moletypes; i++)
        {
            fgets(line,128,fp);//molblock
            fgets(line,128,fp);//moltype              = 0 "Protein_chain_A"
            tik = strtok(line,"\"");
            tik = strtok(NULL,"\"");
            mole_type.strput(tik);
            fgets(line,128,fp);
            tik = strtok(line,"=");
            tik = strtok(NULL,"=");
            mole_num.put(atoi(tik));
            fullmolenum += atoi(tik);
            fgets(line,128,fp);
            fgets(line,128,fp);
        }
        fclose(fp);
    }
    int count = 0;
    int excl_count = 0;
    int excl_offset = 0;
    for (int i = 0; i < moletypes; i++)
    {
        char num[8];
        int tmp_i;
        char dest1[32];
        char dest2[32];
        angle_num.put(0);
        pdih_num.put(0);
        ipdih_num.put(0);
        lj14_num.put(0);
        cons_num.put(0);
        int il;
        int ir;
        
        atom_loc.put(count);
        //printf("1\n");
        if (fp = fopen(filename, "r"))
        {//read bond parameter
            int k = 1;
            int bondNum = 0;
            
            //atoms
            while (k != 0)
		    {
			    fgets(line,128,fp);
                //sscanf :moltype (0)  dest1 for "moltype" and dest2 for "(0)"
                sscanf(line,"%s %s",dest1,dest2);
                if (strcmp(dest1,"moltype") == 0 && strcmp(dest2,"=")!=0)
                {
                    int typeSer = collect_int(line);
                    
                    if (typeSer == i)
                    {
                        break;
                    }
                }
            }

            //arrived the target paragraph
            fgets(line,128,fp);//      name="Protein_chain_A"
            
            tik = strtok(line,"\"");
            tik = strtok(NULL,"\"");
            this->mole_type.strput(tik);
            fgets(line,128,fp);//skip the title line "atoms:"
            fgets(line,128,fp);//atom (23):
            int atomNum = collect_int(line);
            atom_num.put(atomNum);
            for (int j = 0; j < atomNum; j++)
            {
                //atom[     0]={type=  0, typeB=  0, ptype=    Atom, m= 1.40070e+01,
                // q=-9.60000e-01, mB= 1.40070e+01, qB=-9.60000e-01, resind=    0, atomnumber=  7}
                fgets(line,256,fp);
                tik = strtok(line,"=");//atom[     0]
                tik = strtok(NULL,"=");//{type
                tik = strtok(NULL,",");//  0
                atom_type.put(atoi(tik));
                tik = strtok(NULL,"=");// typeB
	            tik = strtok(NULL,"=");//  0, ptype
	            tik = strtok(NULL,"=");//    Atom, m
	            tik = strtok(NULL,",");// 1.40070e+01
	            atom_mass.put(atof(tik));
                tik = strtok(NULL,"=");// q
	            tik = strtok(NULL,",");//-9.60000e-01
                atom_q.put(atof(tik));
                tik = strtok(NULL,"=");// mB
	            tik = strtok(NULL,"=");// 1.40070e+01, qB
	            tik = strtok(NULL,"=");//-9.60000e-01, resind
	            tik = strtok(NULL,",");//    0
                atom_resind.put(atoi(tik));
                count ++;
            }
            
            //excls:
            while (k != 0)
		    {
			    fgets(line,128,fp);
                
                sscanf(line,"%s",dest1);
                if (strcmp(dest1,"excls:") == 0)
                {
                    break;
                }
            }
            fgets(line,128,fp);//nr = 23
            tik = strtok(line,"=");
            tik = strtok(NULL,"=");
            int nr = atoi(tik);
        
            excl_idloc.put(excl_count);
            excl_lsloc.put(excl_offset);
            excl_count += (nr+1);
            fgets(line,128,fp);//nra=237
            for (int ii = 0; ii < nr; ii++)
            {
                fgets(line,128,fp);
                //printf("%s",line);
                //excls[0][0..11]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
                tik = strtok(line,"[");//excls
	            tik = strtok(NULL,"]");//0
                
                tik = strtok(NULL,"[");//NULL,tik eat the right half of string
                //tik = 0..11]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
	            tik = strtok(tik,".");//0
                il = atoi(tik);
                exclist_index.put(il);
                tik = strtok(NULL,".");//NULL,tik ="11]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}"
	            tik = strtok(tik,"]");//11
                ir = atoi(tik);
                tik = strtok(NULL,"{");
                for (int j = il; j < ir+1; j++)
                {
                    tik = strtok(NULL,",");
                    if (strcmp(tik," \n") == 0)
                    {
                        fgets(line,128,fp);
                        if (strchr(line,',')!=NULL)
                        {
                            tik = strtok(line,",");
                        }
                        else
                        {//only one element in the next line
                            tik = strtok(line,"}");
                        }
                    }
                    //printf("%s,",tik);
                    exclist.put(atoi(tik));
                    excl_offset ++;
                }
            }
            exclist_index.put(ir);
            
             //G96Angle
            while (k != 0)
		    {
			    fgets(line,128,fp);
				
                sscanf(line,"%s",dest1);
                if (strcmp(dest1,"G96Angle:") == 0)
                {
                    break;
                }
            }
			
			int atomi,atomj,atomk;
            //arrived the target paragraph
			
            fgets(line,128,fp);//      nr: 148
            
			sscanf(line,"%s %d",dest1,&nr);
			nr = nr / 4;
			fgets(line,128,fp);//         iatoms:
			for (int j = 0; j < nr ; j++)
            {
                //            0 type=169 (UREY_BRADLEY)   1   0   2
				fgets(line,128,fp);
				sscanf(line,"%d %s %s %d %d %d",&tmp_i,dest1,dest2,&atomi,&atomj,&atomk);
                tik = strtok(dest1,"=");
                tik = strtok(NULL,"=");
                funcType_a.put(atoi(tik));
                angleType.put(2);//type 2 g96 angle
                angle_num[i] += 1;
                anglei.put(atomi);
                anglej.put(atomj);
                anglek.put(atomk);
                anglea0.put(0);
                angleaK.put(0);
            }

            //U-B
            while (k != 0)
		    {
			    fgets(line,128,fp);
				
                sscanf(line,"%s",dest1);
                if (strcmp(dest1,"U-B:") == 0)
                {
                    break;
                }
            }
			
            //arrived the target paragraph
			
            fgets(line,128,fp);//      nr: 148
            
			sscanf(line,"%s %d",dest1,&nr);
			nr = nr / 4;
			fgets(line,128,fp);//         iatoms:
			for (int j = 0; j < nr ; j++)
            {
                //            0 type=169 (UREY_BRADLEY)   1   0   2
				fgets(line,128,fp);
				sscanf(line,"%d %s %s %d %d %d",&tmp_i,dest1,dest2,&atomi,&atomj,&atomk);
                tik = strtok(dest1,"=");
                tik = strtok(NULL,"=");
                funcType_a.put(atoi(tik));
                angleType.put(1);
                angle_num[i] += 1;
                anglei.put(atomi);
                anglej.put(atomj);
                anglek.put(atomk);
                anglea0.put(0);
                angleaK.put(0);
            }

            //proper dihedrals
            while (k != 0)
		    {
			    fgets(line,128,fp);
                sscanf(line,"%s",dest1);
                if (strcmp(dest1,"Proper") == 0)
                {
                    break;
                }
            }
			
			int atoml;
            //arrived the target paragraph
			
            fgets(line,128,fp);//      nr: 155
			sscanf(line,"%s %d",dest1,&nr);
			nr /= 5;
            pdih_num[i] += nr;
			fgets(line,128,fp);//         iatoms:
            for (int j = 0; j < nr; j++)
            {
                //0 type=191 (PDIHS)   1   0   3   4
				fgets(line,128,fp);
				sscanf(line,"%d %s %s %d %d %d %d",&tmp_i,dest1,dest2,&atomi,&atomj,&atomk,&atoml);
				
				tik = strtok(dest1,"=");
                tik = strtok(NULL,"=");
				funcType_pdih.put(atoi(tik));
                pdihi.put(atomi);
                pdihj.put(atomj);
                pdihk.put(atomk);
                pdihl.put(atoml);
            }

            //improper dihedrals
            while (k != 0)
		    {
			    fgets(line,128,fp);
                sscanf(line,"%s",dest1);
                if (strcmp(dest1,"Improper") == 0)
                {
                    break;
                }
            }
			
            //arrived the target paragraph
			
            fgets(line,128,fp);//      nr: 155
			sscanf(line,"%s %d",dest1,&nr);
			nr = nr / 5;
            ipdih_num[i] += nr;
			fgets(line,128,fp);//         iatoms:
            for (int j = 0; j < nr; j++)
            {
                //0 type=191 (PDIHS)   1   0   3   4
				fgets(line,128,fp);
				sscanf(line,"%d %s %s %d %d %d %d",&tmp_i,dest1,dest2,&atomi,&atomj,&atomk,&atoml);
				
				tik = strtok(dest1,"=");
                tik = strtok(NULL,"=");
				funcType_ipdih.put(atoi(tik));
                ipdihi.put(atomi);
                ipdihj.put(atomj);
                ipdihk.put(atomk);
                ipdihl.put(atoml);
            }
            //LJ14
            while (k != 0)
		    {
			    fgets(line,128,fp);
                sscanf(line,"%s",dest1);
                if (strcmp(dest1,"LJ-14:") == 0)
                {
                    break;
                }
            }
			
            //arrived the target paragraph
			
            fgets(line,128,fp);//      nr: 155
			sscanf(line,"%s %d",dest1,&nr);
			nr = nr / 3;
            lj14_num[i] += nr;
			fgets(line,128,fp);//         iatoms:
            for (int j = 0; j < nr; j++)
            {
                //0 type=203 (LJ14)   0   6
				fgets(line,128,fp);
				sscanf(line,"%d %s %s %d %d",&tmp_i,dest1,dest2,&atomi,&atomj);
				
				tik = strtok(dest1,"=");
                tik = strtok(NULL,"=");
				funcType_lj14.put(atoi(tik));
                lj14i.put(atomi);
                lj14j.put(atomj);
            }

            //CONSTR
            while (k != 0)
		    {
			    fgets(line,128,fp);
                sscanf(line,"%s",dest1);
                if (strcmp(dest1,"Constraint:") == 0)
                {
                    break;
                }
            }
			
            //arrived the target paragraph
			
            fgets(line,128,fp);//      nr: 155
			sscanf(line,"%s %d",dest1,&nr);
			nr = nr / 3;
            cons_num[i] += nr;
			fgets(line,128,fp);//         iatoms:
            for (int j = 0; j < nr; j++)
            {
                //0 type=225 (CONSTR)   0   1
				fgets(line,128,fp);
				sscanf(line,"%d %s %s %d %d",&tmp_i,dest1,dest2,&atomi,&atomj);
				
				tik = strtok(dest1,"=");
                tik = strtok(NULL,"=");
				funcType_constr.put(atoi(tik));
                constri.put(atomi);
                constrj.put(atomj);
            }

            fclose(fp);   
		}
    }
    excl_idloc.put(excl_count);
    excl_lsloc.put(excl_offset);
}

int t_mole_map::find_loc_byname(char* mole_name_key)
{
    for (int i = 0; i < moletypes; i++)
    {
        strTrim(mole_name_key);
        if (strcmp(mole_name_key,mole_type[i]) == 0)
        {
            return i;
        }
    }
    return -1;
}

void read_functype(array<float> & lj_para, t_ub_para & ub_para, t_g96a_para & g96a_para,
                   t_pdih_para & pdih_para, t_ipdih_para & ipdih_para,
                   t_lj14_para & lj14_para, t_constr_para & constr_para,
                   t_atom_statics & atom_statics
                  )
{
    char dest1[256];
    char tmp[256];
    char dest2[32];
    char filename[64];
    int ntypes = 0;
    ini_get_str("dump_file","in/settings.ini",filename);

    FILE *fp = NULL;
    if(fp = fopen(filename,"r"))
    {
        

        int k = 1;
        
        char line[256];
    
        while (k != 0)
	    {
		    fgets(line,128,fp);
                
            sscanf(line,"%s",dest1);
            if (strcmp(dest1,"ffparams:") == 0)
            {
                break;
            }
        }
        fgets(line,128,fp);//atnr=13
        char *tik = strtok(line,"=");
        tik = strtok(NULL,"=");
        atom_statics.atnr = atoi(tik);
        fgets(line,128,fp);//ntypes=240
        
        tik = strtok(line,"=");
        tik = strtok(NULL,"=");
        ntypes = atoi(tik);
        for (int i = 0; i < ntypes; i++)
        {
            fgets(line,256,fp);
            strcpy(tmp,line);
            tik = strtok(line,"=");
            tik = strtok(NULL,",");
            if (strcmp(tik,"LJ_SR") == 0)
            {
            //functype[0]=LJ_SR, c6= 4.29399917e-03, c12= 5.50861296e-06
                tik = strtok(tmp,"=");//functype[0]
	            tik = strtok(NULL,"=");//LJ_SR, c6
                tik = strtok(NULL,",");//4.29399917e-03
                lj_para.put(atof(tik));
		        tik = strtok(NULL,"=");// c12
		        tik = strtok(NULL,",");//5.50861296e-06
                lj_para.put(atof(tik));
            }
            else if (strcmp(tik,"UREY_BRADLEY") == 0)
            {
            //functype[179]=UREY_BRADLEY, thetaA= 1.21000000e+02
            //, kthetaA= 6.69440002e+02, r13A= 0.00000000e+00, kUBA= 0.00000000e+00
                ub_para.type_index.put(i);
                tik = strtok(tmp,"=");
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
                ub_para.theta0.put(atof(tik));
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
                ub_para.ktheta.put(atof(tik));
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
                ub_para.r130.put(atof(tik));
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
                ub_para.kUB.put(atof(tik));
            }
            else if (strcmp(tik,"PDIHS") == 0)
            {
                pdih_para.type_index.put(i);
                tik = strtok(tmp,"=");
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");//phiA
				pdih_para.phi0.put(atof(tik));
                
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");//cpA
                pdih_para.cp0.put(atof(tik));
                
				tik = strtok(NULL,"=");
				//phiB
				tik = strtok(NULL,"=");
				//cpB
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
                pdih_para.mult.put(atoi(tik));
            }
            else if (strcmp(tik,"IDIHS") == 0)
            {
                ipdih_para.type_index.put(i);
                tik = strtok(tmp,"=");
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");//xiA
				ipdih_para.xi0.put(atof(tik));
                
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");//cxA
                ipdih_para.cx0.put(atof(tik));

            }
            else if (strcmp(tik,"LJ14") == 0)
            {
                lj14_para.type_index.put(i);
                tik = strtok(tmp,"=");
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
				lj14_para.c6.put(atof(tik));
                
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
                lj14_para.c12.put(atof(tik));
            }
            else if (strcmp(tik,"CONSTR") == 0)
            {
                constr_para.type_index.put(i);
                tik = strtok(tmp,"=");
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
				constr_para.d0.put(atof(tik));
            }
            else if (strcmp(tik,"G96ANGLES") == 0)
            {
                g96a_para.type_index.put(i);
                tik = strtok(tmp,"=");
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
				g96a_para.theta0.put(atof(tik));
                
				tik = strtok(NULL,"=");
				tik = strtok(NULL,",");
                g96a_para.ct.put(atof(tik));
            }
            else if (strcmp(tik,"SETTLE") == 0)
            {
                /* code */
                //nothing happened
            }
        }
        
        fclose(fp);
    }
    
}

void gro2csv(char* filename, char* savefile)
{
    char line[GRO_STRLEN];
    int linenum = 0;
    //infomation to be read from .gro file
    char **moleserial;
    char **molename;
    char **atomtype;
    char **atomserial;
    char **rx;
    char **ry;
    char **rz;

    FILE *fp = NULL;
    if (fp = fopen(filename, "r"))
    {
        fgets(line,GRO_STRLEN,fp);//first line: infomation
        fgets(line,GRO_STRLEN,fp);//second line: number of atoms
        linenum = atoi(line);
        printf("linenum = %d\n",linenum);
        moleserial = (char**)malloc(sizeof(char*)*linenum);
        for (int i = 0; i < linenum; i++) moleserial[i] = (char*)malloc(sizeof(char)*5);
        molename = (char**)malloc(sizeof(char*)*linenum);
        for (int i = 0; i < linenum; i++) molename[i] = (char*)malloc(sizeof(char)*5);
        atomtype = (char**)malloc(sizeof(char*)*linenum);
        for (int i = 0; i < linenum; i++) atomtype[i] = (char*)malloc(sizeof(char)*5);
        atomserial = (char**)malloc(sizeof(char*)*linenum);
        for (int i = 0; i < linenum; i++) atomserial[i] = (char*)malloc(sizeof(char)*5);
        rx = (char**)malloc(sizeof(char*)*linenum);
        for (int i = 0; i < linenum; i++) rx[i] = (char*)malloc(sizeof(char)*8);
        ry = (char**)malloc(sizeof(char*)*linenum);
        for (int i = 0; i < linenum; i++) ry[i] = (char*)malloc(sizeof(char)*8);
        rz = (char**)malloc(sizeof(char*)*linenum);
        for (int i = 0; i < linenum; i++) rz[i] = (char*)malloc(sizeof(char)*8);

        for (int i = 0; i < linenum; i++)
        {
            fgets(line,GRO_STRLEN,fp);
            int pos=0;
		// resid
		    strncpy(moleserial[i], line, 5);
		    moleserial[i][5] = '\0';
		    pos += 5;
		// resname
		    strncpy(molename[i], line+pos, 5);
		    molename[i][5] = '\0';
		    pos += 5;
		// atomname
		    strncpy(atomtype[i], line+pos, 5);
		    atomtype[i][5] = '\0';
		    pos += 5;
		// atomid
		    strncpy(atomserial[i], line+pos, 5);
		    atomserial[i][5] = '\0';
		    pos += 5;
		// pos
		// x
		    strncpy(rx[i], line+pos, 8);
		    rx[i][8] = '\0';
		    pos += 8;
		// y
		    strncpy(ry[i], line+pos, 8);
		    ry[i][8] = '\0';
		    pos += 8;
		// z
		    strncpy(rz[i], line+pos, 8);
		    rz[i][8] = '\0';
		    pos += 8;
        }
        fclose(fp);
    }
    else
    {
        perror("can't open gro");
    }

    if (fp = fopen(savefile, "w+"))
    {
        int k = 1;
        int i = 0;
        fprintf(fp,"%s,%d,%s,%s,%s,%s,%s\n",moleserial[i],k,molename[i],atomtype[i],rx[i],ry[i],rz[i]);//first line special
        for (i = 1; i < linenum; i++)
        {
            if (strcmp(moleserial[i],moleserial[i - 1]) == 0)//from 2nd line on
            {//add the inmole atom serial if it is still in the same molecular
                k += 1;
            }
            else
            {
                k = 1;
            }
            fprintf(fp,"%s,%d,%s,%s,%s,%s,%s\n",moleserial[i],k,molename[i],atomtype[i],rx[i],ry[i],rz[i]);
        }
        fclose(fp);
    }
    else
    {
        perror("can't open csv");
    }
    for (int i = 0; i < linenum; i++)
    {
        free(moleserial[i]);
        free(molename[i]);
        free(atomtype[i]);
        free(atomserial[i]);
        free(rx[i]);
        free(ry[i]);
        free(rz[i]);
    }
    free(moleserial);
    free(molename);
    free(atomtype);
    free(atomserial);
    free(rx);
    free(ry);
    free(rz);
    
}

void csv2gro(char* read_file,
               char* save_file
              )
{
    //array<int> mol_serial;
    //array<char*> mol_name;
    //array<>
    FILE *fin = fopen(read_file,"r");
    FILE *fout = fopen(save_file, "w");
    fprintf(fout,"generated gro\n");
    fprintf(fout,"0\n");
    char strline[MAX_STRLEN];
    bool fileNotEnd = true;
    int line_loc = 0;
    while (fileNotEnd)
    {
        if (fgets(strline, MAX_STRLEN, fin)!=NULL)
        {//read one line
            line_loc ++;
            //molecular serial
            char *tik = strtok(strline, ",");
            int moleSer = atoi(tik);//index of mole_array is serial-1
                
            //atom serial
            tik = strtok(NULL, ",");
            
            //molecular name
            tik = strtok(NULL, ",");
            char mole_name[8];
            strTrim(tik);
            strcpy(mole_name,tik);
                
            //atom type 
            tik = strtok(NULL, ",");
            char * atomtype = new char[8];
            strcpy(atomtype,tik);
            atomtype = strLTrim(atomtype);
                
            //rx
            tik = strtok(NULL, ",");
            float rx = atof(tik);

            //ry
            tik = strtok(NULL, ",");
            float ry = atof(tik);

            //rz
            tik = strtok(NULL, ",");
            float rz = atof(tik);

            fprintf(fout, "%5d%-5.5s%5.5s%5d%8.3f%8.3f%8.3f\n", 
						moleSer, mole_name, atomtype, line_loc,
						rx,ry,rz);
        }
        else
        {
            fileNotEnd = false;
        }
	}
}

