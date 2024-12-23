#include "nblist.h"

#define TESTRCUT 1.5

void build_excl_list(array<t_molecule> & molelist,
                    t_mole_map & mole_map,
                     t_excl_list & excl_list,
                     int atom_size
                    )
{
/* refresh when molecule added or loss*/
//clear bonded list
    if (excl_list.bonded_ranks != NULL)
    {
        for (int i = 0; i < excl_list.sizei; i++)
        {
            delete[] excl_list.bonded_ranks[i];
        }
        excl_list.bonded_ranks = NULL;//??
    }
    excl_list.bonded_ranks = new int*[atom_size];
    excl_list.sizei = atom_size;
    for (int i = 0; i < atom_size; i++)
    {
        excl_list.bonded_ranks[i] = new int[EXCL_MAX];
        for (int j = 0; j < EXCL_MAX; j++)
        {
            excl_list.bonded_ranks[i][j] = -1;
        }
        
    }
    
    int atomi;
    int atomj;
    
    int *flag = new int[atom_size]{0};
    for (int i = 0; i < molelist.getsize(); i++)
    {   //in one molecule
        int type = molelist[i].mole_map_rank;
        int jl = mole_map.excl_idloc[type];
        int jr = mole_map.excl_idloc[type+1];
        int offset = mole_map.excl_lsloc[type];
        for (int j = jl; j < jr; j++)
        {
            int kl = mole_map.exclist_index[j];
            int kr = mole_map.exclist_index[j+1];

            for (int k = kl; k < kr; k++)
            {

                atomi = molelist[i].atom_ranks[j - jl];//count from 0 to jr-jl
                atomj = molelist[i].atom_ranks[mole_map.exclist[k+offset]];
                if (atomi>atomj)
                {
                    excl_list.bonded_ranks[atomi][flag[atomi]++] = atomj;
                    excl_list.bonded_ranks[atomj][flag[atomj]++] = atomi;
                }
                
                
                
            }
            
            
            //printf("[%d,%d]\n",atomi,atomj);
        }
    }
    delete[] flag;
}



void build_nblist_img(float *box, 
                        t_nblist & nblist,
                        int atom_size,const t_radius & radius,
                        t_atom_kinematics & atom_kinematics,
                        t_excl_list & excl_list)
{
    float bx = box[XX];
    float by = box[YY];
    float bz = box[ZZ];
    float bxinv = 1.0f / bx;
	float byinv = 1.0f / by;
	float bzinv = 1.0f / bz;
    float rvdw2 = radius.r_vdw * radius.r_vdw;
    float dr2;
    int k;

    if (nblist.va_list != NULL)
    {
        for (int i = 0; i < nblist.size; i++)
        {
            delete[] nblist.va_list[i];
        }
        nblist.va_list = NULL;//??
    }
    nblist.va_list = new int*[atom_size];
    nblist.size = atom_size;
    for (int i = 0; i < atom_size; i++)
    {
        nblist.va_list[i] = new int[atom_size];
        for (int j = 0; j < atom_size; j++)
        {
            nblist.va_list[i][j] = -1;
        }
        
    }
    
    for (int i = 0; i < atom_size; i++)
    {   
        for (int j = i + 1; j < atom_size; j++)//check half of atoms
        {
            float dx = 0;
            float dy = 0;
            float dz = 0;
            int tx = 1;
            int ty = 1;
            int tz = 1;

            float xi = atom_kinematics.rx[i];
            float yi = atom_kinematics.ry[i];
            float zi = atom_kinematics.rz[i];
            float xj = atom_kinematics.rx[j];
            float yj = atom_kinematics.ry[j];
            float zj = atom_kinematics.rz[j];

            for (int t = 0; t < 27; t++)
            {
                tx = t / 9;
                ty = (t % 9) / 3;
                tz = t % 3;
                tx -= 1;
                ty -= 1;
                tz -= 1;
                    
                dx = xi - xj - tx * bx;
                dy = yi - yj - ty * by;
                dz = zi - zj - tz * bz;
                dr2 = dx*dx + dy*dy + dz*dz;
                if (dr2 < rvdw2)
                {
                        nblist.va_list[i][j] = t;
                        nblist.va_list[j][i] = 26 - t;
                        continue;
                }
            }
        }
    }
    
    //delete atoms in excl list
    
    for (int i = 0; i < atom_size; i++)
    {   
        int j = 0;
        while (excl_list.bonded_ranks[i][j] != -1)
        {
            k = excl_list.bonded_ranks[i][j];
            nblist.va_list[i][k] = -1;
            nblist.va_list[k][i] = -1;
            //printf("%d,%d\n",i,k);
            j ++;
        }
    }
    
}
void build_nblist_bin(float *box, 
                      t_nblist & nblist,
                      int atom_size,const t_radius & radius,
                      t_atom_kinematics & atom_kinematics,
                      t_excl_list & excl_list)
{
    float bx = box[XX];
    float by = box[YY];
    float bz = box[ZZ];
    float rx,ry,rz;
    float bxinv = 1.0f / bx;
	float byinv = 1.0f / by;
	float bzinv = 1.0f / bz;
    float rvdw = radius.r_vdw;
    float rvdw2 = radius.r_vdw * radius.r_vdw;
    float dr2;
    int k;

    if (nblist.va_list != NULL)
    {
        for (int i = 0; i < nblist.size; i++)
        {
            delete[] nblist.va_list[i];
        }
        nblist.va_list = NULL;//??
    }
    nblist.va_list = new int*[atom_size];
    nblist.size = atom_size;
    for (int i = 0; i < atom_size; i++)
    {
        nblist.va_list[i] = new int[atom_size];
        for (int j = 0; j < atom_size; j++)
        {
            nblist.va_list[i][j] = -1;
        }
    }
    /*bin_list */
    int nbinx,nbiny,nbinz;//number of bins at each coordinate
    float dbinx,dbiny,dbinz;
    int ibinx,ibiny,ibinz;
    nbinx = int(floor(bx/rvdw));
    nbiny = int(floor(by/rvdw));
    nbinz = int(floor(bz/rvdw));
    int nbin = (nbinx+2) * (nbiny+2) * (nbinz+2);
    printf("nbin = %d\n",nbin);
    
    dbinx = bx / nbinx;
    dbiny = by / nbiny;
    dbinz = bz / nbinz;
    if (nblist.bin_list != NULL)
    {
        for (int i = 0; i < nbin; i++)
        {
            delete[] nblist.bin_list[i];
        }
        nblist.bin_list = NULL;//??
        //delete[] nblist.bin_list;
    }
    nblist.bin_list = new int*[nbin];
    for (int i = 0; i < nbin; i++)
    {
        nblist.bin_list[i] = new int[atom_size + 2];
        nblist.bin_list[i][0] = 0;//used as the number of atoms in a bin
    }
    printf("%d\n",atom_size+1);
    int loc,ser;
    //cycle the atomlist to alloc bin_list
    printf("%f,%f,%f\n",atom_kinematics.rx[0],atom_kinematics.ry[0],atom_kinematics.rz[0]);
    printf("%f,%f,%f\n",dbinx,dbiny,dbinz);
    for (int i = 0; i < atom_size; i++)
    {
        rx = atom_kinematics.rx[i];
        ry = atom_kinematics.ry[i];
        rz = atom_kinematics.rz[i];
        ibinx = int(floor(rx/dbinx)) + 1;
        ibiny = int(floor(ry/dbiny)) + 1;
        ibinz = int(floor(rz/dbinz)) + 1;
        
        loc = nbinx * nbiny * ibinz + nbinx * ibiny + ibinx;
        //printf("(%d)%f,%f,\n",i,ry,dbiny);
        //printf("(%d)%d,%d,\n",i,nbinx,ibiny);
        //printf("(%d)[%d][%d][%d](%d,\n",i,nbinx * nbiny * ibinz,nbinx * ibiny,ibinx,loc);
        nblist.bin_list[loc][0] ++;
        ser = nblist.bin_list[loc][0];
        //printf("%d)",ser);
        
        nblist.bin_list[loc][ser] = i;
    }
    
    /*end bin_list*/
    /*search*/
    printf("search begin\n\n");
    int *nbbinx = new int[nbinx+2];
    int *nbbiny = new int[nbiny+2];
    int *nbbinz = new int[nbinz+2];
    nbbinx[0] = nbinx;
    nbbiny[0] = nbiny;
    nbbinz[0] = nbinz;
    for (int i = 1; i < nbinx + 1; i++)
    {
        nbbinx[i] = i;
    }
    for (int i = 1; i < nbiny + 1; i++)
    {
        nbbiny[i] = i;
    }
    for (int i = 1; i < nbinz + 1; i++)
    {
        nbbinz[i] = i;
    }
    nbbinx[nbinx+1] = 0; 
    nbbiny[nbiny+1] = 0; 
    nbbinz[nbinz+1] = 0;

    for (int i = 0; i < nbin; i++)
    {
        ibinz = i / (nbinx * nbiny);
        ibiny = (i % (nbinx * nbiny)) / nbinx;
        ibinx = i % nbinx;
        printf("bin%d\n",i);
        for (int si = 1; si < nblist.bin_list[i][0]; si++)
        {
            int atomi = nblist.bin_list[i][si];
            //printf("atom %d\n",atomi);
            for (int ii = 0; ii < 3; ii++)
            {
                for (int jj = 0; jj < 3; jj++)
                {
                    for (int kk = 0; kk < 3; kk++)
                    {
                        int jbinz = nbbinz[ibinz + ii];
                        int jbiny = nbbiny[ibiny + jj];
                        int jbinx = nbbinx[ibinx + kk];
                        loc = nbinx * nbiny * jbinz + nbinx * jbiny + jbinx;
                        for (int sj = 1; sj < nblist.bin_list[loc][0]; sj++)
                        {
                            int atomj = nblist.bin_list[loc][sj];
                            if (atomi == atomj)
                            {
                                continue;
                            }
                            
                            float dx = 0;
                            float dy = 0;
                            float dz = 0;
                            int tx = 1;
                            int ty = 1;
                            int tz = 1;
                        

                            float xi = atom_kinematics.rx[atomi];
                            float yi = atom_kinematics.ry[atomi];
                            float zi = atom_kinematics.rz[atomi];
                            float xj = atom_kinematics.rx[atomj];
                            float yj = atom_kinematics.ry[atomj];
                            float zj = atom_kinematics.rz[atomj];

                            for (int t = 0; t < 27; t++)
                            {
                                tx = t / 9;
                                ty = (t % 9) / 3;
                                tz = t % 3;
                                tx -= 1;
                                ty -= 1;
                                tz -= 1;
                    
                                dx = xi - xj - tx * bx;
                                dy = yi - yj - ty * by;
                                dz = zi - zj - tz * bz;
                                dr2 = dx*dx + dy*dy + dz*dz;
                                if (dr2 < rvdw2)
                                {
                                    //printf("%d,%d\n",atomi,atomj);
                                    nblist.va_list[atomi][atomj] = t;
                                    nblist.va_list[atomj][atomi] = 26 - t;
                                    continue;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] nbbinx;
    delete[] nbbiny;
    delete[] nbbinz;
    //delete atoms in excl list
    
    for (int i = 0; i < atom_size; i++)
    {   
        int j = 0;
        while (excl_list.bonded_ranks[i][j] != -1)
        {
            k = excl_list.bonded_ranks[i][j];
            nblist.va_list[i][k] = -1;
            nblist.va_list[k][i] = -1;
            //printf("%d,%d\n",i,k);
            j ++;
        }
    }
    
}
void build_nblist_bin_0(float *box, 
                      t_nblist & nblist,
                      int atom_size,const t_radius & radius,
                      t_atom_kinematics & atom_kinematics,
                      t_excl_list & excl_list)
{
    float bx = box[XX];
    float by = box[YY];
    float bz = box[ZZ];
    float rx,ry,rz;
    float bxinv = 1.0f / bx;
	float byinv = 1.0f / by;
	float bzinv = 1.0f / bz;
    float rvdw = radius.r_vdw;
    float rvdw2 = radius.r_vdw * radius.r_vdw;
    float dr2;
    int k;

    if (nblist.va_list != NULL)
    {
        for (int i = 0; i < nblist.size; i++)
        {
            delete[] nblist.va_list[i];
        }
        nblist.va_list = NULL;//??
    }
    nblist.va_list = new int*[atom_size];
    nblist.size = atom_size;
    for (int i = 0; i < atom_size; i++)
    {
        nblist.va_list[i] = new int[atom_size];
        for (int j = 0; j < atom_size; j++)
        {
            nblist.va_list[i][j] = -1;
        }
    }
    /*bin_list */
    int nbinx,nbiny,nbinz;//number of bins at each coordinate
    float dbinx,dbiny,dbinz;
    int ibinx,ibiny,ibinz;
    nbinx = int(floor(bx/rvdw));
    nbiny = int(floor(by/rvdw));
    nbinz = int(floor(bz/rvdw));
    int nbin = nbinx * nbiny * nbinz;
    printf("nbin = %d\n",nbin);
    
    dbinx = bx / nbinx;
    dbiny = by / nbiny;
    dbinz = bz / nbinz;
    if (nblist.bin_list != NULL)
    {
        for (int i = 0; i < nbin; i++)
        {
            delete[] nblist.bin_list[i];
        }
        nblist.bin_list = NULL;//??
        //delete[] nblist.bin_list;
    }
    nblist.bin_list = new int*[nbin];
    for (int i = 0; i < nbin; i++)
    {
        nblist.bin_list[i] = new int[atom_size + 2];
        nblist.bin_list[i][0] = 0;//used as the number of atoms in a bin
    }
    printf("%d\n",atom_size+1);
    int loc,ser;
    //cycle the atomlist to alloc bin_list
    printf("%f,%f,%f\n",atom_kinematics.rx[0],atom_kinematics.ry[0],atom_kinematics.rz[0]);
    printf("%f,%f,%f\n",dbinx,dbiny,dbinz);
    for (int i = 0; i < atom_size; i++)
    {
        rx = atom_kinematics.rx[i];
        ry = atom_kinematics.ry[i];
        rz = atom_kinematics.rz[i];
        ibinx = int(floor(rx/dbinx));
        ibiny = int(floor(ry/dbiny));
        ibinz = int(floor(rz/dbinz));
        
        loc = nbinx * nbiny * ibinz + nbinx * ibiny + ibinx;
        printf("(%d)%f,%f,\n",i,ry,dbiny);
        printf("(%d)%d,%d,\n",i,nbinx,ibiny);
        printf("(%d)[%d][%d][%d](%d,\n",i,nbinx * nbiny * ibinz,nbinx * ibiny,ibinx,loc);
        nblist.bin_list[loc][0] ++;
        ser = nblist.bin_list[loc][0];
        //printf("%d)",ser);
        
        nblist.bin_list[loc][ser] = i;
    }
    
    /*end bin_list*/
    /*search*/
    int *nbbinx = new int[nbinx+2];
    int *nbbiny = new int[nbiny+2];
    int *nbbinz = new int[nbinz+2];
    nbbinx[0] = nbinx - 1;
    nbbiny[0] = nbiny - 1;
    nbbinz[0] = nbinz - 1;
    for (int i = 1; i < nbinx + 1; i++)
    {
        nbbinx[i] = i - 1;
    }
    for (int i = 1; i < nbiny + 1; i++)
    {
        nbbiny[i] = i - 1;
    }
    for (int i = 1; i < nbinz + 1; i++)
    {
        nbbinz[i] = i - 1;
    }
    nbbinx[nbinx+1] = 0; 
    nbbiny[nbiny+1] = 0; 
    nbbinz[nbinz+1] = 0;

    for (int i = 0; i < nbin; i++)
    {
        ibinz = i / (nbinx * nbiny);
        ibiny = (i % (nbinx * nbiny)) / nbinx;
        ibinx = i % nbinx;
        for (int si = 0; si < nblist.bin_list[i][0]; si++)
        {
            int atomi = nblist.bin_list[i][si];
            for (int ii = 0; ii < 3; ii++)
            {
                for (int jj = 0; jj < 3; jj++)
                {
                    for (int kk = 0; kk < 3; kk++)
                    {
                        int jbinz = nbbinz[ibinz + ii];
                        int jbiny = nbbiny[ibiny + jj];
                        int jbinx = nbbinx[ibinx + kk];
                        loc = nbinx * nbiny * jbinz + nbinx * jbiny + jbinx;
                        for (int sj = 1; sj < nblist.bin_list[loc][0]; sj++)
                        {
                            int atomj = nblist.bin_list[loc][sj];
                            if (atomi == atomj)
                            {
                                continue;
                            }
                            
                            float dx = 0;
                            float dy = 0;
                            float dz = 0;
                            int tx = 1;
                            int ty = 1;
                            int tz = 1;
                        

                            float xi = atom_kinematics.rx[atomi];
                            float yi = atom_kinematics.ry[atomi];
                            float zi = atom_kinematics.rz[atomi];
                            float xj = atom_kinematics.rx[atomj];
                            float yj = atom_kinematics.ry[atomj];
                            float zj = atom_kinematics.rz[atomj];

                            for (int t = 0; t < 27; t++)
                            {
                                tx = t / 9;
                                ty = (t % 9) / 3;
                                tz = t % 3;
                                tx -= 1;
                                ty -= 1;
                                tz -= 1;
                    
                                dx = xi - xj - tx * bx;
                                dy = yi - yj - ty * by;
                                dz = zi - zj - tz * bz;
                                dr2 = dx*dx + dy*dy + dz*dz;
                                if (dr2 < rvdw2)
                                {
                                    //printf("%d,%d\n",atomi,atomj);
                                    nblist.va_list[atomi][atomj] = t;
                                    nblist.va_list[atomj][atomi] = 26 - t;
                                    continue;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] nbbinx;
    delete[] nbbiny;
    delete[] nbbinz;
    //delete atoms in excl list
    
    for (int i = 0; i < atom_size; i++)
    {   
        int j = 0;
        while (excl_list.bonded_ranks[i][j] != -1)
        {
            k = excl_list.bonded_ranks[i][j];
            nblist.va_list[i][k] = -1;
            nblist.va_list[k][i] = -1;
            //printf("%d,%d\n",i,k);
            j ++;
        }
    }
    
}
