// 2021/02/07 
#include "md.h"
int main(int argc, char **argv)
{
	switch (argc)
	{
	case 2:
		//printf("hello, pretreat time\n");
		gro2csv("in/ala.gro","in/ala-ala.csv");
		//csv2gro("in/unstable.csv","out/s16.gro");
		//csv2gro("out/finalstate.csv","out/s2.gro");
		break;
	
	default://default mode, mdrun
	clock_t start_t, end_t;//time conting
	printf("Start\n");
	start_t = clock();
	//initialization
	int step = 0;
	t_force force;//declaring data structures
	t_virial virial;
	t_nstep nstep;
	t_atom_statics atom_statics;
	t_atom_kinematics atom_kinematics;
	t_atom_kinematics old_k;

	array<float> lj_para;
	t_ub_para ub_para;
	t_g96a_para g96a_para;
	t_pdih_para pdih_para;
	t_ipdih_para ipdih_para;
	t_lj14_para lj14_para;
	t_constr_para constr_para;
	t_lincs_data lincs_data;
	
	read_functype(lj_para,ub_para,g96a_para,pdih_para,ipdih_para,lj14_para,constr_para,atom_statics);
	printf("finished to read functype\n");
	t_tcstat tcstat;
	t_tcpara tcpara;
	t_pcstat pcstat;
	t_pcpara pcpara;
	t_nblist nblist;
	
	t_excl_list excl_list;
	t_radius radius;
	t_dsf_para dsf_para;
	t_shift_switch_para sspara;
	t_vcm vcm;
	
	float box[3];
	int tctype = 0;//0 for notc ,1 for Berendsen
	
	tctype = 1 * ini_check_str("tctype", "Berendsen\n", "in/settings.ini")
           + 2 * ini_check_str("tctype", "NoseHoover\n", "in/settings.ini");

	
	//printf("start molemap\n");
	
	
    t_mole_map mole_map;
	printf("end molemap\n");
	
	int mole_size = mole_map.fullmolenum;
	array<t_molecule> mole_array(mole_size);
	for (int i = 0; i < mole_size; i++)
	{
		t_molecule* mp = new t_molecule;
		mole_array.put(*mp);
	}
	printf("molearray prepared\n");
	float dt = ini_get_float("dt","in/settings.ini");
	t_group tcgrp;
	
	init_box(box);
//initialization
	
	int atom_size = init_all(atom_kinematics,
							 old_k,
            				 atom_statics,
							 tcstat,tcpara,
							 pcstat,pcpara,
              				 force,virial,
							 nstep,
							 radius,
							 mole_array,
							 mole_map,
							 dsf_para,
							 vcm
							 );

	if (tcpara.N > 1)
	{
		init_tc_group(tcgrp,tcpara.N);
	}
	else
		init_tc_nogroup(tcgrp,atom_size);
	
	
	init_shift_switch_para(radius,sspara);
	init_lincs_data(mole_array,lincs_data,atom_statics,constr_para);

	

	printf("%f\n",tcpara.nrdf[0]);

	printf("initialized\n");
	printf("tctype = %d\n",tctype);


	//step = 1;//zero step is initialization,first step is now
	//some test code to be merged
	for (int i = 0; i < mole_array.getsize(); i++)
	{
		fix_broken_water(0.2,mole_array[i],atom_kinematics,box);
		
		//check_broken_molecule(0.2,mole_array[i],atom_kinematics);
	}
	kinematics_copy(atom_kinematics,old_k);
	
	printf("atom%d,mole%d\n",atom_size,mole_array.getsize());
	printf("exlist start\n");
	build_excl_list(mole_array, mole_map, excl_list, atom_size);
	printf("nblist start\n");
	/*
	for (int i = 0; i < atom_size; i++)
    {
        for (int j = 0; j < atom_size; j++)
        {
            printf("%d,",excl_list.bonded_ranks[i][j]);
        }
        printf("\n");
    }*/
	
	
	build_nblist_bin(box, nblist, atom_size,radius, atom_kinematics,excl_list);
	/*
	for (int i = 0; i < atom_size; i++)
    {
        for (int j = 0; j < atom_size; j++)
        {
            printf("%d,",nblist.va_list[i][j]);
        }
        printf("\n");
    }*/
	
	printf("nblist finish\n");
	
	
	int fullstep = nstep.full_step;
	if (false)
	{
		gen_vel(atom_kinematics,
				atom_statics,
				mole_array,
				tcpara);
		
	}
	calc_mass_grp(atom_statics,vcm);
	

//main cycle here
	
	while (step++ < fullstep)
	{
		//printf("step %d\n", step);
		one_step(step,nstep,dt,
				box,lj_para,
				dsf_para,
				sspara,ub_para,g96a_para,pdih_para,ipdih_para,lj14_para,
				atom_size,radius,
				nblist,
				atom_kinematics
				,old_k,atom_statics,
				tcstat,tcpara,
				pcstat,pcpara,
				lincs_data,
				force,
				virial,
				tctype,
				excl_list,
				mole_map,
				mole_array,
				vcm,tcgrp);
			
	}

//saving state and conting time
	
	
	end_t = clock();
	float total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
	printf("End, time cost is %fs\n", total_t);


	//*/
	break;
	
	}
	return 0;
}

