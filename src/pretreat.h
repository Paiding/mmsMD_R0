#ifndef _PRETREAT_H
#define _PRETREAT_H
#include "includes.h"
#include "file.h"
#define GRO_STRLEN 128

void gro2csv(char* filename,char* savefile);
void csv2gro(char* read_file,
               char* save_file
              );

class t_atom_map
{
  public:
  //small data structure for pretreat
    int map_size = 0;
    //build in build_atom_map_1
    int moletypes = 1;
    array<char*> atom_type; //string eg OW,HW
    array<char*> atom_name;//eg HW1,HW2
    array<char*> mole_name;//eg SOL
    array<int> atom_loc; //serial in molecular;eg 12345123 12345 for 1st type then 123 for 2nd type of molecular
    //build in build_atom_map_2
    array<float> atom_charge; //electrostatic charge
    array<float> atom_mass;
    array<float> epsilon;//LJ epsilon
    array<float> sigma;//LJ sigma

    void build_atom_map_1();
    void build_atom_map_2();

    t_atom_map(){build_atom_map_1();build_atom_map_2();}

    int find_loc_byname(char *atom_name_key, char *mole_name_key);

};

class t_mole_map
{
    public:
    int moletypes = 1;
    int fullmolenum = 0;
    array<int> mole_num;//number of each type of molecule
    array<int> bond_num;//index is one type of molecular
    array<int> angle_num;//so as bondnum
    array<int> pdih_num;
    array<int> ipdih_num;
    array<int> lj14_num;
    array<int> cons_num;
    array<char*> mole_type;

    array<int> bondi;//index is serial of bond
    array<int> bondj;//one bond one element
    array<float> bondrK;
    array<float> bondr0;
  
    array<int> angleType;
    array<int> funcType_a;
    array<int> anglei;
    array<int> anglej;
    array<int> anglek;
    array<float> anglea0;
    array<float> angleaK;

    array<int> funcType_pdih;
    array<int> pdihi;
    array<int> pdihj;
    array<int> pdihk;
    array<int> pdihl;

    array<int> funcType_ipdih;
    array<int> ipdihi;
    array<int> ipdihj;
    array<int> ipdihk;
    array<int> ipdihl;

    array<int> funcType_lj14;
    array<int> lj14i;
    array<int> lj14j;

    array<int> funcType_constr;
    array<int> constri;
    array<int> constrj;

    //new mole map
    array<int> atom_loc;
    array<int> atom_num;
    array<int> atom_type;
    array<float> atom_mass;
    array<float> atom_q;
    array<int> atom_resind;
    array<int> excl_idloc;
    array<int> excl_lsloc;
    array<int> exclist_index;
    array<int> exclist;

    void build_mole_map();
    t_mole_map(){build_mole_map();}

    int find_loc_byname(char* mole_name_key);
};



//init atom data from csv file, return the lines number of csv file
int csv_get_data(char *filename, t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics,
                t_atom_map & atom_map, t_mole_map & mole_map, array<t_molecule> & mole_array
                );

void read_functype(array<float> & lj_para, t_ub_para & ub_para, t_g96a_para & g96a_para,
                   t_pdih_para & pdih_para, t_ipdih_para & ipdih_para,
                   t_lj14_para & lj14_para, t_constr_para & constr_para,
                   t_atom_statics & atom_statics
                  );




#endif