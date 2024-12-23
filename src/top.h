#ifndef _TOP_H
#define _TOP_H

#include "includes.h"

struct t_atom_kinematics
{
    array<float> rx;
    array<float> ry;
    array<float> rz;
    array<float> vx;
    array<float> vy;
    array<float> vz;
    array<float> ax;
    array<float> ay;
    array<float> az;
};

struct t_atom_statics
{
    array<float> m;//mass
    array<float> m_inv;
    array<float> q;
    array<int> type;//atom type;eg HW1 OW
    array<int> df;//degree of freedom
    array<bool> iswater;
    array<int> mole_serial;//belong to which molecular
    int atnr = 0;

};

struct t_group
{
  int group_num;
  int *index;
  int *list;
};

struct t_force
{
    array<float> fx;
    array<float> fy;
    array<float> fz;

    float esum_vdw = 0.0f;
    float esum_qq = 0.0f;
    float esum_bond = 0.0f;
    float esum_angle = 0.0f;
    float esum_pdih = 0.0f;
    float esum_idih = 0.0f;
};

struct t_nstep
{
  int full_step = 1;
  int n_nblist = 10;
  int n_log = 1000;
};

struct t_radius
{
  float r_vdw;
  float r_dsf;
};

class t_molecule
{
  //private:
  //static int mole_size;
  public:
  int mole_map_rank;//rank of this type in mole_map 
  bool iswater;
  char mole_type[64];
  array<int> atom_ranks;

  //bond angle dih
  int bond_num = 0;
  array<int> bondi;
  array<int> bondj;
  array<float> bondrK;
  array<float> bondr0;
  
  int angle_num = 0;
  array<int> angleType;//0 for normal;1 for ub;2 for g96
  array<int> funcType_a;
  array<int> anglei;
  array<int> anglej;
  array<int> anglek;
  array<float> anglea0;
  array<float> angleaK;

  int pdih_num = 0;
  array<int> pdihType;
  array<int> pdihi;
  array<int> pdihj;
  array<int> pdihk;
  array<int> pdihl;

  int ipdih_num = 0;
  array<int> ipdihType;
  array<int> ipdihi;
  array<int> ipdihj;
  array<int> ipdihk;
  array<int> ipdihl;

  int lj14_num = 0;
  array<int> lj14Type;
  array<int> lj14i;
  array<int> lj14j;

  int cons_num = 0;
  array<int> consType;
  array<int> consi;
  array<int> consj;

  //array<int> excl_index;
  //array<int> excl_list;

  //functions


};

struct t_ub_para
{
  array<int> type_index;//start index of the functype list
  array<float> theta0;
  array<float> ktheta;
  array<float> r130;
  array<float> kUB;

};

struct t_g96a_para
{
  array<int> type_index;//start index of the functype list
  array<float> theta0;
  array<float> ct;
};

struct t_pdih_para
{
  array<int> type_index;//start index of the functype list
  array<float> phi0;
  array<float> cp0;
  array<float> mult;
};

struct t_ipdih_para
{
  array<int> type_index;//start index of the functype list
  array<float> xi0;
  array<float> cx0;
};

struct t_lj14_para
{
  array<int> type_index;//start index of the functype list
  array<float> c6;
  array<float> c12;
};

struct t_constr_para
{
  array<int> type_index;
  array<float> d0;

};


#endif