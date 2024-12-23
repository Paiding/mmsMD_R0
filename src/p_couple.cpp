#include "p_couple.h"

void calc_vir(t_virial & virial,
              t_force & force,
              t_atom_kinematics & atom_kinematics
             )
{
    int size = atom_kinematics.rx.getsize();
    for (int i = 0; i < size; i++)
    {
        virial.vir_x[i] -= 0.5 * atom_kinematics.rx[i] * force.fx[i];
        virial.vir_y[i] -= 0.5 * atom_kinematics.ry[i] * force.fy[i];
        virial.vir_z[i] -= 0.5 * atom_kinematics.rz[i] * force.fz[i];
    }
}

void calc_pres(t_tcstat & tcstat,
               t_pcstat & pcstat,
               t_virial & virial,
               float *box
              )
{
    //int vir_size = virial.vir_x.getsize();
    virial.sum[XX] = virial.vir_x.reduce();
    virial.sum[YY] = virial.vir_y.reduce();
    virial.sum[ZZ] = virial.vir_z.reduce();
    /*
    for (int i = 0; i < vir_size; i++)
    {
        virial.sum[XX] += virial.vir_x[i];
        virial.sum[YY] += virial.vir_y[i];
        virial.sum[ZZ] += virial.vir_y[i];
    }*/
    //virial.sum[XX] /= MILLION;//g/mol m2/s2 to KJ/mol
    //virial.sum[YY] /= MILLION;
    //virial.sum[ZZ] /= MILLION;
    
    //float fac = float(PRESFAC)*2.0f /((box[XX]+0.16) * (box[YY]+0.16) * (box[ZZ]+0.16));
    float fac = float(PRESFAC)*2.0f /((box[XX]) * (box[YY]) * (box[ZZ]));
    pcstat.pres[XX] = fac * (tcstat.ekin_tot[XX] - virial.sum[XX]);
    pcstat.pres[YY] = fac * (tcstat.ekin_tot[YY] - virial.sum[YY]);
    pcstat.pres[ZZ] = fac * (tcstat.ekin_tot[ZZ] - virial.sum[ZZ]);

    pcstat.pres_scalar = (pcstat.pres[XX] + pcstat.pres[YY] + pcstat.pres[ZZ]) / 3.0f;

    pcstat.pres_prev[XX] = pcstat.pres[XX];
    pcstat.pres_prev[YY] = pcstat.pres[YY];
    pcstat.pres_prev[ZZ] = pcstat.pres[ZZ];
}

void Berendsen_pc_mu(t_pcstat & pcstat,
                     t_pcpara & pcpara
                    )
{//computer the mu
    float tmp_pres_scalar = pcstat.pres_scalar;
    pcstat.mu[XXXX] = 1.0f - pcpara.compress[XXXX]*pcpara.dt*(pcpara.p_ref[XXXX]-tmp_pres_scalar) / (pcpara.tau_p * DIM);
    pcstat.mu[YYYY] = 1.0f - pcpara.compress[YYYY]*pcpara.dt*(pcpara.p_ref[YYYY]-tmp_pres_scalar) / (pcpara.tau_p * DIM);
    pcstat.mu[ZZZZ] = 1.0f - pcpara.compress[ZZZZ]*pcpara.dt*(pcpara.p_ref[ZZZZ]-tmp_pres_scalar) / (pcpara.tau_p * DIM);
}

void Berendsen_pc_rx_scale(t_pcstat & pcstat,
                           t_atom_kinematics & atom_kinematics
                          )
{
    int atom_size = atom_kinematics.rx.getsize();
    float mu_xx = pcstat.mu[XXXX];
    float mu_yy = pcstat.mu[YYYY];
    float mu_zz = pcstat.mu[ZZZZ];
    for (int i = 0; i < atom_size; i++)
    {
        atom_kinematics.rx[i] *= mu_xx; 
        atom_kinematics.ry[i] *= mu_yy; 
        atom_kinematics.rz[i] *= mu_zz; 
    }
}

void Berendsen_pc_box_scale(t_pcstat & pcstat,
                           float * box
                           )
{
    float mu_xx = pcstat.mu[XXXX];
    float mu_yy = pcstat.mu[YYYY];
    float mu_zz = pcstat.mu[ZZZZ];
    box[XX] *= mu_xx;
    box[YY] *= mu_yy;
    box[ZZ] *= mu_zz;
}