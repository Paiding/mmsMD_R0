#ifndef _VIRIAL_H_
#define _VIRIAL_H_

#include "includes.h"

typedef struct{
	int Nr;
	float a[3][3]; // 2D array, width: Nr+plus, height:DIM*DIM (9)

	// tmp sum
	int block_size;
	int N_tmp;
	float a_tmp[3][3]; // 2D array, width: N_tmp+plus, height:DIM*DIM (9)
} t_vira;

typedef struct{
	array<float> vir_x;
    array<float> vir_y;
    array<float> vir_z;
	float sum[3] = {0,0,0};
} t_virial;


#endif