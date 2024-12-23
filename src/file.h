#ifndef _FILE_H
#define _FILE_H
#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "top.h"
#define MAX_STRLEN 128



//quest a value by key from ini file
int ini_get_int(char *key, char *filename);
float ini_get_float(char *key, char *filename);
bool ini_check_str(char *key, char *string, char *filename);
void ini_get_str(char *key, char *filename, char *dest);
//trim the space from the end of a string
void strTrim(char *pStr);
//trim the space from the beginning of a string
char* strLTrim(char *pStr);

void csv_print_data(char* filename, t_atom_kinematics & atom_kinematics);
void log_print_data(int step, char* filename, t_atom_kinematics & atom_kinematics);
void log_print_data_2(int step, char* filename, t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics);
//change an int num to a string
char* itoa(int num,char* str,int radix);

//extract 23 from "atom(23)"
int collect_int(const char* str);

void find_float_list(char* filename, char* key,
                     int size,float *dest);


#endif