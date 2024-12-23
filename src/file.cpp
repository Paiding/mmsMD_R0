#include "file.h"

int ini_get_int(char *key, char *filename)
{
    int data;
    FILE *fp = NULL;
    bool fileNotEnd = true;
    if (fp = fopen(filename, "r"))
    {
        int line = 0;
        int len = strlen(key);
        char strline[MAX_STRLEN];
        while (fileNotEnd)
        {
            
            if (fgets(strline, MAX_STRLEN, fp)!=NULL)
            {//reaching the end/
                if (strlen(strline) > len)
                {//longer than the key first
                    if (strncmp(key, strline, len) == 0)
                    {//leftvalue = key?
                        char *tik = strtok(strline, "=");
                        strTrim(tik);
                        if (strcmp(tik,key) == 0)
                        {//case: "keynote" to "key" cannot be read as data
                            tik = strtok(NULL, "=");
                            data = atoi(tik);
                            fileNotEnd = false;
                        }
                    }
                }
            }
            else
            {
                fileNotEnd = false;
            }
        }
        fclose(fp);
    }
    else
    {
        perror("can't open ini");
        return(-1);
    }
    
    
    return data;
}

float ini_get_float(char *key, char *filename)
{
    float data;
    FILE *fp = NULL;
    bool fileNotEnd = true;
    if (fp = fopen(filename, "r"))
    {
        int line = 0;
        int len = strlen(key);
        char strline[MAX_STRLEN];
        while (fileNotEnd)
        {
            
            if (fgets(strline, MAX_STRLEN, fp)!=NULL)
            {//reaching the end/
                if (strlen(strline) > len)
                {//longer than the key first
                    if (strncmp(key, strline, len) == 0)
                    {//leftvalue = key?
                        char *tik = strtok(strline, "=");
                        strTrim(tik);
                        if (strcmp(tik,key) == 0)
                        {//case: "keynote" to "key" cannot be read as data
                            tik = strtok(NULL, "=");
                            data = atof(tik);
                            fileNotEnd = false;
                        }
                    }
                }
            }
            else
            {
                fileNotEnd = false;
            }
        }
    }
    else
    {
        perror("can't open ini");
        return(-1);
    }
    
    return data;
}

void strTrim(char *pStr)
{
        char *pTmp = pStr;  
      
    while (*pStr != '\0')   
    {  
        if (*pStr != ' ')  
        {  
            *pTmp++ = *pStr;  
        }  
        ++pStr;  
    }  
    *pTmp = '\0';  
}

char* strLTrim(char *pStr)
{
    while (*pStr == ' ')
    {
        ++pStr;
    }

    return pStr;
}

bool ini_check_str(char *key,char *string, char *filename)
{
    bool b;
    FILE *fp = NULL;
    bool fileNotEnd = true;
    if (fp = fopen(filename, "r"))
    {
        int line = 0;
        int len = strlen(key);
        char strline[MAX_STRLEN];
        while (fileNotEnd)
        {
            if (fgets(strline, MAX_STRLEN, fp)!=NULL)
            {//reaching the end/
                if (strlen(strline) > len)
                {//longer than the key first
                    if (strncmp(key, strline, len) == 0)
                    {//leftvalue = key?
                        char *tik = strtok(strline, "=");
                        strTrim(tik);
                        if (strlen(tik) == len)
                        {//case: "keynote" to "key" cannot be read as data
                            tik = strtok(NULL, "=");
                            tik = strLTrim(tik);
                            b = strcmp(string,tik);
                            fileNotEnd = false;
                        }
                    }
                }
            }
            else
            {
                fileNotEnd = false;
            }
        }
    }
    else
    {
        perror("can't open ini");
        return(-1);
    }
    
    return !b;
}

void ini_get_str(char *key, char *filename, char *dest)
{
    FILE *fp = NULL;
    bool fileNotEnd = true;
    if (fp = fopen(filename, "r"))
    {
        int line = 0;
        int len = strlen(key);
        char strline[MAX_STRLEN];
        while (fileNotEnd)
        {
            if (fgets(strline, MAX_STRLEN, fp)!=NULL)
            {//reaching the end/
                if (strlen(strline) > len)
                {//longer than the key first
                    if (strncmp(key, strline, len) == 0)
                    {//leftvalue = key?
                        char *tik = strtok(strline, "=");
                        strTrim(tik);
                        if (strlen(tik) == len)
                        {//case: "keynote" to "key" cannot be read as data
                            tik = strtok(NULL, "=");
                            tik = strLTrim(tik);
                            strcpy(dest,tik);
                            if (dest[strlen(dest)-1] == '\n')
                            {//cut the line break character
                                dest[strlen(dest)-1] = '\0';
                            }
                            
                            fileNotEnd = false;
                        }
                    }
                }
            }
            else
            {
                fileNotEnd = false;
            }
        }
    }
    else
    {
        perror("can't open ini");
    }
    
    //return b;
}



void csv_print_data(char* filename, t_atom_kinematics & atom_kinematics)
{
    FILE *fp = NULL;
    int l = atom_kinematics.ax.getsize();
    if (fp = fopen(filename,"wt"))
    {
        for (int i = 0; i < l; i++)
        {
            fprintf(fp,"%d,1,HE,HE,%f,%f,%f\n",i,atom_kinematics.rx[i],atom_kinematics.ry[i],atom_kinematics.rz[i]);
        }
        fclose(fp);
    }
}

void log_print_data(int step, char* filename, t_atom_kinematics & atom_kinematics)
{
    
    FILE *fp = NULL;
    int l = atom_kinematics.ax.getsize();
    if (fp = fopen(filename,"a"))
    {
        fprintf(fp,"step %d\n",step);
        for (int i = 0; i < l; i++)
        {
            fprintf(fp,"atom %d\n",i);
            fprintf(fp,"rx = %f,ry = %f,rz = %f\n",atom_kinematics.rx[i],atom_kinematics.ry[i],atom_kinematics.rz[i]);
            fprintf(fp,"vx = %f,vy = %f,vz = %f\n",atom_kinematics.vx[i],atom_kinematics.vy[i],atom_kinematics.vz[i]);
            fprintf(fp,"ax = %f,ay = %f,az = %f\n",atom_kinematics.ax[i],atom_kinematics.ay[i],atom_kinematics.az[i]);
        }
        fclose(fp);
    }
}

void log_print_data_2(int step, char* filename, t_atom_kinematics & atom_kinematics, t_atom_statics & atom_statics)
{
    
    FILE *fp = NULL;
    int l = atom_kinematics.ax.getsize();
    if (fp = fopen(filename,"a"))
    {
        fprintf(fp,"step %d\n",step);
        for (int i = 0; i < l; i++)
        {
            fprintf(fp,"atom %d\n",i);
            fprintf(fp,"rx = %f,ry = %f,rz = %f\n",atom_kinematics.rx[i],atom_kinematics.ry[i],atom_kinematics.rz[i]);
            fprintf(fp,"vx = %f,vy = %f,vz = %f\n",atom_kinematics.vx[i],atom_kinematics.vy[i],atom_kinematics.vz[i]);
            fprintf(fp,"ax = %f,ay = %f,az = %f\n",atom_kinematics.ax[i],atom_kinematics.ay[i],atom_kinematics.az[i]);
            fprintf(fp,"m = %f,type = %d\n",atom_statics.m[i],atom_statics.type[i]);
        }
        fclose(fp);
    }
}

char* itoa(int num,char* str,int radix)
{//turn a int number to string,radix for the scale system of number/2 8 10 16 ...
    char index[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    unsigned unum;//the absolute value of the int
    int i=0,j,k;//i for the loc of char；k for the start of string;j for the exchange of num
 
    //get the abs
    if(radix==10&&num<0)//minus number
    {
        unum=(unsigned)-num;
        str[i++]='-';
    }
    else unum=(unsigned)num;//> 0
 
    //building string
    do
    {
        str[i++]=index[unum%(unsigned)radix];
        unum/=radix;
 
    }while(unum);
 
    str[i]='\0';
 
    //reverse
    if(str[0]=='-') k=1;//if < 0 don't change the minus character
    else k=0;//>0 reverse all
 
    char temp;
    for(j=k;j<=(i-1)/2;j++)
    {
        temp=str[j];
        str[j]=str[i-1+k-j];
        str[i-1+k-j]=temp;
    }
 
    return str;//turned string
 
}

int collect_int(const char* str)
{//get number from special string,Ex: get 23 from "atom (23):"
    char cpy[128];
    strcpy(cpy,str);
    char* tik = strtok(cpy,"(");
    tik = strtok(NULL,")");
    return atoi(tik);
}

void find_float_list(char* filename, char* key,
                     int size,float *dest)
{
    char dest1[32];
    char* pstr;
    float tmp;
    FILE *fp = NULL;
    bool fileNotEnd = true;
    if (fp = fopen(filename, "r"))
    {
        int line = 0;
        int len = strlen(key);
        char strline[MAX_STRLEN];
        while (fileNotEnd)
        {
            
            if (fgets(strline, MAX_STRLEN, fp)!=NULL)
            {//reaching the end/
                if (strlen(strline) > len)
                {//longer than the key first
                    sscanf(strline,"%s",dest1);
                    if (strncmp(key, dest1, len) == 0)
                    {//leftvalue = key?
                        pstr = strtok(strline,":");
                        pstr = strtok(NULL,":");
                        for (int i = 0; i < size; i++)
                        {
                            while (*pstr == ' ' || *pstr == '\t') pstr++;
                            sscanf(pstr,"%f",&tmp);
                            dest[i] = tmp;
                            pstr++;
                            while ((*pstr >= '0' && *pstr <= '9')||*pstr == '.') pstr++;
                        }
                        
                    }
                }
            }
            else
            {
                fileNotEnd = false;
            }
        }
        fclose(fp);
    }
    else
    {
        perror("can't open ini");
        return;
    }
    
}
