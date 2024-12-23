#include "settle.h"

//check the original file when turn to CUDA!

float mO = 15.9994001;
float mH = 1.00800002;
double wo = 0.88809573957514265;
double wh = 0.055952130212428682;
double wohh = 18.015400171279907;

#define __WATER_SPC__
#ifdef __WATER_SPC__
float rc = 0.0816499963;
float ra = 0.00646074256;
float rb = 0.0512738116;
float rc2 = 0.163299993;
#endif

//#define __WATER_OPLSAA__
#ifdef __WATER_OPLSAA__
float rc = 0.0756950006;
float ra = 0.00655627763;
float rb = 0.0520319976;
float rc2 = 0.151390001;
#endif
float rone = 1.0;


void settle(t_molecule & water,
            t_atom_kinematics & atom_kinematics,
            t_atom_kinematics & old_k,
            t_atom_statics & atom_statics,
			t_virial & virial,
			float invdt)
{//settle algorithm for one water molecule

    //get the location of 3 atoms in atom_kinematics
    int ow1_index = water.atom_ranks[0];
    int hw2_index = water.atom_ranks[1];
    int hw3_index = water.atom_ranks[2];

    float x_ow1 = old_k.rx[ow1_index]; //old xyz of 3 atoms
    float y_ow1 = old_k.ry[ow1_index]; 
    float z_ow1 = old_k.rz[ow1_index];
    float x_hw2 = old_k.rx[hw2_index];
    float y_hw2 = old_k.ry[hw2_index];  
    float z_hw2 = old_k.rz[hw2_index];
    float x_hw3 = old_k.rx[hw3_index];
    float y_hw3 = old_k.ry[hw3_index];  
    float z_hw3 = old_k.rz[hw3_index];

	float b4[9];
	b4[0] = old_k.rx[ow1_index]; //old xyz of 3 atoms
   	b4[1] = old_k.ry[ow1_index]; 
    b4[2] = old_k.rz[ow1_index];
    b4[3] = old_k.rx[hw2_index];
    b4[4] = old_k.ry[hw2_index];  
    b4[5] = old_k.rz[hw2_index];
    b4[6] = old_k.rx[hw3_index];
    b4[7] = old_k.ry[hw3_index];  
    b4[8] = old_k.rz[hw3_index];

    float xb0 = x_hw2 - x_ow1; //old shape of water triangle
    float yb0 = y_hw2 - y_ow1;
    float zb0 = z_hw2 - z_ow1;
    float xc0 = x_hw3 - x_ow1;
    float yc0 = y_hw3 - y_ow1;
    float zc0 = z_hw3 - z_ow1;

    x_ow1 = atom_kinematics.rx[ow1_index];//new xyz of 3 atoms
    y_ow1 = atom_kinematics.ry[ow1_index]; 
    z_ow1 = atom_kinematics.rz[ow1_index];
    x_hw2 = atom_kinematics.rx[hw2_index];
    y_hw2 = atom_kinematics.ry[hw2_index];
    z_hw2 = atom_kinematics.rz[hw2_index];
    x_hw3 = atom_kinematics.rx[hw3_index];
    y_hw3 = atom_kinematics.ry[hw3_index];  
    z_hw3 = atom_kinematics.rz[hw3_index];

    float gama, beta, alpa, xcom, ycom, zcom, al2be2, tmp, tmp2;
	float axlng, aylng, azlng, trns11, trns21, trns31, trns12, trns22, 
				trns32, trns13, trns23, trns33, cosphi, costhe, sinphi, sinthe, 
				cospsi, xaksxd, yaksxd, xakszd, yakszd, zakszd, zaksxd, xaksyd, 
				 xa1;
	float ya1, za1, xb1, yb1;
	float zb1, xc1, yc1, zc1, yaksyd, zaksyd, sinpsi, xa3, ya3, za3, 
				xb3, yb3, zb3, xc3, yc3, zc3, xb0d, yb0d, xc0d, yc0d, 
				za1d, xb1d, yb1d, zb1d, xc1d, yc1d, zc1d, ya2d, xb2d, yb2d, yc2d, 
				xa3d, ya3d, za3d, xb3d, yb3d, zb3d, xc3d, yc3d, zc3d;
	float t1,t2;
	float dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;

    xcom = (x_ow1 * wo + (x_hw2 + x_hw3) * wh);//center of mass of water now
	ycom = (y_ow1 * wo + (y_hw2 + y_hw3) * wh);
	zcom = (z_ow1 * wo + (z_hw2 + z_hw3) * wh);

    xa1 = x_ow1 - xcom;//distance between atom and center of mass
	ya1 = y_ow1 - ycom;
	za1 = z_ow1 - zcom;
	xb1 = x_hw2 - xcom;
	yb1 = y_hw2 - ycom;
	zb1 = z_hw2 - zcom;
	xc1 = x_hw3 - xcom;
	yc1 = y_hw3 - ycom;
	zc1 = z_hw3 - zcom;
    //printf("distance[%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",xa1,ya1,za1,xb1,yb1,zb1,xc1,yc1,zc1);
    xakszd = yb0 * zc0 - zb0 * yc0;//vector product to create new coordinate system
	yakszd = zb0 * xc0 - xb0 * zc0;
	zakszd = xb0 * yc0 - yb0 * xc0;
	xaksxd = ya1 * zakszd - za1 * yakszd;
	yaksxd = za1 * xakszd - xa1 * zakszd;
	zaksxd = xa1 * yakszd - ya1 * xakszd;
	xaksyd = yakszd * zaksxd - zakszd * yaksxd;
	yaksyd = zakszd * xaksxd - xakszd * zaksxd;
	zaksyd = xakszd * yaksxd - yakszd * xaksxd;
	/* 27 flops */
    //compute unit lenth of new coordinate system
	axlng = rsqrtf(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
	aylng = rsqrtf(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
	azlng = rsqrtf(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);
	//printf("{%f,%f,%f}\n",axlng,aylng,azlng);

	trns11 = xaksxd * axlng;
	trns21 = yaksxd * axlng;
	trns31 = zaksxd * axlng;
	trns12 = xaksyd * aylng;
	trns22 = yaksyd * aylng;
	trns32 = zaksyd * aylng;
	trns13 = xakszd * azlng;
	trns23 = yakszd * azlng;
	trns33 = zakszd * azlng;
	/* 24 flops */
//turn vectors to new coordinate
	xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
	yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
	xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
	yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;

	za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;

	xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
	yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
	zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
	xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
	yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
	zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;
	/* 65 flops */
	//printf("[%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",za1d,zb1d,zc1d,xc1,yc1,zc1,trns13,trns23,trns33);
//calculate rotation angle phi psi theta
	sinphi = za1d / ra;
	//printf("sinf = %f\n",sinphi);
	tmp    = rone - sinphi * sinphi;
	//printf("tmp = %f\n",tmp);
	bool dosettle = true;
	if (tmp <= 0.0f) 
	{
		cosphi = 0.0f;
		dosettle = false;
		//printf("[%d,%d,%d]",ow1_index,hw2_index,hw3_index);
	}
	else
		cosphi = tmp*rsqrtf(tmp);

	sinpsi = (zb1d - zc1d) / (rc2 * cosphi);
	tmp2   = rone - sinpsi * sinpsi;
	if (tmp2 <= 0.0f) {
		cospsi = 0.0f;
		dosettle = false;
	}
	else
		cospsi = tmp2*rsqrtf(tmp2);
	/* 46 flops */

	ya2d =  ra * cosphi;
	xb2d = -rc * cospsi;
	t1   = -rb * cosphi;
	t2   =  rc * sinpsi * sinphi;
	yb2d =  t1 - t2;
	yc2d =  t1 + t2;
	/* 7 flops */

	/*     --- Step3  al,be,ga 		      --- */
	alpa   = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
	beta   = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
	gama   = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
	al2be2 = alpa * alpa + beta * beta;
	tmp2   = (al2be2 - gama * gama);
	float tmp3;
	if (tmp2 <= 0.0f)
	{
		tmp3 = 0.0f;
	}
	else
	{
		tmp3 = tmp2*rsqrtf(tmp2);
	}
	
	sinthe = (alpa * gama - beta * tmp3) / al2be2;
	/* 47 flops */

	/*  --- Step4  A3' --- */
	tmp2  = rone - sinthe *sinthe;
	if (tmp2 <= 0.0f)
	{
		costhe = 0.0f;
	}
	else
	{
		costhe = tmp2*rsqrtf(tmp2);
	}
	xa3d = -ya2d * sinthe;
	ya3d = ya2d * costhe;
	za3d = za1d;
	xb3d = xb2d * costhe - yb2d * sinthe;
	yb3d = xb2d * sinthe + yb2d * costhe;
	zb3d = zb1d;
	xc3d = -xb2d * costhe - yc2d * sinthe;
	yc3d = -xb2d * sinthe + yc2d * costhe;
	zc3d = zc1d;
	/* 26 flops */

	/*    --- Step5  A3 --- */

	xa3 = trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
	ya3 = trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
	za3 = trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
	xb3 = trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
	yb3 = trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
	zb3 = trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
	xc3 = trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
	yc3 = trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
	zc3 = trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;
	/* 45 flops */
	//printf("a3[%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",xa3,ya3,za3,xb3,yb3,zb3,xc3,yc3,zc3);
	x_ow1 = xcom + xa3;
	y_ow1 = ycom + ya3;
	z_ow1 = zcom + za3;
	x_hw2 = xcom + xb3;
	y_hw2 = ycom + yb3;
	z_hw2 = zcom + zb3;
	x_hw3 = xcom + xc3;
	y_hw3 = ycom + yc3;
	z_hw3 = zcom + zc3;
	
	if (dosettle)
	{
		atom_kinematics.rx[ow1_index] = x_ow1;
		atom_kinematics.ry[ow1_index] = y_ow1;
    	atom_kinematics.rz[ow1_index] = z_ow1;
    	atom_kinematics.rx[hw2_index] = x_hw2;
    	atom_kinematics.ry[hw2_index] = y_hw2;
    	atom_kinematics.rz[hw2_index] = z_hw2;
    	atom_kinematics.rx[hw3_index] = x_hw3;
    	atom_kinematics.ry[hw3_index] = y_hw3;
    	atom_kinematics.rz[hw3_index] = z_hw3;

		
	}
	else
	{
		//do shake
		
		int iconv;
  		int iatom[3]={0,0,1};
  		int jatom[3]={1,2,2};
  		float rijx,rijy,rijz,tx,ty,tz,im,jm,acor,rp,diff;
 		int i,ll,ii,jj,l3,ix,iy,iz,jx,jy,jz,conv;
		float bond[9];


		float invmass[3];
		invmass[0] = 1.0/mO;
		invmass[1] = 1.0/mH;
		invmass[2] = 1.0/mH;

		float bondsq[3];
		bondsq[0] = 0.01;//dOH2
		bondsq[1] = bondsq[0];
		bondsq[2] = 0.16330 * 0.16330;

		float M2[3];
		M2[0] = 1.0/(2.0*(invmass[0] + invmass[1]));
		M2[1] = M2[0];
		M2[0] = 1.0/(2.0*(invmass[1] + invmass[2]));

		for(ll=0;ll<3;ll++) 
		{
    		l3=3*ll;
    		ix=3*iatom[ll];
    		jx=3*jatom[ll];
    		for(i=0;i<3;i++) 
      			bond[l3+i]= b4[ix+i] - b4[jx+i];
  		}

		float after[9];
		after[0] = atom_kinematics.rx[ow1_index]; //new xyz of 3 atoms
   		after[1] = atom_kinematics.ry[ow1_index]; 
    	after[2] = atom_kinematics.rz[ow1_index];
    	after[3] = atom_kinematics.rx[hw2_index];
    	after[4] = atom_kinematics.ry[hw2_index];  
    	after[5] = atom_kinematics.rz[hw2_index];
    	after[6] = atom_kinematics.rx[hw3_index];
    	after[7] = atom_kinematics.ry[hw3_index];  
    	after[8] = atom_kinematics.rz[hw3_index];

		for(i=0,iconv=0;i<1000 && iconv<3; i++)
  		{
			for(ll=0;ll<3;ll++) 
			{
      			ii = iatom[ll];
      			jj = jatom[ll];
      			l3 = 3*ll;
      			ix = 3*ii;
      			jx = 3*jj;
      			iy = ix+1;
      			jy = jx+1;
      			iz = ix+2;
      			jz = jx+2;

      			rijx = bond[l3];
      			rijy = bond[l3+1];
      			rijz = bond[l3+2];  
      
      			tx   = after[ix]-after[jx];
      			ty   = after[iy]-after[jy];
      			tz   = after[iz]-after[jz];
      
      			rp   = tx*tx+ty*ty+tz*tz;
      			diff = bondsq[ll] - rp;

      			if(fabs(diff)<1e-8) 
				{
					iconv++;
      			} 
				else 
				{
					rp = rijx*tx+rijy*ty+rijz*tz;
				}
				acor = diff*M2[ll]/rp;
				im           = invmass[ii];
				jm           = invmass[jj];
				tx           = rijx*acor;
				ty           = rijy*acor;
				tz           = rijz*acor;
				after[ix] += tx*im;
				after[iy] += ty*im;
				after[iz] += tz*im;
				after[jx] -= tx*jm;
				after[jy] -= ty*jm;
				after[jz] -= tz*jm;
      		}
    	}
		atom_kinematics.rx[ow1_index] = after[0]; //new xyz of 3 atoms
   		atom_kinematics.ry[ow1_index] = after[1];
		atom_kinematics.rz[ow1_index] = after[2];
		atom_kinematics.rx[hw2_index] = after[3];
		atom_kinematics.ry[hw2_index] = after[4];
		atom_kinematics.rz[hw2_index] = after[5];
    	atom_kinematics.rx[hw3_index] = after[6];
		atom_kinematics.ry[hw3_index] = after[7];
		atom_kinematics.rz[hw3_index] = after[8];
  	}
		dax = xa3 - xa1;
		day = ya3 - ya1;
		daz = za3 - za1;
		dbx = xb3 - xb1;
		dby = yb3 - yb1;
		dbz = zb3 - zb1;
		dcx = xc3 - xc1;
		dcy = yc3 - yc1;
		dcz = zc3 - zc1;
	/* 9 flops, counted with the virial */

	if (true)
	{
		atom_kinematics.vx[ow1_index] += dax * invdt;
    	atom_kinematics.vy[ow1_index] += day * invdt;
    	atom_kinematics.vz[ow1_index] += daz * invdt;
    	atom_kinematics.vx[hw2_index] += dbx * invdt;
    	atom_kinematics.vy[hw2_index] += dby * invdt;
    	atom_kinematics.vz[hw2_index] += dbz * invdt;
    	atom_kinematics.vx[hw3_index] += dcx * invdt;
    	atom_kinematics.vy[hw3_index] += dcy * invdt;
    	atom_kinematics.vz[hw3_index] += dcz * invdt;
		
	}
	//cal vir
	if (dosettle)
	{
		float mdax = mO*dax;
        float mday = mO*day;
        float mdaz = mO*daz;
        float mdbx = mH*dbx;
        float mdby = mH*dby;
        float mdbz = mH*dbz;
        float mdcx = mH*dcx;
        float mdcy = mH*dcy;
        float mdcz = mH*dcz;

		float vir_factor = 0.5f * invdt * invdt;
		
		virial.vir_x[ow1_index] -= vir_factor * (b4[0] * mdax);
		virial.vir_y[ow1_index] -= vir_factor * (b4[1] * mday);
		virial.vir_z[ow1_index] -= vir_factor * (b4[2] * mdaz);
		virial.vir_x[hw2_index] -= vir_factor * (b4[3] * mdbx);
		virial.vir_y[hw2_index] -= vir_factor * (b4[4] * mdby);
		virial.vir_z[hw2_index] -= vir_factor * (b4[5] * mdbz);
		virial.vir_x[hw3_index] -= vir_factor * (b4[6] * mdcx);
		virial.vir_y[hw3_index] -= vir_factor * (b4[7] * mdcy);
		virial.vir_z[hw3_index] -= vir_factor * (b4[8] * mdcz);
		
	}
}