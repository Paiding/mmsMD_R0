#ifndef _physics_h
#define _physics_h
#define KILO 1000
#define MILLION 1000000
#define BOLTZMANN	 (1.380658e-23)		/* (J/K)	*/
#define AVOGADRO	 (6.0221367e23)		/* ()		*/
#define RGAS             (8.314511212)   /* (J/(mol K))  */
#define BOLTZ            (RGAS/KILO)            /* (kJ/(mol K)) */
#define BOLTZ_S        ((RGAS * KILO))     /*g/mol * nm^2 / ns^2 */
//#define K_COULUMB  (1.38901566e-4) /*KJ/mol*/
#define K_COULUMB (138.9) //KJ/mol nm e-2
#define PRESFAC (16.6054)
#define RAD2DEG      (180.0/M_PI)                          /* Conversion	*/
#define DEG2RAD      (M_PI/180.0)                          /* id		*/
#endif