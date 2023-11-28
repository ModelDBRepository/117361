/*Header file to accompany tm_util.c*/
/*Written by Dr Robert Stewart for Stewart & Bair, 2009*/

#ifndef INC_TM_UTIL_H
#define INC_TM_UTIL_H
typedef unsigned int uint32; /*synonym for unsigned 32 bit integer*/

typedef struct {
	uint32 n; /*target neuron*/
	float w; 	/*weight*/
	uint32 i;	/*delay interval upgraded to 32 bit integer*/
} synapse;

typedef struct {
	double v;
	double n;
	double m;
	double h;
	double a;
	double b;
	double c;
	double d;
	double I;
	double g_ampa;
	double g_gaba;
	uint32 n_out;
} neuron_tm; 

void tm_derivs(double *, double *, double *);
void tm_first(double **, double **, double *);
void tm_iter(double **, double **, double *, int);

#endif /* INC_TM_UTIL_H */
