/*Header file to accompany iz_util.c*/
/*Written by Dr Robert Stewart for Stewart & Bair, 2009*/

#ifndef INC_IZ_UTIL_H
#define INC_IZ_UTIL_H

#define MAX_IN 20
typedef struct {
	double E; 											/*Electrical elastance 1/C*/
	double vr;											/*Resting mebrane potential*/
	double k;												/*Scaling constant*/
	double l;												/*-k*vt*/
	double a;												/*Recovery variable rate constant*/
	double b;												/*Scaling constant*/
	double v_peak;									/*Peak voltage during spike*/
	double v_reset;									/*Post-spike reset potential*/
	double u_step;									/*Post-spike recovery variable step*/
	double I;												/*Input current*/
	double g_ampa;									/*AMPA conductance*/
	double g_gaba;									/*GABA_A conductance*/
	double E_ampa;									/*AMPA reversal potential*/
	double E_gaba;									/*GABA reversal potential*/
	double v;												/*Membrane voltage*/
	double u;												/*Recovery variable*/
	int n_in;												/*Number of synaptic events in input buffers*/
	double in_t[MAX_IN];						/*Time input buffer*/
	float in_w[MAX_IN];							/*Weight input buffer*/
} neuron_iz;

typedef struct {
	unsigned int n; 								/*Target neuron*/
	float w; 												/*Weight*/
	unsigned int i;									/*Delay interval*/
} synapse;

void iz_derivs(double *, double *, double *);
void iz_first(double **, double **, double *);
void iz_iter(double **, double **, double *, int);

#endif /* INC_IZ_UTIL_H */
