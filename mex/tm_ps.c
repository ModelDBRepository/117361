/* =============================================================
 * tm_ps.c Written by Dr Robert Stewart for Stewart & Bair, 2009
 *
 * MEX-file for solving the TM model neuron described in Brette et al 2007
 * using the Parker-Sochacki method, with the 4th order Runge-Kutta and 
 * Bulirsch-Stoer methods for comparison
 * Compile:   mex tm_ps.c tm_util.c integ_util.c
 * Call as:   [ps_v,rk_v,bs_v,t_cpu] = tm_ps(fp,ip);
 * ============================================================= */
#include "mex.h"
#include "matrix.h"
#include <math.h> 
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "tm_util.h"

#define FP      				prhs[0] /* Input Arguments */
#define IP							prhs[1]
#define PS_V		      	plhs[0] /* Output Arguments */
#define	RK_V						plhs[1]
#define BS_V						plhs[2]
#define T_CPU						plhs[3]

#define NV 25

int tm_bs(double *y,double *y0,double *dydt,double *fp,neuron_tm *nrnp,double dt){
	int flag=0,nv=6,bsk,kmax=50;
	double v,n,m,h,g_ampa,g_gaba,eta[6],tol = fp[17];
 	v = nrnp->v; n = nrnp->n; m = nrnp->m; h = nrnp->h; 
	g_ampa = nrnp->g_ampa; g_gaba = nrnp->g_gaba;
 	eta[0]=tol;	eta[1]=tol;	eta[2]=tol;	eta[3]=tol;	eta[4]=tol;	eta[5]=tol;
	fp[8]=nrnp->I;
  y[0]=nrnp->v; y[1]=nrnp->n; y[2]=nrnp->m; y[3]=nrnp->h;
  y[4]=nrnp->g_ampa; y[5]=nrnp->g_gaba;
  y0[0]=y[0]; y0[1]=y[1]; y0[2]=y[2]; y0[3]=y[3]; y0[4]=y[4]; y0[5]=y[5];
  tm_derivs(y,dydt,fp);
  bsk = bs_step(y,y0,dydt,nv,dt,fp,eta,tm_derivs); if(bsk==kmax)flag = 1;
  nrnp->v=y[0]; nrnp->n=y[1]; nrnp->m=y[2]; nrnp->h=y[3];
  nrnp->g_ampa = y[4]; nrnp->g_gaba = y[5];
	return flag;
}

void tm_rk(double *y,double *y0,double *dydt,double *fp,neuron_tm *nrnp,double dt){
	int nv=6;
	fp[8]=nrnp->I;
	y[0]=nrnp->v; y[1]=nrnp->n; y[2]=nrnp->m; y[3]=nrnp->h;
  y[4]=nrnp->g_ampa; y[5]=nrnp->g_gaba;
	y0[0]=y[0]; y0[1]=y[1]; y0[2]=y[2]; y0[3]=y[3];	y0[4]=y[4]; y0[5]=y[5];
	tm_derivs(y,dydt,fp);
	rk_step(y,y0,dydt,nv,fp,dt,tm_derivs);
	nrnp->v=y[0]; nrnp->n=y[1]; nrnp->m=y[2]; nrnp->h=y[3];
	nrnp->g_ampa = y[4]; nrnp->g_gaba = y[5];
}

int tm_ps(double **yp,double **co,double *yold,double *ynew,neuron_tm *nrnp,double *fp,double dt_full,int order_lim){
	int flag=0,p,nv=6;
	double v,n,m,h,g_ampa,g_gaba,a,b,c,d,co_Na,co_K,n2,n3,n4,m2,m3,m3h,e,f,g,q,r,s;
	double chi,psi,xi,v_alpha_n,v_alpha_m,v_beta_m,alpha_n,alpha_m,beta_m,eta[NV];
	double gNa = 20000, gK = 6000, gL = 10, Vr = -63;
	double E_alpha_n = Vr+15, E_beta_n = Vr+10, E_alpha_m = Vr+13;
	double E_beta_m = Vr+40, E_alpha_h = Vr+17, E_beta_h = Vr+40;
	double co_alpha_n = 0.032, co_alpha_m = 0.32, co_beta_m = 0.28;
	double tol = fp[17], *y, *y0, *dydt;
	double ka = 40.0, kb = 18.0, kc = 5.0, ke = 5.0, kf = 4.0, kg = 5.0; /*slopes*/
	v = nrnp->v; n = nrnp->n; m = nrnp->m; h = nrnp->h; 
	a = nrnp->a; b = nrnp->b;	c = nrnp->c; d = nrnp->d;
	g_ampa = nrnp->g_ampa; g_gaba = nrnp->g_gaba;
	
	/*Get higher power terms*/
	n2=n*n; n3=n2*n; n4=n3*n; m2=m*m; m3=m2*m; m3h=m3*h; co_K=gK*n4; co_Na=gNa*m3h;
	
	/*Tethered and derived values*/
	v_alpha_n = E_alpha_n-v; v_alpha_m = E_alpha_m-v; v_beta_m = v-E_beta_m; 
	c = exp((E_beta_h-v)/kc); d = h/(c+1); 
	e = exp(v_alpha_n/ke); q = v_alpha_n / (e-1); alpha_n = co_alpha_n*q;
 	f = exp(v_alpha_m/kf); r = v_alpha_m / (f-1);	alpha_m = co_alpha_m*r;
 	g = exp(v_beta_m/kg);  s = v_beta_m / (g-1); 	beta_m = co_beta_m*s;
 	chi = -co_K - co_Na - gL - g_ampa - g_gaba;	
 	psi = -(alpha_n+a);	xi = -(alpha_m+beta_m);
	fp[3]=alpha_n; fp[4]=alpha_m; fp[5]=co_K; fp[6]=co_Na; fp[8]=nrnp->I;
 	eta[0]=tol;	eta[1]=tol;	eta[2]=tol;	eta[3]=tol;	eta[4]=tol;	eta[5]=tol;
	
	/*Load variables into solver structure*/
 	yp[0][0] = v; yp[1][0] = n; yp[2][0] = m; yp[3][0] = h;
 	yp[4][0] = g_ampa; yp[5][0] = g_gaba;
 	yp[6][0] = a; yp[7][0] = b; yp[8][0] = c; yp[9][0] = d;
 	yp[10][0] = e; yp[11][0] = f; yp[12][0] = g; 
 	yp[13][0] = q; yp[14][0] = r; yp[15][0] = s;
 	yp[16][0] = chi; yp[17][0] = psi; yp[18][0] = xi;

 	/*Run generic PS solver*/
  p = ps_step(yp,co,yold,ynew,fp,eta,tm_first,tm_iter,0,order_lim,nv+2,nv);
 	
 	if(p<1 || p==order_lim){
 		y = malloc(nv*sizeof(double)); y0 = malloc(nv*sizeof(double));
    dydt = malloc(nv*sizeof(double));
		flag = tm_bs(y,y0,dydt,fp,nrnp,dt_full);
		/*Update subsidiary variables*/
		nrnp->a = 0.5*exp((E_beta_n-nrnp->v)/40);
		nrnp->b = 0.128*exp((E_alpha_h-nrnp->v)/18);
    free(y); free(y0); free(dydt);
	}
	else{		
    nrnp->v=ynew[0]; nrnp->n=ynew[1]; nrnp->m=ynew[2]; nrnp->h=ynew[3];
    nrnp->g_ampa=ynew[4]; nrnp->g_gaba=ynew[5]; nrnp->a=ynew[6]; nrnp->b=ynew[7];
	}
	return flag;
}

/***********************************************************/
void run_sim(double *ps_v,double *rk_v,double *bs_v,double *t_cpu,double *fp_in,int *ip_in){
  int syn_seed=ip_in[0],sim_type=ip_in[1],t_end=ip_in[2];
  int in_seed=ip_in[3],ps_only=ip_in[4],n_nrn=ip_in[5],order_lim=ip_in[99];
  double tol=fp_in[9], dt_rk=fp_in[10], dt_ps=fp_in[11]; 
  double dt,t,t_next,c0,c1,cp,dt_full,fp[100];
	double co_v,co_n,co_m,co_h,co_a,co_b,co_c,co_d,co_e,co_f,co_g,co_K,co_Na;
	neuron_tm *nrn, *nrnp, *nrnx; /*TM neuron pointers*/
	int step,i,p,nv=6,t_ms,flag,n_bs_fails=0,ip[100];  	
	double steps_rk=floor((1.0/dt_rk)+0.5), steps_ps=floor((1.0/dt_ps)+0.5);
  
  /*Cell Parameters*/
	double C = 200; /*pF*/
	double ka = 40.0, kb = 18.0, kc = 5.0, ke = 5.0, kf = 4.0, kg = 5.0;
	
	double tau_ampa_rk = 5.0/dt_rk, co_g_ampa_rk = -1.0/tau_ampa_rk; 
  double tau_gaba_rk = 10.0/dt_rk, co_g_gaba_rk = -1.0/tau_gaba_rk;	
	double tau_ampa_ps = 5.0/dt_ps, co_g_ampa_ps = -1.0/tau_ampa_ps;
	double tau_gaba_ps = 10.0/dt_ps, co_g_gaba_ps = -1.0/tau_gaba_ps;	
  
  double co_alpha_n = 0.032, co_alpha_m = 0.32, co_beta_m = 0.28;
  
  /*Dynamic Data structures for derivs code and generic PS solution*/
  double **co, **yp, *yold, *ynew, *y, *y0, *dydt;
  y = malloc(nv*sizeof(double));  
  y0 = malloc(nv*sizeof(double));
  dydt = malloc(nv*sizeof(double));
  yold = malloc(NV*sizeof(double));  
  ynew = malloc(NV*sizeof(double));
  yp = malloc(NV*sizeof(double *));  
  co = malloc(NV*sizeof(double *));
	for(i=0;i<NV;i++){
		yp[i] = malloc((order_lim+1)*sizeof(double));
		co[i] = malloc((order_lim+1)*sizeof(double));
	}
  nrn = malloc(n_nrn*sizeof(neuron_tm)); nrnx = nrn+n_nrn;
  
  /*Store constant parameters*/
	ip[0]=sim_type; fp[17] = tol;
	
	/*Set variable coefficients*/
	dt = dt_ps;
	co_v = dt/C; co_n = dt;	co_m = dt; co_h = dt; /*time rescaling version*/
	co_a = -1.0/ka;	co_b = -1.0/kb;	co_c = -1.0/kc;
	co_e = -1.0/ke;	co_f = -1.0/kf;	co_g = 1.0/kg;
	
	co[0][0] = co_v; co[1][0] = co_n; co[2][0] = co_m; co[3][0] = co_h;
	co[4][0] = co_g_ampa_ps; co[5][0] = co_g_gaba_ps;	
  co[6][0] = co_a; co[7][0] = co_b;	co[8][0] = co_c; 
  co[10][0] = co_e; co[11][0] = co_f; co[12][0] = co_g;

  /************************************************************/
  /************* Adaptive Parker-Sochacki Method **************/
  /************************************************************/
  mexPrintf("PS. ");
  fp[0] = dt_ps; fp[1] = co_g_ampa_ps; fp[2] = co_g_gaba_ps;
  for(nrnp = nrn; nrnp < nrnx; nrnp++){ /*Initialise neuron structure*/
		nrnp->v = fp_in[0]; nrnp->n = fp_in[1]; nrnp->m = fp_in[2];
		nrnp->h = fp_in[3];	nrnp->a = fp_in[4];	nrnp->b = fp_in[5];
		nrnp->c = fp_in[6];	nrnp->d = fp_in[7];	nrnp->I = fp_in[8];
		nrnp->g_ampa = 0.0; nrnp->g_gaba = 0.0;
	}
	dt_full=1; /*time rescaling*/ fp[7]=dt_full;
	for(p = 1; p < order_lim; p++){ 
		cp = 1.0/(double)(p+1);
  	co[0][p] = co[0][0]*cp; co[1][p] = co[1][0]*cp;	co[2][p] = co[2][0]*cp;
  	co[3][p] = co[3][0]*cp;	co[4][p] = co[4][0]*cp;	co[5][p] = co[5][0]*cp;
  	co[6][p] = co[6][0]*cp; co[7][p] = co[7][0]*cp;	co[8][p] = co[8][0]*cp;
  	co[10][p] = co[10][0]*cp; co[11][p] = co[11][0]*cp; co[12][p] = co[12][0]*cp;
	}
  c0 = (double)clock();
  for(t_ms=0,t=0; t_ms<t_end; t_ms++){
  	ps_v[t_ms] = nrn[0].v;
  	for(step=0; step<steps_ps; step++){
  		t_next = (double)t_ms + (step+1)*dt;/*end of current time step*/
  		for(nrnp = nrn; nrnp < nrnx; nrnp++){ /*loop over neurons*/
        fp[99] = dt_full;
        flag = tm_ps(yp,co,yold,ynew,nrnp,fp,dt_full,order_lim);
		  } /* end loop over neurons*/
		  t=t_next;
  	} /*loop over steps*/
  } /*loop over t_ms*/
  c1 = (double)clock();
  t_cpu[0] = (double)(c1 - c0)/(CLOCKS_PER_SEC);
  mexPrintf("Time = %5.2f. \n",t_cpu[0]); fflush(stdout);
	
	/************************************************************/
  /********************* Bulirsch-Stoer ***********************/
  /************************************************************/
  mexPrintf("BS. ");
  for(nrnp = nrn; nrnp < nrnx; nrnp++){ /*Initialise neuron structure*/
		nrnp->v = fp_in[0]; nrnp->n = fp_in[1]; nrnp->m = fp_in[2]; nrnp->h = fp_in[3];
		nrnp->g_ampa = 0.0; nrnp->g_gaba = 0.0; nrnp->I = fp_in[8];
	}
  c0 = (double)clock();
  for(t_ms=0,t=0; t_ms<t_end; t_ms++){
  	if(ps_only==1) break;
  	bs_v[t_ms] = nrn[0].v;
  	for(step=0; step<steps_ps; step++){
  		t_next = (double)t_ms + (step+1)*dt;/*end of current time step*/
	  	for(nrnp = nrn; nrnp < nrnx; nrnp++){ /*loop over neurons*/
	  		flag = tm_bs(y,y0,dydt,fp,nrnp,1); n_bs_fails+=flag;
		  } /*end loop over neurons*/
		  t=t_next;
  	}/*loop over steps*/
  } /*loop over t_ms*/
  c1 = (double)clock();
  t_cpu[1] = (double)(c1 - c0)/(CLOCKS_PER_SEC);
	mexPrintf("Time = %5.2f. \t",t_cpu[1]);
	mexPrintf("n BS fails = %d.\n",n_bs_fails); fflush(stdout);
	  
  /************************************************************/
  /*************** 4th-order Runge-Kutta Method ***************/
  /************************************************************/
  mexPrintf("RK. ");
	fp[0] = dt_rk; fp[1] = co_g_ampa_rk; fp[2] = co_g_gaba_rk;
  for(nrnp = nrn; nrnp < nrnx; nrnp++){ /*Initialise neuron structure*/
		nrnp->v = fp_in[0]; nrnp->n = fp_in[1]; nrnp->m = fp_in[2]; nrnp->h = fp_in[3];
		nrnp->g_ampa = 0.0; nrnp->g_gaba = 0.0; nrnp->I = fp_in[8];
	}
  c0 = (double)clock();
  for(t_ms=0,t=0; t_ms<t_end; t_ms++){
  	if(ps_only==1) break;
  	rk_v[t_ms] = nrn[0].v;
  	for(step=0; step<steps_rk; step++){
  		t_next = (double)t_ms + (step+1)*dt_rk;/*end of current time step*/
	  	for(nrnp = nrn; nrnp < nrnx; nrnp++){ /*loop over neurons*/
	  		tm_rk(y,y0,dydt,fp,nrnp,1);
		  } /*end loop over neurons*/
		  t=t_next;
  	} /*loop over steps*/
  } /*loop over t_ms*/
  c1 = (double)clock();
  t_cpu[2] = (double)(c1 - c0)/(CLOCKS_PER_SEC);
  mexPrintf("Time = %5.2f. \n",t_cpu[2]); fflush(stdout);
	
	/*Free dynamic arrays*/
  free(nrn); free(yold); free(ynew); free(y); free(y0); free(dydt);
	for(i=0;i<NV;i++){free(yp[i]); free(co[i]);} free(yp); free(co); 
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
  double *fp, *ps_v, *rk_v, *bs_v, *t_cpu, dt;
  int t_end,*ip;  
  
  /*  Check for proper number of arguments. */
  if (nrhs != 2) mexErrMsgTxt("2 inputs required.");
  if (nlhs != 4) mexErrMsgTxt("4 outputs required.");
	
  /* Get the inputs */ 
  fp = (double *)mxGetData(FP);
  ip = (int *)mxGetData(IP); t_end=ip[2];

  /* Create output arguments and assign pointers to them */
  PS_V = mxCreateDoubleMatrix(t_end,1,mxREAL); ps_v = (double *)mxGetData(PS_V);
  RK_V = mxCreateDoubleMatrix(t_end,1,mxREAL); rk_v = (double *)mxGetData(RK_V);
  BS_V = mxCreateDoubleMatrix(t_end,1,mxREAL); bs_v = (double *)mxGetData(BS_V);
  T_CPU = mxCreateDoubleMatrix(3,1,mxREAL); t_cpu = (double *)mxGetData(T_CPU);

  /* Call the C subroutine. */
  run_sim(ps_v,rk_v,bs_v,t_cpu,fp,ip);
}
