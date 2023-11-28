/* =============================================================
 * iz_ps.c Written by Dr Robert Stewart for Stewart & Bair, 2009
 * Benchmark model based on Brette et al 2007 but with Izhikevich neurons
 * Features continuous spike times and synaptic event delivery times
 * Compile:   mex iz_ps.c iz_util.c integ_util.c
 * Call as:   [RK_tf,RK_nrn,PS_tf,PS_nrn,BS_tf,BS_nrn,t_cpu,RK_V,PS_V,BS_V,
 * i_stats,f_stats] = iz_ps(fp,ip);
 * =============================================================*/
#include "mex.h"
#include "matrix.h"
#include <math.h> 
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "iz_util.h"

#define FP      					prhs[0] /* Input Arguments */
#define IP								prhs[1]
#define RK_TF		  				plhs[0] /* Output Arguments */
#define RK_NRN						plhs[1]
#define PS_TF							plhs[2]
#define PS_NRN						plhs[3]
#define BS_TF							plhs[4]
#define BS_NRN						plhs[5]
#define T_CPU							plhs[6] 
#define RK_V							plhs[7]
#define PS_V							plhs[8]
#define BS_V							plhs[9]
#define I_STATS						plhs[10]
#define F_STATS						plhs[11]
#define NR_TOL (1e-20)	/*Newton-Raphson error tolerance*/
#define NV 10

void get_local_inputs(neuron_iz *nrn,int *ne,double *te,synapse **syn,int n_nrn){
  synapse *synp; int i,j,k,n_in;
  for(i = 0; i < n_nrn; i++){nrn[i].n_in = 0;} /*Clear input buffers*/
  for(i = 0; i < n_nrn; i++){
    if(--ne[i]==0){
      for(synp = syn[i]; (synp->n) < n_nrn; synp++){
        k = synp->n;
        n_in = nrn[k].n_in;
        if(n_in<MAX_IN){
	        j=n_in; /*Use insertion sort to maintain ordered synaptic events*/
	        while ((j > 0) && (nrn[k].in_t[j-1] > te[i])){
	          nrn[k].in_t[j] = nrn[k].in_t[j-1]; /*shift*/
	          nrn[k].in_w[j] = nrn[k].in_w[j-1]; j--;
	        }
	        nrn[k].in_t[j] = te[i]; nrn[k].in_w[j] = synp->w; nrn[k].n_in++;
        }
        else printf("Input overflow detected\n");
      }
    }
  }
}

/*Main stepper routine for RK method on Izhikevich neuron - runs full step*/
int iz_rk(double *y, double *y0, double *dydt, double *RK_tf, int *RK_nrn,
	neuron_iz *nrnp, int *ne, double *te, int *ip, double *fp, int *fcount,
	unsigned long long int *icount, FILE *trace_ptr){
  double v,u,g_ampa,g_gaba,E_ampa, E_gaba, co_g_ampa, co_g_gaba, t, start;
	double l,k,a,b,E,I,dv,dx,dx_old,dt_part,dt_full,dt,trace_val;
  int i,flag=0,t_ms,step,sim_type,nrn_ind,steps,nv=4,trace_ind,trace_rec;
	sim_type = ip[0];steps = ip[1];step = ip[2];t_ms = ip[3];nrn_ind = ip[4];
	trace_ind = ip[5]; trace_rec = ip[6];
	co_g_ampa = fp[8];co_g_gaba = fp[9];
	dt = fp[12]; dt_full = fp[13]; t = fp[15]; start = fp[16];
	
  /*extract variables from neuron structure*/
  v = nrnp->v; u = nrnp->u; g_ampa = nrnp->g_ampa; g_gaba = nrnp->g_gaba;
  E_ampa = nrnp->E_ampa; E_gaba = nrnp->E_gaba; k = nrnp->k; l = nrnp->l;
  I = nrnp->I; E = nrnp->E; a = nrnp->a; b = nrnp->b; 

  /*Generic update*/
  fp[0] = I; fp[1] = k; fp[2] = l; fp[3] = E_ampa; fp[4] = E_gaba;
	fp[5] = E; fp[6] = a;fp[7] = b; 
  y0[0] = v; y0[1] = u; y0[2] = g_ampa; y0[3] = g_gaba;
  y[0] = v; y[1] = u; y[2] = g_ampa; y[3] = g_gaba;
  iz_derivs(y, dydt, fp);
  rk_step(y,y0,dydt,nv,fp,dt_full,iz_derivs); (*icount)++;
  if(trace_rec==1){
    /*Full trace output - not collected for current injection sims*/
    if(nrn_ind==trace_ind){
      fwrite(&t,8,1,trace_ptr);fwrite(&v,8,1,trace_ptr);fwrite(&u,8,1,trace_ptr);
      fwrite(&g_ampa,8,1,trace_ptr); fwrite(&g_gaba,8,1,trace_ptr);     
    }
  }

  if (y[0] >= nrnp->v_peak){
    /*Find accurate spike time using Newton-Raphson*/
    dt_part = (nrnp->v_peak-y[0])/dydt[0];	/*First step*/
    dx_old = 100.0;
    for (i = 0; i<20; i++){ /*Up to 20 iterations */
      rk_step(y,y0,dydt,nv,fp,dt_part,iz_derivs); (*icount)++;
      v=y[0];u=y[1];g_ampa=y[2];g_gaba=y[3];
      dv = E*(k*v*v + l*v - u + I - g_ampa*(v-E_ampa)-g_gaba*(v-E_gaba));
      dx = (v-nrnp->v_peak)/dv; dt_part -= dx; if(fabs(dx)<NR_TOL)break;
      if(fabs(dx+dx_old)<NR_TOL)break; dx_old=dx;/*For oscillations*/
    }
    /*if(i==20) mexPrintf("RK NR error tolerance failure. dx: %e\n",dx);*/
    if(dt_part>dt_full || dt_part<0){printf("RK - NR divergence.\n"); flag = 1; dt_part=dt_full/2;}
    /*Record spike and schedule events*/
    RK_tf[*fcount] = t+dt*dt_part; RK_nrn[*fcount] = nrn_ind+1; (*fcount)++;
		if(sim_type==2){ne[nrn_ind]=steps; te[nrn_ind]=start+dt_part;}
    
    if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){
        trace_val = t+dt*dt_part; fwrite(&trace_val,8,1,trace_ptr);
        fwrite(&v,8,1,trace_ptr); fwrite(&u,8,1,trace_ptr);
        fwrite(&g_ampa,8,1,trace_ptr); fwrite(&g_gaba,8,1,trace_ptr);     
      }
    }

    /*Post spike update*/
    y[0] = nrnp->v_reset; y[1] += nrnp->u_step; y[2]=g_ampa; y[3]=g_gaba;
    y0[0] = y[0]; y0[1] = y[1]; y0[2] = y[2]; y0[3] = y[3];
    if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){
        trace_val = t+dt*dt_part; fwrite(&trace_val,8,1,trace_ptr); 
        fwrite(&y[0],8,1,trace_ptr); fwrite(&y[1],8,1,trace_ptr);
        fwrite(&y[2],8,1,trace_ptr); fwrite(&y[3],8,1,trace_ptr);     
      }
    }
    dt_part = dt_full-dt_part;
    iz_derivs(y, dydt, fp);
    rk_step(y,y0,dydt,nv,fp,dt_part,iz_derivs); (*icount)++;
  }
  nrnp->v = y[0];nrnp->u = y[1];nrnp->g_ampa = y[2]; nrnp->g_gaba = y[3];
	return flag;
}

/*Main stepper routine for BS method on Izhikevich neuron - runs full step*/
int iz_bs(double *y, double *y0, double *dydt, double *BS_tf, int *BS_nrn,
	neuron_iz *nrnp, int *ne, double *te, int *ip, double *fp, int *fcount, 
	unsigned long long int *icount, double *mu, int *max_bsk, FILE *trace_ptr){
  double v,u,g_ampa,g_gaba,E_ampa,E_gaba,co_g_ampa,co_g_gaba,t,start,tol;
	double l,k,a,b,E,I,dv,dx,dx_old,dt_part,dt_full,dt,eta[4],trace_val;
  int i,flag=0,bsk,t_ms,step,sim_type,nrn_ind,steps,nv=4,kmax=50,trace_ind,trace_rec;
	sim_type = ip[0];steps = ip[1];step = ip[2];t_ms = ip[3];nrn_ind = ip[4];
	trace_ind = ip[5]; trace_rec = ip[6];
	co_g_ampa = fp[8];co_g_gaba = fp[9];
	dt = fp[12]; dt_full = fp[13]; t = fp[15]; start = fp[16]; tol = fp[17];
  /*extract variables from neuron structure*/
  v = nrnp->v; u = nrnp->u; g_ampa = nrnp->g_ampa; g_gaba = nrnp->g_gaba;
  E_ampa = nrnp->E_ampa; E_gaba = nrnp->E_gaba; k = nrnp->k; l = nrnp->l;
  I = nrnp->I; E = nrnp->E; a = nrnp->a; b = nrnp->b;

  /*Store fp*/
  fp[0] = I; fp[1] = k; fp[2] = l; fp[3] = E_ampa; fp[4] = E_gaba;
  fp[5] = E; fp[6] = a;fp[7] = b;
  y0[0] = v; y0[1] = u; y0[2] = g_ampa; y0[3] = g_gaba;
  
  /*Set error tolerance */
  eta[0] = tol;eta[1] = tol;eta[2] = tol;eta[3] = tol;

  /*Generic update code*/
  y[0] = v; y[1] = u; y[2] = g_ampa; y[3] = g_gaba;
  iz_derivs(y, dydt, fp);
  bsk = bs_step(y,y0,dydt,nv,dt_full,fp,eta,iz_derivs); (*icount)++;
  *mu += ((double)bsk- *mu)/(double)*icount; if(bsk==kmax)flag = 1; 
  if(bsk > *max_bsk)*max_bsk = bsk;
  if(trace_rec==1){
    /*Full trace output - not collected for current injection sims*/
    if(nrn_ind==trace_ind){
      fwrite(&t,8,1,trace_ptr);fwrite(&v,8,1,trace_ptr);fwrite(&u,8,1,trace_ptr);
      fwrite(&g_ampa,8,1,trace_ptr); fwrite(&g_gaba,8,1,trace_ptr); 
      trace_val = (double)bsk; fwrite(&trace_val,8,1,trace_ptr); 
    }
  }

  if (y[0] >= nrnp->v_peak){
    /*Find accurate spike time using Newton-Raphson*/
    dt_part = (nrnp->v_peak-y[0])/dydt[0];	/*First step*/
    dx_old = 100.0;
    for (i = 0; i<20; i++){ /*Up to 20 NR iterations */
      bsk = bs_step(y,y0,dydt,nv,dt_part,fp,eta,iz_derivs); (*icount)++;
      *mu += ((double)bsk- *mu)/(double)*icount; if(bsk==kmax)flag = 1;
      if(bsk > *max_bsk)*max_bsk = bsk;
      v=y[0];u=y[1];g_ampa=y[2];g_gaba=y[3];
      dv = E*(k*v*v + l*v - u + I - g_ampa*(v-E_ampa)-g_gaba*(v-E_gaba));
      dx = (v-nrnp->v_peak)/dv;
      dt_part -= dx;
      if(fabs(dx)<NR_TOL)break;
      if(fabs(dx+dx_old)<NR_TOL)break; dx_old=dx;/*For oscillations*/
    }
    /*if(i==20) mexPrintf("BS NR error tolerance failure. dx: %e\n",dx);*/
    if(dt_part>dt_full || dt_part<0){printf("BS - NR divergence.\n"); flag = 1; dt_part=dt_full/2;}				
    /*Record spike and schedule events*/
    BS_tf[*fcount] = t+dt*dt_part; BS_nrn[*fcount] = nrn_ind+1; (*fcount)++;
		if(sim_type==2){ne[nrn_ind]=steps; te[nrn_ind] = start+dt_part;}
		if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){
        trace_val = t+dt*dt_part; fwrite(&trace_val,8,1,trace_ptr);
        fwrite(&v,8,1,trace_ptr); fwrite(&u,8,1,trace_ptr);
        fwrite(&g_ampa,8,1,trace_ptr); fwrite(&g_gaba,8,1,trace_ptr); 
        trace_val = (double)bsk; fwrite(&trace_val,8,1,trace_ptr); 
      }
    }
    
    /*Post spike update*/
    y[0] = nrnp->v_reset; y[1] += nrnp->u_step; y[2]=g_ampa; y[3]=g_gaba;
    y0[0] = y[0]; y0[1] = y[1]; y0[2] = y[2]; y0[3] = y[3];
    if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){
        trace_val = t+dt*dt_part; fwrite(&trace_val,8,1,trace_ptr);
        fwrite(&y[0],8,1,trace_ptr); fwrite(&y[1],8,1,trace_ptr);
        fwrite(&y[2],8,1,trace_ptr); fwrite(&y[3],8,1,trace_ptr);     
      }
    }
		dt_part=dt_full-dt_part;
    iz_derivs(y, dydt, fp);
    bsk = bs_step(y,y0,dydt,nv,dt_part,fp,eta,iz_derivs); (*icount)++;
    *mu += ((double)bsk- *mu)/(double)*icount; if(bsk==kmax)flag = 1;
    if(bsk > *max_bsk)*max_bsk = bsk;
    if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){trace_val = (double)bsk;fwrite(&trace_val,8,1,trace_ptr);}		
    }
  }
  nrnp->v = y[0];nrnp->u = y[1];nrnp->g_ampa = y[2]; nrnp->g_gaba = y[3];
	return flag;
}

/*Main stepper routine for PS method on Izhikevich neuron - runs full step*/
int iz_ps(double **yp, double **co, double *yold, double *ynew, double *PS_tf, int *PS_nrn,
	neuron_iz *nrnp, int *ne, double *te, int *ip, double *fp, int *fcount,
	unsigned long long int *icount, double *mu, int *max_order, FILE *trace_ptr, int order_lim){
	double v,u,g_ampa,g_gaba,chi,E_ampa,E_gaba,co_g_ampa,co_g_gaba,t,start,tol;
	double l,k,a,b,E,I,vnew,dv,dx,dx_old,dt_part,dt_full,dt,eta[4],trace_val; 
	int ps_order,i,j,flag=0,t_ms,step,sim_type,nrn_ind,steps,trace_ind,trace_rec;
	int err_nv=4, nv=5;
	sim_type = ip[0];steps = ip[1];step = ip[2];t_ms = ip[3];nrn_ind = ip[4];
	trace_ind = ip[5]; trace_rec = ip[6];
	co_g_ampa = fp[8];co_g_gaba = fp[9];/*decay_ampa = fp[10];decay_gaba = fp[11];*/
	dt = fp[12]; dt_full = fp[13]; t = fp[15]; start = fp[16]; tol = fp[17];
		
	/*extract variables from neuron structure*/
	v = nrnp->v; u = nrnp->u; g_ampa = nrnp->g_ampa; g_gaba = nrnp->g_gaba;
	E_ampa = nrnp->E_ampa; E_gaba = nrnp->E_gaba; k = nrnp->k; l = nrnp->l;
	I = nrnp->I; E = nrnp->E; a = nrnp->a; b = nrnp->b;
	chi = k*v - g_ampa - g_gaba + nrnp->l;
	
	/*Set error tolerance */
  eta[0] = tol;eta[1] = tol;eta[2] = tol;eta[3] = tol;

	yp[0][0] = v;yp[1][0] = u;yp[2][0] = g_ampa;yp[3][0] = g_gaba;yp[4][0] = chi;
	fp[0] = I; fp[1] = k; fp[3] = E_ampa; fp[4] = E_gaba;
	fp[5] = E; fp[6] = a;fp[7] = b; fp[99] = dt_full;
  ps_order = ps_step(yp,co,yold,ynew,fp,eta,iz_first,iz_iter,1,order_lim,nv,err_nv); (*icount)++; /*integrated step function*/  

	*mu += ((double)ps_order- *mu)/(double)*icount; 
	if(ps_order > *max_order)*max_order = ps_order;
  vnew=ynew[0]; /*New membrane voltage value*/
	if(trace_rec==1){
    /*Full trace output - not collected for current injection sims*/
    if(nrn_ind==trace_ind){
      fwrite(&t,8,1,trace_ptr); fwrite(&v,8,1,trace_ptr);
      fwrite(&u,8,1,trace_ptr); fwrite(&g_ampa,8,1,trace_ptr); 
      fwrite(&g_gaba,8,1,trace_ptr); trace_val = (double)ps_order;    
      fwrite(&trace_val,8,1,trace_ptr); 
    }
  }
	if (vnew >= nrnp->v_peak){ /*rare*/
		yp[0][0] = v - nrnp->v_peak; /*shifted for root finding*/
		dt_part = -yp[0][0]/yp[0][1];	/*First step*/
    dx_old = 100.0;
		for (i = 0; i<20; i++){ /*Up to 20 NR iterations */
			vnew=yp[0][ps_order]*dt_part+yp[0][ps_order-1];
			dv=yp[0][ps_order];
			for(j=ps_order-2;j>=0;j--){
				dv = vnew + dv*dt_part;
				vnew=yp[0][j]+vnew*dt_part;
			}
			dx = vnew/dv; dt_part -= dx; if(fabs(dx)<NR_TOL)break;
      if(fabs(dx+dx_old)<NR_TOL)break; dx_old=dx;/*For oscillations*/
		}
    /*if(i==20) mexPrintf("PS NR error tolerance failure. dx: %e\n",dx);*/
		if(dt_part>dt_full || dt_part<0){printf("PS - NR divergence.\n"); flag = 1; dt_part=dt_full/2;}
		/*Record spike and schedule events*/
		PS_tf[*fcount] = t+dt*dt_part; PS_nrn[*fcount] = nrn_ind+1; (*fcount)++;
		if(sim_type==2){ne[nrn_ind]=steps; te[nrn_ind]=start+dt_part;} 
		/*Evaluate u, g_ampa, g_gaba at corrected spike time*/
		ps_update(yp,1,ps_order,dt_part,&u);
		ps_update(yp,2,ps_order,dt_part,&g_ampa);
		ps_update(yp,3,ps_order,dt_part,&g_gaba);
		if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){
        trace_val = t+dt*dt_part; fwrite(&trace_val,8,1,trace_ptr);
        fwrite(&v,8,1,trace_ptr); fwrite(&u,8,1,trace_ptr);
        fwrite(&g_ampa,8,1,trace_ptr); fwrite(&g_gaba,8,1,trace_ptr); 
        trace_val = (double)ps_order; fwrite(&trace_val,8,1,trace_ptr); 
      }
    }
		
		/*post spike updates*/
		v = nrnp->v_reset; u += nrnp->u_step; chi = k*v - g_ampa - g_gaba + nrnp->l;
		yp[0][0] = v;yp[1][0] = u;yp[2][0] = g_ampa;yp[3][0] = g_gaba;yp[4][0] = chi;
		if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){
        trace_val = t+dt*dt_part; fwrite(&trace_val,8,1,trace_ptr);
        fwrite(&v,8,1,trace_ptr); fwrite(&u,8,1,trace_ptr);
        fwrite(&g_ampa,8,1,trace_ptr); fwrite(&g_gaba,8,1,trace_ptr); 
      }
    }
		
		dt_part = dt_full-dt_part; fp[99] = dt_part;
    ps_order = ps_step(yp,co,yold,ynew,fp,eta,iz_first,iz_iter,1,order_lim,nv,err_nv); (*icount)++; /*new integrated step function*/  
		*mu += ((double)ps_order- *mu)/(double)*icount;
		if(ps_order > *max_order)*max_order = ps_order;
    if(trace_rec==1){
      /*Full trace output - not collected for current injection sims*/
      if(nrn_ind==trace_ind){trace_val = (double)ps_order;fwrite(&trace_val,8,1,trace_ptr);}
    }
		nrnp->v=ynew[0]; nrnp->u=ynew[1]; nrnp->g_ampa=ynew[2]; nrnp->g_gaba=ynew[3];
	}
	else{
    nrnp->v=ynew[0]; nrnp->u=ynew[1]; nrnp->g_ampa=ynew[2]; nrnp->g_gaba=ynew[3];
	}	return flag;
}


/***********************************************************/
void run_sim(double *RK_tf, int *RK_nrn, double *PS_tf, int *PS_nrn, 
	double *BS_tf, int *BS_nrn, double *t_cpu, double *RK_v, double *PS_v, double *BS_v,
	int *i_stats, double *f_stats, double *fp_in, int *ip_in){
	int n_nrn=ip_in[0], sim_type=ip_in[1], t_end=ip_in[2], syn_seed=ip_in[3];
	int in_seed=ip_in[4], cond=ip_in[5], ps_only=ip_in[6], n_ex = ip_in[7];
	int trace_rec=ip_in[8], order_lim=ip_in[99];
	double inj=fp_in[0], tol=fp_in[1], dt_rk=fp_in[2], dt_ps=fp_in[3]; 
	double C=fp_in[4], vr=fp_in[5], vt=fp_in[6], k=fp_in[7], a=fp_in[8]; 
	double b=fp_in[9], v_reset=fp_in[10], u_step=fp_in[11], v_peak=fp_in[12];
	double tau_ampa=fp_in[13], tau_gaba=fp_in[14], E_ampa=fp_in[15], E_gaba=fp_in[16];
  float w_ampa=(float)fp_in[17]; float w_gaba=(float)fp_in[18];
	double rand_inj_max=fp_in[19], p_connect=fp_in[20];  
	double E=1.0/C, w, t, start; /*Electrical elastance = 1/C*/
	double syn_test, fp[100];
  unsigned int *temp_n; float *temp_w;
  neuron_iz *nrn, *nrnp, *nrnx;/*Izhikevich neuron pointers*/
  synapse **syn; /*Dynamic synapse structure*/
  
	int nrn_ind,fcount,flag,substep,*ne,ip[100],trace_ind=4,nv=4,ps_nv=5;
	int i,j,p,n_syn,step,t_ms,n_in,max_order,max_bsk,trace_n_in=0,n_bs_fails;
	unsigned long long int icount;
	double trace_in_t,trace_in_w,mu_order,mu_bsk;

	clock_t c0_rk, c0_ps, c0_bs;
  double steps_rk=floor((1.0/dt_rk)+0.5), steps_ps=floor((1.0/dt_ps)+0.5);
	
	/*Scale time constants to time step size*/
	double tau_ampa_rk = tau_ampa/dt_rk;
	double tau_gaba_rk = tau_gaba/dt_rk;
	double co_g_ampa_rk = -1.0/tau_ampa_rk, co_g_gaba_rk = -1.0/tau_gaba_rk;	
	
	double tau_ampa_ps = tau_ampa/dt_ps;
	double tau_gaba_ps = tau_gaba/dt_ps;
	double co_g_ampa_ps = -1.0/tau_ampa_ps, co_g_gaba_ps = -1.0/tau_gaba_ps;	
	
	FILE *ps_trace_ptr, *rk_trace_ptr, *bs_trace_ptr, *trace_in_t_ptr, *trace_in_w_ptr;
  char file_stub[100] = "iz_bench";
  char syn_seed_str[10], in_seed_str[10], n_nrn_str[10];
  char cond_str[10];
  char ps_trace_str[100], bs_trace_str[100], rk_trace_str[100];
  char trace_in_t_str[100], trace_in_w_str[100];
  
  /*Dynamic Data structures for derivs code and generic PS solution*/
  double **co, **yp, *yold, *ynew, *y, *y0, *dydt, *te, *rand_inj;
  y = malloc(nv*sizeof(double));  
  y0 = malloc(nv*sizeof(double));
  dydt = malloc(nv*sizeof(double));
  yold = malloc(NV*sizeof(double));  
  ynew = malloc(NV*sizeof(double));
  ne = malloc(n_nrn*sizeof(int));
  te = malloc(n_nrn*sizeof(double));
  rand_inj = malloc(n_nrn*sizeof(double));
  yp = malloc(NV*sizeof(double *));  
  co = malloc(NV*sizeof(double *));
	for(i=0;i<NV;i++){
		yp[i] = malloc((order_lim+1)*sizeof(double));
		co[i] = malloc((order_lim+1)*sizeof(double));
	}
  nrn = malloc(n_nrn*sizeof(neuron_iz)); nrnx = nrn+n_nrn;
  
  /*Store constant parameters*/
	ip[0]=sim_type; ip[5]=trace_ind; ip[6]=trace_rec; fp[17]=tol;
  
	if(sim_type == 0){
		strcat(file_stub, "_inj");
	}
	else if(sim_type == 1){
		strcat(file_stub, "_ext");
	}
	else if(sim_type == 2){
		strcat(file_stub, "_fix");
	}
	sprintf(syn_seed_str,"_%d",syn_seed);
	sprintf(in_seed_str,"_%d",in_seed);
	sprintf(n_nrn_str,"_%d",n_nrn); /*include n neurons in filename*/ 
	strcat(file_stub,n_nrn_str);
	strcat(file_stub,syn_seed_str); strcat(file_stub,in_seed_str);
	sprintf(cond_str,"_%d",cond);
	strcat(file_stub,cond_str);
	printf( "%s\n",file_stub);
	strcpy(ps_trace_str,file_stub);	strcat(ps_trace_str,"_ps_trace");
	strcpy(bs_trace_str,file_stub);	strcat(bs_trace_str,"_bs_trace");
	strcpy(rk_trace_str,file_stub);	strcat(rk_trace_str,"_rk_trace");
	strcpy(trace_in_t_str,file_stub); strcat(trace_in_t_str,"_trace_in_t");
	strcpy(trace_in_w_str,file_stub); strcat(trace_in_w_str,"_trace_in_w");
  if(trace_rec==1){
    if((ps_trace_ptr = fopen(ps_trace_str,"wb"))==NULL) mexErrMsgTxt("Can't open ps_trace file");
    if((bs_trace_ptr = fopen(bs_trace_str,"wb"))==NULL) mexErrMsgTxt("Can't open bs_trace file");
    if((rk_trace_ptr = fopen(rk_trace_str,"wb"))==NULL) mexErrMsgTxt("Can't open rk_trace file");
    if(cond == 4){
      if((trace_in_t_ptr = fopen(trace_in_t_str,"wb"))==NULL) mexErrMsgTxt("Can't open trace_in_t file");
      if((trace_in_w_ptr = fopen(trace_in_w_str,"wb"))==NULL) mexErrMsgTxt("Can't open trace_in_w file");
    }
  }
	
	/*Specify connectivity for Benchmark network from Brette et al 2007*/
	/*srand48((long)(1234567*syn_seed)); /*initialise drand48 rng - non-windows*/
  srand(1234567*syn_seed);/*initialise rand rng - portable*/
	/*synapse *syn[n_nrn];*/
  syn = malloc(n_nrn*sizeof(synapse *));
  temp_n = malloc(n_nrn*sizeof(unsigned int));
  temp_w = malloc(n_nrn*sizeof(float));
	for(i=0;i<n_nrn;i++){
		n_syn = 0;
		for(j=0;j<n_nrn;j++){
			if(i!=j){ /*no auto-synapses*/
				/*syn_test = drand48(); /*drand48 version - non-windows*/
        syn_test = ((double)rand()/((double)RAND_MAX)); /*ANSI rand version*/
				if(syn_test<p_connect){ /*Form a synapse*/
					temp_n[n_syn] = j; /*Record target neuron*/
					if(i<n_ex)temp_w[n_syn] = w_ampa; /*excitatory synapses*/
					else temp_w[n_syn] = w_gaba; /*inhibitory synapses*/
					n_syn++;
				}				
			}
		}
		syn[i] = malloc((n_syn+1)*sizeof(synapse));
		for(j = 0; j < n_syn; j++){
			syn[i][j].n = temp_n[j];
			syn[i][j].w = temp_w[j];
		}
		syn[i][n_syn].n = n_nrn; /*Dummy synapse*/
	}
  free(temp_n);free(temp_w);
	/*random initial injection current values*/
	/*srand48((long)(1234567*in_seed)); /*initialise drand48 rng - non-windows*/
  srand(1234567*in_seed);/*initialise rand rng - portable*/
	/*if(sim_type>0){for(i=0;i<n_nrn;i++)rand_inj[i] = drand48()*rand_inj_max;}*/
  if(sim_type>0){for(i=0;i<n_nrn;i++)rand_inj[i] = rand_inj_max*((double)rand()/((double)RAND_MAX));} /*ANSI rand version*/
	
	/*Initialise neuron structure, with voltages shifted relative to vr*/
	for(nrnp = nrn; nrnp < nrnx; nrnp++){ 
		nrnp->vr = vr; nrnp->k = k; nrnp->l = -k*(vt-vr); nrnp->b = b; 
		nrnp->v_reset = v_reset-vr; nrnp->u_step = u_step; nrnp->v_peak = v_peak-vr; 
		nrnp->E_ampa = E_ampa-vr; nrnp->E_gaba = E_gaba-vr;
 	} 
	
	/*** Run simulations ***/
  /************************************************************/
  /********************** Runge-Kutta 4 ***********************/
  /************************************************************/
  fp[8]=co_g_ampa_rk; fp[9]=co_g_gaba_rk; fp[12]=dt_rk; ip[1]=(int)steps_rk;
  /*Scale time/rate constants such that dt=1 in the equations*/  
 	for(nrnp = nrn; nrnp < nrnx; nrnp++){nrnp->E=E*dt_rk; nrnp->a = a*dt_rk;}
  mexPrintf("RK\n");
  for(nrnp = nrn,i=0; nrnp < nrnx; nrnp++){ 
  	nrnp->v=0; nrnp->u=0; nrnp->g_gaba=0; nrnp->g_ampa=0;
    if(sim_type>0) nrnp->I = rand_inj[i];
    else nrnp->I = inj;
    nrn[i].n_in = 0;ne[i] = -1; i++;
  }
  fcount=0; icount=0; flag=0;  c0_rk = clock();
  for(t_ms=0; t_ms<t_end; t_ms++){
  	if(ps_only==1) break;
  	RK_v[t_ms] = nrn[trace_ind].v;
    if(sim_type>0 && t_ms==50)for(nrnp=nrn;nrnp<nrnx;nrnp++)nrnp->I=0;
  	for(step=0; step<steps_rk; step++){
  		t = (double)t_ms + (double)step*dt_rk;
      if(sim_type == 2)get_local_inputs(nrn,ne,te,syn,n_nrn);  
			for(nrnp=nrn, nrn_ind=0; nrnp<nrnx; nrnp++, nrn_ind++){
				ip[2] = step; ip[3] = t_ms;	ip[4] = nrn_ind; 
        fp[15] = t; /*real time at start of step*/
        /*Work through substeps separated by synaptic events*/
        if(sim_type == 2){
	        n_in = nrnp->n_in;
          start = 0; fp[16] = start; /*start time of substep (in [0 1] of whole step)*/
	        if(n_in){
	          for(substep=0;substep<n_in;substep++){ /*one substep per event*/
	          	fp[13] = nrnp->in_t[substep] - start; 
	          	if(fp[13]>0){
	          		flag = iz_rk(y,y0,dydt,RK_tf,RK_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,rk_trace_ptr);
              	start = nrnp->in_t[substep]; fp[16] = start; fp[15]+=dt_rk*fp[13];
              	/*fp[15] = (double)t_ms + dt_rk*( (double)step + start);*/
	          	}
	          	w = nrnp->in_w[substep];
	          	if(w > 0) nrnp->g_ampa+=w; /*AMPA*/
	          	else nrnp->g_gaba-=w; /*GABA*/
	          }
	          fp[13] = 1-start; /*remainder of time step*/
	          flag = iz_rk(y,y0,dydt,RK_tf,RK_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,rk_trace_ptr);
	        }
	        else{
						fp[13] = 1;  				
		        flag = iz_rk(y,y0,dydt,RK_tf,RK_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,rk_trace_ptr);
	        }
        }
        else{
          fp[13] = 1;  				
	        flag = iz_rk(y,y0,dydt,RK_tf,RK_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,rk_trace_ptr);
					/*fp[13] = 0.5;  				
	        flag = iz_rk(y,y0,dydt,RK_tf,RK_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,rk_trace_ptr);
	        flag = iz_rk(y,y0,dydt,RK_tf,RK_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,rk_trace_ptr);*/
        }
			} /*if(flag==1) break; /*Loop over cells*/
  	} /*if(flag==1) break; /*loop over steps*/
 	} /*loop over t_ms*/
 	t_cpu[0] = (double)(clock() - c0_rk)/(CLOCKS_PER_SEC);
 	mexPrintf("Time = %5.2f. \t",t_cpu[0]); 
 	mexPrintf("Integration steps = %d. \n",icount); fflush(stdout);
  
 	/************************************************************/
  /************ Adaptive Parker-Sochacki method ***************/
  /************************************************************/
  fp[8]=co_g_ampa_ps; fp[9]=co_g_gaba_ps; fp[12]=dt_ps; ip[1]=(int)steps_ps;
  /*Scale time/rate constants such that dt=1 in the equations*/  
 	for(nrnp = nrn; nrnp < nrnx; nrnp++){nrnp->E=E*dt_ps; nrnp->a = a*dt_ps;}
  mexPrintf("PS\n");
  for(p = 1; p < order_lim; p++){  	
  	co[0][p] = nrn[0].E/((double)(p+1));/*assumes all E are the same*/
  	co[1][p] = nrn[0].a/((double)(p+1));/*assumes all a are the same*/
  	co[2][p] = -1.0/(tau_ampa_ps*(double)(p+1));
  	co[3][p] = -1.0/(tau_gaba_ps*(double)(p+1)); 
	}
  for(nrnp = nrn,i=0; nrnp < nrnx; nrnp++){ 
  	nrnp->v=0; nrnp->u=0; nrnp->g_gaba=0; nrnp->g_ampa=0;
    if(sim_type>0) nrnp->I = rand_inj[i];
    else nrnp->I = inj;
    nrn[i].n_in = 0;ne[i] = -1; i++;
  }
  fcount=0; icount=0; flag=0; mu_order = 0.0; max_order = 0; c0_ps = clock(); 
	for(t_ms=0; t_ms<t_end; t_ms++){
		PS_v[t_ms] = nrn[trace_ind].v;
    if(sim_type>0 && t_ms==50)for(nrnp=nrn;nrnp<nrnx;nrnp++)nrnp->I=0;
		for(step=0; step<steps_ps; step++){
			t = (double)t_ms + (double)step*dt_ps;
			if(sim_type == 2)get_local_inputs(nrn,ne,te,syn,n_nrn);
			for(nrnp=nrn, nrn_ind=0; nrnp<nrnx; nrnp++, nrn_ind++){
				ip[2] = step; ip[3] = t_ms;	ip[4] = nrn_ind; 
        fp[15] = t; /*real time at start of step*/
        /*Work through substeps separated by synaptic events*/
        if(sim_type == 2){
	        n_in = nrnp->n_in;
          start = 0; fp[16] = start; /*start time of substep (in [0 1] of whole step)*/
	        if(n_in){
	        	for(substep=0;substep<n_in;substep++){ /*one substep per event*/
	        		if(nrn_ind == trace_ind && cond==4){ /*record inputs to trace neuron*/
	        			trace_in_t = (double)(t+dt_ps*nrnp->in_t[substep]);
                fwrite(&trace_in_t,8,1,trace_in_t_ptr);
	        			trace_in_w = (double)(nrnp->in_w[substep]);
	        			fwrite(&trace_in_w,8,1,trace_in_w_ptr);
	        			trace_n_in++;	
	        		}
	          	fp[13] = nrnp->in_t[substep] - start;
	          	if(fp[13]>0){
	          		flag = iz_ps(yp,co,yold,ynew,PS_tf,PS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_order,&max_order,ps_trace_ptr,order_lim);
              	start = nrnp->in_t[substep]; fp[16] = start; fp[15]+=dt_ps*fp[13];
	          	}
	          	w = nrnp->in_w[substep];
	          	if(w > 0) nrnp->g_ampa+=w; /*AMPA*/
	          	else nrnp->g_gaba-=w; /*GABA*/          	
	          }
            fp[13] = 1-start; /*remainder of time step*/
	          flag = iz_ps(yp,co,yold,ynew,PS_tf,PS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_order,&max_order,ps_trace_ptr,order_lim);
	        }
	        else{
						fp[13] = 1;  				
		        flag = iz_ps(yp,co,yold,ynew,PS_tf,PS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_order,&max_order,ps_trace_ptr,order_lim);
	        }
        }
        else{
          fp[13] = 1;
        	flag = iz_ps(yp,co,yold,ynew,PS_tf,PS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_order,&max_order,ps_trace_ptr,order_lim);
        	/*fp[13] = 0.5;
        	flag = iz_ps(yp,co,yold,ynew,PS_tf,PS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_order,&max_order,ps_trace_ptr,order_lim);
        	flag = iz_ps(yp,co,yold,ynew,PS_tf,PS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_order,&max_order,ps_trace_ptr,order_lim);*/
        }
      } /*if(flag==1) break; /*Loop over cells*/
  	} /*if(flag==1) break; /*loop over steps*/
  } /*loop over t_ms*/
	t_cpu[2] = (double)(clock() - c0_ps)/(CLOCKS_PER_SEC);
	mexPrintf("Time = %5.2f. \t",t_cpu[2]); 
	mexPrintf("Integration steps = %d. \t",icount); 
	mexPrintf("Mean order = %5.10f. \t",mu_order);
	mexPrintf("Max order = %d. \n",max_order);
	fflush(stdout);

	/************************************************************/
  /********************* Bulirsch-Stoer ***********************/
  /************************************************************/
  mexPrintf("BS\n");
  n_bs_fails=0;
  for(nrnp = nrn,i=0; nrnp < nrnx; nrnp++){ 
  	nrnp->v=0; nrnp->u=0; nrnp->g_gaba=0; nrnp->g_ampa=0;
    if(sim_type>0) nrnp->I = rand_inj[i];
    else nrnp->I = inj;
    nrn[i].n_in = 0;ne[i] = -1; i++;
  }
  fcount=0; icount=0; flag=0; mu_bsk = 0.0; max_bsk = 0; c0_bs = clock(); 
  for(t_ms=0; t_ms<t_end; t_ms++){
  	if(ps_only==1) break;
  	BS_v[t_ms] = nrn[trace_ind].v;
    if(sim_type>0 && t_ms==50)for(nrnp=nrn;nrnp<nrnx;nrnp++)nrnp->I=0;
  	for(step=0; step<steps_ps; step++){
  		t = (double)t_ms + (double)step*dt_ps;
      if(sim_type == 2)get_local_inputs(nrn,ne,te,syn,n_nrn);
			for(nrnp=nrn, nrn_ind=0; nrnp<nrnx; nrnp++, nrn_ind++){
				ip[2] = step; ip[3] = t_ms;	ip[4] = nrn_ind; 
        fp[15] = t; /*real time at start of step*/
        /*Work through substeps separated by synaptic events*/
        if(sim_type == 2){
	        n_in = nrnp->n_in;
          start = 0; fp[16] = start; /*start time of substep (in [0 1] of whole step)*/
	        if(n_in){
	        	for(substep=0;substep<n_in;substep++){ /*one substep per event*/
	          	fp[13] = nrnp->in_t[substep] - start; 
	          	if(fp[13]>0){
	          		flag = iz_bs(y,y0,dydt,BS_tf,BS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_bsk,&max_bsk,bs_trace_ptr); n_bs_fails+=flag;
              	start = nrnp->in_t[substep]; fp[16] = start; fp[15]+=dt_ps*fp[13];
	          	}
	          	w = nrnp->in_w[substep];
	          	if(w > 0) nrnp->g_ampa+=w; /*AMPA*/
	          	else nrnp->g_gaba-=w; /*GABA*/
	          }
	          fp[13] = 1-start; /*remainder of time step*/
	          flag = iz_bs(y,y0,dydt,BS_tf,BS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_bsk,&max_bsk,bs_trace_ptr); n_bs_fails+=flag;
	        }
	        else{
						fp[13] = 1;  				
		        flag = iz_bs(y,y0,dydt,BS_tf,BS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_bsk,&max_bsk,bs_trace_ptr); n_bs_fails+=flag;
	        }
        }
        else{
          fp[13] = 1;  				
	        flag = iz_bs(y,y0,dydt,BS_tf,BS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_bsk,&max_bsk,bs_trace_ptr); n_bs_fails+=flag;
        	/*fp[13] = 0.5;  				
	        flag = iz_bs(y,y0,dydt,BS_tf,BS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_bsk,&max_bsk,bs_trace_ptr); n_bs_fails+=flag;
	        flag = iz_bs(y,y0,dydt,BS_tf,BS_nrn,nrnp,ne,te,ip,fp,&fcount,&icount,&mu_bsk,&max_bsk,bs_trace_ptr); n_bs_fails+=flag;*/
        }
      } /*if(flag==1) break; /*Loop over cells*/
  	} /*if(flag==1) break; /*loop over steps*/
 	} /*loop over t_ms*/
 	t_cpu[1] = (double)(clock() - c0_bs)/(CLOCKS_PER_SEC);
 	mexPrintf("Time = %5.2f.\t",t_cpu[1]); mexPrintf("Integration steps = %d. \t",icount);
 	mexPrintf("Mean sub = %5.10f. \t",mu_bsk*2);
 	mexPrintf("Max sub = %d. \t",max_bsk*2);
 	mexPrintf("n BS fails = %d.\n",n_bs_fails); fflush(stdout);
 	
  free(nrn); free(yold); free(ynew); free(y); free(y0); free(dydt);
  free(ne);free(te);free(rand_inj);
 	for(i=0;i<n_nrn;i++)free(syn[i]); free(syn);
 	for(i=0;i<=nv;i++){free(yp[i]); free(co[i]);} free(yp); free(co);
  if(trace_rec==1){
    fclose(ps_trace_ptr); fclose(bs_trace_ptr); fclose(rk_trace_ptr);
    if(cond == 4){fclose(trace_in_t_ptr); fclose(trace_in_w_ptr);}
  }
 	
 	i_stats[0] = max_order; i_stats[1] = max_bsk*2; i_stats[2] = n_bs_fails;
 	f_stats[0] = (double)mu_order; f_stats[1] = (double)mu_bsk*2;
}

/*The gateway routine*/
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
  double *fp, *RK_tf, *PS_tf, *BS_tf, *t_cpu, *RK_v, *PS_v, *BS_v, *f_stats;
  int *ip, *RK_nrn, *PS_nrn, *BS_nrn, *i_stats;
  int n_nrn, t_end;
  
  /*Check for proper number of arguments*/
  if (nrhs != 2) mexErrMsgTxt("Two inputs required.");
  if (nlhs != 12) mexErrMsgTxt("Twelve outputs required.");
	
  /*Get the inputs*/
  fp = (double *)mxGetData(FP);	/*Cell type as integer index*/
  ip = (int *)mxGetData(IP);
  n_nrn = ip[0]; t_end = ip[2]; 	
	
  /*Create output matrices*/
  RK_TF = mxCreateDoubleMatrix(n_nrn*100,1,mxREAL);
  PS_TF = mxCreateDoubleMatrix(n_nrn*100,1,mxREAL);
  BS_TF = mxCreateDoubleMatrix(n_nrn*100,1,mxREAL);
  RK_NRN = mxCreateNumericMatrix(n_nrn*100,1,mxINT32_CLASS,mxREAL);
  PS_NRN = mxCreateNumericMatrix(n_nrn*100,1,mxINT32_CLASS,mxREAL);
  BS_NRN = mxCreateNumericMatrix(n_nrn*100,1,mxINT32_CLASS,mxREAL);
  T_CPU = mxCreateDoubleMatrix(3,1,mxREAL);
  RK_V = mxCreateDoubleMatrix(t_end,1,mxREAL);
  PS_V = mxCreateDoubleMatrix(t_end,1,mxREAL);
  BS_V = mxCreateDoubleMatrix(t_end,1,mxREAL);
  I_STATS = mxCreateNumericMatrix(3,1,mxINT32_CLASS,mxREAL);
  F_STATS = mxCreateDoubleMatrix(2,1,mxREAL);  
  
  /*Set pointers to outputs*/
  RK_tf = mxGetData(RK_TF);
  RK_nrn = mxGetData(RK_NRN);
  PS_tf = mxGetData(PS_TF);
  PS_nrn = mxGetData(PS_NRN);
  BS_tf = mxGetData(BS_TF);
  BS_nrn = mxGetData(BS_NRN);
  t_cpu = mxGetData(T_CPU);
  RK_v = mxGetData(RK_V);
  PS_v = mxGetData(PS_V);
  BS_v = mxGetData(BS_V);
  i_stats = mxGetData(I_STATS);
  f_stats = mxGetData(F_STATS);  

  /*Call the C subroutine*/
  run_sim(RK_tf,RK_nrn,PS_tf,PS_nrn,BS_tf,BS_nrn,t_cpu,RK_v,PS_v,BS_v,i_stats,f_stats,fp,ip);
}
