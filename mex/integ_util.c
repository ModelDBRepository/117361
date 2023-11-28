/*Integration routines for the Parker-Sochacki, Runge-Kutta & Bulirsch-Stoer methods*/

/*Parker-Sochacki stepper - Solve and update in one function*/
int ps_step(double **y,double **co,double *y1,double *ynew,double *fp,
  double eta[],void (*first)(double **,double **,double *), 
	void (*iter)(double **,double **,double *,int),int stop,int ps_limit, int nv, int err_nv){
  int i,p; double dt=fp[99], dt_pow;
	first(y,co,fp); /*Calculate first order terms*/
	if(dt == 1)for(i=0;i<nv;i++)y1[i]=y[i][0]+y[i][1];
  else{for(i=0;i<nv;i++)y1[i]=y[i][0]+dt*y[i][1]; dt_pow=dt*dt;}
  for(p=1;p<(ps_limit-1);p++){/*Iterations*/
    iter(y,co,fp,p);
		if(dt == 1)for(i=0;i<nv;i++)ynew[i]=y1[i]+y[i][p+1]; /*Update*/
    else{for(i=0;i<nv;i++)ynew[i]=y1[i]+y[i][p+1]*dt_pow; dt_pow*=dt;}
    if((fabs(y[0][p+1])>10.0)){p=-1;break;} /*Check for divergence*/
    /*Test error tolerance on variable value change*/
    for(i=0;i<err_nv;i++){if(fabs(ynew[i]-y1[i])>eta[i]) break;}
    if(i==err_nv)break;
    for(i=0;i<nv;i++)y1[i]=ynew[i];
  } p++;
	if(stop==1){if(p==ps_limit)mexErrMsgTxt("PS solve failed.");}
	return p;
}

/*Parker-Sochacki stepper - zero tolerance version*/
int ps_step0(double **y,double **co,double *y1,double *ynew,double *fp,
  void (*first)(double **,double **,double *), 
	void (*iter)(double **,double **,double *,int),int stop,int ps_limit, int nv, int err_nv){
  int i,p; double dt=fp[99], dt_pow;
  first(y,co,fp); /*Calculate first order terms*/
	if(dt == 1)for(i=0;i<nv;i++)y1[i]=y[i][0]+y[i][1];
  else{for(i=0;i<nv;i++)y1[i]=y[i][0]+dt*y[i][1]; dt_pow=dt*dt;}
  for(p=1;p<(ps_limit-1);p++){/*Iterations*/
    iter(y,co,fp,p);
		if(dt == 1)for(i=0;i<nv;i++)ynew[i]=y1[i]+y[i][p+1]; /*Update*/
    else{for(i=0;i<nv;i++)ynew[i]=y1[i]+y[i][p+1]*dt_pow; dt_pow*=dt;}
    if((fabs(y[0][p+1])>10.0)){p=-1;break;} /*Check for divergence*/
    for(i=0;i<err_nv;i++){if(ynew[i]-y1[i]) break;} /*zero tolerance*/
    if(i==err_nv)break;
    for(i=0;i<nv;i++)y1[i]=ynew[i];
  } p++;
	if(stop==1){if(p==ps_limit)mexErrMsgTxt("PS solve failed.");}
	return p;
}

void ps_update(double **y,int var,int ps_order,double dt,double *ynew){
	int i;
	if(dt == 1){for(i=1,*ynew=y[var][0];i<=ps_order;i++)*ynew+=y[var][i];} 
	else{ /*Use Horner's method*/
		*ynew=y[var][ps_order]*dt+y[var][ps_order-1];
		for(i=ps_order-2;i>=0;i--){*ynew=y[var][i]+*ynew*dt;}
	}
}

void cauchy_prod(int p,double *a,double a0,double *b,double b0,double *c){
	/*c is the pth term of the Cauchy product of a and b with zeroth order 
	terms a0, b0 allowing shifted products (e.g (a-1)*(b+1))*/
 	int i;
 	*c = a0*b[p] + b0*a[p];
	for(i = 1; i < p; i++){*c += a[i]*b[p-i];}
}

void cauchy_prod_basic(int p,double *a,double *b,double *c){
	/*c is the pth term of the Cauchy product of a and b*/
 	int i;
 	*c = a[0]*b[p] + b[0]*a[p];
	for(i = 1; i < p; i++){*c += a[i]*b[p-i];}
}

void cauchy_prod_signed(int p,double *a,double a0,int sign_a, 
	double *b,double b0,int sign_b,double *c){
	/*Cauchy product with sign parameters to allow negative series*/
 	int i;
 	if(sign_a == 1){ /*assume sign_b is negative*/
 		*c = -a0*b[p] + b0*a[p];
 		for(i = 1; i < p; i++){*c += -a[i]*b[p-i];}
 	}
 	else if(sign_b == 1){ /*assume sign_a is negative*/
 		*c = a0*b[p] - b0*a[p];
 		for(i = 1; i < p; i++){*c += -a[i]*b[p-i];}
 	}
 	else{ /*assume both signs are negative*/
 		*c = -a0*b[p] - b0*a[p];
 		for(i = 1; i < p; i++){*c += a[i]*b[p-i];}
 	}
}

double cauchy_prod_gain(int p,double *a,double a0,double gain_a, 
	double *b,double b0,double gain_b,double *c){
	/*Cauchy product with gain parameters to allow multiples of series*/
 	int i;
 	*c = a0*gain_b*b[p] + b0*gain_a*a[p];
	for(i = 1; i < p; i++){*c += gain_a*a[i]*gain_b*b[p-i];}
}

void series_div(int p,double a_pp,double *b,double b0,double *c,double c0){
	/*calculates pth term of c = a/b, with zeroth order terms b0, c0*/
	double cb;
	cauchy_prod(p,c,c0,b+1,b[1],&cb);
	c[p+1] = (a_pp - cb)/b0;
}

void series_pow(int p,double *a,double a0,double *b,double b0,double x){
	/*calculates pth coefficient of b = a^x*/
	int i;
	b[p] = x*a[p]*b0; /*i=p special case*/
	for(i=1;i<p;i++){b[p] += ((x+1)*(double)i/(double)p - 1)*a[i]*b[p-i];} 
	b[p] /= a0;
}

/*Non-adaptive 4th order Runge-Kutta stepper*/
void rk_step(double *y,double *y0,double *dydt0,int nv, 
	double *fp,double dt,void (*derivs)(double *,double *,double *)){
    double dt2 = dt/2; double dt6 = dt/6;
		int i; 
    double *rk1,*rk2,*rk3,*dydt;
    rk1 = malloc(nv*sizeof(double)); rk2 = malloc(nv*sizeof(double));
    rk3 = malloc(nv*sizeof(double)); dydt = malloc(nv*sizeof(double));
		for(i=0;i<nv;i++){rk1[i] = dydt0[i]; y[i] = y0[i] + dt2*dydt0[i];} /*1*/
		derivs(y, dydt, fp); /*2*/
		for(i=0;i<nv;i++){rk2[i] = dydt[i]; y[i] = y0[i] + dt2*dydt[i];}		
		derivs(y, dydt, fp); /*3*/
		for(i=0;i<nv;i++){rk3[i] = dydt[i]; y[i] = y0[i] + dt*dydt[i];}	
		derivs(y, dydt, fp); /*4*/
		for(i=0;i<nv;i++)y[i] = y0[i] + dt6*(rk1[i]+dydt[i]+2*(rk2[i]+rk3[i]));
    free(rk1);free(rk2);free(rk3);free(dydt);
}

/*Numerical Recipes modified midpoint method routine - For Bulirsch-Stoer*/
void mmid(double *y,double *dydt,int nv,double htot,int nstep,double *yout,double *fp,void (*derivs)(double *,double *,double *)){
	int n,i;
	double swap,h2,h,*ym,*yn;
	ym = malloc(nv*sizeof(double));	yn = malloc(nv*sizeof(double));
	h = htot/(double)nstep; /*stepsize this trip*/
	for(i=0;i<nv;i++){ym[i]=y[i];yn[i]=y[i]+h*dydt[i];} /*First step*/
	derivs(yn,yout,fp); /*Uses yout for temporary storage of derivatives*/
	h2=2.0*h;
	for(n=2;n<=nstep;n++){ /*General step*/
		for(i=0;i<nv;i++){swap = ym[i] + h2*yout[i];ym[i] = yn[i];yn[i] = swap;}
		derivs(yn,yout,fp);
	}
	for(i=0;i<nv;i++){yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);}/*Last step*/
	free(ym);	free(yn);
}

/*NR polynomial extrapolation function - For Bulirsch-Stoer*/
void pzextr(int iest,double xest,double *yest,double *yz,double *yerr, 
  int nv,double *x,double **d){
/*Evaluate nv functions at x=0 by fitting a	polynomial to a sequence of estimates with 
progressively smaller values x = xest, and corresponding function vectors yest[1..nv]. 
This call is number iest in the sequence of calls. Extrapolated function values are 
output as yz, and their estimated error is output as yerr.*/
	int k1,j;	double q,f2,f1,delta,*c;
  c = malloc(nv*sizeof(double));
  x[iest]=xest; /*save current independent variable*/
  for(j=0;j<nv;j++){yerr[j]=yz[j]=yest[j];}
  if (iest==0){for(j=0;j<nv;j++){d[j][0]=yest[j];}} /*store first element*/
  else{
  	for(j=0;j<nv;j++){c[j]=yest[j];}
  	for (k1=0;k1<iest;k1++){
  		delta=1.0/(x[iest-k1-1]-xest);
  		f1=xest*delta;
  		f2=x[iest-k1-1]*delta;
  		for(j=0;j<nv;j++){ /*propagate tableau 1 diagonal more*/
  			q=d[j][k1]; d[j][k1]=yerr[j];	delta=c[j]-q;
  			yerr[j]=f1*delta;	c[j]=f2*delta; yz[j]+=yerr[j];
  		}
  	}
  	for(j=0;j<nv;j++){d[j][iest]=yerr[j];}
  } free(c);
}

/*NR rational extrapolation function - For Bulirsch-Stoer*/
void rzextr(int iest,double xest,double *yest,double *yz,double *yerr,
  int nv,double *x,double **d){
	int k,j; double yy,v,ddy,c,b1,b,*fx;
  fx = malloc((iest+1)*sizeof(double));
  x[iest]=xest; /*save current independent variable*/
  if (iest == 0){for(j=0;j<nv;j++){yz[j]=yest[j];d[j][0]=yest[j];yerr[j]=yest[j];}}
  else{
  	for(k=0;k<iest;k++){fx[k+1] = x[iest-k-1]/xest;}
  	for(j=0;j<nv;j++){
  		v=d[j][0]; d[j][0]=yy=c=yest[j];
  		for(k=1;k<=iest;k++){
  			b1=fx[k]*v; b=b1-c;
  			if(b){b=(c-v)/b;ddy=c*b;c=b1*b;}
  			else ddy=v;
  			if(k != iest) v=d[j][k];
  			d[j][k]=ddy; yy+=ddy;
  		}
  		yerr[j]=ddy; yz[j]=yy;
  	}
  } free(fx);
}

/*Fixed-step version of Bulirsch-Stoer*/
int bs_step(double *y,double *y0,double *dydt,int nv,double dt,
  double *fp,double *eta,void (*derivs)(double *,double *,double *)){
	/*y holds the state variables; dydt is the derivative at the start of the step*/
	int i,k,n; static int kmax = 50, nseq[100];
	double xest, x[100], *yerr, *yseq, *y1, **d;
  d = malloc(nv*sizeof(double *)); yerr = malloc(nv*sizeof(double));
  yseq = malloc(nv*sizeof(double)); y1 = malloc(nv*sizeof(double));
  if(nv>50)mexErrMsgTxt("BS: Maximum number of variables (50) exceeded");
  for(i=0;i<100;i++)nseq[i] = 2*(i+1);
  for(i=0;i<nv;i++){d[i] = malloc(100*sizeof(double));}
  for(k=0;k<kmax;k++){ /*evaluate the sequence of modified midpoint steps*/
    for(i=0;i<nv;i++)y1[i]=y[i]; /*store previous iteration solution*/
  	n = nseq[k];
 		mmid(y0,dydt,nv,dt,n,yseq,fp,derivs); /*call mmid*/
 		xest=(dt/n)*(dt/n); /*squared, since error series is even*/
 		/*pzextr(k,xest,yseq,y,yerr,nv,x,d); /*polynomial extrapolation*/
 		rzextr(k,xest,yseq,y,yerr,nv,x,d); /*rational extrapolation*/
    /*Solution change-based error tolerance*/
    for(i=0;i<nv;i++){if(fabs(y[i]-y1[i])>eta[i]) break;} if(i==nv)break;
 	}
  for(i=0;i<nv;i++){free(d[i]);} free(d);free(yerr);free(yseq);free(y1);
  return k;
}
