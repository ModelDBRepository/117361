/*Utilities for numerical integration of the modified Traub-Miles neuron*/
/*Written by Dr Robert Stewart for Stewart & Bair, 2009*/

/*Derivs routine*/
void tm_derivs(double *y, double *dydx, double *fp){
	double alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h;
	/*Parameters in nS for cell of 20,000 micometer surface area (2e-4 cm^2)*/
	double gNa = 20000, gK = 6000, gL = 10;
	double C = 200; /*pF*/
 	double ENa = 50; /*Sodium reversal potential in mV*/
	double EK = -90; /*Potassium reversal potential in mV*/
	double EL = -65; /*Leak reversal potential in mV*/
	double E_ampa = 0; /*AMPA reversal potential in mV*/
	double E_gaba = -80; /*GABA reversal potential in mV*/
	double Vr = -63;
	double E_alpha_n = Vr+15;
	double E_beta_n = Vr+10;
	double E_alpha_m = Vr+13;
	double E_beta_m = Vr+40;
	double E_alpha_h = Vr+17;
	double E_beta_h = Vr+40;

	double v = y[0];
	double n = y[1];
	double m = y[2];
	double h = y[3];
	double g_ampa = y[4];
	double g_gaba = y[5];
	double n4 = n*n*n*n; 
	double m3h = m*m*m*h;
	double dt = fp[0];
	double co_g_ampa = fp[1]; 
	double co_g_gaba = fp[2];
	double I_inj = fp[8];

	dydx[0] = dt*(I_inj - gK*n4*(v-EK) - gNa*m3h*(v-ENa) - gL*(v-EL) - g_ampa*(v-E_ampa) - g_gaba*(v-E_gaba))/C;
	if(v == E_alpha_n) alpha_n = 0.032*5; /*protect against div by zero*/
	else alpha_n = 0.032 * (E_alpha_n-v)/(exp((E_alpha_n-v)/5) - 1);
	beta_n = 0.5*exp((E_beta_n-v)/40);      
	dydx[1] = dt*(alpha_n*(1-n) - beta_n*n);
	if(v == E_alpha_m) alpha_m = 0.32*4; /*protect against div by zero*/
	else alpha_m = 0.32 * (E_alpha_m-v)/(exp((E_alpha_m-v)/4) - 1);
	if(v == E_beta_m) beta_m = 0.28*5; /*protect against div by zero*/
	else beta_m = 0.28 * (v-E_beta_m)/(exp((v-E_beta_m)/5) - 1);
	dydx[2] = dt*(alpha_m*(1-m) - beta_m*m);
	alpha_h = 0.128 * exp((E_alpha_h-v)/18); 
	beta_h = 4 / (exp((E_beta_h-v)/5)+1);
	dydx[3] = dt*(alpha_h*(1-h) - beta_h*h);
	dydx[4] = co_g_ampa*g_ampa;
	dydx[5] = co_g_gaba*g_gaba;
}

void tm_first(double **y, double **co, double *fp){
	double v,n,m,h,a,b,c,d,g_ampa,g_gaba,alpha_n,alpha_m,beta_m,n4,m3h,v1;
	double e,f,g,q,r,s,v_alpha_n,v_alpha_m,v_beta_m,co_alpha_n,co_alpha_m,co_beta_m;
	double chi,psi,xi,co_K,co_Na,I_inj;
	double gNa = 20000, gK = 6000, gL = 10;
 	double ENa = 50; /*Sodium reversal potential in mV*/
	double EK = -90; /*Potassium reversal potential in mV*/
	double EL = -65; /*Leak reversal potential in mV*/
	double E_ampa = 0; /*AMPA reversal potential in mV*/
	double E_gaba = -80; /*GABA reversal potential in mV*/
	
	alpha_n = fp[3]; alpha_m = fp[4]; co_K = fp[5]; co_Na = fp[6];
	I_inj = fp[8];
	v=y[0][0]; n=y[1][0]; m=y[2][0]; h=y[3][0]; g_ampa=y[4][0]; g_gaba=y[5][0];
	a=y[6][0]; b=y[7][0]; c=y[8][0]; d=y[9][0]; 
	e=y[10][0]; f=y[11][0];	g=y[12][0]; q=y[13][0]; r=y[14][0]; s = y[15][0];
	chi = y[16][0]; psi = y[17][0]; xi = y[18][0];
	
	v1 = co[0][0]*(v*chi+co_K*EK+co_Na*ENa+gL*EL+g_ampa*E_ampa+g_gaba*E_gaba+I_inj);
	y[0][1] = v1;
	y[1][1] = co[1][0]*(n*psi + alpha_n);/*n*/
	y[2][1] = co[2][0]*(m*xi + alpha_m);/*m*/
	y[3][1] = co[3][0]*(b*(1-h) - 4*d);/*h*/
	y[4][1] = co[4][0]*g_ampa;
	y[5][1] = co[5][0]*g_gaba;
	y[6][1] = co[6][0]*v1*a;/*a*/
	y[7][1] = co[7][0]*v1*b;/*b*/
	y[8][1] = co[8][0]*v1*c;/*c*/
	y[9][1] = (y[3][1] - d*y[8][1])/(c+1); /*d = h/(c+1)*/
	y[10][1] = co[10][0]*v1*e;/*e*/
	y[11][1] = co[11][0]*v1*f;/*f*/
	y[12][1] = co[12][0]*v1*g;/*g*/
	y[13][1] = (-v1 - q*y[10][1])/(e-1);/*q*/
	y[14][1] = (-v1 - r*y[11][1])/(f-1);/*r*/
	y[15][1] = (v1 - s*y[12][1])/(g-1);/*s*/
}

void tm_iter(double **y, double **co, double *fp, int p){
	double v,n,m,h,g_ampa,g_gaba,a,b,c,d,e,f,g,q,r,s,chi,psi,xi;
	double v_chi,alpha_n,alpha_m,beta_m,co_K,co_Na; 
	double n_psi,m_xi,bh,n4,m3h,cd,eq,fr,gs,temp;
	int i,nv=6;
	double gNa = 20000, gK = 6000, gL = 10;
 	double ENa = 50; /*Sodium reversal potential in mV*/
	double EK = -90; /*Potassium reversal potential in mV*/
	double EL = -65; /*Leak reversal potential in mV*/
	double E_ampa = 0; /*AMPA reversal potential in mV*/
	double E_gaba = -80; /*GABA reversal potential in mV*/
	double co_alpha_n = 0.032;
	double co_alpha_m = 0.32;
	double co_beta_m = 0.28;
	
	v=y[0][0]; n=y[1][0]; m=y[2][0]; h=y[3][0]; g_ampa=y[4][0]; g_gaba=y[5][0];
	a=y[6][0]; b=y[7][0]; c=y[8][0]; d=y[9][0]; 
	e=y[10][0]; f=y[11][0];	g=y[12][0]; q=y[13][0]; r=y[14][0]; s = y[15][0];
	chi = y[16][0]; psi = y[17][0]; xi = y[18][0];
	
	y[21][0] = m*m*m; y[23][0] = n*n*n*n; y[19][0] = n*n; y[24][0] = n*n;
  
  /*psi and n times psi*/
  alpha_n = co_alpha_n*y[13][p];    		
  y[17][p] = -(alpha_n+y[6][p]); /*psi*/
  cauchy_prod(p,y[17],psi,y[1],n,&n_psi); /* n*psi */

  /*xi and m times xi*/
  alpha_m = co_alpha_m*y[14][p];
  beta_m = co_beta_m*y[15][p];
  y[18][p] = -(alpha_m+beta_m); /*xi*/
  cauchy_prod(p,y[18],xi,y[2],m,&m_xi); /* m*xi */

  /*powers*/
  cauchy_prod(p, y[1], n, y[1], n, &y[19][p]); /*n^2*/
  cauchy_prod(p, y[2], m, y[2], m, &y[20][p]); /*m^2*/ 
  cauchy_prod(p, y[20], m*m, y[2], m, &y[21][p]); /*m^3*/ 

  /*testing direct power algorithm from Knuth (p526)*/ 
  /*series_pow(p,y[1],n,y[23],n*n*n*n,4.0); /*n^4*/
  /*series_pow(p,y[23],n*n*n*n,y[24],n*n,0.5); /*(n^4)^.5*/
  /*mexPrintf("n2 error: %3.15f\n",y[24][p]-y[19][p]); /*square root works :-)*/

  cauchy_prod(p,y[7],b,y[3],h-1,&bh);    
  cauchy_prod(p, y[19], n*n, y[19], n*n, &n4); 
  cauchy_prod(p, y[21], m*m*m, y[3], h, &m3h);
  co_K = gK*n4; co_Na = gNa*m3h; 
  y[16][p] = -co_K - co_Na - y[4][p] - y[5][p]; /*chi*/
  cauchy_prod(p,y[16],chi,y[0],v, &v_chi); 

  /*v*/
  y[22][p] = co[0][0]*(v_chi + co_K*EK + co_Na*ENa + y[4][p]*E_ampa + y[5][p]*E_gaba);
  y[0][p+1] = y[22][p]/(p+1); /*v*/

  /*n,m,h,g_ampa,g_gaba*/
  y[1][p+1] = co[1][p]*(n_psi + alpha_n); /*n*/
  y[2][p+1] = co[2][p]*(m_xi + alpha_m); /*m*/
  y[3][p+1] = co[3][p]*(-bh - 4*y[9][p]); /*h*/
  y[4][p+1] = co[4][p]*y[4][p]; /*g_ampa*/
  y[5][p+1] = co[5][p]*y[5][p]; /*g_gaba*/

  /*a,b,c,d*/
  cauchy_prod(p,y[22],y[0][1],y[6],a,&y[6][p+1]); y[6][p+1]*=co[6][p];/*a*/ 
  cauchy_prod(p,y[22],y[0][1],y[7],b,&y[7][p+1]); y[7][p+1]*=co[7][p];/*b*/
  cauchy_prod(p,y[22],y[0][1],y[8],c,&y[8][p+1]); y[8][p+1]*=co[8][p];/*c*/
  series_div(p,y[3][p+1],y[8],(c+1),y[9],d); /*d*/

  cauchy_prod(p,y[22],y[0][1],y[10],e,&y[10][p+1]);y[10][p+1]*=co[10][p];/*e*/ 
  series_div(p,-y[0][p+1],y[10],(e-1),y[13],q); /*q*/

  cauchy_prod(p,y[22],y[0][1],y[11],f,&y[11][p+1]);y[11][p+1]*=co[11][p];/*f*/
  series_div(p,-y[0][p+1],y[11],(f-1),y[14],r); /*r*/

  cauchy_prod(p,y[22],y[0][1],y[12],g,&y[12][p+1]);y[12][p+1]*=co[12][p];/*g*/ 
  series_div(p,y[0][p+1],y[12],(g-1),y[15],s); /*g*/
}
