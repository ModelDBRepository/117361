/*Utilities for numerical integration of the Izhikevich model with COBA synapses*/
/*Written by Dr Robert Stewart for Stewart & Bair, 2009*/

/*Calculates the derivatives*/
void iz_derivs(double *y, double *dydx, double *fp){
	double v,u,g_ampa,g_gaba,I,k,l,E_ampa,E_gaba,E,a,b,co_g_ampa,co_g_gaba;
	v = y[0]; u = y[1]; g_ampa = y[2]; g_gaba = y[3]; /*extract variables*/
	I = fp[0]; k = fp[1]; l = fp[2]; E_ampa = fp[3]; E_gaba = fp[4];
	E = fp[5]; a = fp[6]; b = fp[7]; co_g_ampa = fp[8]; co_g_gaba = fp[9];
	dydx[0] = E*(k*v*v + l*v - u + I - g_ampa*(v-E_ampa) - g_gaba*(v-E_gaba));
	dydx[1] = a*(b*v - u);
	dydx[2] = co_g_ampa*g_ampa;
	dydx[3] = co_g_gaba*g_gaba;
}

/*PS - first term*/
void iz_first(double **y, double **co, double *fp){ 
	double v,u,g_ampa,g_gaba,I,k,l,E_ampa,E_gaba,E,a,b,co_g_ampa,co_g_gaba,chi;
	v = y[0][0]; u = y[1][0]; g_ampa = y[2][0]; g_gaba = y[3][0]; chi = y[4][0];
	I = fp[0]; k = fp[1]; E_ampa = fp[3]; E_gaba = fp[4];
	E = fp[5]; a = fp[6]; b = fp[7]; co_g_ampa = fp[8]; co_g_gaba = fp[9];
	y[0][1] = E*(v*chi - u + E_ampa*g_ampa + E_gaba*g_gaba + I);
	y[1][1] = a*(b*v - u);	
	y[2][1] = co_g_ampa*g_ampa;
	y[3][1] = co_g_gaba*g_gaba; 
	y[4][1] = k*y[0][1] - y[2][1] - y[3][1];
}

/*PS - iteration function for higher order terms*/
void iz_iter(double **y, double **co, double *fp, int p){ 
	double v,k,l,E_ampa,E_gaba,b,chi,vchi; int i;
	k = fp[1]; E_ampa = fp[3]; E_gaba = fp[4]; b = fp[7];
	vchi = y[0][0]*y[4][p] + y[4][0]*y[0][p];
	for(i = 1; i < p; i++){vchi += y[0][i]*y[4][p-i];}
	y[0][p+1] = co[0][p]*(vchi - y[1][p] + E_ampa*y[2][p] + E_gaba*y[3][p]);
	y[1][p+1] = co[1][p]*(b*y[0][p] - y[1][p]);
	y[2][p+1] = co[2][p]*y[2][p];
	y[3][p+1] = co[3][p]*y[3][p];
	y[4][p+1] = k*y[0][p+1] - y[2][p+1] - y[3][p+1];
}
