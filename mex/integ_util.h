#ifndef INC_INTEG_UTIL_H
#define INC_INTEG_UTIL_H

extern void ps_update(double **, int, int, double, double *);

extern int ps_step(double **,double **,double *,double *,double *,
  double *,void (*)(double **,double **,double *), 
	void (*)(double **,double **,double *,int), int, int, int, int);

extern int ps_step0(double **,double **,double *,double *,double *,
  void (*)(double **,double **,double *), 
	void (*)(double **,double **,double *,int), int, int, int, int);

extern void cauchy_prod(int, double *, double, double *, double, double *);

extern void cauchy_prod_basic(int, double *, double *, double *);

extern void cauchy_prod_signed(int, double *, double, int, double *, double, int, double *);

extern double cauchy_prod_gain(int, double *, double, double, double *, double, double, double *);

extern void series_div(int, double, double *, double, double *, double);

extern void series_pow(int, double *, double, double *, double, double);

extern void rk_step(double *, double *, double *, int, double *, double,
	void (*)(double *, double *, double *));
  
extern void mmid(double *, double *, int, double, int, double *, double *,
	void (*)(double *, double *, double *));

extern void pzextr(int,double,double *,double *,double *,int,double *,double **);

extern void rzextr(int,double,double *,double *,double *,int,double *,double **);

extern int bs_step(double *, double *, double *, int, double, double *, double *, 
	void (*)(double *, double *, double *));

#endif /* INC_INTEG_UTIL_H */
