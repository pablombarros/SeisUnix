/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.		       */

/* SIRDEPIDMEIC: $Revision: 1.1 $ ; $Date: 2020/04/23 18:58:03 $        */

#include "par.h"
#include "rke.h"

/*********************** self documentation **********************/

char *sdoc[] = {
"								",
" SIRDEPIDMEIC - the SIRD EPIDEMIC model with			",
"								",
"  sirdepidemic > [stdout]					",
"								",
" Required Parameters: none					",
" Optional Parameters:						",
" S0=1000		initial number of susceptibles		",
" I0=1			initial number of infectives		",
" R0=0.0		initial number of removed (should be 0)	",
"	 		(not the basic reproducion rate r0	",
" D0=0			initial number of dead (should be 0)	",
" 								",
" k=.5			transmission rate			",
" b=.3333		removal rate = death + recovery rates	",
" nu=0.01		mortality rate 				",
" 			(0.01 = 1% mortality)			",
" 								",
" h=1			increment in time			",
" tol=1.e-08		error tolerance				",
" stepmax=40		maximum number of steps to compute	",
" mode=SIRD		S, I, R, D successively output		",
"			=S only, =I only, =R only		",
" Notes:							",
" This program is really just a demo showing how to use the 	",
" differential equation solver rke_solve written by Francois 	",
" Pinard, based on a modified form of the 4th order Runge-Kutta ",
" method, which employs the error checking method of R. England ",
" 1969.								",
"								",
" The output consists of unformated C-style binary floats, of	",
" either pairs or triplets as specified by the \"mode\" paramerter.",
"								",
" About the SIRD epidemic model:				",
" r0 = number of new infections per single infected host  	",
"  1 < r0 < 1.5 for influenza (2.2 to 2.7 for Covid-19)		",
"  b, k, S0, and r0 are related via				",
"  k = b*S0/r0							",
"  It is often easier to determine the recovery rate k (in units",
"  of h and to determine reasonable estimate of S0 and of r0 	",
"  and to calculate the infection rate b = k*r0/S0		",
"  nu= mortality rate						",
"								",
" S = total number susceptible to the infection			",
" I = total number of those capable of passing on the infection	",
" R = total number recovered					",
" D = total number dead						",
"								",
" Examples:							",
" Default:							",
" sirdepidemic | xgraph n=40 nplot=4 d1=1 &			",
" 								",
" Hypothetical flu in a school models:				",
"								",
" Without deaths, or with recovered and dead not distinguished:	",
" sirdepidemic h=0.1 stepmax=100 b=.3 k=.9 S0=100 mode=SIRD |	",
"    xgraph n=100 nplot=4 d1=1 style=normal &			",
"								",
" With 10 deaths (10% mortaliy)						",
" sirdepidemic h=0.1 stepmax=100 b=.2 nu=.1 k=.9 S0=100 mode=SIRD |",
"    xgraph n=100 nplot=4 d1=1 style=normal &			",
"								",
" sirdepidemic h=1 stepmax=50 b=0.2 k=0.1 S0=1000 I0=1 mode=SIRD |",
"     xgraph n=50 nplot=4 d1=1 style=normal			",
"								",
" Hong Kong Flu 1968-1969:					", 
" https://services.math.duke.edu/education/ccp/materials/diffcalc/sir/sir1.html",
" Population is S0=7.9 million, r0=1.5, the average period of	",
" infectiveness is  3 days so k=1/3, b=r0*k/S0=6.3e-8, and initial",
" infected is I0=10.						",
" 								",
" Without deaths, or with recovered and dead not distinguished:	",
"  sirdepidemic h=1 stepmax=200 k=.3333 b=6.3e-8 S0=7.9e6 mode=SIRD |",
"      xgraph n=200 nplot=4 d1=1 style=normal			",
"								",
" With deaths and recovered distinguished .1% mortality:	",
"  sirdepidemic h=1 stepmax=200 k=.3323 nu=.001 b=6.3e-8 S0=7.9e6 mode=SIRD |",
"      xgraph n=200 nplot=4 d1=1 style=normal			",
" (you will have to zoom in on the plot to see the 12,000 deaths",
" 								",
" Alternatively, if we normalize S0, I0, and R0 by the total	",
" population, then S0=1.0 the initial population of infectives  ",
" becomes I0=1.27e-6, the average period of infectiveness is	 ",
" 3 days making k=1/3, and the possibility of infective contact ",
" is every 2 days, making b=1/2.  Note that in the S0=1.0 case,",
" r0=1.5 = b/k.							",
"								",
" Without recovered and dead being distinguished		",
" sirdepidemic h=1 stepmax=100 b=.5 k=.333 S0=1.0 I0=1.27e-6 mode=SIRD|",
"   xgraph n=100 nplot=4 d1=1 style=normal			",
" 								",
" With recovered and dead being distinguished 0.1% mortality	",
" sirdepidemic h=1 stepmax=100 b=.5 k=.332 nu=.001 S0=1.0 I0=1.27e-6 mode=SIRD|",
"   xgraph n=100 nplot=4 d1=1 style=normal			",
" 								",
" Thus, it may be easier to determine b and k with S0=1.0, and",
" then try different I0 values					",
NULL};

/*
 * References:
 *
 * Kermack, W. O. and A. G. McKendrick (1927) A contribution to the 
 *  mathematical theory of epidemics, Procedings of the Royal Socieity A.
 *
 * The SRI model describes an epidemic in terms of
 *   S = susceptibles in a population
 *   I = infectives in a population
 *   R = recovered
 *   D = dead
 *
 *   S0 = initial value of the susceptibles
 *   I0 = initial value of the infectives
 *   R0 = initial recovered = 0
 *   D0 = initial dead = 0
 *   
 *  Always S(t) + I(t) + R(t) + D(t) = S0 + I0   
 *   
 *   r0 = b*S0/k  = basic reproduction rate
 *   b = rate of infection
 *   k = rate removal
 *   nu = mortality rate  
 *  
 *   The encounters between susceptibles and the infectives is represented
 *   by the product S*I  
 *	S'(t) =  - b*S*I   (S'(t) is always negative)
 *	I'(t) = b*S*I - k*I - nu*I
 *	R'(t) = k*I 
 *	D'(t) = nu*I
 * S(t)= susceptible members (no births, deaths, immigration, or emigration)
 * I(t)= infective number (includes asymptomatic carriers)
 * R(t)= recovered members
 * D(t)= dead members
 *
 * There is an impiled flow from S(t) -> I(t) -> R(t) or D(t), though infected
 * who are quarantined immediately become part of R(t).
 *
 * The product b*S*I denotes the interaction of the infective population with
 * the susceptible population..
 *
 * Author:  April 2020: John Stockwell
 */

/**************** end self doc ********************************/

/* Prototype of function used internally */
static int
sird_epidemic_equations(double t, double y[4] , double yprime[4]);

/* Define values of imode */
#define SIRD_MODE 0
#define S_MODE 1
#define I_MODE 2
#define R_MODE 3
#define D_MODE 4

int
main(int argc, char **argv)
{
	register int i=0, j=0;		/* counters */
	register int number=4; 		/* the three dependent variables */
	int verbose=0;		/* verbose flag =1 chatty, =0 silent */
	int stepmax=0;		/* maximum number of steps */

	/* initial values of S, I, R */
	double S0=0.0;		/* initial value of susceptible population */
	double I0=0.0;		/* initial value of infectives */
	double R0=0.0;		/* initial value of recovered */
	double D0=0.0;		/* initial value dead */

	double t=0.0;		/* time */
	double h=.001;		/* time increment */
	double tol=0.0;		/* time increment */

	double y[4]={0.0,0.0,0.0,0.0};	/* dependent variable of ODE system */

	rke_variables p;	/* variable used by RKE routines */

	float **yout=NULL; 	/* output array */
	float *tempout=NULL;	/* temporary output array */

	FILE *out_file=stdout;	/* pointer to file that we write out to */

	cwp_String mode="SIRD";	/* output mode of program */
	int imode=SIRD_MODE;	/* integer flag for mode */


	/* Hook up getpar */
	initargs(argc, argv);
	requestdoc(0);

	switch(filestat(STDOUT)) { /* Prevent floats from dumping on screen */
	case BADFILETYPE:
		warn("stdout is illegal filetype");
		pagedoc();
	break;
	case TTY:
		warn("stdout can't be tty");
		pagedoc();
	break; 
	default:			   /* rest are OK */
	break;

	}

	/* Get parameters */
	if (!getparint("stepmax", &stepmax))	stepmax = 40;
	if (!getparint("verbose", &verbose))	verbose = 0;

	/* Initial conditions y[0] = S  y[1]=I  y[2]=R */
	if (!getpardouble("S0", &S0))		S0=1000;
		y[0] = S0;
	if (!getpardouble("I0", &I0))		I0=1.0;
		y[1] = I0;
	if (!getpardouble("R0", &R0))		R0=0.0;
		y[2] = R0;
	if (!getpardouble("D0", &D0))		D0=0.01;
		y[3] = D0;

	if (!getpardouble("h", &h))		h = 1.0;
	if (!getpardouble("tol", &tol))		tol = RKE_ERR_BIAS_INIT;

	/* Get output mode, recall imode initialized to the default FABS */
	if (!getpardouble("h", &h))		h = 1.0;
	getparstring("mode", &mode);
	if (STREQ(mode, "S"))    	imode = S_MODE;
	else if (STREQ(mode, "I"))	imode = I_MODE;
	else if (STREQ(mode, "R"))      imode = R_MODE;
	else if (STREQ(mode, "D"))      imode = D_MODE;
	else if (!STREQ(mode, "SIRD"))
	    err("unknown operation=\"%s\", see self-doc", mode);

	/* allocate space in the output array */
	yout = ealloc2float(4,4*stepmax);

	/* zero out the array */
	memset((void *) yout[0], 0 , 4*stepmax*FSIZE);
	
	/* initialize Runge-Kutta-England routines */
	p = (rke_variables)
		rke_init(4, sird_epidemic_equations);

	/* set tolerance */
	p->error_bias=tol;

	for (i=0; i<stepmax; ++i) {
		double aimed_t;
		t=i*h;
		aimed_t=t+h;
  		if (verbose) {
			warn("using %3d accepted and %3d rejected steps",
	 			p->accepted_steps, p->rejected_steps);
			if (verbose) warn("error tolerance = %10.24f",p->error_bias);
		}

		/* convert doubles in y to floats in yout and write out */
		for(j=0; j<number; ++j) yout[i][j] = (float) y[j];

		/* run the Runge-Kutta-England solver */
  		rke_solve (p, &t, y, aimed_t);
	}

	/* write out according to the mode */
	tempout = ealloc1float(4*stepmax);

	if (imode==S_MODE) {
		for (i=0; i<stepmax; ++i)
			tempout[i] = yout[i][0];
	} else if (imode==I_MODE) {
		for (i=0; i<stepmax; ++i)
			tempout[i] = yout[i][1];
	} else if (imode==R_MODE) {
		for (i=0; i<stepmax; ++i)
			tempout[i] = yout[i][2];
	} else if (imode==D_MODE) {
		for (i=0; i<stepmax; ++i)
			tempout[i] = yout[i][3];
	} else if (imode==SIRD_MODE) {

		for (i=0; i<stepmax; ++i)
			tempout[i] = yout[i][0];
		for (i=0; i<stepmax; ++i)
			tempout[i+stepmax] = yout[i][1];
		for (i=0; i<stepmax; ++i)
			tempout[i+2*stepmax] = yout[i][2];
		for (i=0; i<stepmax; ++i)
			tempout[i+3*stepmax] = yout[i][3];
	}

	if (imode==SIRD_MODE) {
		efwrite(tempout,sizeof(float),4*stepmax,out_file);
	} else {
		efwrite(tempout,sizeof(float),stepmax,out_file);
	}

	/* end the session with rke */
	rke_term(p);

	return EXIT_SUCCESS;
}


static int
sird_epidemic_equations(double t, double y[4] , double yprime[4])
/*********************************************************************
sir_epidemic_equations - the system of ODE's descibing the SIRD epidemic
  model
**********************************************************************
t	independent variable "time"
y 	dependent variable being solved for y(t)
yprime	derivative of dependent variable  y'(t)
**********************************************************************
Notes: This is an example of an autonomous system of ODE's
**********************************************************************/
{
	double b=0.0;	/* infection rate	*/
	double k=0.0;	/* recovery rate	*/
	double nu=0.0;	/* mortality rate	*/

	/* parameters */
	if (!getpardouble("b", &b))		b = 0.5;
	if (!getpardouble("k", &k))		k = 0.333;
	if (!getpardouble("nu", &nu))		nu = 0.0;

	yprime[0] =  -b*y[0]*y[1] ;
	yprime[1] = b*y[0]*y[1]  - k*y[1] -  nu*y[1];
	yprime[2] = k*y[1]; 
	yprime[3] = nu*y[1]; 

    return 1;
}
