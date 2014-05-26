// functions used to derive phase shifts, according to Taylor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include "ptime_time.h"
#include "ptime_align.h"
//#include "nrutil.h"
#define ITMAX 100000  // Maximum allowed number of iterations.
#define EPS 1.0e-16 // Machine double floating-point precision.
//#define EPS 3.0e-8 // Machine floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double pi=3.1415926;

double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A7 in Taylor 92
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < num; j++)
	    {
			A7+=(j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	return A7;
}

double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A9 in Taylor 92
{
	double A9=0.0, sum=0.0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    A9+=a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase);
		    sum+=a_s[i][j]*a_s[i][j];
		    //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	A9=A9/sum;

	return A9;
}

int dft_profiles (int N, double *in, fftw_complex *out)
// dft of profiles
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	//double *in;
	//fftw_complex *out;
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

double zbrent(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
//	Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, a_s, a_p, p_s, p_p, num, nchn),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

int align (int N, double phase, double b, double a, double *real_p, double *real_p_align, double *ima_p, double *ima_p_align, double rotate)
{
	// k is the dimention of amp, N is the dimention of s
	int i;

	// for substraction 
	double amp,cosina,sina;
	for (i=0;i<N/2+1;i++)
	{
		// calculate the sin(phi) and cos(phi) of the profile
		amp=sqrt(real_p[i]*real_p[i]+ima_p[i]*ima_p[i]);
		cosina=real_p[i]/amp;
		sina=ima_p[i]/amp;

		// add phase shift to the profile, phase
		//real_p_align[i]=amp*(cosina)/b;
		//ima_p_align[i]=amp*(sina)/b;
		//real_p_align[i]=amp*(cosina*cos(-i*phase)-sina*sin(-i*phase));
		//ima_p_align[i]=amp*(sina*cos(-i*phase)+cosina*sin(-i*phase));
		//real_p_align[i]=amp*(cosina*cos(-i*phase)-sina*sin(-i*phase))/b;
		//ima_p_align[i]=amp*(sina*cos(-i*phase)+cosina*sin(-i*phase))/b;
		real_p_align[i]=(amp*(cosina*cos(-i*(phase+rotate*pi))-sina*sin(-i*(phase+rotate*pi))))/b;
		ima_p_align[i]=(amp*(sina*cos(-i*(phase+rotate*pi))+cosina*sin(-i*(phase+rotate*pi))))/b;
		
	}

	return 0;
}

int rotate (int N, double *real_p, double *real_p_rotate, double *ima_p, double *ima_p_rotate, double rot)
{
	// k is the dimention of amp, N is the dimention of s
	int i;

	// for substraction 
	double amp,cosina,sina;
	for (i=0;i<N/2+1;i++)
	{
		// calculate the sin(phi) and cos(phi) of the profile
		amp=sqrt(real_p[i]*real_p[i]+ima_p[i]*ima_p[i]);
		cosina=real_p[i]/amp;
		sina=ima_p[i]/amp;

		// rotate profile
		real_p_rotate[i]=amp*(cosina*cos(-i*rot*pi)-sina*sin(-i*rot*pi));
		ima_p_rotate[i]=amp*(sina*cos(-i*rot*pi)+cosina*sin(-i*rot*pi));
		//real_p_rotate[i]=amp*(cosina*cos(-i*pi)-sina*sin(-i*pi));
		//ima_p_rotate[i]=amp*(sina*cos(-i*pi)+cosina*sin(-i*pi));
		
	}

	return 0;
}

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new)
{
	double *dp;
    fftw_plan plan;
	fftw_complex *cp;

    dp = (double *)malloc(sizeof (double) * ncount);
	cp = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp, 0, sizeof (double) * ncount);
	memset(cp, 0, sizeof (fftw_complex) * ncount);

	// initialize the dft...
	double *dp_t;
    fftw_plan plan_t;
	fftw_complex *cp_t;

    dp_t = (double *)malloc(sizeof (double) * ncount);
	cp_t = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp_t, 0, sizeof (double) * ncount);
	memset(cp_t, 0, sizeof (fftw_complex) * ncount);

	int i;
    double real,ima,amp,cosina,sina;

	for (i = 0; i < ncount; i++)
	{
		if (i < ncount/2+1)
		{
            real = real_p[i];
            ima = ima_p[i];
			amp = sqrt(real*real+ima*ima);
			cosina = real/amp;
			sina = ima/amp;

			cp[i][0] = amp*(cosina);
			cp[i][1] = amp*(sina);
			//cp[i][0] = amp*(cosina*cos(-i*3.1415926)-sina*sin(-i*3.1415926));
			//cp[i][1] = amp*(sina*cos(-i*3.1415926)+cosina*sin(-i*3.1415926));
			//cp[i][0]=real_s[i]-real_p[i];
			//cp[i][1]=ima_s[i]-ima_p[i];
			//cp[i][0]=-real_s[i]+real_p[i];
			//cp[i][1]=-ima_s[i]+ima_p[i];
			cp_t[i][0] = real_p[i];
			cp_t[i][1] = ima_p[i];
			//cp[i][0]=real_p[i];
			//cp[i][1]=ima_p[i];
		}
		else
		{
			cp[i][0]=0.0;
			cp[i][1]=0.0;
			cp_t[i][0]=0.0;
			cp_t[i][1]=0.0;
		}
	}

    plan_t = fftw_plan_dft_c2r_1d(ncount, cp_t, dp_t, FFTW_MEASURE);

    fftw_execute(plan_t);

    fftw_destroy_plan(plan_t);

	/////////////////////////////////////////////////////////////////

    plan = fftw_plan_dft_c2r_1d(ncount, cp, dp, FFTW_MEASURE);

    fftw_execute(plan);

    fftw_destroy_plan(plan);

	for (i = 0; i < ncount; i++)
	{
		p_new[i] = dp[i]/ncount;  // normalized by the ncount
		//printf ("%lf\n", p_new[i]);
	}

	return 0;
}

int align_prof (double *s, double *p, double *Q, double *U, double *V, double psrfreq, int nphase, double rot, double frac_off)
// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
{
    //int nphase=1024;
    int nchn=1;

	// dft profile and template
	
	//nchn = n/nphase;
	//printf ("%d\n", nchn);
	int k;  // k=nphase/2

	double real_p[NP], ima_p[NP];
	double amp_s[nchn][NP],amp_p[nchn][NP];  // elements for calculating A7
	double phi_s[nchn][NP],phi_p[nchn][NP];  // the second dim should be NP, which is large enough for different observations

	preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, nphase, nchn, real_p, ima_p);
	//printf ("%d\n", nchn);
	
	// fft s, Q, U, V
	double s_real[NP], s_ima[NP];
	double Q_real[NP], Q_ima[NP];
	double U_real[NP], U_ima[NP];
	double V_real[NP], V_ima[NP];
	preA7_QUV (s, nphase, nchn, s_real, s_ima);
	preA7_QUV (Q, nphase, nchn, Q_real, Q_ima);
	preA7_QUV (U, nphase, nchn, U_real, U_ima);
	preA7_QUV (V, nphase, nchn, V_real, V_ima);

	// initial guess of the phase
    int peak_s, peak_p;	

	find_peak(nphase,s,&peak_s);
	find_peak(nphase,p,&peak_p);

	int d;
	double step;
	double ini_phase,up_phase,low_phase;

	d=corr (s,p,nphase);
	//printf ("%d\n", d);
	//d=peak_p-peak_s;
	//printf ("%d\n", d);
	step=2.0*3.1415926/(10.0*nphase);
	//step=2.0*3.1415926/10240.0;

	if (d>=nphase/2)
	{
		ini_phase=2.0*3.1415926*(nphase-1-d)/nphase;
		//ini_phase=2.0*3.1415926*(1023-d)/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)*A7(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*3.1415926*d/nphase;
		//ini_phase=-2.0*3.1415926*d/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)*A7(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	//printf ("%d %lf %lf %lf\n", d, ini_phase, up_phase, low_phase);

    // calculate phase shift, a and b
    double phase,b;
    phase=zbrent(A7, low_phase, up_phase, 1.0e-16, amp_s, amp_p, phi_s, phi_p, k, nchn);
    //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
    //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
    b=A9(phase, amp_s, amp_p, phi_s, phi_p, k, nchn);

	double a;
    a = (amp_p[0][0]-amp_s[0][0]);
    //a = (amp_p[0][0]-b*amp_s[0][0]);
    //a = (amp_p[0][0]-b*amp_s[0][0])/nphase;
	
	/*
	printf ("phase shift: %.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	printf ("phase shift: %.10lf\n", phase);  // microseconds
	printf ("b: %lf\n", b);  // microseconds
	printf ("a: %lf\n", a);  // microseconds
	printf ("amp_p[0]: %lf\n", amp_p[0][0]);  // microseconds
	printf ("amp_s[0]: %lf\n", amp_s[0][0]);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	*/
	
	// rotate the profile by pi
	double real_s_rotate[nphase/2+1], ima_s_rotate[nphase/2+1];
	rotate (nphase, s_real, real_s_rotate, s_ima, ima_s_rotate, rot);

	// align profile and template
	double real_p_align[nphase/2+1], ima_p_align[nphase/2+1];
	double real_Q_align[nphase/2+1], ima_Q_align[nphase/2+1];
	double real_U_align[nphase/2+1], ima_U_align[nphase/2+1];
	double real_V_align[nphase/2+1], ima_V_align[nphase/2+1];
	align (nphase, phase, b, a, real_p, real_p_align, ima_p, ima_p_align, rot);
	align (nphase, phase, b, a, Q_real, real_Q_align, Q_ima, ima_Q_align, rot);
	align (nphase, phase, b, a, U_real, real_U_align, U_ima, ima_U_align, rot);
	align (nphase, phase, b, a, V_real, real_V_align, V_ima, ima_V_align, rot);

	double s_new[nphase];
	double p_new[nphase];
	double Q_new[nphase];
	double U_new[nphase];
	double V_new[nphase];
	inverse_dft (real_s_rotate, ima_s_rotate, nphase, s_new);
	inverse_dft (real_p_align, ima_p_align, nphase, p_new);
	inverse_dft (real_Q_align, ima_Q_align, nphase, Q_new);
	inverse_dft (real_U_align, ima_U_align, nphase, U_new);
	inverse_dft (real_V_align, ima_V_align, nphase, V_new);

	/*
	// open file to write toa 
	char head[] = "Diff_";
	char channel[] = "_nchn_";
	char pol[] = "_npol_";
	char sub[] = "_nsub_";
	char c1[10], c2[10], c3[10]; 

	char output[100];

	sprintf (c1, "%d", nsub);
	sprintf (c2, "%d", nchan);
	sprintf (c3, "%d", npol);

	strcpy(output, head);
	strcat(output, fname);
	strcat(output, sub);
	strcat(output, c1);
	strcat(output, channel);
	strcat(output, c2);
	strcat(output, pol);
	strcat(output, c3);

	FILE *fp;
	if ((fp = fopen(output, "w+")) == NULL)
	{
        fprintf (stdout, "Can't open file\n");
		exit(1);
	}
	*/

	int index;
	index = def_off_pulse (nphase, p_new, frac_off);

	//double L_temp[nphase];
	double L_new[nphase];
	cal_L (p_new, Q_new, U_new, nphase, index, frac_off, L_new);
	
	// correct V
	double V_corr[nphase];
	corr_V (p_new, V_new, nphase, index, frac_off, V_corr);

	int i;
	for (i = 0; i < nphase; i++)
	{
		//fprintf (fp, "%.3lf %.3lf\n", p_new[i]+a, p[i]);
		printf ("%.6lf %.6lf %.6lf %.6lf %.6lf\n", p_new[i], Q_new[i], U_new[i], V_corr[i], L_new[i]);
		//printf ("%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n", s[i], s_new[i], p_new[i], Q_new[i], U_new[i], V_new[i]);
	}

	double sigma_I, sigma_Q, sigma_U, sigma_V, sigma_L;
	sigma_I = devi (p_new, index, nphase, frac_off);
	sigma_Q = devi (Q_new, index, nphase, frac_off);
	sigma_U = devi (U_new, index, nphase, frac_off);
	sigma_V = devi (V_corr, index, nphase, frac_off);
	sigma_L = devi (L_new, index, nphase, frac_off);

	printf ("%.8lf %.8lf %.8lf %.8lf %.8lf\n", sigma_I, sigma_Q, sigma_U, sigma_V, sigma_L);

	/*
    if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");
	*/

	return 0;
}

int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn, double *real_p, double *ima_p)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[nphase];  // initialize the system, don't know why....

	for (i=0;i<nphase;i++)
	{
		test[i]=s[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	dft_profiles(nphase,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

    fftw_complex *out_s;
	fftw_complex *out_p;
	
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double s_temp[nphase];  // store one template and profile
	double p_temp[nphase];  

	int n;
	double r_s[nphase/2],im_s[nphase/2];
	double r_p[nphase/2],im_p[nphase/2];
	for (i = 0; i < nchn; i++)
	{
	    for (j=0;j<nphase;j++)
	    {
		    s_temp[j]=s[i*nphase + j];
		    p_temp[j]=p[i*nphase + j];
	    }

	    dft_profiles(nphase,s_temp,out_s);
	    //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	    dft_profiles(nphase,p_temp,out_p);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];

		for (j = 0; j < nphase/2+1; j++)                                                  
		{                                                                      
			//real_s[j]=out_s[j][0];                                             
	        //ima_s[j]=out_s[j][1];                                              
			real_p[j]=out_p[j][0];                                             
			ima_p[j]=out_p[j][1];                                              
		}
										
		n = 0;
	    for (j = 0; j <= nphase/2-1; j++)
	    {
		    r_s[j]=out_s[j+1][0];
		    im_s[j]=out_s[j+1][1];
		    r_p[j]=out_p[j+1][0];
		    im_p[j]=out_p[j+1][1];
		    //printf ("%lf %lf\n", r_p[i], im_p[i]);
		    //printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		    n++;
	    }
	    //printf ("%d\n", n);
	    //printf ("%d %d\n", nphase, nchn);

	    for (j = 0; j < n; j++)
	    {
		    amp_s[i][j]=sqrt(r_s[j]*r_s[j]+im_s[j]*im_s[j]);
		    amp_p[i][j]=sqrt(r_p[j]*r_p[j]+im_p[j]*im_p[j]);
		    phi_s[i][j]=atan2(im_s[j],r_s[j]);
		    phi_p[i][j]=atan2(im_p[j],r_p[j]);
		    //printf ("%lf %lf %lf\n", r_s[i], im_s[i], amp_s[i]);
		    //printf ("%lf %lf %lf\n", r_p[i], im_p[i], amp_p[i]);
		    //printf ("%lf\n", amp_s[i]);
		    //printf ("%lf\n", amp_p[i]);
	    }
	}
	(*k)=n;

	fftw_free(out_s); 
	fftw_free(out_p); 
	fftw_free(out_t); 

	return 0;
}

int preA7_QUV (double *p, int nphase, int nchn, double *real_p, double *ima_p)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[nphase];  // initialize the system, don't know why....

	for (i=0;i<nphase;i++)
	{
		test[i]=p[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	dft_profiles(nphase,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

	fftw_complex *out_p;
	
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double p_temp[nphase];  // store one template and profile

	for (i = 0; i < nchn; i++)
	{
	    for (j=0;j<nphase;j++)
	    {
		    p_temp[j]=p[i*nphase + j];
	    }

	    dft_profiles(nphase,p_temp,out_p);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];

		for (j = 0; j < nphase/2+1; j++)                                                  
		{                                                                      
			real_p[j]=out_p[j][0];                                             
			ima_p[j]=out_p[j][1];                                              
		}
										
	}

	fftw_free(out_p); 
	fftw_free(out_t); 

	return 0;
}

