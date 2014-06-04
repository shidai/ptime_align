// calculate the differences between profiles
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "ptime_align.h"

int find_peak (int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}
	peak = temp[n-1];

	for (i = 0; i < n; i++)
	{
		if (fabs(peak-s[i]) < 1.0e-3)
		//if (fabs(peak-fabs(s[i])) < 1.0e-3)
		{
			(*position) = i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		c = (a >= b ? a : b);
		//c = (fabs(a) >= fabs(b) ? a : b);

		temp[i+1] = c;
	}

	return temp[n-1];
}

int corr (double *s, double *p, int nphase)
{
	/*
	FILE *fp1, *fp2;

	if ((fp1 = fopen(argv[1], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	if ((fp2 = fopen(argv[2], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	float x1[1024], x2[1024];
	int i = 0;
	while (fscanf (fp1, "%f", &x1[i]) == 1)
	{
		i++;
	}

	i = 0;
	while (fscanf (fp2, "%f", &x2[i]) == 1)
	{
		i++;
	}
	*/

	double r[nphase];
	int i, j;
	for (j = 0; j < nphase; j++)
	{
		r[j] = 0.0;
		for (i = 0; i < nphase; i++)
		{
			if ((i+j) > (nphase-1))
			{
				r[j] += p[i]*s[i+j-(nphase-1)];
			}
			else
			{
				r[j] += p[i]*s[i+j];
			}
			//printf ("%f %f\n", x1[i], x2[i]);
		}
	}

	int shift;
	find_peak (nphase, r,  &shift);
	/*
	for (j = 0; j < 1024; j++)
	{
		printf ("%f\n", r[j]);
	}
	*/

	return -shift;
}

int on_pulse (int nphase, int peak_position, double *in, double *out, double frac)
// define the on_pulse range, frac_on of the phase
{
	int n = nphase;
	int num = (int)(n*frac);
	int i;
	for (i = 0; i < num; i++)
	{
		if ((peak_position-(int)(num/2)+i) < 0)
		{
			out[i] = in[(n-1)+(peak_position-(int)(num/2)+i)];
		}
		else if ((peak_position-(int)(num/2)+i) > n-1)
		{
			out[i] = in[(peak_position-(int)(num/2)+i)-(n-1)];
		}
		else
		{
			out[i] = in[peak_position-(int)(num/2)+i];
		}
	}

	return 0;
}

int def_off_pulse (int nphase, double *in, double frac_off)
// define the off pulse region based on I, return the starting index of off pulse region
// using frac_off to calculate the off pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i,j;
	double small;
	double temp;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num_off; j++)
			{
				small += (in[j]+30000.0)*(in[j]+30000.0);
			}
			small = sqrt(small/num_off);
		}
			
		temp = 0.0;
		for(j = 0; j < num_off; j++)
		{
			if ((i+j) > n-1)
			{
				temp += (in[(i+j)-(n-1)]+30000.0)*(in[(i+j)-(n-1)]+30000.0);
			}
			else 
			{
				temp += (in[i+j]+30000.0)*(in[i+j]+30000.0);
			}
		}
		temp = sqrt(temp/num_off);

		small = (temp <= small ? temp : small);
		index = (temp <= small ? i : index);
		//printf ("%d %lf \n", index, small);
	}
	//printf ("%d\n",index);

	return index;
}

int off_pulse (int nphase, int index, double *in, double *out, double frac_off)
// get the off_pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i;

	for (i = 0; i < num_off; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int remove_baseline (double *in, int index, double frac_off, int n, double *out)
{
	// define the off_pulse range, frac_off is the fraction of the phase
	// index is the starting point of the off_pulse range
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse (n, index, in, off_0, frac_off);

	int i;
	double baseline = 0.0;
    for (i = 0; i < num_off; i++)
    {
        baseline += off_0[i];
        //average_s += s_off[i];
    }
	baseline = baseline/num_off;

    //printf ("the baseline of std is: %lf \n", baseline);
    //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

int pre_diff (double *s, int nphase, int index, double frac_off, double *s_out)
{
	int n = nphase;
	
	// remove the baseline
	remove_baseline (s, index, frac_off, n, s_out);

    return 0;
}

double devi (double *in, int index, int nphase, double frac_off)
{
	int n = nphase;

	// define the off_pulse range, frac_off is the fraction of the phase
	//double frac_off = 0.6;
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse (n, index, in, off_0, frac_off);

	int i;

	double mean = 0.0;
    for (i = 0; i < num_off; i++)
    {
		mean += off_0[i];
	}
	mean = mean/num_off;

	double sigma = 0.0;
    for (i = 0; i < num_off; i++)
    {
		sigma += (off_0[i]-mean)*(off_0[i]-mean);
	}
	sigma = sqrt(sigma/num_off);
	
	return sigma;
}

int cal_L (double *I, double *Q, double *U, int nphase, int index, double frac_off, double *L)
{
	int i;
	double x;
	double sigma_I;

	sigma_I = devi (I, index, nphase, frac_off);

	for (i = 0; i < nphase; i++)
	{
		//L[i] = sqrt(Q[i]*Q[i] + U[i]*U[i]);
		x = sqrt(Q[i]*Q[i] + U[i]*U[i]);

		if (x/sigma_I >= 1.57)
		{
			L[i] = sigma_I*sqrt((x/sigma_I)*(x/sigma_I) - 1);
		}
		else
		{
			L[i] = 0.0;
		}
	}

	return 0;
}

int corr_V (double *I, double *V, int nphase, int index, double frac_off, double *V_corr)
{
	int i;
	double sigma_I;

	sigma_I = devi (I, index, nphase, frac_off);

	for (i = 0; i < nphase; i++)
	{
		//V_corr[i] = V[i];
		if (fabs(V[i])/sigma_I > 1.0)
		{
			if (V[i] < 0.0)
			{
				V_corr[i] = -pow(pow(V[i],4.0)-pow(sigma_I,4.0),0.25);
			}
			else
			{
				V_corr[i] = pow(pow(V[i],4.0)-pow(sigma_I,4.0),0.25);
			}
		}
		else
		{
			V_corr[i] = 0.0;
		}
	}

	return 0;
}

int align_main (char *rname, char *fname, int nchn_use, double frac_off, double rotate)
//int real_obs (char *fname, char *tname, char *oname, int mode, FILE *fp, double frac_on, double frac_off)
{
	// name of different extension of data files
	char name_data[50]; 
	char name_psrparam[50]; 

	char data[] = "[SUBINT]";
	char psrparam[] = "[PSRPARAM]";

	// read name of reference profile
	char r[50];
	strcpy(r, rname);
	strcat(r, data);

	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	int h, j, p;
	{
		// name of different extension
		strcpy(name_data,fname);
		strcpy(name_psrparam,fname);

		strcat(name_data, data);
		strcat(name_psrparam, psrparam);

		////////////////////////////////////////////////////
		
		double psrfreq;
		psrfreq = read_psrfreq(name_psrparam);
		//printf ("PSR frequency: %.15lf\n", psrfreq);

		////////////////////////////////////////////////////

		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		nchn = get_nchan(name_data);	
		npol = get_npol(name_data);	
		nsub = get_subint(name_data);	
		nphase = get_nphase(name_data);	

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

		// read a reference profile
		double r_prof[nphase];
		double r_new[nphase];

		read_r(r,r_prof,nphase);

		int index_r;
		index_r = def_off_pulse (nphase, r_prof, frac_off);

		pre_diff (r_prof, nphase, index_r, frac_off, r_new);
		////////////////////////////////////////////////////////////////////////////////

		double p_multi[nchn*npol*nphase];
		double I_temp[nphase];
		double Q_temp[nphase];
		double U_temp[nphase];
		double V_temp[nphase];

		// start to calculate shape paremeter for different subint
		for (h = 1; h <= nsub; h++)
		{
			// simulate data

			//SNR = 500.0 + 200.0*i;
			//simulate(n,SNR,s,p_temp);

			// read profiles from data file
			read_prof(name_data,h,p_multi,nphase,npol,nchn);
			//readfile(argv[2],&n,tt,p_multi);

			// start to align profiles
			//for (i = 0; i < npol; i++)
			{
				//for (p = nchn_start-1; p < nchn_end-1; p++)
				p = nchn_use-1;
				{
					for (j = 0; j < nphase; j++)
					{
						I_temp[j] = p_multi[0*nchn*nphase + p*nphase + j];
						//printf ("%lf\n", I_temp[j]);
						Q_temp[j] = p_multi[1*nchn*nphase + p*nphase + j];
						U_temp[j] = p_multi[2*nchn*nphase + p*nphase + j];
						V_temp[j] = p_multi[3*nchn*nphase + p*nphase + j];
					}

					int index;
					index = def_off_pulse (nphase, I_temp, frac_off);

					double I_new[nphase];
					double Q_new[nphase];
					double U_new[nphase];
					double V_new[nphase];
					pre_diff (I_temp, nphase, index, frac_off, I_new);
					pre_diff (Q_temp, nphase, index, frac_off, Q_new);
					pre_diff (U_temp, nphase, index, frac_off, U_new);
					pre_diff (V_temp, nphase, index, frac_off, V_new);
					//printf ("################################\n");

					/*
					for (j = 0; j < nphase; j++)
					{
						printf ("%lf %lf\n", I_new[j], r_new[j]);
					}
					*/

					align_prof (r_new, I_new, Q_new, U_new, V_new, psrfreq, nphase, rotate, frac_off);
				}
			}
		}
	}

	return 0;
}

