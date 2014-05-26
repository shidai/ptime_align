#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "fitsio.h"
//#include "ptime.h"

#define NP 2048
#define PI 3.14159265359

double read_psrfreq ( char *name );

int get_nchan ( char *name );
int get_npol ( char *name );
int get_nphase ( char *name );
int get_subint ( char *name );

int read_scl ( char *name, int subint, double *scl, int nchan, int npol);
int read_offs ( char *name, int subint, double *offs, int nchan, int npol);
int read_value ( char *name, int subint, double *value, int nphase, int nchan, int npol);

int read_prof ( char *name, int subint, double *profile, int nphase, int npol, int nchan);

int read_r (char *name, double *profile, int nphase);

int find_peak (int n, double *s, int *position);

double find_peak_value (int n, double *s);

int corr (double *s, double *p, int nphase);

int on_pulse (int nphase, int peak_position, double *in, double *out, double frac);

int def_off_pulse (int nphase, double *I, double frac_off);

int off_pulse (int nphase, int index, double *in, double *out, double frac_off);

int remove_baseline (double *in, int index, double frac_off, int n, double *out);

int pre_diff (double *s, int nphase, int index, double frac_off, double *s_out);

double devi (double *in, int index, int nphase, double frac_off);

int cal_L (double *I, double *Q, double *U, int nphase, int index, double frac_off, double *L);

int corr_V (double *I, double *V, int nphase, int index, double frac_off, double *V_corr);

int align_main (char *rname, char *fname, int nchn_use, double frac_off, double rot);

/////////////////////////////////////////////////////////////////
//
double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int dft_profiles (int N, double *in, fftw_complex *out);

double zbrent(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn, double *real_p, double *ima_p);

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new);

int align (int N, double phase, double b, double a, double *real_p, double *real_p_align, double *ima_p, double *ima_p_align, double rot);

int align_prof (double *s, double *p, double *Q, double *U, double *V, double psrfreq, int nphase, double rot, double frac_off);

int rotate (int N, double *real_p, double *real_p_rotate, double *ima_p, double *ima_p_rotate, double rot);

int preA7_QUV (double *p, int nphase, int nchn, double *real_p, double *ima_p);


