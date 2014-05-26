// read PSRFITS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "ptime_align.h"
#include "ptime.h"

double read_psrfreq ( char *name )
//int main(int argc, char *argv[])
{
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
	
	char freq[100];
	char F0[100];
	double a[10];

	int colnum = 1;
    int frow;
    int	felem = 1;
    int nelem = 1;
    int	anynull = 0;
    char nval[]="NULL";

	char **line;
	line = (char **)malloc(sizeof(char *));
	line[0] = (char *)malloc(sizeof(char)*1024);

	for (frow = 1; frow < 20; frow++)
	{
		fits_read_col(fptr, TSTRING, colnum, frow, felem, nelem, nval, line, &anynull, &status);           // read the column

		//puts(line[0]);

		int nchar = strlen(line[0]);
		//printf("strlen %d\n", nchar);

		if (line[0][0] == 'F' && line[0][1] == '0')
		{
			int i;
			for (i = 0; i < nchar; i++)
			{
				F0[i] = line[0][i];
			}
			F0[nchar] = '\0';
			//printf("F0 %s\n", F0);

			int l = 0;
			int j = 0;
			for (i = 0; i < nchar+1; i++)
			{
				if( (F0[i] >= '0' && F0[i] <= '9') || F0[i] =='.' ) 
				{ 
					freq[l] = F0[i];
					l++;
				}
				else if(l > 0)
				{
					freq[l] = '\0';
					a[j]=atof(freq);
					j++;
					l=0;
				}
			}
			break;
		}
	}

	double psrfreq;
	psrfreq = a[1];

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file\n " );
    }

	return psrfreq;
}

int get_nchan ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	int nchan;
    if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("number of nchan: %d\n", nchan);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nchan;
}

int get_npol ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
	int npol;
    if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("number of npol: %d\n", npol);

	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return npol;
}

int get_nphase ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
	int nbin;
    if ( fits_read_key(fptr, TINT, (char *)"NBIN", &nbin, NULL, &status) )           // get the row number
    {
        printf( "error while getting the nbin number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("number of nbin: %d\n", nbin);

	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nbin;
}

int get_subint ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("number of subint: %ld\n", nrows);
    
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nrows;
}

int read_scl ( char *name, int subint, double *scl, int nchan, int npol)
//int main (int argc, char *argv[] )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status) )           // get the colnum number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = nchan*npol;
    int null = 0;
    int	anynull = 0;

	//int subint = 1;
	//int nchan = 8;
	//double wts[nchan];
	frow = subint;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, scl, &anynull, &status);           // read the column

	//int i;
    //for (i = 0; i < nchan; i++)                             // print the results
	//{
		//puts(line[0]);
    //    printf("%lf\n", wts[i]);
		//fprintf (fp, "%s\n", line[0]);
	//}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return 0;
}

int read_offs ( char *name, int subint, double *offs, int nchan, int npol)
//int main (int argc, char *argv[] )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status) )           // get the colnum number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = nchan*npol;
    int null = 0;
    int	anynull = 0;

	//int subint = 1;
	//int nchan = 8;
	//double wts[nchan];
	frow = subint;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, offs, &anynull, &status);           // read the column

	//int i;
    //for (i = 0; i < nchan; i++)                             // print the results
	//{
		//puts(line[0]);
    //    printf("%lf\n", wts[i]);
		//fprintf (fp, "%s\n", line[0]);
	//}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return 0;
}

//int main ( int argc, char *argv[] )
int read_value ( char *name, int subint, double *value, int nphase, int nchan, int npol)
{  
//double *read_arrival_time( char *input, long *nrows )
    //int subint = 1;
	//double profile[8*1024];
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    //if ( fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status) )           // get the row number
    if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

	///////////////////////////////////////////////////////////////////////////

	int nbin;
    int frow;
    int felem;
    int nelem;
    int null;
    int anynull;
    //double *profile;     // the array to store the profile   

	nbin = nphase;
    //profile = ( double *)malloc( (nchan*npol*nbin) * sizeof( double ) );               // allocate space for column value
    frow = subint;
    felem = 1;
    nelem = nbin*nchan*npol;
    //nelem = 1024;
    null = 0;
    anynull = 0;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, value, &anynull, &status);           // read the column

	//int i;
    //for (i = 0; i < (nchan*npol*nbin); i++)                             // print the results
    //    printf("%d %lf \n", i, profile[i]);

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return 0;
}

//int main ( int argc, char *argv[] )
int read_prof ( char *name, int subint, double *profile, int nphase, int npol, int nchan)
{  
	double value[nphase*npol*nchan];
	double scl[npol*nchan];
	double offs[npol*nchan];

	read_value (name, subint, value, nphase, nchan, npol);
	read_offs (name, subint, offs, nchan, npol);
	read_scl (name, subint, scl, nchan, npol);

	int i,j,h;
    for (i = 0; i < npol; i++)                             // print the results
    //for (i = 0; i < nchan; i++)                             // print the results
	{
		for (h = 0; h < nchan; h++)                             // print the results
		//for (h = 0; h < npol; h++)                             // print the results
		{
			for (j = 0; j < nphase; j++)                             // print the results
			{
				//profile[i*nchan*nphase+h*nphase+j] = value[i*nchan*nphase+h*nphase+j];
				profile[i*nchan*nphase+h*nphase+j] = value[i*nchan*nphase+h*nphase+j]*scl[i*nchan+h] + offs[i*nchan+h];
				//profile[i*npol*nphase+h*nphase+j] = value[i*npol*nphase+h*nphase+j]*scl[i*npol+h] + offs[i*npol+h];
				//profile[i*nphase+j] = value[i*nphase+j];
				//printf("%lf \n", profile[i*nchan*nphase+h*nphase+j]);
			}
		}
	}

    return 0;
}

int read_r (char *name, double *profile, int nphase)
{  
	int j;
	{
		fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
		int status;

		status = 0;

		//if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
		if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
		{
			printf( "error while openning file\n" );
		}

		//////////////////////////////////////////////////////////////////////////
		int npola;
		if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npola, NULL, &status) )           // get the pol number
		{
			printf( "error while getting the npol number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", npol);

		int nchan;
		if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the chan number
		{
			printf( "error while getting the nchan number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", nchan);
		///////////////////////////////////////////////////////////////////////////

		if ( fits_close_file(fptr, &status) )
		{
			printf( " error while closing the file\n" );
		}

		int subint = 1;
		// read the data
		double temp[npola*nphase*nchan];
		read_prof (name, subint, temp, nphase, npola, nchan);

		for (j = 0; j < nphase; j++)
		{
			profile[j] = temp[j];
		}
	}

	return 0;
}

