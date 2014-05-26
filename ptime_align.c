#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptime_align.h"

int main (int argc, char *argv[])
{
	//int h,i,j,k;

	//////////////////////////////////////////////////////
	char fname[128];   // name of data file
	char rname[128];   // name of reference profile

	int nchn_use;  // sub-channel index
	double frac_off, rotate;

	int i;

	for (i = 1; i < argc; i++)
    {
		if (strcmp(argv[i],"-r") == 0)
		{
			strcpy(rname,argv[++i]);
			//printf ("Profiles will aligned to %s\n", rname);
		}
		else if (strcmp(argv[i],"-f") == 0)
		{
			strcpy(fname,argv[++i]);
			nchn_use = atof(argv[++i]);
			//printf ("The sub-channel %d of %s will be aligned and scaled to %s\n", nchn_use, fname, rname);
		}
		else if (strcmp(argv[i],"-frac") == 0)
		{
			frac_off = atof(argv[++i]);
		}
		else if (strcmp(argv[i],"-rotate") == 0)
		{
			rotate = atof(argv[++i]);
		}
		else 
		{
			printf ("Wrong options!!!\n");
			printf ("Usage: ptime_align.out -r rname -f fname nchn -frac on off -rotate rot\n"
					"       align and scale profiles.\n"
					"       -r rname: the reference profile;\n" 
					"       -f fname: data file;\n"
					"       -frac off: off pulse fraction;\n"
					"       -rotate rot: fraction of pi to rotate profile;\n");
			exit (0);
		}
    }
	//printf ("%d\n", smode);

	// start to deal with different data file
	//
	// open file to write toa 
	//FILE *fp;
	//if ((fp = fopen(oname, "w+")) == NULL)
	//{
    //	fprintf (stdout, "Can't open file\n");
	//	exit(1);
	//}
    //fprintf (fp, "S0    S    err\n");
	/////////////////////////////////////////////////////////
	
	align_main (rname, fname, nchn_use, frac_off, rotate);

    //if (fclose (fp) != 0)
	//	fprintf (stderr, "Error closing\n");

	return 0;
}
