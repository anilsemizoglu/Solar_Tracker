#include "dependencies.h"       /* must precede my other headers */
#include "sofajpl.h"            /* SOFA and JPL functions */
#include "sofajpl_sup1.h"       /* supplementary routines */
#include "hipparcos.h"          /* Hipparcos star catalog routines */

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <math.h>
using namespace std;

const double PI =3.14159265359 ;
const double R_SUN = 696000.;           /* Sun radius, km */
const double R_EARTH = 6378.14;         /* Earth radius */
const char *ephem = "binp2000.406";

void 
psh_dsex_u( 
const int dmsf[],       /* deg (or hr), min, sec, dec. sec */ 
int skip, 
const int precision, 
const char *str)        /* terminating string: "  ", "\n", etc. */

/*
prints out degrees to a certain precision
*/

{
	/* The "skip" argument lets you omit fields on the left if you know
	they are always 0. */

	if (skip < 0 || skip > 3)       /* illegal value */
		skip = 0;

	/* The "precision" argument has the same meaning as in the SOFA
	routines that generate sexagesimal (A2AF etc.):
	... -1 = 10" precision, 0 = 1" precision, 1 = .1" precision, ...
	Here it's used to supress meaningless fields, e.g., seconds when
	precision is only 1 minute. It also controls the formatting of
	decimal seconds. */

	switch (skip) {
	case 0:
		printf("%d", dmsf[0]);          /* deg or hours */
		if (precision <= -4)    /* 1 deg precision */
			break;
		printf(" ");
	case 1:
		printf("%02d", dmsf[1]);        /* min */
		if (precision <= -2)    /* 1 min precision */
			break;
		printf(" ");
	case 2:
		printf("%02d", dmsf[2]);        /* sec */
		if (precision <= 0)     /* 1 sec precision */
			break;
	case 3:
		printf(".%0*d", precision, dmsf[3]);    /* dec. sec */
	}
	printf(str);
}


void
psh_dsex_s(
const char sign,
const int dmsf[],
const int skip,
const int precision,
const char *str)
/* Like the preceding function, except input is sexagesimal components and
a sign, and output is displayed with that sign. */
{
	printf("%c", sign);
	psh_dsex_u(dmsf, skip, precision, str);
}


void
psh_dangle_u(
double angle,           /* radians */
const int skip,
const int precision,
const char *str)
/* Convert radians to unsigned sexagesimal angle, print. The angle is put into
the range 0 to 2pi before conversion to degrees.

The "skip" argument is the number of fields to omit on the left. This eliminates
unnecessary leading zeros when it's known the angle is always small. Use 0 to
display everything, 1 to omit degrees, 2 to omit degrees and minutes, 3 to omit
everything except the decimal part of seconds.

Argument "precision" controls rounding of the displayed angle:
... -3 = 10', -2 = 1', -1 = 10", 0 = 1", 1 = .1" ...

The "str" is a string to print after the angle. E.g., "\n" or "  ". */
{
	int dmsf[4];
	char sign;

	/* put angle in 0 to 2pi range */
	angle = iau_anp(angle);

	iau_a2af(precision, angle, &sign, dmsf);
	psh_dsex_u(dmsf, skip, precision, str);
}
void
psh_dangle_s(
double angle,
const int skip,
const int precision,
const char *str)
/* Like psh_dangle_u(), but its input is put into the range -pi to pi before
conversion to degrees, and the signed result printed. */
{
	int dmsf[4];
	char sign;

	/* put angle in -pi to +pi range */
	angle = iau_anpm(angle);

	iau_a2af(precision, angle, &sign, dmsf);
	psh_dsex_s(sign, dmsf, skip, precision, str);
}


void
psh_dtime_u(
const double angle,
const int skip,
const int precision,
const char *str)
/* Like psh_dangle_u(), but angle is converted to unsigned sexagesimal time. */
{
	int hmsf[4];    /* hour, min, sec, decimal sec */
	char sign;      /* ignored, but the SOFA routine needs it */

	/* convert angle to sexagesimal time, result in hmsf[] */
	iau_a2tf(precision, angle, &sign, hmsf);

	/* print result */
	psh_dsex_u(hmsf, skip, precision, str);
}


void sun_equatorial(double &ra, double &dec, int y, int m, int d, int h, int min, int s)
/* Sun geocentric equatorial apparent RA, dec, HP, SD (true equinox of date */
{

	/* observer GCRS position & velocity; set to 0 to generate geocentric
	coordinates */
	const IAU_PVVEC obs_pv = { { 0., 0., 0. }, { 0., 0., 0. } };


	/* no. of decimal places in the preceding values */
	const int ndp_t = 5;    	
	const int ndp_d = 5;    	
	const int ndp_a = 2;  
	
	/* calculated values */
	double ra_out, dec_out, distance_out;

	/* temporaries */
	double tt[2];           /* 2-part Julian date (TT) */
	IAU_PVEC p_vec;         /* position vector */
	/* constants */

	/* desired accuracy in Sun apparent position = .1" (in radians) */
	const double ACCURACY = psh_mseca2rad(.1);
	const int SUN = jpl_sun;


//	printf("\nSun geocentric apparent RA, dec, HP, SD:\n");

	/* Julian Date in TT */
	double frac = h/24.0 + min/1440.0 + s/86400.0;

	psh_tdate2jd(y, m, d,frac , tt);

	/* put geocentric apparent position of Sun in p_vec */
	if (psh_jopen(ephem)) {
		printf("sun_equatorial(): can't open JPL ephemeris\n");
		return;
	}
	if (psh_jplanet(tt, SUN, 2, ACCURACY, obs_pv, p_vec)) {
		printf("sun_equatorial(): error from JPL ephemeris\n");
		return;
	}

	/* transform to true equator & equinox of date using IAU 2000B
	nutation */
	psh_ccrs2eq(tt, p_vec, 'b', 0, p_vec);


	/* Convert to spherical (RA & dec), print. */
	psh_cvec2sph(p_vec, &ra_out, &dec_out, &distance_out);
/*
	cout << y << "/"<< m << "/" <<d <<" | "<< h <<":"<< min <<":"<< s << endl;
	psh_dtime_u(ra_out, 0, ndp_t, "\n");
	cout << "ra_out" << ": " << ra_out*180/PI << endl;
	cout << endl;
	psh_dangle_s(dec_out, 0, ndp_d, "\n");
	cout << "dec_out" << ": " << dec_out*180/PI << endl;
*/
ra = ra_out*180/PI;
dec = dec_out*180/PI;
		
	/* Convert distance to km, obtaining the conversion factor from the
	ephemeris. */
//	distance_out *= psh_jconst("AU", &status);
	//if (status) {
//		printf("sun_equatorial(): AU constant not found\n");
//		return;
//	}

}
