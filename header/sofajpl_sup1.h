/*
sofajpl_sup1.h
Supplementary routines for my SOFA/JPL/Hipparcos DLL.


Copyright 2006 by Paul S. Hirose. You may use and distribute this software
(as-is or modified), but not for profit. Unless you are the sole user, I must
receive credit in a place reasonably obvious to users (e.g., by clicking
Help, About).

Furthermore, you're bound by the SOFA Software License, which I have included
in one of the files of this software. It's also available at the IAU SOFA web
site.

Furthermore, you should give credit to:
1) NASA Jet Propulsion Laboratory, for their Planetary and Lunar Ephemerides
2) ESA Hipparcos Space Astrometry Mission, for their Hipparcos Catalog
3) SIMBAD database, operated at CDS, Strasbourg, France.


See the SOFA and JPL documentation at:

http://www.iau-sofa.rl.ac.uk/
http://ssd.jpl.nasa.gov/eph_info.html

Some of these functions aren't fully understandable unless you have some
knowledge of the IAU 2003 resolutions. U.S. Naval Observatory Circular 179 is a
readable explanation:
http://aa.usno.navy.mil/publications/docs/Circular_179.html

More detailed but less readable are these documents:
http://maia.usno.navy.mil/conv2003.html


All functions in this header begin with psh_ followed by a letter indicating the
purpose of the function:
	c	coordinate conversion
	j	JPL planetary/lunar ephemeris
	m	mathematical conversions and constants
	t	time


2005-09-28  Began work.
2006-07-10  Initial release.
2006-07-25  Minor correction to comment in psh_cstar_ap().
2007-05-02  Minor clarification in psh_tdate2jd() comment.
2007-06-07  Correction to psh_tutc_tt() comment.
*/


/*############################# Test Functions ###############################*/

/*########################## coordinate conversion ###########################*/


EXTERN DLLFUNC void STDCALL
psh_ccrs2eq(
	const double tt[2],	/* 2-part Julian Date (TT) */
	constIAU_PVEC in_vec,	/* input vector in the ICRS */
	const char frame,	/* 'M' (mean) or 'A' (2000A) or 'B' (2000B) */
	const int reverse,	/* flag for reverse conversion */
	IAU_PVEC out_vec);	/* output vector */
/* Rotate the coordinate frame from the ICRS to equator & equinox of date, mean
or true, using IAU 2000 frame bias and precession, and IAU 2000A or 2000B
nutation. (Nutation is not applicable to mean coordinates.) Depending the
"frame" argument, the SOFA PNM00A, PNM00B, or BP00 routine is used.

Terrestrial Time is a 2-part Julian Date in the SOFA manner. It is automatically
normalized (the components are re-apportioned) for maximum accuracy. 

The "reverse" flag, if nonzero, reverses the conversion: from the selected
frame in in_vec to the ICRS in out_vec.

It is allowable to use the same argument as the input and output vector. Output
vector is zero if the letter for "frame" is illegal (not case sensitive). */


EXTERN DLLFUNC void STDCALL
psh_ccrs2ec(
	const double tt[2],	/* 2-part Julian Date (TT) */
	constIAU_PVEC in_vec,	/* input vector in the ICRS */
	const char frame,	/* 'M' (mean) or 'A' (2000A) or 'B' (2000B) */
	const int reverse,	/* flag for reverse conversion */
	IAU_PVEC out_vec);	/* output vector */
/* Like psh_ccrs2eq(), but ecliptical coordinates are output. I.e., converting
the output vector to polar form gives ecliptical latitude and longitude.

The "frame" argument selects the equinox: mean (IAU 2000), or the true equinox
computed with the IAU 2000A or 2000B nutation model. */


EXTERN DLLFUNC void STDCALL
psh_ccrs2ci(
	const double tt[2],	/* 2-part Julian Date (TT) */
	constIAU_PVEC in_vec,	/* input vector in the ICRS */
	const char frame,	/* 'A' (2000A) or 'B' (2000B) */
	const int reverse,	/* flag for reverse conversion */
	IAU_PVEC out_vec);	/* output vector */
/* Like psh_ccrs2eq(), but output frame is the celestial intermediate system.
There is no "mean" option in this system; it always uses the true equator. */


EXTERN DLLFUNC void STDCALL
psh_ccrs2eq_m(
	const double tt[2],	/* 2-part Julian Date (TT) */
	const char frame,	/* 'M' (mean) or 'A' (2000A) or 'B' (2000B) */
	const int reverse,	/* flag for reverse conversion */
	IAU_RMAT rot_mat);	/* output: rotation matrix */

EXTERN DLLFUNC void STDCALL
psh_ccrs2ec_m(
	const double tt[2],	/* 2-part Julian Date (TT) */
	const char frame,	/* 'M' (mean) or 'A' (2000A) or 'B' (2000B) */
	const int reverse,	/* flag for reverse conversion */
	IAU_RMAT rot_mat);	/* output: rotation matrix */

EXTERN DLLFUNC void STDCALL
psh_ccrs2ci_m(
	const double tt[2],	/* 2-part Julian Date (TT) */
	const char frame,	/* 'A' (2000A) or 'B' (2000B) */
	const int reverse,	/* flag for reverse conversion */
	IAU_RMAT rot_mat);	/* output: rotation matrix */

/* Siblings of the preceding set of functions. Instead of rotating a coordinate
frame, they output the corresponding rotation matrix. This can be more efficient
when repeated coordinate transformations have to be done and the output frame
maintains essentially the same orientation throughout. */


EXTERN DLLFUNC void STDCALL
psh_cnut00(
	const double tt[2],	/* 2-part Julian Date (Terrestrial Time) */
	const char model,	/* 'A' or 'B' */
	double *nutlo,		/* nutation in longitude */
	double *nutob);		/* nutation in obliquity */
/* Returns nutation angles in radians, given Terrestrial Time as a 2-part JD in
the SOFA format. Uses the IAU 2000A or 2000B nutation model, depending on
argument "model" (it is not case sensitive). If the letter is not one of the
allowed letters, exactly 1.0 in returned in both nutations. */


EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_cob00(
	const double jd[2],	/* 2-part Julian Date (Terrestrial Time) */
	const char model);	/* 'M', 'A', or 'B' */
/* Returns the obliquity of the ecliptic (radians), given Terrestrial Time as a
2-part Julian Date in the SOFA format. The "model" argument selects the
obliquity model: IAU 2000 mean, IAU 2000A true, or IAU 2000B true. (The letter
is not case sensitive.) If the letter is illegal, 1.0 is returned. */


EXTERN DLLFUNC void STDCALL
psh_cobsvec(
	const double lat,	/* geodetic latitude, radians */
	const double lon,	/* geodetic longitude, radians */
	const double h,		/* height above ellipsoid */
	const double a,		/* ellipsoid equatorial radius */
	const double f_inv,	/* ellipsoid 1/f */
	IAU_PVEC xyz);		/* rectangular coordinates out */
/* Converts geodetic coordinates from spherical to rectangular. Any unit of
measure may be used for arguments h and a, but it must be the same for both.
The output uses the same unit. */


EXTERN DLLFUNC void STDCALL
psh_cobspv(
	constIAU_PVEC r_obs,	/* observer rectangular coords, meters */
	constIAU_RMAT c2t_mat,	/* celestial to terrestrial matrix */
	IAU_PVVEC pv_obs);	/* output GCRS vectors */
/* Computes position and velocity of the observer in the Geocentric Celestial
Reference System with respect to the geocenter. These values are needed to get
the topocentric position of a body.

Vector r_obs is the observer's rectangular coordinates in the ITRS (or to
practical accuracy, WGS-84 and similar geodetic datums), in meters. A separate
function will compute this vector from latitude, longitude, and height.

Rotation matrix c2t_mat is the celestial-to-terrestrial matrix that is output
from the SOFA C2T family of routines. It represents the difference between the
GCRS and the ITRF.

The output position and velocity are in AU and AU/day. */


EXTERN DLLFUNC double STDCALL
psh_csta_press1(
	const double h,		/* meters above sea level */
	const double p);	/* altimeter setting */
/* Given altimeter setting, returns station pressure at h meters above sea
level, for h in the troposphere (less than 11000 m). Any units of measure may be
used for p; the returned value uses the same units.

In my experience, the "barometric pressure" reported by news media and even the
National Weather Service (unless you know where to look) is really altimeter
setting.

To convert inches of mercury to millibars, use the relationship 1013.25 mb =
760 mm Hg. Therefore, 1 inch Hg = 33.86 mb.

The output is always correct at the point where the altimeter setting is
measured. It's also correct above or below that point if the intervening
atmosphere is at standard temperature (15 C at sea level, decreasing 6.5 C
per km). */


EXTERN DLLFUNC double STDCALL
psh_csta_press2(
	const double h,		/* station's height above sea level, meters */
	const double p,		/* altimeter setting */
	const double t,		/* temperature at station, C */
	const double hm);	/* height where p is measured, meters */
/* More accurate version of psh_csta_press1(), for instances when the observing
station is at significantly different altitude from the point where altimeter
setting is measured. Corrects for temperature offset from standard, but assumes
standard 6.5 C per km lapse rate. */


EXTERN DLLFUNC double STDCALL
psh_csta_press3(
	const double h,		/* station's height above sea level, meters */
	const double p,		/* altimeter setting */
	const double t,		/* temperature at station, C */
	const double hm,	/* height where p is measured, meters */
	const double tm);	/* temperature where p is measured, C */
/* Still more accurate. Corrects for nonstandard lapse rate too. Assumes linear
temperature change between h and hm. */


EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_cref_u(
	double alt,		/* UNREFRACTED altitude, radians */
	const double p,		/* air pressure, millibars */
	const double t,		/* temperature at observer, C */
	const double acc);	/* desired accuracy, radians */
/* Returns refraction in radians at the given unrefracted altitude. I.e., add
the returned value to unrefracted altitude to get apparent altitude.

Argument p is "station pressure": the actual air pressure (not corrected to sea
level) at the observer.

Uses the refraction formulas from the Astronomical Almanac (p. B78 in the 2006
edition). It switches between the high and low altitude formulas at 20 degrees;
both formulas agree at that point. An accuracy argument is necessary because the
solution is iterative below 20 degrees. Outside that zone it is ignored. */


EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_cref_r(
	double alt,		/* REFRACTED altitude, radians */
	const double p,		/* air pressure, millibars */
	const double t,		/* temperature at observer, C */
	const double acc);	/* desired accuracy, radians */
/* Like the preceding function, except that the input altitude is the refracted
value. Subtract the returned value from refracted altitude to obtain unrefracted
altitude.

Above 20 degrees the solution is iterative, using argument "acc" to decide when
to stop. */


EXTERN DLLFUNC void STDCALL
psh_cstar_pmg(
	const double t0[2],	/* initial epoch, Julian Date */
	const double t1[2],	/* final epoch, Julian Date */
	double a0,		/* RA at t0, radians */
	double d0,		/* dec. at t0, radians */
	double mu_a,		/* proper motion in RA at t0, radians/year */
	double mu_d,		/* proper motion in dec at t0, radians/year */
	double pi,		/* parallax at t0, radians */
	double v,		/* radial velocity at t0, km/s */
	IAU_PVEC out_vec);	/* output: vector to star */
/* Applies proper motion and radial velocity to update star coordinates from t0
to t1. Output is rectangular coordinates in the same frame as the initial
coordinates. Uses the rigorous procedure on p. B64 in the Astronomical Almanac
2006. Annual parallax is not applied; there is a separate function for that.

Times are supplied as 2-part Julian Dates in the SOFA manner. I.e., initial
epoch = t0[0]+t0[1]. The high resolution of this format is unnecessary for star
positions, but it's used for consistency with other DLL functions.

Proper motions are radians per Julian year (365.25 days exactly). RA proper
motion is the great circle rate. That is, the unit has constant size regardless
of declination. This is the same convention used by the Hipparcos catalog. For
RA proper motion measured in coordinate rate, a different function is provided.

Radial velocity (positive if receding) is significant only for a few nearby
stars with rapid proper motion and high space velocity. Normally, it may be set
to zero.

In this function, parallax is used only to scale the radial velocity. If either
radial velocity or parallax is zero, the other argument has no effect. Zero or
negative parallax is allowed; the latter can occur due to measurement error, and
has the same effect as reversing the sign of radial velocity.

Output vector is very close to a unit vector, differing by space motion
(relative to distance at t0[]) over the time period. This makes it easier to
apply a formally correct annual parallax in one of my other functions. */


EXTERN DLLFUNC void STDCALL
psh_cstar_pmc(
	const double t0[2],	/* initial epoch, Julian Date */
	const double t1[2],	/* final epoch, Julian Date */
	double a0,		/* RA at t0, radians */
	double d0,		/* dec. at t0, radians */
	double mu_a,		/* proper motion in RA at t0, radians/year */
	double mu_d,		/* proper motion in dec at t0, radians/year */
	double pi,		/* parallax at t0, radians */
	double v,		/* radial velocity at t0, km/s */
	IAU_PVEC out_vec);	/* output: vector to star */
/* Identical to preceding function, except that proper motion in RA must be the
coordinate rate, rather than the great circle rate. I.e., the RA proper motions
used here will tend to be larger for stars near the poles. */


EXTERN DLLFUNC int STDCALL
psh_cstar_ap(
	const double jd[2],	/* pointer to 2-part Julian Date */
	constIAU_PVEC in_vec,	/* ICRS vector to star */
	const double pi,	/* parallax (radians) */
	const int place,	/* 1=astrometric, 2=apparent */
	constIAU_PVEC obs_v,	/* observer GCRS velocity */
	IAU_PVEC out_vec);	/* output GCRS vector */
/* Converts barycentric star position in the ICRS to its apparent geocentric or
topocentric position in the GCRS.

Time is supplied as 2-part Julian Date (jd[0] and jd[1]) in the SOFA manner.
The high resolution of this format is unnecessary for star positions, but it's
used for consistency with other DLL functions. Formally, the time scale is TT or
TDB. In practice, UTC could probably be used without noticeable effect.

Proper motion must be applied before calling this function. For strictly correct
results, in_vec must be the sum of a unit vector to the star, and the small
vector representing the proper motion correction. The psh_cstar_pmg() and
psh_cstar_pmc() functions generate such output.

Zero or negative parallax is allowed. Measurement error can cause the latter
when parallax is very small. It will causes the parallax correction to be
applied with reversed sign.

Argument "place": 0 = geometric place (parallax only), 1 = astrometric place
(same as 0), 2 = apparent place (parallax, light deflection, and aberration are
applied), 3 = same as 2, except no relativistic deflection of light by the Sun's
gravity. That last option is provided mainly to satisfy your curiousity.

Argument obs_v is the observer's velocity (AU/day) in the GCRS. A function is
provided elsewhere to compute this. To generate geocentric coordinates, or if
the dinural aberration (.3 arc seconds maximum) is negligible for your purposes,
set obs_v to {0.,0.,0.}.

The algorithm is the rigorous procedure on p. B64 in the 2006 Astronomical
Almanac, except that the observer's GCRS velocity is added to Earth's velocity,
causing dinural aberration to automatically appear in the result.

Geocentric parallax is not applied. For an Earthbound observer it is negligible,
even at the accuracy of the Hipparcos catalog.

The same variable may be used for in_vec and out_vec. That is, you can overwrite
the input with the output.

Parallax, aberration, and light deflection computations will use the JPL DE/LE
ephemeris if a file covering the desired time is open. Otherwise, the SOFA
ephemeris routine EPV00 is used. That routine is accurate. In a Monte Carlo test
spanning 1900 to 2100, error was .005 arc second root mean square, relative to
DE406. Even with much lower accuracy there should be no cause for concern, since
a star apparent position is relatively insensitive to Earth's position and
velocity.

Returned value is 0 if ok. That should always be the case now, unless the
function has a bug. (It used to depend on the JPL ephemeris, and nonzero
meant that you forgot to open an ephemeris file.) */


EXTERN DLLFUNC void STDCALL
psh_ctrs2hor(
	constIAU_PVEC in_vec,		/* ITRS coordinates */
	const double lon,		/* radians; longitude first */
	const double lat,
	const double xi,		/* radians; deflection of vertical */
	const double eta,		/* radians; deflection of vertical */
	const int reverse,		/* flag to reverse the conversion */
	IAU_PVEC out_vec);		/* may be same address as in_vec */
/* Rotate a vector's coordinate frame to the observer's horizontal frame, where
x, y, and z point to east, north, and zenith, respectively.

Deflection of the vertical can be applied; xi and eta are positive when the
plumb line intersects the celestial sphere north and east (respectively) of the
ellipsoidal normal.

The frame of in_vec must be oriented parallel to the ITRS (the Earth-centered
Earth-fixed frame that has +x directed to the longitude origin, +y to longitude
90 east, and +z to the north pole).

To reverse the conversion (i.e., in_vec is in the observer's horizontal frame),
pass a nonzero value in argument "reverse".

It is permissible to pass the same argument as in_vec and out_vec; in that case
output overwrites the input.

If there will be repeated calls with the same observer coordinates, consider
using psh_ctrs2hor_m(). */


EXTERN DLLFUNC void STDCALL
psh_ctrs2hor_m(
	const double lon,		/* radians; longitude first */
	const double lat,		/* latitude */
	const double xi,		/* radians; deflection of vertical */
	const double eta,		/* radians; deflection of vertical */
	IAU_RMAT rot_mat);		/* rotation matrix */
/* Sibling to the preceding function. Output is the matrix that rotates the
coordinate frame to the observer's horizontal frame. Since that matrix is
constant if the observer remains in the same place, it can be computed once and
passed to iau_rxp__() or iau__trxp() for repeated rotations. */


EXTERN DLLFUNC void STDCALL
psh_cvec2sph(
	constIAU_PVEC vec,	/* xyz coordinates */
	double *lon,		/* 0 to 2pi */
	double *lat,		/* -pi/2 to +pi/2 */
	double *r);		/* radius */
/* Convert rectangular to spherical coordinates. Identical to the SOFA P2S
routine, except P2S returns longitude in range -pi to +pi. */


EXTERN DLLFUNC void STDCALL
psh_cvec2azel(
	constIAU_PVEC vec,	/* x = east, y = north */
	double *az,		/* 0 to 2pi */
	double *el,		/* -pi/2 to +pi/2 */
	double *r);		/* radius */
/* Convert rectangular coordinates to azimuth and elevation, such that azimuth
has the conventional sense: zero at north, increasing east. (If you use a
standard rectangular to spherical conversion, azimuth is zero at east and
increases north.) */


EXTERN DLLFUNC void STDCALL
psh_cazel2vec(
	const double az,	/* radians; north = 0, east = pi/2 */
	const double el,	/* -pi/2 to +pi/2 */
	const double r,		/* radius */
	IAU_PVEC vec);		/* output vector */
/* Convert azimuth and elevation to rectangular. Inverse of psh_cvec2azel(). */


EXTERN DLLFUNC void STDCALL
psh_cpasep2vec(
	constIAU_PVEC vec_in,	/* not necessarily a unit vector */
	const double pa,	/* position angle, radians */
	const double sep,	/* separation angle, radians */
	IAU_PVEC vec_out);	/* unit vector */
/* Computes the point at pa and sep, relative to vec_in. This is the inverse
of computing the position angle and separation angle of two vectors. The
orientation of sep is such that pa = 0 will put vec_out north of vec_in, and
pa = pi/2 will put it east. */


/*######################## JPL Ephemeris Functions ###########################*/

EXTERN DLLFUNC int STDCALL
psh_jopen(
	const char *filename);
/* Opens the named JPL ephemeris, returns 0 if ok. This is a front end for
jpl_opneph(), simpler and less prone to accidental misuse. */


EXTERN DLLFUNC int STDCALL
psh_jspan(
	double *start,
	double *stop);
/* Return start and stop Julian Date of the time span covered by the open JPL
ephemeris. These JDs are NOT 2-part values in the SOFA manner; each JD is
in one variable. Returns 0 if ok, nonzero if no ephemeris file is open. */


EXTERN DLLFUNC double STDCALL
psh_jconst(
	const char *cname,	/* name of constant, CASE SENSITIVE */
	int *status);		/* 0 = ok (constant found) */
/* Return value of the named constant from the open JPL ephemeris. */


EXTERN DLLFUNC int STDCALL
psh_jplanet(
	const double t[2],	/* time (TT) */
	const int target,	/* planet or Moon or Sun */
	const int place,	/* 0=geometric, 1=astrometric, 2=apparent */
	const double accuracy,	/* desired accuracy, radians */
	constIAU_PVVEC obs_pv,	/* observer GCRS position and velocity */
	IAU_PVEC out_vec);	/* output GCRS vector */
/* Computes the Geocentric Celestial Reference System rectangular coordinates of
the Moon, the Sun, or a planet, with respect to the geocenter or topocenter.

Input time scale is Terrestrial Time in a 2-part SOFA-style Julian Date
(t[0] and t[1]). For maximum resolution these parts should be apportioned in
the SOFA "date/time" format, i.e., 0 <= t[1] < 1.

The "target" argument is the same code for the body that the JPL ephemeris
routines use. Enumerations for these are defined in sofajpl.h. (Earth and the
barycenters are not legal bodies for this function.)

Argument "place": 0 = geometric place, 1 = astrometric place (includes light
time), 2 = apparent place (light time, aberration, relativistic light
deflection by the Sun's gravity), 3 = same as 2, except NO light deflection.
That last option is provided mainly to satisfy your curiousity.

If geocentric geometric place is wanted, consider calling either jpl_pleph__()
or jpl_dpleph__(); that may be less bother than using this function.

The "accuracy" argument is necessary because the correction for light time is
computed iteratively until the solution converges to the requested accuracy.
Light time is not applicable to geometric place; in that case "accuracy" is
ignored.

Argument obs_pv is the observer's GCRS position (AU) and velocity (AU/day) in
rectangular form. To generate geocentric coordinates, use
{{0.,0.,0.},{0.,0.,0.}}.

In this function, except for the Sun or Moon, the returned point is the
center of mass of the planet and its satellites. "Therefore, the positions, when
converted to geocentric apparent places -- angular coordinates as seen from
Earth -- do not precisely indicate the center of the apparent planetary disk.
Displacements can amount to a few tens of milliarcseconds for Jupiter and
Saturn, a few milliarcseconds for Uranus and Neptune, and about 0.1 arcsecond
[sic] for Pluto." (USNO Circular 179)

The algorithm is the rigorous procedure on p. B60 (B64 for the Sun) in the 2006
Astronomical Almanac, except that the observer's position and velocity are added
to the values for the geocenter. In this way, geocentric parallax and dinural
aberration appear in the result without any extra effort.

To obtain full accuracy, a JPL DE/LE ephemeris file covering the desired time
must be open. A function is provided to do that. Any of the JPL DE/LE
ephemerides is acceptable; the same routine can handle them all.

If a JPL ephemeris is not available, the function automatically uses the self
contained SOFA ephemeris. In that case, Pluto and the Moon are not allowed.
Accuracy depends on the body. In a Monte Carlo test spanning 1900 to 2100, the
geocentric vector to the Sun had .005 arc second root mean square error vs.
DE406. Planetary coordinates were less accurate; I got these results (arc
seconds RMS): Mercury = .4, Venus = 1.5, Mars = 7.4, Jupiter = 20.3,
Saturn = 29.9.

The modulus of the output vector is geometric distance at time t, in AU.

Function return value is 0 if ok. There is no indication if the JPL ephemeris
fails; the function quietly switches to the SOFA ephemeris. (This can be avoided
by ensuring a JPL file is open and that it covers the time period.) A nonzero
value means the SOFA routine has failed to converge or "target" is illegal.
PLAN94 warns if date is outside 1000-3000 AD, but that's ignored here. */


/*################# Mathematical Conversions & Constants #####################*/


EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mdeg2rad(double degrees);
/* degrees to radians */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mrad2deg(double radians);
/* radians to degrees */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mmina2rad(double minutes);
/* minutes of arc to radians */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mrad2mina(double radians);
/* radians to minutes of arc */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mseca2rad(double seconds);
/* seconds of arc to radians */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mrad2seca(double radians);
/* radians to seconds of arc */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mhou2rad(double hours);
/* hours to radians */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mrad2hou(double radians);
/* radians to hours */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mmint2rad(double minutes);
/* minutes of time to radians */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mrad2mint(double radians);
/* radians to minutes of time */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_msect2rad(double seconds);
/* seconds of time to radians */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mrad2sect(double radians);
/* radians to seconds of time */

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
psh_mpi();
/* return pi */


/*############################# Time Functions ###############################*/

EXTERN DLLFUNC int STDCALL
psh_tdate2jd(
	const int y, const int m, const int d,
	const double fracday,	/* time of day, as a fractional day */
	double jd[2]);		/* returned 2-part Julian Date */
/* Convert Gregorian calendar date and time to a SOFA-style 2-part Julian Date.
Basically, this is just a shell to provide a more convenient interface to the
SOFA function iau_cal2jd__().

On return, jd[0] is midnight (ends in .5) and 0 <= jd[1] < 1. This is true
even if fracday is outside the [0, 1) interval. Output is valid if function
returns 0. */


EXTERN DLLFUNC void STDCALL
psh_tjdcalf(
	const int ndp,		/* no. of decimal places in ymdf[3] */
	const double dj1,	/* Julian date part 1 */
	const double dj2,	/* Julian date part 2 */
	int iymdf[4],		/* year, month, day, fractional day */
	int *j);		/* status */
/* Debugged version of SOFA routine IAU_JDCALF. The original function fails to
initialize *j to 0, so a garbage value is returned there when there is no error.
This fixes the bug. */


EXTERN DLLFUNC void STDCALL
psh_tjcal2jd(
	int y,		/* -1 = 2 BC, 0 = 1 BC, 1 = 1 AD, etc. */
	int m,		/* 1 = Jan, 2 = Feb, etc. */
	int d,
	double *jd0, double *jd1,
	int *error);
/* Like the SOFA cal2jd subroutine, but converts dates in the JULIAN calendar
to Julian Date. The result is returned in 2 parts: *jd0 == 2400000.5 (always)
and *jd1 == Modified Julian Date at 0h UT on the input date.

Unlike the SOFA routines, inputs are passed by value rather than by pointer.

Status is returned in *error:
	 0 == ok
	-1 == year rejected to avoid overflow
	-2 == no such month
	-3 == d is outside the normal range
If -1 or -2, the function returns with *jd0 and *jd1 unchanged.
If -3, there is not necessarily a true error. For example, in astromomical usage
January -4 or December 37 may be encountered. So, a -3 return value is simply an
advisory. The JD is computed correctly with any d, as long as the equivalent
date in normal format is within the allowed range.

The algorithm itself has unlimited range. However, to avoid overflowing a
32-bit integer, the year must be limited to -5879488 < y < 5879489.
*/


EXTERN DLLFUNC void STDCALL
psh_tjd2jcal(
	double jd0, double jd1,
	int *yea,		/* -1 = 2 BC, 0 = 1 BC, 1 = 1 AD, etc. */
	int *mon,		/* 1 = Jan, 2 = Feb, etc. */
	int *day,
	double *fracday,	/* fractional days */
	int *status);		/* 0 = ok; -1 = JD bad; -2 = internal error */
/* Convert a SOFA-style 2-part Julian Date to date and fractional days in the
JULIAN calendar. Any date in years -5879487 to 5879488 (inclusive) is OK.

This is a Julian calendar counterpart to SOFA routine jd2cal. Note that unlike
SOFA, the JD inputs are passed by value instead of by pointer, since that's
the natural way to do it in C/C++. */


/* A 2-part Julian date can apportion its components in various ways. To use the
SOFA terminology:
	2453996.674672     0.           JD method
	2451545.        2451.674672     J2000 method
	2400000.5      53996.174672     MJD method
	2453996.5           .174672     date & time method
The following functions will normalize (re-apportion) Julian date components to
the last three formats with minimum loss of accuracy. The first component of
the result will be an exact value. */

EXTERN DLLFUNC void STDCALL
psh_tnorm_datetime(double jd[2]);
/* Re-apportion the components of a 2-part Julian Date to the SOFA date/time
format: jd[0] is at midnight (ends in .5) and 0 <= jd[1] < 1 and the sum of the
components is unchanged. */

EXTERN DLLFUNC void STDCALL
psh_tnorm_j2000(double jd[2]);
/* Normalize a 2-part Julian Date so jd[0] = 2451545 (J2000 base epoch) and the
sum of the components is unchanged. */

EXTERN DLLFUNC void STDCALL
psh_tnorm_mjd(double jd[2]);
/* Normalize a 2-part Julian Date so jd[0] = 2400000.5 and the sum of the
components is unchanged. I.e., jd[1] = Modified Julian Date. */


EXTERN DLLFUNC double STDCALL
psh_ttdb_tt(
	const double jd[2]);
/* Return TDB-TT timescale offset in days. The topocentric component is
ignored. Julian Date is input as a 2-element array of doubles whose sum is the
JD. Formally, this JD should be TDB, but if it's TT there's no practical
difference. */


/*
UTC FUNCTIONS:

The following UTC functions use a table independent of the SOFA function
iau_dat(). Both tables are immediately available at DLL startup, but the SOFA
table is hard-coded and cannot be altered except by re-compiling. The other
table can be modified by the user.

Both tables are updated to the latest leap second announcement before I release
the DLL.
*/


EXTERN DLLFUNC int STDCALL
psh_tutc_read(
	const char *filename);
/* Replace the current UTC to TAI table with a user-provided table. The file
must be a text file in the format of the U.S. Naval Observatory leap second
table at ftp://maia.usno.navy.mil/ser7/tai-utc.dat
(The DLL comes pre-loaded with the latest version.)

The data structure automatically expands as needed, so for practical purposes
there's no upper bound to the table size.

The returned value is 0 if ok, 1 if the file cannot be opened, 2 if no valid
data are found. */


EXTERN DLLFUNC int STDCALL
psh_tutc_write(
	const char *filename);
/* Write the UTC to TAI table to a text file in USNO format. Returns 0 if ok,
nonzero if file write error. */


EXTERN DLLFUNC double STDCALL
psh_tutc_begin();
/* Returns Julian Date of the first entry in the UTC to TAI table (0 if the
table is empty. */


EXTERN DLLFUNC double STDCALL
psh_tutc_end();
/* Returns Julian Date of the last entry in the UTC to TAI table (0 if the
table is empty. */


EXTERN DLLFUNC int STDCALL
psh_tutc_size();
/* Returns number of entries in the leap second table. */


EXTERN DLLFUNC int STDCALL
psh_tutc2tt(	
	const int y, const int m, const int d,
	double fracday,		/* time of day, as a fractional day */
	double jd[2]);		/* returned 2-part Julian Date */
/* Convert Gregorian calendar date and time (UTC) to a SOFA-style 2-part Julian
Date in the Terrestrial Time scale.

On return, jd[0] is midnight (ends in .5) and 0 <= jd[1] < 1. This is true
even if fracday >= 1 (which occurs during a positive leap second, and yes, the
algorithm does work correctly in such a situation).

To be more precise, the applicable entry in the leap second table is determined
*only* by y, m, d. Once the correct entry is found, fracday is taken into
account. In the "rubber second" era before 1972, the TAI-UTC offset changes
continuously, so it's affected slightly by fracday. In 1972 and later, fracday
has no effect on the TAI-UTC offset; it simply adds directly to the output JD.

Result is valid if function returns 0. */


EXTERN DLLFUNC double STDCALL
psh_tutc_tt(
	const double tt[2],	/* 2-part Julian Date (TT) */
	double *step,		/* seconds; 0. except during leap second */
	int *status);
/* Returns the UTC-TT offset in *days* (always a negative number) at a given
Julian Date (TT), except that during a positive step adjustment, the returned
value is what the offset will be after the step adjustment ends.

Example output at a positive leap second (for clarity, the returned value is
expressed as seconds, and TT as date and time, to the nearest tenth second):

                   TT   returned value   step
2006 01 01 00:01:03.7            -64.2    0.0  (23:59:59.5 UTC)
                :04.2            -65.2    1.0  (23:59:60.0 UTC)
                :04.7            -65.2    1.0  (23:59:60.5 UTC)
                :05.2            -65.2    0.0  (00:00:00.0 UTC)

If the leap second had been negative:

2006 01 01 00:01:02.7            -64.2    0.0  (23:59:58.5 UTC)
                :03.2            -63.2    0.0  (00:00:00.0 UTC)
                :03.7            -63.2    0.0  (00:00:00.5 UTC)

In both cases, TT plus the returned value is a quasi-UTC, which can be
converted to sexagesimal components by routines which know nothing about leap
seconds. Then, adding "step" to the seconds component produces UTC time. Or
"step" can be ignored if a backward clock jump is acceptable. Another
alternative is to use "step" as a flag to freeze the seconds count during a
leap second.

The algorithm does not limit the value returned in "step" to 0 or 1, though in
practice that will be the case in 1972 and later. Fractional second steps will
occur in the 1960s.

The value in status:
-1 = definite bad date (predates first entry in leap second table)
 1 = possible bad date (after last table entry)
 0 = OK

A date before the beginning of the conversion table is a fatal error, and the
largest possible value is returned from the function to make the error obvious.

On the other hand, a date after last table entry may be OK. Depending on when
the next leap second occurs, the table may be still be valid. (The table with
the 1998 leap second as its last entry would have been valid until 2005.)
Therefore, returning 1 in *status is merely an advisory. To some extent it can
be suppressed by keeping the table updated to the latest BIH announcements. (To
show that a leap second will NOT be added before some future date, create a
table with that date as its last line, but the same TAI-UTC offset as the
previous line.)

Beware of the midnight pitfall. Due to the inexactness of floating point math,
an input time which is nominally 00:00:00 may be represented as 23:59:59.999...
If this happens at the changeover point between adjacent leap second table
entries, the wrong UTC offset will be returned. You can prevent that by adding
a small bias to TT before calling this function. */


EXTERN DLLFUNC double STDCALL
psh_ttai_utc(
	const int y, const int m, const int d,		/* UTC date */
	const double fd,	/* fractional days (UTC) */
	int *status);		/* 0 if ok */
/* Returns delta AT (== TAI-UTC) in seconds at a given epoch of UTC. Similar to
the SOFA DAT routine. Implements the conversion table available from the U.S.
Naval Observatory at ftp://maia.usno.navy.mil/ser7/tai-utc.dat

The value in *status:
	 1 = possible bad date (after end of table)
	 0 = OK
	-1 = definite bad date (before beginning of table)
	-2 = illegal date (no such month, etc.)

Argument d may exceed the normal range, e.g., Dec 32 is ok. When this occurs,
the computation uses the equivalent "normal" date.

A date before the beginning of the conversion table is a fatal error, and the
largest possible value is returned from the function to make the error obvious.
The same thing occurs if month is not 1 <= m <= 12.

On the other hand, a date after the end of the table is not necessarily bad.
Depending on when the next leap second occurs, the table may be still be valid.
(After the 1998 leap second, there were none until 2005.) Therefore, returning
1 in *status is merely an advisory.

Unlike the SOFA routine, fd (fractional day) is not restricted to 0 <= fd < 1.
The formula in the USNO text file is applied literally regardless. If fd is
outside the normal range, the date used to enter the delta AT lookup table does
not change; that date is controlled by (y, m, d) only.

For dates in the leap second era (1972 and later), fd is ignored. */
