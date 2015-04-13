/*
sofajpl.h
C++ shells for the SOFA and JPL Fortran routines in my DLL.


Copyright 2006 by Paul S. Hirose. You may use and distribute this software
(as-is or modified), but not for profit. Unless you are the sole user, I must
receive credit in a place reasonably obvious to users (e.g., by clicking
Help, About).

Furthermore, you're bound by the SOFA Software License. A copy is at the end of
this file. It's also available at the IAU SOFA web site.

Furthermore, you should give credit to:
1) NASA Jet Propulsion Laboratory, for their Planetary and Lunar Ephemerides
2) ESA Hipparcos Space Astrometry Mission, for their Hipparcos Catalog
3) SIMBAD database, operated at CDS, Strasbourg, France.


It is necessary to read the SOFA and JPL documentation to understand these
functions. See the SOFA web site at
http://www.iau-sofa.rl.ac.uk/
and download README.txt from the JPL FTP site at
ftp://ssd.jpl.nasa.gov//pub/eph/export/


The SOFA copyright notice is at the end of this file.


2006-07-02  Created file.
2006-07-07  Ready for initial release.
*/


/* Typedefs for defining SOFA arrays */

typedef double IAU_PVEC[3];	/* position vector */
typedef double IAU_PVVEC[2][3];	/* position/velocity vector */
typedef double IAU_RMAT[3][3];	/* rotation matrix */
typedef double IAU_RVEC[3];	/* rotation vector */

/* Typedefs for declaring const SOFA array arguments in function declarations */

typedef const double constIAU_PVEC[3];		/* position vector */
typedef const double constIAU_PVVEC[2][3];	/* position/velocity vector */
typedef const double constIAU_RMAT[3][3];	/* rotation matrix */
typedef const double constIAU_RVEC[3];		/* rotation vector */

/* NOTE: Fortran and C have opposite conventions for storing matrices. See the
SOFA documentation. This can be a problem if you create a rotation matrix in
C and pass it to a SOFA function! You must transpose the matrix first (a SOFA
function can do that). */


/* JPL ephemeris codes for the planets, the solar system barycenter, and the
Earth-Moon barycenter. */

enum {jpl_mercury = 1, jpl_venus, jpl_earth, jpl_mars, jpl_jupiter, jpl_saturn,
jpl_uranus, jpl_neptune, jpl_pluto, jpl_moon, jpl_sun, jpl_ssbary, jpl_embary};


/*############################# SOFA Functions ###############################*/


EXTERN DLLFUNC void STDCALL
iau_a2af(const int ndp, const double angle, char *sign, int *idmsf);

EXTERN DLLFUNC void STDCALL
iau_a2tf(const int ndp, const double angle, char *sign, int *ihmsf);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_anp(const double a);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_anpm(const double a);

EXTERN DLLFUNC void STDCALL
iau_bi00(double *dpsibi, double *depsbi, double *dra);

EXTERN DLLFUNC void STDCALL
iau_bp00(const double date1, const double date2, IAU_RMAT rb, IAU_RMAT rp,
 IAU_RMAT rbp);

EXTERN DLLFUNC void STDCALL
iau_bpn2xy(constIAU_RMAT rbpn, double *x, double *y);

EXTERN DLLFUNC void STDCALL
iau_c2i00a(const double date1, const double date2, IAU_RMAT rc2i);

EXTERN DLLFUNC void STDCALL
iau_c2i00b(const double date1, const double date2, IAU_RMAT rc2i);

EXTERN DLLFUNC void STDCALL
iau_c2ibpn(const double date1, const double date2, constIAU_RMAT rbpn,
 IAU_RMAT rc2i);

EXTERN DLLFUNC void STDCALL
iau_c2ixy(const double date1, const double date2, const double x,
 const double y, IAU_RMAT rc2i);

EXTERN DLLFUNC void STDCALL
iau_c2ixys(const double x, const double y, const double s, IAU_RMAT rc2i);

EXTERN DLLFUNC void STDCALL
iau_c2s(constIAU_PVEC p, double *theta,
 double *phi);

EXTERN DLLFUNC void STDCALL
iau_c2t00a(const double tta, const double ttb, const double uta,
 const double utb, const double xp, const double yp, IAU_RMAT rc2t);

EXTERN DLLFUNC void STDCALL
iau_c2t00b(const double tta, const double ttb, const double uta,
 const double utb, const double xp, const double yp, IAU_RMAT rc2t);

EXTERN DLLFUNC void STDCALL
iau_c2tceo(constIAU_RMAT rc2i, const double era, constIAU_RMAT rpom,
 IAU_RMAT rc2t);

EXTERN DLLFUNC void STDCALL
iau_c2teqx(constIAU_RMAT rbpn, const double gst, constIAU_RMAT rpom,
 IAU_RMAT rc2t);

EXTERN DLLFUNC void STDCALL
iau_c2tpe(const double tta, const double ttb, const double uta,
 const double utb, const double dpsi, const double deps, const double xp,
 const double yp, IAU_RMAT rc2t);

EXTERN DLLFUNC void STDCALL
iau_c2txy(const double tta, const double ttb, const double uta,
 const double utb, const double x, const double y, const double xp,
 const double yp, IAU_RMAT rc2t);

EXTERN DLLFUNC void STDCALL
iau_cal2jd(const int iy, const int im, const int id, double *djm0, double *djm,
 int *j);

EXTERN DLLFUNC void STDCALL
iau_cp(constIAU_PVEC p, IAU_PVEC c);

EXTERN DLLFUNC void STDCALL
iau_cpv(constIAU_PVVEC pv, IAU_PVVEC c);

EXTERN DLLFUNC void STDCALL
iau_cr(constIAU_RMAT r, IAU_RMAT c);

EXTERN DLLFUNC void STDCALL
iau_d2tf(const int ndp, const double days, char *sign, int *ihmsf);

EXTERN DLLFUNC void STDCALL
iau_dat(const int iy, const int im, const int id, const double fd,
 double *deltat, int *j);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_dtdb(const double epoch1, const double epoch2, const double ut,
 const double elong, const double u, const double v);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_ee00(const double date1, const double date2, const double epsa,
 const double dpsi);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_ee00a(const double date1, const double date2);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_ee00b(const double date1, const double date2);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_eect00(const double date1, const double date2);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_epb(const double dj1, const double dj2);

EXTERN DLLFUNC void STDCALL
iau_epb2jd(const double epb, double *djm0, double *djm);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_epj(const double dj1, const double dj2);

EXTERN DLLFUNC void STDCALL
iau_epj2jd(const double epj, double *djm0, double *djm);

EXTERN DLLFUNC void STDCALL
iau_epv00(const double epoch1, const double epoch2, IAU_PVVEC pvh,
 IAU_PVVEC pvb, int *jstat);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_eqeq94(const double epoch1, const double epoch2);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_era00(const double dj1, const double dj2);

EXTERN DLLFUNC void STDCALL
iau_fk52h(const double r5, const double d5, const double dr5, const double dd5,
 const double px5, const double rv5, double *rh, double *dh, double *drh,
 double *ddh, double *pxh, double *rvh);

EXTERN DLLFUNC void STDCALL
iau_fk5hip(IAU_RMAT r5h, IAU_RVEC s5h);

EXTERN DLLFUNC void STDCALL
iau_fk5hz(const double r5, const double d5, const double epoch1,
 const double epoch2, double *rh, double *dh);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_gmst00(const double uta, const double utb, const double tta,
 const double ttb);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_gmst82(const double dj1, const double dj2);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_gst00a(const double uta, const double utb, const double tta,
 const double ttb);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_gst00b(const double uta, const double utb);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_gst94(const double uta, const double utb);

EXTERN DLLFUNC void STDCALL
iau_h2fk5(const double rh, const double dh, const double drh, const double ddh,
 const double pxh, const double rvh, double *r5, double *d5, double *dr5,
 double *dd5, double *px5, double *rv5);

EXTERN DLLFUNC void STDCALL
iau_hfk5z(const double rh, const double dh, const double epoch1,
 const double epoch2, double *r5, double *d5, double *dr5, double *dd5);

EXTERN DLLFUNC void STDCALL
iau_ir(IAU_RMAT r);

EXTERN DLLFUNC void STDCALL
iau_jd2cal(const double dj1, const double dj2, int *iy, int *im, int *id,
 double *fd, int *j);


#if 0

/* CAUTION - the next function has a bug. Status variable *j is not initialized,
so it returns garbage when no error has occured. My function psh_jdcalf__() is
identical except for the bug. (Its prototype is in another file.) */

EXTERN void iau_jdcalf(const int ndp, const double dj1, const double dj2,
 int *iymdf, int *j);

#endif


EXTERN DLLFUNC void STDCALL
iau_num00a(const double date1, const double date2, IAU_RMAT rmatn);

EXTERN DLLFUNC void STDCALL
iau_num00b(const double date1, const double date2, IAU_RMAT rmatn);

EXTERN DLLFUNC void STDCALL
iau_numat(const double epsa, const double dpsi, const double deps,
 IAU_RMAT rmatn);

EXTERN DLLFUNC void STDCALL
iau_nut00a(const double date1, const double date2, double *dpsi, double *deps);

EXTERN DLLFUNC void STDCALL
iau_nut00b(const double date1, const double date2, double *dpsi, double *deps);

EXTERN DLLFUNC void STDCALL
iau_nut80(const double epoch1, const double epoch2, double *dpsi, double *deps);

EXTERN DLLFUNC void STDCALL
iau_nutm80(const double epoch1, const double epoch2, IAU_RMAT rmatn);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_obl80(const double date1, const double date2);

EXTERN DLLFUNC void STDCALL
iau_p2pv(constIAU_PVEC p, IAU_PVVEC pv);

EXTERN DLLFUNC void STDCALL
iau_p2s(constIAU_PVEC p, double *theta, double *phi, double *r);

EXTERN DLLFUNC void STDCALL
iau_pap(constIAU_PVEC a, constIAU_PVEC b, double *theta);

EXTERN DLLFUNC void STDCALL
iau_pas(const double al, const double ap, const double bl, const double bp,
 double *theta);

EXTERN DLLFUNC void STDCALL
iau_pdp(constIAU_PVEC a, constIAU_PVEC b, double *adb);

EXTERN DLLFUNC void STDCALL
iau_plan94(const double epoch1, const double epoch2, const int np, IAU_PVVEC pv,
 int *j);

EXTERN DLLFUNC void STDCALL
iau_pm(constIAU_PVEC p, double *r);

EXTERN DLLFUNC void STDCALL
iau_pmat00(const double date1, const double date2, IAU_RMAT rbp);

EXTERN DLLFUNC void STDCALL
iau_pmat76(const double epoch1, const double epoch2, IAU_RMAT rmatp);

EXTERN DLLFUNC void STDCALL
iau_pmp(constIAU_PVEC a, constIAU_PVEC b, IAU_PVEC amb);

EXTERN DLLFUNC void STDCALL
iau_pn(constIAU_PVEC p, double *r, IAU_PVEC u);

EXTERN DLLFUNC void STDCALL
iau_pn00(const double date1, const double date2, const double dpsi,
 const double deps, double *epsa, IAU_RMAT rb, IAU_RMAT rp, IAU_RMAT rbp,
 IAU_RMAT rn, IAU_RMAT rbpn);

EXTERN DLLFUNC void STDCALL
iau_pn00a(const double date1, const double date2, double *dpsi, double *deps,
 double *epsa, IAU_RMAT rb, IAU_RMAT rp, IAU_RMAT rbp, IAU_RMAT rn,
 IAU_RMAT rbpn);

EXTERN DLLFUNC void STDCALL
iau_pn00b(const double date1, const double date2, double *dpsi, double *deps,
 double *epsa, IAU_RMAT rb, IAU_RMAT rp, IAU_RMAT rbp, IAU_RMAT rn,
 IAU_RMAT rbpn);

EXTERN DLLFUNC void STDCALL
iau_pnm00a(const double date1, const double date2, IAU_RMAT rbpn);

EXTERN DLLFUNC void STDCALL
iau_pnm00b(const double date1, const double date2, IAU_RMAT rbpn);

EXTERN DLLFUNC void STDCALL
iau_pnm80(const double epoch1, const double epoch2, IAU_RMAT rmatpn);

EXTERN DLLFUNC void STDCALL
iau_pom00(const double xp, const double yp, const double sp, IAU_RMAT rpom);

EXTERN DLLFUNC void STDCALL
iau_ppp(constIAU_PVEC a, constIAU_PVEC b, IAU_PVEC apb);

EXTERN DLLFUNC void STDCALL
iau_ppsp(constIAU_PVEC a, const double s, constIAU_PVEC b, IAU_PVEC apsb);

EXTERN DLLFUNC void STDCALL
iau_pr00(const double date1, const double date2, double *dpsipr,
 double *depspr);

EXTERN DLLFUNC void STDCALL
iau_prec76(const double ep01, const double ep02, const double ep11,
 const double ep12, double *zeta, double *z, double *theta);

EXTERN DLLFUNC void STDCALL
iau_pv2p(constIAU_PVVEC pv, IAU_PVEC p);

EXTERN DLLFUNC void STDCALL
iau_pv2s(constIAU_PVVEC pv, double *theta, double *phi, double *r, double *td,
 double *pd, double *rd);

EXTERN DLLFUNC void STDCALL
iau_pvdpv(constIAU_PVVEC a, constIAU_PVVEC b, double *adb);

EXTERN DLLFUNC void STDCALL
iau_pvm(constIAU_PVVEC pv, double *r, double *s);

EXTERN DLLFUNC void STDCALL
iau_pvmpv(constIAU_PVVEC a, constIAU_PVVEC b, IAU_PVVEC amb);

EXTERN DLLFUNC void STDCALL
iau_pvppv(constIAU_PVVEC a, constIAU_PVVEC b, IAU_PVVEC apb);

EXTERN DLLFUNC void STDCALL
iau_pvstar(constIAU_PVVEC pv, double *ra, double *dec, double *pmr, double *pmd,
 double *px, double *rv, int *j);

EXTERN DLLFUNC void STDCALL
iau_pvu(const double dt, constIAU_PVVEC pv, IAU_PVVEC upv);

EXTERN DLLFUNC void STDCALL
iau_pvup(const double dt, constIAU_PVVEC pv, IAU_PVEC p);

EXTERN DLLFUNC void STDCALL
iau_pvxpv(constIAU_PVVEC a, constIAU_PVVEC b, IAU_PVVEC axb);

EXTERN DLLFUNC void STDCALL
iau_pxp(constIAU_PVEC a, constIAU_PVEC b, IAU_PVEC axb);

EXTERN DLLFUNC void STDCALL
iau_rm2v(constIAU_RMAT r, IAU_RVEC w);

EXTERN DLLFUNC void STDCALL
iau_rv2m(constIAU_RVEC w, IAU_RMAT r);

EXTERN DLLFUNC void STDCALL
iau_rx(const double phi, IAU_RMAT r);

EXTERN DLLFUNC void STDCALL
iau_rxp(constIAU_RMAT r, constIAU_PVEC p, IAU_PVEC rp);

EXTERN DLLFUNC void STDCALL
iau_rxpv(constIAU_RMAT r, constIAU_PVVEC pv, IAU_PVVEC rpv);

EXTERN DLLFUNC void STDCALL
iau_rxr(constIAU_RMAT a, constIAU_RMAT b, IAU_RMAT atb);

EXTERN DLLFUNC void STDCALL
iau_ry(const double theta, IAU_RMAT r);

EXTERN DLLFUNC void STDCALL
iau_rz(const double psi, IAU_RMAT r);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_s00(const double date1, const double date2, const double x, const double y);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_s00a(const double date1, const double date2);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_s00b(const double date1, const double date2);

EXTERN DLLFUNC void STDCALL
iau_s2c(const double theta, const double phi, IAU_PVEC c);

EXTERN DLLFUNC void STDCALL
iau_s2p(const double theta, const double phi, const double r, IAU_PVEC p);

EXTERN DLLFUNC void STDCALL
iau_s2pv(const double theta, const double phi, const double r, const double td,
 const double pd, const double rd, IAU_PVVEC pv);

EXTERN DLLFUNC void STDCALL
iau_s2xpv(const double s1, const double s2, constIAU_PVVEC pv, IAU_PVVEC spv);

EXTERN DLLFUNC void STDCALL
iau_sepp(constIAU_PVEC a, constIAU_PVEC b, double *s);

EXTERN DLLFUNC void STDCALL
iau_seps(const double al, const double ap, const double bl, const double bp,
 double *s);

EXTERN DLLFUNC double STDCALL ATTRIBUTE((const))
iau_sp00(const double date1, const double date2);

EXTERN DLLFUNC void STDCALL
iau_starpm(const double ra1, const double dec1, const double pmr1,
 const double pmd1, const double px1, const double rv1, const double ep1a,
 const double ep1b, const double ep2a, const double ep2b, double *ra2,
 double *dec2, double *pmr2, double *pmd2, double *px2, double *rv2, int *j);

EXTERN DLLFUNC void STDCALL
iau_starpv(const double ra, const double dec, const double pmr,
 const double pmd, const double px, const double rv, IAU_PVVEC pv, int *j);

EXTERN DLLFUNC void STDCALL
iau_sxp(const double s, constIAU_PVEC p, IAU_PVEC sp);

EXTERN DLLFUNC void STDCALL
iau_sxpv(const double s, constIAU_PVVEC pv, IAU_PVVEC spv);

EXTERN DLLFUNC void STDCALL
iau_tr(constIAU_RMAT r, IAU_RMAT rt);

EXTERN DLLFUNC void STDCALL
iau_trxp(constIAU_RMAT r, constIAU_PVEC p, IAU_PVEC trp);

EXTERN DLLFUNC void STDCALL
iau_trxpv(constIAU_RMAT r, constIAU_PVVEC pv, IAU_PVVEC trpv);

EXTERN DLLFUNC void STDCALL
iau_xys00a(const double date1, const double date2, double *x, double *y,
 double *s);

EXTERN DLLFUNC void STDCALL
iau_xys00b(const double date1, const double date2, double *x, double *y,
 double *s);

EXTERN DLLFUNC void STDCALL
iau_zp(IAU_PVEC p);

EXTERN DLLFUNC void STDCALL
iau_zpv(IAU_PVVEC pv);

EXTERN DLLFUNC void STDCALL
iau_zr(IAU_RMAT r);


/*######################## JPL Ephemeris Functions ###########################*/


EXTERN DLLFUNC int STDCALL
jpl_pleph(
	const double et,		/* Julian date (TT) */
	const int ntarg,		/* target body */
	const int ncent,		/* center body */
	IAU_PVVEC rrd);			/* position & velocity */
/* Input time is Julian date (Terrestrial Time) in a double precision variable.

Output is a SOFA pv-vector, a [2][3] array of doubles. Position is AU, and
velocity AU/day.

Rewritten from the original JPL code so the ephemeris file is chosen at runtime,
and may be changed during the run as often as desired. Another difference is
that it runs silently and return a status code:
	0  ok
	1  date outside ephemeris range
	2  no nutations or librations in this ephemeris
	3  file read error
In other respects this is the JPL routine; see the comments in the Fortran
source code. */


EXTERN DLLFUNC int STDCALL
jpl_dpleph(
	const double et2z[2],
	const int ntarg,
	const int ncent,
	IAU_PVVEC rrd);
/* Like jpl_pleph(), except time is a 2-part SOFA-style Julian Date
(et2z[0] and et2z[1]). For maximum interpolation accuracy these parts should be
apportioned in the SOFA "date/time" format, i.e., 0 <= t[1] < 1. */


EXTERN DLLFUNC int STDCALL
jpl_const(
	char nam[],
	double val[],
	double sss[3],
	int *n);
/* Loads nam[], val[], and sss[3] with constants from the header of the open
ephemeris file. The number of names in nam[] and doubles in val[] (same number
for both) is returned in *n. Each name is 6 chars long, padded on the right with
spaces if necessary.

A buffer overrun will occur if nam[] or val[] is too short. (The Fortran code
allocates enough space for *n = 400.)

For the status code returned by the function, see jpl_pleph() comments. */


EXTERN DLLFUNC int STDCALL
jpl_a2eph(
	const char *head,	/* name of ASCII header file */
	const char *ephdat,	/* name of ephemeris ASCII data file */
	const char *ephout,	/* name of binary output file */
	char *errmsg,		/* error message */
	const int errmsg_len);	/* size of errmsg[] buffer */
/* Convert JPL DE/LE ephemeris files from ASCII distribution format to the
binary format in which they are used.

It is an error if the output file already exists.

Returns 0 if OK; if not, a null-terminated message will be in errmsg[]. If the
message is longer than errmsg[], it will be truncated to fit. Contents of
errmsg[] are undefined if the function returned 0.

This is a modification of JPL's ASC2EPH program; see their FTP site for source
code of the original routine, input files, and documentation. */


EXTERN DLLFUNC int STDCALL
jpl_bmerge(
	const char *ineph1,	/* name of input file #1 */
	const char *ineph2,	/* name of input file #2 */
	const char *outeph,	/* output file */
	char *errmsg,		/* error message */
	const int errmsg_len);	/* size of errmsg[] buffer */
/* Combines binary JPL ephemerides ineph1 and ineph2 into single ephemeris
outeph. The time spans of ineph1 and ineph2 must abut or overlap.

It is an error if the output file already exists.

Returns 0 if OK; if not, a null-terminated message will be in errmsg[]. If the
message is longer than errmsg[], it will be truncated to fit. Contents of
errmsg[] are undefined if the function returned 0.

This is a modification of JPL's BINMERGE program; see their FTP site for source
code of the original routine, and documentation. */


EXTERN DLLFUNC int STDCALL
jpl_bshort(
	const char *ineph,	/* name of input file */
	const double t1,	/* output file start time, 1-part JD */
	const double t2,	/* output file end time, 1-part JD */
	const char *outeph,	/* name of output file */
	char *errmsg,		/* error message */
	const int errmsg_len);	/* size of errmsg[] buffer */
/* Make a JPL DE/LE ephemeris binary file with a shorter time span than the
input file.

Note that t1 and t2 each point to a Julian date in stored in one double
precision variable, not a 2-part Julian date.

It is an error if the output file already exists.

Returns 0 if OK; if not, a null-terminated message will be in errmsg[]. If the
message is longer than errmsg[], it will be truncated to fit. Contents of
errmsg[] are undefined if the function returned 0.

This is a modification of JPL's BINSHORT program; see their FTP site for
documentation and source code of the original routine. */



/* SOFA copyright notice */
#if 0

copyr.lis                                                 2005 January 1


COPYRIGHT NOTICE


Text like the following appears at the end of every SOFA routine.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2005
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. The Software is made available free of charge for use by:
*
*     a) private individuals for non-profit research; and
*
*     b) non-profit educational, academic and research institutions.
*
*  3. Commercial use of the Software is specifically excluded from the
*     terms and conditions of this license.  Commercial use of the
*     Software is subject to the prior written agreement of the Board on
*     terms to be agreed.
*
*  4. The provision of any version of the Software under the terms and
*     conditions specified herein does not imply that future versions
*     will also be made available under the same terms and conditions.
*
*  5. The user may modify the Software for his/her own purposes.  The
*     user may distribute the modified software provided that the Board
*     is informed and that a copy of the modified software is made
*     available to the Board on request.  All modifications made by the
*     user shall be clearly identified to show how the modified software
*     differs from the original Software, and the name(s) of the
*     affected routine(s) shall be changed.  The original SOFA Software
*     License text must be present.
*
*  6. In any published work produced by the user and which includes
*     results achieved by using the Software, the user shall acknowledge
*     that the Software was used in producing the information contained
*     in such publication.
*
*  7. The user may incorporate or embed the Software into other software
*     products which he/she may then give away free of charge but not
*     sell provided the user makes due acknowledgement of the use which
*     he/she has made of the Software in creating such software
*     products.  Any redistribution of the Software in this way shall be
*     made under the same terms and conditions under which the user
*     received it from the SOFA Center.
*
*  8. The user shall not cause the Software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  9. The Software is provided to the user "as is" and the Board makes
*     no warranty as to its use or performance.   The Board does not and
*     cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Board makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*
*-----------------------------------------------------------------------

      END


/* end SOFA copyright notice */
#endif
