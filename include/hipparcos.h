/*
hipparcos.h
Hipparcos routines for my sofajpl DLL.


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


See the Hipparcos web site:
http://www.rssd.esa.int/Hipparcos/



COMMA SEPARATED VALUE FORMAT

Some of these functions read or write text files in a CSV (comma separated
value) format compatible with Microsoft Excel and other spreadsheet and database
software. In this format, each line is one record (e.g., one star in the
Hipparcos catalog). Commas separate the fields of a record. Surround a field
with quotation marks if it contains a comma, space, or quotation mark. A field
with none of those special characters requires no quotes.

Inside a quoted field, double quotation marks ("") convert to a single quotation
mark. A single quotation mark terminates the quoted field. A comma or space
is considered part of the data.

Outside a quoted field, spaces are ignored and commas are field separators.

Empty fields are allowed. That is, a comma may be immediately followed by
another comma.


2006-07-07  File ready for release.
2006-07-20  Fixed errors in comments for psh_hstar().
*/



EXTERN DLLFUNC int STDCALL
psh_hcat_read(
	const char *filename);
/* Loads the star catalog with data from a text file in comma-separated value
format (see description at beginning of this file). Each line must contain
the same 78 fields as a Hipparcos machine readable catalog entry, although some
of these fields may be empty. (Every Hipparcos catalog entry I've seen contains
empty fields.)

The stars in the input file do not need to be sorted in any particular order.

The catalog is not cleared before loading a file, so you can load any number of
files in succession. During loading, stars are indexed by their Hipparcos
catalog (HIP) identifier (field 1 of each entry). If a HIP is already in the
catalog, the new data replaces the old.

Storage is acquired from the operating system as needed, so there is no specific
maximum number of stars.

Return value is 0 if ok. A nonzero value indicates incorrect file format or
memory exhausted.

NOTE: This function does not permanently affect the star catalog. At startup
it's initialized with a default load which has the 150 brightest stars and all
the navigational stars. */


EXTERN DLLFUNC void STDCALL
psh_hcat_clear();
/* Deletes all stars from the catalog. This only affects the catalog for the
duration of the run; the default catalog is always present at startup. */


EXTERN DLLFUNC int STDCALL
psh_hcat_size();
/* Returns number of stars in the catalog. */


EXTERN DLLFUNC int STDCALL
psh_hcat_write(
	const char *filename);
/* Outputs the star catalog as a comma-separated text file with 78 fields per
line. See the comments for the catalog read function. The stars will be sorted
in ascending Hipparcos catalog identifier order. The catalog itself is
unchanged. Return value is 0 if ok. */


EXTERN DLLFUNC char * STDCALL
psh_hfield(
	const int hip,		/* Hipparcos catalog number of star */
	const int field,	/* selects the data field to return (0 - 77) */
	char *str);	/* output destination; must be [13] or longer */
/* Returns one data field for a star. The data fields are numbered 0 to 77 in
the same sequence as the machine-readable Hipparcos catalog. The output
string will have 0 to 12 chars, followed by one null terminator. Return value
from function is str if ok, or 0 if error (no such star or field). */


EXTERN DLLFUNC int STDCALL
psh_hbary(
	const int hip,	/* Hipparcos catalog number of star */
	const double *jd,	/* 2-part Julian Date (TT) */
	IAU_PVEC output);
/* Returns the star's barycentric rectangular coordinates in the ICRS, adjusted
for proper motion to the epoch specified by jd. Returns 0 if ok, nonzero if
error (probably because the HIP number was not found in the catalog).

Julian Date is supplied in 2 parts (jd[0] and jd[1] in the SOFA manner.
Formally, the scale is Terrestrial Time, but the slowness of proper motion
makes this unimportant in practice.

The computational method uses the standard model of stellar motion described in
Hipparcos catalog Section 1.2.8 (online at the Hipparcos site), including the
radial velocities in Table 1.2.3. Specifically, Equation 1.2.16 is used, except
that normalization to a unit vector is omitted. That allows parallax to be
applied by simply adding the last term of Equation 1.2.18. */


EXTERN DLLFUNC int STDCALL
psh_hstar(
	const double *jd,	/* pointer to 2-part Julian date (TT) */
	const int hip,		/* Hipparcos catalog number of star */
	const int place,	/* 1=astrometric, 2=apparent */
	constIAU_PVEC v_obs,	/* observer velocity */
	IAU_PVEC out_vec);	/* output */
/* Computes the topocentric or geocentric coordinates of a star. The star is
specified by its Hipparcos catalog number. It must be in the internal star
catalog; see the functions that manipulate that catalog. The returned vector is
in the GCRS frame, and is adjusted for parallax, proper motion, and (if desired)
aberration.

Time is supplied as 2-part Julian Date (jd[0] and jd[1]) in the SOFA manner.
The high resolution of this format is unnecessary for star positions, but it's
used for consistency with other functions. Formally, the time scale is TT or
TDB. In practice, UTC could probably be used without noticeable effect.

Argument "place": 0 = geometric place (parallax only), 1 = astrometric place
(same as 0), 2 = apparent place (parallax, light deflection, and aberration are
applied), 3 = same as 2, except no relativistic deflection of light by the Sun's
gravity. That last option is provided mainly to satisfy your curiousity.

Argument v_obs is the observer's velocity (AU/day) in the GCRS. A function is
provided elsewhere to compute this. To generate geocentric coordinates, or if
the dinural aberration (.3 arc seconds maximum) is negligible for your purposes,
set v_obs to {0.,0.,0.}.

Geocentric parallax is not applied. For an Earthbound observer it is negligible,
even at the accuracy of the Hipparcos catalog.

Parallax, aberration, and light deflection computations will use the JPL DE/LE
ephemeris if a file covering the desired time is open. Otherwise, the SOFA
ephemeris routine EPV00 is used. That routine is accurate. In a Monte Carlo test
spanning 1900 to 2100, error was .005 arc second root mean square, relative to
DE406. Even with much lower accuracy there should be no cause for concern, since
a star apparent position is relatively insensitive to Earth's position and
velocity.

This function implements the standard model of stellar motion described in
Hipparcos catalog Section 1.2.8 (online at the Hipparcos site), including the
radial velocities in Table 1.2.3. Specifically, Equation 1.2.16 is used.

Parallax, aberration, and relativistic light deflection are computed with the
rigorous procedure on p. B64 in the 2006 Astronomical Almanac. The observer's
GCRS velocity is added to Earth's velocity, causing dinural aberration to
automatically appear in the result.

The returned value is 0 if ok, nonzero if the HIP number is not in the catalog
(in which case the contents of out_vec are undefined). */


EXTERN DLLFUNC int STDCALL
psh_hname_hip(
	const char *name);
/* Looks up the name (not case sensitive) in the star name table, returns the
Hipparcos catalog number (0 if not found). The name must NOT be enclosed in
quotes, even if it contains spaces. */


EXTERN DLLFUNC int STDCALL
psh_hname_read(
	const char *filename);
/* Loads a file into the star name table that cross references a star
designation (name) to its Hipparcos catalog number. Note that this table is
independent of the star catalog.

The input file is a text file with 2 comma separated values per line. (For a
description of the comma separated value format, see the beginning of this
file.) The first value on each line is a star designation of arbitrary length.

The second value is the Hipparcos catalog number. The same number may appear
multiple times in the input file. (On separate lines, NOT the same line!) For
example, there may be several designations for Sirius:

"NAME SIRIUS A",32349
"* ALF CMA A",32349
"* 9 CMA",32349

The different designations for the same star need not be together. Lines may be
arranged in any order in the file.

There is no definite upper limit on the number of table entries. The data
structure acquires storage from the system and expands dynamically as needed.

Data from the file augment the data already in the table. Any number of files
may be loaded. If a duplicate name occurs, new data overwrite the old.

Returns 0 if file loaded ok.

This function does not permanently affect the star name table. At startup the
table is automatically initialized with a default data set which contains the
proper name (sometimes more than one name or spelling), Bayer designation,
Flamsteed number, and navigation number, for the 150 brightest stars and all the
navigational stars. (Not all these designations are present for every star.) */


EXTERN DLLFUNC int STDCALL
psh_hname_write(
	const char *filename);
/* Write the star name table to a text file of comma separated values. The
format is compatible with the read function. Returns 0 if successful. */


EXTERN DLLFUNC int STDCALL
psh_hname_size();
/* returns number of names in the star name table */


EXTERN DLLFUNC void STDCALL
psh_hname_clear();
/* Deletes all data from the star name table. This only affects the table for
the duration of the run; the default table is always present at startup. The
empty table may be reloaded with the read function. */
