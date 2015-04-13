// ephemeris headers for sofeJPL libraries
#include "dependencies.h"       /* must precede my other headers */
#include "sofajpl.h"            /* SOFA and JPL functions */
#include "sofajpl_sup1.h"       /* supplementary routines */
#include "hipparcos.h"          /* Hipparcos star catalog routines */

void psh_dsex_u( const int dmsf[], int skip, const int precision, const char *str);
void psh_dsex_s( const char sign, const int dmsf[], const int skip, const int precision, const char *str);
void psh_dangle_u( double angle,const int skip,const int precision,const char *str);
void psh_dangle_s( double angle, const int skip, const int precision, const char *str);
void psh_dtime_u( const double angle, const int skip, const int precision, const char *str);
void sun_equatorial(double &ra, double &dec, int y, int m, int d, int h, int min, int s);

