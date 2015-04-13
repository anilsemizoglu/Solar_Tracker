/*
dependencies.h
Paul S. Hirose

2006-06-30  Began work.
2006-07-07  Ready for initial release.
*/


/* Comment this out if not using GNU compiler. Should be the only change
needed in this file.
#define GNU
*/


#ifdef __cplusplus
#define EXTERN extern "C" 
#else
#define EXTERN extern
#endif


/* Declare imported DLL functions to help compiler generate better code, but
control exported DLL functions with a .def file, not the "dllexport" keyword. */

#ifdef DEFINE_DLL	/* true inside the DLL source code */
#define DLLFUNC
#else			/* DLL users see this */
#define DLLFUNC __declspec(dllimport)
#endif


#ifdef GNU	/* using GNU compiler */

#define ATTRIBUTE(a) __attribute__(a)
#define STDCALL __attribute__((stdcall))

#else		/* using Microsoft compiler */

#define ATTRIBUTE(a)
#define STDCALL __stdcall

#endif		/* end compiler dependencies */
