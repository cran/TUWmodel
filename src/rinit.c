#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(hbvmodel)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hbvmodel_dual)(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"hbvmodel",      (DL_FUNC) &F77_NAME(hbvmodel),      9},
    {"hbvmodel_dual", (DL_FUNC) &F77_NAME(hbvmodel_dual), 9},
    {NULL, NULL, 0}
};

void R_init_TUWmodel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
