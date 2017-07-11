#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(contract_mst_fortran)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(generate_neighbour_list)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gen_manycuts)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gen_progind_from_adjlst)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"contract_mst_fortran",    (DL_FUNC) &F77_NAME(contract_mst_fortran),     7},
    {"generate_neighbour_list", (DL_FUNC) &F77_NAME(generate_neighbour_list), 19},
    {"gen_manycuts",            (DL_FUNC) &F77_NAME(gen_manycuts),            10},
    {"gen_progind_from_adjlst", (DL_FUNC) &F77_NAME(gen_progind_from_adjlst), 13},
    {NULL, NULL, 0}
};

void R_init_CampaRi(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
