#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void ARMARECURSION1 (void *, void *, void *, void *, void *, void *, void *, void *);
extern void LGARCHSIM (void *, void *, void *, void *, void *, void *);
extern void VARMARECURSION1 (void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* describe routines */
static const R_CMethodDef CEntries[] = {
  {"ARMARECURSION1", (DL_FUNC) &ARMARECURSION1, 8},
  {"LGARCHSIM", (DL_FUNC) &LGARCHSIM, 6},
  {"VARMARECURSION1", (DL_FUNC) &VARMARECURSION1, 9},
  {NULL, NULL, 0}
};

/* register routines */
void R_init_lgarch(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
