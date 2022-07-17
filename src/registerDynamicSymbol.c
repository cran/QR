// RegisteringDynamic Symbols

#include <R.h>
#include <R_ext/Complex.h>
#include <R_ext/Rdynload.h>

extern void qrc(int *MX, int *NX, double *X, int *LDX, double *tau, int outlwork);
extern void qrcomplexc(int *MX, int *NX, Rcomplex *X, int *LDX, Rcomplex *tau, int outlwork);

static const R_CMethodDef CEntries[] = {
  {"qrc",      ( DL_FUNC ) &qrc,         6},
  {"qrcomplexc",      ( DL_FUNC ) &qrcomplexc,         6},
  {NULL, NULL, 0}
};

void R_init_QR(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
