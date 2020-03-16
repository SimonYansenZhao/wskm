#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "wskm.h"

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

static const R_CMethodDef cMethods[] = {
    CALLDEF(ewkm, 15),
	CALLDEF(fgkm, 22),
	CALLDEF(twkm, 22),
    CALLDEF(parseGroup, 3),
    NULL
};


void R_init_wskm(DllInfo *dll) {
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
