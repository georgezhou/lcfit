%module oblatenessfast
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}
%apply(double * INPLACE_ARRAY1, int DIM1){(double *phi, int np), (double *deficitFlux, int nf)};
%{
#include "oblateness.h"
%}

%include "oblateness.h"
