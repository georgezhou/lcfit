%module elliptic
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY1, int DIM1){(double *n, int ln), (double *k, int lk), (double *nk, int lnk),(double *ee, int lee),(double *kk, int lkk)};

%{
#include <iostream.h>
#include "elliptic.h"
%}

%include "elliptic.h"
