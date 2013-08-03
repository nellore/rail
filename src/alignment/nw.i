%module nw

%{
    #define SWIG_FILE_WITH_INIT
    #include "nw.h"
%}

%include "numpy.i"


%init %{
  import_array();
%}
/*type maps for input arrays and strings*/
%apply (int* INPLACE_ARRAY2, int DIM1, int DIM2) {(int* D, int Dx, int Dy)}
%apply (int* IN_ARRAY2, int DIM1, int DIM2) {(int* mat, int mx, int my)}
%apply (char *STRING, int LENGTH) {(char *xstr, int xL),(char *ystr, int yL)}
//%apply (char* IN_ARRAY, int DIM1){(char *xstr, int xL),(char *ystr, int yL)}

%include "nw.h"
