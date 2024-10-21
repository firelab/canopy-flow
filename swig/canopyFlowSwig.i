%module canopyFlowSwig  // Ensure the module name is specified correctly

%{
#include "canopyFlowSwig.h"  // Include the main header for your classes and functions
%}


// Include the necessary standard library for using std::vector
%include "std_vector.i" 
%include "std_string.i"

%template(DoubleVector) std::vector<double>;

%include "canopyFlowSwig.h"  