%module pev_bindings

%{
#include "wrap.hpp"
%}

%include "std_string.i"   // <-- add this if wrap.hpp uses std::string
%include "std_vector.i"

/* tell SWIG how to convert std::vector<int> */
%include "std_vector.i"
namespace std {
  %template(IntVector) vector<int>;
}

/* wrap pileup-events */
%include "wrap.hpp"
