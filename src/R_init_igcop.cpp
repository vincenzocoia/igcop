#include <Rcpp.h>
using namespace Rcpp;

// As per Rcpp GitHub Issue #636
// https://github.com/RcppCore/Rcpp/issues/636#issuecomment-280985661
// Attempt to avoid the message upon a package Check:
//-----
//   Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’
//
// It is good practice to use registered native symbols and to disable
//   symbol search.
//-----
void R_init_igcop(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
