#include <Rcpp.h>


#include "mshrmassintensmap.h"

// template<> SEXP Rcpp::wrap(Rcpp::DataFrame(msms)
// template<> operator HrMassIntensMap ( Rcpp::DataFrame & );


Rcpp::DataFrame toDF(const HrMassIntensMap &);

HrMassIntensMap toMap(const Rcpp::DataFrame & );
