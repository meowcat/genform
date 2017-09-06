#include <Rcpp.h>
using namespace Rcpp;

#include "wrap.h"



Rcpp::DataFrame toDF(const HrMassIntensMap & MIM)
{
  int rows = MIM.size();
  NumericVector _mz(rows);
  NumericVector _int(rows);
  NumericVector::iterator i_mz = _mz.begin();
  NumericVector::iterator i_int = _int.begin();
  for(HrMassIntensMap::const_iterator it=MIM.begin();it!=MIM.end();it++)
  {
    *i_mz = it->first;
    *i_int = it->second;
    i_mz++;
    i_int++;
  }
  return(DataFrame::create(Named("mz") = _mz,
                                Named("int")= _int));
}

HrMassIntensMap toMap(const Rcpp::DataFrame & df )
{
    NumericVector _mz = df["mz"];
    NumericVector _int = df["int"];
    int rows = df.nrows();

    HrMassIntensMap MIM;

    for(int row = 0; row < rows; row++)
    {
      MIM[_mz[row]] = _int[row];
    }
    return(MIM);

}
