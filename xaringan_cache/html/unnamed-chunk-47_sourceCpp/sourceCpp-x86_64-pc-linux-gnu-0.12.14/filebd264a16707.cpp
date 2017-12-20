#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix PairGRM(const NumericMatrix &G,
                      const NumericVector &p) {

  // In our algorithms, individuals are stored in columns 
  // and SNPs are stored in rows.
   
  int nSNP = G.nrow(); // number of SNPs
  int nIND = G.ncol(); // number of individuals
  NumericMatrix GRM(nIND, nIND); // Genetic Relationship Matrix

  for (int i = 0; i < nIND; i++) {
    for (int j = 0; j < nIND; j++) {

      // value = GRM(i, j)
      double value = 0;
      double tmp = 0; 
      // number of missing values for each pair 
      // of individuals (i, j)
      // nbmv = \sum_{k = 1}^{nSNP} \delta_{ik} \delta{jk}
      int nbmv = 0;

      // Loop over the SNPs to compute the dot product
      for (int k = 0; k < nSNP; k++) {
        if ((!NumericVector::is_na(G(k, i))) && 
            (!NumericVector::is_na(G(k, j)))) {
          tmp = (G(k, i) - 2 * p[k]) * (G(k, j) - 2 * p[k]);
          value +=  tmp / (2 * p[k] * (1 - p[k]));
        } else {
          nbmv++;
        }
      }
      // Divide by the number of non-missing values for (i, j)
      GRM(i, j) = value / (nSNP - nbmv);
    }
  }
  return GRM;
}
