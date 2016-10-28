#include <Rcpp.h>
using namespace Rcpp;

// Pearson product-moment correlation coefficient
// 
// [[Rcpp::export]]
double myPearson(NumericVector x, NumericVector y) {
  
  // http://see-programming.blogspot.de/2013/07/c-program-to-calculate-correlation.html

  int n = x.length();
  double  xsum, ysum, xysum, xsqr_sum, ysqr_sum;
  double num, deno, coeff;

  xsum = ysum = xysum = xsqr_sum = ysqr_sum = 0;

  // calculate sums
  for ( int i = 0; i < n; i++) {
    xsum += x[i];
    ysum += y[i];
    xysum += x[i] * y[i]; 
    xsqr_sum += x[i] * x[i]; 
    ysqr_sum += y[i] * y[i]; 
  }
  
  // calculate correlation coefficient 
  num = (n * xysum) - (xsum * ysum);
  deno = (n * xsqr_sum - xsum * xsum)* (n * ysqr_sum - ysum * ysum);
  coeff = num / sqrt(deno);
  
  return coeff;

}
// calculate Pearson correlation of each data to clusters
// and return dataXcluster correlation matrix
//'@export
// [[Rcpp::export]]
NumericMatrix clusterCor_c(NumericMatrix data, NumericMatrix clusters) {

  int N = data.nrow();
  int d = data.ncol();
  int C = clusters.nrow();
  NumericMatrix clCor(N,C); // data X median correlation matrix
  
  for ( int i=0; i<N; i++ ) {
    NumericVector x(d);
    NumericVector y(d);
    for ( int k=0; k<d; k++ )  
      x(k) = data(i,k);
    // find max. correlating cluster median
    for ( int c=0; c<C; c++ ) {
      for ( int k=0; k<d; k++ )  
	y(k) = clusters(c,k); 
      clCor(i,c) = myPearson(x,y); //cor(x,y);
    }   
  }
  return clCor;
}
// calculate Pearson correlation of each data to clusters
// and report the (first!) cluster which had max correlation
//'@export
// [[Rcpp::export]]
NumericVector clusterMaxCor_c(NumericMatrix data, NumericMatrix clusters, float mincor=0.0, int warn=0) {

  int N = data.nrow();
  int d = data.ncol();
  int C = clusters.nrow();
  NumericVector cls(N); // should be cluster assignment by highest cor
  
  for ( int i=0; i<N; i++ ) {
    NumericVector x(d);
    NumericVector y(d);
    for ( int k=0; k<d; k++ )  
      x(k) = data(i,k);
    // find max. correlating cluster median
    NumericVector clscor(C);
    for ( int c=0; c<C; c++ ) {
      for ( int k=0; k<d; k++ )  
	y(k) = clusters(c,k); 
      clscor[c] = myPearson(x,y); //cor(x,y);
    }   
    // find max clscor
    // NOTE: not unique!
    double mxcor = max(clscor);
    if ( mxcor < mincor ) {
      cls[i] = 0;
    } else {
      // check multiple?
      if ( warn== 1) {
	int cnt=0;
	for ( int m=0; m<C; m++ ) 
	  if ( clscor[m] == mxcor )
	    cnt++;
        if ( cnt>1 ) 
          Rcpp::warning("%d mxcor at data row %d; mxcor=%f", cnt, i+1, mxcor);
      }
      int m = C; 
      while( m-- >= 0 )  // find first - not unique!
	if ( clscor[m] == mxcor )  {
    	  cls[i] = m+1;
	  break;
	}
    }
  }
  return cls;
  
}
