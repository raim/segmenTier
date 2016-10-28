#include <Rcpp.h>
using namespace Rcpp;

/// INDIVIDUAL SCORING FUNCTIONS
// only scorecls_c is actually used in algorithm, others for testing

//' Scoring Function "icor" - Test
//' @details  Scoring function "icor" calculates the sum of similarities of
//' positions k:i to cluster c.
//' The similarities are calculated e.g., as a (Pearson) correlation between
//' the individual positions and the tested cluster c center. 
//' This function is for testing only, the dynamic programming algorithm 
//' uses function ccSMicor.
//' NOTE: individual scoring functions are only used for testing
//' the algorithm uses the matrix-filling functions below
//' @param k start position for score calculation
//' @param i end position for score calculation
//' @param c the cluster to which similarities are to be calculated
//' @param seq the cluster sequence (where positions k:i are considered);
//' notably this is not required here, but used as an argument for
//' consistency with other scoring functions.
//' @param M initial penalty parameter, used to set minimal segment sizes
//' @param csim position-cluster similarity matrix, where the rows
//' are the positions in the sequence \code{seq} and columns are the
//' the clusters
//' @export
// [[Rcpp::export]]
double scoreicor_c(int k, int i, int c, NumericVector seq,
		   int M, NumericMatrix csim) {

  if ( c == 0 ) return 0; // nuissance cluster 

  // sum of similarities of positions k:i to cluster c
  double scr = -M;
  for ( ; k<= i; k++ ) 
    scr += csim( k-1, c-1 ); 
  return scr;
}

//' Scoring Function "cor" - Test
//' @details  Scoring function "cor" calculates the sum of similarities
//' between the clusters at positions k:i to cluster c. Note the difference
//' to "icor" where real data from positions are comapred to clusters, while
//' here two clusters are compared.
//' NOTE: This function is for testing only, the dynamic programming algorithm 
//' uses function ccSMcor.
//' @param k start position for score calculation
//' @param i end position for score calculation
//' @param c the cluster to which similarities are to be calculated
//' @param seq the cluster sequence (where clusters at positions k:i are
//' considered)
//' @param M initial penalty parameter, used to set minimal segment sizes
//' @param csim cluster-cluster similarity matrix
//' @export
// [[Rcpp::export]]
double scorecor_c(int k, int i, int c, NumericVector seq,
		  int M, NumericMatrix csim) {

  if ( c == 0 ) return 0; // nuissance cluster 

  // sum of similarities of clusters at positions k:i to cluster c
  double scr = -M;
  for ( ; k<= i; k++ ) 
    if ( seq[k-1]>0 ) 
      scr += csim( seq[k-1]-1, c-1 ); 
  return scr;
}

//' Scoring Function "cls"
//' @details  Scoring function "cls" merely counts the number of
//' of clusters in sequence k:i that are identical to the tested
//' cluster \code{c}, and sub-tracts a minimal size penality
//' and a penalty the for the count of non-identical clusters.
//' NOTE: This function is used in the matrix function ccSMcls.
//' @param k start position for score calculation
//' @param i end position for score calculation
//' @param c the cluster to which similarities are to be calculated
//' @param seq the cluster sequence (where clusters at positions k:i are
//' considered)
//' @param M initial penalty parameter, used to set minimal segment sizes
//' @param a the penalty for non-matching clusters
//' @export
// [[Rcpp::export]]
double scorecls_c(int k, int i, int c, NumericVector seq, int M, int a) {

  // NOTE: this function is used in ccSMcls!!

  if ( c == 0 ) return 0; // nuissance cluster 
  
  int k0=k;
  LogicalVector isC(i-k+1);
  for ( ; k<= i; k++ ) 
    isC[k-k0] = seq[k-1] == c;
  
  double scr = -M + sum(isC) - a*sum(!isC); // counting c and c'
  return scr;
}


/// MATRIX SCORING FUNCTIONS
// these are called in the algorithm and calculate the scoring function
// matrices for individual scoring functions

//' Scoring Function Matrix "icor"
//' @details  Scoring function "icor" calculates the sum of similarities of
//' positions k:i to cluster c over all k and i.
//' The similarities are calculated e.g., as a (Pearson) correlation between
//' the individual positions and the tested cluster c center.
//' Note the difference to "ccor" where the cluster centers are compared
//' instead of original data at positions k and i with a cluster.
//' @param seq the cluster sequence (where positions k:i are considered);
//' notably this argument not required here, but only used for
//' consistency with other scoring functions
//' @param c the cluster to which similarities are to be calculated; note, 
//' that c=1 is the nuissance cluster
//' @param M initial penalty parameter, used to set minimal segment sizes
//' @param Mn penalty parameter for nuissance clusters, should be Mn<M
//' @param csim position-cluster similarity matrix, where the rows
//' are the positions in the sequence \code{seq} and columns are the
//' the clusters
//' @export
// [[Rcpp::export]]
NumericMatrix ccSMicor(NumericVector seq, int c, int M, int Mn,
		       NumericMatrix csim) {

  if ( c==1 ) M = Mn; // nuissance cluster - lower M!
  int nrow = seq.length();
  NumericMatrix SM(nrow,nrow);
  std::fill( SM.begin(), SM.end(), NumericVector::get_na() ) ;
  for (int i = 0; i < nrow; i++) {
    // sum of similarities of positions k:i to cluster c
    SM(i,i) = -M + csim( i, c-1 ); 
    for ( int k = i-1; k >= 0; k-- ) 
      SM(k,i) =  SM(k+1,i) + csim( k, c-1 ); 
  }
  return SM;
}

//' Scoring Function Matrix "ccor" 
//' @details  Scoring function "ccor" calculates the sum of similarities
//' between the clusters at positions k:i to cluster c over all k and i.
//' Note the difference to "icor" where real data from positions are
//' compared to clusters, while here two clusters are compared.
//' @param seq the cluster sequence (where clusters at positions k:i are
//' considered)
//' @param c the cluster to which similarities are to be calculated
//' @param M initial penalty parameter, used to set minimal segment sizes
//' @param Mn penalty parameter for nuissance clusters, should be Mn<M
//' @param csim cluster-cluster similarity matrix
//' @export
// [[Rcpp::export]]
NumericMatrix ccSMccor(NumericVector seq, int c, int M, int Mn, 
		       NumericMatrix csim) {
  
  if ( c==1 ) M = Mn; // nuissance cluster - lower M!
  int nrow = seq.length();
  NumericMatrix SM(nrow,nrow);
  std::fill( SM.begin(), SM.end(), NumericVector::get_na() ) ;
  for (int i = 0; i < nrow; i++) {
    // sum of similarities of clusters at positions k:i to cluster c
    SM(i,i) = -M + csim( seq[i]-1, c-1 );
    for ( int k = i-1; k >= 0; k-- ) 
      SM(k,i) =  SM(k+1,i) + csim( seq[k]-1, c-1 ); 
  }
  return SM;
}


//' like ccSMccor, but also calculates scores for the nuissance cluster 
//' @inheritParams ccSMccor
//' @export
// [[Rcpp::export]]
NumericMatrix ccSMncor(NumericVector seq, int c, int M, NumericMatrix csim) {

  int nrow = seq.length();
  NumericMatrix SM(nrow,nrow);
  std::fill( SM.begin(), SM.end(), NumericVector::get_na() ) ;
  for (int i = 0; i < nrow; i++) {
    // sum of similarities of clusters at positions k:i to cluster c
    SM(i,i) = -M + csim( seq[i]-1, c-1 );
    for ( int k = i-1; k >= 0; k-- ) 
      if ( seq[k]==0 ) SM(k,i) =  SM(k+1,i);
      else SM(k,i) =  SM(k+1,i) + csim( seq[k]-1, c-1 ); 
  }
  return SM;
}


//' Scoring Function Matrix "cls"
//' @details  Scoring function "cls" merely counts the number of
//' of clusters in sequence k:i, over all k and i, that are identical
//' to the tested cluster \code{c}, and sub-tracts a minimal size penality
//' and a penalty the for the count of non-identical clusters.
//' Note: this function used in the scoring unction scorecls_c for individual
//' calculations.
//' @param seq the cluster sequence (where clusters at positions k:i are
//' considered)
//' @param c the cluster to which similarities are to be calculated
//' @param M initial penalty parameter, used to set minimal segment sizes
//' @param Mn not use here, present just for consistency between
//' scoring matrix function calls
//' @param csim integer, the penalty for non-matching clusters
//' @export
// [[Rcpp::export]]
NumericMatrix ccSMcls(NumericVector seq, int c, int M, int Mn, int csim) {

  // note: "Mn" not used, required for consistency with ccSM functions
  // to allow call from R wrapper

  //int a = csim; // note: "csim" required for call from R wrapper (as Mn)

  int nrow = seq.length();
  NumericMatrix SM(nrow,nrow);
  std::fill( SM.begin(), SM.end(), NumericVector::get_na() ) ;
  for (int i = 0; i < nrow; i++) 
    for ( int k = 0; k <= i; k++) 
      SM(k,i) = scorecls_c(k+1, i+1, c, seq, M, csim);
  return SM;
}

//' dynamic programming routine 
//' @details: This is \code{\link{segmenTier}}'s core dynamic programing
//' routine. It takes the scoring function matrices for all clusters,
//'  and dynamically constructs the total score matrix S(i,c).
//' @param seq the cluster sequence (where clusters at positions k:i are
//' considered)
//' @param C the list of clusters
//' @param SM list of scoring function matrices
//' @param multi if multiple \code{k} are found which return the same maximal
//' score, should the "max" (shortest distance) or "min" \code{k} be used?
//' This has little effect on real-life large data sets, since the situation
//' will rarely occur. Default is "max".
//' @param verb level of verbosity, currently not used (TODO: rm?)
//' @export
// [[Rcpp::export]]
List calculateTotalScore(NumericVector seq, NumericVector C, 
			 List SM, String multi="max", int verb=1) {
  // list of scoring function matrices for each cluster 
  List SML(SM);
  
  // result matrices S(i,c) and K(i,c)
  int N = seq.length();
  int M = C.length();
  NumericMatrix S(N,M); // total score matrix
  NumericMatrix K(N,M); // for backtracing
  
  // initialize matrix to 0 and first seq cluster to 1
  // S(0,-1) = 0  wins over S(0,c) = -Inf; 
  std::fill( S.begin(), S.end(), 0.0 ) ;
  std::fill( K.begin(), K.end(), 1 ) ;

  // initialize second position: score(1,2,c)
  for ( int c=0; c<M; c++ ) {
    NumericMatrix SM = SML[c]; // scoring function matrix!
    S(1,c) = SM(0,1);
    K(1,c) = 1;
  }

  // go through sequence of clusters
  // start at position 3, since 1-2 were initialized already
  for ( int i=2; i<N; i++ ) {

    for ( int c=0; c<M; c++ ) {
      
      int kmax = i-1; // 
      
      NumericMatrix SM = SML[c]; // scoring function matrix!
      NumericVector scr(kmax);

      // S(k-1,c') + score(k,i,c)
      // fill from k=0..i
      for ( int k=0; k<kmax; k++ ) {
	// score(k,i,c)
	scr[k] = SM(k,i);  // scorecor_c(k,i,c,seq,M,a);
	// + max_c' S(k-1,c')
	if ( k > 0 ) {
	  // TODO: is -Inf dangerous? smarter solution?
	  double mxsk = - std::numeric_limits<double>::infinity();
	  for ( int cp=0; cp<M; cp++ ) 
	    if ( cp!=c ) 
	      if ( S(k-1,cp) > mxsk ) mxsk = S(k-1,cp);
	  scr[k] += mxsk;
	}
      }
      // max_k ( max_c' S(k-1,c') + score(k,i,c) )
      float mxsc = max( scr );
      S(i,c) = mxsc;

      // store which k was used for back-tracing
      if ( multi=="max" ) {
	NumericVector rcs = clone<NumericVector>(scr);
	std::reverse(rcs.begin(), rcs.end());
	K(i,c) = kmax - which_max( rcs );
      } else {
	K(i,c) = which_max( scr ) + 1;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("S") = S,
			    Rcpp::Named("K") = K);
}

