#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//
// C++ main function `rcppOSS`
//
// R function `OSS` core code, `Rcpp`-version by `the package itself`.
//
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// Scale a matrix
//
// @param X A matrix.
//
// @return Scaled matrix.
//
// NumericMatrix ScaleMatrix(NumericMatrix X){
//   int p = X.cols();
//   for(int j = 0; j < p; j++){
//     X(_,j) = (X(_,j) - mean(X(_,j))) / sd(X(_,j));
//   }
//   return X;
// }

//' Get L2 norm
//'
//' Get L2 norm of a matrix or data frame.
//'
//' @param X A matrix or data.frame.
//'
//' @return L2 norm of `X`(every row).
// [[Rcpp::export]]
NumericVector L2norm(NumericMatrix X){
  int N = X.rows();
  NumericVector norm(N);
  for(int i=0; i<N; i++){
    norm[i] = sum(X(i,_) * X(i,_));
    // norm[i] = std::inner_product(X(i,_).begin(), X(i,_).end(), X(i,_).begin(),0);
  }
  return norm;
}

//' Find t smallest index of a vector
//'
//' @param loss A vector.
//' @param t A int.
//'
//' @return The index of the t smallest element of the vector.
// [[Rcpp::export]]
IntegerVector bottom_t_index(NumericVector loss, int t){
  IntegerVector remain(t);
  NumericVector losscopy = clone(loss);
  std::nth_element(losscopy.begin(), losscopy.begin() + t - 1, losscopy.end());
  for(int i=0, ii=0; i<loss.length(); i++){
    if (loss[i] <= losscopy[t-1]) remain[ii++] = i;
  }
  return(remain);
}

//' Compute loss function for OSS
//'
//' @param candi The index of the candidate set.
//' @param last_index The index of the seleted point in last iteration.
//' @param X The whole data.
//' @param norm Norm of the whole data.
//'
//' @return Loss of every point in candidate set.
// [[Rcpp::export]]
NumericVector ComputeLoss(IntegerVector candi, int last_index, NumericMatrix X, NumericVector norm){
  int p = X.cols();
  int k = candi.length();
  NumericVector loss(k);
  for(int i=0; i<k; i++){
    int delta = sum(sign(X(candi[i],_)) == sign(X(last_index,_)));
    loss[i] = pow(p - norm[candi[i]]/2 - norm[last_index]/2 + delta, 2);
  }
  return loss;
}

//' OSS `Rcpp`-version by `the package itself` (`OSS` core code)
//'
//' @param X A matrix.
//' @param n Subsample size.
//'
//' @return Subsample index.
// [[Rcpp::export]]
IntegerVector rcppOSS(NumericMatrix X, int n) {
  // X = ScaleMatrix(X);  // Standardize in R
  int N = X.rows();

  IntegerVector index(n);
  IntegerVector candi = seq_len(N)-1;
  NumericVector norm = L2norm(X);
  double r = log(N) / log(n);

  index[0] = which_max(norm);
  candi.erase(index[0]);
  NumericVector loss = ComputeLoss(candi, index[0], X, norm);
  // Rcout << "i = 0" << "\n";
  // Rcout << "index: " << index[0] << " \n";
  // Rcout << "length of candi: " << candi.length() << " \n";
  // Rcout << "candi: " << candi << "\n";
  // Rcout << "length of loss: " << loss.length() << " \n";
  // Rcout << "loss: " << loss << "\n" << "\n";

  for(int i=1; i<n; i++){

    int tmp = which_min(loss);
    index[i] = candi[tmp];
    // Rcout << "i = " << i << "\n";
    // Rcout << "tmp: " << tmp << "\n";
    // Rcout << "index: " << index[i] << "\n";

    candi.erase(tmp);
    loss.erase(tmp);
    // Rcout << "length of candi: " << candi.length() << "\n";
    // Rcout << "candi: " << candi << "\n";
    // Rcout << "length of loss: " << loss.length() << "\n";
    // Rcout << "loss: " << loss << "\n";

    double t = 0;
    if (N > pow(n,2)) {
      t = (N / (i+1)) * 3;
    } else {
      t = (N / pow(i+1, r-1)) * 3;
    }
    // Rcout << "t = " << t << "\n";
    if (candi.length() > t) {
      IntegerVector remain = bottom_t_index(loss, floor(t));
      //   Rcout << "length of remain: " << remain.length() << "\n";
      //   Rcout << "remain:" << remain << "\n";
      candi = candi[remain];
      loss = loss[remain];
    }

    // Rcout << "After eliminate:" << "\n";
    // Rcout << "length of candi: " << candi.length() << "\n";
    // Rcout << "candi: " << candi << "\n";
    // Rcout << "length of loss: " << loss.length() << "\n";
    // Rcout << "loss: " << loss << "\n";
    loss = loss + ComputeLoss(candi, index[i], X, norm); // loss can't located.
    // Rcout << "After Update loss:" << "\n";
    // Rcout << "length of loss: " << loss.length() << "\n";
    // Rcout << "loss: " << loss << "\n" << "\n";
  }
  return index + 1;
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//
// C++ main function `armaOSS`
//
// R function `armaOSS` core code, `RcppArmadillo`-version by `Zhu`.
//
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

//' Find t smallest index of a vector (RcppArmadillo-version)
//'
//' @param x A vector.
//' @param k A int.
//'
//' @return The index of the t smallest element of the vector.
// [[Rcpp::export]]
arma::vec armabottom_k(arma::vec x, unsigned int k) {
  arma::vec x2 = x; // save a copy of x
  arma::vec ind(k); // save the indexes of the smallest k numbers
  std::nth_element(x.begin(), x.begin() + k - 1, x.end()); // std::greater<double>());
  for(int ii=0, i=0; i<int(x.n_elem) && ii<int(k); i++){
    if(x2[i] <= x[k-1])  ind[ii++] = i;  // +1 for R
  }
  return ind;
}

//' Scale a matrix (RcppArmadillo-version)
//'
//' @param X A matrix.
//'
//' @return Scaled matrix.
// [[Rcpp::export]]
arma::mat armaScaleMatrix(arma::mat X){
  // new add to scale (GOSS do not have)
  int p = X.n_cols;
  for(int j = 0; j < p; j++){
    X.col(j) = (X.col(j) - mean(X.col(j))) / stddev(X.col(j));
  }
  return X;
}

//' Compute loss function for OSS (RcppArmadillo-version)
//'
//' @param X Matrix of the candidate set.
//' @param xa Norm of the candidate set.
//' @param y A vector. The point which be selected last iteration.
//' @param ya Norm of `y`.
//' @param tPow The power of the loss function.
//'
//' @return Loss of the candidate set.
// [[Rcpp::export]]
arma::vec armaComputeLoss(arma::mat X, arma::vec xa, arma::mat y, double ya, double tPow) {
  int n=X.n_rows;
  int p=X.n_cols;
  arma::vec B = zeros<vec>(n);
  for(int i=0; i<n; i++){
    B(i) = pow(accu(X.row(i)==y)+p-xa(i)/2-ya/2,tPow);
  }
  return B;
}

//' OSS `RcppArmadillo`-version by `Zhu` (`myArma_OSS` core code)
//' @param x A matrix.
//' @param k Subsample size.
//' @param tPow The power of the loss function.
//'
//' @return Subsample index.
// [[Rcpp::export]]
arma::uvec armaOSS(arma::mat x, int k, double tPow=2){
  // x = armaScaleMatrix(x); // Standardize in R
  int n=x.n_rows;
  arma::uvec candi=linspace<uvec>(1,n,n);
  arma::uvec ind=linspace<uvec>(1,k,k);
  arma::vec L=sum(pow(x,2),1);
  arma::vec xa=L;
  uword mm=L.index_max();
  ind(0)=candi(mm);
  candi.shed_row(mm);
  L.shed_row(mm);

  // Rcout << "i = 0" << "\n";
  // Rcout << "index: " << ind(0)-1 << "\n";
  // Rcout << "length of candi: " << candi.n_elem << "\n";
  // Rcout << "candi: " << candi.t()-1 << "\n";

  arma::mat sx=sign(x);

  /* GOSS original:
   double r=log(n/k)/log(k);
   */
  double r=log(n)/log(k); // modified
  for(int i=1; i<k; i++){
    if(i==1)
      L = armaComputeLoss(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);
    else
      L = L + armaComputeLoss(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);

    mm=L.index_min();
    ind(i)=candi(mm);
    candi.shed_row(mm);
    L.shed_row(mm);

    // Rcout << "i = " << i << "\n";
    // Rcout << "index: " << ind(i)-1 << "\n";
    // Rcout << "length of candi: " << candi.n_elem << "\n";
    // Rcout << "candi: " << candi.t() - 1 << "\n";
    // Rcout << "length of loss: " << L.n_elem << "\n";
    // Rcout << "loss: " << L.t() << "\n";

    /* GOSS original:
     int nc=floor(n/pow(i,r));
     */
    double nc = n/pow(i+1,r-1); // modified

    /* GOSS original:
     if((i>1) & (L.n_elem>double(nc))){
     arma::uvec tt=arma::conv_to<arma::uvec>::from(bottom_k(L,nc));
     L=L.elem(tt);
     candi=candi.elem(tt);
     */
    if( L.n_elem > nc ){ //modified
      arma::uvec tt=arma::conv_to<arma::uvec>::from(armabottom_k(L,floor(nc)));
      L=L.elem(tt);
      candi=candi.elem(tt);

    // Rcout << "t = " << nc << "\n";
    // Rcout << "length of remain: " << tt.n_elem << "\n";
    // Rcout << "remain: " << tt.t() << "\n";
    // Rcout << "After eliminatie" << "\n";
    // Rcout << "length of candi: " << candi.n_elem << "\n";
    // Rcout << "candi: " << candi.t()-1 << "\n";
    // Rcout << "length of loss: " << L.n_elem << "\n";
    // Rcout << "loss: " << L.t() << "\n";

    }
  }
  return ind;
}

