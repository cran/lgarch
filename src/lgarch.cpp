#include <R.h>

extern "C" {
  void ARMARECURSION1 ( int * iStart, int * iEnd, double * phi1, double * theta1, double * yzeroadj, double * innov, double * lny2adj, double * uadj);
}

void ARMARECURSION1 ( int * iStart, int * iEnd, double * phi1, double * theta1, double * yzeroadj, double * innov, double * lny2adj, double * uadj){

  for (int i=*iStart; i < *iEnd; i++) {
    if(yzeroadj[i]==0){
      lny2adj[i] = innov[i] + *phi1*lny2adj[i-1] + *theta1*uadj[i-1];
      uadj[i] = 0;
    }else{
      uadj[i] = lny2adj[i] - innov[i] - *phi1*lny2adj[i-1] - *theta1*uadj[i-1];
    }
  }

}


extern "C" {
  void LGARCHSIM ( int * maxpq, int * nmaxpq, double * lnsigma2, double * phi, double * phisum, double * innov);
}

void LGARCHSIM ( int * maxpq, int * nmaxpq, double * lnsigma2, double * phi, double * phisum, double * innov) {

  for (int i=*maxpq; i < *nmaxpq; i++) {
    for (int j=0; j < *maxpq; j++) {
      phisum[i] = phisum[i] + phi[j] * lnsigma2[i-j-1];
    }
    lnsigma2[i] = innov[i] + phisum[i];
  }

}


extern "C" {
  void VARMARECURSION1 ( int * iStart, int * n, int * m, double * mU, double * mY, double * mInnov, double * PHI, double * THETA, double * mYiszeroadj);
}

void VARMARECURSION1 ( int * iStart, int * n, int * m, double * mU, double * mY, double * mInnov, double * PHI, double * THETA, double * mYiszeroadj) {

  double PHIfit; PHIfit=0;
  double THETAfit; THETAfit=0;
  double mYfit; mYfit=0;

  /*observation loop:*/
  for(int i=*iStart; i < *n; i++){

    /*dimension loop:*/
    for(int j=0; j < *m; j++){

      /*AR and MA loop:*/
      PHIfit = 0;
      THETAfit = 0;
      for(int k=0; k < *m; k++) {
        PHIfit = PHIfit + PHI[j+k * *m] * mY[i-1+k * *n];
        THETAfit = THETAfit + THETA[j+k * *m] * mU[i-1+k * *n];
      }

      /*sum the fits:*/
      mYfit = mInnov[i+j * *n] + PHIfit + THETAfit;

      /*check for zero:*/
      if(mYiszeroadj[i+j * *n]==1){
        mY[i+j * *n] = mYfit;
      }

      /*recursion (equation-by-equation):*/
      mU[i+j* *n] = mY[i+j * *n] - mYfit;

    } /*end dimension loop*/
  } /*end observation loop*/

}
