#include <R.h>

extern "C" {
  void LGARCHRECURSION1 ( int * iStart, int * iEnd, double * phi1, double * theta1, double * yzeroadj, double * innov, double * lny2adj, double * uadj);
}

void LGARCHRECURSION1 ( int * iStart, int * iEnd, double * phi1, double * theta1, double * yzeroadj, double * innov, double * lny2adj, double * uadj) {

  for (int i=*iStart; i < *iEnd; i++) {
    if(yzeroadj[i]==0){
      lny2adj[i] = innov[i] + *phi1*lny2adj[i-1] + *theta1*uadj[i-1];
      uadj[i] = 0;
    }else{
      uadj[i] = lny2adj[i] - innov[i] - *phi1*lny2adj[i-1] - *theta1*uadj[i-1];
    }
  }

}
