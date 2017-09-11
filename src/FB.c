#include <R.h>
#include <Rinternals.h>

#include <stdio.h>  /* directives au préprocesseur */
//#include <iostream>
#include <math.h>
#include <stdlib.h>
//#include <fstream>
//#include <string>

//#include <gsl/gsl_sf_log.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_blas.h>

//#include "vecmatfunc.h"


// Convert a vector into the matrix
void VectToMat(double *Vect,double **Mat, int nbcol,int nbrow)
{
  int x,y;
  for(x=0; x<nbrow; x++){
    for(y=0; y<nbcol; y++){
      Mat[x][y] = Vect[x+y*nbrow];
    }
  }
}


// Convert the matrix into a vector
void MatToVect(double **Mat,double *Vect, int nbcol,int nbrow)
{
  int x,y;
  for(x=0; x<nbrow; x++){
    for(y=0; y<nbcol; y++){
      Vect[x+y*nbrow]=Mat[x][y];
    }
  }
}

// production for two matrix
double matProd(double **A,double **B, int m, int q, int lg)
{
  int x;
  double res;
  res=0;

  for(x=0; x<lg; x++){
    res=res+(A[m][x] * B[x][q]);
  }
  return res;
}

// colSums for the matrix
double colSums(double **mat, int m,int lg)
{
  int x;
  double res;
  res=0;

  for(x=0; x<lg; x++){
    res=res+mat[m][x];
  }
  return res;
}


/******************************************************************************/
		/*                    Forward.c	                              */
/******************************************************************************/

// F0
void CalculF0(double *init,double **emiss,double **F, int Q)
{
	//Rprintf("CalculF0\n");

   int x, y;
   double NormF0 = 0;

   for(x=0; x<Q; x++){
     NormF0 = NormF0 +exp(log(init[x]) + log(emiss[0][x]));
   }

   for(y=0; y<Q; y++){
            F[0][y] =  exp(log(init[y])+log(emiss[0][y])-log(NormF0));
   }
}

// Normalise Ft
double NormF(double **trans,double **emiss,double **F, int Q, int t)
{
	//Rprintf("Normalise\n");

    int y;
    double tot;
    tot=0;

    for(y=0; y<Q; y++){
          tot = tot + exp(log(emiss[t][y])+log(matProd(F, trans, t-1, y, Q)));
    }

    return tot;
}


//  Compute  Ft
void CalculF(double **trans,double **emiss,double **F, int Q,int lg)
{
     int t,y;
     double constF = 0;

     for(t=1; t<lg; t++){
       constF =  NormF(trans,emiss,F,Q,t);
          for(y=0; y<Q; y++){
              F[t][y] =  exp(log(emiss[t][y])+ log(matProd(F, trans, t-1, y, Q)) - log(constF));
          }
        //printf("%d ",NormF);
     }
}


//    loglikelihood
void CalculVrais(double **emiss,double *init,double **trans, int Q,int lg,double *Lambda)
{
	//Rprintf("CalculF\n");
int x,y,i;
double **A = (double**)malloc(lg*sizeof(double*));
double **Apast = (double**)malloc(lg*sizeof(double*));

for(i=0;i<lg;i++){
  A[i] = (double*)malloc(Q*sizeof(double));
}

for(i=0;i<lg;i++){
  Apast[i] = (double*)malloc(Q*sizeof(double));
}

for(x=0; x<(lg); x++){
	for(y=0; y<Q; y++){
			Apast[x][y]=0;
			A[x][y]=0;
	}
}

for(y=0; y<Q; y++){
	Apast[0][y] = init[y]*emiss[0][y];
		//printf("%d ",emiss[0][y]);
  		// printf("%d ",init[y]);
		//printf("%d ",A[0][y]);
	}

Lambda[0] = 1/colSums(Apast,0,Q);

for(y=0; y<Q; y++){
Apast[0][y] = Apast[0][y]*Lambda[0];
}

for(x=0; x<(lg-1); x++){
   for(y=0; y<Q; y++){
       A[x+1][y] = matProd(Apast,trans,x,y,Q) * emiss[x+1][y];
   }
    Lambda[x+1] = 1/colSums(A,x+1,Q);

   for(y=0; y<Q; y++){
			Apast[x+1][y] = A[x+1][y]*Lambda[x+1];
		}
}


for(i=0;i<lg;i++){
    free(A[i]) ;
    }

free(A);

for(i=0;i<lg;i++){
    free(Apast[i]) ;
    }
free(Apast);
}



/******************************************************************************/
/*       Forward step          */
/******************************************************************************/

void Forward(double *emissVec, double *init, double *transVec, int *lg, int *Q, double *FVect,double *Lambda)
{
  int i;

  double **emiss = (double**)malloc(*lg*sizeof(double*));
  for(i=0;i<*lg;i++){
      emiss[i] = (double*)malloc(*Q*sizeof(double));
  }

  double **F = (double**)malloc(*lg*sizeof(double*));
  for(i=0;i<*lg;i++){
      F[i] = (double*)malloc(*Q*sizeof(double));
  }

  double **trans = (double**)malloc(*Q*sizeof(double*));
  for(i=0;i<*Q;i++){
      trans[i] = (double*)malloc(*Q*sizeof(double));
      }

  VectToMat(emissVec,emiss, *Q,*lg);
  VectToMat(transVec,trans, *Q,*Q);
  VectToMat(FVect,F, *Q,*lg);

  // Compute F0
  CalculF0(init,emiss,F,*Q);

  // Compute Ft
  CalculF(trans,emiss,F,*Q,*lg);

    // Likelihood calculation
  CalculVrais(emiss,init,trans,*Q,*lg,Lambda);

   MatToVect(emiss,emissVec,*Q,*lg);
   MatToVect(F,FVect,*Q,*lg);
   MatToVect(trans,transVec,*Q,*Q);


 for(i=0;i<*lg;i++){ free(emiss[i]);}

  free(emiss);

 for(i=0;i<*Q;i++){free(trans[i]);}

  free(trans);

 for(i=0;i<*lg;i++){free(F[i]) ;}

  free(F);
}



/******************************************************************************/
		/*                    Backward.c	                              */
/******************************************************************************/

// trans * tau / G
double sumTau(double **trans,double **Tau,double **G, int Q,int j,int t)
{

   int x;
   double res;
   res=0;

   for(x=0; x<Q; x++){
             res=res+ trans[j][x] * Tau[t][x] / G[t][x];
   }

   return res;
}

// Normalise G
void NormG(double **trans, double **F, int Q,int lg,double **G)
{
  int x,y;

  for(x=0; x<(lg-1); x++){
      for(y=0; y<Q; y++){
        G[x+1][y] = matProd(F, trans, x, y, Q);
        G[0][y] = 0;
      }
  }
}


// Tau at times N-1, N-2, ..., 1
void CalculTau(double **trans,double **F, int Q,int lg,double **G,double **Tau)
{

     int t,y;

     NormG(trans,F,Q,lg,G);

     for(t=(lg-2); t>=0; t--){
        for(y=0; y<Q; y++){
          Tau[t][y] =  exp(log(F[t][y])+log(sumTau(trans,Tau,G,Q,y,t+1)));
        }
     }
}

// Tau at times N
void CalculTauN(double **F,int lg, int Q, double **Tau)
{
    int y;
    for(y=0; y<Q; y++){
    Tau[lg-1][y]=F[lg-1][y];
    }
}

/******************************************************************************/
/*   Backward step                             */
/******************************************************************************/

void Backward(double *FVect, double *transVec,int *lg, int *Q, double *TauVect, double *GVect)
{
   int i;

  double **F = (double**)malloc(*lg*sizeof(double*));
  for(i=0;i<*lg;i++){
      F[i] = (double*)malloc(*Q*sizeof(double));
  }

  double **trans = (double**)malloc(*Q*sizeof(double*));
  for(i=0;i<*Q;i++){
       trans[i] = (double*)malloc(*Q*sizeof(double));
  }

  double **Tau = (double**)malloc(*lg*sizeof(double*));
  for(i=0;i<*lg;i++){
       Tau[i] = (double*)malloc(*Q*sizeof(double));
  }

  double **G = (double**)malloc(*lg*sizeof(double*));
  for(i=0;i<*lg;i++){
        G[i] = (double*)malloc(*Q*sizeof(double));
  }

      VectToMat(transVec,trans, *Q,*Q);
      VectToMat(FVect,F, *Q,*lg);
      VectToMat(TauVect,Tau, *Q,*lg);
      VectToMat(GVect,G, *Q,*lg);


    CalculTauN(F,*lg,*Q,Tau);

    CalculTau(trans,F,*Q,*lg,G,Tau);

    MatToVect(F,FVect,*Q,*lg);
    MatToVect(trans,transVec,*Q,*Q);
    MatToVect(Tau,TauVect,*Q,*lg);
    MatToVect(G,GVect,*Q,*lg);


for(i=0;i<*lg;i++){
                  free(F[i]) ;
           }
   free(F);
for(i=0;i<*lg;i++){
                  free(Tau[i]) ;
           }
   free(Tau);
for(i=0;i<*lg;i++){
                  free(G[i]) ;
           }
   free(G);
for(i=0;i<*Q;i++){
                  free(trans[i]) ;
           }
   free(trans);
}


