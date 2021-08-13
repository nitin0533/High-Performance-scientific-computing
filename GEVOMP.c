#include<time.h>
#include<stdio.h>
#include <stdio.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
int main(){
double start;
double end;
start = omp_get_wtime();

int d,s,h,l,i,j,m,n,k,L=250;
#define n 500
float** K=malloc(sizeof(float*)*n);
float** M=malloc(sizeof(float*)*n);

for(j=0;j<n;j++){
        K[j]=malloc(sizeof(float)*n);
        M[j]=malloc(sizeof(float)*n);

}
//double K[3][3]={{4.0,1.0,-3.0},{1.0,5.0,-4.0},{-3.0,-4.0,7.0}};
//double M[3][3]={{1.3,1.4,1.4},{1.4,6.9,6.3},{1.4,6.3,7.4}};

//Import matries

        FILE *file;
        file=fopen("KG.txt","r");
        if(file==NULL){
            printf("file doesnt exist");
            return 0;
        }
        while(!feof(file)){
                for(i=0;i<n;i++){
                    for(j=0;j<n;j++){
                       fscanf(file, "%f", &K[i][j]);
                }
                }
         }

        file=fopen("MG.txt","r");
         if(file==NULL){
            printf("file doesnt exist");
            return 0;
         }
         while(!feof(file)){
                for(i=0;i<n;i++){
                    for(j=0;j<n;j++){
                       fscanf(file, "%f", &M[i][j]);
                    }
        }
        }
double nr;

//****************Dynamic Memory Allocations**************

double** A=malloc(sizeof(double*)*n);
double** q=malloc(sizeof(double*)*n);
double** temp=malloc(sizeof(double*)*n);
double** CR=malloc(sizeof(double*)*n);
double** Q=malloc(sizeof(double*)*n);
double** Temp=malloc(sizeof(double*)*n);
double** P=malloc(sizeof(double*)*n);
double** phiBt=malloc(sizeof(double*)*n);
//double** Minv=malloc(sizeof(double*)*n);

double err=2*n;
for(j=0;j<n;j++){
   A[j]=malloc(sizeof(double)*n);
        q[j]=malloc(sizeof(double)*n);
        temp[j]=malloc(sizeof(double)*n);
        CR[j]=malloc(sizeof(double)*n);
        Q[j]=malloc(sizeof(double)*n);
        Temp[j]=malloc(sizeof(double)*n);
        P[j]=malloc(sizeof(double)*n);
      phiBt[j]=malloc(sizeof(double)*n);
    //  Minv[j]=malloc(sizeof(double)*n);
}

#pragma omp parallel for private(j)
for(i=0;i<n;i++){
for(j=0;j<n;j++){
        A[i][j]=M[i][j];
        temp[i][j]=M[i][j];
}
}
 printf("A allocated \n");

for(d=0;d<2;d++){
int l=0;
double error;
while(l<L){   //iterations of QR algorithm...error value checks whether matrix is completely transformed

//********QR decomposition********

//for(int l=0;l<L;l++){
for(i=0;i<n;i++){
        if(i>0){
                for(k=1;k<=i;k++){
                        nr=0.0;
                        #pragma omp parallel for reduction(+:nr)
                        for(j=0;j<n;j++){
                        nr+=A[j][i]*q[j][k-1];
                        }
                        #pragma omp parallel for
                        for(j=0;j<n;j++){
                       temp[j][i]=temp[j][i]-nr*q[j][k-1];
                        }
                }
                nr=0.0;
                #pragma omp parallel for reduction(+:nr)
                for(j=0;j<n;j++){
                nr+=temp[j][i]*temp[j][i];
                }
                #pragma omp parallel for
                for(j=0;j<n;j++){
                q[j][i]=temp[j][i]/sqrt(nr);
                }
        }
        if(i==0){
                nr=0.0;
                #pragma omp parallel for reduction(+:nr)
                for(j=0;j<n;j++){
                nr+=A[j][0]*A[j][0];
                }
                #pragma omp parallel for
                for(j=0;j<n;j++){
                q[j][0]=A[j][0]/sqrt(nr);
                }
        }
}

//**************QR decomposition end*******

//**********New A calculation QtAQ**********
#pragma omp parallel for private(j,k)
for(h=0;h<n;h++){
        for(j=0;j<n;j++){
        CR[h][j]=0;
//#pragma omp parallel for
                for(k=0;k<n;k++){
                CR[h][j]+=q[k][h]*A[k][j];
                }
        }
}
#pragma omp parallel for private(j,k)
for(h=0;h<n;h++){
        for(j=0;j<n;j++){
        Temp[h][j]=0.0;
                for(k=0;k<n;k++){
                Temp[h][j]+=CR[h][k]*q[k][j];
                }
        A[h][j]=Temp[h][j];
        temp[h][j]=A[h][j];
        }
}//***********Accumulate Q=Qnew*Qold********
if(l==0){
        #pragma omp parallel for private(j)
        for(h=0;h<n;h++){
                for(j=0;j<n;j++){
                Q[h][j]=q[h][j];
                }
        }
}

if(l>0){
        #pragma omp parallel for private(j,k)
        for( h=0;h<n;h++){
                for(j=0;j<n;j++){
                P[h][j]=0.0;
                        for(k=0;k<n;k++){
                        P[h][j]+=Q[h][k]*q[k][j];
                        }
        }
}
        #pragma omp parallel for private(j)
        for(h=0;h<n;h++){
                for(j=0;j<n;j++){
                Q[h][j]=P[h][j];
                }
        }

}

//***********checking convergence***********
double err=0.0;
for(h=0;h<n;h++){
        #pragma omp parallel for reduction(+:err)
        for(j=0;j<n;j++){
        if(h!=j){
        err+=A[h][j]*A[h][j];
        }       }
}
error=sqrt(err);
l++;
if(error<1e-4){
break;
}
printf("ml=%d \n",l);
}

if(d==0){
        //********************Calculate phiBt using eq 38 , phiB is identity matrix *****************
        for(i=0;i<n;i++){
        #pragma omp parallel for private(j)
        for(j=0;j<n;j++){
                phiBt[i][j]=0.0;
                phiBt[i][j]=Q[i][j]/sqrt(A[j][j]);
        }
        }

        //************************Calculate transformed K matrix (Eq 39)*****************************
        for(h=0;h<n;h++){
        #pragma omp parallel for private(j,k)
                for(j=0;j<n;j++){
                CR[h][j]=0;
                        for(k=0;k<n;k++){
                        CR[h][j]+=phiBt[k][h]*K[k][j];
                        }
                }
        }
        for(h=0;h<n;h++){
            #pragma omp parallel for private(j,k)
                for(j=0;j<n;j++){
                Temp[h][j]=0.0;
                        for(k=0;k<n;k++){
                        Temp[h][j]+=CR[h][k]*phiBt[k][j];
                        }
                A[h][j]=Temp[h][j];
                temp[h][j]=A[h][j];
                }
        }
}

}


/***************GEVP solved**************/


//***********Eigen values of transformed K is same as generalized eigen values Eq 46*********
printf("Eigen value are \n");
for(i=n-10;i<n;i++){
        printf("%f \n", A[i][i]);
        //printf("\n");
}


//**************Calculate Eigen vectors using Eq 43********************
#pragma omp parallel for private(j,k)
for(h=0;h<n;h++){
for(j=0;j<n;j++){
CR[h][j]=0;
for(k=0;k<n;k++){
CR[h][j]+=phiBt[h][k]*Q[k][j];
}
}
}

printf("EV are \n");
for(i=n-10;i<n;i++){
for(j=n-10;j<n;j++){
        printf("%f ", CR[i][j]);
}
        printf("\n");
}

/*#####implementation of mdm#####*/
        int modes=10;
        double w=20000.0,dt=0.00001; //frequency
        double phi[n][modes], f[n],fmat[modes],tmat[modes][80],wr[modes],u[n][80];

        //*********setting EV and ev and f ************
        #pragma omp parallel for private(j)     
        for(i=0;i<modes;i++){
        for(j=0;j<n;j++){
        phi[j][i]=CR[j][n-i-1];
        f[j]=0.0;
        }
        wr[i]=sqrt(A[n-1-i][n-1-i]);
        }
        f[89]=5*pow(10,10);
        f[440]=5*pow(10,10);

        //************cal fmat******
        #pragma omp parallel for private(k)
        for(h=0;h<modes;h++){
        fmat[h]=0;
        for(k=0;k<n;k++){
        fmat[h]+=phi[k][h]*f[k];
        }
        }
        //******q cal********
        #pragma omp parallel for private(j)
        for(i=0;i<modes;i++){
        for(j=0;j<80;j++){
          tmat[i][j]=fmat[i]*(1/(pow(wr[i],2)-pow(w,2)))*(sin(w*j*dt)-((w/wr[i])*sin(wr[i]*j*dt)));
        }
        }
        //**********u cal*********
        #pragma omp parallel for private(j,k)
        for(h=0;h<n;h++){
                for(j=0;j<80;j++){
                u[h][j]=0.0;
                        for(k=0;k<modes;k++){
                        u[h][j]+=phi[h][k]*tmat[k][j];
                        }
        }
        }
        printf("response\n");
        for(h=0;h<80;h++){
        printf("%8.8f ",u[89][h]);
        }

end = omp_get_wtime();
printf("\n");
printf("Time elapsed %f mins\n", (end - start)/60);
return 0;
}
          

