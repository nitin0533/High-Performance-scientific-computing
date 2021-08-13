#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 500 

// Declared  new variable type as array
typedef double arr[N];

/******** Code for elementwise matrix multiplication for dot product*********************/

__global__ void dotmat(arr *A, arr *Temp, arr *q, int k,arr *temp) {

    int row=blockIdx.x;
    int col=blockIdx.y;
    if( col<k){
        temp[row][col]=A[row][k]*q[row][col];
     }
}

/******summing up above matrix column wise to get k dot products******************/

__global__ void dot(arr *temp,arr *r,arr *Temp, arr *q, int k) {

    int j;
    int x=threadIdx.x;
    r[0][x]=0.0;
    if(x<k){
        for(j=0;j<N;j++){
          r[0][x]+=temp[j][x];
         }
      }
}

/**********suming over all dot products and respective column of Q**************/
__global__ void submat(arr *LR, arr *Temp, arr *q, int k,arr *r) {

    int row=blockIdx.x;
    int col=blockIdx.y;
     if( col<k){
        LR[row][col]=r[0][col]*q[row][col];
     }
}
/**********Add all the columns element wise into single vector and substract from A[k] and calculate q[k]********************/

__global__ void Qcal(arr *LR, arr *Temp, arr *q, int k) {

  __shared__ double p[N];
  __shared__ double nr;
    int  j;
    int x=threadIdx.x;
    p[x]=0.0;
    for(j=0;j<k;j++){
          p[x]+=LR[x][j];
    }
   __syncthreads();
    Temp[x][k]=Temp[x][k]-p[x];
   __syncthreads();

      if(x==0){
        nr=0.0;
        for(j=0;j<N;j++){
        nr+=Temp[j][k]*Temp[j][k];
        }
        }
   __syncthreads();
        q[x][k]=Temp[x][k]/sqrt(nr);
}

/**********Matrix multiplication AtB***************************************/

__global__ void matmultt(arr *l,arr *m, arr *n)
{
    int x=blockIdx.x;
    int y=blockIdx.y;
    __shared__ double p[N];

    int i;
    int k=threadIdx.x;
    n[x][y]=0;
    p[k]=0;
    p[k]=l[k][x]*m[k][y];

  __syncthreads();
    if(k==0){
     for(i=0;i<N;i++){
        n[x][y]=n[x][y]+p[i];
      }
   }
}

/**************matrix multiplication AB************/
__global__ void matmult(arr *l,arr *m, arr *n)
{
    int x=blockIdx.x;
    int y=blockIdx.y;
    __shared__ double p[N];
    int i;
    int k=threadIdx.x;

    n[x][y]=0;
    p[k]=0;
    p[k]=l[x][k]*m[k][y];
  __syncthreads();
    if(k==0){
       for(i=0;i<N;i++){
          n[x][y]=n[x][y]+p[i];
        }
     }
}

int main(int argc, char** argv) {
    int i,j,k,l,L=250;
    double time_spent = 0.0;
    clock_t begin = clock();
    size_t bytes = N * N * sizeof(double);

    // Allocate memory for our matrices
        arr *A, *q,*Temp,*temp,*MR,*KR,*FR,*LR,*r,*CR,*Q,*phiBt,*KK;
        cudaMallocManaged(&A, bytes);
        cudaMallocManaged(&q, bytes);
        cudaMallocManaged(&Temp, bytes);
        cudaMallocManaged(&temp, bytes);
        cudaMallocManaged(&r, bytes);
        cudaMallocManaged(&CR, bytes);
        cudaMallocManaged(&Q, bytes);
        cudaMallocManaged(&KK, bytes);
        cudaMallocManaged(&phiBt, bytes);
        cudaMallocManaged(&FR, bytes);
        cudaMallocManaged(&LR, bytes);
        cudaMallocManaged(&KR, bytes);
        cudaMallocManaged(&MR, bytes);

/****************Import matries**********************/
        float K[N][N];
        float M[N][N];
        FILE *filek;
        filek=fopen("KG.txt","r");
        if(filek==NULL){
            printf("file doesnt exist");
            return 0;
        }
        while(!feof(filek)){
        printf("entered1");
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){

        fscanf(filek, "%f", &K[i][j]);
        }}}

        FILE *file;
        file=fopen("MG.txt","r");
        if(file==NULL){
            printf("file doesnt exist");
            return 0;
        }
        while(!feof(file)){
        printf("entered2");
         for(i=0;i<N;i++){
         for(j=0;j<N;j++){
               fscanf(file, "%f", &M[i][j]);
         }}}

printf("K of A is:\n");
for(i = N-5; i < N; i++) {
for (j = N-5; j < N; j++) {
           printf("%f ",K[i][j]);
   } printf("\n");
   }

printf("M is:\n");
for(i = 0; i <10; i++) {
for(j = 0; j < 10; j++) {
           printf("%f ",M[i][j]);
   } printf("\n");
   }


//initialize q and Temp
      for(i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
             A[i][j] = M[i][j];
             q[i][j] = 0.0;
             Temp[i][j] = M[i][j];
            KK[i][j] = K[i][j];
        }
    }
//******************setting grid parameters*****************

 dim3 grid(N,N);

for(int x = 0;x< 2; x++) {
//*************QR decompose*********************
   printf("A for %dis:\n",x);
    for(i = 0; i < 10; i++) {
    for(j = 0; j < 10; j++) {
           printf("%f ",A[i][j]);
           } printf("\n");
            }

for(l = 0;l< L; l++) {
      printf("l is %d for x \n",l);

      for(i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        q[i][j] = 0.0;
        MR[i][j] = 0.0;
        Temp[i][j] = A[i][j];
        CR[i][j] = A[i][j];

        }
        }
//Calculate q[0]
        double nr=0.0;
        for(j=0;j<N;j++){
        nr+=A[j][0]*A[j][0];
        }
        for(j=0;j<N;j++){
        q[j][0]=A[j][0]/sqrt(nr);
        }

for(k = 1;k< N; k++) {
      for(i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
           temp[i][j] = 0.0;
           r[0][i] = 0.0;
         }
        }
//**************kernel function for dot products**************

dotmat<<<grid,1>>>(A,Temp,q,k,temp);
cudaDeviceSynchronize();

dot<<<1,N>>>(temp,r,Temp,q,k);
cudaDeviceSynchronize();

submat<<<grid,1>>>(LR,Temp,q,k,r);
cudaDeviceSynchronize();

Qcal<<<1,N>>>(LR,Temp,q,k);
cudaDeviceSynchronize();

/*

        for(j=0;j<N;j++){
        for(i=0;i<N;i++){
        r[0][j]+=temp[i][j];
        }
        }

        for(j=0;j<k;j++){
        for(i=0;i<N;i++){
        Temp[i][k]=Temp[i][k]-r[0][j]*q[i][j];
        }
        }

double  nr=0.0;
        for(j=0;j<N;j++){
        nr+=Temp[j][k]*Temp[j][k];
//      nr+=r[1][j];
=0;j<N;j++){
        q[j][k]=Temp[j][k]/sqrt(nr);
        }
*/
}

//****************QRD ends****************
//*************New A cal******************
      matmultt<<<grid, N>>>(q, A, MR);
      cudaDeviceSynchronize();

      matmult<<<grid,N>>>(MR, q, A);
      cudaDeviceSynchronize();

//***************New Q cal*****************
    if(l<1){
        for(i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
           Q[i][j] =q[i][j];
        }
        }
    }
    if(l>0){
        matmult<<<grid, N>>>(Q, q, KR);
        cudaDeviceSynchronize();
        for(i = 0; i < N; i++) {
        for (j = 0; j < N; j++){
          Q[i][j] =KR[i][j];
        }
        }
    }
//*****************************************

} //l loop

if(x<1){
//*******Cal phiBt**********
for(i = 0; i < N; i++) {
for (j = 0; j < N; j++) {
        LR[i][j] =0.00;
        if(i==j){
        LR[i][i] =1/sqrt(A[i][i]);
        }
     } }

matmult<<<grid,N>>>(Q, LR,phiBt);
cudaDeviceSynchronize();

//***********************transformed K******************
matmultt<<<grid, N>>>(phiBt,KK,FR);
cudaDeviceSynchronize();

matmult<<<grid, N>>>(FR, phiBt,A);
cudaDeviceSynchronize();

}

} //xloop
//**************Eigen vectors********
matmult<<<grid, N>>>(phiBt, Q,CR);
cudaDeviceSynchronize();

printf("Eigenvalues of A are\n");
for(i = N-5; i < N; i++) {
for (j =N-5; j < N; j++) {
           printf("%f ",A[i][j]);
   } printf("\n");
   }

printf("Eigenvectors of A is:\n");
for(i = N-5; i < N; i++) {
for (j = N-5; j < N; j++) {
           printf("%f ",CR[i][j]);
   } printf("\n");
   }


/*#####implementation of mdm#####*/
        int h,modes=10;
        double w=20000.0,dt=0.00001; //frequency
        double phi[N][modes], f[N],fmat[modes],tmat[modes][80],wr[modes],u[N][80];

        //*********setting EV and ev and f ************
        for(i=0;i<modes;i++){
        for(j=0;j<N;j++){
        phi[j][i]=CR[j][N-i-1];
        f[j]=0.0;
        }
        wr[i]=sqrt(A[N-1-i][N-1-i]);
        }
        f[89]=5*pow(10,10);
        f[440]=5*pow(10,10);

        //************cal fmat******
       for(h=0;h<modes;h++){
        fmat[h]=0;
        for(k=0;k<N;k++){
        fmat[h]+=phi[k][h]*f[k];
        }
        }
        //******q cal********
        for(i=0;i<modes;i++){
        for(j=0;j<80;j++){
          tmat[i][j]=fmat[i]*(1/(pow(wr[i],2)-pow(w,2)))*(sin(w*j*dt)-((w/wr[i])*sin(wr[i]*j*dt)));
        }
        }
        //**********u cal*********
        for(h=0;h<N;h++){
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

    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time elpased is %f min\n\n", time_spent/60);

return 0;
}

