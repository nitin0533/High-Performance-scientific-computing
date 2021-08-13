#include <stdio.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"


int main(int argc, char** argv){

int slot, rem,my_PE_num;
int core;

float f, start, end;
MPI_Status status;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);
MPI_Comm_size( MPI_COMM_WORLD, &core );
start=MPI_Wtime();

int i,j,m,n,k,L=250;
#define n 500
float** K=malloc(sizeof(float*)*n);
float** M=malloc(sizeof(float*)*n);

for(j=0;j<n;j++){
        K[j]=malloc(sizeof(float)*n);
        M[j]=malloc(sizeof(float)*n);
}


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

//double K[3][3]={{4.0,1.0,-3.0},{1.0,5.0,-4.0},{-3.0,-4.0,7.0}};
//double M[3][3]={{1.3,1.4,1.4},{1.4,6.9,6.3},{1.4,6.3,7.4}};


double nr;

//**********************************************************************************************************//

if (my_PE_num==0){

double** A=malloc(sizeof(double*)*n);
double** q=malloc(sizeof(double*)*n);
double** temp=malloc(sizeof(double*)*n);
double** CR=malloc(sizeof(double*)*n);
double** Q=malloc(sizeof(double*)*n);
double** Temp=malloc(sizeof(double*)*n);
double** P=malloc(sizeof(double*)*n);
double** PJ=malloc(sizeof(double*)*n);
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
        PJ[j]=malloc(sizeof(double)*n);
      phiBt[j]=malloc(sizeof(double)*n);
    //  Minv[j]=malloc(sizeof(double)*n);
}

//Import mass matrix
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
//double M[5][5]={{10,1.4,1.4,1.8,7},{1.4,30,6.3,2,5},{1.4,6.3,7.4,2.8,3},{1.8,2,2.8,4,2},{7,5,3,2,30}};

for(i=0;i<n;i++){
for(j=0;j<n;j++){
        A[i][j]=M[i][j];
        Temp[i][j]=A[i][j];
}}
printf("A allocated \n");
//printf("A is\n"); 
//      for(i=0;i<n;i++){
//      for(j=0;j<n;j++){
//              printf("%f ", A[i][j]);
//      }

int h,x;
for(x=0;x<2;x++){

//***************QR algorithm***********************
int l=0;
double error;
int s,d;
for(l=0;l<L;l++){
printf("l %d is %d\n",x,l);
                for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                        q[i][j]=0.0;
                        temp[i][j]=0.0;
                }}
                for(j=1;j<core;j++){
                for(s=0;s<n;s++){
                        MPI_Send( &(A[s][0]), n, MPI_DOUBLE,j, 10000*l+s, MPI_COMM_WORLD);
                }}
//Calculate q[0]
                double nr=0.0;
                int f;
                for(j=0;j<n;j++){
                nr+=A[j][0]*A[j][0];
                }
                for(j=0;j<n;j++){
                q[j][0]=A[j][0]/sqrt(nr);
                }
//send q0
                for(j=1;j<core;j++){
                for(s=0;s<n;s++){
                        MPI_Send( &(q[s][0]), 1, MPI_DOUBLE,j, 10000*l+s, MPI_COMM_WORLD);
                }}

//Cal all other columns of q
        for(i=1;i<n;i++){
                slot=floor(i/core);
                int rem=i-slot*(core-1);
                for(d=0;d<n;d++){
                for(j=0;j<n;j++){
                        temp[d][j]=0.0;
                }}

                for(k=1;k<=rem;k++){
                if(k<=i){
                        nr=0.0;
                for(j=0;j<n;j++){
                        nr+=A[j][i]*q[j][k-1];
                        }
                        for(j=0;j<n;j++){
                        temp[j][0]+=nr*q[j][k-1];
                        }}
                }
        //recv temp
                for(j=1;j<core;j++){
                for(s=0;s<n;s++){
                MPI_Recv( &(temp[s][j]), 1, MPI_DOUBLE,j,10000*l+s, MPI_COMM_WORLD,&status);
                }}
        //cal new Q
                for(j=0;j<core;j++){
                for(s=0;s<n;s++){
                Temp[s][i]=Temp[s][i]-temp[s][j];
                }}
                double nr=0.0;
                int f;
                for(j=0;j<n;j++){
                nr+=Temp[j][i]*Temp[j][i];
                }

                for(j=0;j<n;j++){
                q[j][i]=Temp[j][i]/sqrt(nr);
                }
        //send Q
                for(j=1;j<core;j++){
                for(s=0;s<n;s++){
                        MPI_Send( &(q[s][i]), 1, MPI_DOUBLE,j, 10000*l+s, MPI_COMM_WORLD);
                }}}//i loop

//New A calculation

        double slots=floor(n/core);
        //CR cal 
                for(i=(core-1)*slots;i<n;i++){
                for(j=0;j<n;j++){
                CR[i][j]=0;
                for(k=0;k<n;k++){
                CR[i][j]+=q[k][i]*A[k][j];
                }}}
        //CR receive
                for(k=1;k<core;k++){
                for(j=(k-1)*slots;j<(k*slots);j++){
                MPI_Recv( &(CR[j][0]), n, MPI_DOUBLE,k, j, MPI_COMM_WORLD,&status);
                }}
        //A cal 
                for(i=(core-1)*slots;i<n;i++){
                for(j=0;j<n;j++){
                P[i][j]=0;
                for(k=0;k<n;k++){
                P[i][j]+=CR[i][k]*q[k][j];
                }}}
        //A receive
                for(k=1;k<core;k++){
                for(j=(k-1)*slots;j<(k*slots);j++){
                MPI_Recv( &(P[j][0]), n, MPI_DOUBLE,k, j, MPI_COMM_WORLD,&status);
                }}

                for(k=0;k<n;k++){
                for(j=0;j<n;j++){
                A[k][j]=P[k][j];
                Temp[k][j]=P[k][j];
                }}

//Accumulate Q
        if(l==0){
                for(h=0;h<n;h++){
                        for(j=0;j<n;j++){
                        Q[h][j]=q[h][j];
                        }}
}
        if(l>0){
        //PJ cal 
                for(i=(core-1)*slots;i<n;i++){
                for(j=0;j<n;j++){
                PJ[i][j]=0;
                for(k=0;k<n;k++){
                PJ[i][j]+=Q[i][k]*q[k][j];
                }}}

                for(h=0;h<n;h++){
                        for(j=0;j<n;j++){
                        Q[h][j]=PJ[h][j];
                        }}
        //PJ receive
                for(k=1;k<core;k++){
                for(j=(k-1)*slots;j<(k*slots);j++){
                MPI_Recv( &(Q[j][0]), n, MPI_DOUBLE,k, j, MPI_COMM_WORLD,&status);
                }}}//l if
}//l loop
//      printf("A from %d is\n",my_PE_num); 
//              for(k=n-5;k<n;k++){
//              for(j=n-5;j<n;j++){
//                      printf("%f ", A[k][j]);
//              }
//                      printf("\n");
//              }

//      printf("Q from %d is\n",my_PE_num); 
//              for(k=0;k<n;k++){
//              for(j=0 ;j<n;j++){
//                      printf("%f ", Q[k][j]);
//              }printf("\n");
//              }

if(x<1){
        //********************Calculate phiBt using eq 38 , phiB is identity matrix *****************
        for(i=0;i<n;i++){
        for(j=0;j<n;j++){
                phiBt[i][j]=0.0;
                phiBt[i][j]=Q[i][j]/sqrt(A[j][j]);
        }
        }

                for(j=1;j<core;j++){
                for(s=0;s<n;s++){
                        MPI_Send( &(phiBt[s][0]), n, MPI_DOUBLE,j, 10000*l+s, MPI_COMM_WORLD);
                }}

        //************************Calculate transformed K matrix (Eq 39)*****************************


//New K calculation

        double slots=floor(n/core);
        //CR cal 
                for(i=(core-1)*slots;i<n;i++){
                for(j=0;j<n;j++){
                CR[i][j]=0;
                for(k=0;k<n;k++){
                CR[i][j]+=phiBt[k][i]*K[k][j];
                }}}
        //CR receive
                for(k=1;k<core;k++){
                for(j=(k-1)*slots;j<(k*slots);j++){
                MPI_Recv( &(CR[j][0]), n, MPI_DOUBLE,k, j, MPI_COMM_WORLD,&status);
                }}
        //K cal 
                for(i=(core-1)*slots;i<n;i++){
                for(j=0;j<n;j++){
                P[i][j]=0;
                for(k=0;k<n;k++){
                P[i][j]+=CR[i][k]*phiBt[k][j];
                }}}
        //K receive
                for(k=1;k<core;k++){
                for(j=(k-1)*slots;j<(k*slots);j++){
                MPI_Recv( &(P[j][0]), n, MPI_DOUBLE,k, j, MPI_COMM_WORLD,&status);
                }}

                for(k=0;k<n;k++){
                for(j=0;j<n;j++){
                A[k][j]=P[k][j];
                Temp[k][j]=A[k][j];
                }}
} // if x
//      printf("A2 from %d is\n",my_PE_num); 
//              for(k=n-5;k<n;k++){
//              for(j=n-5;j<n;j++){
//                      printf("%f ", A[k][j]);
//              }
//                      printf("\n");
//              }

}// x loop
        for(h=0;h<n;h++){
        for(j=0;j<n;j++){
        CR[h][j]=0;
        for(k=0;k<n;k++){
        CR[h][j]+=phiBt[h][k]*Q[k][j];
        }
        }
}

        printf("ev from %d is\n",my_PE_num);
                for(k=n-10;k<n;k++){
                for(j=n-10;j<n;j++){
                        printf("%f ", A[k][j]);
                }
                        printf("\n");
                }
/*
        printf("EV from %d is at 500\n",my_PE_num); 
                for(j=0;j<50;j++){
                        printf("%f\n", CR[j][500]);
                }

        printf("EV from %d is at 0\n",my_PE_num); 
                for(j=0;j<50;j++){
                        printf("%f\n ", CR[j][0]);
                }
*/
///////////////*********MDM************************************
        int modes=10;
        double w=20000.0,dt=0.00001; //frequency
        double phi[n][modes], f[n],fmat[modes],tmat[modes][80],wr[modes],u[n][80];

        //*********setting EV and ev and f ************
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
        for(h=0;h<modes;h++){
        fmat[h]=0;
        for(k=0;k<n;k++){
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

        FILE *fp;
        fp = fopen ("Vector.txt","w");
        for(i = 0; i < n;i++){
        for(j = 0; j < n;j++){
            fprintf (fp, "%f ",CR[i][j]);
        }  fprintf (fp,"\n");
        }
        fclose (fp);
        FILE *fv;
        fv = fopen ("Values.txt","w");
        for(i = 0; i < n;i++){
          fprintf (fv, "%f\n",A[i][i]);
        }
        fclose (fv);
//***********************MDM********************

}//if loop
//*************************************************************************************************************//
if (my_PE_num>0){

        double** LA=malloc(sizeof(double*)*n);
        double** Lq=malloc(sizeof(double*)*n);
        double** LQ=malloc(sizeof(double*)*n);
        double** LCR=malloc(sizeof(double*)*n);
        double** LP=malloc(sizeof(double*)*n);
        double** LPJ=malloc(sizeof(double*)*n);
        double** LphiBt=malloc(sizeof(double*)*n);
        double** tempq=malloc(sizeof(double*)*n);
        double** LK=malloc(sizeof(double*)*n);
        for(j=0;j<n;j++){
                LA[j]=malloc(sizeof(double)*n);
                Lq[j]=malloc(sizeof(double)*n);
                LQ[j]=malloc(sizeof(double)*n);
                LCR[j]=malloc(sizeof(double)*n);
                LP[j]=malloc(sizeof(double)*n);
                LPJ[j]=malloc(sizeof(double)*n);
                tempq[j]=malloc(sizeof(double)*1);
                LphiBt[j]=malloc(sizeof(double)*n);
                LK[j]=malloc(sizeof(double)*n);
        }
//*********************************QR algorithm***********************************************
int x;
for(x=0;x<2;x++){
int l=0,s,d,h;
for(l=0;l<L;l++){
        for(s=0;s<n;s++){
        MPI_Recv( &(LA[s][0]), n, MPI_DOUBLE,0,10000*l+s, MPI_COMM_WORLD,&status);
        }
        for(s=0;s<n;s++){
        MPI_Recv( &(Lq[s][0]), 1, MPI_DOUBLE,0,10000*l+s, MPI_COMM_WORLD,&status);
        }

//****************QRD*******************************
        for(i=1;i<n;i++){
                for(d=0;d<n;d++){
                        tempq[d][0]=0.0;
                }
                slot=floor(i/core);
                int rem=i-slot*(core-1);
                for(k=rem+(my_PE_num-1)*slot+1;k<=rem+(my_PE_num-1)*slot+slot;k++){
                if(k<=i){
                        nr=0.0;
                for(j=0;j<n;j++){
                        nr+=LA[j][i]*Lq[j][k-1];
                        }
                for(j=0;j<n;j++){
                tempq[j][0]+=nr*Lq[j][k-1];
                }}}

                for(s=0;s<n;s++){
                        MPI_Send( &(tempq[s][0]), 1, MPI_DOUBLE,0, 10000*l+s, MPI_COMM_WORLD);
                }
                for(s=0;s<n;s++){
                MPI_Recv( &(Lq[s][i]), 1, MPI_DOUBLE,0,10000*l+s, MPI_COMM_WORLD,&status);
                }
        }
//****************QRD ends**************************
//***************New A calculation******************

                double slots=floor(n/core);
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                for(j=0;j<n;j++){
                LCR[i][j]=0;
                for(k=0;k<n;k++){
                LCR[i][j]+=Lq[k][i]*LA[k][j];
                }}}
        //C send
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                MPI_Send( &(LCR[i][0]), n, MPI_DOUBLE,0, i, MPI_COMM_WORLD);
                }
        //A cal
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                for(j=0;j<n;j++){
                LP[i][j]=0;
                for(k=0;k<n;k++){
                LP[i][j]+=LCR[i][k]*Lq[k][j];
                }}}
        //A send
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                MPI_Send( &(LP[i][0]), n, MPI_DOUBLE,0, i, MPI_COMM_WORLD);
                }
//*********************Accumulate Q*********************
        if(l==0){
                for(h=0;h<n;h++){
                        for(j=0;j<n;j++){
                        LQ[h][j]=Lq[h][j];
                        }}
        }

        if(l>0){
        //LPJ cal
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                for(j=0;j<n;j++){
                LPJ[i][j]=0;
                for(k=0;k<n;k++){
                LPJ[i][j]+=LQ[i][k]*Lq[k][j];
                }}}

                for(h=0;h<n;h++){
                        for(j=0;j<n;j++){
                        LQ[h][j]=LPJ[h][j];
                        }}
        //LPJ send
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                MPI_Send( &(LQ[i][0]), n, MPI_DOUBLE,0, i, MPI_COMM_WORLD);
                }
        }// l if loop
}//l loop
if(x==0){
//***************New K calculation******************
        for(s=0;s<n;s++){
        MPI_Recv( &(LphiBt[s][0]), n, MPI_DOUBLE,0,10000*l+s, MPI_COMM_WORLD,&status);
        }

        double slots=floor(n/core);
        for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
        for(j=0;j<n;j++){
        LCR[i][j]=0;
        for(k=0;k<n;k++){
        LCR[i][j]+=LphiBt[k][i]*K[k][j];
        }}}
        //C send
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                MPI_Send( &(LCR[i][0]), n, MPI_DOUBLE,0, i, MPI_COMM_WORLD);
                }
        //A cal
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                for(j=0;j<n;j++){
                LP[i][j]=0;
                for(k=0;k<n;k++){
               LP[i][j]+=LCR[i][k]*LphiBt[k][j];
                }}}
        //A send
                for(i=(my_PE_num-1)*slots;i<(my_PE_num)*slots;i++){
                MPI_Send( &(LP[i][0]), n, MPI_DOUBLE,0, i, MPI_COMM_WORLD);
                }

} // x if
} // x loop
} //if loop

end=MPI_Wtime();
printf("\n Time taken by process %d is %f \n",my_PE_num,(end - start)/60);
MPI_Finalize();
return 0;

}

