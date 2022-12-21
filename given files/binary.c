#include<stdio.h>
#include<stdlib.h>


int main(){

   int dim1, dim2, dim3;
   dim1=3;
   dim2=3;
   dim3=1;
   double *** array = (double ***)malloc(dim1*sizeof(double**));
   for (int i = 0; i< dim1; i++) {
   array[i] = (double **) malloc(dim2*sizeof(double *));
      for (int j = 0; j < dim2; j++) {
         array[i][j] = (double *)malloc(dim3*sizeof(double));
      }
   }
  
   
   for (int i = 0; i < 3 ; i++) {
      for (int j = 0; j < 3 ; j++){
            for(int k = 0; k <1 ; k++){
               array[i][j][k] = i+j+k;
              
            }
      }
   }
   printf("%lf ",array[1][1][0]);
   printf("%lf ",array[1][2][0]);
   printf("--------------------ARR-matrice--------------\n");
   for (int i = 0; i < 3 ; i++) {
      for (int j = 0; j < 3 ; j++){
            for(int k = 0; k <1 ; k++){
               printf("%lf ",array[i][j][k]);
         }
      }  
      printf("\n");
   }
}
