#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


const int dim1, dim2, dim3;  /* Global variables, dimension*/
#define VALUES_IN_C_INTERPOLED(i,j,k) (array_in_c_interpoled[dim2*dim3*i + dim3*j + k])
#define VALUES_IN_RHO_INTERPOLED(i,j,k) (array_in_rho_interpoled[dim2*dim3*i + dim3*j + k])


float interpolate1D(float v1, float v2, float x){
    return v1*(1-x) + v2*x;
}

float interpolate2D(float v1, float v2, float v3, float v4, float x, float y){
    float s = interpolate1D(v1, v2, x);
    float t = interpolate1D(v3, v4, x);
    return interpolate1D(s, t, y);
}
float interpolate3D(float v1, float v2, float v3, float v4,float v5, float v6, float v7, float v8, float x, float y, float z){
    float s = interpolate2D(v1, v2, v3, v4, x, y);
    float t = interpolate2D(v5, v6, v7, v8, x, z);
    return interpolate1D(s, t, z);
}

int main(int argc, char *argv[]){
    int tailleMatX = 5;
    int tailleMatY = 5;
    int tailleMatZ = 5;

    double * array_in_c_interpoled = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));
    double * array_in_rho_interpoled = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));

    for (int i = 0; i <tailleMatX ; i++) {
            for (int j = 0; j <tailleMatY ; j++){
                for(int k = 0; k <tailleMatZ ; k++){
                    VALUES_IN_C_INTERPOLED(i,j,k) = interpolate3D(340,340,340,340,340,340,340,340,i,j,k); 
                    printf("VALUES_IN_C(%d,%d,%d) : %lf",i,j,k,VALUES_IN_C_INTERPOLED(i,j,k)); 
                }
            }
    }
}