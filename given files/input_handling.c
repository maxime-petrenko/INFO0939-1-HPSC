#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


const int dim1, dim2, dim3;  /* Global variables, dimension*/

#define VALUES_IN_C(i,j,k) (array_in_c[dim2*dim3*i + dim3*j + k])
#define VALUES_IN_RHO(i,j,k) (array_in_rho[dim2*dim3*i + dim3*j + k])

/* Our structure for the first terms of the file */
struct rec{
    int nx;
    int ny;
    int nz;
    char xmin[8];
    char xmax[8];
    char ymin[8];
    char ymax[8];
    char zmin[8];
    char zmax[8];
};

int main()
{
    FILE *ptr_myfile;
    struct rec my_record_in_c;
    ptr_myfile=fopen("in_c.dat","rb");
    if (!ptr_myfile){
        printf("Unable to open file!");
        return 1;
    }
    fread(&my_record_in_c,sizeof(struct rec),1,ptr_myfile);

    int nx_in_c = my_record_in_c.nx;
    int ny_in_c = my_record_in_c.ny;
    int nz_in_c = my_record_in_c.nz;
    double xmin_in_c = *((double*)my_record_in_c.xmin);
    double xmax_in_c = *((double*)my_record_in_c.xmax);
    double ymin_in_c = *((double*)my_record_in_c.ymin);
    double ymax_in_c = *((double*)my_record_in_c.ymax);
    double zmin_in_c = *((double*)my_record_in_c.zmin);
    double zmax_in_c = *((double*)my_record_in_c.zmax);

   
    double * array_in_c = (double *)malloc(nx_in_c*ny_in_c*nz_in_c*sizeof(double)); 
    if(array_in_c == NULL){
        printf("Memory not allocated");
        return 2;
    }else{    
        for (int i = 0; i < nz_in_c; i++) {
            for (int j = 0; j < ny_in_c; j++){
                for(int k = 0; k < nx_in_c; k++){
                    char bytes[8];
                    fread(&bytes,sizeof(double),1,ptr_myfile) ;
                    VALUES_IN_C(k,j,i) = *((double*)bytes);
                    printf("VALUES_IN_C(%d,%d,%d) : %lf",k,j,i,VALUES_IN_C(k,j,i));
                    printf("\n");
                }
            }
        }
    }
    fclose(ptr_myfile);
    
    FILE *ptr_myfile2;
    struct rec my_record_in_rho;
    ptr_myfile2=fopen("in_rho.dat","rb");
    if (!ptr_myfile2){
        printf("Unable to open file!");
        return 1;
    }
    fread(&my_record_in_rho,sizeof(struct rec),1,ptr_myfile2);

    int nx_in_rho = my_record_in_rho.nx;
    int ny_in_rho = my_record_in_rho.ny;
    int nz_in_rho = my_record_in_rho.nz;
    double xmin_in_rho = *((double*)my_record_in_rho.xmin);
    double xmax_in_rho = *((double*)my_record_in_rho.xmax);
    double ymin_in_rho = *((double*)my_record_in_rho.ymin);
    double ymax_in_rho = *((double*)my_record_in_rho.ymax);
    double zmin_in_rho = *((double*)my_record_in_rho.zmin);
    double zmax_in_rho = *((double*)my_record_in_rho.zmax);
    double * array_in_rho = (double *)malloc(nx_in_rho*ny_in_rho*nz_in_rho*sizeof(double)); 
    if(array_in_rho == NULL){
        printf("Memory not allocated");
        return 2;
    }else{    
        for (int i = 0; i < nz_in_rho; i++) {
            for (int j = 0; j < ny_in_rho; j++){
                for(int k = 0; k < nx_in_rho; k++){
                    char bytes[8];
                    fread(&bytes,sizeof(double),1,ptr_myfile2) ;
                    VALUES_IN_RHO(k,j,i) = *((double*)bytes);
                    printf("VALUES_IN_RHO(%d,%d,%d) : %lf",k,j,i,VALUES_IN_RHO(k,j,i));
                    printf("\n");
                }
            }
        }
    }
    fclose(ptr_myfile2);
    return 0;
}
