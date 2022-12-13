/**
 * @file main.c
 * @author oussama, maxime
 * @brief we re trying to solve the accoustic pde using FDTD method
 * @version 0.1
 * @date 2022-10-16
 * chaque 50 timestep il faut créer un fichier output (nombre de fichier output est 81)
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// 1 à quoi nous servenet le xmin et xmax dans les différents fichier binaire on résout le système
//2 c'est quoi la différence entre le domaine et le grid de résolution 
// 3 c'est quoi la taille des matrices p vx vy vz est ce que c'est nx ny nz


//used to organize the values of the input 
const int dim1, dim2, dim3;  /* Global variables, dimension*/

#define VALUES_IN_C(i,j,k) (array_in_c[dim2*dim3*i + dim3*j + k])
#define VALUES_IN_C_INTERPOLED(i,j,k) (array_in_c_interpoled[dim2*dim3*i + dim3*j + k])
#define VALUES_IN_RHO(i,j,k) (array_in_rho[dim2*dim3*i + dim3*j + k])
#define VALUES_IN_RHO_INTERPOLED(i,j,k) (array_in_rho_interpoled[dim2*dim3*i + dim3*j + k])
#define VALUES_OUT_P(i,j,k) (array_out_p[dim2*dim3*i + dim3*j + k])
#define VALUES_OUT_Vx(i,j,k) (array_out_vx[dim2*dim3*i + dim3*j + k])
#define VALUES_OUT_Vy(i,j,k) (array_out_vy[dim2*dim3*i + dim3*j + k])
#define VALUES_OUT_Vz(i,j,k) (array_out_vz[dim2*dim3*i + dim3*j + k])

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



// file handling functions 
bool file_exists(char *filename)
{
    FILE *fp = fopen(filename, "r");
    bool is_exist = false;
    if (fp != NULL)
    {
        is_exist = true;
        fclose(fp); // close the file
    }
    return is_exist;
}

int write_int(int param, char *filename){
    FILE *fp = fopen( filename , "ab" );
    fwrite(&param , 1 , sizeof(int) , fp );
    fclose(fp);
    return(0);
}

int write_double(double param, char *filename){
    FILE *fp = fopen( filename , "ab" );
    fwrite(&param , 1 , sizeof(double) , fp );
    fclose(fp);
    return(0);
}

int rem_if_exists(char *filename){
    if (remove(filename) == 0){
        return 0;
	}
    else{
        return 1;
	}
}
// end of file handling functions
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
// functions for the name of file 
char* concat(const char *s1, const char *s2){
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}
int count(int output_file_nb){
    int digits_count = 0;

    do {
        output_file_nb /= 10;
        ++digits_count;
        
    } while (output_file_nb != 0);
     return digits_count;
}


int main(int argc, char *argv[]){
	/*get and affect the values of the input files */
    /*getting the values of the param file  */
    /*read the parameter ascii file  */
	int idx = 0;
	FILE *ifp;
	char *filename;

    
	// Check if a filename has been specified in the command
	if (argc < 2){
			printf("Missing Filename\n");
			return(1);
	}
	else{
			filename = argv[1];
			printf("Filename : %s\n", filename);
	}
	ifp = fopen(filename,"r");
	if (ifp == NULL){
		printf("Error opening file line 128!\n");
		exit(1);
	}


    double delta,delta_t,max_t;
    int sampling_rate;
    char source_type[50],input_speed_filename[50], input_density_filename[50], output_pressure_base_filename[50],output_velocity_x_base_filename[50],output_velocity_y_base_filename[50], output_velocity_z_base_filename[50];
	fscanf(ifp,"%lf\n%lf\n%lf\n%d\n%s\n%s\n%s\n%s\n%s\n%s\n%s",&delta,&delta_t,&max_t,&sampling_rate,source_type,input_speed_filename,input_density_filename,output_pressure_base_filename,output_velocity_x_base_filename,output_velocity_y_base_filename,output_velocity_z_base_filename);
	fclose(ifp);
	

   
    printf("sampling rate :%d\n",sampling_rate);
    /*get the parameters of the file in_c*/
	FILE *ptr_myfile;
    struct rec my_record_in_c;
    printf("input density filename : %s\n ",input_density_filename);
    ptr_myfile=fopen(input_speed_filename,"rb");
    if (!ptr_myfile){
        printf("Unable to open file input speed filename!");
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
    //printf("nz_in_c : |%d| zmin_in_c :|%lf| zmax_in_c :|%lf|\n",nz_in_c,zmin_in_c,zmax_in_c);
    double * array_in_c = (double *)malloc(nx_in_c*ny_in_c*nz_in_c*sizeof(double)); 
    if(array_in_c == NULL){
        printf("Memory not allocated");
    }else{    
        for (int i = 0; i < nx_in_c; i++) {
            for (int j = 0; j < ny_in_c; j++){
                for(int k = 0; k < nz_in_c; k++){
                    char bytes[8];
                    fread(&bytes,sizeof(double),1,ptr_myfile) ;
                    VALUES_IN_C(i,j,k) = *((double*)bytes);
                    //printf("VALUES_IN_C(%d,%d,%d) : %lf",i,j,k,VALUES_IN_C(i,j,k));
                    //printf("\n");
                }
            }
        }
    }
    fclose(ptr_myfile);
    /*get the parameters of the file in_rho*/
    FILE *ptr_myfile2;
    struct rec my_record_in_rho;
    ptr_myfile2=fopen(input_density_filename,"rb");
    if (!ptr_myfile2){
        printf("Unable to open file input density filename!");
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
    }else{    
        for (int i = 0; i < nx_in_rho; i++) {
            for (int j = 0; j < ny_in_rho; j++){
                for(int k = 0; k < nz_in_rho; k++){
                    char bytes[8];
                    fread(&bytes,sizeof(double),1,ptr_myfile2) ;
                    VALUES_IN_RHO(i,j,k) = *((double*)bytes);
                    //printf("VALUES_IN_RHO(%d,%d,%d) : %lf",i,j,k,VALUES_IN_RHO(i,j,k));
                    //printf("\n");
                }
            }
        }
    }
    fclose(ptr_myfile2);
    /*first the value of c and rho should be interpoled so they can be used in the grid of the resolution*/

    int tailleMatX =(xmax_in_c-xmin_in_c)/delta;
    int tailleMatY =(ymax_in_c-ymin_in_c)/delta;
    int tailleMatZ =(zmax_in_c-zmin_in_c)/delta+1;
    printf("taille mat z : %d",tailleMatZ);


    double * array_in_c_interpoled = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));
    double * array_in_rho_interpoled = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));
    

    for (int i = 0; i <tailleMatX ; i++)
    {
        for (int j = 0; j <tailleMatY ; j++){
            for(int k = 0; k <tailleMatZ ; k++){
                VALUES_IN_C_INTERPOLED(i,j,k)  = interpolate2D(VALUES_IN_C(0,0,0),VALUES_IN_C(0,1,0),VALUES_IN_C(1,0,0),VALUES_IN_C(1,1,0),i,j) ;
                //printf("VALUES_IN_C(%d,%d,%d) : %lf",i,j,k,VALUES_IN_C_INTERPOLED(i,j,k)); 
                VALUES_IN_RHO_INTERPOLED(i,j,k)  = interpolate2D(VALUES_IN_RHO(0,0,0),VALUES_IN_RHO(0,1,0),VALUES_IN_RHO(1,0,0),VALUES_IN_RHO(1,1,0),i,j);
                //later we ll use the interpolate 3d so we can do the 3d version
            }
        }
    }
    
    
    
    /* initialisation des différentes valeur de vélocité et pression pour les différents points à t0 */
    double * array_out_p = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));
    double * array_out_vx = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));
    double * array_out_vy = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));
    double * array_out_vz = (double *)malloc(tailleMatX*tailleMatY*tailleMatZ*sizeof(double));
    //printf("tailleMatriceX est : %d\n",tailleMatX);


    printf("declaration of the arrays done\n");

    if(array_out_p == NULL ||array_out_vx == NULL ||array_out_vy == NULL ||array_out_vz == NULL ){
        printf("Memory not allocated\n");
    }else{    
        for (int i = 0; i <tailleMatX ; i++) {
            for (int j = 0; j <tailleMatY ; j++){
                for(int k = 0; k <tailleMatZ ; k++){
                    VALUES_OUT_P(i,j,k)  = 0.0 ;
                    VALUES_OUT_Vx(i,j,k) = 0.0 ;
                    VALUES_OUT_Vy(i,j,k) = 0.0 ;
                    VALUES_OUT_Vz(i,j,k) = 0.0 ;
                }
            }
        }

    }


	printf("values of the different output files were initialized \n");

	int n;
	int output_file_nb=0;
	for(n=0;n<=max_t/delta_t;n++){


        double freq;
        if(strcmp(source_type,"point_source_middle_3400")==0){
             freq = 3400;
        } 
        //UpdateP();
        double A,B,C;

        printf("x:%d, y:%d, z:%d",tailleMatX, tailleMatY,  tailleMatZ);
        int xmid = tailleMatX/2,ymid = tailleMatY/2,zmid = tailleMatZ/2;  

                
        VALUES_OUT_P(xmid,ymid,zmid)=sin(2*M_PI*n*delta*freq); /// PROBLEM : que des zeros !!!!!!!!!!!!!!!!!!

        printf("VALUES_out_p(%d,%d,%d) : %lf",xmid,ymid,zmid,VALUES_OUT_P(xmid,ymid,zmid));
        printf("\n");

        for (int i = 0; i < tailleMatX ; i++) {

            for (int j = 0; j <= tailleMatY ; j++){


                for(int k = 0; k <tailleMatZ ; k++){

                    


                    // if(i==0){A= VALUES_OUT_Vx(i,j,k); 
                    // }else{A= VALUES_OUT_Vx(i-1,j,k);}
                    // if(j==0){B = VALUES_OUT_Vy(i,j,k); 
                    // }else{B = VALUES_OUT_Vy(i,j-1,k);}
                    // if(k==0){C= VALUES_OUT_Vx(i,j,k); 
                    // }else{C= VALUES_OUT_Vx(i,j,k-1);}


                    //PROBLEM : il ne pas faut pas update la source das la boucle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    VALUES_OUT_P(i,j,k) +=  -( VALUES_IN_RHO(i,j,k) * pow(VALUES_IN_C(i,j,k),2) * delta_t / delta)
                        * ( ( VALUES_OUT_Vx(i,j,k) - A  ) + ( VALUES_OUT_Vy(i,j,k) - B ) + (VALUES_OUT_Vz(i,j,k) - C));  

                    
                }
            }
        }

        
       
        // condition on the center
        
        //condition on the boundaries not done yet TO DO

		//UpdateV();
        double D,E;
        for (int i = 0; i <tailleMatX ; i++) {
            for (int j = 0; j <tailleMatY ; j++){
                for(int k = 0; k <tailleMatZ ; k++){
                    // if(i=tailleMatX-1){D= 0; 
                    // }else{D= VALUES_OUT_P(i+1,j,k);}
                    // if(j=tailleMatY-1){E = 0;
                    // }else{E = VALUES_OUT_P(i,j+1,k);}
                    VALUES_OUT_Vx(i,j,k) += - delta_t / (VALUES_IN_RHO(i,j,k) * delta) * ( D - VALUES_OUT_P(i,j,k) );
                    VALUES_OUT_Vy(i,j,k) += - delta_t / (VALUES_IN_RHO(i,j,k) * delta) * ( E - VALUES_OUT_P(i,j,k) );
                    //VALUES_OUT_Vz(i,j,k) += - delta_t / (VALUES_IN_RHO(i,j,k) * delta) * ( VALUES_OUT_P(i,j,k+1) - VALUES_OUT_P(i,j,k) );
                }
            }
	    }
      
        
        // creating the output files for the pressure values using the sampling rate     
        // change to  100 to sampling rate after
        if (n%sampling_rate == 0){
            char* filename;
            int digits_count = count(output_file_nb);
            printf("digits_count(%d)\n ", digits_count);
            char str_output_file_nb[20];
            sprintf(str_output_file_nb, "%06d", output_file_nb);
            char *filenames[4]= {output_pressure_base_filename,output_velocity_x_base_filename,output_velocity_y_base_filename,output_velocity_z_base_filename};
            for(int l_filenames=0;l_filenames<4;l_filenames++){
            filename = concat(filenames[l_filenames],str_output_file_nb);
            printf("filename written into: %s\n",filename );
            // printf("%s \n",filename)


    


            double zero= 0;
            double max5=99;
			write_int(tailleMatX,filename);
			write_int(tailleMatY,filename);
			write_int(tailleMatZ,filename);
			write_double(zero,filename);
			write_double(max5,filename);
			write_double(zero,filename);
			write_double(max5,filename);
			write_double(zero,filename);
			write_double(max5,filename);
            for (int i = 0; i <tailleMatX ; i++) {
                for (int j = 0; j <tailleMatY ; j++){
                    for(int k = 0; k <tailleMatZ ; k++){
                    switch(l_filenames){
                            case 0:
                                write_double(VALUES_OUT_P(i,j,k),filename);
                                // printf("VALUES_out_p(%d,%d,%d) : %lf",i,j,k,VALUES_OUT_P(i,j,k));
                                // printf("\n");
                                break;
                            case 1:
                                write_double(VALUES_OUT_Vx(i,j,k),filename);
                                break;
                            case 2:
                                write_double(VALUES_OUT_Vy(i,j,k),filename);
                                break;
                            case 3:
                                write_double(VALUES_OUT_Vz(i,j,k),filename);
                                break;
                    }}
			    }

			}
             // deallocate the string                       
            }
            output_file_nb = output_file_nb +1 ;
            free(filename);
        }
	}
	return(0);
}
