/**
 * @file main.c
 * @author oussama, maxime
 * @brief we re trying to solve the accoustic pde using FDTD method
 * @version 0.1
 * @date 2022-10-16

 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define  NX 5				/* X [pixels] 100*/
#define  NY 5				/*  Y [pixels] 100*/

#define dx 0.01				/*  [m] */
#define dt 20.0e-6			/*  [s] */

#define Nstep 10			/* Nstep*dt=durée de la simulation 0.01 s fait en 10000 pas*/

#define freq 3.4e3			/*  [Hz] */

#define rho 1.225 			/* [kg/m^3] */
#define c 340 /*célerité en m *s**(-1)*/



double Vx[NX+1][NY];		/* [m/s] */
double Vy[NX][NY+1];		/*  [m/s] */
double P[NX][NY];			/* [Pa] */

bool file_exists(const char *filename)
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

FILE *fp;
int write_out_p_int(int param){
   fp = fopen( "out_p.dat" , "ab" );
   fwrite(&param , 1 , sizeof(int) , fp );
   fclose(fp);
   return(0);
}
int write_out_p_double(double param){
   fp = fopen( "out_p.dat" , "ab" );
   fwrite(&param , 1 , sizeof(double) , fp );
   fclose(fp);
   return(0);
}

void UpdateV(){

	int i,j;
	/*we only change the valuex of vx from i=1,j=0 to i=Nx-1,j=Ny-1 so that the values of the velocity 
	 in the boundaries dont get changed and stay equal to 0 same principle is used for Vy*/
	for(i=1;i<NX+1;i++){
		for(j=0;j<NY;j++){
			Vx[i][j] += - dt / (rho * dx) * ( P[i][j] - P[i-1][j] );
		}
	}

	for(i=0;i<NX;i++){
		for(j=1;j<NY+1;j++){
			Vy[i][j] += - dt / (rho * dx) * ( P[i][j] - P[i][j-1] );
		}
	}
}
void UpdateP(){
	int i,j;
	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
			P[i][j] += - ( rho* c * c * dt / dx )
			            * ( ( Vx[i+1][j] - Vx[i][j] ) + ( Vy[i][j+1] - Vy[i][j] ) );
		}
    }
	/*
	for(i=0;i<=NX;i++){
		P[i][0] = P[i][1];
		P[i][NY]	= P[i][NY-1];
	}
	for(j=1;j<=NY;j++){
		P[0][j] = P[1][j];
		P[NX][j] = P[NX-1][j];	
	}
	*/
}

int main(void){
	/*read the parameter ascii file  */

	char line[1000] = "";  // assume each line has at most 999 chars (1 space for NUL character)
	char *lines[1000] = { NULL }; // assume max number of lines is 1000
	int idx = 0;
	FILE *ifp;

	ifp = fopen("param.txt", "r");
	if (ifp == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	while(fgets(line, sizeof(line), ifp) != NULL)
	{
		lines[idx] = strdup(line);
		/*printf("%s", lines[idx]);*/
		idx++;
	}


	char delta_str[1000],delta_t_str[1000],max_t_str[1000],sampling_rate_str[1000],source_type[1000],input_speed_filename[1000], input_density_filename[1000], output_pressure_base_filename[1000],output_velocity_x_base_filename[1000],output_velocity_y_base_filename[1000], output_velocuty_z_base_filename[1000];
	strcpy(delta_str,lines[0]);
	strcpy(delta_t_str,lines[1]);
	strcpy(max_t_str,lines[2]);
	strcpy(sampling_rate_str,lines[3]);
	strcpy(source_type,lines[4]);
	strcpy(input_speed_filename,lines[5]);
	strcpy(input_density_filename,lines[6]);
	strcpy(output_pressure_base_filename,lines[7]);
	strcpy(output_velocity_x_base_filename,lines[8]);
	strcpy(output_velocity_y_base_filename,lines[9]);
	strcpy(output_velocity_y_base_filename,lines[10]);

	double delta =strtod(delta_str, NULL),delta_t =  strtod(delta_t_str,NULL) ,max_t = strtod(max_t_str,NULL);
	int sampling_rate = (int) sampling_rate_str;

	/*writing in the out_p file*/

	char *filename = "out_p.dat";
    if (file_exists(filename)){
    	remove("out_p.dat");
	    printf("The file is deleted successfully.\n");
    } else {
        printf("The file is not deleted.\n");
    }

	write_out_p_int(NX);
    write_out_p_int(NY);
	int NZ=0;
	write_out_p_int(NZ);
	write_out_p_double(0);
	write_out_p_double(NX);
	write_out_p_double(0);
	write_out_p_double(NY);
	write_out_p_double(0);
	write_out_p_double(0);
	int i,j;


	int n;
	/* initialisation des différentes valeur de vélocité et pression pour les différents points à 0 */
	for(i=0;i<NX+1;i++){
		for(j=0;j<NY;j++){
			Vx[i][j] = 0.0;
		}
	}
	for(i=0;i<NX;i++){
		for(j=0;j<NY+1;j++){
			Vy[i][j] = 0.0;
		}
	}
	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
			P[i][j]  = 0.0;
		}
	}
	for(n=0;n<=Nstep;n++){
		//CI
        P[(int)(NX/2)][(int)(NY/2)] = sin(2*M_PI*freq*n*dt);

        UpdateV();
		UpdateP();
	 	for(int j = NY-1 ;j >=0 ; j--){
			for(int i=0 ; i<NX;i++){
				printf("%lf  ,  ",P[i][j]);
			}
			printf("\n");
		}

		for(int j = 0 ; j<NY ;j++){
			for(int i=0 ; i<NX ;i++){
				write_out_p_double(P[i][j]);
			}
		}

	}
	return(0);
}
