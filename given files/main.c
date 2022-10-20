/**
 * @file main.c
 * @author oussama, maxime
 * @brief we re trying to solve the accoustic pde using FDTD method
 * @version 0.1
 * @date 2022-10-16

 */
#include <stdio.h>
#include <math.h>

#define NX 5				/* X [pixels] 100*/
#define NY 5				/*  Y [pixels] 100*/

#define dx 0.01				/*  [m] */
#define dt 20.0e-6			/*  [s] */

#define Nstep 10			/* Nstep*dt=durée de la simulation 0.01 s fait en 10000 pas*/

#define freq 1.0e3			/*  [Hz] */

#define rho 1.3				/* [kg/m^3] */
#define c 340 /*célerité en m *s**(-1)*/



double Vx[NX+1][NY];		/* [m/s] */
double Vy[NX][NY+1];		/*  [m/s] */
double P[NX][NY];			/* [Pa] */


void UpdateV(){

	int i,j;
	for(i=1;i<NX;i++){
		for(j=0;j<NY;j++){
			Vx[i][j] += - dt / (rho * dx) * ( P[i][j] - P[i-1][j] );
		}
	}

	for(i=0;i<NX;i++){
		for(j=1;j<NY;j++){
			Vy[i][j] += - dt / (rho * dx) * ( P[i][j] - P[i][j-1] );
		}
	}
}
void UpdateP(){
	int i,j;
	for(i=1;i<NX;i++){
		for(j=1;j<NY;j++){
			P[i][j] += - ( rho* c * c * dt / dx )
			            * ( ( Vx[i+1][j] - Vx[i][j] ) + ( Vy[i][j+1] - Vy[i][j] ) );
		}
    }
	for(i=0;i<=NX;i++){
		P[i][0] = P[i][1];
		P[i][NY]	= P[i][NY-1];
	}
	for(j=1;j<=NY;j++){
		P[0][j] = P[1][j];
		P[NX][j] = P[NX-1][j];	
	}
}

int main(void){
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
	 	for(int j = 0;j <= NY; j++){
			for(int i=0 ; i<=NX;i++){
				printf("%lf",P[i][j]);
			}
			printf("\n");
		}
		printf("\n\n");
		/* initialisation de la source
		if( n < (1.0/freq)/dt ){
			sig = (1.0-cos(2.0*M_PI*freq*n*dt))/2.0 * sin(2.0*M_PI*freq*n*dt);
		}
		else{
			sig = 0.0;
		}
		P[(int)(NX/2)][(int)(NY/2)] = sig;
        */

	}
	return(0);
}
