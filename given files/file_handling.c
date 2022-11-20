#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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

int write_out_p_int(int param, char *filename){
    FILE *fp = fopen( filename , "ab" );
    fwrite(&param , 1 , sizeof(int) , fp );
    fclose(fp);
    return(0);
}

int write_out_p_double(double param, char *filename){
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

int main(void){
    int Nstep = 20 ;
    int num_file = 0;
    char filename[14];
    for(int i=0;i<Nstep;i++){
        if(i%5==0){
            num_file++;
            snprintf(filename, 14 , "out_p_%d.dat", num_file);
            rem_if_exists(filename);
            write_out_p_int(3,filename);
            write_out_p_int(4,filename);
            write_out_p_int(6,filename);
        }
    }

}