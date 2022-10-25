#include<stdio.h>

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

int main(){
    
    write_out_p_int(-1);
    write_out_p_double(7);

}
