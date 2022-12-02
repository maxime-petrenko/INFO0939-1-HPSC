#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {

    


    char tab[]="SalutSalut";

    int lenght = strlen(tab); 

    //printf(" lenght : %d",lenght);


    char delta_str[1000];
    int lt_delta_str = strlen(delta_str); 
    printf(" delta_str : %d",lt_delta_str);

    strcpy(delta_str,tab);

    int lenght_delta_str = strlen(delta_str); 

    printf(" newdelta_str : %d",lenght_delta_str);





    FILE *fp = fopen("param.txt", "r");
    if(fp == NULL) {
        perror("Unable to open file!");
        exit(1);
    }
    

    // Read lines using POSIX function getline
    // This code won't work on Windows
    char *line = NULL;
    size_t len = 0;

    while(getline(&line, &len, fp) != -1) {

        printf(".");
        //printf("line length: %zd\n", strlen(line));

    }

    printf("\n\nMax line size: %zd\n", len);

    fclose(fp);
    free(line);     // getline will resize the input buffer as necessary
                    // the user needs to free the memory when not needed!

    

}



//  gcc main.c -o main -lmn

//  cc input.c input
//  ./input