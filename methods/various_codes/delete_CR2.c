/*   delete all new line characters and add brackets
 *
 *   Paul Rosin
 *   March 2021
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 50000

int main(int argc, char *argv[])
{
    int l;
    FILE *fp1,*fp2;
    char buf[MAX_LINE];
    int flag = 0;
    
    if (argc != 3) {
        fprintf(stderr,"usage: %s infile outfile\n",argv[0]);
        exit(-1);
    }

    if ((fp1=fopen(argv[1],"r")) == NULL){
        fprintf(stderr,"cant open %s\n",argv[1]);
        exit(-1);
    }

    if ((fp2=fopen(argv[2],"w")) == NULL){
        fprintf(stderr,"cant open %s\n",argv[2]);
        exit(-1);
    }

    while (fgets(buf, MAX_LINE, fp1) != NULL) {
        l = strlen(buf);
            buf[l-1] = '\0';
        //fprintf(fp2,"%s",buf);
        fprintf(fp2,"(%s) ",buf);
    }
	fprintf(fp2,"\n");

    fclose(fp1);
    fclose(fp2);
}

