/*
 *    simplified version of lines.c
 *    performs recursive subdivision until fixed error threshold on
 *    maximum deviation
 *
 *    effectively assumes open curve
 *
 *    Paul Rosin
 *    February 1996
 *    Paul.Rosin@brunel.ac.uk
 *
 *    delete duplicate endpoint
 *    May 2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef FALSE
# define FALSE 0
# define TRUE (!FALSE)
#endif

#define SIZE 100000    /* used to be 2048 - not enough for racer image */
#define NINETY 1.570796327

#define OK 0
#define NULL_BRKPT -3
#define NULL_SIG -9999

#define SQR(x) ((x)*(x))
/* our version of abs for all number representations */
#define ABS(x) (x)<0?(-(x)):(x)

/* one list of pixels */
float x_c[SIZE],y_c[SIZE];
int flags[SIZE];
int number_pixels;
int list_no;

/* x_c2/y_c2 transformed and rotated */
float x_c2[SIZE],y_c2[SIZE];
int end2;        /* end of segments in x_c2/y_c2 */

FILE *fp, *fp_out, *fp_i;

float min_dev;    /*  minimum deviation for a line */

void chop_pixels(int, int*);
void segment(int, int);
void transform(int, int);
void deviation(int*, float*, int*);
void read_pixels(int*, int*, int*);
void ise();

int main(int argc, char *argv[])
{
    int i,j;
    char *infile,*outfile,*i_file;
    int end_of_file;
    char file_type[50];
    long no_pixels,no_lines_1;
    int printed_endf = FALSE;
    int print_indices = FALSE;
    int chop = 0;
    int num = -1;

    min_dev = 2;
    i_file = infile = outfile = NULL;

    if (argc == 1) {
        printf("            %s\n",argv[0]);
        printf("Usage: line -i file -o file [options]\n");
        printf("options:\n");
        printf(" -d  set deviation threshold (default: %.1f)\n",min_dev);
        printf(" -p  print indices\n");
        printf(" -P  print int\n");
        printf(" -c  chop off some data\n");
        exit(-1);
    }

    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'P':
                    i++;
                    num = atoi(argv[i]);
                    i_file = argv[i];
                    break;
                case 'p':
                    i++;
                    print_indices = TRUE;
                    i_file = argv[i];
                    break;
                case 'i':
                    i++;
                    infile = argv[i];
                    break;
                case 'o':
                    i++;
                    outfile = argv[i];
                    break;
                case 'c':
                    i++;
                    chop = atoi(argv[i]);
                    break;
                case 'd':
                    i++;
                    min_dev = atof(argv[i]);
                    break;
                default:
                    printf("unknown option %s\n",argv[i]);
                    exit(-1);
            }
        }
        else {
            printf("unknown option %s\n",argv[i]);
            exit(-1);
        }
    }

    if ((infile == NULL) || (outfile == NULL)) {
        fprintf(stderr,"need input & output files\n");
        exit(-1);
    }

    /*
    printf("deviation threshold: %.1f\n",min_dev);
    */

    if ((fp = fopen(infile,"r")) == NULL) {
        printf("cant open %s\n",outfile);
        exit(-1);
    }

    /* read magic word for format of file */
    fscanf(fp,"%s\n",file_type);
    if (strcmp(file_type,"pixel") != 0) {
        fprintf(stderr,"not pixel data file - aborting\n");
        exit(-1);
    }

    if ((fp_out = fopen(outfile,"w")) == NULL) {
        fprintf(stderr,"cant open %s\n",outfile);
        exit(-1);
    }

    /* write magic word at top of file */
    fprintf(fp_out,"super\n");

    /* variables for statistics */
    no_pixels = 0;
    no_lines_1 = 0;

    do {
        read_pixels(&number_pixels,&end_of_file,&list_no);

        // added - PLR May 2013 to avoid error caused by initial line from start to end
        // degenerating to a point
        /* eliminate duplicated end points */
        if ((x_c[1] == x_c[number_pixels]) && (y_c[1] == y_c[number_pixels]))
            number_pixels--;

        /* remove some of the data */
        if (chop != 0)
            chop_pixels(chop,&number_pixels);

        if (number_pixels < 2) {
            printf("WARNING: skipping pixel list with only %d pixels\n",number_pixels);
            continue;
        }

        /* count number of pixels */
        no_pixels += number_pixels;

        for (i = 1; i<= number_pixels; i++) {
            flags[i] = NULL_BRKPT;
        }

        segment(1,number_pixels);

        /* count number of lines after first stage */
        for (i=1;i<number_pixels;i++)
            if (flags[i] != NULL_BRKPT)
                no_lines_1++;

        /* calculate integral square error of lines */
        ise();

        /* super data output format */
        fprintf(fp_out,"list: %d\n",list_no);

        fprintf(fp_out,"line: 0 %6.0f %6.0f ",x_c[1],y_c[1]);

        if (print_indices) {
            if ((fp_i = fopen(i_file,"w")) == NULL) {
                printf("cant open %s\n",i_file);
                exit(-1);
            }
            if (chop == 0) {
                if (num != -1)
                    fprintf(fp_i,"1 %d\n",num);
                else
                    fprintf(fp_i,"1 %.2f\n",min_dev);
            }
            else
                fprintf(fp_i,"%d %d\n",1+chop,chop);
        }

        for (j = 2; j <= number_pixels - 1; j++) {
            if (flags[j] != NULL_BRKPT) {
                fprintf(fp_out,"%6.0f %6.0f\n",x_c[j],y_c[j]);
                fprintf(fp_out,"line: 0 %6.0f %6.0f ",x_c[j],y_c[j]);

                if (print_indices) {
                    if (chop == 0) {
                        if (num != -1)
                            fprintf(fp_i,"%d %d\n",j,num);
                        else
                            fprintf(fp_i,"%d %.2f\n",j,min_dev);
                    }
                    else
                        fprintf(fp_i,"%d %d\n",j+chop,chop);
                }
            }
        }
        fprintf(fp_out,"%6.0f %6.0f\n",x_c[number_pixels],y_c[number_pixels]);
        fprintf(fp_out,"endl:\n");

        if (print_indices) {
            if (chop == 0) {
                if (num != -1)
                    fprintf(fp_i,"%d %d\n",number_pixels,num);
                else
                    fprintf(fp_i,"%d %.2f\n",number_pixels,min_dev);
            }
            else
                fprintf(fp_i,"%d %d\n",number_pixels+chop,chop);
            fclose(fp_i);
        }

        if (end_of_file) {
            fprintf(fp_out,"endf:\n");
            printed_endf = TRUE;
        }
    } while (end_of_file == 0);

    if (!printed_endf)
        fprintf(fp_out,"endf:\n");

    fclose(fp_out);

    /* print out statistics */
    /*
    printf("number of pixels processed:   %d\n",no_pixels);
    printf("number of lines:              %d\n",no_lines_1);
    */
}

void chop_pixels(int chop, int *number_pixels)
{
    int i;

    for (i = 1; i <= *number_pixels-chop; i++) {
        x_c[i] = x_c[i+chop];
        y_c[i] = y_c[i+chop];
    }
    *number_pixels -= chop;
}

void segment(int start_in, int finish_in)
{
    int pos;
    float dev;
    int ok;

    transform(start_in,finish_in);
    deviation(&pos,&dev,&ok);
    pos = pos + start_in - 1;
/***
ok = FALSE;
printf("pos = %d\n",pos);
***/
    if (((finish_in - start_in) < 3) || (dev < min_dev) || (ok == FALSE)) {
        /* terminate recursion */
        flags[start_in] = flags[finish_in] = OK;
    }
    else{
        /* recurse to next level down */
        segment(start_in,pos);
        segment(pos,finish_in);
    }
}

void transform(int start, int finish)
{
    int i,j;
    float x_offset,y_offset,x_end,y_end;
    double angle,sine,cosine;
    double temp;

    x_offset = x_c[start];
    y_offset = y_c[start];
    x_end = x_c[finish];
    y_end = y_c[finish];
    if ((x_end - x_offset) == 0.0) {
        if ((y_end - y_offset) > 0.0)
            angle = -NINETY;
        else
            angle = NINETY;
    }
    else{
        temp = ((float)(y_end-y_offset) / (float)(x_end-x_offset));
        angle = -atan(temp);
    }
    cosine = cos(angle);
    sine = sin(angle);
    j = 0;
    for (i=start;i<=finish;i++) {
        j++;
        x_c2[j] = x_c[i] - x_offset;
        y_c2[j] = y_c[i] - y_offset;
        temp = (float)(cosine * x_c2[j]) - (float)(sine * y_c2[j]);
        y_c2[j] = (float)(sine * x_c2[j]) + (float)(cosine * y_c2[j]);
        x_c2[j] = temp;
    }
    end2 = j;
}

void deviation(int *pos, float *dev, int *ok)
{
    int i;
    int pos1;
    float max1,temp;  /* temp used for abs deviation - dont change!! */

    pos1 = 0;
    max1 = 0.0;
    for (i = 1;i <= end2;i++) {
       temp = ABS(y_c2[i]);
       if (temp > max1) {
          max1 = temp;
          pos1 = i;
       }
    }
    /* if no peak found - signal with ok */
    if (max1 == 0.0)
       *ok = FALSE;
    else
        *ok = TRUE;
    *pos = pos1;
    *dev = max1;
}

void read_pixels(int *number_pixels, int *end_of_file, int *list_no)
{
    char dumstring[50];
    int j;
    int tx,ty;

    fscanf(fp,"%s %d\n",dumstring,list_no);
    /* printf("list: %d ",*list_no); */
     /* read integer data */
     j = 0;
     do {
         j++;
         fscanf(fp,"%d %d\n",&tx,&ty);
         x_c[j] = tx;
         y_c[j] = ty;
         if (j > SIZE) {
             printf("ERROR: too many pixels\n");
             exit(-1);
         }
     } while (x_c[j] != -1);
     if (y_c[j] == -1)
            *end_of_file = 1;
     else
            *end_of_file = 0;
     *number_pixels = --j;
     /* printf("pixels: %d\n",*number_pixels); */
}

void ise()
{
    int i,j,st,fi;
    float dist,sum = 0;
    int no_dps,dps[10000];

    dps[0] = 1;
    no_dps = 1;
    for (j = 2; j <= number_pixels - 1; j++) {
        if (flags[j] != NULL_BRKPT)
            dps[no_dps++] = j;
    }
    dps[no_dps++] = number_pixels;

    for (i = 0; i < no_dps-1; i++) {
        st = dps[i]; fi = dps[i+1];
        for (j = st; j < fi; j++) {
            dist = (y_c[st] - y_c[fi]) * x_c[j] -
                   (x_c[st] - x_c[fi]) * y_c[j] -
                    x_c[fi] * y_c[st] +
                    x_c[st] * y_c[fi];
            dist = SQR(dist) /
                 (double)(SQR(x_c[st] - x_c[fi]) + SQR(y_c[st] - y_c[fi]));
            sum += dist;
        }
    }
    //printf("LINES %d ISE = %08.2f\n",no_dps-1,sum);
}
