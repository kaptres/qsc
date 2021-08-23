/* 
 * complexity measure:
 * Page, D. L., Koschan, A. F., Sukumar, S. R., Roui-Abidi, B., & Abidi, M. A.
 * Shape analysis algorithm based on information theory
 * Int. Conf. on Image Processing, 2003
 *
 * I use Sturges' rule to get bin size - maybe still not good enough
 *
 * I don't think this method is as good as complexity.c
 *
 * Paul Rosin
 * February 2021
 */

/*
 * measure convexity of a region
 *
 * ratio of area of region to area of convex hull
 * 1 = perfect convexity; 0 = terrible convexity
 *
 * input is a pixel list of boundary coordinates
 *
 * Paul Rosin
 * December 1998
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define MAX_POINTS 100000

#define INTEGER 0

#define PI 3.1415926

#define ABS(x)    (((x)<0.0)? (-(x)): (x))
#define SQR(x)    ((x) * (x))

/* temporary arrays for reading in pixel data */
int no_pixels,no_pixels2;
int x[MAX_POINTS];
int y[MAX_POINTS];
int x2[MAX_POINTS];
int y2[MAX_POINTS];

double hist[1000];
int no_bins;

int xt[MAX_POINTS], yt[MAX_POINTS];

/* file i/o stuff */
char file_type[50];
int file_data_type;
FILE *fp1;

int shift = 1;

void read_data(FILE*, int*);
void read_integer_data(FILE*, int*, int*);
void make_ccw();

int main(int argc, char *argv[])
{
    int endoffile;
    int i,j,jp,jm;
    double dx1,dy1,dx2,dy2;
    double mag1,mag2;
    double dot_product,cross_product;
    double angle;
    int subsample = 4;
    double sum,entropy;

    if (argc == 3) {
        subsample = atoi(argv[2]);
    }
    else if (argc != 2) {
       fprintf(stderr,"%s infile [subsample-interval]\n",argv[0]);
       exit(-1);
    }

    if ((fp1=fopen(argv[1],"r")) == NULL) {
        fprintf(stderr,"cant open %s\n",argv[1]);
        exit(-1);
    }

    /* read magic word for format of file */
    fscanf(fp1,"%s\n",file_type);
    if (strcmp(file_type,"pixel") == 0) {
        file_data_type = INTEGER;
    }
    else {
        fprintf(stderr,"input file not pixel type - aborting\n");
        exit(-1);
    }

    /* read pixel data */
    read_data(fp1,&endoffile);

    if (subsample > 1) {
        no_pixels2 = 0;
        for (i = 0; i < no_pixels; i++) {
            if ((i % subsample) == 0) {
                x[no_pixels2] = x[i];
                y[no_pixels2] = y[i];
                no_pixels2++;
            }
        }
        no_pixels = no_pixels2;
    }

    // make curve clockwise - to make sign of curvature standard for concavity test
    //make_ccw();

    // count number of concave points
    /* calculate subtended angles using dot product AND cross product of vectors */

    // Sturges rule
    no_bins = 1 + ceil(log2(no_pixels));

    for (i = 0; i < no_bins; i++)
        hist[i] = 0;
    for (i = 0; i < no_pixels; i++) {
        jp = (i+shift) % no_pixels;
        jm = (i-shift+no_pixels) % no_pixels;

        dx1 = x[jm] - x[i];
        dy1 = y[jm] - y[i];

        dx2 = x[jp] - x[i];
        dy2 = y[jp] - y[i];

        mag1 = sqrt(SQR(dx1) + SQR(dy1));
        mag2 = sqrt(SQR(dx2) + SQR(dy2));

        // don't need to normalise since the product of lengths will cancel out
        dot_product = (dx1 * dx2 + dy1 * dy2)/(mag1*mag2);
        cross_product = (dx1 * dy2 - dx2 * dy1)/(mag1*mag2);
        angle = atan2(cross_product, dot_product);
        angle *= 180 / PI;
        if (angle < 0) angle += 360;
        //printf("(%d %d)  (%d %d)  (%d %d)\n",x[jm],y[jm],x[i],y[i],x[jp],y[jp]);
        //printf("-> %f\n",angle);
        j = angle * no_bins / 360;
        hist[j]++;
    }

    // normalise histogram
    sum = 0;
    for (i = 0; i < no_bins; i++)
        sum += hist[i];
    for (i = 0; i < no_bins; i++)
        hist[i] /= sum;

    entropy = 0;
    for (i = 0; i < no_bins; i++) {
        //printf("%2d %f\n",i,hist[i]);
        if (hist[i] > 0)
            entropy += hist[i] * log2(hist[i]);
    }
    entropy = -entropy;

    printf("complexity = %f\n",entropy);

    fclose(fp1);
}

void read_data(FILE *fp, int *endoffile)
{
    int list_no;

    read_integer_data(fp,endoffile,&list_no);

    if (no_pixels >= MAX_POINTS) {
        fprintf(stderr,"ERROR: too many pixels\n");
        exit(-1);
    }
}

void read_integer_data(FILE *fp, int* endoffile, int* list_no)
{
    char dumstring[50];
    int j;
    int tx,ty;

    fscanf(fp,"%s %d\n",dumstring,list_no);

    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&tx,&ty);
       x[j] = tx; y[j] = ty;
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;

    /* eliminate duplicated end points */
    if ((x[0] == x[no_pixels-1]) && (y[0] == y[no_pixels-1]))
        no_pixels--;
}

/* ############################################################ */

void make_ccw()
{
    int i,j;
    double area;

    if ((x[0] != x[no_pixels - 1]) && (y[0] != y[no_pixels - 1])) {
        x[no_pixels] = x[0];
        y[no_pixels] = y[0];
        no_pixels++;
    }

    /* determine sense of polygon by calculating "signed" area */
    area = 0;
    for (i = 0; i < no_pixels; i++) {
        j = (i + 1) % no_pixels;
        area += x[i] * y[j];
        area -= y[i] * x[j];
    }

    if ((x[0] == x[no_pixels - 1]) && (y[0] == y[no_pixels - 1]))
        no_pixels--;

    /* reverse list */
    if (area > 0) {
        printf("reversing list\n");
        for (i = 0; i < no_pixels; i++) {
            xt[i] = x[no_pixels - 1 - i];
            yt[i] = y[no_pixels - 1 - i];
        }
        for (i = 0; i < no_pixels; i++) {
            x[i] = xt[i];
            y[i] = yt[i];
        }
    }
}
