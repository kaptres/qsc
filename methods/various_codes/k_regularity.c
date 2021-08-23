/*
 *    measure k-regularity (i.e. wiggliness or fractal dimension) of curves
 *    Vasselle & Giraudon, 2-D digital curve analysis: a regularity measure
 *    4th ICCV, 1993, pp. 556-561
 *
 *    Paul Rosin
 *    Joint Research Centre
 *    Ispra, Italy
 *    November 1994
 *    paul.rosin@jrc.it
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef FALSE
# define FALSE 0
# define TRUE (!FALSE)
#endif

#define MAX_PIXELS 100000
#define MAX_SIZE 1500

int no_pixels;
int xdata[MAX_PIXELS],ydata[MAX_PIXELS];

unsigned char image_char[MAX_SIZE][MAX_SIZE];
float image_float[MAX_SIZE][MAX_SIZE];
int width,height;

double max_measure,min_measure;

void analyse(int, int);
double r(int, int, int);
double dist(int, int);
int read_link_data(FILE*, int*);

int main(int argc, char *argv[])
{
    FILE *fp1;
    char file_type[50];
    int j;
    int endoffile;
    int list;
    double measure,factor;
    int x,y;
    int s = 2;
    int k = 3;

    width = height = 256;

    /*
    if (argc == 5) {
        s = atoi(argv[2]);
        k = atoi(argv[4]);
    }
    else if (argc != 3) {
        printf("usage: %s input_file output_file [s k]\n",argv[0]);
        exit(-1);
    }
    */
    if (argc == 4) {
        s = atoi(argv[2]);
        k = atoi(argv[3]);
    }
    else if (argc != 2) {
        printf("usage: %s input_file [s k]\n",argv[0]);
        exit(-1);
    }

    //printf("S = %d;  K = %d\n",s,k);

    if ((fp1=fopen(argv[1],"r")) == NULL) {
        printf("cant open %s\n",argv[1]);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp1,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0) {
        printf("not link data file - aborting\n");
        exit(-1);
    }

    max_measure = 0;
    min_measure = 999;
#if 0
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            image_float[x][y] = 0;
#endif

    do {
        list = read_link_data(fp1,&endoffile);
        analyse(s,k);
    } while (!endoffile);
    fclose(fp1);

    printf("k-irregularity = %f\n",max_measure);

#if 0
    printf("measure = [%4.2f, %4.2f]\n",min_measure,max_measure);
    if (max_measure < 1) max_measure = 1;
    if (min_measure > 0.6) min_measure = 0.6;

    factor = 255.0 / (max_measure - min_measure);

    /* rescale */
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            if (image_float[x][y] > 0)
                /*
                image_char[x][y] = (image_float[x][y] - min_measure) * factor;
                */
                image_char[x][y] = image_float[x][y] * 255;
            else
                image_char[x][y] = 0;

    write_image(image_char,argv[2]);
#endif
}

void analyse(int s, int k)
{
    double r();
    int i;
    double measure = 0;

    if (no_pixels < (k*s+1)) return;

    for (i = 0; i < no_pixels-k*s; i++) {
        measure += r(s,k,i);
    }
    measure /= (double)(no_pixels-k*s);

    if (measure > max_measure)
        max_measure = measure;
    if (measure < min_measure)
        min_measure = measure;

#if 0
    /* store in image */
    for (i = 0; i < no_pixels; i++) {
        image_float[xdata[i]][ydata[i]] = measure;
    }
#endif
}

double r(int s, int k, int i)
{
    int j;
    double d1,d2;
    double dist();

    d1 = dist(i+k*s,i);

    d2 = 0;
    for (j = 1; j <= k; j++)
        d2 += dist(i+j*s,i+(j-1)*s);

    return(d1/d2);
}

double dist(int a, int b)
{
    double dx,dy;

    dx = xdata[a] - xdata[b];
    dy = ydata[a] - ydata[b];
    return(sqrt(dx*dx + dy*dy));
}

int read_link_data(FILE *fp, int *endoffile)
{
    char dumstring[50];
    int j;
    int list;

    fscanf(fp,"%s %d\n",dumstring,&list);
    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&xdata[j],&ydata[j]);
    } while(xdata[j] != -1);
    *endoffile = (ydata[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectlydata terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
    if (no_pixels >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        exit(-1);
    }
    return(list);
}
