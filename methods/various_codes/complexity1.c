/*
 *    calculate complexity
 *    cumulative residual entropy (see Wang, Vemuri, ICCV 2003) of angles
 *
 *    calculates angles using cross product of vectors (cf. pix_concavities.c)
 *
 *    Paul Rosin
 *    September 2004
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PIXELS 10000

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define PI 3.1415926

#define SQR(a)    ((a)*(a))

int no_pixels;
int x[MAX_PIXELS], y[MAX_PIXELS];

double cre(double*, int);
void sort_double(int, double*);
void read_link_data(FILE*, int*);

int main(int argc, char *argv[])
{
    FILE *fp1;
    char file_type[50];
    int i,j;
    int endoffile;
    double dx1,dy1,dx2,dy2;
    double a[MAX_PIXELS];
    double mag1,mag2;
    int step = 1;
    double tmp;

    if (argc != 3) {
        printf("usage: %s input_file step-size\n",argv[0]);
        return 1;
    }

    step = atoi(argv[2]);

    if ((fp1=fopen(argv[1],"r")) == NULL) {
        printf("cant open %s\n",argv[1]);
        return 1;
    }
    /* read magic word for format of file */
    fscanf(fp1,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0) {
        printf("not link data file - aborting\n");
        return 1;
    }

    read_link_data(fp1,&endoffile);
    fclose(fp1);

    if (step > no_pixels/2)
        return 1;
    
    /* calculate subtended angles using cross product of vectors */
    for (i = 0; i < no_pixels; i++) {
        dx1 = x[(i+step)%no_pixels] - x[i];
        dy1 = y[(i+step)%no_pixels] - y[i];
        dx2 = x[i] - x[(i-step+no_pixels)%no_pixels];
        dy2 = y[i] - y[(i-step+no_pixels)%no_pixels];
        mag1 = sqrt(SQR(dx1) + SQR(dy1));
        mag2 = sqrt(SQR(dx2) + SQR(dy2));
        tmp = (dx1 * dy2 - dx2 * dy1)/(mag1*mag2);
        if (tmp > 1) {
            printf("resetting %f -> 1\n",tmp);
            tmp = 1;
        }
        if (tmp < -1) {
            printf("resetting %f -> -1\n",tmp);
            tmp = -1;
        }
        /* HACK - THIS ONLY HAPPENS IF THERE ARE DUPLICATED POINTS IN THE DATA */
        if ((mag1*mag2) == 0)
            a[i] = a[i-1];
        else
            a[i] = asin(tmp);

        /***
        printf("%3d) %7.1f\n",i,a[i]*180/PI);
        ***/
    }

    printf("STEP %d CRE %f\n",step,cre(a,no_pixels));

    /* plot CDF */
    /*
    printf("newgraph\n");
    printf("newcurve\n");
    printf("linetype solid\n");
    printf("pts\n");
    for (i = 1; i <= no_pixels; i++) {
        printf("%.3f %.3f\n",a[i]*180/PI,(i-1)/(double)(no_pixels-1));
    }
    */
}

double cre(double data[],int n)
{
    int i;
    double I,a,b,w;
    double entropy;

    /* shift data up to start at 1 (not 0) */
    for (i = n; i > 0; i--)
        data[i] = data[i-1];

    sort_double (n, data);

    /*
    for (i = 0; i < n; i++)
        printf("%.3f ",data[i+1]);
    printf("\n");
    */

    entropy = 0;
    for (i = 0; i < n-1; i++) {
        I = i + 1;
        a = (double)i/(n-1);
        b = (double)(i+1)/(n-1);
        w = data[(int)I+1] - data[(int)I];

        /* this is the limit as a->0 */
        if (a == 0)
            entropy -= (b*w*(-1 + 2*log(b)))/4.;
        else
            entropy -=
               (w*(-SQR(a) + SQR(b) + 2*SQR(a)*log(a) - 
               2*SQR(b)*log(b)))/(4.*(a - b));
    }

    return entropy;
}

/* from numerical recipes in C */
void sort_double(int n, double ra[])
{
    int l,j,ir,i;
    double rra;

    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1)
            rra=ra[--l];
        else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                return;
            }
        }
        i=l;
        j=l << 1;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j+1]) ++j;
            if (rra < ra[j]) {
                ra[i]=ra[j];
                j += (i=j);
            }
            else j=ir+1;
        }
        ra[i]=rra;
    }
}

void read_link_data(FILE *fp, int *endoffile)
{
    char dumstring[50];
    int j;

    fscanf(fp,"%s %d\n",dumstring,&j);
    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&x[j],&y[j]);
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
    if (no_pixels >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        return;
    }
    /* eliminate duplicated end points */
    if ((x[0] == x[no_pixels-1]) && (y[0] == y[no_pixels-1]))
        no_pixels--;
}
