/*
 *    calculate complexity
 *    mean absolute curvature
 *    NOTE: not scale invariant, and so the shape should be normalised first using normalise_shape_mom3
 *
 *    adapted from 3D version:
 *    Quantification of “complexity” in curved surface shape using total absolute curvature✩
 *    Taishi Matsumotoa, Koichiro Satob, Yoshiyuki Matsuokaa, Takeo Katoa
 *    Computers & Graphics, 2019
 *
 *    calculates curvatures using Lowe's code
 *
 *    Paul Rosin
 *    March 2021
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PIXELS 10000
#define MAX_WIN 3000

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define PI 3.1415926

#define SQR(a)    ((a)*(a))
#define ABS(x)    (((x)<0.0)? (-(x)): (x))
#define sqr(a)    ((a)*(a))
#define cube(a)    ((a)*(a)*(a))

int no_pixels;
double x[MAX_PIXELS],y[MAX_PIXELS];
/* smoothed pixel data */
double xsmooth[MAX_PIXELS],ysmooth[MAX_PIXELS];
double curvature[MAX_PIXELS];

double G[MAX_WIN];         /* convolution mask for Gaussian smoothing */
double G1[MAX_WIN];        /* convolution mask for first derivative of Gaussian */
double G2[MAX_WIN];        /* convolution mask for second derivative of Gaussian */

double cor_table[197];

int masksize,halfmask;

void stats(double[], int , double *, double *);
double gaussian(double, int);
double gaussian1(double, int);
double gaussian2(double, int);
void smooth_data (double);
void make_cor_table(double);
void wraparound_data(int);
void G_G2inv_pair(double, double*, double*);
double shrink_correct(double, double, double);
void read_link_data (FILE*, int*);

int main(int argc, char *argv[])
{
    FILE *fp1;
    char file_type[50];
    int i,j;
    int endoffile;
    double sigma;
    double sd,mean;

    if (argc != 3) {
        printf("usage: %s input_file sigma\n",argv[0]);
        exit(-1);
    }

    sigma = atof(argv[2]);

    if ((fp1=fopen(argv[1],"r")) == NULL) {
        printf("cant open %s\n",argv[1]);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp1,"%s\n",file_type);
    j = strcmp(file_type,"pixel_float");
    if (j != 0) {
        printf("not link data file - aborting\n");
        exit(-1);
    }

    read_link_data(fp1,&endoffile);
    fclose(fp1);

    if (sigma > (no_pixels / 5.0)) {
        fprintf(stderr,"ERROR: not enough pixels for requested sigma\n");
        exit(-1);
    }

    smooth_data(sigma);

    /*
    for (i = 0; i < no_pixels-masksize; i++)
        printf("%f %f\n",xsmooth[i+halfmask],ysmooth[i+halfmask]);
    exit(-1);
    */

    /* shift values back to original position (starting at 0) */
    for (i = 0; i < no_pixels-masksize; i++)
        curvature[i] = ABS(curvature[i+halfmask]);
    no_pixels -= masksize;

    stats(curvature,no_pixels,&sd,&mean);

    printf("SIGMA %.3f mean absolute curvature %f\n",sigma,mean);
}

void stats(double data[], int n, double *sd, double *mean)
{
    double sum,sum_sqr;
    int i;

    /* calculate mean & sd */
    sum = sum_sqr = 0;
    for (i = 0; i < n; i++) {
        sum += data[i];
        sum_sqr += SQR(data[i]);
    }
    *mean = sum / n;
    *sd = sqrt(sum_sqr/n - SQR(sum)/SQR(n));
}


double gaussian(double sigma, int t)
{
    double g;

    g = exp(-(t*t)/(2*sigma*sigma)) / (sigma * sqrt(2*PI));

    return(g);
}

/* first derivative */
double gaussian1(double sigma, int t)
{
    double g;

    g = exp(-(t*t)/(2*sigma*sigma)) * (-t)
        / (sigma * sigma * sigma * sqrt(2*PI));

    return(g);
}

/* second derivative */
double gaussian2(double sigma, int t)
{
    double g;

    g = exp(-(t*t)/(2*sigma*sigma)) * (t*t/(sigma*sigma)-1)
        / (sigma * sigma * sigma * sqrt(2*PI));

    return(g);
}

void smooth_data(double sigma)
{
    int j,t;
    double sum = 0.0;
    double sum_p = 0.0;
    double sum_n = 0.0;
    int point;
    double xG,yG,xG1,yG1,xG2,yG2;
    int xpoint,ypoint;
    double dist;
    double shrink_correct();

    masksize = 5 * sigma;
    if ((masksize % 2) == 0) masksize++;
    if (masksize >= MAX_WIN) masksize = MAX_WIN-1;
    if (masksize >= 2*no_pixels) {
        masksize = 2*no_pixels-1;
        printf("truncated masksize due to lack of pixel data!!\n");
    }
    halfmask = masksize / 2;

    wraparound_data(halfmask);

    /* generate masks for Gaussians & derivatives */
    for (t = -halfmask; t <= halfmask; t++)
        G[halfmask+t] = gaussian(sigma,t);
    for (t = 0; t < masksize; t++) sum += G[t];
    for (t = 0; t < masksize; t++) G[t] /= sum;

    for (t = -halfmask; t <= halfmask; t++)
        G1[halfmask+t] = gaussian1(sigma,t);

    for (t = -halfmask; t <= halfmask; t++)
        G2[halfmask+t] = gaussian2(sigma,t);

    for (t = 0; t < masksize; t++)
        if (G2[t] > 0)
            sum_p += G2[t];
        else
            sum_n -= G2[t];
    for (t = 0; t < masksize; t++)
        if (G2[t] < 0)
            G2[t] = G2[t]*sum_p/sum_n;

    /* generate shrinkage correction table */
    make_cor_table(sigma);

    for (t = 0; t < no_pixels-masksize; t++) {
        xG = yG = xG1 = yG1 = xG2 = yG2 = 0;
        for(j = 0; j < masksize; j++) {
            xpoint = x[t+j];
            ypoint = y[t+j];
            xG += (G[j]*xpoint);
            yG += (G[j]*ypoint);
            xG1 += (G1[j]*xpoint);
            yG1 += (G1[j]*ypoint);
            xG2 += (G2[j]*xpoint);
            yG2 += (G2[j]*ypoint);
        }
        dist = sqrt(sqr(xG1)+sqr(yG1));
        point = t + halfmask;
        xsmooth[point] = shrink_correct(sigma,xG,xG2);
        ysmooth[point] = shrink_correct(sigma,yG,yG2);
        curvature[point] = ((xG1*yG2)-(yG1*xG2))/cube(dist);
    }
}

void wraparound_data(int length)
{
    int s,t;

/***
    if (length > no_pixels) {
        printf("length bigger than no_pixels!!!!\n");
    }
***/

    /* shift up the pixel data */
    for (t = no_pixels-1; t >= 0; t--) {
       x[t+length] = x[t];
       y[t+length] = y[t];
    }

    /* copy data to start of pixel list */
    no_pixels += length;
    for (s = length-1,t = no_pixels-1; t > no_pixels-1-length; s--,t--) {
       x[s] = x[t];
       y[s] = y[t];
    }

    /* copy data to end of pixel list */
    for (t = length; t < length*2; t++) {
       x[no_pixels] = x[t];
       y[no_pixels] = y[t];
       no_pixels++;
    }
}

/* ##### ALL CODE BELOW IS TRANSLATED FROM DAVID LOWE'S LISP VERSION ##### */

/*********************************************************************
Produce a table giving shrinkage correction values as a function of
the second derivative convolution.
The table goes from 1/G2 = sigma to 1/G2 = 50*sigma in intervals of
0.25*sigma (197 values).  Good results could be achieved with a smaller
table, but play it safe.
*********************************************************************/
void make_cor_table(double sigma)
{
    double curRadius = 50.0;
    double curG,curG2i;
    double prevG,prevG2i;
    double desG2i;
    int i;
    int first=TRUE;

    G_G2inv_pair(sigma*curRadius,&curG,&curG2i);
    for (i=0;i<197;i++){
        desG2i=sigma*(50-i*0.25);
        /* lower radius until we bracket desired value of G2 */
        while (first || (curG2i > desG2i)) {
            first = FALSE;
            prevG = curG;
            prevG2i = curG2i;
            curRadius -= 0.25;
            G_G2inv_pair(sigma*curRadius,&curG,&curG2i);
            /*  if G2i stops shrinking then fill rest of table
                with current value of G */
            if (curG2i>prevG2i) {
                curG2i = 0;
                curG = prevG;
            }
        }
        /* place interpolated shrinkage error in correction table */
        cor_table[196-i] = curG +
            ((desG2i-curG2i)*(prevG-curG))/(prevG2i-curG2i);
    }
}

/*********************************************************************
We compute the amount of shrinkage that occurs for a given value of
the G2 convolution by actually convolving the current masks with a
path-length parameterized circle of the given radius.  While this is
more time-consuming than simply calculating the formulas given in
my paper, it is more accurate because it compensates for the sampling
and truncation errors in the current convolution masks (these errors
can amount to 10% or more for the current mask size of 5*sigma).
Since this is done during precomputation, there is little need to worry
about efficiency.
*********************************************************************/
void G_G2inv_pair(double radius, double* curG, double* curG2i)
{
    double s;    /* current path length parameter value */
    double x;    /* corresponding x coordinate value */
    double g,g2;
    int i;

    g = g2 = 0;
    for (i=0;i<=masksize;i++){
        s = i - masksize/2;
        x = radius * (1 - cos(s/radius));
        g += (x*G[i]);
        g2 += (x*G2[i]);
    }
    *curG = g;
    *curG2i = 1.0/g2;
}

/*********************************************************************
Correct the Gaussian convolution value g according to the measure of
curvature g2.  Correction is computed using table lookup and linear
interpolation.
Efficiency is very important here, because this is executed for each
point that is smoothed.
*********************************************************************/
double shrink_correct(double sigma, double g, double g2)
{
    double ratio,index,sign,frac;
    int int1;

    if (fabs(g2*sigma) < 0.02) return(g);
    ratio = fabs(1.0/(g2*sigma));
    index = 4*(ratio-1);
    if (index > 195) {
        index = 195;
    }
    if (index < 0) {
        index = 0;
    }
    if (g2>0) sign = 1; else sign = -1;
    int1 = (int)index;
    frac = index - int1;
    return(g-sign*((1-frac)*(cor_table[int1])+(frac*cor_table[1+int1])));
}


void read_link_data(FILE *fp, int *endoffile)
{
    char dumstring[50];
    int j;

    fscanf(fp,"%s %d\n",dumstring,&j);
    j = -1;
    do{
       j++;
       fscanf(fp,"%lf %lf\n",&x[j],&y[j]);
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
    if (no_pixels >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        exit(-1);
    }
    /* eliminate duplicated end points */
    if ((x[0] == x[no_pixels-1]) && (y[0] == y[no_pixels-1]))
        no_pixels--;
}
