/*
 *    calculate fractal dimension of a curve
 *    uses the hybrid method in:
 *        Hayward, J, Orford, JD, Whalley, WB.
 *        Three implementations of fractal 40 analysis of particle outlines.
 *        Computers & Geosciences 1989;15(2):199-207.
 *
 *    inputs pixel_float lists (i.e. floating point data)
 *
 *    Paul Rosin
 *    September 2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PIXELS 200000
#define MAX_RULERS 2000

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define ABS(x)    (((x)<0.0)? (-(x)): (x))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define SQR(a)    ((a)*(a))

int no_pixels;
double x[MAX_PIXELS], y[MAX_PIXELS];

double dist2(int, int);
void read_link_data(FILE*, int*);
void linear_regression(double*, double*, int);

int main(int argc, char *argv[])
{
    FILE *fp1;
    char file_type[50];
    int i,ii,j,n;
    int endoffile;
    double count[MAX_RULERS];
    double save_count[MAX_RULERS];
    double ruler2,ruler[MAX_RULERS],r;
    double dist,err1,err2;
    int plot = FALSE;
    int start_ruler = 4;

    if (argc != 2) {
        printf("usage: %s input_file\n", argv[0]);
        exit(-1);
    }

    if ((fp1 = fopen(argv[1], "r")) == NULL) {
        printf("cant open %s\n", argv[1]);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp1, "%s\n", file_type);
    j = strcmp(file_type, "pixel_float");
    if (j != 0) {
        printf("not link data file - aborting\n");
        exit(-1);
    }

    do {
        read_link_data(fp1, &endoffile);

        // dummy value to ensure large error
        x[no_pixels] = y[no_pixels] = -9999999;

        n = 0;
        save_count[0] = 9999;

        for (r = start_ruler; (save_count[n] > 5); r *= 2) {
        //for (r = 8; r < no_pixels / 10; r *= sqrt(2.0)) {
        //for (r = 8; r < no_pixels / 10; r += 1) {
            ruler2 = SQR(r);
            ii = i = 0;
            count[n] = 0;
            if (plot)
                printf("pixel_float\nlist: 0\n%f %f\n",x[0],y[0]);
            do {
                // find end of ruler
                j = i;
                do {
                    j++;
                } while (dist2(i,j) < ruler2);
                count[n]++;

                // choose the pixel with smallest error wrt ruler length
                err1 = ABS(dist2(i,j-1) - ruler2);
                err2 = ABS(dist2(i,j) - ruler2);
                if (err1 < err2) {
                    dist = sqrt(dist2(i,j-1));
                    i = j-1;
                }
                else {
                    dist = sqrt(dist2(i,j));
                    i = j;
                }

                if (plot) printf("%f %f\n",x[i],y[i]);
                
                if ((j < no_pixels) && (i == ii)) {
                    fprintf(stderr,"stuck in loop - need to start with larger ruler\n");
                    exit(-1);
                }
                else
                    ii = i;
            //} while (i < no_pixels);
            } while (j < no_pixels);

            // modify last step to make it fractional
            count[n] += dist/r - 1;

            //printf("RAW ruler + count: %f  %f\n",r,count[n]);

            save_count[n+1] = count[n];
            count[n] *= r;
    
            ruler[n] = log(r);
            count[n] = log(count[n]);
            printf("%f  %f\n",ruler[n],count[n]);

            n++;
            if (n > MAX_RULERS) {
                fprintf(stderr,"ERROR: Too many rulers\n");
                exit(-1);
            }
        }
        
        linear_regression(ruler,count,n);

    } while (!endoffile);
    fclose(fp1);
}

/* squared straight distance between i & j points of curve */
double dist2(int i, int j)
{
    double dx,dy;

    if (j >= no_pixels)
        return 9e9;
    else {
        dx = x[i] - x[j];
        dy = y[i] - y[j];
        return(SQR(dx) + SQR(dy));
    }
}

void read_link_data(FILE *fp, int *endoffile)
{
    char dumstring[50];
    int j;

    fscanf(fp, "%s %d\n", dumstring, &j);
    j = -1;
    do {
        j++;
        fscanf(fp, "%lf %lf\n", &x[j], &y[j]);
    } while (x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr, "Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
    if (no_pixels >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        exit(-1);
    }
}

/* 
Takes in an array of data points (x,y) and performs a LMS fit to the data.
Minimises the squared error of the y value of each point from the regression line.

from Kreysig 2nd edition 1968
*/
void linear_regression(double *x, double *y, int no_points)
{
    int loop1;
    double m,c;
    double sigma_x,sigma_x_2,sigma_y_2;
    double sigma_y,sigma_xy;
    double temp1,temp2;
    double x_mean,y_mean;

    /* clear variables */
    sigma_x = sigma_x_2 = sigma_y = sigma_y_2 = sigma_xy = 0.0;
    m = c = 0;

    for (loop1 = 0; loop1 < no_points; loop1++) {
        sigma_x += x[loop1];
        sigma_y = sigma_y + y[loop1];
        sigma_x_2 += (x[loop1] * x[loop1]);
        sigma_xy = sigma_xy + (y[loop1] * x[loop1]);
    }
    x_mean = sigma_x / no_points;
    y_mean = sigma_y / no_points;
    temp1 = sigma_x_2 - ((sigma_x * sigma_x)/no_points);
    temp2 = sigma_xy - ((sigma_x * sigma_y)/no_points);
    m = temp2 / temp1;
    c = - x_mean * m + y_mean;

    printf("slope = %f\n",m);
    printf("fractal dimension = %f\n",1-m);

    /* 
    determine coords of start and end of line 
    take x coords of line as truthful
    */
    printf("endpoints of fitted line\n");
    printf("%f %f\n",x[0],x[0] * m + c);
    printf("%f %f\n",x[no_points-1],x[no_points-1] * m + c);
}
