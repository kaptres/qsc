/*
 *    calculate fractal dimension of a curve
 *
 *    uses the averaged mass dimension method:
 *            Borkowski W., 1999
 *            Fractal dimension based features are useful descriptors of leaf complexity and shape
 *            Can.J.For.Res., 29, 1301-1310.
 *
 *    probably only works properly if the shape is a closed figure
 *
 *    inputs pixel_float lists (i.e. floating point data)
 *
 *    Paul Rosin
 *    October 2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PIXELS 10000
#define MAX_RULERS 200

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define ABS(x)    (((x)<0.0)? (-(x)): (x))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define SQR(a)    ((a)*(a))

int no_pixels;
double x[MAX_PIXELS], y[MAX_PIXELS];

double dist2();
double kselect();

FILE *fp_out = NULL;

int robust = FALSE;

double dist2(int, int);
void read_link_data(FILE*, int*);
void linear_regression(double*, double*, int);
void fitline2(int, double[], double[], double*, double*);
int solve_line(int, int, double*, double*, double[], double[]);
double euclidean_dist(double, double, double, double, double, double);
double kselect(unsigned long, unsigned long, double[]);
void options(char*);

int main(int argc, char *argv[])
{
    FILE *fp1;
    char file_type[50];
    int i,ii,j,n;
    int endoffile;
    static double distances2[MAX_PIXELS][MAX_PIXELS];
    double count[MAX_RULERS];
    double save_count[MAX_RULERS];
    double radius2,radius[MAX_RULERS],r;
    int start_radius = 1;
    int nbrs;
    int ignore = -1;
    double m,c;

    char *infile = NULL, *outfile = NULL;

    /* parse command line */
    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'i':
                    i++;
                    infile = argv[i];
                    break;
                case 'o':
                    i++;
                    outfile = argv[i];
                    break;
                case 'r':
                    robust = TRUE;
                    break;
                case 'I':
                    i++;
                    ignore = atoi(argv[i]);
                    break;
                default:
                    printf("unknown option %s\n",argv[i]);
                    options(argv[0]);
            }
        }
        else {
            printf("unknown option %s\n",argv[i]);
            options(argv[0]);
        }
    }

    if ((infile == NULL) || (outfile == NULL))
        options(argv[0]);

    if ((fp1 = fopen(infile, "r")) == NULL) {
        printf("cant open %s\n", infile);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp1, "%s\n", file_type);
    j = strcmp(file_type, "pixel_float");
    if (j != 0) {
        printf("not link data file - aborting\n");
        exit(-1);
    }

    if ((fp_out = fopen(outfile, "w")) == NULL) {
        printf("cant open %s\n", outfile);
        exit(-1);
    }

    do {
        read_link_data(fp1, &endoffile);

        fprintf(fp_out,"newgraph\n");
        fprintf(fp_out,"xaxis label fontsize 14 : log (radius)\n");
        fprintf(fp_out,"yaxis label fontsize 14 : log (mean neighbours)\n");
        fprintf(fp_out,"newcurve\n");
        //fprintf(fp_out,"linetype solid\n");
        fprintf(fp_out,"pts\n");

        for (i = 0; i < no_pixels; i++)
            for (j = 0; j < no_pixels; j++)
                distances2[i][j] = dist2(i,j);

        n = 0;
        save_count[0] = 0;

        //for (r = start_radius; save_count[n] < SQR(no_pixels); r *= 2) {
        for (r = start_radius; save_count[n] < SQR(no_pixels); r *= sqrt(2.0)) {
        //for (r = start_radius; (save_count[n] > 5); r *= 2) {
        //for (r = 8; r < no_pixels / 10; r *= sqrt(2.0)) {
        //for (r = 8; r < no_pixels / 10; r += 1) {
            radius2 = SQR(r);

            nbrs = 0;
            for (i = 0; i < no_pixels; i++) {
                for (j = 0; j < no_pixels; j++) {
                    if (distances2[i][j] <= radius2)
                        nbrs++;
                }
            }
            count[n] = (double)nbrs / no_pixels;

            //printf("RAW radius + count: %f  %f\n",r,count[n]);

            save_count[n+1] = nbrs;

            radius[n] = log(r);
            count[n] = log(count[n]);
            fprintf(fp_out,"%f  %f\n",radius[n],count[n]);

            n++;
            if (n > MAX_RULERS) {
                fprintf(stderr,"ERROR: Too many radii\n");
                exit(-1);
            }
        }

        // ignore last points (?)
        n -= 1;

        // shift data down overwritting data at start
        if (ignore > 0) {
            for (i = 0; i < n - ignore; i++) {
                radius[i] = radius[i+ignore];
                count[i] = count[i+ignore];
            }
            n -= ignore;
        }

        if (robust)
            fitline2(n,radius,count,&m,&c);
        else
            linear_regression(radius,count,n);

    } while (!endoffile);
    fclose(fp1);
    fclose(fp_out);
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

    fprintf(fp_out,"title fontsize 16 : fractal dimension = %f\n",m);

    /*
    determine coords of start and end of line
    take x coords of line as truthful
    */
    fprintf(fp_out,"newcurve\n");
    fprintf(fp_out,"marktype none\n");
    fprintf(fp_out,"linetype solid\n");
    fprintf(fp_out,"pts\n");
    fprintf(fp_out,"%f %f\n",x[0],x[0] * m + c);
    fprintf(fp_out,"%f %f\n",x[no_points-1],x[no_points-1] * m + c);
}

/* robust LMedA (absolute not squares) straight line fit */
void fitline2(int no_points, double datax[], double datay[], double *median_m, double *median_c)
{
    int p1,p2;
    int i,j,k;
    double median_error,best_error=9e9;
    static double error[MAX_PIXELS];
    double euclidean_dist();
    double m,c;

    for (i = 0; i < no_points-1; i++) {
        for (j = i+1; j < no_points; j++) {
            if ((datax[i] == datax[j]) && (datay[i] == datay[j]))
                continue;

            /* for each possible pair of points test the line
             * calculate residuals over all the points
             */
            for (k = 0; k < no_points; k++)
                error[k] = euclidean_dist(datax[i],datay[i],
                                          datax[j],datay[j],
                                          datax[k],datay[k]);

            // in case NR in C program starts array at 1
            error[no_points] = error[0];
            median_error = kselect(no_points/2+1,no_points,error);

            if (median_error < best_error) {
                best_error = median_error;
                p1 = i;
                p2 = j;
            }
        }
    }
 
    solve_line(p1,p2,&m,&c,datax,datay);

    *median_m = m;
    *median_c = c;
 
    fprintf(fp_out,"title fontsize 16 : fractal dimension = %f\n",m);
    fprintf(fp_out,"newcurve\n");
    fprintf(fp_out,"marktype none\n");
    fprintf(fp_out,"linetype solid\n");
    fprintf(fp_out,"pts\n");
    fprintf(fp_out,"%f %f\n",datax[0],datax[0] * m + c);
    fprintf(fp_out,"%f %f\n",datax[no_points-1],datax[no_points-1] * m + c);
}

/* equation of line through points indexed by i & j */
int solve_line(int i, int j, double *m, double *c, double datax[], double datay[])
{
    double x1,y1,x2,y2;

    x1 = datax[i];
    y1 = datay[i];
    x2 = datax[j];
    y2 = datay[j];

    if (x2 == x1) {
        return FALSE;
    }

    *c = (x2*y1-x1*y2) / (x2-x1);
    *m = (y2-y1) / (x2-x1);

    return TRUE;
}

/* distance of (x3,y3) from the line (x1,y1) -> (x2,y2) */
double euclidean_dist(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double factor_squared;
    double dist;

    dist = (y2 - y1) * x3 - (x2 - x1) * y3 -
            x1 * y2 + x2 * y1;
    dist = ABS(dist);

    factor_squared = SQR(x1 - x2) + SQR(y1 - y2);

    if (factor_squared == 0)
        dist = 0;
    else
        dist /= sqrt(factor_squared);

    return(dist);
}

/* ========================================================== */
/* select k'th smallest value - from "Numerical recipes in C" */
/* REMEMBER: calling this function alters (reorders) the data */
/* ========================================================== */

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

double kselect(unsigned long k, unsigned long n, double arr[])
{
    unsigned long i,ir,j,l,mid;
    double a,temp;

    /* hack to ensure that there are no problems with starting from 0 or 1
     * PLR February 2007
     * BUT IT IS WRONG IF kselect IS CALLED MORE THAN ONCE!!!
     */
    arr[n] = arr[0];

    l=1;
    ir=n;
    for (;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[ir] < arr[l]) {
                SWAP(arr[l],arr[ir])
            }
            return arr[k];
        } else {
            mid=(l+ir) >> 1;
            SWAP(arr[mid],arr[l+1])
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1],arr[ir])
            }
            if (arr[l] > arr[ir]) {
                SWAP(arr[l],arr[ir])
            }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l])
            }
            i=l+1;
            j=ir;
            a=arr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j])
            }
            arr[l]=arr[j];
            arr[j]=a;
            if (j >= k) ir=j-1;
            if (j <= k) l=i;
        }
    }
}
#undef SWAP

void options(char *progname)
{
    printf("usage: %s [options]\n",progname);
    printf("     -i file   input pixel list\n");
    printf("     -o file   output jgraph\n");
    printf("     -r        robust line fitting\n");
    printf("     -I int    ignore 1st n point during fitting\n");
    exit(-1);
}
