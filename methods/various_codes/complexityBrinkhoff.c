/* 
 * modified to produce complexity measure:
 *    Brinkhoff, T., Kriegel, H. P., Schneider, R., & Braun, A.
 *    Measuring the Complexity of Polygonal Objects.
 *    Proc. 3rd ACM International Workshop on Advances in GIS, 1995.
 *
 * seems a bit problematic regarding straight/nearly straight sections & sampling
 *
 * testing on the data from the paper I get the grouping as in figure 8
 * if I subsample the data first (e.g. every 4th point)
 *
 * Paul Rosin
 * January 2021
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

int xt[MAX_POINTS], yt[MAX_POINTS];

/* file i/o stuff */
char file_type[50];
int file_data_type;
FILE *fp1;

void read_data(FILE*, int*);
void read_integer_data(FILE*, int*, int*);
double poly_length(int[], int[], int);
int polyCentroid(int[], int[], int n, double*, double*, double*);
void convex_hull(double*, int, int*, int*, char*, double);
void centroid(double*, int, double*, double*);
void make_rotation_matrix(double, double, double*);
void make_ccw();
int isColinear(int, int, int);
int area(int, int, int);
int isConcave(int, int, int);

int shift = 1;

int main(int argc, char *argv[])
{
    double p[MAX_POINTS*2];
    int n;
    int ch[MAX_POINTS];
    int nch;
    static char directions[] = "norm";
    double epsilon = 0.001;
    int endoffile;
    double l_ch,l_orig;
    double complexity,convexity,amplitude,notches,frequency;
    double xc,yc,region_area,ch_area;
    int flag;
    int i,jp,jm;
    double dx1,dy1,dx2,dy2;
    double mag1,mag2;
    double dot_product,cross_product;
    double angle;
    int count_concavity;
    int subsample = 4;

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
    make_ccw();

    flag = polyCentroid(x,y,no_pixels,&xc,&yc,&region_area);

    for (i = 0; i < no_pixels; i++) {
        p[i*2] = x[i];
        p[i*2+1] = y[i];
    }
    n = no_pixels;

    convex_hull(p,n,ch,&nch,directions,epsilon);

    for (i = 0; i < nch; i++) {
        n = ch[i];
        x2[i] = x[n];
        y2[i] = y[n];
    }

    flag = polyCentroid(x2,y2,nch,&xc,&yc,&ch_area);

    /*
    printf("AREA %f %f\n",region_area,ch_area);
    */
    convexity = region_area / ch_area;
    //printf("area based convexity %f\n",convexity);

    l_ch = poly_length(x2,y2,nch);
    l_orig = poly_length(x,y,no_pixels);
    convexity = l_ch / l_orig;
    //printf("length based convexity %f\n",convexity);
/*
printf("AREA   ORIG  CH%f %f\n",region_area,ch_area);
printf("LENGTH CH   ORIG   %f %f\n",l_ch,l_orig);
*/
    // count number of concave points
    /* calculate subtended angles using dot product AND cross product of vectors */

    count_concavity = 0;
    for (i = 0; i < no_pixels; i++) {
        jp = (i+shift) % no_pixels;
        jm = (i-shift+no_pixels) % no_pixels;

        dx1 = x[jm] - x[i];
        dy1 = y[jm] - y[i];

        dx2 = x[jp] - x[i];
        dy2 = y[jp] - y[i];

        mag1 = sqrt(SQR(dx1) + SQR(dy1));
        mag2 = sqrt(SQR(dx2) + SQR(dy2));

        dot_product = (dx1 * dx2 + dy1 * dy2)/(mag1*mag2);
        cross_product = (dx1 * dy2 - dx2 * dy1)/(mag1*mag2);
        angle = atan2(cross_product, dot_product);
        //printf("(%d %d)  (%d %d)  (%d %d)\n",x[jm],y[jm],x[i],y[i],x[jp],y[jp]);
        //printf("%f\n",angle*180/PI);
        if (angle < 0)
            count_concavity++;
        /*
        if (isConcave(jm,i,jp))
            printf("*** isConcave\n");
        if (isColinear(jm,i,jp))
            printf("ccccc\n");
        */
    }

    notches = (double)count_concavity / (no_pixels - 3);
    //printf("NOTCHES %f   %d %d\n",notches,count_concavity,no_pixels);
    frequency = 16 * pow(notches-0.5,4.0) - 8 * pow(notches-0.5,2.0) + 1;

    amplitude = (l_orig - l_ch) / l_orig;
    convexity = (ch_area - region_area) / region_area;

    complexity = 0.8 * amplitude * frequency + 0.2 * convexity;
    printf("complexity = %f\n",complexity);

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

void read_integer_data(FILE *fp, int *endoffile, int *list_no)
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

double poly_length(int x[], int y[], int n)
{
     int i;
     double dx,dy,sum = 0;

     for (i = 0; i < n; i++) {
         dx = x[i] - x[(i+1)%n];
         dy = y[i] - y[(i+1)%n];
         sum += sqrt(SQR(dx) + SQR(dy));
     }
     return sum;
}

/* ############################################################ */

/*********************************************************************
polyCentroid: Calculates the centroid (xCentroid, yCentroid) and area
of a polygon, given its vertices (x[0], y[0]) ... (x[n-1], y[n-1]). It
is assumed that the contour is closed, i.e., that the vertex following
(x[n-1], y[n-1]) is (x[0], y[0]).  The algebraic sign of the area is
positive for counterclockwise ordering of vertices in x-y plane;
otherwise negative.

Returned values:  0 for normal execution;  1 if the polygon is
degenerate (number of vertices < 3);  and 2 if area = 0 (and the
centroid is undefined).
**********************************************************************/
int polyCentroid(int x[], int y[], int n,
         double *xCentroid, double *yCentroid, double *area)
{
     register int i, j;
     double ai, atmp = 0, xtmp = 0, ytmp = 0;

     if (n < 3) return 1;
     for (i = n-1, j = 0; j < n; i = j, j++) {
         ai = x[i] * y[j] - x[j] * y[i];
         atmp += ai;
         xtmp += (x[j] + x[i]) * ai;
         ytmp += (y[j] + y[i]) * ai;
     }
     *area = atmp / 2;

     /* !!!!! make area unsigned !!!!! */
     *area = ABS(*area);

     if (atmp != 0) {
         *xCentroid = xtmp / (3 * atmp);
         *yCentroid = ytmp / (3 * atmp);
         return 0;
     }
     return 2;
}

/* ############################################################ */

/* 
convexhull.c

Find the convex hull of a set of points.

The algorithm used is as describe in 

Shamos, Michael,  "Problems
in Computational Geometry"  May, 1975 (a set of photocopies -- 
QA 447.S52 1983 in Carlson).  

It originally appeared in a more complicated form in 

Jarvis, R.A., "On the Identification of the Convex Hull of a 
Finite Set of Points in the Plane", Info. Proc. Letters 2(1973),
pp. 18-21.

The algorithm is of complexity k*n, where n is the number of points in the
input and k the number of points in the convex hull.

usage:
convex_hull(p,n,ch,&nch,directions);
where p is an n*2 array of doubles containing the set of points,
n is the number of points,
ch is an array of size n integers to hold the list of points in the 
    convex hull, numbered 0 to nch-1;
In nch the number of points in the convex hull is returned.
directions is either "full" or "norm".  If directions="full" all the
possible points on the convex hull are returned.  If directions="norm"
a minimal set of points to describe the convex hull is returned.

epsilon is the angle tolerance in radians.  When two angles are closer than
PI radians they are considered equal.
*/

/* ### CONTENTS OF INCLUDE FILE INSERTED HERE ### */

/* incl.h */

#define boolean int
#define true 1
#define false 0

#define ERROR -1
#define OK 0

#define EPSILON 1.e-6
#define BIGNUM 1.e20

#define sqr(x) ((x)*(x))

/* ### -------------------------------------- ### */

static int debug=0;

void convex_hull(double *p, int n, int *ch, int *p_nch, char *directions, double epsilon)
{
    double phi, max_phi, dist, max_cen_dist, min_dist, max_dist, 
        x, y, cen_x, cen_y, xp, yp, 
    xx, yy, m[2][2];
    int i, i_keep, vertex, furthest, ch_vert;
    boolean full;

    if (!strcmp(directions,"full")) {
    full = true;
    }
    else if (!strcmp(directions,"norm")) {
    full = false;
    }
    else {
    fprintf(stderr,"convex_hull: invalid argument \"%s\"\n",directions);
    exit(ERROR);
    }
    centroid(p,n,&cen_x,&cen_y);
    /* find the point furthest from the centroid */
    max_cen_dist = 0.;
    for (i=0; i<n; i++) {
    x = p[2*i];
    y = p[2*i+1];
    dist = sqrt(sqr(x-cen_x) + sqr(y-cen_y));
    if (dist>max_cen_dist) {
        max_cen_dist = dist;
        furthest = i;
    }
    }

    /* 
    Determine rotation matrix so that coordinate system for determining
    the orientations of line segments is wrt the line from the point
    under consideration to the centroid.  Then all angles will be 
    < 90 degrees.  To maintain the strict inequality the furthest 
    point along the extreme ray will be used as the next point on 
    the convex hull.
    */

    make_rotation_matrix((cen_x-p[furthest*2])/max_cen_dist,
             (cen_y-p[furthest*2+1])/max_cen_dist,(double *)m);
    ch_vert = 0;
    vertex = furthest;

    /* HACK: to avoid loop with problem data - PLR */
    ch[1] = -999;

    do {
        ch[ch_vert] = vertex;
        /* Find the ray with the maximum and minimum angle in the new 
           coordinate system */
        if (debug) printf("vertex %d\n",vertex);
        max_phi = - BIGNUM;
        min_dist = BIGNUM;
        max_dist = 0.;
        for (i=0; i<n; i++) {
            /* calculate phi, the angle in the new coordinate system */
            x = p[2*i] - p[2*vertex];
            y = p[2*i+1] - p[2*vertex+1];
            xp = m[0][0]*x + m[0][1]*y;
            yp = m[1][0]*x + m[1][1]*y;
            if (debug) printf("\ti %d x %f y %f xp %f yp %f\n", i,x,y,xp,yp);
            if ((dist=sqrt(sqr(x)+sqr(y))) > EPSILON) {
                phi = atan2(yp,xp);
                if (debug) printf("\tphi %f\n", phi);
                if (phi > max_phi-epsilon) {
                    if (full) {
                        /* use the closest point */
                        if (phi > max_phi + epsilon || dist < min_dist) {
                            min_dist = dist;
                            max_phi = phi;
                            i_keep = i;
                        }
                    }
                    else {
                        /* use the furthest point */
                        if (phi > max_phi + epsilon || dist > max_dist) {
                            max_dist = dist;
                            max_phi = phi;
                            i_keep = i;
                        }
                    }
                    if (debug) printf("\t\tmax_phi %f i_keep %d\n", max_phi,i_keep);
                }
            }
        }
        vertex = i_keep;
        xx = cen_x - p[vertex*2];
        yy = cen_y - p[vertex*2+1];
        dist = sqrt(sqr(xx) + sqr(yy));
        make_rotation_matrix(xx/dist,yy/dist,(double *)m);
        ch_vert++;
    /*
    } while (vertex != ch[0]);
    */
    /* HACK: to avoid loop with problem data - PLR */
    } while ((vertex != ch[0]) && (vertex != ch[1]));

    *p_nch = ch_vert;
}

void centroid(double *pts, int n, double *p_cen_x, double *p_cen_y)
/* 
Determines the centroid of the set of n points in the n*2 array of
doubles pts and returns its x-coordinate and y-coordinate in p_cen_x
and p_cen_y respectively.
*/
{
    double sumx, sumy;
    int i;

    sumx = sumy = 0;
    for (i=0; i<n; i++) {
        sumx += pts[2*i];
        sumy += pts[2*i+1];
    }
    *p_cen_x = sumx / n;
    *p_cen_y = sumy / n;
}

void make_rotation_matrix(double cos_theta, double sin_theta, double *m)
/* 
Given a unit vector (vx,vy) system, finds the matrix m that
corresponds to rotating the coordinate system so that its x-axis 
lies in the direction defined by (vx,vy).
*/
{
    *m = cos_theta;
    *(m+1) = sin_theta;
    *(m+2) = - sin_theta;
    *(m+3) = cos_theta;
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
        //printf("reversing list\n");
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

int isColinear(int p1, int p2, int p3)
{
    if (area(p1,p2,p3)==0)
        return TRUE;
    else
        return FALSE;
}

/*
        Returns twice the signed area of the triangle determined by a,b,c,
*/
int area(int a, int b, int c)
{
        int t = 
                x[a] * y[b] - y[a] * x[b] +
                y[a] * x[c] - x[a] * y[c] +
                x[b] * y[c] - x[c] * y[b];

        return ABS(t);
}

int isConcave(int p1, int p2, int p3)
{
    int zcrossproduct;
    int dx1,dy1,dx2,dy2;

    dx1 = x[p2] - x[p1];
    dy1 = y[p2] - y[p1];
    dx2 = x[p3] - x[p2];
    dy2 = y[p3] - y[p2];
    zcrossproduct = dx1*dy2 - dy1*dx2;

    return(zcrossproduct > 0);
}
