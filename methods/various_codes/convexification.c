/*
 *    convexification of polygon
 *    I have used it to compute a convexity measure
 *    modified from my attempt to approximate a solution to the potato peeling problem
 *    (potato6.c)
 *
 *    each iteration the edge list with max deviation from the CH is chosen for flipping
 *
 *    there seems to be some instability as translating the data gives
 *    different results! (although not too drastic)
 *    - it affects the number of iterations more than the area (i.e. the actual polygon
 *      looks much the same)
 *
 *    Paul Rosin
 *    March 2004
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FALSE 0
#define TRUE (!FALSE)

#define MAX_ITERATIONS 2000

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define ABS(x)    (((x)<0.0)? (-(x)): (x))
#define SQR(x)    ((x)*(x))

#define MAX_PIXELS 10000

int no_pixels;
float x[MAX_PIXELS],y[MAX_PIXELS];       // input data
float xx[MAX_PIXELS],yy[MAX_PIXELS];     // data after one iteration of convexification
float tmpx[MAX_PIXELS],tmpy[MAX_PIXELS]; // store section of polygon pre-reversal */

int ch[MAX_PIXELS];
int nch;

FILE *fp1,*fp2;

int max_iterations = MAX_ITERATIONS;
int reverse = FALSE;
int do_random = FALSE;

void convexification();
void reflect(float, float, float, float, float, float, float*, float*);
void read_link_data(FILE*, int*);
void convex_hull(double*, int, int*, int*, char*, double);
void centroid(double*, int, double*, double*);
void make_rotation_matrix(double, double, double*);
int polyCentroid(float[], float[], int n, double*, double*, double*);
int my_random(int);
void options(char*);

int main(int argc, char *argv[])
{
    char file_type[50];
    int i,j;
    int endoffile;
    int dx,dy;
    char *infile,*outfile;

    infile = outfile = NULL;

    /* parse command line */
    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'm':
                    i++;
                    max_iterations = atoi(argv[i]);
                    break;
                case 'R':
                    do_random = TRUE;
                    break;
                case 'r':
                    reverse = TRUE;
                    break;
                case 'i':
                    i++;
                    infile = argv[i];
                    break;
                case 'o':
                    i++;
                    outfile = argv[i];
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

    if (reverse)
        printf("using REVERSE flip\n");
    else
        printf("using REGULAR flip\n");

    if (do_random)
        printf("using RANDOM flip\n");
    else
        printf("using MAX-DEV flip\n");

    if ((fp1=fopen(infile,"r")) == NULL) {
        fprintf(stderr,"can't open %s\n",infile);
        exit(-1);
    }
    if ((fp2=fopen(outfile,"w")) == NULL) {
        fprintf(stderr,"can't open %s\n",outfile);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp1,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0) {
        printf("not link data file - aborting\n");
        exit(-1);
    }

    fprintf(fp2,"pixel\n");
    do {
        read_link_data(fp1,&endoffile);

        // delete duplicate end point
        dx = x[0] - x[no_pixels-1];
        dy = y[0] - y[no_pixels-1];
        if ((dx == 0) && (dy == 0))
            no_pixels--;

        fprintf(fp2,"list: 0\n");

        convexification();

        if (!endoffile)
            fprintf(fp2,"-1 0\n");
        else
            fprintf(fp2,"-1 -1\n");
    } while (!endoffile);
    fclose(fp1);
    fclose(fp2);
}

/* perform convexification of a polygon */
void convexification()
{
    double xc,yc,area_orig,area_final;
    double p[MAX_PIXELS*2];
    int n1,n2;
    static char directions[] = "norm";
    double epsilon = 0.001;
    int flag;
    float xo,yo;
    int count,count2;
    int i,j;
    double d,max_val,max_dist[MAX_PIXELS];
    int max_pos;
    int no_iterations = 0;
    double save_distances[MAX_PIXELS];
    float xst,yst,xfi,yfi;
    double dx,dy,theta;
    float xm,ym,xn,yn,xoo,yoo;
    int start_count;

    flag = polyCentroid(x,y,no_pixels,&xc,&yc,&area_orig);

    do {
        /* convex hull of figure */
        for (i = 0; i < no_pixels; i++) {
            p[i*2] = x[i];
            p[i*2+1] = y[i];
        }
        n1 = no_pixels;
        convex_hull(p,n1,ch,&nch,directions,epsilon);
    
        /* sort CH indices into ascending order */
        for (i = 0; i < nch-1; i++)
            for (j = i+1; j < nch; j++)
                if (ch[i] > ch[j]) {
                    n1 = ch[i];
                    ch[i] = ch[j];
                    ch[j] = n1;
                }

        /* measure maximum distance for reflection of point set across for each CH edge */
        for (i = 0; i < nch; i++) {
            n1 = ch[i];
            n2 = ch[(i+1)%nch];
            max_dist[i] = 0;
            for (j = (n1+1)%no_pixels; j != n2; j = (j+1)%no_pixels) {
                xst = x[n1]; yst = y[n1];
                xfi = x[n2]; yfi = y[n2];

                reflect(x[j],y[j],xst,yst,xfi,yfi,&xo,&yo);

                d = sqrt( SQR(x[j]-xo) + SQR(y[j]-yo) );
                max_dist[i] = MAX(max_dist[i],d);
            }
        }

        max_pos = 0;
        max_val = max_dist[0];
        for (i = 0; i < nch; i++) {
            if (max_dist[i] > max_val) {
                max_val = max_dist[i];
                max_pos = i;
            }
        }
        //printf("%2d) MAX DIST %f at (%.1f %.1f)\n",no_iterations,max_val,x[ch[max_pos]],y[ch[max_pos]]);
        save_distances[no_iterations] = max_val;

        /* random selection of pocket to flip
         * actually, the CH edge may not be at a pocket, but do it anyway
         */
        if (do_random) max_pos = my_random(nch);

        /* reflect selected points across corresponding convex hull edges */
        count = 0;
        for (i = 0; i < nch; i++) {
            n1 = ch[i];
            n2 = ch[(i+1)%nch];
            xx[count] = x[n1]; yy[count] = y[n1];
            count++;
            /* reflect points */
            count2 = 0;
            if (i == max_pos) {
                start_count = count;
                for (j = (n1+1)%no_pixels; j != n2; j = (j+1)%no_pixels) {
                    xst = x[n1]; yst = y[n1];
                    xfi = x[n2]; yfi = y[n2];
    
                    reflect(x[j],y[j],xst,yst,xfi,yfi,&xo,&yo);
    
                    if (reverse) {
                        /* midpoint of CH edge */
                        xm = (x[n1] + x[n2]) / 2.0;
                        ym = (y[n1] + y[n2]) / 2.0;
        
                        dx = x[n1] - x[n2];
                        dy = y[n1] - y[n2];
                        theta = atan2(-dx,dy);
                        xn = 1000 * cos(theta) + xm;
                        yn = 1000 * sin(theta) + ym;
        
                        /* reflect again to get reverse flip */
                        reflect(xo,yo,xm,ym,xn,yn,&xoo,&yoo);

                        /* store flipped points */
                        tmpx[count2] = xoo; tmpy[count2] = yoo;
                        count2++;
                        xo = xoo; yo = yoo;
                    }
                    xx[count] = xo; yy[count] = yo;
                    count++;
                }

                if (reverse) {
                    /* reinsert flipped points in reverse order */
                    for (j = (n1+1)%no_pixels; j != n2; j = (j+1)%no_pixels) {
                        count2--;
                        xx[start_count] = tmpx[count2]; yy[start_count] = tmpy[count2];
                        start_count++;
                    }
                }
            }
            /* leave points unchanged */
            else {
                for (j = (n1+1)%no_pixels; j != n2; j = (j+1)%no_pixels) {
                    xx[count] = x[j]; yy[count] = y[j];
                    count++;
                }
            }
        }

        /* copy back: ready to flip again */
        for (i = 0; i < count; i++) {
            x[i] = xx[i];
            y[i] = yy[i];
        }

        /***
        for (i = 0; i < count; i++)
            printf("%d %d\n",xx[i],yy[i]);
        printf("%d %d\n",xx[0],yy[0]);
        printf("-1 -1\n");
        ***/

        no_iterations++;
    //} while (max_val > 0);
    } while ((max_val > 0) && (no_iterations < max_iterations));

    /* didn't converge properly
     * so go back to first instance of the current error value was generated
     */
    if (no_iterations == MAX_ITERATIONS) {
        i = no_iterations - 1;
        do {
            i--;
        } while (ABS(save_distances[i]-save_distances[no_iterations-1]) < 0.1);
        no_iterations = i+2;
        //printf("RESETTING iterations to %d\n",no_iterations-1);
    }

    for (i = 0; i < count; i++)
        fprintf(fp2,"%.0f %.0f\n",xx[i],yy[i]);
    fprintf(fp2,"%.0f %.0f\n",xx[0],yy[0]);

    flag = polyCentroid(x,y,no_pixels,&xc,&yc,&area_final);

    /*
    printf("convexification (iterations) %f\n",
        (double)(no_iterations-1)/(double)(SQR(no_pixels) - 4*no_pixels + 1));
    printf("convexification (normalised iterations) %f\n",1.0/((no_iterations-1)+1));
    */
    printf("convexification (iterations) %d\n",no_iterations-1);
    printf("convexification (area) %f\n",area_orig/area_final);
}

/* reflect a point (xi,yi) across the corresponding
 * CH edge (xst,yst)->(xfi,yfi)
 */
void reflect(float xi, float yi, float xst, float yst, float xfi, float yfi, float *xo, float *yo)
{
    double m,offset;
    double dx,dy,xx,yy;

    if ((xst == xfi) && (yst == yfi)) {
        fprintf(stderr,"ERROR: cannot invert point\n");
        exit(-1);
    }
    else if (xst == xfi) {
        offset = xst - xi;
        *xo = (int)(xi + 2 * offset);
        *yo = yi;
    }
    else if (yst == yfi) {
        offset = yst - yi;
        *xo = xi;
        *yo = (int)(yi + 2 * offset);
    }
    else {
        m = (double)(yst-yfi) / (double)(xst-xfi);
        xx = (xi + m*(m*xst + yi - yst))/(1.0 + SQR(m));
        yy = m * (xx - xst) + yst;
        dx = xx - xi; dy = yy - yi;
        *xo = xi + 2 * dx;
        *yo = yi + 2 * dy;
    }
}

void read_link_data(FILE *fp, int* endoffile)
{
    char dumstring[50];
    int j;
    int tx,ty;

    fscanf(fp,"%s %d\n",dumstring,&j);
    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&tx,&ty);
       x[j] = tx;
       y[j] = ty;
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
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
#define PI 3.1415927

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

/*
 * ANSI C code from the article
 * "Centroid of a Polygon"
 * by Gerard Bashein and Paul R. Detmer,
    (gb@locke.hs.washington.edu, pdetmer@u.washington.edu)
 * in "Graphics Gems IV", Academic Press, 1994
 */

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
int polyCentroid(float x[], float y[], int n, double *xCentroid, double *yCentroid, double *area)
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

int my_random(int range)
{
    int i;

    i = random();
    i = i % range;
    return(i);
}

void options(char *progname)
{
    printf("usage: %s [options]\n",progname);
    printf("     -i file   input\n");
    printf("     -o file   output\n");
    printf("     -r        reverse flip\n");
    printf("     -R        random pocket selection\n");
    printf("     -m int    maximum number of iterations\n");
    exit(-1);
}
