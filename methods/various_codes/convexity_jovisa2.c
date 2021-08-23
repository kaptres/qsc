/*
 *    compute CONVEXITY of curve
 *    Jovisa Zunic's method:
 *        minimise ratio of perimeter of axis aligned bounding rectangle to
 *        L1 perimeter of shape
 *    sensitive to noise, so it's best to simplify the curve first
 *
 *    READS LINE SUPERDATA FORMAT (merges multiple lists into one)
 *
 *    I think it can give odd results (> 1) if the shape is not a closed
 *
 *    still numerically goes through all orientations brute-force
 *    (I've updated convexity_jovisa.c)
 *
 *    Paul Rosin
 *    February 2003
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_POINTS 10000

#define PI 3.1415926

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define LISTS 1
#define LINES 2
#define ARCS 3
#define ENDL 4
#define ENDF 5

#define ABS(x)       (((x)<0.0)? (-(x)): (x))
#define SQR(a)       ((a)*(a))
#define MAX(a,b)     (((a) > (b)) ? (a) : (b))
#define MIN(a,b)     (((a) < (b)) ? (a) : (b))

int no_points;
int x[MAX_POINTS],y[MAX_POINTS];
double xo[MAX_POINTS],yo[MAX_POINTS];

double convexity();

FILE *fp_in,*fp_out;
void output_pixels();
double convexity();
void transform(double, double, double);
void read_super_data(int*, int*);
int curve_type(char[]);
void usage(char*);

int main(int argc, char *argv[])
{
    char *infile,*outfile,file_type[50];
    int i,j;
    int endoffile;
    double theta,c;
    double min_c,min_t;
    int close = FALSE;

    infile = outfile = NULL;

    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'c':
                    close = TRUE;
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
                    fprintf(stderr,"ERROR: unknown option %s\n",argv[i]);
                    usage(argv[0]);
                    exit(-1);
            }
        }
        else {
            fprintf(stderr,"ERROR: unknown option %s\n",argv[i]);
            usage(argv[0]);
            exit(-1);
        }
    }

    if (infile == NULL) {
        fprintf(stderr,"need input & output files\n");
        usage(argv[0]);
        exit(-1);
    }

    if ((fp_in=fopen(infile,"r")) == NULL) {
        printf("cant open %s\n",infile);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp_in,"%s\n",file_type);
    j = strcmp(file_type,"super");
    if (j != 0) {
        printf("not super data file - aborting\n");
        exit(-1);
    }

    if (outfile != NULL) {
        if ((fp_out=fopen(outfile,"w")) == NULL) {
            printf("cant open %s\n",outfile);
            exit(-1);
        }
        fprintf(fp_out,"pixel\n");
    }

    /* read in all line data first */
    do {
        read_super_data(&no_points,&endoffile);
    } while (!endoffile);
    fclose(fp_in);

    /* now process it as single list */

    /* add line to close contour */
    if (close) {
        x[no_points] = x[no_points-1];
        y[no_points] = y[no_points-1];
        x[no_points+1] = x[0];
        y[no_points+1] = y[0];
        no_points += 2;
    }

    min_c = 9e9;
    for (theta = 0; theta < PI/2; theta += 0.001) {
        transform(0.0,0.0,theta);
        c = convexity();
        /*
        printf("%f %f\n",theta*180/PI,c);
        */
        if (c < min_c) {
            min_c = c;
            min_t = theta;
        }
    }
    printf("min convexity at theta = %.3f rad (= %.1f deg)\n",min_t,min_t*180.0/PI);

    transform(0.0,0.0,min_t);
    c = min_c;

    if (outfile != NULL)
        output_pixels();

    printf("convexity = %f\n",c);
    /*
    printf("convexity = %.3f\n",c);
    */

    if (outfile != NULL)
        fclose(fp_out);
}

void output_pixels()
{
    int i;

    fprintf(fp_out,"list: 0\n");
    for (i = 0; i < no_points; i++)
        fprintf(fp_out,"%.0f %.0f\n",xo[i],yo[i]);
    fprintf(fp_out,"-1 0\n");
}

double convexity()
{
    int i;
    double dx,dy,px,py;
    double sum_projections;
    double c;
    double min_x,max_x,min_y,max_y;
    double perim_rectangle;

    max_x = min_x = xo[0];
    max_y = min_y = yo[0];
    px = py = 0;
    for (i = 0; i < no_points; i += 2) {
        max_x = MAX(max_x,xo[i]);
        min_x = MIN(min_x,xo[i]);
        max_y = MAX(max_y,yo[i]);
        min_y = MIN(min_y,yo[i]);

        max_x = MAX(max_x,xo[i+1]);
        min_x = MIN(min_x,xo[i+1]);
        max_y = MAX(max_y,yo[i+1]);
        min_y = MIN(min_y,yo[i+1]);

        dx = xo[i+1] - xo[i];
        dy = yo[i+1] - yo[i];
        px += ABS(dx);
        py += ABS(dy);
    }

    sum_projections = px + py;
    perim_rectangle = 2 * (max_x - min_x) + 2 * (max_y - min_y);

    c = perim_rectangle / sum_projections;

    return c;
}

void transform(double xshift, double yshift, double theta)
{
    int i;
    float temp,xt,yt;
    float sine,cosine;

    /* translate */
    for (i = 0; i < no_points; i++) {
        xo[i] = x[i] + xshift;
        yo[i] = y[i] + yshift;
    }

    /* rotate */
    if (theta != 0) {
        sine = sin(-theta);
        cosine = cos(-theta);
        for (i = 0; i < no_points; i++) {
            xt = xo[i];
            yt = yo[i];
            temp = xt*cosine - yt*sine;
            xo[i] = temp;
            temp = xt*sine + yt*cosine;
            yo[i] = temp;
        }
    }
}

/* reads lists and merges them into a single list */
void read_super_data(int *no_segs, int *end_of_file)
{
    int i,j,l,t;
    float f;
    char data_type[40],ch;
    /* equals 0 1st call only */
    static int k = 0;

    /***
    k = 0;
    ***/
    do {
        i = -1;
        do {
            fscanf(fp_in,"%c",&ch);
            data_type[++i] = ch;
        } while (ch != ':');
        data_type[++i] = '\0';
        j = curve_type(data_type);

        if (j == LISTS) {
            fscanf(fp_in,"%d\n",&l);
        }
        else if (j == LINES) {
            fscanf(fp_in,"%f %d %d %d %d\n",
                &f,&x[k],&y[k],&x[k+1],&y[k+1]);
            k += 2;
         if (k > MAX_POINTS) {
            printf("ERROR: too many points&n");
            exit(-1);
         }
        }
        else if (j == ARCS) {
            printf("warning! ARC in data - ignoring\n");
            fscanf(fp_in,"%f %d %d %d %d %d %d %d %d\n",
                &f,&t,&t,&t,&t,&t,&t,&t,&t);
        }
        else if (j == ENDL) {  /* read to end of line */
            fscanf(fp_in,"\n");
        }
    } while ((j != ENDL) && (j != ENDF));
    if (j == ENDF)
        *end_of_file = TRUE;
    else
        *end_of_file = FALSE;
    *no_segs = k;
}

int curve_type(char array[])
{
    int i,j;

    if ((j = strcmp(array,"list:")) == 0)
        i = 1;
    else if ((j = strcmp(array,"line:")) == 0)
        i = 2;
    else if ((j = strcmp(array,"arc:")) == 0)
        i = 3;
    else if ((j = strcmp(array,"endl:")) == 0)
        i = 4;
    else if ((j = strcmp(array,"endf:")) == 0)
        i = 5;
    return(i);
}

void usage(char *progname)
{
    printf("usage: %s -i infile (lines)\n",progname);
    printf("options:\n");
    printf("         -c        close contour - only good for single contour data\n");
    printf("         -o  file  output normalised curve\n");
    exit(-1);
}
