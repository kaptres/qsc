/*
 *    normalise pixel region using moments (area, centroid & orientation)
 *    also normalise the number of pixel samples
 *
 *    Paul Rosin
 *    March 2021
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PIXELS 200000

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define ABS(x)    (((x)<0.0)? (-(x)): (x))
#define SQR(a)    ((a)*(a))
#define MAX(a,b)     (((a) > (b)) ? (a) : (b))
#define MIN(a,b)     (((a) < (b)) ? (a) : (b))

double area = 40000;
int no_pixels1;
int no_pixels2 = 1000;
int x[MAX_PIXELS], y[MAX_PIXELS];
double xdata1[MAX_PIXELS];
double ydata1[MAX_PIXELS];
double xdata2[MAX_PIXELS];
double ydata2[MAX_PIXELS];

int use_gap = TRUE;
double gap;

void read_link_data(FILE*, int*);
void store_link_data(FILE*, int);
double distance(double, double, double, double);
void interp();
void options(char*);

int main(int argc, char *argv[])
{
    FILE *fp_in,*fp_out;
    char file_type[50];
    int i, j;
    int endoffile;
    double a,dx,dy;
    double m00,m11,m10,m01,m20,m02;
    double xc,yc,tx,ty,xx,yy;
    double cosine,sine,theta,scale;

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
                case 'a':
                    i++;
                    area = atof(argv[i]);
                    break;
                case 'g':
                    i++;
                    gap = atof(argv[i]);
                    break;
                case 'n':
                    i++;
                    no_pixels2 = atof(argv[i]);
                    use_gap = FALSE;
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

    if ((fp_in = fopen(infile, "r")) == NULL) {
        printf("cant open %s\n", argv[1]);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp_in, "%s\n", file_type);
    j = strcmp(file_type, "pixel");
    if (j != 0) {
        printf("not link data file - aborting\n");
        exit(-1);
    }

    if ((fp_out = fopen(outfile, "w")) == NULL) {
        printf("cant open %s\n", outfile);
        exit(-1);
    }

   fprintf(fp_out,"pixel_float\n");
   do {
        read_link_data(fp_in, &endoffile);

        m00 = m10 = m01 = m11 = m20 = m02 = 0;
        for (i = 0; i < no_pixels1; i++) {
            j = (i-1+no_pixels1) % no_pixels1;

            dx = x[i] - x[j];
            dy = y[i] - y[j];
            a = x[i]*dy - y[i]*dx;

            m00 += a;
            m10 += a * (x[i] - dx/2);
            m01 += a * (y[i] - dy/2);
            m11 += a * (x[i]*y[i] - x[i]*dy/2 - y[i]*dx/2 + dx*dy/3);
            m20 += a * (x[i]*x[i] - x[i]*dx + dx*dx/3);
            m02 += a * (y[i]*y[i] - y[i]*dy + dy*dy/3);
        }
        m00 /= 2.0;
        m10 /= 3.0;
        m01 /= 3.0;
        m11 /= 4.0;
        m20 /= 4.0;
        m02 /= 4.0;

        xc = m10/m00;
        yc = m01/m00;

        /*
        printf("m00 %f\n",m00);
        printf("m10 %f\n",m10);
        printf("m01 %f\n",m01);
        printf("m11 %f\n",m11);
        printf("m20 %f\n",m20);
        printf("m02 %f\n",m02);
        printf("xc %f\n",m10/m00);
        printf("yc %f\n",m01/m00);
        */

        theta = atan2(2*(m00*m11-m10*m01) , (m00*m20 - m10*m10)-(m00*m02 - m01*m01) ) / 2;
        //printf("theta %f\n",theta*58);

        scale = area / ABS(m00);
        /* since scaling both x & y */
        scale = sqrt(scale);
        //printf("S %f\n",scale);

        cosine = cos(-theta);
        sine = sin(-theta);

        for (i = 0; i < no_pixels1; i++) {
            tx = x[i] - xc;
            ty = y[i] - yc;
            xx = cosine*tx - sine*ty;
            yy = sine*tx + cosine*ty;
            // hack - interpolation code starts indexing at 1
            xdata1[i+1] = xx * scale + 500;
            ydata1[i+1] = yy * scale + 500;
        }

        interp();
        store_link_data(fp_out,endoffile);
    } while (!endoffile);
    fclose(fp_in);
    fclose(fp_out);
}

void read_link_data(FILE *fp, int *endoffile)
{
    char dumstring[50];
    int j;

    fscanf(fp, "%s %d\n", dumstring, &j);
    j = -1;
    do {
        j++;
        fscanf(fp, "%d %d\n", &x[j], &y[j]);
    } while (x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr, "Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels1 = j;
    if (no_pixels1 >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        exit(-1);
    }
}

void store_link_data(FILE *fp, int endoffile)
{
    int j;

    fprintf(fp,"list: 0\n");
    // HACK - last value sometimes ends up as zero! - PLR
    //for (j = 1; j <= no_pixels2; j++)
    for (j = 1; j < no_pixels2; j++)
       fprintf(fp,"%f %f\n",xdata2[j],ydata2[j]);
    if (endoffile)
        fprintf(fp,"-1 -1\n");
    else
        fprintf(fp,"-1 0\n");
}

double distance(double x1,double y1,double x2,double y2)
{
    double t1,t2;

    t1 = SQR(x1-x2);
    t2 = SQR(y1-y2);
    return sqrt(t1+t2);
}

void interp()
{
    double position,inc,length;
    double temp1,temp2;
    int i,j;
    int no_out;
    double leng[10000];
    
    /* determine the total distance along the curve */
    length = 0.0;
    leng[1] = 0.0;
    for (i=2;i<=no_pixels1;i++) {
        length  = length + distance(xdata1[i],ydata1[i],xdata1[i-1],ydata1[i-1]);
        leng[i] = length;
        /* printf("length: %6.3f leng[%d]: %6.3f\n",length,i,leng[i]); */
    }
    if (use_gap)
        no_pixels2 = length / gap;

    // hack? - PLR
    //inc = length / no_pixels2;
    inc = length / (no_pixels2-1);

    printf("number of new pixels: %6d\n",no_pixels2);
    printf("new pixel spacing %f\n",inc);
    xdata2[1] = xdata1[1]; ydata2[1] = ydata1[1];
    no_out = 1;
    position = 0.0;
    for (i=1;i<= no_pixels2;i++) {
        position += inc;
        j = 0;
        do {
           j++;
           if ((leng[j] <= position) && (position < leng[j+1])) {
               no_out++;
               /* printf("no: %6d 1st: %6.3f pos: %6.3f 2nd: %6.3f\n",
                      no_out,leng[j],position,leng[j+1]); */
               /* determine accurate coords of interpolated point */
               temp1 = distance(xdata1[j],ydata1[j],xdata1[j+1],ydata1[j+1]);
               temp2 = position - leng[j];
               xdata2[no_out] = xdata1[j]+(xdata1[j+1]-xdata1[j])*temp2/temp1;
               ydata2[no_out] = ydata1[j]+(ydata1[j+1]-ydata1[j])*temp2/temp1;
           }
        }while(j < no_pixels1);
    }
    // hack - PLR
    xdata2[no_pixels2] = xdata1[no_pixels1]; ydata2[no_pixels2] = ydata1[no_pixels1];
}

void options(char *progname)
{
    printf("usage: %s [options]\n",progname);
    printf("     -i file   input pixel list\n");
    printf("     -o file   output normalised pixel list\n");
    printf("     -n int    number of output pixels (default: %d)\n",no_pixels2);
    printf("     -g float  intersample spacing of output pixels (default: instead, use number of output pixels)\n");
    printf("     -a float  area (default: %.3f)\n",area);
    exit(-1);
}
