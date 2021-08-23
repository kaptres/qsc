#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQR(x)       ((x)*(x))

int main(int argc, char *argv[])
{
    float f;
    double sd,mean;
    double sum,sum_sqr;
    int count = 0;

    sum_sqr = sum = 0;
    while (scanf("%f",&f) == 1) {
        sum_sqr += SQR(f);
        sum += f;

        count++;
    }
    mean = sum / count;
    sd = sqrt(sum_sqr/count - SQR(sum)/SQR(count));

    //printf("%d items; MEAN %.4f  SD %.4f  SUM %f\n",count,mean,sd,sum);
    printf("%f\n",mean);
}
