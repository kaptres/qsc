// computes average of each set of 12 numbers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQR(x)       ((x)*(x))

int main(int argc, char *argv[])
{
    float f;
    double mean;
    double sum;
    int count = 0;
	int i;

    sum = 0;
    while (scanf("%f",&f) == 1) {

		if ((count % 12 == 0) && (count > 0)) {
			//printf("%f\n",sum/12);
			printf("(%f) ",sum/12);
			sum = 0;
		}

        sum += f;
        count++;
    }

	if ((count % 12 == 0) && (count > 0)) {
		//printf("%f\n",sum/12);
		printf("(%f)\n",sum/12);
		sum = 0;
	}
}
