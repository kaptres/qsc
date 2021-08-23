// Author Murat Genctav

#include "mex.h"
#include "matrix.h"

void buildRHS(double *b,int numShapeNodes,int numScales,double *edt,int *shapeNodeIndices,int numNodes,double scaleRes) {
	double s;
	int t,k;
	int cId;
	
	for(t=0;t<numScales;t++) {
		// calculate current scale
		s = t*scaleRes;
		for(k=0;k<numNodes;k++) {
			cId = shapeNodeIndices[k];
			if(cId>0) { // current node is inside the shape
				if(edt[k]>s) { // belongs to positive region
					b[t*numShapeNodes+cId-1] = 1.0;
				} else { // belongs to negative region
					b[t*numShapeNodes+cId-1] = -1.0;
				}
			}
		}
	}
}

void buildSystemMatrix(double *aval,mwIndex *arow,mwIndex *acol,int numShapeNodes,const mwSize *dims,int *shapeNodeIndices,double alpha) {
	int i,j,k;
	int I,J,K,jI,kIJ;
	int cId;
	int numNonZeros;
	double centerValue;
	
	centerValue = -6-alpha;
	
	numNonZeros = 0;
	I = dims[0];
	J = dims[1];
	K = dims[2];
	for(k=1;k<K-1;k++) {
		kIJ = k*I*J;
		for(j=1;j<J-1;j++) {
			jI = j*I;
			for(i=1;i<I-1;i++) {
				cId = shapeNodeIndices[i+jI+kIJ];
				if (cId > 0) { // a shape node
					acol[cId-1] = numNonZeros;
					// add neighboring elements before
					if (shapeNodeIndices[i + jI + (k - 1)*I*J] > 0) {
						aval[numNonZeros] = 1.0;
						arow[numNonZeros] = shapeNodeIndices[i + jI + (k - 1)*I*J] - 1;
						numNonZeros++;
					}
					if (shapeNodeIndices[i + (j - 1)*I + kIJ] > 0) {
						aval[numNonZeros] = 1.0;
						arow[numNonZeros] = shapeNodeIndices[i + (j - 1)*I + kIJ] - 1;
						numNonZeros++;
					}
					if (shapeNodeIndices[(i - 1) + jI + kIJ] > 0) {
						aval[numNonZeros] = 1.0;
						arow[numNonZeros] = shapeNodeIndices[(i - 1) + jI + kIJ] - 1;
						numNonZeros++;
					}
					// add diagonal element
					
					aval[numNonZeros] = centerValue;
					arow[numNonZeros] = cId - 1;
					numNonZeros++;
					// add neighboring elements after
					if (shapeNodeIndices[(i + 1) + jI + kIJ] > 0) {
						aval[numNonZeros] = 1.0;
						arow[numNonZeros] = shapeNodeIndices[(i + 1) + jI + kIJ] - 1;
						numNonZeros++;
					}
					if (shapeNodeIndices[i + (j + 1)*I + kIJ] > 0) {
						aval[numNonZeros] = 1.0;
						arow[numNonZeros] = shapeNodeIndices[i + (j + 1)*I + kIJ] - 1;
						numNonZeros++;
					}
					if (shapeNodeIndices[i + jI + (k + 1)*I*J] > 0) {
						aval[numNonZeros] = 1.0;
						arow[numNonZeros] = shapeNodeIndices[i + jI + (k + 1)*I*J] - 1;
						numNonZeros++;
					}
				}
			}
		}
	}
	acol[numShapeNodes] = numNonZeros;
}

void assignShapeNodeIndices(int *shapeNodeIndices,int *numShapeNodes,double *edt,int numNodes) {
	int i;
	
	// counting and indexing shape nodes
	*numShapeNodes = 0;
	for(i=0;i<numNodes;i++) {
		if(edt[i]>0) {
			(*numShapeNodes)++;
			shapeNodeIndices[i] = *numShapeNodes; // Attention! Indices are 1 based!
		}
	}
}

/* Usage:
 * [A,B] = buildScreenedPoissonSystemMex(EDT,alpha,scaleRes,numScales);
 */

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *edt;
	double alpha;
	double scaleRes;
	int numScales;
	const mwSize *dims;
	int *shapeNodeIndices;
	int numNodes;
	int numShapeNodes;
	double *aval;
	mwIndex *arow;
	mwIndex *acol;
	double *b;
	
	/*
	printf("Setting up the screened Poisson solver\n");
	printf("Issuing input arguments...\n");
	/* */
	
	/* Handle input arguments */
	/* arg0: distance transform */
	edt = mxGetPr(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	numNodes = dims[0]*dims[1]*dims[2];
	/* arg1: screening parameter */
	alpha = mxGetScalar(prhs[1]);
	/* arg2: scale resolution */
	scaleRes = mxGetScalar(prhs[2]);
	/* arg3: # of scales = # of right hand sides */
	numScales = (int) mxGetScalar(prhs[3]);
	
	/*
	printf("\tGrid dimensions: %dx%dx%d\n",dims[0],dims[1],dims[2]);
	printf("\talpha = %.5e\n",alpha);
	printf("\tScale resolution (s step size) = %lf\n",scaleRes);
	printf("\tNumber of scales to solve for = %d\n",numScales);
	printf("Done\n");
	/* */
	
	/*
	printf("Assigning index numbers to grid points inside the shape...");
	/* */

	/* Grid numbering */
	// allocate memory for index array
	shapeNodeIndices = (int *) calloc(numNodes,sizeof(int));
	assignShapeNodeIndices(shapeNodeIndices,&numShapeNodes,edt,numNodes);

	/*
	printf("Done\n");
	/* */
	
	/*
	printf("Allocating memory for output (System matrix and RHS)...");
	/* */

	/* Allocate memory for output arguments */
	/* larg0: system matrix */
	plhs[0] = mxCreateSparse(numShapeNodes,numShapeNodes,numShapeNodes*7,mxREAL);
	aval = mxGetPr(plhs[0]);
	arow = mxGetIr(plhs[0]);
	acol = mxGetJc(plhs[0]);
	/* larg1: multiple right hand sides */
	plhs[1] = mxCreateDoubleMatrix(numShapeNodes,numScales,mxREAL);
	b = mxGetPr(plhs[1]);

	/*
	printf("Done\n");
	/* */
	
	/*
	printf("Building system matrix...");
	/* */

	/* Build system matrix */
	buildSystemMatrix(aval,arow,acol,numShapeNodes,dims,shapeNodeIndices,alpha);

	/*
	printf("Done\n");
	/* */

	/*
	printf("Building multiple RHSs...");
	/* */
	
	/* Build multiple right hand sides */
	buildRHS(b,numShapeNodes,numScales,edt,shapeNodeIndices,numNodes,scaleRes);

	/*
	printf("Done\n");
	/* */
	
	/* Clean up allocated memory for internal processing */
	free(shapeNodeIndices);
}