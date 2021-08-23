// Author Murat Genctav muratgenctav@gmail.com
#include "mex.h"
#include "matrix.h"

void buildSystemMatrixRHS(double *aval,mwIndex *arow,mwIndex *acol,int numShapeNodes,const mwSize *dims,int *shapeNodeIndices,double alpha,double *b) {
	int i,j,k;
	int I,J,K,jI,kIJ;
	int cId;
	int numNonZeros;
	double centerValue;
	
	centerValue = -4-alpha;
	
	numNonZeros = 0;
	I = dims[0];
	J = dims[1];
	for(j=1;j<J-1;j++) {
		jI = j*I;
		for(i=1;i<I-1;i++) {
			cId = shapeNodeIndices[i+jI];
			if (cId > 0) { // a shape node
                b[cId-1] = 0;
				acol[cId-1] = numNonZeros;
				// add neighboring elements before
				if (shapeNodeIndices[i + (j - 1)*I] > 0) {
					aval[numNonZeros] = 1.0;
					arow[numNonZeros] = shapeNodeIndices[i + (j - 1)*I] - 1;
					numNonZeros++;
				}
                else {
                    b[cId-1] = b[cId-1]-1;
                }
				if (shapeNodeIndices[(i - 1) + jI] > 0) {
					aval[numNonZeros] = 1.0;
					arow[numNonZeros] = shapeNodeIndices[(i - 1) + jI] - 1;
					numNonZeros++;
				}
                else {
                    b[cId-1] = b[cId-1]-1;
                }
				// add diagonal element
				aval[numNonZeros] = centerValue;
				arow[numNonZeros] = cId - 1;
				numNonZeros++;
				// add neighboring elements after
				if (shapeNodeIndices[(i + 1) + jI] > 0) {
					aval[numNonZeros] = 1.0;
					arow[numNonZeros] = shapeNodeIndices[(i + 1) + jI] - 1;
					numNonZeros++;
				}
                else {
                    b[cId-1] = b[cId-1]-1;
                }
				if (shapeNodeIndices[i + (j + 1)*I] > 0) {
					aval[numNonZeros] = 1.0;
					arow[numNonZeros] = shapeNodeIndices[i + (j + 1)*I] - 1;
					numNonZeros++;
				}
                else {
                    b[cId-1] = b[cId-1]-1;
                }
			}
		}
	}
	acol[numShapeNodes] = numNonZeros;
}

void assignShapeNodeIndices(int *shapeNodeIndices,int *numShapeNodes,double *mask,int numNodes) {
	int i;	
	// counting and indexing shape nodes
	*numShapeNodes = 0;
	for(i=0;i<numNodes;i++) {
		if(mask[i]>0) {
			(*numShapeNodes)++;
			shapeNodeIndices[i] = *numShapeNodes; // Attention! Indices are 1 based!
		}
	}
}

/* Usage:
 * [A,b] = buildScreenedPoissonSystemMex(mask,alpha);
 */

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *mask;
	double alpha;
	const mwSize *dims;
	int *shapeNodeIndices;
	int numNodes;
	int numShapeNodes;
	double *aval;
	mwIndex *arow;
	mwIndex *acol;
	double *b;
	
	/* arg0: distance transform */
	mask = mxGetPr(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	numNodes = dims[0]*dims[1];
	/* arg1: screening parameter */
	alpha = mxGetScalar(prhs[1]);
	
	/* Grid numbering */
	// allocate memory for index array
	shapeNodeIndices = (int *) calloc(numNodes,sizeof(int));
	assignShapeNodeIndices(shapeNodeIndices,&numShapeNodes,mask,numNodes);

	/* Allocate memory for output arguments */
	/* larg0: system matrix */
	plhs[0] = mxCreateSparse(numShapeNodes,numShapeNodes,numShapeNodes*5,mxREAL);
	aval = mxGetPr(plhs[0]);
	arow = mxGetIr(plhs[0]);
	acol = mxGetJc(plhs[0]);
	/* larg1: right hand side */
	plhs[1] = mxCreateDoubleMatrix(numShapeNodes,1,mxREAL);
	b = mxGetPr(plhs[1]);

	/* Build system matrix */
	buildSystemMatrixRHS(aval,arow,acol,numShapeNodes,dims,shapeNodeIndices,alpha,b);
	
	/* Clean up allocated memory for internal processing */
	free(shapeNodeIndices);
}