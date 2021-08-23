// author: mferhata
#include "mex.h"
#include <math.h>

// compute v function of a 3D shape
// 
// Usage: v = spe_Linfty_iterative_mex_3D 
//                  (S, rho, delta_time, tolerance, iter_limit, v0)
// 
// INPUTS
// S                - 3D shape voxel (M x N x K)
//                    true at shape pixels, false at remaining pixels
// rho              - v function smoothness parameter (scalar)
// delta_time       - time step
//                    in range (0,1/6] for stability
// tolerance        - convergence criterion (scalar)
//                    Stop iterations if maximum of absolute v function
//                    difference between consecutive iterations is smaller
//                    than tolerance.
// iter_limit       - maximum iteration number (scalar)
// v0               - (optional) initial v function
//                    If not supplied, initial v function is equal to 0.
// OUTPUT
// v                - 3D v function of shape (M x N x K)
//                    v function is equal to 0 on pixels not included in 
//                    shape.
// 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    #define SHAPE_VOX_IN        prhs[0]
    #define RHO_IN              prhs[1]
    #define DELTA_TIME_IN       prhs[2]
    #define ITER_STOP_TOL_IN    prhs[3]
    #define N_MAX_ITERS_IN      prhs[4]
    #define V0_IN               prhs[5]
    #define V_OUT               plhs[0]
    
    int dim1;
    int dim2;
    int dim3;
    bool *S;
    int boundary_value;
    float rho;
    float delta_time;
    float tolerance;
    int iter_limit;
    float *v0;
    float *v;
    
    int voxel_cnt;
    int i, j, k;
    int di, dj, dk;
    int *dim1_inds;
    int *dim2_inds;
    int *dim3_inds;
    int count;
    
    int     n;
    float val;
    float term;
    float temp;
    float max, min;
    float *v_change;
    float change;
    float max_change;
    
    const mwSize *dims  = mxGetDimensions (SHAPE_VOX_IN);
    dim1                = dims[0];
    dim2                = dims[1];
    dim3                = dims[2];
    S                   = (bool *) mxGetData (SHAPE_VOX_IN);
    rho                 = mxGetScalar (RHO_IN);
    delta_time          = mxGetScalar (DELTA_TIME_IN);
    tolerance           = mxGetScalar (ITER_STOP_TOL_IN);
    iter_limit          = mxGetScalar (N_MAX_ITERS_IN);
    
    printf ("size:       %dx%dx%d\n",   dim1, dim2, dim3);
    //printf ("rho:        %f\n",         rho);
    //printf ("delta_time: %f\n",         delta_time);
    //printf ("tolerance:  %e\n",         tolerance);
    //printf ("iter_limit: %d\n",         iter_limit);
    
    V_OUT   = mxCreateNumericArray (3, dims, mxSINGLE_CLASS, mxREAL);
    v       = (float *) mxGetData (V_OUT);
    
    // v initialization
    if (nrhs == 6)
    {
        v0 = (float *) mxGetData (V0_IN);
        for (i = 0; i < dim1*dim2*dim3; i++)
        {
            if (S[i])
                v[i] = v0[i];
            else
                v[i] = 0;
        }
    }
    else
    {
        for (i = 0; i < dim1*dim2*dim3; i++)
            v[i] = 0;
    }

    // count number of shape pixels
    voxel_cnt = 0;
    for (i = 0; i < dim1*dim2*dim3; i++)
    {
        if (S[i])
            voxel_cnt++;
    }
    printf ("voxel_cnt:  %d\n", voxel_cnt);
    
    // get indices of shape pixels
    dim1_inds   = mxMalloc (sizeof(int) * voxel_cnt);
    dim2_inds   = mxMalloc (sizeof(int) * voxel_cnt);
    dim3_inds   = mxMalloc (sizeof(int) * voxel_cnt);
    count = 0;
    for (k=0; k<dim3; k++)
    {
        for (j=0; j<dim2; j++)
        {
            for (i=0; i<dim1; i++)
            {
                if (S[k*dim1*dim2 + j*dim1 + i])
                {
                    dim1_inds[count] = i;
                    dim2_inds[count] = j;
                    dim3_inds[count] = k;
                    count++;
                }
            }
        }
    }
    
    v_change    = mxMalloc (sizeof(float) * voxel_cnt);
    term        = 1 / (rho * rho);
    for (n=0; n<iter_limit; n++)
    {
        // compute the change
        for (count=0; count<voxel_cnt; count++)
        {
            i   = dim1_inds[count];
            j   = dim2_inds[count];
            k   = dim3_inds[count];
            
            max= INT_MIN;
            min= INT_MAX;
            // The strel is hard-coded as a 3x3 cube
            for (dk=-1; dk<2; dk++){
                for (dj=-1; dj<2; dj++){
                    for (di=-1; di<2; di++){
                        temp = v[(k+dk)*dim1*dim2 + (j+dj)*dim1 + i+di];
                        if (temp > max)
                            max = temp;
                        if (temp < min)
                            min = temp;
                    }
                }
            }

            val = v[k*dim1*dim2 + j*dim1 + i];
            //              |    Inf_Laplacian    | Screening  |
            //              | v v v v v v v v v v | v v v v v  |
            v_change[count] = max + min - 2 * val - term * val + 1;
        }
        
        // compute maximum of absolute v difference between consecutive iterations
        // update v
        max_change = INT_MIN;
        for (count=0; count<voxel_cnt; count++)
        {
            i       = dim1_inds[count];
            j       = dim2_inds[count];
            k       = dim3_inds[count];
            
            change  = fabs (delta_time * v_change[count]);
            if (change > max_change)
                max_change = change;
            
            v[k*dim1*dim2 + j*dim1 + i] += delta_time * v_change[count];
        }
        if (n % 50 == 49)
            printf ("%d max change is %.5f\n", n+1, max_change);
        
        if (max_change < tolerance)
            break;
    }
    
    printf ("number of iterations: %d max_change: %e\n", n, max_change);
    
    mxFree (dim1_inds);
    mxFree (dim2_inds);
    mxFree (dim3_inds);
    mxFree (v_change);
    
    return;
}
