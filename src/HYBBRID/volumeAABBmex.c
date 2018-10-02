/* VOLUMEAABB Computation of the volume of the minimal AABB enclosing a dataset.
 *
 *    VOLUMEAABB(data, R) computes the volume of the minimal axis-aligned 
 *      bounding box enclosing the set of points described in data (one point
 *      per row), after a rotation defined by the matrix R.
 *
 * output #1 : Volume of the AABB enclosing data*R'.
 *
 * REMARKS:
 *    This is a MEX-function corresponding to the Matlab function volumeAABB.m.
 *
 * SEE ALSO:
 *    HYBBRID
 *
 * AUTHORS:
 *    Chia-Tche Chang <cchang.uclouvain@gmail.com>
 *    Bastien Gorissen <bastien@panopticgame.com>
 *    Samuel Melchior <samuel.melchior@epfl.ch>
 *
 * REFERENCES:
 *    - C.-T.Chang, B.Gorissen, S.Melchior, 
 *      "Fast oriented bounding box optimization on the rotation group SO(3, R)",
 *      submitted to ACM Transactions on Graphics, 2011.
 *
 * This file is part of HYBBRID.
 * Copyright © 2010 Chia-Tche Chang, Bastien Gorissen, Samuel Melchior
 *
 * HYBBRID is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HYBBRID is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HYBBRID. If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include "mex.h"

#define	DATA_IN	prhs[0]
#define	R_IN	prhs[1]
#define	VOL_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

static void volumeAABBmex(double* vol, double* data, double* r, unsigned int nbData) {
    unsigned int i, j, s;
    double r1, r2, r3;
    double min, max, val;
    *vol = 1.0;
    
    for(i = 0 ; i < 3 ; i++) {
        r1 = r[i];
        r2 = r[i+3];
        r3 = r[i+6];
        min = data[0]*r1 + data[nbData]*r2 + data[2*nbData]*r3;
        max = min;
        for(j = 1 ; j < nbData ; j++) {
            val = data[j]*r1 + data[nbData+j]*r2 + data[2*nbData+j]*r3;
            min = MIN(min, val);
            max = MAX(max, val);
        }
        *vol *= max - min;
    }
    
    return;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 2)
        return;

    /* Create a matrix for the return argument */
    VOL_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    /* Assign pointers to the various parameters */
    double* vol = mxGetPr(VOL_OUT);
    double* data = mxGetPr(DATA_IN);
    double* r = mxGetPr(R_IN);
    unsigned int nbData = mxGetM(DATA_IN);
    
    /* Do the actual computations in a subroutine */
    volumeAABBmex(vol, data, r, nbData);
    return;
}
