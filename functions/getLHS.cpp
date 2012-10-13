#include <mex.h>
#include <math.h>
#include <matrix.h>

// Matlab mex wrapper
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

  //declare variables
  mxArray *m_in, *G_in;
  const mwSize *dims, *mdim;
  double *m, *G, *LHS;
  int i, j, k, nx, nz, ny;
  long int mcell;

  //unpack MEX wrapper input
  m_in = mxDuplicateArray(prhs[0]);
  G_in = mxDuplicateArray(prhs[1]);
  nx = (int) mxGetScalar(prhs[2]);
  ny = (int) mxGetScalar(prhs[3]);
  nz = (int) mxGetScalar(prhs[4]);

  //figure out dimensions of forward Operator
  dims = mxGetDimensions(G_in);
  mdim = mxGetDimensions(m_in);
  // Check input dimensions
  //dims[0] is should be data dimension. dims[1] is model dimension
  if( (int) dims[0] > (int) dims[1] ) {
    mexPrintf("G matrix should be nxm where n < m\n");
    return;
  }
  if ( (int) dims[1] != (int) mdim[0] ) {
    mexPrintf("Model dimensions should match dimension n in G if G is nxm \n");
    return;
  }

  // Create MEX wrapper output array
  plhs[0] = mxCreateDoubleMatrix(dims[1], 1, mxREAL);

  //associate pointers with MEX Input
  m = mxGetPr(m_in); // model
  G = mxGetPr(G_in); // Forward operator
  //associate pointers with MEX output
  LHS = mxGetPr(plhs[0]);

  ///// Calculate Derivatives ///////////////////
  mcell = nx * ny * nz;
  double dmx[mcell], dmy[mcell], dmz[mcell];
  ///// DX ////////
  for(i = 0; i < mcell; i++) {
    if (i < nz) {
      dmx[i] = m[i] - m[i + nz];
    }
    else if ( (i >= nz) && (i < (mcell - nx * nz)) ) {
      dmx[i] = -m[i - nz] + 2 * m[i] - m[i + nz];
    }
    else if (  (i >= (mcell - nx * nz))  && (i < (mcell - (nx - 1) * nz)) ) {
      dmx[i] = -m[i - nz] + m[i];
    }
    else {
      dmx[i] = 0;
    }
  }

  ///// DY ////////
  for(i = 0; i < mcell; i++) {
    if (i < nz * nx) {
      dmy[i] = m[i] - m[i + nz * nx];
    }
    else if ( (i >= nz * nx) && (i < (mcell - nx * nz)) ) {
      dmy[i] = -m[i - nz * nx] + 2 * m[i] - m[i + nz * nx];
    }
    else {
        dmy[i]= -m[i - nz * nx] + m[i];
    }
  }

  ///// DZ ////////
  for(i = 0; i < mcell; i++) {
    if ( (i % nz) == 0) {
       dmz[i] = m[i] - m[i+1];
    }
    else if ( (i % nz) == (nz - 1) ) {
      dmz[i] = -m[i - 1] + 2 * m[i];
    }
    else {
      dmz[i] = -m[i - 1] + 2 * m[i] - m[i + 1];
    }
  }


  ////// Calculate GtGm /////////////
  double GtG_row[dims[1]];

  for(i = 0; i < dims[1]; i++) {
    LHS[i] = 0;
    for(j = 1; j < dims[1]; j++) {
      GtG_row[j] = 0;
      for(k = 0; k < dims[0]; k++) {
            GtG_row[j] = GtG_row[j] + G[k + dims[0] * i] * G[k + dims[0] * j];
      }
      LHS[i] = LHS[i] + GtG_row[j] * m[j];
    }
    LHS[i] = LHS[i] + dmx[i] + dmy[i] + dmz[i];
  }
}
