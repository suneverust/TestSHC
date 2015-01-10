/***************************************************************************
  * Modified from "test_s2_rotate_fftw_wrap.c"
  * in package soft
  **************************************************************************
***************************************************************************
  **************************************************************************
  
  SOFT: SO(3) Fourier Transforms
  Version 2.0

  Copyright (c) 2003, 2004, 2007 Peter Kostelec, Dan Rockmore
  
  This file is part of SOFT.

  SOFT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  SOFT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/

/*

 a somewhat memory-friendly test routine to rotate a spherical function
 by massaging its S^2 Fourier coefficients with Wigner-D functions

 bwIn = bandwidth of input signal
 bwOut = bandwidth of output signal (can up- or down-sample)
 degOut = max degree of spherical harmonic you want to use ( < bwOut )
 alpha, beta, gamma -> the three Euler angles

             0 <= alpha, gamma < 2*pi
             0 <= beta <= pi

 inputSamples -> filename of input samples
 outputSamples -> filename of output (rotated) samples

 isReal = 1: samples are strictly real (so no imaginary parts)
 isReal = 0: samples are complex (so in interleaved format)


 Here are order of rotation events:
  1) rotate by gamma about the z-axis
  2) rotate by beta about the y-axis
  3) rotate by alpha about the z-axis.

 example: test_s2_rotate_fftw_wrap bw alpha beta gamma inputSamples outputSamples isReal

 example: test_s2_rotate_fftw_wrap 8 0.37 2.32 4.37 randomS2sig_bw8.dat yyy.dat 0


 NOTE: Sometimes there is a segmentation fault *after* all the rotating and
 writing out of the output file is complete. I haven't tracked this down yet,
 but I believe it has to do with freeing up the memory associated with doing
 the S^2 transforms ... my array of double pointers are not pointing in the
 right places when I try to free memory. However, the rotation itself is
 correct.

 ************************************************************************
 *
 * Modified by Bo Sun
 * From TAMS, University of Hamburg
 * "tams_s2_rotate_fftw" is modified from "test_s2_rotate_fftw_wrap.c"
 * in package soft.
 * Author: Bo Sun
 * Afflication: TAMS, University of Hamburg
 * E-Mail: bosun@informatik.uni-hamburg.de
 *         user_mail@QQ.com
 *
 */

#include <eigen3/Eigen/Dense>
#include <vector>

extern "C" {
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "wrap_fftw.h"
#include "csecond.h"

/**************************************************************/
/**************************************************************/


void tams_s2_rotate_fftw (const Eigen::VectorXf TAMS_S2_real,
                          int TAMS_bandwidth_,
                          double alpha, double beta, double gamma,
                          Eigen::VectorXf &Output_S2)
{
    int bw, i, isReal ;
    double *sigIn, *sigOut ;

    bw = TAMS_bandwidth_;
    isReal = 1;

    if ( isReal )
    {
        sigIn = (double *) malloc(sizeof(double)*(4*bw*bw));
        sigOut = (double *) malloc(sizeof(double)*(4*bw*bw));
    }

    if ( (sigIn == NULL ) || (sigOut == NULL ) )
    {
        perror("Error in allocating memory");
        exit( 1 ) ;
    }

    fprintf(stdout,"reading in signal ...\n");

    /* read in signal */
    if ( isReal )
    {
        for ( i = 0 ; i < (4*bw*bw) ; i ++ )
        {
            sigIn[i] = TAMS_S2_real(i);
        }
    }

    fprintf(stdout,"about to rotate ...\n");

    s2RotateFFTW( bw,
                  sigIn,
                  sigOut,
                  alpha, beta, gamma,
                  isReal );

    fprintf(stdout,"finished rotating ...\n");

    /* write out rotated signal */
    if ( isReal )
    {
        for ( i = 0 ; i < (4*bw*bw) ; i ++ )
        {
            Output_S2(i)= sigOut[i];
        }
    }

    fprintf(stdout,"finished writing ...\n");

    free(sigOut);
    free(sigIn);


}


} // end of extension "C"
