/*********************************************
 * Author: Bo Sun                            *
 * Afflication: TAMS, University of Hamburg  *
 * E-Mail: bosun@informatik.uni-hamburg.de   *
 *         user_mail@QQ.com                  *
 * Date: Dec 02, 2014                        *
 * Licensing: GNU GPL license.               *
 *********************************************/

#ifndef TESTSHC_H_
#define TESTSHC_H_

#include <cmath>
#include <vector>
#include <algorithm>
#include <eigen3/Eigen/Dense>

// Types
#ifndef TYPE_DEFINITION_
#define TYPE_DEFINITION_
typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
#endif /*TYPE_DEFINITION_*/

/** \brief cart2sph(X,Y,Z,azimuth,polar) transforms Cartesian coordinates stored in
  * corresponding elements of X, Y, and Z into spherical coordinates.
  * azimuth and polar are angular displacements in radians.
  * azimuth(longitudinal) is the counterclockwise angle in the x-y plane measured from the positive x-axis.
  * polar(colatitudianl) is the polar angle measured from the positive z axis.
  * 0 < azimuth < 2*M_PI; 0 < polar < M_PI
  */
void tams_cart2sph (float x, float y, float z,
                    float& azimuth, float& polar);

/** \brief tams_vector_normalization normalize the input vector
  * Parameters:
  * \param[in]     tams_vector   the input vector
  * \param[out]    tams_vector   the normalized vector (values are in range [0,1])
  */
void tams_vector_normalization (std::vector<float> &tams_vector);

/** \brief tams_vector2entropy compute the entropy of a vector
  * Parameters:
  * [in]   tams_vector     the input vector
  * [in]   hist_bin        the size of histogram in entropy computation
  * [out]  entropy         the resultant entropy
  */
void tams_vector2entropy( const std::vector<float> tams_vector,
                          const size_t hist_bin,
                          float& entropy);

/** \brief computeSEI calculate the SEI of the input point cloud
  * Parameters:
  * [in]     cloud         the input vector
  * [in]     sei_dim       the dimension of SEI(sei_dimXsei_dim)
  * [in]     hist_bin      the size of histogram in entropy computation
  * [out]    entropy       the resulstant SEI stored in a (sei_dim X sei_dim) matrix
  */
void computeSEI (const PointCloudT cloud,
                 size_t sei_dim,
                 size_t hist_bin,
                 Eigen::MatrixXf &entropy);

#endif /*TESTSHC_H_*/
