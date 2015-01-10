/*********************************************
 * Author: Bo Sun                            *
 * Afflication: TAMS, University of Hamburg  *
 * E-Mail: bosun@informatik.uni-hamburg.de   *
 *         user_mail@QQ.com                  *
 * Date: Dec 02, 2014                        *
 * Licensing: GNU GPL license.               *
 *********************************************/

#ifndef TESTSHC_IMPL_H_
#define TESTSHC_IMPL_H_

#include "TestSHC.h"

void tams_cart2sph(const float x, const float y, const float z,
                   float& azimuth, float& polar)
{
    polar = atan2(hypot(x,y),z);
    azimuth = atan2(y,x);
    if (azimuth < 0)
        azimuth = azimuth + 2*M_PI;
}

void tams_vector_normalization (std::vector<float> &tams_vector)
{
    float max_element = (*std::max_element(tams_vector.begin(),tams_vector.end()));
    float min_element = (*std::min_element(tams_vector.begin(),tams_vector.end()));

    for (std::vector<float>::iterator itr = tams_vector.begin();
         itr != tams_vector.end(); itr ++)
    {
        // save memory but dangerous!
        (*itr)=((*itr)-min_element)/(max_element-min_element);
    }
}

void tams_vector2entropy (const std::vector<float> tams_vector,
                          const size_t hist_bin,
                          float &entropy)
{
    std::vector<float> temp_hist(hist_bin+1);
    for (std::vector<float>::const_iterator itr = tams_vector.begin();
         itr !=tams_vector.end(); itr++)
    {
        temp_hist[floor((*itr)*hist_bin)]++;
    }
    temp_hist[hist_bin-1]++;
    temp_hist.pop_back();
    if(temp_hist.size()!=hist_bin)
        pcl::console::print_warn(
                    "Warning: something maybe wrong in computing Histogram!\n");

    // Parzen Window: [0.05, 0.25, 0.40, 0.25, 0.05]
    std::vector<float> temp_hist_pad;
    temp_hist_pad.push_back(0.0);
    temp_hist_pad.push_back(0.0);
    for(std::vector<float>::iterator itr = temp_hist.begin();
        itr != temp_hist.end(); itr++)
    {
        temp_hist_pad.push_back(*itr);
    }
    temp_hist_pad.push_back(0.0);
    temp_hist_pad.push_back(0.0);
    std::vector<float>().swap(temp_hist);

    std::vector<float> tams_hist;
    for(std::vector<float>::iterator itr=temp_hist_pad.begin()+2;
        itr !=temp_hist_pad.end()-2; itr++)
    {
        tams_hist.push_back( (*(itr-2))*0.05
                             +(*(itr-1))*0.25
                             +(*(itr  ))*0.40
                             +(*(itr+1))*0.25
                             +(*(itr+2))*0.05);
    }
    if(tams_hist.size()!=hist_bin)
    {
        pcl::console::print_error("Error: Histogram Parzen Window Failed\n");
        return;
    }
    std::vector<float>().swap(temp_hist_pad);
    
    entropy = 0.0;
    for (std::vector<float>::iterator itr = tams_hist.begin();
         itr !=tams_hist.end(); itr++)
    {
        if ((*itr)>0)
            entropy += -(*itr)*log((*itr));
    }
    std::vector<float>().swap(tams_hist);
}

void computeSEI (const PointCloudT cloud,
                 size_t sei_dim,
                 size_t hist_bin,
                 Eigen::MatrixXf &entropy)
{
    float sei_azimuth_spa = 2*M_PI/(2*sei_dim);
    float sei_polar_spa   = M_PI/(2*sei_dim);

    // Point Division
    Eigen::Array<std::vector<float>, Eigen::Dynamic, Eigen::Dynamic>
            sei_points_divi(2*sei_dim, 2*sei_dim);

    float temp_az , temp_polar;
    size_t temp_sei_azth, temp_sei_polarth;
    double dist;
    for (PointCloudT::const_iterator itr=cloud.begin();
         itr!=cloud.end(); itr++)
    {
        dist = sqrt((*itr).x*(*itr).x+
                    (*itr).y*(*itr).y+
                    (*itr).z*(*itr).z);
        tams_cart2sph((*itr).x, (*itr).y, (*itr).z,
                      temp_az, temp_polar);
        if (temp_az < sei_azimuth_spa/2 ||
                temp_az >= 2*M_PI-sei_azimuth_spa/2)
            temp_sei_azth = 0;
        else
            temp_sei_azth = floor((temp_az-sei_azimuth_spa/2)/sei_azimuth_spa)+1;
        temp_sei_polarth = floor(temp_polar/sei_polar_spa);
        if (temp_sei_azth > 2*sei_dim)
        {
            std::cout << "temp_az:" <<temp_az << std::endl;
        }
        if (temp_sei_polarth > 2*sei_dim)
        {
            std::cout << "temp_polar:"<<temp_polar << std::endl;
        }
        sei_points_divi(temp_sei_azth, temp_sei_polarth).push_back(dist);
    }

    // compute entropy
    for(temp_sei_azth = 0; temp_sei_azth < 2*sei_dim; temp_sei_azth++)
    {
        for(temp_sei_polarth = 0; temp_sei_polarth < 2*sei_dim; temp_sei_polarth++)
        {
            if (sei_points_divi(temp_sei_azth, temp_sei_polarth).size()<5)
            {
                entropy(temp_sei_azth, temp_sei_polarth) = 0;
                continue;
            }
            if (    (*std::max_element(sei_points_divi(temp_sei_azth, temp_sei_polarth).begin(),
                                       sei_points_divi(temp_sei_azth, temp_sei_polarth).end()))
                    ==
                    (*std::min_element(sei_points_divi(temp_sei_azth, temp_sei_polarth).begin(),
                                       sei_points_divi(temp_sei_azth, temp_sei_polarth).end())))
            {
                entropy(temp_sei_azth, temp_sei_polarth) = 0;
                continue;
            }

            tams_vector_normalization(sei_points_divi(temp_sei_azth, temp_sei_polarth));

            tams_vector2entropy(sei_points_divi(temp_sei_azth, temp_sei_polarth),
                                hist_bin,
                                entropy(temp_sei_azth, temp_sei_polarth));
        }
    }
}

#endif /*TESTSHC_IMPL_H_*/
