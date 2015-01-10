/*********************************************
 * Author: Bo Sun                            *
 * Afflication: TAMS, University of Hamburg  *
 * E-Mail: bosun@informatik.uni-hamburg.de   *
 *         user_mail@QQ.com                  *
 * Date: Dec 02, 2014                        *
 * Licensing: GNU GPL license.               *
 *********************************************/

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/filter.h>
#include <pcl/console/print.h>
#include <pcl/io/pcd_io.h>

#include <fstream>

#include "TestSHC.h"
#include "TestSHC.hpp"

extern "C"{
#include "tams_s2_rotate_fftw.h"
#include "tams_s2_semi_memo_for.h"
}

// Types
#ifndef TYPE_DEFINITION_
#define TYPE_DEFINITION_
typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
#endif /*TYPE_DEFINITION_*/

int 
main(int argc, char **argv)
{
    // Initiate Point Cloud
    PointCloudT::Ptr object(new PointCloudT);

    pcl::console::print_highlight("Loading point clouds...\n");
    if (pcl::io::loadPCDFile<PointT>(argv[1],*object)<0)
    {
        pcl::console::print_error("Error loading object file!\n");
    }

    // Remove the NaN points in object if any
    pcl::console::print_highlight("Removing the NaN points if any...\n");
    std::vector<int> dummy_indices;
    pcl::removeNaNFromPointCloud(*object, *object, dummy_indices);
    std::vector<int>().swap(dummy_indices);

    size_t sei_dim = 16;
    size_t hist_bin = 12;
    size_t tams_bandwidth_ = sei_dim;

    ///////////////////////////////////////////////////////////////////////////
    // compute the SEI of object, SEI is a spherical function
    Eigen::MatrixXf sei_object = Eigen::MatrixXf::Zero(2*sei_dim, 2*sei_dim);
    computeSEI(*object, sei_dim, hist_bin, sei_object);
    (*object).resize(0);

    // compute the SHC of SEI
    Eigen::VectorXf sei_object_real = Eigen::VectorXf::Zero(2*sei_dim*2*sei_dim);
    sei_object.resize(2*sei_dim*2*sei_dim,1);
    sei_object_real << sei_object;
    sei_object.resize(2*sei_dim, 2*sei_dim);
    pcl::console::print_highlight("Computing the SHC of SEI (object)...\n");
    std::vector<double> object_sh_real;
    std::vector<double> object_sh_imag;
    tams_s2_semi_memo_for(sei_object_real,
                          tams_bandwidth_,
                          object_sh_real,
                          object_sh_imag);

    if (object_sh_real.size()!=object_sh_imag.size())
    {
        pcl::console::print_error("Error compute SHC...\n");
    }

    // compute the magnitude of SHC
    pcl::console::print_highlight("Computing the magnitude of SHC (object)...\n");
    std::vector<float> object_sh_mag(object_sh_real.size());

    for(int i=0; i<object_sh_mag.size(); i++)
    {
        object_sh_mag[i] = sqrt(object_sh_real[i]*object_sh_real[i]+
                                object_sh_imag[i]*object_sh_imag[i]);
    }

    // store the magnitude of SHC (code-order)
    std::ofstream file1("object_mag.txt");
    file1 << "Here is the magnitude of SHC (object)..\n";
    for(std::vector<float>::iterator itr = object_sh_mag.begin();
        itr != object_sh_mag.end();itr++)
    {
        file1 << (*itr) << '\n';
    }
    file1.close();
    std::vector<double>().swap(object_sh_real);
    std::vector<double>().swap(object_sh_imag);
    std::vector<float>().swap(object_sh_mag);

    ////////////////////////////////////////////////////////////////////////////////
    // rotate the SEI by soft2.0
    pcl::console::print_highlight("Rotating the SEI (by soft2.0)...\n");
    Eigen::VectorXf sei_object_rot1_real = Eigen::VectorXf::Zero(2*sei_dim*2*sei_dim);
    float alpha = M_PI/4;
    float beta =  0.0;
    float gamma = 0.0;
    tams_s2_rotate_fftw(sei_object_real, tams_bandwidth_,
                        alpha, beta, gamma,
                        sei_object_rot1_real);

    // compute the SHC of Rotated SEI
    pcl::console::print_highlight("Computing the SHC of SEI (object_rot by library soft2.0)...\n");
    std::vector<double> object_rot1_sh_real;
    std::vector<double> object_rot1_sh_imag;
    tams_s2_semi_memo_for(sei_object_rot1_real,
                          tams_bandwidth_,
                          object_rot1_sh_real,
                          object_rot1_sh_imag);
    if(object_rot1_sh_real.size()!=object_rot1_sh_imag.size())
        pcl::console::print_error("Error compute SHC...\n");

    // compute the magnitude of SHC
    pcl::console::print_highlight("Computing the magnitude of SHC (object_rot by library soft2.0)...\n");
    std::vector<float> object_rot1_sh_mag(object_rot1_sh_real.size());
    for (int i=0; i<object_rot1_sh_real.size(); i++)
        object_rot1_sh_mag[i] = sqrt(object_rot1_sh_real.at(i)*object_rot1_sh_real.at(i)+
                                  object_rot1_sh_imag.at(i)*object_rot1_sh_imag.at(i));


    std::ofstream file3("object_rot_soft_mag.txt");
    file3 << "Here is the magnitude of SHC (object_rot by library soft2.0)...\n";
    for(std::vector<float>::iterator itr = object_rot1_sh_mag.begin();
        itr != object_rot1_sh_mag.end();itr++)
    {
        file3 << (*itr) << '\n';
    }
    file3.close();
    std::vector<double>().swap(object_rot1_sh_real);
    std::vector<double>().swap(object_rot1_sh_imag);
    std::vector<float>().swap(object_rot1_sh_mag);

    ////////////////////////////////////////////////////////////////////////////
    // rotate the SEI by hand
    pcl::console::print_highlight("Rotating the SEI (by hand)...\n");
    int width_offset = 20;
    Eigen::MatrixXf sei_object_rot = Eigen::MatrixXf::Zero(2*sei_dim, 2*sei_dim);
    for (int i= 0; i<2*sei_dim; i++)
    {
        for (int j=0; j<2*sei_dim; j++)
        {
            if (i-width_offset<0)
                sei_object_rot(i,j) = sei_object(i-width_offset+2*sei_dim,j);
            else if (i-width_offset>=2*sei_dim)
                sei_object_rot(i,j) = sei_object(i-width_offset-2*sei_dim,j);
            else
                sei_object_rot(i,j) = sei_object(i-width_offset,j);
        }
    }

    // compute the SHC of Rotated SEI
    Eigen::VectorXf sei_object_rot_real = Eigen::VectorXf::Zero(2*sei_dim*2*sei_dim);
    sei_object_rot.resize(2*sei_dim*2*sei_dim,1);
    sei_object_rot_real << sei_object_rot;
    pcl::console::print_highlight("Computing the SHC of SEI (object_rot by hand)...\n");
    std::vector<double> object_rot_sh_real;
    std::vector<double> object_rot_sh_imag;
    tams_s2_semi_memo_for(sei_object_rot_real,
                          tams_bandwidth_,
                          object_rot_sh_real,
                          object_rot_sh_imag);

    if (object_rot_sh_real.size()!=object_rot_sh_imag.size())
    {
        pcl::console::print_error("Error compute SHC...\n");
    }

    // compute the magnitude of SHC
    pcl::console::print_highlight("Computing the magnitude of SHC (object_rot by hand)...\n");
    std::vector<float> object_rot_sh_mag(object_rot_sh_real.size());
    for(int i=0; i<object_rot_sh_real.size(); i++)
    {
        object_rot_sh_mag.at(i) = sqrt(object_rot_sh_real.at(i)*object_rot_sh_real.at(i)+
                                   object_rot_sh_imag.at(i)*object_rot_sh_imag.at(i));
    }

    // store the magnitude of SHC (code-order)
    std::ofstream file2("object_rot_hand_mag.txt");
    file2 << "Here is the magnitude of SHC (object_rot by hand)..\n";
    for(std::vector<float>::iterator itr = object_rot_sh_mag.begin();
        itr != object_rot_sh_mag.end();itr++)
    {
        file2 << (*itr) << '\n';
    }
    file2.close();
    std::vector<double>().swap(object_rot_sh_real);
    std::vector<double>().swap(object_rot_sh_imag);
    std::vector<float>().swap(object_rot_sh_mag);
}
