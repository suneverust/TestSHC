# Afflication: TAMS, University of Hamburg
# E-Mail: bosun@informatik.uni-hamburg.de

# We gather the ESSENTIAL files for compute the spherical harmonics 
# in s2kit10, and write this CMakeLists.txt in order to let CMake to 
# configure and compile the project. 
 
# We use "FindFFTW3.cmake" to let CMake find the FFTW3_INCLUDE_DIR 
# and FFTW3_LIBRARIES automatically, please make sure 
# the "FindFFTW3.cmake" at the right place.


cmake_minimum_required (VERSION 2.8)

project(s2rotate)

find_package (FFTW3 REQUIRED)

include_directories (${FFTW3_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/s2rotate/include)

link_directories (${FFTW3_LIBRARIES})

add_library (s2rotate 	csecond.c 
			cospmls.c
			FFTcode.c
			fft_grids.c
			FST_semi_memo.c
			indextables.c
			legendreTransforms.c
                        makeweights.c
                        makeWigner.c
			newFCT.c
			oddweights.c
			OURmods.c
			OURperms.c
			permroots.c
                        rotate_so3_fftw.c
                        rotate_so3_utils.c
                        s2_cospmls.c
                        s2_legendreTransforms.c
                        s2_primitive.c
                        s2_semi_fly.c
                        s2_semi_memo.c
			so3_correlate_fftw.c
			soft_fftw.c
                        utils_so3.c
                        utils_vec_cx.c
			weights.c
                        wignerTransforms_fftw.c
                        wrap_s2_rotate_fftw.c
			wrap_soft_fftw.c
			wrap_soft_fftw_cor2.c)

target_link_libraries (s2rotate ${FFTW3_LIBRARIES} -lm)
