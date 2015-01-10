# TestSHC

TestSHC is an example package to show that: 
the magnitude of Spherical Harmonic Coefficients (SHC) of a spherical function
does not change when that spherical funciton is rotated aroud Z axis. 

1. To produce the spherical function, this package receive a point cloud file 
   and compute its Spherical Entropy Image (SEI) which is a spherical function
2. We use two methods to rotate the spherical function                                                                   1) rotate the spherical function using the 's2RotateFFTW' in library s2kit
      (we modified the original function and stored the modified file in subdirectory 's2rotate')                        2) rotate the spherical function by hand
3. We compute the SHC based on the library s2kit
   (we modified the original function and stored the modified file in subdirectory 's2kit')

-------
The package contains two subdirectories: 
  --'s2kit' which contains the files used for computing SHC of spherical functions, 
    called by 'tams_s2_semi_memo_for'
  --'s2rotate' which contains the files used for rotate spherical functions, 
    called by 'tams_s2_rotate_fftw'
  

