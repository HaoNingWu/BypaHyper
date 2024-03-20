# BypaHyper
This is a collection of MATLAB codes to reproduce all the figures in our paper "Bypassing the quadrature exactness assumption of hyperinterpolation on the sphere," which has been published in [Journal of Complexity]([https://arxiv.org/abs/2202.13691](https://www.sciencedirect.com/science/article/pii/S0885064X23000584)) (vol. 80, paper no. 101789) in 2024.

# Reproducing figures
* As the names of these M.files indicate, each M.file corresponding to the reproduction of the indicated figure in the paper.
* For illustrating hyperinterpolants rather than error curves, please visit our another repository [MZHyper](https://github.com/HaoNingWu/MZHyper/). More point sets are available in this repository.  
* Please download the sphere_approx_toolbox_v3.0 (available [here](https://1drv.ms/u/s!AmzdJkQhNBOrhlhZ7TNzdUYOb7X1?e=rfGNGn)) and add it onto path before running the codes.
* Please go to /sphere_approx_toolbox_v3.0/utilities/ and change the bold part of the path '**/Users/haoningwu/Documents/MATLAB/BypaHyper**/sphere_approx_toolbox_v3.0/data/xx' in **loadMD.m**, **loadME.m**, and **loadStd.m** to your own path storing the sphere_approx_toolbox_v3.0. Otherwise, MATLAB would report error:
  >The file '/Users/haoningwu/Documents/MATLAB/BypaHyper/sphere_approx_toolbox_v3.0/data/xx/xxxx' could not be opened because: No such file or
directory
