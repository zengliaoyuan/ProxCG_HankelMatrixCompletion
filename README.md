# ProxCG_HankelMatrixCompletion
This is a Matlab package that implements a single-loop proximal-conditional-gradient penalty method for solving Hankel matrix completion problems with Laplacian noise. More details about the application can be found in our paper: [A single-loop proximal-conditional-gradient penalty method](https://arxiv.org/abs/2409.14957).


### The Matlab source codes are:

**proxCG_L1.m**: the main function of the proximal-conditional-gradient type method (proxCG) for solving the Hankel MC problem with Laplacian noise  <br />

  - **Hkxad_fast.m** : a subroutine used in proxCG_L1 for computing the adjoint of the Hankel operator $H^*(Y)$; the code is a part of the HSGD.m from [https://github.com/caesarcai/HSGD](https://github.com/caesarcai/HSGD) <br />
  - **yHx_hdl.m** : a subroutine used in proxCG_L1; the function handle of $Y-H(x)$ with $Y=UDV^T$
    - **fhmvmultiply_1D.m** : a subroutine used in yHx_hdl.m; a fast multiplication of a Hankel matrix $H(x)$ and a vector $w$; downloaded from [https://github.com/caesarcai/HSGD](https://github.com/caesarcai/HSGD) <br />
  - **update_svd_thin.m** : a subroutine used in proxCG_L1.m; update the svd representation by svd_rank1_update_qr.m and truncate the singular values by thinSVD.m <br />
    - **svd_rank1_update_qr.m** : a subroutine used in update_svd_thin.m; a fast svd of a rank-1 updated matrix by QR decomposition <br />
    - **thinSVD.m** : a subrountine used in update_svd_thin.m <br />



**completeADM.m** : the main function of ADMM <br />

  - **Hkx.m** and **Hkxad.m**: subrountines used in completeADM.m; compute the Hankel operator $H(x)$ and the adjoint of Hankel $H^*(Y)$, respectively; downloaded from [https://www.polyu.edu.hk/ama/profile/pong/Hankel_Final_Codes/](https://www.polyu.edu.hk/ama/profile/pong/Hankel_Final_Codes/) <br />
  - **project_sigVal_to_simplex.m** : a subroutine used in completeADM.m; project the singular values vector (computed by svd) to the simplex <br />


**runProxCG.m** : the runcode comparing proxCG and ADMM <br />

  - **generate_signal.m** : a subroutine used in runProxCG.m; generate random Hankel matrix completion problems <br />


