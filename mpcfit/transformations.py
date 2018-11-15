#/mpcfit/mpcfit/transformations.py
"""

--------------------------------------------------------------
 
 Oct 2018
 
 Payne

 Functions for transformation between coordinate frames
 
 To include covariance transformation capabilities 
 
 See Bernstein's orbit/transforms.c code for comparison
 
 See B&K (2000) and/or Farnocchia et al (2015) for math

 --------------------------------------------------------------
 
"""

#Import third-party packages
# --------------------------------------------------------------
import numpy as np
import collections

#Import neighboring packages
#--------------------------------------------------------------
#import mpcutilities.phys_const as PHYS


#Fitting functions & classes
#--------------------------------------------------------------


class TRANSFORMS:
    '''
    To hold functions for transformation between coordinate frames
        
    '''
    
    
    def __init__(self, ):
        '''
        ...
        '''
        self.inputCov =np.array( [ [1,2,3,4], [5,6,7,8], [9,10,11,12] ,[13,14,15,16] ])
        self.partial  =np.array( [ [1,2,], [3,4], [5,6] , [7,8]])

    def covarianceRemap(self, inputCov, partial):
        '''
        Remap covariance matrix from one basis to another, using
        the partial deriv matrix 
        
        \Gamma_{out} = A \Gamma_{in} A^T
        A = \frac{\partial X_{out}}{\partial X_{in}}
        X_{in}, X_{out} = input & output coordinate systems
        
        Compare to Bernstein orbit/orbfit1.c/covar_map
        
        kin = dimension on input, kout = dimension on output.
                
        '''

        # Need to have numpy arrays / matricees
        assert isinstance(inputCov, np.ndarray)
        assert isinstance(partial, np.ndarray)
        
        # Do the matrix multiplication in the correct order
        print("inputCov", inputCov)
        print("partial", partial)
        pT = np.transpose(partial)
        print("pT",pT )
        tmp = np.matmul(inputCov, pT)
        print("tmp", tmp)
        outputCov = np.matmul(partial, tmp)
        print("outputCov", outputCov)
        
        
        # Check on mat-mult calc:
        '''
        for i in range(
        for (i=1; i<=kout; i++)
        for (j=1; j<=kout; j++) {
                sum = 0.;
                for (m=1; m<=kin; m++)
                for (n=1; n<=kin; n++)
                        sum += derivs[i][m]*derivs[j][n]*covar_in[m][n];
                        covar_out[i][j] = sum;
                }
        '''

        
        return outputCov


