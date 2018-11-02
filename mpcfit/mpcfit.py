# /mpcfit/mpcfit/mpcfit.py
"""
# --------------------------------------------------------------
# Oct 2018
# Payne
#
# Generalized class(es) for fitting
# Currently just sketching out framework
# 
# Want to be able to fit using various orbit-advancement methods 
# - Which will be called from mpcadvancer
# (i) keplerian; (ii) n-body; 
# 
# Want to be able to do single/simple "best-fit" ...
# ... as well as partial-derivatives / tangent-eqns
#
#
# --------------------------------------------------------------
"""


# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import collections

# Import neighboring packages
# --------------------------------------------------------------
#import mpcutilities.phys_const as PHYS


# Fitting functions & classes
# --------------------------------------------------------------


class ORBIT_FIT:
    '''
    To hold full results of orbit fit
    
    MJP:    At present am using this as a scratch-pad for
            ideas as I wade through Milani / Farnocchia / etc
            
        
    '''
    
    
    def __init__(self, ):
        '''
        Quantities which *may* be needed ...
        ... Following Farnocchia for now ...
        '''
        self.x_Vec      = ...       # Vector of orbital parameters
        self.xBest_Vec  = ...       # Vector of best-fit orbital parameters (nominal soln)
        self.nu_Vec     = ...       # Vector of residuals (O-C)
        self.W_Mat      = ...       # Weight matrix (from observational uncertainties)
        self.Q_Mat      = ...       # Cost function (Q = \nu^T W \nu ) <<-- This is chiSq
        self.B_Mat      = ...       # Design matrix (B = \partial \nu / \partial x
        self.C_Mat      = ...       # Normal Matrix (C = B^T W B )
        self.G_Mat      = ...       # Covariance matrix, \Gamma = C^{-1}


class FIT:
    '''
    Provides access to methods to fit OBSERVATIONS-class object
    with an orbit of some kind
    
    Will presumably utilize MPCAdvancer methods for the 
    orbit propagation part 
    
    Need to think about uniform-versus-nonuniform target-times in MPCAdvancer:
    - fitting observations will require non-uniform
    - for n-body, suspect passing to cheby/poly-func will be easiest

    '''


    def __init__(self, observations):
        '''
        Must initialize with a valid OBSERVATIONS object
        
        Parameters
        ----------
        observations   :   OBSERVATIONS-class
        
        '''
        
        # Add assertions to assure contents
        #  - after development-phase, could disable assertions
        # assert hasattr(observations, '???')






    def get_residuals(x_Vec , t_Vec, d_Vec, params):
        '''
        Calculation of residual vector:
            (Observed Data) - (Calculated Values)
            
        This is *NOT* iterated.
        Will call MPCAdvancer function to get heliocentric
        Will call MPCSky function to get sky-plane
        
        Parameters
        ----------
        x_Vec       :       ndarray
            Vector of function parameters (orbit spec)
        t_Vec       :       ndarray
            Vector of times of observations
            Independent variable
        d_Vec       :       ndarray
            Vector of observations
            Dependent variable(s)
        params      :       dict
            A means to pass various specifiction parameters. 
            E.g. 
             - The method of orbital advance/evaluation
             - ...
        
        Returns
        -------
        nu_Vec      :       ndarray
            Vector of residuals (O-C)
        
        Examples
        --------
        >>> ...
        
        '''
    
        # Use MPCAdvancer to get heliocentric cartesian positions at times t_Vec
        
        # Convert heliocentric cartesians to topocentric RA,Dec
        
        # Return vector of residuals
        pass
        

    def get_best_fit_LAZY():
        '''
        Iterative method to find best-fit orbit parameters 
        for a set of observations
        
        I think that a lazy / black-box approach is to use 
        scipy.optimize.least_squares /
        scipy.optimize.minimize / 
        scipy.optimize.curve_fit
         - Experimentation & preparation is in the notebook experimentation_with_leastsquares.ipynb
        
        Parameters
        ----------
        
        Returns
        -------
        
        Examples
        --------
        >>> ...
        
        '''
        # Conceptually we could get a lest-squares fit as follows ...
        # (this is not yet using jacobian approach) 
        res_lsq = least_squares(get_residuals, x_VecInit, args=(t_Vec, d_Vec, params))
        
        pass


    def get_best_fit_EXPLICIT():
        '''
        Iterative method to find best-fit orbit parameters 
        for a set of observations
        
        Explicit evaluation of design matrix, 
        application of differential correction, etc,
        following Milani / Farnocchia / etc
        
        Parameters
        ----------
        
        Returns
        -------
        
        Examples
        --------
        >>> ...
        
        '''
        pass
        


    def function_to_do_outlier_rejection_and-or_deal_with_non_convergence_and_hence_do_multiple_iterations_of_get_best_fit_XXX():
        '''
        Obviously the name will change
        
        Parameters
        ----------
        
        Returns
        -------
        
        Examples
        --------
        >>> ...
        
        '''
        pass




















