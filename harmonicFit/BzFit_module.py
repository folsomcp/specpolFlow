'''
Module with tools for doing MCMC fitting of Bz data with MCMC
'''

import numpy as np
import scipy
import emcee

class data:
    '''
    A helper class to hold a set of values (and associated error bars) as a function of JD. 
    '''
    
    def __init__(self, jd, val, val_err):
        self.jd = jd
        self.val = val
        self.val_err = val_err
    
    def __getitem__(self, key):
        """
        Returns an data object with only the values at the specified index(s).

        :param key: the index or slice being checked

        :rtype: data
        """
        jd_s = self.jd[key]
        val_s = self.val[key]
        val_err_s = self.val_err[key]
        
        slice = data(jd_s, val_s, val_err_s)
        return slice

    def __setitem__(self, key, newval):
        """
        Sets all values of the bzdata object at the specified location equal to the
        input object's values.

        :param key: the index or slice being overwritten
        :param newval: LSD object used to replace the overwritten values
        """
        if not(isinstance(newval, data)):
            raise TypeError()
        else:
            self.jd[key] = newval.jd
            self.val[key] = newval.val
            self.val_err[key] = newval.val_err

    def foldrot(self, P, jd0):
        '''
        Returns the Rotation.Phase that corresponds to self.jd,
        for the specified period and JD0.

        Note that Rotation can be negative, if JD0 is after the JD. 
        '''
        return (self.jd - jd0)/P
    
    def foldphase(self, P, jd0):
        '''
        Return the phase that corresponds to the self.jd,
        for the specified period and JD0. 
        '''
        phi = self.foldrot(P, jd0)
        phi = phi - np.fix(phi)
        phi[phi<0] = 1+phi[phi<0]

        return phi

class dipole_model:
    '''
    Helper function to define a Bpole model for Bz
    '''

    def __init__(self, i, beta, Bpole, u):
        self.i=i
        self.beta=beta
        self.Bpole=Bpole
        self.u=u

        const = self.Bpole/20*(15+u)/(3-u)
        A = const*np.sin(self.i)*np.sin(self.beta)
        mean = const*np.cos(self.i)*np.cos(self.beta)
        self.cos_model = cos_model(A, mean, 0.0)

    def to_array(self):
        '''
        Helper function to return an array of parameters
        to use with minimizing functions
        '''
        return [self.i, self.beta, self.Bpole, self.u]
    
    def to_dict(self):
        d = {
            'i':self.i,
            'beta':self.beta,
            'Bpole':self.Bpole,
            'u':self.u,
        }
        return d

    def get_model(self, phi):
        '''
        Function to return the Bz value for a given phi for this dipolar model

        :param phi: The rotational phase in radiants
        '''
        Bz = self.cos_model.get_model(phi)
        return Bz

    def get_cosalpha(self, phi):
        '''
        Function to return cos(alpha), where alpha is the angle between the magnetic axis and the line of sight.
        cos(alpha) = cos(beta)cos(i) + sin(beta)sin(i)cos(phi)

        :param phi: The rotational phase in radiants
        '''
        cos_alpha = np.cos(self.i)*np.cos(self.beta) + np.sin(self.i)*np.sin(self.beta)*np.cos(phi)
        return cos_alpha
    
def dipole_model_from_df(df):
    '''
    Create a dipole_model object from a dictionary
     or from a single row of a pandas dataframe. 

    :param 
    '''
    return dipole_model(df["i"], df["beta"], df["Bpole"], df["u"])

class cos_model:
    '''
    Helper function to define a cos model for Bz:
    Model = A*cos(phi+phi0) + mean
    '''

    def __init__(self, A, mean, phi0):
        self.A=A
        self.mean=mean
        self.phi0=phi0

        self.max = mean + A
        self.min = mean - A
        self.r = self.min / self.max

    def to_array(self):
        '''
        Helper function to return an array of parameters
        to use with minimizing functions
        '''
        return [self.A, self.mean, self.phi0]
    
    def to_dict(self):
        d = {
            'A':self.A,
            'mean':self.mean,
            'phi0':self.phi0,
            'max':self.max,
            'min':self.min,
            'r':self.r
        }
        return d

    def get_model(self, phi):
        '''
        Function to return value of the cos_model for a given phi
        '''
        return self.A*np.cos(phi+self.phi0) + self.mean


    def get_i(self, beta):
        '''
        Function to return the value of i implied by given value(s) of beta
        '''

        # using the arctan2 to avoir dividing by zero
        i = np.arctan2( (1-self.r), (1+self.r)*np.tan(beta) )
        
        return i
    
    def get_beta(self, i):
        '''
        Function to return the value of beta implied by given value(s) of i
        '''

        # using the arctan2 to avoir dividing by zero
        beta = np.arctan2( (1-self.r), (1+self.r)*np.tan(i) )
        
        return beta
    
def cos_model_from_df(df):
    '''
    Create a cos_model object from a dictionary object or
    from a single row of a pandas dataframe.
    Note: only A, mean, and phi0 are used, because min, max, and r are computed 
    at the object initialization. 
    '''
    return cos_model(df["A"], df["mean"], df["phi0"])

class harmonic_model:
    '''
    Helper function to define a hamonics model:
    Model = a0 + a1cos(phi) + b1sin(phi) + a2 cos(2phi) + b2sin(2phi) + ...
    '''

    def __init__(self, coeff):
        '''
        Constructor for the harmonic_model. 

        :param coeff: an array with the coefficients, in this order: a0, a1, b1, a2, b2, ...
        '''
        self.coeff = coeff

        self.n_coeff = len(coeff)
        self.order = int(self.n_coeff / 2)
        ia = [0] + list(range(1,self.n_coeff,2))
        self.a = self.coeff[ia]
        self.b = self.coeff[range(2,self.n_coeff,2)]

    def get_model_loop(self,phi):
        '''
        Function to return value of the harmonic_model for a given phi
        This is a testing funciton written for clarity instead of efficiency
        to check the matrix construction for the inversion. 
        '''
        model = np.zeros(phi.size)
        for i,a in enumerate(self.a):
            model = model + a*np.cos(i*phi)
        for i,b in enumerate(self.b):
            model = model + b*np.sin((i+1)*phi)
        return model
    
    def get_model_matrix(self,phi):
        '''
        Funciton to get the model for a given phi array
        '''

        MM = self.get_phi_matrix(phi)
        model = MM @ self.coeff

        return model

    
    def get_phi_matrix(self, phi_arr):
        '''
        Helper function to get the matrix that is used to fit the harmonics model

        :param phi_arr: Array of phases for the data (in radians)
        '''

        ia_orders = np.arange(1, self.a.size)
        ib_orders = np.arange(1, self.b.size+1)
        #print('ia: ',ia_orders)
        #print('ib: ',ib_orders)
        
        MM = np.zeros((phi_arr.size, self.n_coeff))
        MM[:,0] = 1.0
        
        #phi_arr = np.ones(phi_arr.size) # debugging

        #MM[:,range(1,self.n_coeff,2)] = 99 # debugging
        #MM[:,range(1,self.n_coeff,2)] = np.outer(phi_arr, ia_orders) # debugging
        MM[:,range(1,self.n_coeff,2)] = np.cos(np.outer(phi_arr, ia_orders))
        
        #MM[:,range(2,self.n_coeff,2)] = 55 # debugging
        #MM[:,range(2,self.n_coeff,2)] = np.outer(phi_arr, ib_orders) # debugging
        MM[:,range(2,self.n_coeff,2)] = np.sin(np.outer(phi_arr, ib_orders))
        
        return MM

def harmonic_fit(data, P, jd0, n):
    '''
    Compute the coefficient of a harmonics fit 
    using least-square fitting. 

    :param data: A time_series object.
    :param P: rotation period to fold the data
    :param jd0: the zeroth of the ephemeris
    :param n: the number of coefficients to use [a0,a1,b1,a2,b2,...]
    '''
    ####
    # The solution 
    # Val = a0 + a1 cos(phi) + b1 sin(phi) + a2 cos(2phi) + b2 sin(2phi)+..
    #
    # We can change this to a matrix problem:
    # |val1|   | 1  cos(phi1) sin(phi1) cos(2*phi1) sin(2*phi1) ...| |a0|
    # |val2| = | 1  cos(phi2) sin(phi2) cos(2*phi2) sin(2*phi2) ...| |a1|
    # |val3|   | 1  cos(phi3) sin(phi3) cos(2*phi3) sin(2*phi3) ...| |b1|
    #
    # V = MM * C
    #
    # The least-square solution is 
    # (M^T S^2 M)^-1 (M^T S^2 V)
    # where S is a diagonal matrix with 1/val_err
    # We can save some memory by splitting the S^2 into two:
    # (M^T S) (S M).
    # With a diagonal matrix, 
    # a1 = np.array([1,2,3])
    # a2 = np.array([[1,0,0],[0,2,0],[0,0,3]])
    # np.dot(a2,b) gives the same results as b*a1[:, np.newaxis]

    # get the phase of the observations
    # and change into radians
    phi = data.foldphase(P, jd0) * 2*np.pi

    # Create a harmonic_model object with zero coefficients
    empty_model = harmonic_model(np.zeros(n))

    # Use the class functions to get the design matrix
    # (doesn't depends on the coefficients)
    MM = empty_model.get_phi_matrix(phi)
    #print(MM)

    MM /= data.val_err[:, np.newaxis]
    #print(MM)
    #Evaluate the solution to the linear least squares problem
    valsOverSig = data.val/data.val_err
    #fitCoeff = np.linalg.inv(MM.T @ MM) @ (MM.T @ valsOverSig)
    covar = np.linalg.inv(MM.T @ MM)
    auto = (MM.T @ valsOverSig)
    fitCoeff = covar @ auto
    fitErr = np.sqrt(np.diag(covar))

    #Get the reduced chi^2 for this optimal fit
    chi2 = (valsOverSig - MM@fitCoeff).T @ (valsOverSig - MM@fitCoeff)
    chi2nu = chi2/(data.val.size - fitCoeff.size)
    
    return harmonic_model(fitCoeff), fitErr, chi2, chi2nu    



def get_logLH(model, data, P, jd0, b=1):
    '''
    Function that calculates the natural logarith of the likelihood 
    between a data object and a model object, 
    with the data folder to the given Period and JD0. 
    The parameter `b` encodes a noise scaling factor into the likelihood calculation
    such that sigma_i = error_bar_i / d to account for potential under/over estimation
    of the uncertainty, and for features in the data the model cannot explain. 
    '''
    phi = data.foldphase(P, jd0)*2*np.pi
    model_values = model.get_model(phi)
    logLH = np.sum( np.log(b/(data.val_err*(2*np.pi)**0.5)) - 0.5*b**2/data.val_err**2 * (data.val-model_values)**2 )
    return logLH

def minimize(data, P, jd0, model_func, theta0=[2, 2, 0]):

    def function_to_minimize(theta, data, P, jd0):
        model = model_func(*theta)
        return -1*get_logLH(model, data, P, jd0)
    
    #theta0 = [2, 2, 0]
    bounds = ((0,10),(-10,10),(-np.pi/2,np.pi/2))

    print(function_to_minimize(theta0, data, P, jd0))   

    result = scipy.optimize.minimize(function_to_minimize, 
                    theta0, args=(data, P, jd0), 
                    bounds=bounds,
                    method='Nelder-Mead',
                    options={'return_all':True,'disp':True,'adaptive':True,'xatol':1e-7,'fatol':1e-7})

    print(result.message)

    return model_func(*result.x)


def log_prior_flat(theta, bounds):
    
    isIn = True
    for i in range(0,theta.size):
        # is the ith parameter in the bounds?
        # if yes, then isIn remains True. 
        # if not, then isIn becomes False
        isIn = isIn and bounds[i][0] < theta[i] < bounds[i][1]

    if isIn == True:
        return 0.0
    return -np.inf

def log_probability(theta, data, P, jd0, bounds):
    # get the flat prior in the bounds
    lp = log_prior_flat(theta, bounds)
    if not np.isfinite(lp):
        return -np.inf
    
    model = cos_model(*theta)
    return lp + get_logLH(model, data, P, jd0, b=1)

def set_mcmc(data, P, jd0, filename='test.h5'):

    nwalkers = 32
    ndim = 3
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    pos = np.zeros((nwalkers,ndim))

    #bounds = ((0,10),(-10,10),(0.0,2*np.pi),(0.1,2))
    bounds = ((0,10),(-10,10),(-np.pi/2,np.pi/2))

    # setting up the initial positions of walkers
    # to be uniformly distributed in the bounds
    for i in range(0,ndim):
        pos[:,i] = np.random.uniform(bounds[i][0], bounds[i][1],nwalkers)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(data, P, jd0, bounds), backend=backend)
    return sampler, pos
    #sampler.run_mcmc(pos, 20000, progress=True)
    
  






