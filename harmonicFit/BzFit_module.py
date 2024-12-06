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

class cos_model:
    '''
    Helper function to define a cos model for Bz:
    Bz = A*cos(phi+phi0) + mean
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
        Function to return the Bz value for a given phi for this cos model
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
    Create a cos_model object from a single row of a dataframe.
    Note: only A, mean, and phi0 are used, because min, max, and r are computed 
    at the object initialization. 
    '''
    return cos_model(df["A"], df["mean"], df["phi0"])

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

def minimize(data, P, jd0, theta0=[2, 2, 0]):

    def function_to_minimize(theta, data, P, jd0):
        model = cos_model(*theta)
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

    return cos_model(*result.x)


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
    
  






