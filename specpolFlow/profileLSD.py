## @module profileLSD.py
"""
Tools for generating and manipulating LSD profiles.  Includes calculating
longitudinal magnetic fields and radial velocities.  Also includes a wrapper
around LSDpy.
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as specialf
from scipy.optimize import curve_fit

from .obsSpec import Spectrum

###################################
###################################
class LSD:
    """
    Holds the LSD profile data.

    This usually includes numpy arrays:
    
    * vel - velocity grid for the LSD profile
    * specI - the Stokes I profile
    * specSigI - the uncertainties on Stokes I
    * specV - the polarization profile, most often Stokes V
    * specSigI - the uncertainties on the polarization profile
    * specN1 - the Null1 profile
    * specSigN1 - the uncertainties on specN1
    * specN2 - the Null2 profile, if it exists, otherwise zeros
    * specSigN2 - the uncertainties on Null2 if they exist
    * header - an optional text header for the LSD profile file

    And integers:
    
    * numParam - The number of Stokes I V & Null profiles used (usually 1, 3, or 4)
    * npix - the number of pixels (points in velocity) in the LSD profile
    """
    def __init__(self, vel, specI, specSigI, specV=None, specSigV=None,
                 specN1=None, specSigN1=None, specN2=None, specSigN2=None,
                 header=None):
        """
        Initialize an LSD profile, using data from an existing profile. 

        :param vel: velocity grid for the LSD profile
        :param specI: the Stokes I profile
        :param specSigI: the uncertainties on Stokes I
        :param specV: the polarization profile, usually Stokes V
        :param specSigV: the uncertainties on the polarization profile 
        :param specN1: the Null1 profile
        :param specSigN1: the uncertainties on Null1
        :param specN2: the Null2 profile, if it exists, otherwise it is all 0s (optional)
        :param specSigN2: the uncertainties on Null2 if there are any, otherwise all 0s (optional)
        :param header: the header of the LSD profile file, if it exists (optional) 
        """
        self.vel = vel
        self.specI = specI
        self.specSigI = specSigI
        self.npix = vel.size
        self.numParam = 1
        
        if(specV is not None and specSigV is not None):
            self.specV = specV
            self.specSigV = specSigV
            self.numParam = 2
        else:
            self.specV = np.zeros(self.npix)
            self.specSigV = np.zeros(self.npix)
            
        if(specN1 is not None and specSigN1 is not None):
            self.specN1 = specN1
            self.specSigN1 = specSigN1
            self.numParam = 3
        else:
            self.specN1 = np.zeros(self.npix)
            self.specSigN1 = np.zeros(self.npix)
        
        if(specN2 is not None and specSigN2 is not None):
            self.specN2 = specN2
            self.specSigN2 = specSigN2
            self.numParam = 4
        else:
            self.specN2 = np.zeros(self.npix)
            self.specSigN2 = np.zeros(self.npix)
        
        if(header is not None):
            self.header = header
        else:
            self.header = ""
    
    def __len__(self):
        """
        Return the length of each array in an LSD profile.
        They should all be the same length, so it just returns the length of
        the velocity array. 

        :rtype: int
        """
        return len(self.vel)

    def __getitem__(self, key):
        """
        Returns an LSD with only the values at the specified index(s).

        :param key: the index or slice being checked

        :rtype: LSD
        """
        vel_s = self.vel[key]
        specI_s = self.specI[key]
        specSigI_s = self.specSigI[key]
        specV_s = self.specV[key]
        specSigV_s = self.specSigV[key]
        specN1_s = self.specN1[key]
        specSigN1_s = self.specSigN1[key]
        specN2_s = self.specN2[key]
        specSigN2_s = self.specSigN2[key]
        slice_prof = LSD(vel_s, specI_s, specSigI_s, specV_s, specSigV_s,
                         specN1_s, specSigN1_s, specN2_s, specSigN2_s, self.header)
        slice_prof.numParam = self.numParam
        return slice_prof

    def __setitem__(self, key, newval):
        """
        Sets all values of the LSD at the specified location equal to the
        input profile's values.

        :param key: the index or slice being overwritten
        :param newval: LSD object used to replace the overwritten values
        """
        if not(isinstance(newval, LSD)):
            raise TypeError()
        else:
            self.vel[key] = newval.vel
            self.specI[key] = newval.specI
            self.specSigI[key] = newval.specSigI
            self.specV[key] = newval.specV
            self.specSigV[key] = newval.specSigV
            self.specN1[key] = newval.specN1
            self.specSigN1[key] = newval.specSigN1
            self.specN2[key] = newval.specN2
            self.specSigN2[key] = newval.specSigN2

    def save(self, fname):
        """
        Save the LSD profile to a file.
        
        This saves to a text file in Donati's format.
        
        :param fname: the name of the file the LSD profile to save to.
        """
        
        oFile = open(fname, 'w')
        if self.header != None:
            if self.header == "":
                #For blank headers include some placeholder text
                oFile.write('#LSD profile\n')
            else:
                oFile.write(self.header)
                #Make sure there is a line break after the first header text
                if self.header[-1] != '\n': oFile.write('\n')
            if self.numParam <= 3:
                oFile.write(' {:d} 6\n'.format(self.npix))
            elif self.numParam == 4:  #if Null1 & Null2 exist
                oFile.write(' {:d} 8\n'.format(self.npix))

        if self.numParam <= 3:
            for i in range(self.npix):
                oFile.write('{:>12.6f} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e}\n'.format(
                    self.vel[i], self.specI[i], self.specSigI[i], self.specV[i],
                    self.specSigV[i], self.specN1[i], self.specSigN1[i]))
        elif self.numParam == 4:  #if Null1 & Null2 exist
            for i in range(self.npix):
                oFile.write('{:>12.6f} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e}\n'.format(
                    self.vel[i], self.specI[i], self.specSigI[i], self.specV[i],
                    self.specSigV[i], self.specN1[i], self.specSigN1[i],
                    self.specN2[i], self.specSigN2[i]))
        oFile.close()
        return

    def norm(self, normValue):
        """
        Return a renormalize an LSD profile, divide the I, V, and null 
        profiles by a value.
        
        :param normValue: the value to renormalize (divide) the LSD profile by
        :rtype: LSD
        """
 
        new = LSD(self.vel, 
                  self.specI/normValue, self.specSigI/normValue, 
                  self.specV/normValue, self.specSigV/normValue,
                  self.specN1/normValue, self.specSigN1/normValue, 
                  self.specN2/normValue, self.specSigN2/normValue, 
                  self.header)
        new.numParam = self.numParam

        return new

    def vshift(self, velShift):
        """
        Return a LSD profile with a shift to the velocities of the LSD profile
        
        :param velShift: the change in velocity to be added
        :rtype: LSD
        """
        new = LSD(self.vel-velShift, 
                  self.specI, self.specSigI, 
                  self.specV, self.specSigV,
                  self.specN1, self.specSigN1, 
                  self.specN2, self.specSigN2, 
                  self.header)
        new.numParam = self.numParam

        return new

    def scale(self, scale_int, scale_pol):
        """
        Return a LSD profile with rescaled amplitudes of the LSD profile
        (see also set_weights())
        
        :param scale_int: scale the intensity profile by this
        :param scale_pol: scale the polarization and null profiles by this
        :rtype: LSD
        """

        new = LSD(self.vel, 
                  1.0 - ((1.0-self.specI)*scale_int), self.specSigI*scale_int, 
                  self.specV*scale_pol, self.specSigV*scale_pol,
                  self.specN1*scale_pol, self.specSigN1*scale_pol, 
                  self.specN2*scale_pol, self.specSigN2*scale_pol, 
                  self.header)
        new.numParam = self.numParam

        return new

    def set_weights(self, wint_old, wpol_old, wint_new, wpol_new):
        """
        Return a LSD profile with different mask normalization/scaling weights
        (see also scale())
        
        :param wint_old: The current intensity weight (d)
        :param wpol_old: The current polarization weight (g*d*lambda)
        :param wint_new: The new intensity weight (d)
        :param wpol_new: The new polarization weight (g*d*lambda)
        :rtype: LSD
        """
        return self.scale(wint_new/wint_old, wpol_new/wpol_old)
   
    def plot(self, figsize=(10,10), sameYRange=True, plotZeroLevel=True,
             fig=None, **kwargs):
        """
        Plot the LSD profile

        This generates a matplotlib figure. To display the figure, after
        running this function, use 'import matplotlib.pyplot as plt' and
        'plt.show()'. To save the figure you can use 
        'fig.savefig("fileName.pdf")' on the new figure object.
        
        :param figsize: the size of the figure being created
        :param sameYRange: optionally set the Stokes V and Null plots to have
                           the same y scale
        :param plotZeroLevel: optionally include a dashed line at 0 for
                              the V and Null plots
        :param fig: optionally, a matplotlib figure used for plotting the LSD 
                    profile. By default (or if None) a new figure will be
                    generated.
        :return: a matplotlib figure object, and an axes object,
                 containing the plot. 
        """

        #Chose y-axis limits for the V and N plots so that they can have 
        #the same scale.  Use 1.05 time biggest divergence from zero
        plotVLims = np.max([np.abs(1.05*np.min(self.specV)),
                            np.abs(1.05*np.max(self.specV)),
                            np.abs(1.05*np.min(self.specN1)),
                            np.abs(1.05*np.max(self.specN1)),
                            np.abs(1.05*np.min(self.specN2)),
                            np.abs(1.05*np.max(self.specN2)),
                            1.05*np.max(self.specSigV)])
        
        #Get the number of Stokes parameters to plot
        #(Check the number of non-zero data field)
        #Support just Stokes I; Stokes I, V, Null1; Stokes I, V, Null1, Null2
        npar = self.numParam

        firstPlt = True
        if fig == None:
            fig, ax = plt.subplots(npar, 1, figsize=figsize, sharex=True)
        else:
            ax = fig.axes
            firstPlt = False

        #set up a custom color cycle, with black then the default colors
        if firstPlt:
            prop_cycle = plt.rcParams['axes.prop_cycle']
            colorList = ['#000000'] + prop_cycle.by_key()['color']

        #If there is only 1 plot (only a Stokes I spectrum),
        #wrap the matplotlib axis in a list to make it iterable
        if not (isinstance(ax, np.ndarray) or isinstance(ax, list)): ax = [ax]
        #protect against the number of axes not matching self.numParam
        for i, axi in enumerate(ax):
            if firstPlt: axi.set_prop_cycle(color=colorList)
            if i == len(ax) - 1: #Stokes I (always the last/bottom plot)
                axi.errorbar(self.vel, self.specI, yerr=self.specSigI, 
                             xerr=None, fmt='o', ms=3, ecolor='0.5', **kwargs)
                axi.set_xlim(xmin=np.min(self.vel),xmax=np.max(self.vel))
                axi.set_ylabel('I')
                axi.set_xlabel('Velocity (km/s)')
            elif (i == 0) and (npar > 1): #Stokes V
                axi.errorbar(self.vel, self.specV, yerr=self.specSigV, 
                             xerr=None, fmt='o', ms=3, ecolor='0.5', **kwargs)
                if plotZeroLevel: axi.axhline(y=0.0, ls='--', c='k', alpha=0.4)
                if sameYRange: axi.set_ylim(ymin=-plotVLims,ymax=plotVLims)
                axi.set_ylabel('V')
            elif (i == 1) and (npar > 2): #Null1
                axi.errorbar(self.vel, self.specN1, yerr=self.specSigN1, 
                             xerr=None, fmt='o', ms=3, ecolor='0.5', **kwargs)
                if plotZeroLevel: axi.axhline(y=0.0, ls='--', c='k', alpha=0.4)
                if sameYRange: axi.set_ylim(ymin=-plotVLims,ymax=plotVLims)
                axi.set_ylabel('N1')
            elif (i == 2) and (npar > 3): #Null2
                axi.errorbar(self.vel, self.specN2, yerr=self.specSigN2, 
                             xerr=None, fmt='o', ms=3, ecolor='0.5', **kwargs)
                if plotZeroLevel: axi.axhline(y=0.0, ls='--', c='k', alpha=0.4)
                if sameYRange: axi.set_ylim(ymin=-plotVLims,ymax=plotVLims)
                axi.set_ylabel('N2')
        
        plt.subplots_adjust(hspace=.0)
        return fig, ax

    def fit_gaussian_rv(self, velrange=None, plotFit=False, fullOutput=False,
                        scaleErrs=False):
        """
        Fit a Gaussian function to this LSD profile, to determine a
        radial velocity.
        
        :param velrange: Range of velocity to fit include in the fit
                         (a tuple or list with 2 elements)
                         If not given, the whole range will be used
        :param plotFit: If True, Plot the fit to the LSD profile with matplotlib
        :param fullOutput: If True, Return all the Gaussian best fit parameters,
                           rather than just the RV and RV error
        :param scaleErrs: If True scale the errors on the fitting parameters by
                          the square root of the reduced chi^2
        :return: Tuple of best fit RV and RV error.
                 If fullOutput is True then also includes continuum & error,
                 amplitude & error, and width & error
        """

        if velrange == None:
            velrange = (self.vel[1], self.vel[-2])
        elif not(isinstance(velrange, list) or isinstance(velrange, tuple)):
            raise TypeError('velrange in fitGaussianRV should be a list or '
                            +'tuple, with two elements')
        indVelUse = (self.vel >= velrange[0]) & (self.vel <= velrange[1])

        #Set the initial estimates for the fitting parameters
        indMinI = np.argmin(self.specI[indVelUse])
        initVel = self.vel[indVelUse][indMinI]
        initCont = 1.0
        initAmpl = 0.1
        initWidth = 10.
        param0 = [initCont, initAmpl, initVel, initWidth]

        #run the fitting routine
        fitVal, covarVal = curve_fit(_gaussProf, self.vel[indVelUse],
                                     self.specI[indVelUse], p0=param0,
                                     sigma=self.specSigI[indVelUse],
                                     absolute_sigma=(not scaleErrs))
        fitSigma = np.sqrt(np.diag(covarVal))
        
        fitCont = fitVal[0]
        fitAmpl = fitVal[1]
        fitVel = fitVal[2]
        fitWidth = fitVal[3]
        fitContErr = fitSigma[0]
        fitAmplErr = fitSigma[1]
        fitVelErr = fitSigma[2]
        fitWidthErr = fitSigma[3]

        #Optionally plot the fit to the observed LSD profile
        if plotFit == True:
            plt.errorbar(self.vel, self.specI,  yerr=self.specSigI, c='grey')
            plt.errorbar(self.vel[indVelUse], self.specI[indVelUse],
                          yerr=self.specSigI[indVelUse], c='k')
            plt.plot(self.vel[indVelUse], _gaussProf(
                self.vel[indVelUse], fitCont, fitAmpl, fitVel, fitWidth), 'b.')
            plotSynVels = np.linspace(velrange[0], velrange[1], 1000)
            plt.plot(plotSynVels, _gaussProf(
                plotSynVels, fitCont, fitAmpl, fitVel, fitWidth), 'b')
            plt.xlabel('Velocity (km/s)')
            plt.ylabel('I')
            plt.show()

        #Optionally return all fit parameters
        if fullOutput == True:
            return (fitVel, fitVelErr, fitCont, fitContErr, fitAmpl, 
                    fitAmplErr, fitWidth, fitWidthErr)
        return fitVel, fitVelErr

    def cog_I(self, Ic=1.0):
        '''
        Helper function to return the center of gravity of Stokes I for
        this LSD profile, for a given continuum level

        :param Ic: the continnum level to use the the COG calculation
                   (float, default=1.0)
        :return: the velocity of the center of gravity
        '''
        nominator = np.trapz(self.vel * (Ic-self.specI), x=self.vel)
        denominator = np.trapz( Ic-self.specI, x=self.vel )
        return nominator/denominator

    def cog_IV(self, Ic=1.0):
        '''
        Helper function to return the center of gravity of (Stokes I)*(Stoke V)
        for this LSD profile, for a given continuum level

        :param Ic: the continnum level to use the the COG calculation
                   (float, default=1.0)
        :return: the velocity of the center of gravity
        '''
        nominator = np.trapz(self.vel * np.abs( self.specV * (Ic-self.specI) ),
                             x=self.vel )
        denominator = np.trapz(np.abs( self.specV * (Ic-self.specI) ),
                               x=self.vel )
        return nominator/denominator

    def cog_V(self):
        '''
        Helper function to return the center of gravity of this Stokes V
        LSD profile
        
        :return: the velocity of the center of gravity
        '''
        nominator = np.trapz(self.vel * np.abs(self.specV), x=self.vel )
        denominator = np.trapz( np.abs(self.specV), x=self.vel )
        return(nominator/denominator)        

    def cog_min(self):
        '''
        Helper function to return the velocity of the minimum of this
        Stokes I LSD profile
        
        :return: the velocity of the Stokes I minumum
        '''
        
        cog_min = self.vel[self.specI.argmin()]
        if cog_min.size > 1:
            cog_min = cog_min[0]
        return cog_min 

    def calc_fap(self):
        '''Helper function that calculates the FAP for Stokes V, null1,
        and null2, for this LSD profile.

        The False Alarm Probability (FAP) is the probability that the observed
        data are consistent with the null hypothesis of no magnetic field.
        In this case, the null hypothesis is a flat line in Stokes V
        (an offset from 0 is allowed to account for continuum polarization).
        The probability is evaluated from chi^2.
        
        If you would like a specific range in velocity, simply slice the LSD
        object beforehand.  Note that the calcBz function also returns the
        FAP inside the spectral line, over the same velocity range as used
        in the Bz calculation. 

        :return: the values FAP V, FAP N1, FAP N2. 
        '''

        approxDOF = (self.npix-1.)

        #'fitting' the flat line (essentially an average weighted by 1/sigma^2)
        contV = np.sum(self.specV/self.specSigV**2) / np.sum(1./self.specSigV**2)
        chi2V = np.sum(((self.specV - contV)/self.specSigV)**2)
        probV = 1.0-specialf.gammainc(approxDOF/2., chi2V/2.)
        #repeat for the Null1 and Null2 profiles (if they exist)
        probN1 = 0.
        probN2 = 0.
        
        if self.numParam > 2:
            contN1 = np.sum(self.specN1/self.specSigN1**2) / np.sum(1./self.specSigN1**2)
            chi2N1 = np.sum(((self.specN1 - contN1)/self.specSigN1)**2)
            probN1 = 1.0-specialf.gammainc(approxDOF/2., chi2N1/2.)
        if self.numParam > 3:
            contN2 = np.sum(self.specN2/self.specSigN2**2) / np.sum(1./self.specSigN2**2)
            chi2N2 = np.sum(((self.specN2 - contN2)/self.specSigN2)**2)
            probN2 = 1.0-specialf.gammainc(approxDOF/2., chi2N2/2.)
        
        return(probV, probN1, probN2)

    def calc_bz(self, cog='I', norm='auto', lambda0=500., geff=1.2,
                velrange=None, bzwidth=None, plot=True):
        '''Calculate the Bz of an LSD profile
        
        :param cog: The value, or calculation method for the center of gravity.
                    The choices are:
                    'I': center of gravity of I,
                    'V': center of gravity of V,
                    'IV': center of gravity of I*V,
                    or a float: a user defined value in km/s.
        :param norm: calculation method for the continuum. The choices are:
                    'auto': the median of I outside of velrange (if defined)
                    or the full range (if velrange is not defined),
                    or float: a user defined value to use for Ic.
        :param lambda0: wavelength of the transition in nanometers (default=500).
                    For an LSD profile, this is the lambda value the LSD
                    profile shape was scaled with.
        :param geff: effective Lande factor of the transition.
                    For an LSD profile, this is the geff value the LSD profile
                    shape was scaled with.
        :param velrange: range of velocity to use for the determination of the
                    line center and the continuum. If not defined, the whole
                    range will be used. If bzwidth is not defined, this range
                    will also be used for the line Bz calculation.
        :param bzwidth: distance from the line center for the Bz calculation.
                    One element = same on each side of line center.
                    Two elements, left and right of line center.
                    Not defined: using velrange.
        :param plot: whether or not a graph is generated, and returned.
        :return: a dictionary with Bz (in G) and FAP calculations,
                 optionally also a matplotlib figure.
        '''
        
        # Velrange is used to identify the position of the line,
        # for calculating the cog, and for calculating the position
        # of the continuum.
        # If Velrange is not defined, it will use the whole range.
        # The range for calculating Bz itself is controlled by bzwidth below
        if velrange != None:
            inside = np.logical_and(self.vel>=velrange[0], self.vel<=velrange[1])
            lsd_in = self[inside]
            lsd_out = self[np.logical_not(inside)]
        else:
            lsd_in=copy.deepcopy(self)
        
        # Check if norm is a string.
        if isinstance(norm, str):
            print('using AUTO method for the normalization')
            if velrange != None:
                print('  using the median of the continuum outside of the line')
                norm_val = np.median(lsd_out.specI)
            else:
                print('  no range in velocity given, using the median of the whole specI to determine continuum')
                norm_val = np.median(lsd_in.specI)
        else:
            norm_val = copy.copy(norm)
            print('using given norm value')

        # Check if cog is a string.
        if isinstance(cog, str):
            # calculate the cog for the chosen method
            if cog == 'I':
                cog_val = lsd_in.cog_I(norm_val)
            elif cog == 'min':
                cog_val = lsd_in.cog_min()
            elif cog == 'IV':
                cog_val = lsd_in.cog_IV(norm_val)
            elif cog == 'V':
                cog_val = lsd_in.cog_V()
            else:
                raise ValueError('calcBz got unrecognized value for cog: {:}'.format(cog))
        else:
            cog_val=copy.copy(cog)
        
        # Now we define the position of the line for the Bz calculation itself.
        # If the keyword bzwidth is defined, we use that range from the chosen cog.
        if bzwidth == None:
            # No bzwidth. using vrange if defined
            if velrange != None:
                print('no bzwidth defined, using velrange to calculate Bz')
                lsd_bz = copy.copy(lsd_in)
                # saving the range for plotting later.
                p_bzwidth = np.copy(velrange)
            else:
                print('no bzwidth nor velrange defined, using full velocity range to calculate Bz')
                p_bzwidth = [self.vel.min(), self.vel.max()]
                lsd_bz = copy.copy(self)
        else:
            # Check whether it is a numpy array
            if isinstance(bzwidth, list) or isinstance(bzwidth, tuple):
                if len(bzwidth) == 1:
                    #print('list with one element')
                    # keeping the actual bz calculation range for plotting later.
                    p_bzwidth = [cog_val-bzwidth, cog_val+bzwidth]
                    lsd_bz = self[ np.logical_and(self.vel >= p_bzwidth[0], self.vel <= p_bzwidth[1]) ]
                elif len(bzwidth) == 2:
                    #print('list with two elements')
                    p_bzwidth = [cog_val-bzwidth[0], cog_val+bzwidth[1]]
                    lsd_bz = self[ np.logical_and(self.vel >= p_bzwidth[0], self.vel <= p_bzwidth[1]) ]
                else:
                    print('bzwidth has too many elements (need one or two)')
                    raise ValueError('bzwidth has too many elements {:} (need 1 or 2)'.format(len(bzwidth)))
            else:
                p_bzwidth = [cog_val-bzwidth, cog_val+bzwidth]
                lsd_bz = self[ np.logical_and(self.vel >= p_bzwidth[0], self.vel <= p_bzwidth[1]) ]

        # Actual calculation of the Bz:
        # Call the integration function for each V, Null1, Null2 parameter
        blv, blvSig = _integrate_bz(lsd_bz.vel, lsd_bz.specV, lsd_bz.specSigV, 
                    geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
        if self.numParam > 2:
            bln1, bln1Sig = _integrate_bz(lsd_bz.vel, lsd_bz.specN1, lsd_bz.specSigN1,
                    geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
        else: bln1, bln1Sig = (0.0, 0.0)
        if self.numParam > 3:
            bln2, bln2Sig = _integrate_bz(lsd_bz.vel, lsd_bz.specN2, lsd_bz.specSigN2,
                    geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
        else: bln2, bln2Sig = (0.0, 0.0)
        
        # Get the FAP in the same range as the one used for Bz
        FAP_V, FAP_N1, FAP_N2 = lsd_bz.calc_fap()
        
        result = {
                'Ic': norm_val,
                'cog': cog_val,
                'Bzwidth min': p_bzwidth[0],
                'Bzwidth max': p_bzwidth[1],
                'V bz (G)': blv,
                'V bz sig (G)': blvSig,
                'V FAP': FAP_V,
                'N1 bz (G)': bln1,
                'N1 bz sig (G)': bln1Sig,
                'N1 FAP': FAP_N1,
                'N2 bz (G)': bln2,
                'N2 bz sig (G)': bln2Sig,
                'N2 FAP': FAP_N2
                }

        if plot:
            fig  = _plot_bz_calc(self, lsd_in, lsd_bz, velrange,
                            p_bzwidth, norm_val, cog_val, cog)
            return result,fig
        else:
            return result

##################
##################

def _gaussProf(x, cont, ampl, center, width):
    """
    Define a simple Gaussian line profile (for fitting with)
    """
    y = cont - ampl*np.exp(-(x-center)**2/(2.*width**2))
    return y


def _integrate_bz(vel, spec, specSig, geff, lambda0, cog_val,
                  specI, specSigI, norm_val):
    '''
    Helper function that is used in the Bz calculations. Internal.
    to integrate one Stokes parameter. 

    :param vel: (numpy array) the velocity array
    :param spec: (numpy array) the associated Stokes parameter array
    :param specSig: (numpy array) the error bar on the Stokes parameter
    :param geff: the effective lande factor of the profile
    :param lambda0: the nominal wavelength of the profile
    :cog_val: the chose center of gravity value (km/s)
    :specI: (numpy array) the associated Stokes I parameter
    :specSigI: (numpy array) the associated error bar on Stokes I
    :return: Bz and its error
    '''
    #Evaluate the integral equation for Bz.  Only intended for use in calcBz.
    
    # This is the constant for the Zeeman splitting
    # Lambda_B = constant * lambda0**2 B
    # lambda_B_constant = e/(4*pi*m_e*c) = 4.668644778304102e-05 in cgs units
    lambda_B_constant = 4.668644778304102e-12 #to match/cancel lambda0 in nm
    cvel = 2.99792458e5 #c in km/s to match/cancel the LSD profile velocities
    
    # set the velocity step for error propagation
    deltav = (vel[1] - vel[0]) # This is in km/s
    
    # Calculation of the integral in the numerator of the Bz function
    # with a trapezoidal numerical integral
    # For the error calculation, we propagate like we would for
    # summation in a numerical integral.
    fnum = np.trapz( (vel - cog_val) * spec, x=vel-cog_val ) #in (km/s)^2
    sfnum = np.sqrt(np.sum( (vel - cog_val )**2 * specSig**2 ) * deltav**2)
    
    # Calculation of the integral in the denominator of the Bz function
    # with a trapezoidal numerical integral
    # For the square error calculation, we propagate like we would for
    # summation in a numerical integral.
    ri0v = np.trapz(norm_val-specI, x=vel ) # in km/s
    si0v = np.sqrt(np.sum(specSigI**2 )* deltav**2) # in km/s
    
    # Make the actual Bz calculation.
    # for the units to work out, lambda0 should be in nm
    bl = -1*fnum / ( ri0v*geff*lambda0*cvel*lambda_B_constant)
    blSig = np.abs(bl * np.sqrt( (sfnum/fnum)**2 + (si0v/ri0v)**2 ))
    return bl, blSig


def _plot_bz_calc(lsd, lsd_in, lsd_bz, velrange, p_bzwidth, norm_val, cog_val, cog):
    """
    Generate a plot showing the center of gravity and integration ranges used
    in the calculation of Bz from and LSD profile.  Called by the calc_bz
    function. Mostly just for internal use.

    :param lsd: the full LSD profile to plot
    :param lsd_in: slice of an LSD profile with the range used for COG calculation
    :param lsd_bz: slice of an LSD profile with the range used for integration of Bz
    :param velrange: the velocity range used for COG calculation
    :param p_bzwidth: the velocity range used for integration of Bz
    :param norm_val: the continuum level used for normalization
    :param cog_val: the final COG used for Bz calculation
    :param cog: the input COG flag/value given by the user
    :return: a matplotlib figure object
    """
    #This function relies on the plot method of the LSD profile class
    fig, ax = lsd.plot(sameYRange=False)
    
    for item in ax:
        if velrange != None:
            item.axvline(x=velrange[0], ls='--', label='velrange')
            item.axvline(x=velrange[1], ls='--')
        item.axvline(x=p_bzwidth[0], ls='dotted', label='bzwidth')
        item.axvline(x=p_bzwidth[1], ls='dotted')
        
    ax[-1].axhline(y=norm_val, ls='--', c='pink', label='Ic')
    
    # for the plot, calculate and display all of the possible methods
    # for calculating the cog.
    for item in ax:
        item.axvline(x=lsd_in.cog_min(), label='cog min I', lw=3, alpha=0.5, c='blue')
        item.axvline(x=lsd_in.cog_I(norm_val), label='cog I',lw=3, alpha=0.5, c='red')
        item.axvline(x=lsd_in.cog_IV(norm_val), label='cog I*V',lw=3, alpha=0.5, c='orange')
        item.axvline(x=lsd_in.cog_V(), label='cog V',lw=3, alpha=0.5, c='green')
        item.axvline(x=cog_val, label='chosen cog: {}'.format(cog), ls='--', c='k')
       
    ax[-1].legend(loc=0)
    
    red = lsd_bz.vel > cog_val
    blue = lsd_bz.vel < cog_val
    ax[0].fill_between(lsd_bz.vel[red], lsd_bz.specV[red], step='mid', color='red')
    ax[0].fill_between(lsd_bz.vel[blue], lsd_bz.specV[blue], step='mid', color='blue')
    if(len(ax) > 2):
        ax[1].fill_between(lsd_bz.vel[red], lsd_bz.specN1[red], step='mid', color='red')
        ax[1].fill_between(lsd_bz.vel[blue], lsd_bz.specN1[blue], step='mid', color='blue')
    if(len(ax) > 3):
        ax[2].fill_between(lsd_bz.vel[red], lsd_bz.specN2[red], step='mid', color='red')
        ax[2].fill_between(lsd_bz.vel[blue], lsd_bz.specN2[blue], step='mid', color='blue')
    
    return fig


def read_lsd(fname):
    """
    Read in an LSD profile and return an instance of the LSD class
    
    The LSD profiles are in Donati's text format.
    However, the two lines of header in Donati's format are optional.
    
    :param fname: the name of the file containing the LSD profile
    :rtype: LSD
    """
    #check if this file has a header
    fcheck = open(fname, 'r')
    head1txt = fcheck.readline()
    head2txt = fcheck.readline()
    head3txt = fcheck.readline()
    fcheck.close()
    head1 = head1txt.split()
    head2 = head2txt.split()
    nskip = 2
    try: #Test if the first two lines look like LSD profile numbers
        float(head1[0])
        if(len(head2) > 2 and len(head1) == len(head2)):
            nskip = 0
    except(ValueError):
        nskip = 2

    header = None
    if(nskip > 0): header = head1txt

    #Check the number of columns in the LSD profile
    ncols = len(head3txt.split())
        
    #Read the profile, skipping the header
    __prof = np.loadtxt(fname, skiprows=nskip, unpack=True)
    
    if ncols == 9:
        vel = __prof[0,:]
        specI = __prof[1,:]
        specSigI = __prof[2,:]
        specV = __prof[3,:]
        specSigV = __prof[4,:]
        specN1 = __prof[5,:]
        specSigN1 = __prof[6,:]
        specN2 = __prof[7,:]
        specSigN2 = __prof[8,:]
        #Check for profiles with placeholder columns of 0 (LSDpy may do this!)
        if (np.all(specN2 == 0.0) and np.all(specSigN2 == 0.0)):
            ncols = 7
        elif (np.all(specV == 0.0) and np.all(specSigV == 0.0)
              and np.all(specN1 == 0.0) and np.all(specSigN1 == 0.0)):
            ncols = 3
        else:
            prof = LSD(vel, specI, specSigI, specV, specSigV,
                            specN1, specSigN1, specN2, specSigN2, header=header)
    if ncols == 7:
        vel = __prof[0,:]
        specI = __prof[1,:]
        specSigI = __prof[2,:]
        specV = __prof[3,:]
        specSigV = __prof[4,:]
        specN1 = __prof[5,:]
        specSigN1 = __prof[6,:]
        #Check for profiles with placeholder columns of 0 (LSDpy may do this!)
        if (np.all(specV == 0.0) and np.all(specSigV == 0.0)
              and np.all(specN1 == 0.0) and np.all(specSigN1 == 0.0)):
            ncols = 3
        else:
            prof = LSD(vel, specI, specSigI, specV, specSigV,
                        specN1, specSigN1, header=header)
    
    if ncols == 3:
        vel = __prof[0,:]
        specI = __prof[1,:]
        specSigI = __prof[2,:]
        prof = LSD(vel, specI, specSigI, header=header)

    #For unsupported formats or numbers of columns
    if ncols != 3 and ncols != 7 and ncols != 9:
        raise ValueError("Read an unexpected number of columns from "
                         +"{:}, can't read as an LSD profile.".format(fname))
    return prof

###################################

def run_lsdpy(obs=None, mask=None, outName='prof.dat',
              velStart=None, velEnd=None, velPixel=None, 
              normDepth=0.2, normLande=1.2, normWave=500.0,
              removeContPol=1, trimMask=1, sigmaClipIter=0,
              sigmaClip=500, interpMode=1, outModelName='',
              fLSDPlotImg=1, fSavePlotImg=0, outPlotImgName='figProf.pdf'):
    """
    Run the LSDpy code and return an LSD object  
    (a convenience wrapper around the lsdpy.main() function)

    This requires LSDpy to be installed and in your Python path.
    (see https://github.com/folsomcp/LSDpy )
    
    Any arguments not specified will be read from the file inlsd.dat.
    The file inlsd.dat is optional, but if the file does not exist and 
    any arguments are 'None', the program will error and halt.
    Some arguments have default values, which will be used if they are not
    explicitly specified and if the inlsd.dat file is missing.
    
    Arguments are: 
    
    :param obs:           name of the input observed spectrum file
    :param mask:          name of the input LSD mask file
    :param outName:       name of the output LSD profile (Default = 'prof.dat')
    :param velStart:      float, starting velocity for the LSD profile (km/s)
    :param velEnd:        float, ending  velocity (km/s)
    :param velPixel:      float, velocity pixel size (km/s)
    :param normDepth:     float, normalizing line depth
    :param normLande:     float, normalizing effective Lande factor
    :param normWave:      float, normalizing wavelength
    :param removeContPol: int, flag for whether continuum polarization is 
                          subtracted from the LSD profile (0=no, 1=yes)
                          (Default = 1)
    :param trimMask:      int, flag for whether very closely spaced lines 
                          should be removed from the line mask (0=no, 1=yes)
                          (Default = 1)
    :param sigmaClipIter: int, number of iterations for sigma clipping, 
                          rejecting possible bad pixels based on the fit to
                          Stokes I. Set to 0 for no sigma clipping.
                          (Default = 0, no sigma clipping)
    :param sigmaClip:     float, if sigma clipping, reject pixels where the
                          observation differs from the model by more than this
                          number of sigma.  Should be a large value so only very
                          bad pixels are rejected.
                          (Default = 500.)
    :param interpMode:    int, mode for interpolating the model on to the
                          observation during LSD 0 = nearest neighbour,
                          1 = linear interpolation.
                          (Default = 1)
    :param outModelName:  name of the file for the output model spectrum 
                          (if saved). If this is '' a model 
                          will be generated but not saved to file.
                          (Default = '')
    :param fLSDPlotImg:   int, flag for whether to plot the LSD profile
                          (using matplotlib) (0=no, 1=yes)
                          (Default = 1)
    :param fSavePlotImg:  int, flag for whether to save the plot of the 
                          LSD profile (0=no, 1=yes)
                          (Default = 0)
    :param outPlotImgName: name of the plotted figure of the LSD profile 
                          (if saved) (Default = 'figProf.pdf')
    :return: an LSD profile object, and a 'Spectrum' object containing the model
                          spectrum (the fit to the observation,
                          i.e. the convolution of the line mask and LSD profile)
    """
    import LSDpy

    vel, sI, sIerr, sV, sVerr, sN1, sN1err, headerTxt, specList = LSDpy.lsd(
        obs, mask, outName, velStart, velEnd, velPixel, 
        normDepth, normLande, normWave, removeContPol, trimMask, 
        sigmaClipIter, sigmaClip, interpMode, 1, outModelName, 
        fLSDPlotImg, fSavePlotImg, outPlotImgName)

    prof = LSD(vel, sI, sIerr, sV, sVerr, sN1, sN1err, header=headerTxt)
    modelSpec = Spectrum(specList[0], specList[1], specList[2], specList[3],
                            np.zeros_like(specList[0]), np.zeros_like(specList[0]))
    return prof, modelSpec

###################################
