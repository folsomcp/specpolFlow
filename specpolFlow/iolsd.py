## @module iolsd.py
# Documentation for iolsd.py
#
# Tools for reading and writing files, related to calculating
# and analyzing LSD profiles.

import numpy as np
import matplotlib.pyplot as plt
import copy
import scipy.special as specialf


##################
##################

def integrate_bz(vel, spec, specSig, geff, lambda0, cog_val, specI, specSigI, norm_val):
    '''
    Helper function that is used in the Bz calculations. Internal.
    to integrate one Stokes parameter. 

    :param vel: (numpy array) the velocity array
    :param spec: (numpy array) the associated Stokes parameter array
    :param specSig: (numpy array) the associated error bar on the Stokes parameter
    :param geff: the effective lande factor of the profile
    :param lambda0: the nominal wavelength of the profile
    :cog_val: the chose center of gravity value (km/s)
    :specI: (numpy array) the associated Stokes I parameter
    :specSigI: (numpy array) the associated error bar on Stokes I
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

def plot_bz_calc(lsd, lsd_in, lsd_bz, velrange, p_bzwidth, norm_val, cog_val, cog):
    """
    Generate a plot showing the center of gravity and integration ranges used in the calculation of Bz from and LSD profile.
    Called by the main calcBz function. Internal

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

    And integers:
    
    * numParam - The number of Stokes I V & Null profiles used (usually 1, 3, or 4)
    * npix - the number of pixels (points in velocity) in the LSD profile
    """
    def __init__(self, vel, specI, specSigI, specV=[], specSigV=[],
                 specN1=[], specSigN1=[], specN2=[], specSigN2=[], header=None):
        """
        Initialize an LSD profile, using data from an existing profile. 

        :param self: LSD being created
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
        
        if(specV != [] and specSigV != []):
            self.specV = specV
            self.specSigV = specSigV
            self.numParam = 2
        else:
            self.specV = np.zeros(self.npix)
            self.specSigV = np.zeros(self.npix)
            
        if(specN1 != [] and specSigN1 != []):
            self.specN1 = specN1
            self.specSigN1 = specSigN1
            self.numParam = 3
        else:
            self.specN1 = np.zeros(self.npix)
            self.specSigN1 = np.zeros(self.npix)
        
        if(specN2 != [] and specSigN2 != []):
            self.specN2 = specN2
            self.specSigN2 = specSigN2
            self.numParam = 4
        else:
            self.specN2 = np.zeros(self.npix)
            self.specSigN2 = np.zeros(self.npix)
        
        if(header != None):
            self.header = header
        else:
            self.header = ""
    
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

    #overloaded functions
    def __len__(self):
        """
        Return the length of each array in an LSD profile. They should all be the same length, so it just returns the length of the velocity array. 

        :param self: LSD whose length is being checked

        :rtype: int
        """
        return len(self.vel)

    def __getitem__(self, key):
        """Overloaded getitem function. Returns an LSD with only the values at the specified index(s).

        :param self: LSD being queried
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
        Overloaded setitem function. Sets all values of the LSD at the specified location equal to the input profile's values.

        :param self: LSD object being edited
        :param key: the index or slice being overwritten
        :param newval: LSD whose values are to replace the overwritten ones
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

    def norm(self, normValue):
        """
        Return a renormalize an LSD profile, divide the I, V, and null profiles by a value.
        
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
        # VERO: old def that changes the original object. 
        # Can be removed in cleanup later.
        #self.vel = self.vel + velShift
        #return self
        new = LSD(self.vel-velShift, 
                        self.specI, self.specSigI, 
                        self.specV, self.specSigV,
                        self.specN1, self.specSigN1, 
                        self.specN2, self.specSigN2, 
                        self.header)
        new.numParam = self.numParam

        return new

    def scale(self, scale_int, scale_pol):
        '''Return a LSD profile with rescaled amplitudes of the LSD profile (see also set_weights())
        
        :param scale_int: scale the intensity profile by this
        :param scale_pol: scale the polarization and null profiles by this
        :rtype: lsd object
        '''

        new = LSD(self.vel, 
                        1.0 - ((1.0-self.specI) * scale_int), self.specSigI * scale_int, 
                        self.specV*scale_pol, self.specSigV*scale_pol,
                        self.specN1*scale_pol, self.specSigN1*scale_pol, 
                        self.specN2*scale_pol, self.specSigN2*scale_pol, 
                        self.header)
        new.numParam = self.numParam

        return new

    def set_weights(self, wint_old, wpol_old, wint_new, wpol_new):
        '''Change the weight of the LSD profile (see also scale())
        
        :param wint_old: The current intensity weight (d)
        :param wpol_old: The current polarization weight (g*d*lambda)
        :param wint_new: The new intensity weight (d)
        :param wpol_new: The new polarization weight (g*d*lambda)
        :rtype: lsd object
        '''
        return(self.scale(wint_new/wint_old, wpol_new/wpol_old))
   
    def plot(self, figsize=(10,10), sameYRange=True, plotZeroLevel=True, **kwargs):
        '''Plot the LSD profile
        
        :param self: the lsd object to plot 
        :param figsize: the size of the figure being created
        :param sameYRange: optionally set the Stokes V and Null plots to have the same y scale
        :param plotZeroLevel: optionally include a dashed line at 0 for the V and Null plots.
        :rtype: returns an fig, ax matplotlib container. 
        '''

        #Chose y-axis limits for the V and N plots so that they can have the same scale
        #use 1.05 time biggest divergence from zero
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
        nplots = self.numParam
        
        fig, ax = plt.subplots(nplots,1,figsize=figsize,sharex=True)
        #If there is only 1 plot (only a Stokes I spectrum),
        #wrap the matplotlib axis in a list to make it iterable
        if not isinstance(ax, np.ndarray): ax = [ax]
        for i, axi in enumerate(ax):
            if i == nplots - 1: #Stokes I
                axi.errorbar(self.vel, self.specI, yerr=self.specSigI, xerr=None,
                            fmt='o', ms=3,ecolor='0.5',c='k',**kwargs)
                axi.set_xlim(xmin=np.min(self.vel),xmax=np.max(self.vel))
                axi.set_ylabel('I')
                axi.set_xlabel('Velocity (km/s)')
            elif i == 0: #Stokes V (if nplots > 1)
                axi.errorbar(self.vel, self.specV, yerr=self.specSigV, xerr=None,
                             fmt='o', ms=3, ecolor='0.5',c='k',**kwargs)
                if plotZeroLevel: axi.axhline(y=0.0, ls='--', c='k', alpha=0.4)
                if sameYRange: axi.set_ylim(ymin=-plotVLims,ymax=plotVLims)
                axi.set_ylabel('V')
            elif i == 1: #Null1
                axi.errorbar(self.vel, self.specN1, yerr=self.specSigN1, xerr=None,
                             fmt='o', ms=3,ecolor='0.5',c='k',**kwargs)
                if plotZeroLevel: axi.axhline(y=0.0, ls='--', c='k', alpha=0.4)
                if sameYRange: axi.set_ylim(ymin=-plotVLims,ymax=plotVLims)
                axi.set_ylabel('N1')
            elif i == 2: #Null2
                axi.errorbar(self.vel, self.specN2, yerr=self.specSigN2, xerr=None,
                             fmt='o', ms=3,ecolor='0.5',c='k',**kwargs)
                if plotZeroLevel: axi.axhline(y=0.0, ls='--', c='k', alpha=0.4)
                if sameYRange: axi.set_ylim(ymin=-plotVLims,ymax=plotVLims)
                axi.set_ylabel('N2')
        
        plt.subplots_adjust(hspace=.0)
        return(fig, ax)

    def rvfit():
        return

    def cog_I(self, Ic=1.0):
        '''
        Helper function to return the center of gravity of a Stokes I LSD profile, 
        for a given continuum level (default =1.0)

        :param self: LSD object
        :param Ic: (float, default=1.0) the continnum level to use the the COG calculation
        '''
        nominator = np.trapz(self.vel * (Ic-self.specI), x=self.vel)
        denominator = np.trapz( Ic-self.specI, x=self.vel )
        return nominator/denominator

    def cog_IV(self, Ic=1.0):
        '''
        Helper function to return the center of gravity of a (Stokes I)*(Stoke V) for a LSD profile, 
        for a given continuum level (default =1.0)

        :param self: LSD object
        :param Ic: (float, default=1.0) the continnum level to use the the COG calculation
        '''
        nominator = np.trapz(self.vel * np.abs( self.specV * (Ic-self.specI) ), x=self.vel )
        denominator = np.trapz( np.abs( self.specV * (Ic-self.specI) ), x=self.vel )
        return nominator/denominator

    def cog_V(self):
        '''
        Helper function to return the center of gravity of a Stokes V LSD profile, 

        :param self: LSD object
        '''
        nominator = np.trapz(self.vel * np.abs(self.specV), x=self.vel )
        denominator = np.trapz( np.abs(self.specV), x=self.vel )
        return(nominator/denominator)        

    def cog_min(self):
        '''
        Helper function to return the velocity of the minimum of a Stokes I profile
        
        :param self: a LSD object
        '''
        
        cog_min = self.vel[self.specI.argmin()]
        if cog_min.size > 1:
            cog_min = cog_min[0]
        return cog_min 

    def calc_fap(self):
        '''Helper function to return the V, null1, and null2 FAP for a given LSD object

        The False Alarm Probability (FAP) is the probability that the observed
        data are consistent with the null hypothesis of no magnetic field.
        In this case, the null hypothesis is a flat line in Stokes V
        (an offset from 0 is allowed to account for continuum polarization).
        The probability is evaluated from chi^2.
        
        If you would like a specific range in velocity, simply slice the LSD object beforehand. 
        Note that the calcBz function also returns the FAP inside the spectral line, 
        over the same velocity range as used in the Bz calculation. 

        :param self: LSD object (input)

        :return: FAP V, FAP N1, FAP N2. 
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

    def calc_bz(self, cog='I', norm='auto', lambda0=500., geff=1.2, velrange=None, bzwidth=None, plot=True):
        '''Calculate the Bz of an LSD profile
        
        :param self: lsd object (input). It is assumed that the lsd.vel is in km/s.
        :param cog: The value, or calculation method for the center of gravity. The choices are:
                    'I': center of gravity of I,
                    'V': center of gravity of V,
                    'IV': center of gravity of I*V,
                    a float: a user defined value in km/s.
        :param norm: calculation method for the continuum. The choices are:
                        'auto': the median of I outside of velrange (if defined) or the full range (if velrange is not defined)
                        float: a user defined value to use for Ic.
        :param lambda0: wavelength of the transition in nanometer (default=500).
                        For an LSD profile, this is the lambda value the LSD profile shape was scaled with.
        :param geff: effective Lande factor of the transition.
                    For an LSD profile, this is the geff value the LSD profile shape was scaled with.
        :param velrange: range of velocity to use for the determination of the
                        line center and the continuum. If not defined, the whole range
                        will be used. If bzwidth is not defined, this range will also be
                        used for the line Bz calculation.
        :param bzwidth: distance from the line center to use in the Bz calculation.
                        One element = same on each side of line center.
                        Two elements, left and right of line center.
                        Not defined: using velrange.
        :param plot: whether or not a graph is generated, and returned.
        :return: a dictionary with Bz (in G) and FAP calculations,
                optionally a matplotlib figure.
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
        blv, blvSig = integrate_bz(lsd_bz.vel, lsd_bz.specV, lsd_bz.specSigV, 
                    geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
        if self.numParam > 2:
            bln1, bln1Sig = integrate_bz(lsd_bz.vel, lsd_bz.specN1, lsd_bz.specSigN1,
                    geff, lambda0, cog_val, lsd_bz.specI, lsd_bz.specSigI, norm_val)
        else: bln1, bln1Sig = (0.0, 0.0)
        if self.numParam > 3:
            bln2, bln2Sig = integrate_bz(lsd_bz.vel, lsd_bz.specN2, lsd_bz.specSigN2,
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
            fig  = plot_bz_calc(self, lsd_in, lsd_bz, velrange,
                            p_bzwidth, norm_val, cog_val, cog)
            return result,fig
        else:
            return result

def read_lsd(fname):
    """
    function that reads in a LSD profile.
    
    The LSD profiles are in Donati's text format.
    The two lines of header in Donati's format is optional.
    
    :param fname: the name of the file containing the LSD profile
    :rtype: returns an instance of the LSD class, defined in this module
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
###################################

def run_lsdpy(obs=None, mask=None, outName='prof.dat',
         velStart=None, velEnd=None, velPixel=None, 
         normDepth=None, normLande=None, normWave=None,
         removeContPol=None, trimMask=None, sigmaClipIter=None, sigmaClip=None, 
         interpMode=None, outModelName='',
         fLSDPlotImg=None, fSavePlotImg=None, outPlotImgName=None):
    """Run the LSDpy code and return an LSD object.
    (A convenience wrapper around the lsdpy.main() function.)
    
    Any arguments not specified will be read from the file inlsd.dat.
    The file inlsd.dat is optional, but if the file dose not exist and 
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
                          (Default = 0)
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
    :rtype: Returns an LSD profile object.  Also the model spectrum (the fit
                          to the observation, i.e. the convolution of the line
                          mask and LSD profile), inside an 'observation' object.
    """
    import LSDpy

    vel, sI, sIerr, sV, sVerr, sN1, sN1err, headerTxt, specList = LSDpy.lsdpy.main(
        obs, mask, outName, velStart, velEnd, velPixel, 
        normDepth, normLande, normWave, removeContPol, trimMask, 
        sigmaClipIter, sigmaClip, interpMode, 1, outModelName, 
        fLSDPlotImg, fSavePlotImg, outPlotImgName)

    prof = LSD(vel, sI, sIerr, sV, sVerr, sN1, sN1err, header=headerTxt)
    modelSpec = observation(specList[0], specList[1], specList[2], specList[3],
                            np.zeros_like(specList[0]), np.zeros_like(specList[0]))
    return prof, modelSpec

###################################
###################################

class Mask:
    """
    The data for the LSD line Mask.

    This usually contains arrays:
    
    * wl - wavelengths of lines
    * element - the element+ion code for the line
    * depth - the depth of the line
    * lande - the effective Lande factor of the line
    * iuse - integer flag for whether the line is used
    """
    def __init__(self, wl, element, depth, excite, lande, iuse):
        """
        Generate a Mask object
        """
        self.wl = wl
        self.element = element
        self.depth = depth
        self.excite = excite
        self.lande = lande
        self.iuse = iuse

    def __getitem__(self, key):
        """Overloaded getitem function. Returns a Mask object with only the values at the specified index(s).

        :param key: the index or slice being checked
        :rtype: Mask
        """
        wl_s = self.wl[key]
        element_s = self.element[key]
        depth_s = self.depth[key]
        excite_s = self.excite[key]
        lande_s = self.lande[key]
        iuse_s = self.iuse[key]
        return Mask(wl_s, element_s, depth_s, excite_s, lande_s, iuse_s)

    def __setitem__(self, key, newval):
        """
        Overloaded setitem function. Sets all values of the mask at the specified location equal to the input mask's values.

        :param key: the index or slice being overwritten
        :param newval: Mask whose values are to replace the overwritten ones
        """
        if not(isinstance(newval, Mask)):
            raise TypeError()
        else:
            self.wl[key] = newval.wl
            self.element[key] = newval.element
            self.depth[key] = newval.depth
            self.excite[key] = newval.excite
            self.lande[key] = newval.lande
            self.iuse[key] = newval.iuse

    def __len__(self):
        '''Overloaded len function that returns the number of lines in the mask'''
        return len(self.wl)
    
    def prune(self):
        """
        Return a Mask object with unused lines removed from the Mask.
        
        Remove lines if iuse index is set to 0,
        restricting the Mask to only lines used in LSD.
        """
        #Restrict the Mask to only lines flagged to be used
        ind2 = np.where(self.iuse != 0)
        return Mask(self.wl[ind2],
                        self.element[ind2],
                        self.depth[ind2],
                        self.excite[ind2],
                        self.lande[ind2],
                        self.iuse[ind2]
                        )
        
    def get_weights(self, normDepth, normWave, normLande):
        """
        Returns the calculated the LSD weights of all the lines in the mask
        (no matter is iuse is 0 or 1).
        
        This assumes the Stokes I lines are weighted as depth,
        and Stokes V is weighted as depth*wavelength*Lande factor
        
        :param normDepth: the normalizing line depth for the mask/LSD profile
        :param normWave: the normalizing wavelength (nm) for the mask/LSD profile
        :param normLande: the normalizing effective Lande factor for the mask/LSD profile
        :return: weightI, weightV
        """
        weightI = self.depth / normDepth
        weightV = self.depth*self.wl*self.lande / (normDepth*normWave*normLande)
        return weightI, weightV
    
    def save(self, fname):
        """
        Save the line Mask to a text file, in Donati's and LSDpy format.
        
        :param fname: the file name to output the Mask to.
        """
        
        nlines = self.wl.shape[0]
        with open(fname, 'w') as oFile:
            oFile.write('{:d}\n'.format(nlines))
            for i in range(nlines):
                oFile.write('{:9.4f} {:6.2f} {:6.3f} {:6.3f} {:6.3f} {:2d}\n'.format(
                    self.wl[i], self.element[i], self.depth[i],
                    self.excite[i], self.lande[i], self.iuse[i]))

        return
    
    def clean(self, regions):
        '''
        Use a ExcludeMaskRegions object and returns a Mask Object in which
        set the iuse column is set to zero for
        for spectral lines with wavelengths in these regions. 
        
        :param regions: An ExcludeMaskRegions object
        '''
        # Making a copy of the 
        mask_clean = copy.deepcopy(self)
        nregions = len(regions)
        for i in range(0,nregions):
            is_in = np.logical_and( (self.wl>=regions[i].start), (self.wl<=regions[i].stop) )
            mask_clean.iuse[is_in] = 0

        return mask_clean

def read_mask(fname):
    """
    Read in an LSD line mask file and returns a Mask object.

    The mask file should one line of header and columns of:
    * Wavelength (nm)
    * Atomic number + (ionization state)*0.01
    * Line depth
    * Excitation potential of the lower level (eV)
    * Effective Lande factor
    * Flag for whether the line is used (1=use, 0=skip).

    :param fname: the name of the file to read.
    :rtype: Mask
    """
    tmpMask = np.loadtxt(fname, skiprows=1, unpack=True)
    
    #Sort the line mask so wavelength is always increasing
    ind = np.argsort(tmpMask[0,:])
    
    wl = tmpMask[0, ind]
    element = tmpMask[1, ind]
    depth = tmpMask[2, ind]
    excite = tmpMask[3, ind]
    lande = tmpMask[4, ind]
    iuse = tmpMask[5, ind].astype(int)
    
    return Mask(wl, element, depth, excite, lande, iuse)

###################################
###################################

class ExcludeMaskRegions:
    '''
    Class for a region object that records spectral regions to exclude from a Mask.
    '''

    def __init__(self, start, stop, type):
        self.start = start
        self.stop = stop
        self.type = type

    def __getitem__(self, key):
        """Overloaded getitem function. Returns a region object with only the values at the specified index(s).

        :param key: the index or slice being checked
        :rtype: ExcludeMaskRegions
        """
        return ExcludeMaskRegions(self.start[key], self.stop[key], self.type[key])

    def __setitem__(self, key, newval):
        """
        Overloaded setitem function. Sets all values of the Mask at the specified location equal to the input mask's values.

        :param key: the index or slice being overwritten
        :param newval: Mask whose values are to replace the overwritten ones
        """ 
        if not(isinstance(newval, ExcludeMaskRegions)):
            raise TypeError()
        else:
            self.start[key] = newval.start
            self.stop[key] = newval.stop
            self.type[key] = newval.type

    def __len__(self):
        '''Overloaded len function that returns the number of lines in the mask'''
        return len(self.start)

    def __add__(self, other):
        '''
        Overloaded addition function to concatenate two ExcludeMaskRegions objects
        '''
        start = np.concatenate([self.start, other.start])
        stop = np.concatenate([self.stop, other.stop])
        type = np.concatenate([self.type, other.type])
        return ExcludeMaskRegions(start, stop, type)

    def save(self, fname):
        '''
        Save the ExcludeMaskRegions object to a text file.
        
        :param fname: the file path/name
        """
        '''
        with open(fname, 'w') as ofile:
            for item in self:
                ofile.write('{} {} {}\n'.format(item.start, item.stop, item.type))
        return

    def to_dict(self):
        '''
        Function to return the ExcludeMaskRegion as a dictionary. 
        This is useful to transform ExcludeMaskRegions objects to Panda dataframes
        '''
        return({'start':self.start, 'stop':self.stop,'type':self.type})

def read_exclude_mask_regions(fname):
    '''
    Read in a ExcludeMaskRegion file into an ExcludeMaskRegion object. 

    :param fname: the path/name of the file
    :rtype ExcludeMaskRegion: 
    '''
    start = []
    stop = []
    type = []
    with open(fname, 'r') as f:
        nLines = 0
        for txtLine in f:
            split = txtLine.split() # split by WS
            start.append(split[0]) 
            stop.append(split[1]) 
            type.append(' '.join(split[2:])) # join all the rest with ' ' in case there was multiple words
    
    start = np.asarray(start, dtype=float)
    stop = np.asarray(stop, dtype=float)
    
    return ExcludeMaskRegions(start, stop, np.array(type))

def get_Balmer_regions_default(velrange=500):
    '''
    Returns a ExcludeMaskRegion object with regions around Balmer H-lines (alpha to epsilon) 
    up to a given radial velocity, and a region that exclude the Balmer jump (from 360-392 nm)
    
    :param velrange: (default 500 km/s) velocity range around the H-line to be excluded.
    :rtype: ExcludeMaskRegion
    '''

    c = 299792.458 # Speed of light in km/s
    # Mask should be in nm
    wavelengths = [
        656.281,
        486.14,
        434.05,
        410.17,
        397.01
    ]
    types = [
            'Halpha',
            'Hbeta',
            'Hgamma',
            'Hdelta',
            'Hepsilon',
            'Hjump'
    ]
    start = []
    stop = []
    for w in wavelengths:
        start.append(-1*velrange/c*w + w)
        stop.append(velrange/c*w + w)

    # Adding the Balmer jump
    start.append(360)
    stop.append(392)

    return ExcludeMaskRegions(np.array(start), np.array(stop), np.array(types,dtype=object))

def get_telluric_regions_default():
    '''
    Returns a ExcludeMaskRegions object with regions with heavy telluric regions in the optical
    '''
    start = np.array([587.5,627.5,684.0,717.0,757.0,790.0,809.0])  # nm
    stop   = np.array([592.0,632.5,705.3,735.0,771.0,795.0,990.0])  # nm

    return ExcludeMaskRegions(start, stop, np.array(['telluric']*len(start),dtype=object))

###################################
###################################

class observation:
    """
    Contains an observed spectrum, usually spectropolarimetric data.

    Usually contains arrays:
    
    * wl - wavelengths
    * specI - Stokes I spectrum
    * specV - polarized spectrum, usually Stokes V
    * specN1 - the first polarimetric null spectrum
    * specN2 - the second polarimetric null spectrum
    * specSig - the formal uncertainties, which apply to the other spectra
    """

    def __init__(self,wl, specI, specV, specN1, specN2, specSig, header=None):
        
        self.header=header
        self.wl = wl
        self.specI=specI
        self.specV=specV
        self.specN1=specN1
        self.specN2=specN2
        self.specSig=specSig
        
    def __getitem__(self, key):
        """Overloaded getitem function. Returns an observation object with only the values at the specified index(s).

        :param self: observation being queried
        :param key: the index or slice being checked

        :rtype: observation
        """
        wl_s = self.wl[key]
        specI_s = self.specI[key]
        specSig_s = self.specSig[key]
        specV_s = self.specV[key]
        specN1_s = self.specN1[key]
        specN2_s = self.specN2[key]
        #The header may be None but that is ok.
        slice_obs = observation(wl_s, specI_s, specV_s, specN1_s, specN2_s, specSig_s, header=self.header)
        return slice_obs

    def __setitem__(self, key, newval):
        """
        Overloaded setitem function. Sets all values of the observation at the specified location equal to the input observation's values.

        :param self: observation object being edited
        :param key: the index or slice being overwritten
        :param newval: observation whose values are to replace the overwritten ones
        """
        if not(isinstance(newval, observation)):
            raise TypeError()
        else:
            self.wl[key] = newval.wl
            self.specI[key] = newval.specI
            self.specSig[key] = newval.specSig
            self.specV[key] = newval.specV
            self.specN1[key] = newval.specN1
            self.specN2[key] = newval.specN2

    def __len__(self):
        return len(self.wl)

    def save(self, fname, saveHeader=True):
        '''
        Write the observation into a .s LibreESPRIT style format
        Optionally skip writing the two lines of header

        :param saveHeader: optional flag to skip writing the header if False
        '''

        #Note, the LibreESPRIT .s format header counts the number of columns
        #data columns not counting the first wavelength column
        ncols = 5
        #Support 3 column (intensity only) spectra
        if np.all(self.specV == 0.):
            ncols = 2
            #and 2 column (intensity only, no errorbars) spectra
            if np.all(self.specSig == 0):
                ncols = 1
        
        with open(fname, 'w') as f:
            #Optionaly write 2 lines of header            
            if saveHeader:
                if self.header is None:
                    f.write('*** Spectrum of\n')
                else:
                    f.write(self.header)
                    #Make sure there is a line break after the first header text
                    if self.header[-1] != '\n': f.write('\n')
                f.write('{:7n} {:1n}\n'.format(int(self.wl.size), ncols))
            
            if ncols == 5:
                for i in range(self.wl.size):
                    f.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(
                        self.wl[i], self.specI[i], self.specV[i],
                        self.specN1[i], self.specN2[i], self.specSig[i]))
            elif ncols == 2:
                for i in range(self.wl.size):
                    f.write('{:10.4f} {:11.4e} {:11.4e}\n'.format(
                        self.wl[i], self.specI[i], self.specSig[i]))
            elif ncols == 1:
                for i in range(self.wl.size):
                    f.write('{:10.4f} {:11.4e}\n'.format(
                        self.wl[i], self.specI[i]))
        return

def read_spectrum(fname, trimBadPix=False, sortByWavelength=False):
    """
    Read in the observed spectrum and save it.
    
    This follows the .s format from Donati's LibreESPRIT,
    files can either have two lines of header or no header.
    This supports 6 column spectropolarimetric files
    (wavelength, I, V|Q|U, null1, null2, errors),
    and also 3 column spectra (wavelength, I, errors).

    :param fname: the name of the file to read.
    :param sortByWavelength: reorder the points in the spectrum to always
                             increase in wavelength, if set to True.
    :param trimBadPix: optionally remove the more obviously bad pixels if True.
                       Removes pixels with negative flux or error bars of zero,
                       pixels within 3sigma of zero (large error bars),
                       and pixels with extremely large values.
    :return: a Spectrum object containing the observation.
    """
    # Reading manually is often faster than np.loadtxt for a large files
    fObs = open(fname, 'r')
    #Check if the file starts with data or a header (assume 2 lines of header)
    line1 = fObs.readline()
    line2 = fObs.readline()
    line3 = fObs.readline()
    #assume at least the 3rd line contains real data
    ncolumns = len(line3.split())
    if ncolumns != 6 and ncolumns != 3 and ncolumns != 2:
        print('{:} column spectrum: unknown format!\n'.format(ncolumns))
        raise ValueError(('Reading {:} as an {:} column spectrum: '
                          'unknown format!').format(fname, ncolumns))
    
    if len(line1.split()) == ncolumns and len(line2.split()) == ncolumns:
        #If the column counts are consistent there may be no header
        try:
            #and if the first line starts with numbers, probably no header
            float(line1.split()[0])
            float(line1.split()[1])
            float(line2.split()[0])
            float(line2.split()[1])
            obs_header = None
            nHeader = 0
        except ValueError:
            #Otherwise assume there are two lines of header,
            #one line of with comment and a second with file dimensions,
            #but we figure that by reading the file.
            obs_header = line1
            nHeader = 2
    else:
        obs_header = line1
        nHeader = 2

    #Get the number of lines of data in the file
    nLines = 3 - nHeader #we've alread read 3 lines
    for line in fObs:
        words = line.split()
        if len(words) == ncolumns:
            nLines += 1
        else:
            print('ERROR: reading observation, '
                  +'line {:}, {:} columns :\n{:}'.format(
                      nLines, len(words), line))

    obs = observation(np.zeros(nLines), np.zeros(nLines), np.zeros(nLines),
                      np.zeros(nLines), np.zeros(nLines), np.zeros(nLines),
                      header=obs_header)
    
    #Rewind to start then advance the file pointer 2 lines
    fObs.seek(0)
    if obs_header is not None:
        fObs.readline()
        fObs.readline()
    #Then read the actual data of the file
    for i, line in enumerate(fObs):
        words = line.split()
        if (len(words) == ncolumns and ncolumns == 6):
            obs.wl[i] = float(words[0])
            obs.specI[i] = float(words[1])
            obs.specV[i] = float(words[2])
            obs.specN1[i] = float(words[3])
            obs.specN2[i] = float(words[4])
            obs.specSig[i] = float(words[5])
        elif (len(words) == ncolumns and ncolumns == 3):
            obs.wl[i] = float(words[0])
            obs.specI[i] = float(words[1])
            obs.specSig[i] = float(words[2])
        elif (len(words) == ncolumns and ncolumns == 2):
            obs.wl[i] = float(words[0])
            obs.specI[i] = float(words[1])
            
    fObs.close()

    #Optionally, remove bad pixels.
    if trimBadPix:
        use_flag = (obs.specI > 0.)
        if not np.all(obs.specSig == 0):  #(all zeros => no error column)
            use_flag = use_flag & (obs.specSig > 0.)
        use_flag = use_flag & (obs.specI > 3*obs.specSig)
        use_flag = use_flag & (obs.specI < 10*np.percentile(obs.specI, 99.9))
        obs = obs[use_flag]
    
    #Optionally, sort the observation so wavelength is always increasing
    if sortByWavelength:
        obs_ind = np.argsort(obs.wl)
        obs = obs[obs_ind]
    
    return obs

###################################
###################################

class line_list:
    """
    Container for a set of spectral line data, usually from VALD.

    This usually contains:
    
    * nLines - number of lines in the line list
    * ion - list of species identifiers (element or molecule and ionization)
    * wl - array of wavelengths
    * loggf - array of oscillator strengths (log gf)
    * Elo - array of excitation potentials for the lower level in the transition (in eV)
    * Jlo - array of J quantum numbers for the lower level
    * Eu - array of excitation potentials for the upper level in the transition (in eV)
    * Jup - array of J quantum numbers for the upper level
    * landeLo - array of Lande factors for the lower level
    * landeUp - array of Lande factors for the upper level
    * landeEff - array of effective Lande factors for the transition
    * rad - array of radiative damping coefficients
    * stark - array of quadratic Stark damping coefficients
    * waals - array of van der Waals damping coefficients
    * depth - depth at the centre of the spectral line, as estimated by VALD
    * configLo - list of strings with the electron configuration and term symbols for the lower level
    * configUp - list of strings with the electron configuration and term symbols for the upper level
    * refs - list of references for the sources of the line data (optional)
    """
    def __init__(self, ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo, landeUp,
                 landeEff, rad, stark, waals, depth, configLo, configUp, refs):
        self.ion      = ion
        self.wl       = wl
        self.loggf    = loggf
        self.Elo      = Elo
        self.Jlo      = Jlo
        self.Eup      = Eup
        self.Jup      = Jup
        self.landeLo  = landeLo
        self.landeUp  = landeUp
        self.landeEff = landeEff
        self.rad      = rad
        self.stark    = stark
        self.waals    = waals
        self.depth    = depth
        self.configLo = configLo
        self.configUp = configUp
        self.refs     = refs
        self.nLines   = self.wl.size

    def __getitem__(self, key):
        """Overloaded getitem function. Returns a line_list with only the values at the specified index(s).

        :param key: the index or slice being checked
        :rtype: line_list
        """
        ion      = self.ion[key]
        wl       = self.wl[key]
        loggf    = self.loggf[key]
        Elo      = self.Elo[key]
        Jlo      = self.Jlo[key]
        Eup      = self.Eup[key]
        Jup      = self.Jup[key]
        landeLo  = self.landeLo[key]
        landeUp  = self.landeUp[key]
        landeEff = self.landeEff[key]
        rad      = self.rad[key]
        stark    = self.stark[key]
        waals    = self.waals[key]
        depth    = self.depth[key]
        configLo = self.configLo[key]
        configUp = self.configUp[key]
        refs     = self.refs[key]
        lList =  line_list(ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo,
                           landeUp, landeEff, rad, stark, waals, depth,
                           configLo, configUp, refs)
        return lList

    def __setitem__(self, key, newval):
        """
        Overloaded setitem function. Sets all values of the line_list at the specified location equal to the input line_list values.

        :param key: the index or slice being overwritten
        :param newval: line_list used to replace the values given by key
        """
        if not(isinstance(newval, line_list)):
            raise TypeError()
        else:
            self.ion[key]      = newval.ion
            self.wl[key]       = newval.wl
            self.loggf[key]    = newval.loggf
            self.Elo[key]      = newval.Elo
            self.Jlo[key]      = newval.Jlo
            self.Eup[key]      = newval.Eup
            self.Jup[key]      = newval.Jup
            self.landeLo[key]  = newval.landeLo
            self.landeUp[key]  = newval.landeUp
            self.landeEff[key] = newval.landeEff
            self.rad[key]      = newval.rad
            self.stark[key]    = newval.stark
            self.waals[key]    = newval.waals
            self.depth[key]    = newval.depth
            self.configLo[key] = newval.configLo
            self.configUp[key] = newval.configUp
            self.refs[key]     = newval.refs
            self.nLines   = self.wl.size

    def __len__(self):
        return self.nLines

    def __str__(self):
        """
        Generate a nicely formatted string of line data for printing
        """
        strList = []
        outStr = '\n'
        for line in self:
            fmt = ("'{:s}',{:17.4f},{:7.3f},{:8.4f},{:5.1f},{:8.4f},{:5.1f},"
                   "{:7.3f},{:7.3f},{:7.3f},{:6.3f},{:6.3f},{:8.3f},{:6.3f},")
            if len(line.ion) == 3: fmt = fmt[:7] + ' ' + fmt[7:]
            strList.append(fmt.format(line.ion, line.wl, line.loggf,
                                line.Elo, line.Jlo, line.Eup, line.Jup,
                                line.landeLo, line.landeUp, line.landeEff,
                                line.rad, line.stark, line.waals, line.depth))
            strList.append("'  {:>88}'".format(line.configLo))
            strList.append("'  {:>88}'".format(line.configUp))
            strList.append("'{:}'".format(line.refs))
        return outStr.join(strList)

    def write_VALD(self, fname):
        """
        Write a line list to a text file.
        This outputs using the VALD version 3 'extract stellar' 'long' format.
        
        A few details (e.g. references) are omitted since they are not saved
        in the line_list class.
        
        :param fname: the file name to save the output to
        """
        
        fOut = open(fname, 'w')
        fOut.write(("{:11.5f},{:12.5f},{:5d},{:7d},{:4.1f}, Wavelength "
                    "region, lines selected, lines processed, Vmicro\n"
                    ).format(self.wl[0], self.wl[-1], self.nLines, 999999, 0.))
        fOut.write("                                                 "
                   "                    Lande factors       "
                   "Damping parameters   Central\n")
        fOut.write("Spec Ion       WL_air(A)  log gf* E_low(eV) J lo "
                   "E_up(eV)  J up  lower   upper    mean   Rad.  "
                   "Stark   Waals   depth\n")
        fOut.write(str(self)+'\n')
        fOut.close()
        return
    
def line_list_zeros(nLines):
    """
    Generate a line list of zeros and blank text.
    
    Used by read_VALD (It can be a bit faster to allocate all the array space at once)
    
    :param nLines: the number of lines in the line_list of zeros
    :rtype: line_list
    """
    ion      = np.tile(np.array([''], dtype='U6'), nLines)
    wl       = np.zeros(nLines)
    loggf    = np.zeros(nLines)
    Elo      = np.zeros(nLines)
    Jlo      = np.zeros(nLines)
    Eup      = np.zeros(nLines)
    Jup      = np.zeros(nLines)
    landeLo  = np.zeros(nLines)
    landeUp  = np.zeros(nLines)
    landeEff = np.zeros(nLines)
    rad      = np.zeros(nLines)
    stark    = np.zeros(nLines)
    waals    = np.zeros(nLines)
    depth    = np.zeros(nLines)
    configLo = np.tile(np.array([''], dtype='U128'), nLines)
    configUp = np.tile(np.array([''], dtype='U128'), nLines)
    refs = np.tile(np.array(['_          unknown source'],dtype='U180'), nLines)
    lList = line_list(ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo,
                      landeUp, landeEff, rad, stark, waals, depth,
                      configLo, configUp, refs)
    return lList

def read_VALD(fname):
    """
    Read a list of spectral line data from VALD and return a line_list

    This expects VALD version 3 line list, in an 'extract stellar' 'long' format.

    :param fname: the file name for the VALD line list.
    :rtype: line_list object, containing arrays of line data
    """
    fVald = open(fname, 'r')
    i = 0
    j = 0
    nLines = 0
    for txtLine in fVald:
        if i == 0:
            nLines = int(txtLine.split(',')[2])
            llist = line_list_zeros(nLines)
        #There should be 3 lines of header,
        #then 4 file lines for each set of spectra line data.
        if i > 2 and i < nLines*4+3:
            ii = (i-3)%4
            if ii == 0:
                #The line data are: Spec Ion, Wl(A), log gf, E_lower(eV), 
                #J_lower, E_upper(eV), J_upper, Lande factor lower, 
                #Lande factor upper, effective Lande factor, Rad. damping, 
                #Stark damping, van der Waals damping, central depth
                vals = txtLine.split(',')
                llist.ion[j]      = vals[0].strip('\'')
                llist.wl[j]       = float(vals[1])
                llist.loggf[j]    = float(vals[2])
                llist.Elo[j]      = float(vals[3])
                llist.Jlo[j]      = float(vals[4])
                llist.Eup[j]      = float(vals[5])
                llist.Jup[j]      = float(vals[6])
                llist.landeLo[j]  = float(vals[7])
                llist.landeUp[j]  = float(vals[8])
                llist.landeEff[j] = float(vals[9])
                llist.rad[j]      = float(vals[10])
                llist.stark[j]    = float(vals[11])
                llist.waals[j]    = float(vals[12])
                llist.depth[j]    = float(vals[13])
            if ii == 1:
                llist.configLo[j] = txtLine.strip(' \n\'')
            if ii == 2:
                llist.configUp[j] = txtLine.strip(' \n\'')
            if ii == 3:
                llist.refs[j] = txtLine.strip(' \n\'')
                j += 1
        i += 1
    
    fVald.close()
    return llist
