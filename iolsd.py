## @module iolsd.py
# Documentation for iolsd.py
#
# Tools for reading and writing files, related to calculating
# and analyzing LSD profiles.

import numpy as np
import matplotlib.pyplot as plt

class lsd_prof:
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

        :param self: lsd_prof being created
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
        :param header: optionally, one line of header text for the output file.  This text string should end with a newline.
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
            oFile.write(' {:d} 6\n'.format(self.npix))
        
        for i in range(self.npix):
            oFile.write('{:>12.6f} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e}\n'.format(
                self.vel[i], self.specI[i], self.specSigI[i], self.specV[i],
                self.specSigV[i], self.specN1[i], self.specSigN1[i]))
        oFile.close()
        return

    #overloaded functions
    def __len__(self):
        """
        Return the length of each array in an LSD profile. They should all be the same length, so it just returns the length of the velocity array. 

        :param self: lsd_prof whose length is being checked

        :rtype: int
        """
        return len(self.vel)

    def __getitem__(self, key):
        """Overloaded getitem function. Returns an lsd_prof with only the values at the specified index(s).

        :param self: lsd_prof being queried
        :param key: the index or slice being checked

        :rtype: lsd_prof
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
        slice_prof = lsd_prof(vel_s, specI_s, specSigI_s, specV_s, specSigV_s,
                              specN1_s, specSigN1_s, specN2_s, specSigN2_s, self.header)
        slice_prof.numParam = self.numParam
        return slice_prof

    def __setitem__(self, key, newval):
        """
        Overloaded setitem function. Sets all values of the lsd_prof at the specified location equal to the input profile's values.

        :param self: lsd_prof object being edited
        :param key: the index or slice being overwritten
        :param newval: lsd_prof whose values are to replace the overwritten ones
        """
        if not(isinstance(newval, lsd_prof)):
            raise TypeError()
        else:
            self.vel[key] = newval.vel[:]
            self.specI[key] = newval.specI[:]
            self.specSigI[key] = newval.specSigI[:]
            self.specV[key] = newval.specV[:]
            self.specSigV[key] = newval.specSigV[:]
            self.specN1[key] = newval.specN1[:]
            self.specSigN1[key] = newval.specSigN1[:]
            self.specN2[key] = newval.specN2[:]
            self.specSigN2[key] = newval.specSigN2[:]

    def __mul__(self, other):
        """
        Overloaded multiplication function. Allows you to do lsd * n and multiply all profile values in the lsd, other than the velocity, by n.

        :param self: lsd_prof being scaled
        :param other: number to multiply by

        :rtype: lsd_prof
        """
        print('* for LSD profiles is to be depreciated, use norm() instead')
        self.specI = np.multiply(self.specI, other)
        self.specSigI = np.multiply(self.specSigI, other)
        self.specV = np.multiply(self.specV, other)
        self.specSigV = np.multiply(self.specSigV, other)
        self.specN1 = np.multiply(self.specN1, other)
        self.specSigN1 = np.multiply(self.specSigN1, other)
        self.specN2 = np.multiply(self.specN2, other)
        self.specSigN2 = np.multiply(self.specSigN2, other)
        return self

    def __rmul__(self, other):
        #Overloaded reverse multiplication function, for n * lsd.
        return self*other
    
    
    def __add__(self, other):
        """
        Overloaded addition function. Allows you to do lsd + n and add n to all values in the velocity array. 
        :param self: lsd_prof being added to
        :param other: 
        
        :rtype: lsd_prof
        """
        print('+ for LSD profiles is to be depreciated, use shift() instead')
        self.vel = self.vel + other
        return self
    
    def __radd__(self, other):
        #Overloaded reverse addition function, for n + lsd.
        return self + other
    
    def __sub__(self, other):
        """
        Overloaded subtraction function. Allows you to do lsd - n and subtract n from all values in the velocity array. 
        :param self: lsd_prof being subtracted from
        :param other: 
        
        :rtype: lsd_prof
        """
        print('- for LSD profiles is to be depreciated, use shift() instead')
        self.vel = self.vel - other
        return self
    
    def __rsub__(self, other):
        #Overloaded reverse subtraction function, for n + lsd.
        print('- for LSD profiles is to be depreciated, use shift() instead')
        self.vel = other - self.vel
        return self
    
    #def __floordiv__(self, other):  #Currently just performs regular division not 'floor division' //
    #    self.specI = np.divide(self.specI, other)
    #    self.specSigI = np.divide(self.specSigI, other)
    #    self.specV = np.divide(self.specV, other)
    #    self.specSigV = np.divide(self.specSigV, other)
    #    self.specN1 = np.divide(self.specN1, other)
    #    self.specSigN1 = np.divide(self.specSigN1, other)
    #    self.specN2 = np.divide(self.specN2, other)
    #    self.specSigN2 = np.divide(self.specSigN2, other)
    #    return self
    
    def __truediv__(self, other):
        print('/ for LSD profiles is to be depreciated, use norm() instead')
        self.norm(other)
        return self

    def norm(self, normValue):
        """
        Return a renormalize an LSD profile, divide the I, V, and null profiles by a value.
        
        :param normValue: the value to renormalize (divide) the LSD profile by
        :rtype: lsd_prof
        """
        # FROM VERO:
        # It could look like this:
        #new = copy.deepcopy(self) 
        #new.specI = np.divide(new.specI, normValue)
        #new.specSigI = np.divide(new.specSigI, normValue)
        #new.specV = np.divide(new.specV, normValue)
        #new.specSigV = np.divide(new.specSigV, normValue)
        #new.specN1 = np.divide(new.specN1, normValue)
        #new.specSigN1 = np.divide(new.specSigN1, normValue)
        #new.specN2 = np.divide(new.specN2, normValue)
        #new.specSigN2 = np.divide(new.specSigN2, normValue)
        # But maybe using the constructor directly saves a package. 
        new = lsd_prof(self.vel, 
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
        :rtype: lsd_prof
        """
        # VERO: old def that changes the original object. 
        # Can be removed in cleanup later.
        #self.vel = self.vel + velShift
        #return self
        new = lsd_prof(self.vel-velShift, 
                        self.specI, self.specSigI, 
                        self.specV, self.specSigV,
                        self.specN1, self.specSigN1, 
                        self.specN2, self.specSigN2, 
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
        # VERO: old def that changes the original object. 
        # Can be removed in cleanup later.
        #self.scale(wint_new/wint_old, wpol_new/wpol_old)
        return(self.scale(wint_new/wint_old, wpol_new/wpol_old))

    def scale(self, scale_int, scale_pol):
        '''Return a LSD profile with rescaled amplitudes of the LSD profile (see also set_weights())
        
        :param scale_int: scale the intensity profile by this
        :param scale_pol: scale the polarization and null profiles by this
        :rtype: lsd object
        '''
        
        # VERO: old def that changes the original object. 
        # Can be removed in cleanup later.
        #self.specI = 1.0 - ((1.0-self.specI) * scale_int)
        #self.specSigI = self.specSigI * scale_int
        #self.specV = self.specV * scale_pol
        #self.specSigV = self.specSigV * scale_pol
        #self.specN1 = self.specN1 * scale_pol
        #self.specSigN1 = self.specSigN1 * scale_pol
        #self.specN2 = self.specN2 * scale_pol
        #self.specSigN2 = self.specSigN2 * scale_pol

        new = lsd_prof(self.vel, 
                        1.0 - ((1.0-self.specI) * scale_int), self.specSigI * scale_int, 
                        self.specV*scale_pol, self.specSigV*scale_pol,
                        self.specN1*scale_pol, self.specSigN1*scale_pol, 
                        self.specN2*scale_pol, self.specSigN2*scale_pol, 
                        self.header)
        new.numParam = self.numParam

        return new
    
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

def read_lsd(fname):
    """
    function that reads in a LSD profile.
    
    The LSD profiles are in Donati's text format.
    The two lines of header in Donati's format is optional.
    
    :param fname: the name of the file containing the LSD profile
    :rtype: returns an instance of the lsd_prof class, defined in this module
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
            prof = lsd_prof(vel, specI, specSigI, specV, specSigV,
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
            prof = lsd_prof(vel, specI, specSigI, specV, specSigV,
                        specN1, specSigN1, header=header)
    
    if ncols == 3:
        vel = __prof[0,:]
        specI = __prof[1,:]
        specSigI = __prof[2,:]
        prof = lsd_prof(vel, specI, specSigI, header=header)

    #For unsupported formats or numbers of columns
    if ncols != 3 and ncols != 7 and ncols != 9:
        raise ValueError("Read an unexpected number of columns from "
                         +"{:}, can't read as an LSD profile.".format(fname))
    return prof


def run_lsdpy(obs=None, mask=None, outName='prof.dat',
         velStart=None, velEnd=None, velPixel=None, 
         normDepth=None, normLande=None, normWave=None,
         removeContPol=None, trimMask=None, sigmaClipIter=None, sigmaClip=None, 
         interpMode=None, outModelName='',
         fLSDPlotImg=None, fSavePlotImg=None, outPlotImgName=None):
    """Run the LSDpy code and return an lsd_prof object.
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

    prof = lsd_prof(vel, sI, sIerr, sV, sVerr, sN1, sN1err, header=headerTxt)
    modelSpec = observation(specList[0], specList[1], specList[2], specList[3],
                            np.zeros_like(specList[0]), np.zeros_like(specList[0]))
    return prof, modelSpec


class mask:
    """
    The data for the LSD line mask.

    This usually contains arrays:
    
    * wl - wavelengths of lines
    * element - the element+ion code for the line
    * depth - the depth of the line
    * lande - the effective Lande factor of the line
    * iuse - integer flag for whether the line is used
    """
    def __init__(self, wl, element, depth, excite, lande, iuse):
        """
        Generate a mask object
        """
        self.wl = wl
        self.element = element
        self.depth = depth
        self.excite = excite
        self.lande = lande
        self.iuse = iuse

    def __getitem__(self, key):
        """Overloaded getitem function. Returns a mask object with only the values at the specified index(s).

        :param key: the index or slice being checked
        :rtype: mask
        """
        wl_s = self.wl[key]
        element_s = self.element[key]
        depth_s = self.depth[key]
        excite_s = self.excite[key]
        lande_s = self.lande[key]
        iuse_s = self.iuse[key]
        return mask(wl_s, element_s, depth_s, excite_s, lande_s, iuse_s)

    def __setitem__(self, key, newval):
        """
        Overloaded setitem function. Sets all values of the mask at the specified location equal to the input mask's values.

        :param key: the index or slice being overwritten
        :param newval: mask whose values are to replace the overwritten ones
        """
        if not(isinstance(newval, mask)):
            raise TypeError()
        else:
            self.wl[key] = newval.wl[:]
            self.element[key] = newval.element[:]
            self.depth[key] = newval.depth[:]
            self.excite[key] = newval.excite[:]
            self.lande[key] = newval.lande[:]
            self.iuse[key] = newval.iuse[:]

    def __len__(self):
        return self.wl.size
    
    def prune(self):
        """
        Remove unused lines from the mask.
        
        Remove lines if iuse index is set to 0,
        restricting the mask to only lines used in LSD.
        This deletes the unused lines, but may be convenient
        for more efficient processing of the mask later.
        """
        #Restrict the mask to only lines flagged to be used
        ind2 = np.where(self.iuse != 0)
        self.wl = self.wl[ind2]
        self.element = self.element[ind2]
        self.depth = self.depth[ind2]
        self.excite = self.excite[ind2]
        self.lande = self.lande[ind2]
        self.iuse = self.iuse[ind2]
        return
        
    def set_weights(self, normDepth, normWave, normLande):
        """
        Calculate the weights of the lines used for LSD calculations.
        
        This assumes the Stokes I lines are weighted as depth,
        and Stokes V is weighted as depth*wavelength*Lande factor
        
        :param normDepth: the normalizing line depth for the mask/LSD profile
        :param normWave: the normalizing wavelength (nm) for the mask/LSD profile
        :param normLande: the normalizing effective Lande factor for the mask/LSD profile
        """
        self.weightI = self.depth / normDepth
        self.weightV = self.depth*self.wl*self.lande / (normDepth*normWave*normLande)
        return
    
    def save(self, fname):
        """
        Save the line mask to a text file, in Donati's and LSDpy format.
        
        :param fname: the file name to output the mask to.
        """
        
        nlines = self.wl.shape[0]
        oFile = open(fname, 'w')
        oFile.write('{:d}\n'.format(nlines))
        
        for i in range(nlines):
            oFile.write('{:9.4f} {:6.2f} {:6.3f} {:6.3f} {:6.3f} {:2d}\n'.format(
                self.wl[i], self.element[i], self.depth[i],
                self.excite[i], self.lande[i], self.iuse[i]))
        return


def read_mask(fname):
    """
    Read in an LSD line mask and return a mask object.

    The mask file should one line of header and columns of:
    * Wavelength (nm)
    * Atomic number + (ionization state)*0.01
    * Line depth
    * Excitation potential of the lower level (eV)
    * Effective Lande factor
    * Flag for whether the line is used (1=use, 0=skip).

    :param fname: the name of the file to read.
    :rtype: mask
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
    
    return mask(wl, element, depth, excite, lande, iuse)


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
            self.wl[key] = newval.wl[:]
            self.specI[key] = newval.specI[:]
            self.specSig[key] = newval.specSig[:]
            self.specV[key] = newval.specV[:]
            self.specN1[key] = newval.specN1[:]
            self.specN2[key] = newval.specN2[:]

    def __len__(self):
        return self.wl.size

    def write_s(self, fname, noHeader=False):
        '''
        Write the observation into a .s LibreESPRIT style format
        Optionally skip writing the two lines of header

        :param noHeader: flag to skip writing the header if True
        '''

        #Note, the LibreESPRIT .s format header counts the number of columns
        #not counting the first wavelength column
        ncols = 5
        #Support 3 column (intensity only) spectra
        if np.all(self.specV == 0.): ncols = 2
        
        with open(fname, 'w') as f:
            #Optionaly write 2 lines of header            
            if noHeader == False:
                if self.header == None:
                    f.write('*** Spectrum of\n')
                else:
                    f.write(self.header)
                f.write('{:7i} {:1i}\n'.format(int(self.wl.size)), ncols)
            
            if ncols == 5:
                for i in range(0,self.wl.size):
                    f.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(
                        self.wl[i], self.specI[i], self.specV[i],
                        self.specN1[i], self.specN2[i], self.specSig[i]))
            elif ncols == 2:
                for i in range(0,self.wl.size):
                    f.write('{:10.4f} {:11.4e} {:11.4e}\n'.format(
                        self.wl[i], self.specI[i], self.specSig[i]))
        return


def read_spectrum(fname, sortByWavelength=False):
    """
    Read in the observed spectrum and save it.
    
    This follows the .s format from Donati's LibreESPRIT,
    files can either have two lines of header or no header.
    This supports 6 column spectropolarimetric files
    (wavelength, I, V|Q|U, null1, null2, errors),
    and also 3 column spectra (wavelength, I, errors).

    :param fname: the name of the file to read.
    :param sortByWavelength: reorder the points in the spectrum to always increase in wavelength, if set to True.
    """
    ## Reading manually is ~4 times faster than np.loadtxt for a large files
    fObs = open(fname, 'r')
    #Check if the file starts with data or a header (assume it is two lines)
    line = fObs.readline()
    words = line.split()
    try:
        float(words[0])
        float(words[1])
        float(words[2])
        line2 = fObs.readline()
        words2 = line2.split()
        if(len(line2) > 2 and len(words) == len(words2)):
            #If the first line behaves like spectrum data,
            #and the second line is similar to the first line
            obs_header = None
            fObs.seek(0)
        else:
            #Otherwise assume there are two lines of header,
            #the first one being a comment the second should be
            #dimensions of the file but we figure that by reading the file.
            obs_header = line
    except ValueError:
        obs_header = line
        fObs.readline()

    #Get the number of lines of data in the file
    nLines = 0
    for line in fObs:
        words = line.split()
        if(nLines == 0):
            ncolumns = len(words)
            if (ncolumns != 6):
                if(ncolumns == 3):
                    print('Apparent Stokes I only spectrum')
                    print('Generating place holder V and N columns')
                else:
                    print('{:} column spectrum: unknown format!'.format(ncolumns))
                    import sys
                    sys.exit()
        if len(words) == ncolumns:
            if ncolumns == 6:
                if(float(words[1]) > 0. and float(words[5]) > 0.):
                    nLines += 1
            elif ncolumns == 3:
                if(float(words[1]) > 0. and float(words[2]) > 0.):
                    nLines += 1
        else:
            print('ERROR: reading observation, line {:}, {:} columns :\n{:}'.format(nLines, len(words), line))

    obs_wl = np.zeros(nLines)
    obs_specI = np.zeros(nLines)
    obs_specV = np.zeros(nLines)
    obs_specN1 = np.zeros(nLines)
    obs_specN2 = np.zeros(nLines)
    obs_specSig = np.zeros(nLines)
    
    i = 0
    #Rewind to start then advance the file pointer 2 lines
    fObs.seek(0)
    if obs_header != None:
        fObs.readline()
        fObs.readline()
    #Then read the actual data of the file
    for line in fObs:
        words = line.split()
        if (len(words) == ncolumns and ncolumns == 6):
            if(float(words[1]) > 0. and float(words[5]) > 0.):
                obs_wl[i] = float(words[0])
                obs_specI[i] = float(words[1])
                obs_specV[i] = float(words[2])
                obs_specN1[i] = float(words[3])
                obs_specN2[i] = float(words[4])
                obs_specSig[i] = float(words[5])
                i += 1
        elif (len(words) == ncolumns and ncolumns == 3):
            if(float(words[1]) > 0. and float(words[2]) > 0.):
                obs_wl[i] = float(words[0])
                obs_specI[i] = float(words[1])
                obs_specSig[i] = float(words[2])
                obs_specV[i] = 0.
                obs_specN1[i] = 0.
                obs_specN2[i] = 0.
                i += 1
            
    fObs.close()
    
    #Optionally, sort the observation so wavelength is always increasing
    if sortByWavelength:
        obs_ind = np.argsort(obs_wl)
        
        obs_wl = obs_wl[obs_ind]
        obs_specI = obs_specI[obs_ind]
        obs_specV = obs_specV[obs_ind]
        obs_specN1 = obs_specN1[obs_ind]
        obs_specN2 = obs_specN2[obs_ind]
        obs_specSig = obs_specSig[obs_ind]
    
    return(observation(obs_wl, obs_specI, obs_specV, obs_specN1, obs_specN2,
                       obs_specSig, header=obs_header))


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
        self.configLo = configLo #['' for i in range(nLines)]
        self.configUp = configUp #['' for i in range(nLines)]
        self.refs     = refs #["'_          unknown source'" for i in range(nLines)]
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
            self.ion[key]      = newval.ion[:]
            self.wl[key]       = newval.wl[:]
            self.loggf[key]    = newval.loggf[:]
            self.Elo[key]      = newval.Elo[:]
            self.Jlo[key]      = newval.Jlo[:]
            self.Eup[key]      = newval.Eup[:]
            self.Jup[key]      = newval.Jup[:]
            self.landeLo[key]  = newval.landeLo[:]
            self.landeUp[key]  = newval.landeUp[:]
            self.landeEff[key] = newval.landeEff[:]
            self.rad[key]      = newval.rad[:]
            self.stark[key]    = newval.stark[:]
            self.waals[key]    = newval.waals[:]
            self.depth[key]    = newval.depth[:]
            self.configLo[key] = newval.configLo[:]
            self.configUp[key] = newval.configUp[:]
            self.refs[key]     = newval.refs[:]
            self.nLines   = self.wl.size

    def __len__(self):
        return self.nLines

    def write_VALD(self, fname):
        """
        Write a line list to a text file.
        This outputs using the VALD version 3 'extract stellar' 'long' format.
        
        A few details (e.g. references) are omitted since they are not saved
        in the line_list class.
        
        :param fname: the file name to save the output to
        """
        
        fOut = open(fname, 'w')
        fOut.write("{:11.5f},{:11.5f},{:5d},{:7d},{:4.1f}, Wavelength region, lines selected, lines processed, Vmicro\n".format(
            self.wl[0], self.wl[-1], self.nLines, 999999, 0.0))
        fOut.write("                                                                     Lande factors       Damping parameters  Central\n")
        fOut.write("Spec Ion       WL_air(A)  log gf* E_low(eV) J lo E_up(eV)  J up  lower   upper    mean   Rad.   Stark  Waals  depth\n")
        
        for i, line in enumerate(self):
            fmt = "'{:s}',{:16.4f},{:8.3f},{:8.4f},{:5.1f},{:8.4f},{:5.1f},{:7.3f},{:7.3f},{:7.3f},{:6.3f},{:6.3f},{:6.3f},{:6.3f},\n"
            if len(line.ion) == 3: fmt = fmt[:7] + ' ' + fmt[7:]
            fOut.write(fmt.format(
                line.ion, line.wl, line.loggf, line.Elo, line.Jlo,
                line.Eup, line.Jup, line.landeLo, line.landeUp, line.landeEff,
                line.rad, line.stark, line.waals, line.depth))
            fOut.write("'  {:}'\n".format(line.configLo))
            fOut.write("'  {:}'\n".format(line.configUp))
            fOut.write("'{:}'\n".format(line.refs))
        fOut.close()
        return
    
def line_list_zeros(nLines):
    """
    Generate a line list of zeros and blank text.
    
    Used by read_VALD (It can be a bit faster to allocate all the array space at once)
    
    :param nLines: the number of lines in the line_list of zeros
    :rtype: line_list
    """
    ion      = ['' for i in range(nLines)]
    wl      = np.zeros(nLines)
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
    configLo = ['' for i in range(nLines)]
    configUp = ['' for i in range(nLines)]
    refs = ["'_          unknown source'" for i in range(nLines)]
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
                llist.wl[j]      = float(vals[1])
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
