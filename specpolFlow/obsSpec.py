"""
Tools for manipulating spectra, typically spectropolarimetric observations.
"""

import copy
import warnings
import numpy as np

# use a version dependent name for trapezoidal integration,
# since Numpy 2.0 changed names (this could also be a try except)
if np.lib.NumpyVersion(np.__version__) >= '2.0.0':
    from numpy import trapezoid as _trapezoid
else:
    from numpy import trapz as _trapezoid

###################################

class Spectrum:
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

    def __init__(self, wl, specI, specV, specN1, specN2, specSig, header=None):
        
        self.header=header
        self.wl = wl
        self.specI=specI
        self.specV=specV
        self.specN1=specN1
        self.specN2=specN2
        self.specSig=specSig
        
    def __getitem__(self, key):
        """
        Returns a Spectrum object with only the values at the specified index(s)

        :param key: the index or slice being checked
        :rtype: Spectrum
        """
        wl_s = self.wl[key]
        specI_s = self.specI[key]
        specSig_s = self.specSig[key]
        specV_s = self.specV[key]
        specN1_s = self.specN1[key]
        specN2_s = self.specN2[key]
        #The header may be None but that is ok.
        slice_spec = Spectrum(wl_s, specI_s, specV_s, specN1_s, specN2_s,
                              specSig_s, header=self.header)
        return slice_spec

    def __setitem__(self, key, newval):
        """
        Sets all values of the Spectrum at the specified location equal
        to the input Spectrum's values.

        :param key: the index or slice being overwritten
        :param newval: Spectrum whose values are to replace the overwritten ones
        """
        if not(isinstance(newval, Spectrum)):
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

    def concatenate(self, *args):
        """
        Combine this Spectrum with other Spectrum objects, passed as arguments,
        concatenating them in order, into a new Spectrum object.

        :param args: other Spectrum objects (or a list or tuple of them)
                     to concatenate with this one
        :rtype: Spectrum
        """
        cat_header = self.header
        cat_wl = self.wl
        cat_specI = self.specI
        cat_specV = self.specV
        cat_specN1 = self.specN1
        cat_specN2 = self.specN2
        cat_specSig = self.specSig
        args2 = args
        #Check if the user passed a list or tuple of spectra (unwrap it)
        if len(args) > 0:
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                if isinstance(args[0][0], Spectrum):
                    args2 = args[0]
        for arg in args2:
            if not isinstance(arg, Spectrum):
                raise ValueError('concatenate can only use Spectrum objects!')
            cat_wl = np.concatenate((cat_wl, arg.wl))
            cat_specI = np.concatenate((cat_specI, arg.specI))
            cat_specV = np.concatenate((cat_specV, arg.specV))
            cat_specN1 = np.concatenate((cat_specN1, arg.specN1))
            cat_specN2 = np.concatenate((cat_specN2, arg.specN2))
            cat_specSig = np.concatenate((cat_specSig, arg.specSig))
        cat_spec = Spectrum(cat_wl, cat_specI, cat_specV, cat_specN1,
                            cat_specN2, cat_specSig, header=cat_header)
        return cat_spec

    def coadd(self, *args, byOrders=True, mergeOrders='none'):
        """
        coadd this spectrum with other spectra

        Uses the the wavelength grid of this spectrum for the coadded spectrum.
        Other spectra are interpolated onto this spectrum's wavelengths before
        coadding.  Coadding here essentially averages spectra weighted
        by 1/sigma**2. This assumes the spectra are continuum normalized,
        and have reliable uncertainties.

        Spectra can be split into individual spectral orders and then have each
        order coadded (when byOrders=True) or the orders in the spectra can
        be merged before coadding (when mergeOrders='trim' or ='coadd').
        If neither option is used, and there are regions with order overlap,
        or where wavelength goes backwards, will produce incorrect results
        with this routine. 
        
        :param args: other Spectrum objects (or a list or tuple of them)
                     to coadd with this one
        :param byOrders: flag for whether individual orders are coadded (True),
                     or whole spectra are coadded at once assuming no overlap
                     (False).  If this is False, orders should already be merged,
                     or the mergeOrders flag should be set to 'trim' or 'coadd'.  
        :param mergeOrders: optional method for merging spectral orders before
                     coadding.  Options are 'trim', 'coadd', or 'none'.
                     See the merge_orders function for details.  (If this is
                     not 'trim' or 'coadd' no merging will be done.)
                     If both this flag and byOrders=True are set, then orders
                     are merged before coadding, effectively overriding
                     byOrders=True.
        :rtype: Spectrum
        """
        spec = copy.deepcopy(self) #work on a copy (not self!)
        specList = list(args)
        # Check if the user passed a list or tuple of spectra (unwrap it)
        if len(args) > 0:
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                if isinstance(args[0][0], Spectrum):
                    specList = list(args[0])
        # Type checking
        for arg in specList:
            if not isinstance(arg, Spectrum):
                raise ValueError('coadd can only use Spectrum objects!')
            if np.any(arg.specSig <= 0.0):
                raise ValueError("coadd requires uncertainties "
                                 "(>0) for all points!")

        # Merge orders if requested
        if mergeOrders == 'trim' or mergeOrders == 'coadd':
            spec = spec.merge_orders(mode=mergeOrders)
            for i, arg in enumerate(specList):
                specList[i] = arg.merge_orders(mode=mergeOrders)

        # Split orders if requested (ignore gaps, since only overlaps
        # cause problems for the interpolation routine)
        if byOrders:
            specOrds = spec.get_orders(ignoreGaps=True)
            norders = len(specOrds)
            specListOrds = []
            for arg in specList:
                argOrds = arg.get_orders(ignoreGaps=True)
                if len(argOrds) != norders:
                    raise ValueError('in Spectrum.coadd(), when coadding by '
                                     'orders, all spectra must have the same '
                                     'number of orders!  ')
                specListOrds += [argOrds]
        else: # Wrap spectra in lists (1 element lists with the whole spectrum)
            specOrds = [spec]
            specListOrds = [[arg] for arg in specList]

        # loop over orders adding spectra for each order
        for i, spec in enumerate(specOrds):
            #Set up the weighted average using self as the first entry
            #weight by 1/sigma**2, and save the sum of the weights
            weight = 1./spec.specSig**2
            spec.specI = spec.specI*weight
            spec.specV = spec.specV*weight
            spec.specN1 = spec.specN1*weight
            spec.specN2 = spec.specN2*weight
            totalWeight = weight.copy()
            
            for argOrds in specListOrds:
                arg = argOrds[i] #get this order of this spectrum
                if not np.all(np.diff(arg.wl) > 0):
                    warnings.warn('\nin Spectrum.coadd: Using spectra '
                        'with order overlap or a decreasing wavelength.\n'
                        'Doing this will cause incorrect results in those '
                        'regions!\n(You can try using the mergeOrders="trim" '
                        'or  byOrders=True options to fix this.)',
                        stacklevel=2)

                #Interpolate this spectrum onto the reference wavelengths
                bspecI   = np.interp(spec.wl, arg.wl, arg.specI)
                bspecV   = np.interp(spec.wl, arg.wl, arg.specV)
                bspecN1  = np.interp(spec.wl, arg.wl, arg.specN1)
                bspecN2  = np.interp(spec.wl, arg.wl, arg.specN2)
                bspecSig = np.interp(spec.wl, arg.wl, arg.specSig)

                #Add this spectrum to the weighted sums
                weight = 1./bspecSig**2
                spec.specI += bspecI*weight
                spec.specV += bspecV*weight
                spec.specN1 += bspecN1*weight
                spec.specN2 += bspecN2*weight
                totalWeight += weight

            #Once all spectra have been added, complete the weighted average
            spec.specI /= totalWeight
            spec.specV /= totalWeight
            spec.specN1 /= totalWeight
            spec.specN2 /= totalWeight
            spec.specSig = np.sqrt(1/totalWeight)

        #concatenate orders if necessary back into one Spectrum object
        specOut = specOrds[0]
        if len(specOrds) > 1:
            specOut = specOrds[0].concatenate(specOrds[1:])
        return specOut

    def doppler_shift(self, velocity):
        '''
        Doppler shift the spectrum according to an input radial velocity.

        :param velocity: the radial velocity in km/s
        :rtype: Spectrum
        '''
        c = 299792.458  #speed of light in km/s
        spec = copy.deepcopy(self) #work on a copy (not self!)
        spec.wl = spec.wl + spec.wl*velocity/c
        return spec

    def individual_line(self, lambda0, lwidth):
        '''
        Select an individual line in the spectrum and return it
        as an LSD profile object

        :param lambda0: wavelength of the line (same units as self.wl)
        :param lwidth: distance from the line center, in wavelength,
                       for the wavelength window used for the line profile.
                       One element: same distance on each side of line center.
                       Two elements: distance to the left and right of
                       line center.
        :rtype: LSD
        '''
        #This is nearly a circular import, since profileLSD imports obsSpec
        #Maybe move this to a stand alone function in profileLSD?
        from .profileLSD import LSD
        
        # Select observed line
        if isinstance(lwidth, list) or isinstance(lwidth, tuple):
            if len(lwidth) == 1:
                p_lwidth = [lambda0 - lwidth[0], lambda0 + lwidth[0]]
            elif len(lwidth) == 2:
                p_lwidth = [lambda0 - lwidth[0], lambda0 + lwidth[1]]
            else:
                print('lwidth has too many elements (need one or two)')
                raise ValueError('lwidth has too many elements: '
                                 '{:} (need 1 or 2)'.format(len(lwidth)))
        else:
            p_lwidth = [lambda0 - lwidth, lambda0 + lwidth]
        
        obs_line = self[(self.wl >= p_lwidth[0]) & (self.wl <= p_lwidth[1])]

        # Now we convert wavelengths to velocity space
        c = 299792.458  #speed of light in km/s
        vel = c*(obs_line.wl-lambda0)/lambda0

        prof = LSD(vel, obs_line.specI, obs_line.specSig, obs_line.specV, 
                   obs_line.specSig, obs_line.specN1, obs_line.specSig,
                   header=obs_line.header)
        return prof 

    def save(self, fname, saveHeader=True):
        '''
        Write the Spectrum into a text .s file in a LibreESPRIT style format.
        Optionally skip writing the two lines of header.

        :param saveHeader: optional flag, skip writing the header if False
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
            #Optionally write 2 lines of header            
            if saveHeader:
                if self.header is None:
                    f.write('*** Spectrum of\n')
                else:
                    f.write(self.header)
                    #Make sure there is a line break after the first header text
                    if self.header[-1] != '\n': f.write('\n')
                f.write('{:7d} {:1d}\n'.format(int(self.wl.size), ncols))
            
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

    def get_orders(self, ignoreGaps=False):
        """
        Split an observed echelle spectrum into individual echelle orders.
        
        This splits a concatenated spectrum into echelle orders, based on overlap
        (places where wavelength decreases), and optionally gaps in wavelength
        between orders.

        :param ignoreGaps: if True then regions separated by gaps in wavelength
                           are treated as one order, otherwise gaps are used as
                           order edges
        :return: a list of Spectrum objects, one for each identified order
        """
        #Define spectral order edges by a step backwards in wavelength (velocity)
        #or a step forward in velocity more than 10x the average velocity step size.
        #(use velocity rather than wavelength since pixel size in velocity is
        # more consistent across a spectrum)
        c = 299792.458  #speed of light in km/s
        gapSize = 20.
        velSteps = (self.wl[1:-1] - self.wl[0:-2])/self.wl[1:-1]*c
        meanVelStep = np.mean(velSteps)
        if ignoreGaps:
            orderEdges = velSteps < 0.
        else:
            orderEdges = np.logical_or(velSteps < 0., velSteps > gapSize*meanVelStep)
        indOrderEdges =  np.nonzero(orderEdges)[0] #last point in a spectral order
        numOrders = indOrderEdges.shape[0] + 1

        #Slice the orders out of this spectrum and put them in a new list
        scp = copy.deepcopy(self) #start with a copy (don't just modify self!)
        orders = []
        indLast = 0
        for i in range(numOrders - 1):
            orders.append(scp[indLast : indOrderEdges[i]+1])
            indLast = indOrderEdges[i]+1
        orders.append(scp[indLast :])
        
        return orders

    def get_orders_in_range(self, wl1, wl2=None, ignoreGaps=False):
        """
        Split an observed echelle spectrum into individual orders,
        and return orders that include the specified wavelength range.

        This splits a concatenated spectrum into echelle orders, based on overlap
        (places where wavelength decreases), and optionally gaps in wavelength
        between orders.  Returning the order, or orders, that include any part
        of the specified wavelength range (wl1 and wl2 can be identical).

        :param wl1: the start of the wavelength range to include
        :param wl2: the end of the wavelength range to include
                    (if not provided, orders including just the point in
                    wavelength wl1 will be returned)
        :param ignoreGaps: if True then regions separated by gaps in wavelength
                           are treated as one order, otherwise gaps are used as
                           order edges
        :return: a list of Spectrum objects, one for each relevant order
        """
        if wl2 is None:  #If we just want to check one point in wavelength
            wl2 = wl1
        if wl1 > wl2:
            raise ValueError('Start wavelength must be smaller than end wavelength'
                             ' got {:} {:}'.format(wl1, wl2))
        
        allOrders = self.get_orders(ignoreGaps)
        retOrders = []
        for order in allOrders:
            #if the order starts before the desired range ends,
            #and the order ends after the desired range starts,
            #then some (or all) of it should be in range
            if(wl2 > order.wl[0] and wl1 < order.wl[-1]): 
                retOrders.append(order)

        return retOrders

    def merge_orders(self, mode='trim', midpoint=0.5):
        """
        Simple merging of spectral orders, for echelle spectra.

        The default mode 'trim' simply uses one order up to the midpoint
        of an overlap region, and then uses the next order past the midpoint.
        This is more robust against continuum normalization errors but
        doesn't optimize the total signal-to-noise.
        
        The mode 'coadd' attempts to coadd orders in the overlap region.
        It interpolates the second order onto the wavelength grid of the
        first order, then averages the spectra weighted by 1/sigma**2.
        This optimizes the total signal-to-noise but is vulnerable 
        to continuum normalization errors at the edges of orders.
        This mode requires reliable uncertainties.

        :param mode: choice of 'trim', 'coadd'
        :param midpoint: for mode 'trim', the fraction of the way through
                         an overlap region treated as the midpoint
        :rtype: Spectrum
        """

        orders = self.get_orders(ignoreGaps=True)
        numOrders = len(orders)

        if mode.lower() == 'trim':
            #Merge by splitting orders at the midpoint of their overlap
            wlStartMid = 0.0
            for i in range(numOrders):
                #find the midpoint of the overlap between this and the next order
                if i == numOrders - 1: #for the last order include to the end
                    wlEndMid = orders[i].wl[-1]
                else:
                    wlEndMid = (orders[i].wl[-1]*midpoint
                                + orders[i+1].wl[0]*(1. - midpoint))
                #use this order between this midpoint and the previous midpoint
                induse = (orders[i].wl > wlStartMid) & (orders[i].wl <= wlEndMid)

                #concatenate the used parts of the orders
                if i == 0: #(or just use this order if it's the 1st)
                    specM = orders[i][induse]
                else:
                    specM = specM.concatenate(orders[i][induse])
                wlStartMid = wlEndMid
            
        elif mode.lower() == 'coadd':
            #Merge by coadding orders inside the overlap region
            if np.any(self.specSig <= 0.0):
                raise ValueError("merge_orders in mode 'coadd' requires uncertainties "
                                 "(>0) for all points!")
            for i in range(numOrders):
                #use this order's wavelength grid in the overlap with the next order
                #(assuming all wavelength solutions are equally good)
                if i < numOrders - 1:
                    #identify the overlap range between orders
                    ord1 = orders[i]
                    ord2 = orders[i+1]
                    indO1 = ord1.wl >= ord2.wl[0]
                    indO2 = ord2.wl <= ord1.wl[-1]
                    #Include one extra point from the next order for interpolation
                    indO2[np.nonzero(np.logical_not(indO2))[0][0]] = True
                    
                    #coadd the region in the overlap
                    ord1[indO1] = ord1[indO1].coadd(ord2[indO2], byOrders=False,
                                                    mergeOrders='none')
                    
                #Save this order with coadded values,
                #except for the overlap with the previous order
                if i == 0: #(just use this order if it's the 1st)
                    specM = orders[i]
                else:
                    induse = orders[i].wl > orders[i-1].wl[-1]
                    specM = specM.concatenate(orders[i][induse])            
        else:
            raise ValueError("in merge_orders unrecognized mode '{:}'!".format(mode))
        return specM

    def convolveR(self, R):
        """
        Convolve the spectrum with a Gaussian instrumental profile
        corresponding to a resolution R.  (R is the FWHM of the Gaussian)

        Note: Uncertainties are not propagated in this routine.
        The convolution operation introduces strong correlations into
        the uncertainties for nearby pixels, which is requires careful
        treatment. This routine should not be used for cases where statistics
        or uncertainties are important.  
        
        :param R: the instrumental resolution, unit-less
                  (in the form lambda/delta_lambda)
        :rtype: Spectrum
        """
        if np.all(self.wl[1:] >= self.wl[:-1]):
            specS = self
        else:
            print('Warning: order overlap in spectrum to be convolved!\n'
                  'This will cause incorrect results around overlap regions.\n'
                  '(Sorting pixels by wavelength before convolving.)\n'
                  'Consider using merge_orders before running convolveR.')
            indSort = np.argsort(self.wl)
            specS = self[indSort]

        # check if the polarization spectrum exists
        doPol = (np.any(self.specV != 0.0) | np.any(self.specN1 != 0.0)
                 | np.any(self.specN2 != 0.0))
        
        specC = copy.deepcopy(specS) #output a copy (not self!)
        
        sigmaIntRange = 4.0
        # loop over the spectrum
        for i in range(len(specS)):
            # generate the instrumental profile in wavelength units,
            # for this pixel's wavelength
            fwhm = specS.wl[i]/R
            sigma = fwhm/2.3548200450309493 # FWHM/(2*sqrt(2*log(2)))
            # get the portion of the spectrum relevant for convolution
            # for this pixel, i.e. within n*sigma of this pixel
            ind0, ind1 = np.searchsorted(specS.wl,
                                         [specS.wl[i] - sigmaIntRange*sigma,
                                          specS.wl[i] + sigmaIntRange*sigma])
            specT = specS[ind0:ind1]
            specG = 1./(sigma*np.sqrt(2.*np.pi))*np.exp(-(specT.wl - specS.wl[i])**2/(2.*sigma**2))
            #normG = _trapezoid(specG, wlG)
            #in benchmarking this seems to be ~25% faster than numpy's trapz function
            normG = np.sum((specG[:-1] + specG[1:])*0.5*(specT.wl[1:] - specT.wl[:-1]))
            specG /= normG
            
            # evaluate the convolution for this point
            #specC.specI[i] = _trapezoid(specT.specI*specG, specT.wl)
            prod = specT.specI*specG
            specC.specI[i] = np.sum((prod[:-1] + prod[1:])*0.5*(specT.wl[1:] - specT.wl[:-1]))
            if doPol:
                prod = specT.specV*specG
                specC.specV[i] = np.sum((prod[:-1] + prod[1:])*0.5*(specT.wl[1:] - specT.wl[:-1]))
                prod = specT.specN1*specG
                specC.specN1[i] = np.sum((prod[:-1] + prod[1:])*0.5*(specT.wl[1:] - specT.wl[:-1]))
                prod = specT.specN2*specG
                specC.specN2[i] = np.sum((prod[:-1] + prod[1:])*0.5*(specT.wl[1:] - specT.wl[:-1]))
        return specC


def read_spectrum(fname, trimBadPix=False, sortByWavelength=False):
    """
    Read in the observed spectrum and save it.
    
    This follows the .s format from Donati's LibreESPRIT.
    Files can either have two lines of header or no header.
    This supports 6 column spectropolarimetric files
    (wavelength, I, V|Q|U, null1, null2, errors),
    and also 3 column spectra (wavelength, I, errors).

    :param fname: the name of the file to read.
    :param trimBadPix: optionally remove the more obviously bad pixels if True.
                       Removes pixels with negative flux or error bars of zero,
                       pixels with flux within 3sigma of zero (large errors),
                       and pixels with extremely large values.
    :param sortByWavelength: reorder the points in the spectrum to always
                             increase in wavelength, if set to True.
    :rtype: Spectrum
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
    nLines = 3 - nHeader #we've already read 3 lines
    for line in fObs:
        words = line.split()
        if len(words) == ncolumns:
            nLines += 1
        else:
            print('ERROR: reading observation, '
                  +'line {:}, {:} columns :\n{:}'.format(
                      nLines, len(words), line))

    obs = Spectrum(np.zeros(nLines), np.zeros(nLines), np.zeros(nLines),
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
