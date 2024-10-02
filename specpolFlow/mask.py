"""
Tools for generating and manipulating line masks for LSD.
Includes creating a line mask from a VALD line list, cleaning problem 
lines from the mask, managing regions to exclude from the mask,
and interacively cleaning and tweaking.
"""

import numpy as np
import warnings
#from . import lineList as lineListLib #(moved inside some functions)

###################################
###################################

class Mask:
    """
    The data for the LSD line Mask.

    This usually contains the arrays: 

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
        """
        Returns a Mask object with only the values at the specified index(s).

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
        Sets all values of the mask at the specified location
        equal to the input mask's values.

        :param key: the index or slice being overwritten
        :param newval: Mask object used to replace the overwritten values
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
        """
        Returns the number of lines in the mask.
        """
        return len(self.wl)
    
    def prune(self):
        """
        Return a Mask object with unused lines removed from the Mask.
        
        Remove lines if iuse index is set to 0,
        restricting the Mask to only lines used in LSD.
        
        :rtype: Mask
        """
        trimmedMask = self[self.iuse != 0]
        return trimmedMask
        
    def get_weights(self, normDepth, normWave, normLande):
        """
        Returns the calculated LSD weights of all the lines in the mask
        (includes lines with both the iuse flag 0 and 1).
        
        This assumes the Stokes I lines are weighted as depth,
        and Stokes V is weighted as depth*wavelength*Lande factor.
        
        :param normDepth: the normalizing line depth for the mask/LSD profile
        :param normWave: the normalizing wavelength (nm) for
                         the mask/LSD profile
        :param normLande: the normalizing effective Lande factor for
                          the mask/LSD profile
        :return: weightI, weightV (as arrays)
        """
        weightI = self.depth / normDepth
        weightV = self.depth*self.wl*self.lande / (normDepth*normWave*normLande)
        return weightI, weightV
    
    def save(self, fname):
        """
        Save the line Mask to a text file, in Donati's and LSDpy format.
        
        :param fname: the file name to output the Mask to
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
        Remove lines inside a set of exclude regions from a mask.
        
        Takes an ExcludeMaskRegions object, and returns a Mask object with
        the iuse values set to zero for spectral lines with wavelengths
        inside these regions. 
        
        :param regions: An ExcludeMaskRegions object
        :rtype: Mask
        '''
        # Making a copy of the mask (or use copy.deepcopy(self) )
        mask_clean = Mask(self.wl.copy(), self.element.copy(),
                          self.depth.copy(), self.excite.copy(),
                          self.lande.copy(), self.iuse.copy())
        nregions = len(regions)
        for i in range(0,nregions):
            is_in = np.logical_and( (self.wl>=regions[i].start),
                                    (self.wl<=regions[i].stop) )
            mask_clean.iuse[is_in] = 0

        return mask_clean

def read_mask(fname):
    """
    Read in an LSD line mask file and return a Mask object.

    The mask file will have one line of header and columns of:  

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
    ind = np.argsort(tmpMask[0,:], kind='stable')
    
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
    """
    Class for a region object that records spectral regions to exclude from a Mask.

    Usually contrains arrays of:  
    
    * start - starting wavelengths for the regions
    * stop - ending wavelengths for the regions
    * type - optionally, text comments for the type of region \
             (inside a numpy array with dtype=object)
    """

    def __init__(self, start, stop, type):
        self.start = start
        self.stop = stop
        self.type = type

    def __getitem__(self, key):
        """
        Returns a region object with only the values at the specified index(s).

        :param key: the index or slice being checked
        :rtype: ExcludeMaskRegions
        """
        return ExcludeMaskRegions(self.start[key], self.stop[key],
                                  self.type[key])

    def __setitem__(self, key, newval):
        """
        Sets all values of the Mask at the specified location
        equal to the input mask's values.

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
        """
        Returns the number of lines in the mask
        """
        return len(self.start)

    def __add__(self, other):
        """
        Concatenate two ExcludeMaskRegions objects.
        
        :param other: an ExcludeMaskRegions object to combine with this one
        :rtype: ExcludeMaskRegions
        """
        start = np.concatenate([self.start, other.start])
        stop = np.concatenate([self.stop, other.stop])
        type = np.concatenate([self.type, other.type])
        return ExcludeMaskRegions(start, stop, type)

    def save(self, fname):
        """
        Save the ExcludeMaskRegions object to a text file.

        The file contains one line for each region, with the start wavelength,
        end wavelength, and a type comment.
        
        :param fname: the file path/name
        """
        with open(fname, 'w') as ofile:
            for item in self:
                ofile.write('{:.4f} {:.4f} {}\n'.format(item.start, item.stop,
                                                        item.type))
        return

    def to_dict(self):
        """
        Return the ExcludeMaskRegions as a dictionary. 
        This is useful to transform ExcludeMaskRegions objects
        to Panda dataframes.
        """
        return {'start':self.start, 'stop':self.stop,'type':self.type}

def read_exclude_mask_regions(fname):
    """
    Read in an ExcludeMaskRegions file into an ExcludeMaskRegions object. 

    The file should be a text file with one line for each region,
    containing the start wavelength, end wavelength, and the type comment.
    
    :param fname: the path/name of the file
    :rtype: ExcludeMaskRegions
    """
    start = []
    stop = []
    type = []
    with open(fname, 'r') as f:
        nLines = 0
        for txtLine in f:
            split = txtLine.split() # split by WS
            start.append(split[0]) 
            stop.append(split[1])
            # join all the rest with ' ' in case there was multiple words
            type.append(' '.join(split[2:]))
    
    start = np.asarray(start, dtype=float)
    stop = np.asarray(stop, dtype=float)
    
    return ExcludeMaskRegions(start, stop, np.array(type, dtype=object))

def get_Balmer_regions_default(velrange=500):
    '''
    Generate an ExcludeMaskRegions object with regions around Balmer H-lines
    (alpha to epsilon) up to a given radial velocity,
    and a region that excludes the Balmer jump (from 360-392 nm).
    
    :param velrange: velocity range around the H-line to be excluded
                     (default 500 km/s)
    :rtype: ExcludeMaskRegions
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
    start.append(360.0)
    stop.append(392.0)

    return ExcludeMaskRegions(np.array(start), np.array(stop),
                              np.array(types,dtype=object))

def get_telluric_regions_default():
    '''
    Generate an ExcludeMaskRegions object with regions containing
    heavy telluric line contamination in the optical.
    
    :rtype: ExcludeMaskRegions
    '''
    start = np.array([587.5,627.5,684.0,717.0,757.0,790.0,809.0])  # nm
    stop   = np.array([592.0,632.5,705.3,735.0,771.0,795.0,990.0])  # nm

    return ExcludeMaskRegions(start, stop,
                              np.array(['telluric']*len(start),dtype=object))


###################################
###################################

def make_mask(lineListFile, outMaskName=None, depthCutoff=0.0, 
              wlStart=None, wlEnd=None, landeStart=None, landeEnd=None,
              elementsUsed=[], elementsExclude=[],
              atomsOnly=True, includeNoLande=False, defaultLande=1.0):
    """
    Generate an LSD line mask from a VALD line list and save it.
    This is the main interface for generating new line masks with Python scrips.
    
    This generates an LSD line mask from atomic lines in the list,
    and estimates Lande factors for lines where possible.  Lines for which 
    no Lande factor can be estimated are omitted by default.
    This assumes input line list wavelengths in A, and returns a mask in nm.
    
    :param lineListFile: The name of the file containing the line list
                         (in the VALD3 'extract stellar' 'long' format)
    :param outMaskName: The name of the to write the mask to
                        (set this to None to avoid saving a file, default)
    :param depthCutoff: Only include lines in the mask deeper than this value
                        (defaults to 0, all lines included)
    :param wlStart: Optionally, only use lines with wavelengths above this
                    (note: value in nm!)
    :param wlEnd: Optionally, only use lines with wavelengths below this
                  (note: value in nm!)
    :param landeStart: Optionally, only use lines with effective Lande factors
                       above this
    :param landeEnd: Optionally, only use lines with effective Lande factors
                     below this
    :param elementsUsed: Optionally, provide a list of elements to include
                         in the line mask. This should be a list of strings
                         with the element symbols, e.g. ['C', 'O', 'Si', 'Fe']
                         (an empty list or a list starting with 'all' will
                         include all elements)
    :param elementsExclude: Optionally, provide a list of elements to exclude
                            from the line mask. Also a list of strings with
                            the element symbols.
    :param atomsOnly: Only include atomic lines (no molecular lines) and
                      exclude H lines (defaults to True)
    :param includeNoLande: Include lines in the mask even if no Lande factor
                           can be estimated for them (defaults to False)
    :param defaultLande: The value of the effective Lande factor to use if
                         no value can be estimated (only used if
                         includeNoLande = True)
    :rtype: Mask
    """

    from . import lineList as lineListLib
    lineList = lineListLib.read_VALD(lineListFile)
    
    mask = convert_list_to_mask(lineList, depthCutoff=depthCutoff,
                                atomsOnly=atomsOnly,
                                includeNoLande=includeNoLande,
                                defaultLande=defaultLande)
    
    mask = filter_mask(mask, depthCutoff=depthCutoff,
                       wlStart=wlStart, wlEnd=wlEnd,
                       landeStart=landeStart, landeEnd=landeEnd,
                       elementsUsed=elementsUsed,
                       elementsExclude=elementsExclude)

    if outMaskName is not None:
        mask.save(outMaskName)
    return mask

    
def convert_list_to_mask(lineList, depthCutoff=0.0, atomsOnly = True,
                         includeNoLande = False, defaultLande=1.0):
    """
    Convert a lineList to an LSD line mask.
    
    This estimates Lande factors when VALD doesn't have a known value.
    This can also automatically reject any H lines and any molecular lines.
    Note, this assumes an input line list in A, and returns a mask in nm.
    See also the make_mask() function.
    
    :param lineList: The line list object to be converted to a mask
    :param depthCutoff: Only include lines deeper than this depth cutoff
                        (defaults to 0, all lines included)
    :param atomsOnly: Flag to use only atomic lines (no molecular lines)
                      and exclude H lines (defaults to True)
    :param includeNoLande: Flag to include lines in the mask even if no Lande
                           factor can be estimated for them (defaults to False)
    :param defaultLande: The value of the effective Lande factor to use for
                         lines where no value can be estimated
                         (only used if includeNoLande = True)
    :rtype: Mask
    """
    from . import lineList as lineListLib
    elements=('H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
              'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
              'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
              'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
              'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
              'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
              'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
              'At','Rn','Fr','Ra','Ac','Th','Pa','U' )
    nLines = lineList.nLines
    mask = Mask(np.zeros(nLines), np.zeros(nLines), np.zeros(nLines),
                np.zeros(nLines), np.zeros(nLines),
                np.zeros(nLines,dtype=int))

    skippedIon = []
    missingLandeIon = []
    missingLandeNum = 0
    j = 0
    #Loop through the line list making sure to use good data
    for i in range(lineList.nLines):
        element = lineList.ion[i].split()[0]
        ion = lineList.ion[i].split()[1]

        if (lineList.depth[i] >= depthCutoff) and (
            (atomsOnly == False) or (element in elements and element != 'H')):
            #Estimate an effective Lande factor if necessary, VALD flags 
            #levels without a known Lande factor using the value 99.0
            landeEff = lineList.landeEff[i]
            if lineList.landeEff[i] >= 90.:
                landeLo = lineList.landeLo[i]
                if landeLo >= 90.:
                    landeLo = lineListLib.estimateLande(
                        lineList.Jlo[i], lineList.configLo[i])
                landeUp = lineList.landeUp[i]
                if landeUp >= 90.:
                    landeUp = lineListLib.estimateLande(
                        lineList.Jup[i], lineList.configUp[i])
                #If we have managed to estimate Lande factors
                if landeLo < 90. and landeUp < 90.:
                    landeEff = lineListLib.getEffectiveLande(
                        landeLo, landeUp, lineList.Jlo[i], lineList.Jup[i])

            #If we can't calculate a Lande factor but still want to use the line
            if landeEff >= 90. and includeNoLande == True:
                landeEff = defaultLande

            #For lines with a Lande factor
            if landeEff < 90.:
                #Save the data for good lines
                mask.wl[j]     = lineList.wl[i]
                mask.depth[j]  = lineList.depth[i]
                mask.excite[j] = lineList.Elo[i]
                mask.lande[j]  = landeEff
                mask.iuse[j]   = 1
                #The element code here is the atomic number plus 
                #the ionization state/100, so Fe II becomes 26.01
                try:
                    mask.element[j] = float(elements.index(element)+1) \
                        + (float(ion)-1)*0.01
                except ValueError:
                    mask.element[j] = 0.00
                
                j += 1
            else:
                #Track lines not used due to missing Lande factors
                missingLandeNum += 1
                if not lineList.ion[i] in missingLandeIon:
                    missingLandeIon += [lineList.ion[i]]
        elif (lineList.depth[i] >= depthCutoff):
            #Track lines not used due unsupported species
            if not lineList.ion[i] in skippedIon:
                skippedIon += [lineList.ion[i]]

    if missingLandeNum > 0:
        print('missing Lande factors for {:} lines (skipped) from:'.format(
            missingLandeNum))
        print(missingLandeIon)
    if len(skippedIon) > 0:
        print('skipped all lines for species:')
        print(skippedIon)
    
    #Convert wavelengths from A to values in nm
    mask.wl = mask.wl*0.1
    #Remove blank entries in the line mask
    mask = mask.prune()
    
    return mask


def filter_mask(mask, depthCutoff = 0.0, wlStart = None, wlEnd = None,
               landeStart = None, landeEnd = None,
               elementsUsed = [], elementsExclude = []):
    """
    Remove lines from an LSD line mask, based on their depth, wavelength,
    effective Lande factor, or element.

    Usually, only one of elementsUsed and elementsExclude should be used.
    Both can be used, in that case only elements in elementsUsed are included,
    and then any elements in elementsExclude are removed from that.
    (Filtering by element currently only works for atomic lines,
    not molecular lines)

    :param mask: The line mask object to filter lines out of
    :param depthCutoff: Only include lines in the mask deeper than this value
                        (defaults to 0, all lines included)
    :param wlStart: Optionally, only use lines with wavelengths above this
                    (note: value in nm!)
    :param wlEnd: Optionally, only use lines with wavelengths below this
                  (note: value in nm!)
    :param landeStart: Optionally, only use lines with effective Lande factors
                       above this
    :param landeEnd: Optionally, only use lines with effective Lande factors
                     below this
    :param elementsUsed: Optionally, provide a list of elements to include
                         in the line mask. This should be a list of strings
                         with the element symbols, e.g. ['C', 'O', 'Si', 'Fe']
                         (an empty list or a list starting with 'all' will
                         include all elements)
    :param elementsExclude: Optionally, provide a list of elements to exclude
                            from the line mask. Also a list of strings with
                            the element symbols.
    :rtype: Mask
    """
    
    elements=('H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
              'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
              'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
              'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
              'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
              'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
              'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
              'At','Rn','Fr','Ra','Ac','Th','Pa','U' )

    if wlStart is None: wlStart = 0.0
    if wlEnd is None: wlEnd = 1e10
    if landeStart is None: landeStart = -1e10
    if landeEnd is None: landeEnd = 1e10,
    
    #Remove lines that are outside the requested range in wavelength, Lande factor, or depth.
    indRemove = ((mask.depth < depthCutoff)
                 | (mask.wl < wlStart) | (mask.wl > wlEnd) 
                 | (mask.lande < landeStart) | (mask.lande > landeEnd))
    
    #If there is also a list of specific elements to use
    if type(elementsUsed) != type([]): print('ERROR: list of elements to include not of type list')
    if elementsUsed != [] and elementsUsed[0].lower() != 'all':
        indKeep = np.zeros_like(mask.iuse, dtype=bool)
        #Set a flag to only keep elements in the requested list
        for elUse in elementsUsed:
            try: #get the atomic number for this element
                elNumber = elements.index(elUse.strip())+1
            except ValueError:
                raise ValueError('Unrecognized element {:} in elementsUsed'.format(elUse))
            indKeep = indKeep | (mask.element.astype(int) == elNumber)
        indRemove = indRemove | np.logical_not(indKeep)
    
    #Or if there is also a list of specific elements to exclude from the mask
    if type(elementsExclude) != type([]): print('ERROR: list of elements to exclude not of type list')
    if elementsExclude != [] :
        #if elementsUsed != [] and elementsUsed[0].lower() != 'all':
        #    print('Warning: got both a list of elements to include and to exclude. ')
        indExclude = np.zeros_like(mask.iuse, dtype=bool)
        #Set a flag to remove elements in the exclude list
        for elCut in elementsExclude:
            try: #get the atomic number for this element
                elNumber = elements.index(elCut.strip())+1
            except ValueError:
                raise ValueError('Unrecognized element {:} in elementsExclude'.format(elUse))
            indExclude = indExclude | (mask.element.astype(int) == elNumber)
        indRemove = indRemove | indExclude
    
    #Apply the filtering
    mask.iuse[indRemove] = 0
    
    #Remove blank entries in the line mask
    mask = mask.prune()
    
    return mask
