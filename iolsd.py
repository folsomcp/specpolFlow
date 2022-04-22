## @module iolsd.py
# Documentation for iolsd.py
#
# More information can come here.

import numpy as np

class lsd_prof:
    def __init__(self, velStart=None, velEnd=None, pixVel=None):
        """
        Initialize an empty LSD profile.
        
        Optionally this can set up the grid in velocity for the profile.

        :param velStart: the velocity for the start (blue side) of the LSD profile
        :param velEnd: the velocity for the end (red side) of the the LSD profile
        :param pixVel: the step in velocity between pixels in the LSD profile
        """
        if velStart != None and velEnd != None and pixVel != None:
            self.vel = np.arange(velStart, velEnd+pixVel, pixVel)
            self.npix = self.vel.shape[0]
        else:
            self.npix = 0
            self.vel = np.zeros(self.npix)
        self.specI = np.ones(self.npix)
        self.specSigI = np.zeros(self.npix)
        self.specV = np.zeros(self.npix)
        self.specSigV = np.zeros(self.npix)
        self.specN1 = np.zeros(self.npix)
        self.specSigN1 = np.zeros(self.npix)
        self.specN2 = np.zeros(self.npix)
        self.specSigN2 = np.zeros(self.npix)
        self.header = None

    def save(self, fname):
        """
        Save the LSD profile to a file.
        
        This saves to a text file in Donati's format.
        
        :param fname: the name of the file the LSD profile to save to.
        :param header: optionally, one line of header text for the output file.  This text string should end with \n.
        """
        
        oFile = open(fname, 'w')
        if self.header != None:
            oFile.write(self.header)
            oFile.write(' {:d} 6\n'.format(self.npix))
        
        for i in range(self.npix):
            oFile.write('{:>12.6f} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e} {:>13.6e}\n'.format(
                self.vel[i], self.specI[i], self.specSigI[i], self.specV[i],
                self.specSigV[i], self.specN1[i], self.specSigN1[i]))
        oFile.close()
        return


def read_lsd(fname):
    """
    function that reads in a LSD profile.
    
    The LSD profiles are in Donati's text format.
    The two lines of header in Donati's format is optional
    
    :param fname: the name of the file containing the LSD profile
    :rtype: returns an instance of the lsd_prof class, defined in this module.
    """
    #check if this file has a header
    fcheck = open(fname, 'r')
    head1txt = fcheck.readline()
    head2txt = fcheck.readline()
    fcheck.close()
    head1 = head1txt.split()
    head2 = head2txt.split()
    nskip = 2
    try:
        float(head1[0])
        if(len(head2) > 2 and len(head1) == len(head2)):
            nskip = 0
    except(ValueError):
        nskip = 2
    
    #Read the LSD profile, skipping header lines
    __prof = np.loadtxt(fname, skiprows=nskip, unpack=True)

    prof = lsd_prof()
    prof.vel       = __prof[0,:]
    prof.specI     = __prof[1,:]
    prof.specSigI  = __prof[2,:]
    prof.specV     = __prof[3,:]
    prof.specSigV  = __prof[4,:]
    prof.specN1    = __prof[5,:]
    prof.specSigN1 = __prof[6,:]
    
    prof.npix = prof.vel.shape[0]
    prof.specN2 = np.zeros(prof.npix)
    prof.specSigN2 = np.zeros(prof.npix)
    if (nskip > 0): prof.header = head1txt
    
    return prof


class mask:
    def __init__(self, fname=None):
        """
        Read in an LSD mask file and save it to an instance of the mask class.
        
        The mask files are in the format for Donati's LSD and LSDpy.
        
        :param fname: the name of the file containing the mask
        """
        if fname == None:
            self.wl = np.zeros(0)
            self.element = np.zeros(0)
            self.depth = np.zeros(0)
            self.lande = np.zeros(0)
            self.iuse = np.zeros(0)
            
        else:
            #Columns are: wavelength (nm),
            #             element + ionization*0.01,
            #             line depth,
            #             excitation potential of the lower level,
            #             effective Lande factor,
            #             flag for whether the line is used (1=use, 0=skip).
            tmpMask = np.loadtxt(fname, skiprows=1, unpack=True)
            
            #Sort the line mask so wavelength is always increasing
            self.ind = np.argsort(tmpMask[0,:])
            
            self.wl = tmpMask[0, self.ind]
            self.element = tmpMask[1, self.ind]
            self.depth = tmpMask[2, self.ind]
            self.excite = tmpMask[3, self.ind]
            self.lande = tmpMask[4, self.ind]
            self.iuse = tmpMask[5, self.ind].astype(int)

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


class observation:
    def __init__(self, fname, sortByWavelength=False):
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
                self.header = None
                fObs.seek(0)
            else:
                #Otherwise assume there are two lines of header,
                #the first one being a comment the second should be
                #dimensions of the file but we figure that by reading the file.
                self.header = line
        except ValueError:
            self.header = line
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

        self.wl = np.zeros(nLines)
        self.specI = np.zeros(nLines)
        self.specV = np.zeros(nLines)
        self.specN1 = np.zeros(nLines)
        self.specN2 = np.zeros(nLines)
        self.specSig = np.zeros(nLines)
        
        i = 0
        #Rewind to start then advance the file pointer 2 lines
        fObs.seek(0)
        if self.header != None:
            fObs.readline()
            fObs.readline()
        #Then read the actual data of the file
        for line in fObs:
            words = line.split()
            if (len(words) == ncolumns and ncolumns == 6):
                if(float(words[1]) > 0. and float(words[5]) > 0.):
                    self.wl[i] = float(words[0])
                    self.specI[i] = float(words[1])
                    self.specV[i] = float(words[2])
                    self.specN1[i] = float(words[3])
                    self.specN2[i] = float(words[4])
                    self.specSig[i] = float(words[5])
                    i += 1
            elif (len(words) == ncolumns and ncolumns == 3):
                if(float(words[1]) > 0. and float(words[2]) > 0.):
                    self.wl[i] = float(words[0])
                    self.specI[i] = float(words[1])
                    self.specSig[i] = float(words[2])
                    self.specV[i] = 0.
                    self.specN1[i] = 0.
                    self.specN2[i] = 0.
                    i += 1
                
        fObs.close()
        
        #Optionaly, sort the observation so wavelength is always increasing
        if sortByWavelength:
            self.ind = np.argsort(self.wl)
            
            self.wl = self.wl[self.ind]
            self.specI = self.specI[self.ind]
            self.specV = self.specV[self.ind]
            self.specN1 = self.specN1[self.ind]
            self.specN2 = self.specN2[self.ind]
            self.specSig = self.specSig[self.ind]
        
        return
