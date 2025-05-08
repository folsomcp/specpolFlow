"""
Read a set of FITS files, in the ESPaDOnS spectra format from the CADC
Can be modified for other spectra formats,
if you know the structure of the records in the FITS file.
"""

import numpy as np
import astropy.io.fits as fits
from .. import obsSpec as spf

def espadons(flist, flistout=None, ftype=None, writeSpecHeader=False):
    """
    Convert a list of .fits files in the CADC ESPaDOnS format
    into text .s files.  This function can process both intensity only i.fits
    spectra and intensity+polarization p.fits spectra.

    By default the type of file will be inferred based on the name, if it ends
    with either i.fits or p.fits.  Alternatively, the ftype function parameter
    can be set to either 'i' or 'p' for those respective file types.
    
    This function generates two files in a .s format:  
    
    * The UPENA normalized spectrum (n.s), with automated radial velocity corrections from the telluric lines.
    * The UPENA unnormalized spectrum (u.s), using the automated radial velocity correction from the normalized spectrum.
    
    This is done starting from the unnormalized spectrum without the automated 
    radial velocity correction, to which we apply the radial velocity
    correction determined from the normalized spectrum. The reason behind this
    is that the UPENA automated radial velocity determination performed on 
    unnormalized spectra does not produce consistently reliable results. 
    The content of the fits header is also saved in a .out text file.
    If flistout=None (default), the files are written at the same path as
    the fits-format data, with the '.fits' stripped, 
    and 'n.s' and 'u.s' appended to the filename root. 
    If flistout is a list of paths ending with a rootname, the files will be
    saved at that path with that rootname, and 'n.s' and 'u.s' appended to
    the filename root.

    :param flist: (list of strings) a list of ESPaDOnS filenames
    :param flistout: (list of strings) optional, list of output file rootnames.
                     If this is not given output file names will be based on
                     the input names.
    :param ftype: optional, set the type of .fits file to process.  
                  By default, if this is not set, the type will be inferred
                  from the file name, if the name ends in i.fits or p.fits. 
                  This can be set to the text string 'i' for intensity
                  (no polarization) spectra, or to 'p' for spectropolarimetric
                  (intensity + polarization) spectra.
    :param writeSpecHeader: Flag for whether to write two lines of header in
                            the .s text file (True = include header)
    :return: list of normalized Spectrum objects,
             and a list of unnormalized Spectrum objects
    """

    if isinstance(flist, str): flist = [flist,]
    #Some basic error checking
    if not isinstance(flist, list):
        raise ValueError('in espadons(), the flist argument must be a Python '
                         'list of input file names (or a single file name)')
    if not (flistout is None):
        if isinstance(flistout, str): flistout = [flistout,]
        if not isinstance(flistout, list):
            raise ValueError('in espadons(), the flistout argument must be a '
                             'Python list of output file names')
        if len(flist) != len(flistout):
            raise ValueError('in espadons(), the list of input file names and '
                             'list of output file names must have the same length')

    speclist_n = []
    speclist_u = []
    for i, fname in enumerate(flist):
        #If an output file name exists use it
        if flistout is None:
            fnameOut = None
        else:
            fnameOut = flistout[i]
            
        #If no file type is given, try guessing from the file name
        _ftype = ftype
        if ftype is None:
            if len(fname) > 6:
                if fname[-6:] == 'p.fits':
                    _ftype = 'p'
                elif fname[-6:] == 'i.fits':
                    _ftype = 'i'
                else:
                    print('Warning: failed to infer file type from name:', fname)

        if _ftype == 'p':
            spec_n, spec_u = espadons_p(fname, fnameOut,
                                        writeSpecHeader=writeSpecHeader)
        elif _ftype == 'i':
            spec_n, spec_u = espadons_i(fname, fnameOut,
                                        writeSpecHeader=writeSpecHeader)
        else:
            raise ValueError('in espadons(), unrecognized ftype: '
                             '{:} (only can use: i or p)'.format(ftype))
           
        speclist_n += spec_n
        speclist_u += spec_u
    return speclist_n, speclist_u

    
def espadons_p(flist, flistout=None, writeSpecHeader=False):
    """
    Convert a list of p.fits files in the CADC ESPaDOnS format 
    for spectropolarimetric observations into text .s files. 

    This function generates two files in a .s format:  
    
    * The UPENA normalized spectrum (n.s), with automated radial velocity corrections from the telluric lines.
    * The UPENA unnormalized spectrum (u.s), using the automated radial velocity correction from the normalized spectrum.
    
    This is done starting from the unnormalized spectrum without the automated 
    radial velocity correction, to which we apply the radial velocity
    correction determined from the normalized spectrum. The reason behind this
    is that the UPENA automated radial velocity determination performed on 
    unnormalized spectra does not produce consistently reliable results. 
    The content of the fits header is also saved in a .out text file.
    If flistout=None (default), the files are written at the same path as
    the fits-format data, with the '.fits' stripped, 
    and 'n.s' and 'u.s' appended to the filename root. 
    If flistout is a list of paths ending with a rootname, the files will be
    saved at that path with that rootname, and 'n.s' and 'u.s' appended to
    the filename root.

    :param flist: (list of strings) a list of ESPaDOnS filenames
    :param flistout: (list of strings) optional, list of output file rootnames
    :param writeSpecHeader: Flag for whether to write two lines of header in the
                         .s text file (True = include header)
    :return: list of normalized Spectrum objects,
             and a list of unnormalized Spectrum objects
    """

    if isinstance(flist, str): flist = [flist,]
    #Some basic error checking
    if not isinstance(flist, list):
        raise ValueError('in espadons_p(), the flist argument must be a Python '
                         'list of input file names (or a single file name)')
    if not (flistout is None):
        if isinstance(flistout, str): flistout = [flistout,]
        if not isinstance(flistout, list):
            raise ValueError('in espadons_p(), the flistout argument must be a '
                             'Python list of output file names')
        if len(flist) != len(flistout):
            raise ValueError('in espadons_p(), the list of input file names and '
                             'list of output file names must have the same length')
    
    speclist_n = []
    speclist_u = []
    for i, fname in enumerate(flist):
        print('converting ', fname.strip())
        if flistout is None:
            # striping of white spaces
            fnameOut = fname.strip()
            # removing the '.fits' from the end of the string, 
            # to create the root name for the output files. 
            if len(fnameOut) > 5:
                if fnameOut[-5:] == '.fits' or fnameOut[-5:] == '.FITS':
                    fnameOut = fnameOut[:-5]
        else:
            fnameOut = flistout[i]
        # open the fits file with astropy
        fitsSpec = fits.open(fname.strip())

        header = fitsSpec[0].header

        # Useful for debugging
        #print('File format info')
        #print(fitsSpec.info())
        #print('Header information')
        #print(repr(header))

        # We need to extract the radial velocity correction 
        # that was determined from the normalized spectrum,
        # so that we can apply it to the unnormalized spectrum

        # the radial velocity correction is in the comment section
        # of the header. Using the 'COMMENT' keyword returns a 
        # list of strings
        comments = header['COMMENT']
        # searching for the two strings in the comments that contains
        # the radial velocity correction:
        matches = [match for match in comments if
                   "radial velocity correction from telluric lines" in match]
        # the one calculated for the normalized spectrum is the first instance
        radvel = float(matches[0].split(':')[1])
        
        target = header['OBJECT']
        dateUTC = header['DATE-OB1']
        timeUTC = header['UTIME1']
        # extracting the table of data (24 columns)
        specTab = fitsSpec[0].data

        # The normalized spectrum with radial velocity correction from telluric lines
        # is the first 6 columns
        wln = specTab[0]
        specIn = specTab[1]
        specVn = specTab[2]
        specN1n = specTab[3]
        specN2n = specTab[4]
        specSign = specTab[5]
    
        # The unnormalized spectrum *without* the radial velocity correction
        # from the telluric lines is the last of 4 data blocks
        # (Use this since the radial velocity correction for the unnormalized 
        #  spectrum is unreliable)
        wlu = specTab[18]
        specIu = specTab[19]
        specVu = specTab[20]
        specN1u = specTab[21]
        specN2u = specTab[22]
        specSigu = specTab[23]
    
        fitsSpec.close()

        # Now we apply the velocity correction to the unnormalized data
        c = 299792.458  #speed of light in km/s, since radvel is in km/s
        wlu = wlu + wlu*radvel/c

        spec_n = spf.Spectrum(wln, specIn, specVn, specN1n, specN2n, specSign)
        spec_n.header = '***Reduced spectrum {:} {:} {:}\n'.format(
            target, dateUTC, timeUTC)
        spec_n.save(fnameOut+'n.s', saveHeader=writeSpecHeader)

        spec_u = spf.Spectrum(wlu, specIu, specVu, specN1u, specN2u, specSigu)
        spec_u.header = '***Reduced spectrum {:} {:} {:}\n'.format(
            target, dateUTC, timeUTC)
        spec_u.save(fnameOut+'u.s', saveHeader=writeSpecHeader)
        
        outHeader = open(fnameOut+'.out','w')
        outHeader.write(repr(header))
        outHeader.close()
        
        speclist_n += [spec_n]
        speclist_u += [spec_u]
    return speclist_n, speclist_u


def espadons_i(flist, flistout=None, writeSpecHeader=False):
    """
    Convert a list of i.fits files in the CADC ESPaDOnS format
    for intensity only (no polarization) spectra, into text .s files.
    
    This function generates two files in a .s format:  
    
    * The UPENA normalized spectrum (n.s), with automated radial velocity corrections from the telluric lines.
    * The UPENA unnormalized spectrum (u.s), using the automated radial velocity correction from the normalized spectrum.
    
    This is done starting from the unnormalized spectrum without the automated 
    radial velocity correction, to which we apply the radial velocity
    correction determined from the normalized spectrum. The reason behind this
    is that the UPENA automated radial velocity determination performed on 
    unnormalized spectra does not produce consistently reliable results. 
    The content of the fits header is also saved in a .out text file.
    If flistout=None (default), the files are written at the same path as
    the fits-format data, with the '.fits' stripped, 
    and 'n.s' and 'u.s' appended to the filename root. 
    If flistout is a list of paths ending with a rootname, the files will be
    saved at that path with that rootname, and 'n.s' and 'u.s' appended to
    the filename root.

    :param flist: (list of strings) a list of ESPaDOnS filenames
    :param flistout: (list of strings) optional, list of output file rootnames
    :param writeSpecHeader: Flag for whether to write two lines of header in the
                         .s text file (True = include header)
    :return: list of Spectrum objects
    """

    if isinstance(flist, str): flist = [flist,]
    #Some basic error checking
    if not isinstance(flist, list):
        raise ValueError('in espadons_i(), the flist argument must be a Python '
                         'list of input file names (or a single file name)')
    if not (flistout is None):
        if isinstance(flistout, str): flistout = [flistout,]
        if not isinstance(flistout, list):
            raise ValueError('in espadons_i(), the flistout argument must be a '
                             'Python list of output file names')
        if len(flist) != len(flistout):
            raise ValueError('in espadons_i(), the list of input file names and '
                             'list of output file names must have the same length')
    
    speclist_n = []
    speclist_u = []
    for i, fname in enumerate(flist):
        print('converting ', fname.strip())
        if flistout is None:
            # striping of white spaces
            fnameOut = fname.strip()
            # removing the '.fits' from the end of the string, 
            # to create the root name for the output files.
            if len(fnameOut) > 5:
                if fnameOut[-5:] == '.fits' or fnameOut[-5:] == '.FITS':
                    fnameOut = fnameOut[:-5]
        else:
            fnameOut = flistout[i]
        # open the fits file with astropy
        fitsSpec = fits.open(fname.strip())

        header = fitsSpec[0].header

        # Useful for debugging
        #print('File format info')
        #print(fitsSpec.info())
        #print('Header information')
        #print(repr(header))

        # We need to extract the radial velocity correction 
        # that was determined from the normalized spectrum,
        # so that we can apply it to the unnormalized spectrum

        # the radial velocity correction is in the comment section
        # of the header. Using the 'COMMENT' keyword returns a 
        # list of strings
        comments = header['COMMENT']
        # searching for the two strings in the comments that contains
        # the radial velocity correction:
        matches = [match for match in comments if
                   "radial velocity correction from telluric lines" in match]
        # the one calculated for the normalized spectrum is the first instance
        radvel = float(matches[0].split(':')[1])

        target = header['OBJECT']
        dateUTC = header['DATE-OBS']
        timeUTC = header['UTC-OBS']
        # extracting the table of data (12 columns)
        specTab = fitsSpec[0].data
        
        # The normalized spectrum with radial velocity correction from telluric lines
        # is the first 3 columns
        wln = specTab[0]
        specIn = specTab[1]
        specSign = specTab[2]
    
        # The unnormalized spectrum *without* the radial velocity correction
        # from the telluric lines is the last of 4 data blocks
        # (Use this since the radial velocity correction for the unnormalized 
        #  spectrum is unreliable)
        wlu = specTab[9]
        specIu = specTab[10]
        specSigu = specTab[11]
    
        fitsSpec.close()

        # Now we apply the velocity correction to the unnormalized data
        c = 299792.458  #speed of light in km/s, since radvel is in km/s
        wlu = wlu + wlu*radvel/c

        zeros = np.zeros_like(specIn)
        spec_n = spf.Spectrum(wln, specIn, zeros, zeros, zeros, specSign)
        spec_n.header = '***Reduced spectrum {:} {:} {:}\n'.format(
            target, dateUTC, timeUTC)
        spec_n.save(fnameOut+'n.s', saveHeader=writeSpecHeader)

        spec_u = spf.Spectrum(wlu, specIu, zeros, zeros, zeros, specSigu)
        spec_u.header = '***Reduced spectrum {:} {:} {:}\n'.format(
            target, dateUTC, timeUTC)
        spec_u.save(fnameOut+'u.s', saveHeader=writeSpecHeader)
        
        outHeader = open(fnameOut+'.out','w')
        outHeader.write(repr(header))
        outHeader.close()
        
        speclist_n += [spec_n]
        speclist_u += [spec_u]
    return speclist_n, speclist_u
