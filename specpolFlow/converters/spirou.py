"""
Converters that read a set of FITS files, in SPIRou APERO spectrum formats,
and save them to text .s formats.
"""

import numpy as np
import astropy.io.fits as fits
from .. import obsSpec as spf

def spirou(flist, flistout=None, ftype=None, nanTreatment='replace',
           removeSegmentSize=200, allowGapSize=3,
           saveFitsHeader=True, writeSpecHeader=False):
    '''
    Convert a .fits, file or list of files, from SPIRou APERO formats into
    text .s files.

    This supports intensity spectrum e.fits files, telluric corrected intensity
    spectrum t.fits files, and polarized spectrum p.fits files.
    This calls the spirou_e, spirou_t, and spirou_p
    functions.

    The .fits files usually contain nan values for pixels where the telluric
    correction failed, or where the spectrum extraction failed.
    This function provides options: to remove those pixels (and optionally 
    remove some of the small fragments of spectrum that get left behind),
    to replace those pixels with placeholder values, or to keep the
    nan values. (Keeping the nan values will cause problems for other
    SpecpolFlow functions!)

    For e.fits and t.fits files: this function estimates uncertainties
    on the spectrum simply as sqrt(N), applies a blaze correction to
    the spectrum, and applies a barycentric velocity correction.

    For p.fits files: this function converts the input polarized spectrum from
    being normalized by the intensity spectrum (e.g. V/I) to being normalized
    by the continuum (e.g. V/Ic, like I/Ic).

    :param flist: .fits file name or list filenames (a string or list of strings)
    :param flistout: optional, an output filename or list of filenames.
                     If this is not provided output names will be generated
                     from the input filenames, replacing .fits with .s
    :param ftype: the type of input .fits file to process.
                  'e' is e.fits files, with intensity spectra, still with
                  separate spectral orders.
                  't' is t.fits files, with telluric corrected spectra,
                  still with separate spectral orders.
                  'p' is p.fits files, with polarized spectra.
                  If this is not provided, the file type will be estimated
                  from the input filename, but the files must end with
                  e.fits, t.fits, or p.fits.
    :param nanTreatment: Flag for how to deal with nan values in the input file,
                         can be: 'remove', 'replace', 'keep'. 
                         nan values can be from places where the spectrum 
                         extraction failed or the telluric correction was too 
                         uncertain.  'remove' deletes pixels with nans.
                         'replace' replaces the nan values with 0, and sets the
                         uncertainties to 100. 
                         'keep' retains the nan values in the spectrum.
    :param removeSegmentSize: If nanTreatment = 'remove', this can leave a lot
                         of little fragments of spectrum in the output, mostly
                         in regions with strong telluric lines.  
                         Fragments that are shorter than removeSegmentSize 
                         pixels will be removed from the output. (Set this to
                         0 to keep all the spectrum)
    :param allowGapSize: If nanTreatment = 'remove', and small fragments
                         of spectrum are also being removed, there can be small
                         gaps of a few pixels in the spectrum. Setting
                         allowGapSize ignores gaps smaller than this value when
                         calculating the size of segments for removeSegmentSize.
    :param saveFitsHeader: Flag for whether to save the .fits header information
                         to a text file, named using the output filename
                         .out (True = save to file).
    :param writeSpecHeader: Flag for whether to write two lines of header in the
                         .s text file (True = include header)
    :return: list of Spectrum objects
    '''
    if isinstance(flist, str): flist = [flist,]
    #Some basic error checking
    if not (isinstance(flist, list) or isinstance(flist, tuple)):
        raise ValueError('in spirou(), flist must be a string '
                         'or list of strings')
    if not flistout is None:
        if isinstance(flistout, str): flistout = [flistout,]
        if not (isinstance(flistout, list) or isinstance(flistout, tuple)):
            raise ValueError('in spirou(), flistout must be a string '
                             'or list of strings with output file names')
        if len(flist) != len(flistout):
            raise ValueError('in spirou(), the list of input file names and '
                             'list of output file names must have the same length')
    
    speclist = []
    for i, fname in enumerate(flist):
        #If an output file name exists use it
        if flistout is None:
            outname = None
        else:
            outname = flistout[i]
        #If no file type is given, try guessing from the file name
        if ftype is None:
            if len(fname) > 6:
                if fname[-6:] == 'p.fits':
                    _ftype = 'p'
                elif fname[-6:] == 'e.fits':
                    _ftype = 'e'
                elif fname[-6:] == 't.fits':
                    _ftype = 't'
        else:
            _ftype = ftype
        
        if _ftype == 'p':
            spec = spirou_p(fname, outname, nanTreatment=nanTreatment,
                            removeSegmentSize=removeSegmentSize,
                            allowGapSize=allowGapSize,
                            saveFitsHeader=saveFitsHeader,
                            writeSpecHeader=writeSpecHeader)
        elif _ftype == 'e':
            spec = spirou_e(fname, outname, nanTreatment=nanTreatment,
                            removeSegmentSize=removeSegmentSize,
                            allowGapSize=allowGapSize,
                            saveFitsHeader=saveFitsHeader,
                            writeSpecHeader=writeSpecHeader)
        elif _ftype == 't':
            spec = spirou_t(fname, outname, nanTreatment=nanTreatment,
                            removeSegmentSize=removeSegmentSize,
                            allowGapSize=allowGapSize,
                            saveFitsHeader=saveFitsHeader,
                            writeSpecHeader=writeSpecHeader)
        else:
            raise ValueError('in spirou(), unrecognized ftype: '
                             '{:} (only can use: p t e)'.format(ftype))

        speclist += [spec]
    return speclist


def spirou_p(fname, outname=None, nanTreatment='replace',
             removeSegmentSize=200, allowGapSize=3,
             saveFitsHeader=True, writeSpecHeader=False):
    '''
    Function to read in a SPIRou p.fits file containing a polarized spectrum,
    save it to a .s file, and return a Spectrum object

    The p.fits files usually contain the polarization spectrum normalized by
    the full Stokes I spectrum, not the continuum of Stokes I
    (e.g. V/I not V/Ic).  This function converts the polarized spectrum to be
    normalized by the continuum (outputting e.g. V/Ic).

    The p.fits files usually contain nan values for pixels where the telluric
    correction failed, or where the spectrum extraction failed.
    This function provides options: to remove those pixels (and optionally 
    remove some of the small fragments of spectrum that get left behind),
    to replace those pixels with placeholder values, or to keep the
    nan values. (Keeping the nan values will cause problems for other
    SpecpolFlow functions!)
    APERO p.fits files also have nan values for wavelengths, so in the 
    'replace' mode wavelengths in nan regions are only estimates, 
    based on the surrounding wavelength values.
    
    :param fname: name of the .fits file to read
    :param outname: name of the file to save the spectrum to.  If not provided
                    the function will default to using the input fname and
                    replacing .fits with .s
    :param nanTreatment: Flag for how to deal with nan values in the input file,
                         can be: 'remove', 'replace', 'keep'. 
                         nan values can be from places where the spectrum 
                         extraction failed or the telluric correction was too 
                         uncertain. 'remove' deletes pixels with nans.
                         'replace' replaces the nan values with 0, and sets the
                         uncertainties to 100. 
                         'keep' retains the nan values in the spectrum.
    :param removeSegmentSize: If nanTreatment = 'remove', this can leave a lot
                         of little fragments of spectrum in the output, mostly
                         in regions with strong telluric lines.  
                         Fragments that are shorter than removeSegmentSize 
                         pixels will be removed from the output. (Set this to
                         0 to keep all the spectrum)
    :param allowGapSize: If nanTreatment = 'remove', and small fragments
                         of spectrum are also being removed, there can be small
                         gaps of a few pixels in the spectrum. Setting
                         allowGapSize ignores gaps smaller than this value when
                         calculating the size of segments for removeSegmentSize.
    :param saveFitsHeader: Flag for whether to save the .fits header information
                         to a text file, named using the output filename
                         .out (True = save to file).
    :param writeSpecHeader: Flag for whether to write two lines of header in the
                         .s text file (True = include header)
    :rtype: Spectrum
    '''
    if isinstance(fname, list) or isinstance(fname, tuple):
        raise  ValueError('spirou_p() can only process one file'
                          ' at a time, not a list of files!')
    elif not isinstance(fname, str):
        raise  ValueError('in spirou_p(), file name must be a string!')
    if outname is None:
        outname = fname+'.s'
        if outname[-7:-2].lower() == '.fits': outname = outname[:-7]+'.s'
    print('converting', fname, 'to', outname)
    if (nanTreatment != 'remove' and nanTreatment != 'replace'
        and nanTreatment != 'keep'):
        raise ValueError('in spirou_p(), unrecognized nanTreatment '
                         'value: {:}'.format(nanTreatment))
    # read the fits file
    fitsSpec = fits.open(fname)
    header = fitsSpec[0].header
    header2 = fitsSpec[1].header
    
    #print('File format info')
    #print(fitsSpec.info())
    #print('Header information')
    #print(repr(header))

    target = header['OBJECT']
    dateUTC = header['DATE-OBS']
    timeUTC = header['UTIME']
    # get the 2D data arrays, with one row for each order
    polar2D    = fitsSpec[1].data
    polarErr2D = fitsSpec[2].data
    specI2D    = fitsSpec[3].data
    specIerr2D = fitsSpec[4].data
    null12D    = fitsSpec[5].data
    null22D    = fitsSpec[6].data
    wl2D       = fitsSpec[7].data
    blaze2D    = fitsSpec[8].data
    fitsSpec.close()

    # trim nans of the start and end of each order
    # (orders are often extracted past the end of the usable flux, and the data
    # are padded with nans so the rows for the orders all have the same length)
    spec = None
    nOrd, nPix = wl2D.shape
    for i in range(nOrd):
        # save this order into a Spectrum object (trimming nans off the edges)
        nuse = np.nonzero(np.isfinite(specI2D[i,:]))[0]
        # deal with orders that are entirely nan
        # this should be rare (e.g. telluric correction issues, very low S/N)
        if nuse.size < 1:
            if nanTreatment == 'keep':
                nuse = np.arange(nPix, dtype=int)
            else:
                continue
        
        ord = spf.Spectrum(wl2D[i,:][nuse[0]:nuse[-1]+1],
                           specI2D[i,:][nuse[0]:nuse[-1]+1],
                           polar2D[i,:][nuse[0]:nuse[-1]+1],
                           null12D[i,:][nuse[0]:nuse[-1]+1],
                           null22D[i,:][nuse[0]:nuse[-1]+1],
                           polarErr2D[i,:][nuse[0]:nuse[-1]+1])
        ordBlaze = blaze2D[i,:][nuse[0]:nuse[-1]+1]

        # APERO p.fits files have a continuum normalized I/Ic,
        # but *not* continuum normalized polarization spectra (like V/I, not V/Ic).
        # So convert the polarization spectrum to be continuum normalized,
        # like the intensity spectrum.
        ord.specV = ord.specV*ord.specI
        ord.specN1 = ord.specN1*ord.specI
        ord.specN2 = ord.specN2*ord.specI
        ord.specSig = ord.specSig*ord.specI
        
        # get remaining nan values in the order
        nuse = np.isfinite(ord.wl) & np.isfinite(ord.specI) & np.isfinite(ord.specV)
        nan = np.logical_not(nuse)
        
        if nanTreatment == 'replace':
            # replace NaNs with zeros (for the normalization routine)
            ord.specI[nan]   = 0.0
            ord.specV[nan]   = 0.0
            ord.specN1[nan]  = 0.0
            ord.specN2[nan]  = 0.0
            ord.specSig[nan] = 100.0

            # make up some wavelengths where they don't exist
            # (this version can be a bit slow, assigning one wl value at a time)
            indfix = np.nonzero(np.isnan(ord.wl))[0]
            indok = np.nonzero(np.isfinite(ord.wl))[0]
            for j in indfix:
                lastOk = indok[indok < j][-1]
                nextOk = indok[indok > j][0]
                wlStep = (ord.wl[nextOk] - ord.wl[lastOk])/float(nextOk - lastOk)
                wlGuess = ord.wl[lastOk] + wlStep*(j - lastOk)
                ord.wl[j] = wlGuess

        elif nanTreatment == 'remove':
            # cut the NaNs, and also small fragments of spectrum between nans
            _flag_spec_fragments(nuse, removeSegmentSize, allowGapSize)
            ord = ord[nuse]

        # add this order to the full spectrum
        if spec is None:
            spec = ord
        else:
            spec = spec.concatenate(ord)

    # add a header and save the result
    spec.header = '***Reduced spectrum {:} {:} {:}\n'.format(target, dateUTC, timeUTC)
    spec.save(outname, saveHeader=writeSpecHeader)
    
    if saveFitsHeader:
        outname2 = outname+'.out'
        if outname2[-6:-4] == '.s': outname2 = outname2[:-6]+'.out'
        fout = open(outname2, 'w')
        fout.write('# Main FITS header\n')
        fout.write(repr(header)+'\n')
        fout.write('# Additional spectrum FITS header\n')
        # avoid repeating redundant header information
        for card in header2.cards:
            if (card.keyword not in header) or (card.keyword == 'HISTORY'
                                            or card.keyword == 'COMMENT'):
                fout.write(str(card)+'\n')
        fout.close()
    return spec


def spirou_e(fname, outname=None, nanTreatment='replace',
             removeSegmentSize=100, allowGapSize=3, bervCorr=True,
             saveFitsHeader=True, writeSpecHeader=False):
    '''
    Function to read in a SPIRou e.fits file containing an intensity spectrum,
    save it to a .s file, and return a Spectrum object

    This function estimates uncertainties on the spectrum simply as the 
    square root of the flux, then applies a blaze correction to the spectrum.
    
    The e.fits files usually contain nan values for pixels where the spectrum
    extraction failed (most commonly near order edges where the flux is low).
    This function provides options: to remove those pixels (and optionally
    remove some of the small fragments of spectrum that get left behind), 
    to replace those pixels with placeholder values, or to keep the
    nan values. (Keeping the nan values will cause problems for other
    SpecpolFlow functions!)

    :param fname: name of the .fits file to read
    :param outname: name of the file to save the spectrum to.  If not provided
                    the function will default to using the input fname and
                    replacing .fits with .s
    :param nanTreatment: Flag for how to deal with nan values in the input file,
                         can be: 'remove', 'replace', 'keep'. 
                         'remove' deletes pixels with any nans.
                         'replace' replaces the nan values with 0, and sets the
                         uncertainties to 100. 
                         'keep' retains the nan values in the spectrum.
    :param removeSegmentSize: If nanTreatment='remove', this can leave
                         little fragments of spectrum in the output.
                         Fragments that are shorter than removeSegmentSize
                         pixels will be removed from the output. (Set this to
                         0 to keep all the spectrum)
    :param allowGapSize: If nanTreatment = 'remove', and small fragments
                         of spectrum are also being removed, there can be small
                         gaps of a few pixels in the spectrum. Setting
                         allowGapSize ignores gaps smaller than this value when
                         calculating the size of segments for removeSegmentSize.
    :param bervCorr:     Flag for whether to apply a barycentric radial velocity
                         correction to the spectrum (using the BERV from the
                         fits header, True = apply correction).  Without this
                         correction the observation is in the observer's rest
                         frame, rather than the solar system barycentric rest
                         frame.
    :param saveFitsHeader: Flag for whether to save the .fits header information
                         to a text file, named using the output filename
                         .out (True = save to file).
    :param writeSpecHeader: Flag for whether to write two lines of header in the
                         .s text file (True = include header)
    :rtype: Spectrum
    '''
    if isinstance(fname, list) or isinstance(fname, tuple):
        raise  ValueError('spirou_e() can only process one file'
                          ' at a time, not a list of files!')
    elif not isinstance(fname, str):
        raise  ValueError('in spirou_e(), file name must be a string!')
    if outname is None:
        outname = fname+'.s'
        if outname[-7:-2].lower() == '.fits': outname = outname[:-7]+'.s'
    print('converting', fname, 'to', outname)
    if (nanTreatment != 'remove' and nanTreatment != 'replace'
        and nanTreatment != 'keep'):
        raise ValueError('in spirou_e(), unrecognized nanTreatment '
                         'value: {:}'.format(nanTreatment))
    # read the fits file
    fitsSpec = fits.open(fname)
    header = fitsSpec[0].header
    header2 = fitsSpec[1].header
    header3 = fitsSpec[2].header
    
    #print('File format info')
    #print(fitsSpec.info())
    #print('Header information')
    #print(repr(header))

    target = header['OBJECT']
    dateUTC = header['DATE-OBS']
    timeUTC = header['UTIME']
    # get the 2D data arrays, with one row for each order
    specI2D = fitsSpec[1].data
    wl2D    = fitsSpec[5].data
    blaze2D = fitsSpec[9].data
    fitsSpec.close()
    
    # trim nans of the start and end of each order
    # (orders are often extracted past the end of the usable flux, and the data
    # are padded with nans so the rows for the orders all have the same length)
    spec = None
    nOrd, nPix = wl2D.shape
    for i in range(nOrd):
        # save this order into a Spectrum object (trimming nans off the edges)
        nuse = np.nonzero(np.isfinite(specI2D[i,:]))[0]
        # deal with orders that are entirely nan
        # this should be rare (e.g. telluric correction issues, very low S/N)
        if nuse.size < 1:
            if nanTreatment == 'keep' or nanTreatment == 'replace':
                nuse = np.nonzero(np.isfinite(wl2D[i,:]))[0]
                if nuse.size < 1: nuse = np.arange(nPix, dtype=int)
            else:
                continue

        ord = spf.Spectrum(wl2D[i,:][nuse[0]:nuse[-1]+1],
                           specI2D[i,:][nuse[0]:nuse[-1]+1],
                           np.zeros(nuse[-1] - nuse[0] + 1),
                           np.zeros(nuse[-1] - nuse[0] + 1),
                           np.zeros(nuse[-1] - nuse[0] + 1),
                           100.*np.ones(nuse[-1] - nuse[0] + 1))
        ordBlaze = blaze2D[i,:][nuse[0]:nuse[-1]+1]

        # get remaining nan values in the order
        nuse = np.isfinite(ord.wl) & np.isfinite(ord.specI) & np.isfinite(ordBlaze)
        nan = np.logical_not(nuse)
        
        # estimate uncertainties
        iok = (nuse) & (ord.specI > 0.0)
        ord.specSig[iok] = np.sqrt(ord.specI[iok])
        
        # apply blaze correction
        ord.specI = ord.specI/ordBlaze
        ord.specSig = ord.specSig/ordBlaze
        # and set uncertainty estimates for pixels with negative flux to be large
        if nanTreatment == 'keep':
            ord.specSig[np.logical_not(iok)] = np.nan
        else:
            if np.count_nonzero(iok) > 10:
                #ord.specSig[np.logical_not(iok)] = np.percentile(ord.specSig[iok], 99.9)
                ord.specSig[np.logical_not(iok)] = 0.2*np.percentile(ord.specI[iok], 99.0)

        if nanTreatment == 'replace':
            # replace NaNs with zeros (for the normalization routine)
            ord.specI[nan] = 0.0
            ord.specSig[nan] = 100.0
        
        elif nanTreatment == 'remove':
            # cut the NaNs, and also small fragments of spectrum between nans
            _flag_spec_fragments(nuse, removeSegmentSize, allowGapSize)
            ord = ord[nuse]

        # add this order to the full spectrum
        if spec is None:
            spec = ord
        else:
            spec = spec.concatenate(ord)

    # The e.fits spectra default to being in the observer's rest frame rather
    # than the solar system barycentric rest frame. Optionally, correct for this
    # velocity shift (p.fits files seem to already have this correction).
    if bervCorr:
        try:
            berv = header2['BERV']
            bjd = header2['BJD']
            print("applying barycentric radial velocity correction "
                  "{:8.4f} km/s, Barycentric Julian date {:.6f}".format(berv, bjd))
        except:
            berv = 0.0
            print("Warning: Failed to find barycentric radial velocity, "
                  "no correction applied!")
        c = 299792.458
        spec.wl = spec.wl + spec.wl*berv/c

    # add a header and save the result
    spec.header = '***Reduced spectrum {:} {:} {:}\n'.format(target, dateUTC, timeUTC)
    spec.save(outname, saveHeader=writeSpecHeader)

    if saveFitsHeader:
        outname2 = outname+'.out'
        if outname2[-6:-4] == '.s': outname2 = outname2[:-6]+'.out'
        fout = open(outname2, 'w')
        fout.write('# Main FITS header\n')
        fout.write(repr(header)+'\n')
        fout.write('# Additional spectrum FITS header\n')
        # avoid repeating redundant header information
        for card in header2.cards:
            if (card.keyword not in header) or (card.keyword == 'HISTORY'
                                            or card.keyword == 'COMMENT'):
                fout.write(str(card)+'\n')
        fout.close()
    return spec


def spirou_t(fname, outname=None, nanTreatment='replace',
             removeSegmentSize=200, allowGapSize=3, bervCorr=True,
             saveFitsHeader=True, writeSpecHeader=False):
    '''
    Function to read in a SPIRou t.fits file containing a telluric corrected
    intensity spectrum, save it to a .s file, and return a Spectrum object

    This function estimates uncertainties on the spectrum simply as the 
    square root of the flux, then applies a blaze correction to the spectrum.
    
    The t.fits files usually contain the telluric corrected spectrum
    and a copy of the telluric line spectrum removed from the original
    observation.  This function saves both to two different files
    (with the telluric spectrum in [outname].telluric).

    The t.fits files usually contain nan values for pixels where the telluric
    correction failed, or where the spectrum extraction failed.
    This function provides options: to remove those pixels (and optionally
    remove some of the small fragments of spectrum that get left behind),
    to replace those pixels with placeholder values, or to keep the
    nan values. (Keeping the nan values will cause problems for other
    SpecpolFlow functions!)
    
    :param fname: name of the .fits file to read
    :param outname: name of the file to save the spectrum to.  If not provided
                    the function will default to using the input fname and
                    replacing .fits with .s
    :param nanTreatment: Flag for how to deal with nan values in the input file,
                         can be: 'remove', 'replace', 'keep'. 
                         nan values can be from bad pixels or from places where
                         the telluric correction was too uncertain.
                         'remove' deletes pixels with any nans.
                         'replace' replaces the nan values with 0, and sets the
                         uncertainties to 100. 
                         'keep' retains the nan values in the spectrum.
    :param removeSegmentSize: If nanTreatment = 'remove', this can leave a lot
                         of little fragments of spectrum in the output, mostly
                         in regions with strong telluric lines.  
                         Fragments that are shorter than removeSegmentSize 
                         pixels will be removed from the output. (Set this to
                         0 to keep all the spectrum)
    :param allowGapSize: If nanTreatment = 'remove', and small fragments
                         of spectrum are also being removed, there can be small
                         gaps of a few pixels in the spectrum. Setting
                         allowGapSize ignores gaps smaller than this value when
                         calculating the size of segments for removeSegmentSize.
    :param bervCorr:     Flag for whether to apply a barycentric radial velocity
                         correction to the spectrum (using the BERV from the
                         fits header, True = apply correction).  Without this
                         correction the observation is in the observer's rest
                         frame, rather than the solar system barycentric rest
                         frame.
    :param saveFitsHeader: Flag for whether to save the .fits header information
                         to a text file, named using the output filename
                         .out (True = save to file).
    :param writeSpecHeader: Flag for whether to write two lines of header in the
                         .s text file (True = include header)
    :rtype: Spectrum
    '''
    if isinstance(fname, list) or isinstance(fname, tuple):
        raise  ValueError('spirou_t() can only process one file'
                          ' at a time, not a list of files!')
    elif not isinstance(fname, str):
        raise  ValueError('in spirou_t(), file name must be a string!')
    if outname is None:
        outname = fname+'.s'
        if outname[-7:-2].lower() == '.fits': outname = outname[:-7]+'.s'
    print('converting', fname, 'to', outname, 'and', outname+'.telluric')
    if (nanTreatment != 'remove' and nanTreatment != 'replace'
        and nanTreatment != 'keep'):
        raise ValueError('in spirou_t(), unrecognized nanTreatment '
                         'value: {:}'.format(nanTreatment))
    # read the fits file
    fitsSpec = fits.open(fname)
    header = fitsSpec[0].header
    header2 = fitsSpec[1].header
    
    #print('File format info')
    #print(fitsSpec.info())
    #print('Header information')
    #print(repr(header))

    target = header['OBJECT']
    dateUTC = header['DATE-OBS']
    timeUTC = header['UTIME']
    # get the 2D data arrays, with one row for each order
    specI2D = fitsSpec[1].data
    wl2D    = fitsSpec[2].data
    blaze2D = fitsSpec[3].data
    tellurc2D = fitsSpec[4].data
    # This extra OH emission spectrum could be useful,
    # but I'm not sure how to interpret or use it.
    OHline2D = fitsSpec[5].data
    fitsSpec.close()
    
    # trim nans of the start and end of each order
    # (orders are often extracted past the end of the usable flux, and the data
    # are padded with nans so the rows for the orders all have the same length)
    spec = specTell = None
    nOrd, nPix = wl2D.shape
    for i in range(nOrd):
        # save this order into a Spectrum object (trimming nans off the edges)
        nuse = np.nonzero(np.isfinite(specI2D[i,:]))[0]
        # deal with orders that are entirely nan
        # this should be rare (e.g. telluric correction issues, very low S/N)
        if nuse.size < 1:
            if nanTreatment == 'keep' or nanTreatment == 'replace':
                nuse = np.nonzero(np.isfinite(wl2D[i,:]))[0]
                if nuse.size < 1: nuse = np.arange(nPix, dtype=int)
            else:
                continue

        ord = spf.Spectrum(wl2D[i,:][nuse[0]:nuse[-1]+1],
                           specI2D[i,:][nuse[0]:nuse[-1]+1],
                           np.zeros(nuse[-1] - nuse[0] + 1),
                           np.zeros(nuse[-1] - nuse[0] + 1),
                           np.zeros(nuse[-1] - nuse[0] + 1),
                           100.*np.ones(nuse[-1] - nuse[0] + 1))
        ordTell = spf.Spectrum(wl2D[i,:][nuse[0]:nuse[-1]+1],
                               tellurc2D[i,:][nuse[0]:nuse[-1]+1],
                               np.zeros(nuse[-1] - nuse[0] + 1),
                               np.zeros(nuse[-1] - nuse[0] + 1),
                               np.zeros(nuse[-1] - nuse[0] + 1),
                               np.zeros(nuse[-1] - nuse[0] + 1))
        ordBlaze = blaze2D[i,:][nuse[0]:nuse[-1]+1]

        # get remaining nan values in the order
        nuse = np.isfinite(ord.wl) & np.isfinite(ord.specI) & np.isfinite(ordBlaze)
        nan = np.logical_not(nuse)
        
        # estimate uncertainties
        iok = (nuse) & (ord.specI > 0.0)
        ord.specSig[iok] = np.sqrt(ord.specI[iok])
        
        # apply blaze correction
        ord.specI = ord.specI/ordBlaze
        ord.specSig = ord.specSig/ordBlaze
        # and set uncertainty estimates for pixels with negative flux to be large
        if nanTreatment == 'keep':
            ord.specSig[np.logical_not(iok)] = np.nan
        else:
            if np.count_nonzero(iok) > 10:
                #ord.specSig[np.logical_not(iok)] = np.percentile(ord.specSig[iok], 99.9)
                ord.specSig[np.logical_not(iok)] = 0.2*np.percentile(ord.specI[iok], 99.0)

        if nanTreatment == 'replace':
            # replace NaNs with zeros (for the normalization routine)
            ord.specI[nan] = 0.0
            ord.specSig[nan] = 100.0
            ordTell.specI[nan] = 0.0
        
        elif nanTreatment == 'remove':
            # cut the NaNs, and also small fragments of spectrum between nans
            _flag_spec_fragments(nuse, removeSegmentSize, allowGapSize)
            ord = ord[nuse]
            ordTell = ordTell[nuse]

        # add this order to the full spectrum
        if spec is None:
            spec = ord
            specTell = ordTell
        else:
            spec = spec.concatenate(ord)
            specTell = specTell.concatenate(ordTell)

    # The t.fits spectra default to being in the observer's rest frame rather
    # than the solar system barycentric rest frame. Optionally, correct for this
    # velocity shift (p.fits files seem to already have this correction).
    if bervCorr:
        try:
            berv = header2['BERV']
            bjd = header2['BJD']
            print("applying barycentric radial velocity correction "
                  "{:8.4f} km/s, Barycentric Julian date {:.6f}".format(berv, bjd))
        except:
            berv = 0.0
            print("Warning: Failed to find barycentric radial velocity, "
                  "no correction applied!")
        c = 299792.458
        spec.wl = spec.wl + spec.wl*berv/c
        specTell.wl = specTell.wl + specTell.wl*berv/c

    # add a header and save the result
    spec.header = '***telluric corrected spectrum {:} {:} {:}\n'.format(target, dateUTC, timeUTC)
    spec.save(outname, saveHeader=writeSpecHeader)
    specTell.header = '***telluric spectrum for {:} {:} {:}\n'.format(target, dateUTC, timeUTC)
    specTell.save(outname+'.telluric', saveHeader=writeSpecHeader)

    if saveFitsHeader:
        outname2 = outname+'.out'
        if outname2[-6:-4] == '.s': outname2 = outname2[:-6]+'.out'
        fout = open(outname2, 'w')
        fout.write('# Main FITS header\n')
        fout.write(repr(header)+'\n')
        fout.write('# Additional spectrum FITS header\n')
        # avoid repeating redundant header information
        for card in header2.cards:
            if (card.keyword not in header) or (card.keyword == 'HISTORY'
                                            or card.keyword == 'COMMENT'):
                fout.write(str(card)+'\n')
        fout.close()

    return spec


def _flag_spec_fragments(nuse, removeSegmentSize, allowGapSize):
    '''
    Helper function for converters.  Takes an array of bool flags for whether
    pixels are used.  Look for small fragments of spectrum that are used,
    and flag them to not be used.  Modifies the array of flags in place.
    '''
    istartOk = 0
    sizeOk = 0
    sizeNaN = 0
    # go through each pixel and flag short fragments of spectrum to remove.
    for i in range(nuse.size):
        if nuse[i] == False:
            sizeNaN += 1
            # If we reach the end of an ok region, but that region is too
            # small, don't use it.
            if sizeOk > 0 and sizeOk < removeSegmentSize:
                # But wait until we have found a large enough gap
                # before removing the small ok region. (if the NaN gap
                # is too small, keep adding to the current ok region)
                if sizeNaN > allowGapSize:
                    nuse[istartOk:i] = False
            # If we have found too large a NaN gap, reset the ok region
            if sizeNaN > allowGapSize and sizeOk > 0:
                sizeOk = 0
        else:
            if sizeOk == 0: istartOk = i
            sizeOk += 1
            sizeNaN = 0
    # If we haven't closed the last region to remove
    if sizeNaN > 0 and sizeOk > 0:
        nuse[istartOk:i+1] = False
    # If we have a fragment as the last region
    elif sizeOk > 0 and  sizeOk < removeSegmentSize:
        nuse[istartOk:i+1] = False        
    return
