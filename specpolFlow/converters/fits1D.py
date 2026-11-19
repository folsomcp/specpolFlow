"""
Reads a set of FITS spectrum files and convert to s format.
"""

import numpy as np
import astropy.io.fits as fits
from .. import obsSpec as spf

def fits1D(flist,flistout=None,writeSpecHeader=False,isEchelle=False,wlUnit=None):
    """
    Convert a list of 1D .fits files for intensity only (no polarization) spectra, into text .s files. 
    No radial velocity or heliocentric velocity correction is applied by this function. The error is set to be zero, the user must manual change this.
    The headers of these files MUST contain the following keywords:

    * CRPIX1 - the pixel where the reference wavelength (CRVAL1) is given
    * CRVAL1 - the wavelength value of the reference pixel
    * CDELT1 - the dispersion in units of wavelength/pixel
    * NAXIS1 - the length of the data array
    
    Additionally if writeSpecHeader=True, OBJNAME (the object name) 
    and DATE-OBS (the UTC date and time of the observation) must also be specified.

    If wlUnit=None, CUNIT1 must be specified and the value must be (`Angstrom` or `nm`)

    This function generates a single file in a .s format:  
    
    * The UPENA spectrum (.s) 

    The content of the fits header is also saved in a .out text file.
    If flistout=None (default), the files are written at the same path as each entry of flist
    and '.s' is appended to the filename root. If isEchelle is True, then the output file is the same as the first entry of flist. 
    If flistout is a list of paths ending with a rootname, each file will be
    saved at that path with that rootname, and '.s' is appended to
    the filename root. If isEchelle is True then only the first entry of flistout will be used.

    :param flist: (list of strings) a list of fits filenames
    :param flistout: (list of strings) optional, output file rootname
    :param writeSpecHeader: Flag for whether to write the reduced fits header in the
                         .s text file (True = include header)
    :param isEchelle: Flag for whether each file in flist is 1 order of an echelle spectra. The output will then be a single .s file with all orders.
    :param wlUnit: Must be None, 'Angstrom', or 'nm'. Specifies the wavelength unit of the input files. If None is specified then each fits file must contain the CUNIT1 keyword
    :return: List of spectrum objects. If isEchelle is True then a single spectrum object is returned
    """

    if isinstance(flist, str): flist = [flist,]
    #Some basic error checking
    if not isinstance(flist, list):
        raise ValueError('in fits1D(), the flist argument must be a Python '
                        'list of input file names (or a single file name)')
    if not (flistout is None):
        if not isinstance(flist, list):
            raise ValueError('in fits1D(), the flistout argument must be a '
                            'Python list of output file names')
    
    #set default output file name
    if flistout is None:
        flistout=flist


    wl_arr=[]
    I_arr=[]
    Ierr_arr=[]

    target_arr=[]
    dateUTC_arr=[]

    for i, fname in enumerate(flist):

        hdul = fits.open(fname) #Opens the fits file

        ## checks for required keywords
        try:
            hdul[0].header['CRPIX1']
        except:
            raise ValueError('in fits1D(), the file {} is missing the `CRPIX1` entry from the header'.format(fname))
        try:
            hdul[0].header['CRVAL1']
        except:
            raise ValueError('in fits1D(), the file {} is missing the `CRVAL1` entry from the header'.format(fname))
        try:
            hdul[0].header['CDELT1']
        except:
            raise ValueError('in fits1D(), the file {} is missing the `CDELT1` entry from the header'.format(fname))
        try:
            hdul[0].header['NAXIS1']
        except:
            raise ValueError('in fits1D(), the file {} is missing the `NAXIS1` entry from the header'.format(fname))

        ##checks for required wavelength keyword if wlUnit is not specified
        if wlUnit==None:
            try:
                hdul[0].header['CUNIT1']
            except:
                raise ValueError('in fits1D(), `wlUnit` is `None` but the file {} is missing the `CUNIT1` entry from the header'.format(fname))
        if (wlUnit!=None) and (wlUnit!='Angstrom') and (wlUnit!='nm'):
            raise ValueError('in fits1D(), `wlUnit` must be `None`, `Anstrom`, or `nm`')

        ref= hdul[0].header['CRPIX1'] #location of the reference pixel for the wavelength array
        p1val=hdul[0].header['CRVAL1'] #wavelength value of the reference pixel
        dp1=hdul[0].header['CDELT1'] #step size of the wavelength array
        naxis1=hdul[0].header['NAXIS1'] #length of the wavelength array

        #Setting up the wavelength array
        wl=np.arange(p1val+(1-ref)*dp1,p1val+(1-ref)*dp1+dp1*naxis1-dp1/2,dp1) 

        #Grabbing the spectum data
        I=hdul[0].data

        #Since the spectrum could be normalized or unnormalized we will set the Ierr as 10% of I. The user can change this after the fact.
        Ierr=I*0.0

        #convert the wavelength array to nanometers
        if wlUnit==None:
            if hdul[0].header['CUNIT1']=='Angstrom':
                wl_arr.append(wl/10)
            
            if hdul[0].header['CUNIT1']=='nm':
                wl_arr.append(wl)

            #error catching
            if hdul[0].header['CUNIT1'] != 'Angstrom' and hdul[0].header['CUNIT1'] != 'nm':
                raise ValueError('in fits1D(), the file {} must have `Angstrom` or `nm` specified in `CUNIT1` entry of the header'.format(fname)  )  

        if wlUnit=='Angstrom':
            wl_arr.append(wl/10)
        if wlUnit=='nm':
            wl_arr.append(wl)
        
        
        I_arr.append(I)
        Ierr_arr.append(Ierr)

        # write the header if specified
        flag=0
        if writeSpecHeader==True:
            flag=1
            try:
                hdul[0].header['OBJNAME']
            except:
                print('WARNING: in fits1D(), `writeSpecHeader` is `True` but the file {} is missing the `OBJNAME` entry from the header'.format(fname))
                flag=0
            try:
                hdul[0].header['DATE-OBS']
            except:
                print('WARNING: in fits1D(), `writeSpecHeader` is `True` but the file {} is missing the `DATE-OBS` entry from the header'.format(fname))
                flag=0

        if flag==1:
            target_arr.append(hdul[0].header['OBJNAME'])
            dateUTC_arr.append(hdul[0].header['DATE-OBS'])
        if flag==0: 
            target_arr.append('N/A')
            dateUTC_arr.append('N/A')

    # if flist is a list of echelle orders, combine them into a single output spectrum
    if isEchelle==True:


        zeros = np.zeros_like(np.concatenate(I_arr))
        spec_u=spf.Spectrum(np.concatenate(wl_arr), np.concatenate(I_arr), zeros, zeros, zeros, np.concatenate(Ierr_arr))
        spec_u.header = '***Reduced spectrum {:} {:}\n'.format(
                target_arr[0], dateUTC_arr[0])
        
        spec_u.save(flistout[0]+'.s', saveHeader=writeSpecHeader)

        outHeader = open(flistout[0]+'.out','w')
        outHeader.write(repr(hdul[0].header))
        outHeader.close()

        return spec_u

    # if flist is a list of spectra, output a list of spectra
    if isEchelle==False:
        specs=[]
        for i, fname in enumerate(flist):
            zeros = np.zeros_like(I_arr[i])
            spec_u=spf.Spectrum(wl_arr[i], I_arr[i], zeros, zeros, zeros, Ierr_arr[i])
            spec_u.header = '***Reduced spectrum {:} {:}\n'.format(
                    target_arr[i], dateUTC_arr[i])
            spec_u.save(flistout[i]+'.s', saveHeader=writeSpecHeader)

            outHeader = open(flistout[i]+'.out','w')
            outHeader.write(repr(hdul[0].header))
            outHeader.close()

            specs.append(spec_u)

        return specs





