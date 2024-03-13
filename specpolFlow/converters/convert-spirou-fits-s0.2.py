#!/usr/bin/python
# Convert SPIRou reduced spectra in fits files to a text column
# format similar to LibreESPRIT.

import numpy as np
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(description='Convert SPIRou data reduction pipeline .fits files into plain text .s files')
parser.add_argument('file')
parser.add_argument("-d", "--head", action="store_true", help='Write two lines of headder in the LibreESPRIT format.')

args = parser.parse_args()
fname = args.file
headder = args.head

print('reading {:}'.format(fname))
pfile = fits.open(fname)

#print 'headder'
#print pfile[0].header
#print pol.header
nextend = pfile[0].header['NEXTEND']
target = pfile[0].header['OBJECT']
dateUTC = pfile[0].header['DATE-OBS']
timeUTC = pfile[0].header['UTIME']
print('observed {:} on {:} at {:}'.format(target, dateUTC, timeUTC))

if fname[-5:] == '.fits':
    outName = fname[:-5]+'.s'
    print('writing {:}'.format(outName))
else:
    outName = fname+'.s'
outFile = open(outName,'w')
#outFile.write('***spec. of {:} on {:} at {:} UTC\n'.format(target, dateUTC, timeUTC))


#For polarization spectra "*p.fits", expect 8 extensions in the fits file
if nextend == 8:
    if fname[-6:] != 'p.fits':
        print('Warning: unexpected file extension {:}, treating file as a polarization spectrum'.format(fname[-6:]))

    pol = pfile[1]
    polerr = pfile[2]
    specI = pfile[3]
    specIerr = pfile[4]
    null1 = pfile[5]
    null2 = pfile[6]
    wl = pfile[7]
    blaze = pfile[8]
    
    #calculate output file length
    flen = 0
    for i in range(len(pol.data)):
        jblaze = np.where(np.isfinite(blaze.data[i]))[0]
        flen += len(jblaze)
    if headder:
        outFile.write('***Reduced spectrum {:} {:} {:}\n'.format(target, dateUTC, timeUTC))
        outFile.write('{:} 6\n'.format(flen))

    print('wl, specI, pol, null1, null2, polerr, specIerr (blaze corrected)')
    for i in range(len(pol.data)):
        #Assume the polarization parameter is X/I (not X or X/Ic)
        #Apply blaze correction
        specIb = specI.data[i]/blaze.data[i]
        specIerrb = specIerr.data[i]/blaze.data[i]
        #The blaze function should have divided out in the X/I calculation
        #Convert to just X rather than X/I,
        #so that the normalization code can divide by Ic to get X/Ic
        #(this should make the two uncertainty values nearly identical)
        polb = pol.data[i]*specIb
        polerrb = np.abs(polerr.data[i]*specIb)
        null1b = null1.data[i]*specIb
        null2b = null2.data[i]*specIb
        #print '{:} {:e} {:e} {:e} {:e} {:e}'.format(i, specI.data[i][100], specIerr.data[i][100], pol.data[i][100], polerr.data[i][100], blaze.data[i][100] )
        #Check for bad values (nan) to exclude/treat
        jblaze = np.where(np.isfinite(blaze.data[i]))[0]
        jOk = np.isfinite(polb) & np.isfinite(specIb) & np.isfinite(polerrb) & np.isfinite(specIerrb)
        for j in jblaze:
            if jOk[j]:
                outFile.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(wl.data[i][j], specIb[j], polb[j], null1b[j], null2b[j], polerrb[j], specIerrb[j]))
            else:
                #replace nan's with zeros (rather than just discarding,
                # helps for finding order edges)
                outFile.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(wl.data[i][j], 0.0, 0.0, 0.0, 0.0, 100., 100.))


#For merged 1D intensity "*s.fits", expect 2 extensions in the fits file, use constant in velocity pixels
if nextend == 2:
    if fname[-6:] != 's.fits':
        print('Warning: unexpected file extension {:}, treating file as a merged intensity spectrum'.format(fname[-6:]))

    wl = pfile[2].data.field(0)
    specI = pfile[2].data.field(1)
    specIerr = pfile[2].data.field(2)
    specIcorr = pfile[2].data.field(9)
    specIcorrerr = pfile[2].data.field(10)

    jOk = np.isfinite(specI) & np.isfinite(specIerr) & np.isfinite(specIcorr) & np.isfinite(specIcorrerr)

    if headder:
        outFile.write('***Reduced spectrum {:} {:} {:}\n'.format(target, dateUTC, timeUTC))
        outFile.write('{:} 4\n'.format(len(wl)))

    print('wl, specI, specIerr, specI_tcorr, specI_tcorrerr (merged orders)')
    for j in range(len(wl)):
        if jOk[j]:
            outFile.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(wl[j], specI[j], specIerr[j], specIcorr[j], specIcorrerr[j]))
        else:
            outFile.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(wl[j], 0.0, 100., 0.0, 100.))
    
    
#For telluric corrected intensity spectra "*t.fits", expect 4 extensions in the fits file
if nextend == 4:
    if fname[-6:] != 't.fits':
        print('Warning: unexpected file extension {:}, treating file as a telluric corrected intensity spectrum'.format(fname[-6:]))

    specI = pfile[1]
    wl = pfile[2]
    blaze = pfile[3]
    telluric = pfile[4]

    #calculate output file length
    flen = 0
    for i in range(len(blaze.data)):
        jblaze = np.where(np.isfinite(blaze.data[i]))[0]
        flen += len(jblaze)
    if headder:
        outFile.write('***Reduced spectrum {:} {:} {:}\n'.format(target, dateUTC, timeUTC))
        outFile.write('{:} 2\n'.format(flen))

    print('wl, specI_tcorr, telluric (blaze corrected)')
    for i in range(len(specI.data)):
        specIb = specI.data[i]/blaze.data[i]
        jblaze = np.where(np.isfinite(blaze.data[i]))[0]
        jOk = np.isfinite(specIb) & np.isfinite(telluric.data[i])
        for j in jblaze:
            if jOk[j]:
                outFile.write('{:10.4f} {:11.4e} {:11.4e}\n'.format(wl.data[i][j], specIb[j], telluric.data[i][j]))
            else:
                #replace nan's with zeros (rather than just discarding,
                # helps for finding order edges)
                outFile.write('{:10.4f} {:11.4e} {:11.4e}\n'.format(wl.data[i][j], 0.0, 0.0))
    #print telluric.header


#For intensity spectra of each beam "*e.fits", expect 12 extensions in the fits file
if nextend == 12:
    if fname[-6:] != 'e.fits':
        print('Warning: unexpected file extension {:}, treating file as an intensity spectrum with beams'.format(fname[-6:]))

    specI = pfile[1]
    wl = pfile[5]
    blaze = pfile[9]
    specIA = pfile[2]
    wlA = pfile[6]
    blazeA = pfile[10]
    specIB = pfile[3]
    wlB = pfile[7]
    blazeB = pfile[11]

    #calculate output file length
    flen = 0
    for i in range(len(blaze.data)):
        jblaze = np.where(np.isfinite(blaze.data[i]))[0]
        flen += len(jblaze)
    if headder:
        outFile.write('***Reduced spectrum {:} {:} {:}\n'.format(target, dateUTC, timeUTC))
        outFile.write('{:} 1\n'.format(flen))
    
    print('wl, specI (blaze corrected)')
    for i in range(len(specI.data)):
        specIb = specI.data[i]/blaze.data[i]
        jblaze = np.where(np.isfinite(blaze.data[i]))[0]
        jOk = np.isfinite(specIb)
        for j in jblaze:
            if jOk[j]:
                outFile.write('{:10.4f} {:11.4e}\n'.format(wl.data[i][j], specIb[j]))
            else:
                #replace nan's with zeros (rather than just discarding,
                # helps for finding order edges)
                outFile.write('{:10.4f} {:11.4e}\n'.format(wl.data[i][j], 0.0))

outFile.close()
    
#import matplotlib.pyplot as plt
#for i in range(len(wl.data)):
#    plt.plot(wl.data[i], specI.data[i])
#    plt.plot(wl.data[i], blaze.data[i])
#    plt.plot(wl.data[i], telluric.data[i])
#plt.show()

