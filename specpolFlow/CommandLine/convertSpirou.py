#!/usr/bin/python3
#
# Convert SPIRou .fits spectrum files from the APERO format to a text .s format.

import argparse
import matplotlib.pyplot as plt
try:
    import specpolFlow as pol
except ModuleNotFoundError:
    #If specpolFlow is not installed or in the python path
    #try guessing the location and adding it temporarily
    import sys, os
    loc0 = os.path.dirname(__file__) 
    sys.path.insert(0, os.path.join(loc0, '..', '..'))
    import specpolFlow as pol


def spirou_cli():
    """
    Main function for running the SPIRou converters as a terminal program.
    
    Takes no arguments, but uses command line parameters instead.
    """
    import argparse
    parser = argparse.ArgumentParser(description='Convert FITS file spectra '
                'from the SPIRou APERO format to text .s files. Output files '
                'as [filename].s and [filename].out for header information. \n'
                'Supports e.fits (extracted spectrum in separate orders), '
                't.fits (telluric corrected spectrum in separate orders), '
                'and p.fits (polarized spectrum) files. \n'
                'The FITS files contain nan values where the spectrum '
                'extraction or telluric correction failed.  There are options '
                'to either remove or replace those nan values.')
    parser.add_argument("observation", nargs='*',
                        help='file or list of FITS files to process.')
    parser.add_argument('-t', '--type', default=None,
                        help='optional, specify the type of the fits file to '
                        'process.  Can be e, t, or p.  If not specified the '
                        'type will be inferred from the file name '
                        '(if the name ends in e.fits, t.fits or p.fits)')
    parser.add_argument('-n', '--nan', metavar='TREATMENT', default='replace',
                        help='specify how to treat nan values, can be '
                        '"replace" (replace nan with 0 and uncertainties '
                        'with 100), "remove" (delete pixels with nan and '
                        'some small fragments), or "keep".')
    parser.add_argument('-s', '--removeSegment', metavar='NUMPIX',
                        type=int, default=200,
                        help='optional, only used with "-n remove". Removing '
                        'nan values can leave small fragments of spectrum, '
                        'which are hard to normalize or analyze. '
                        'Fragments smaller than this number of pixels will '
                        'also be removed. Set this to 0 to include all pixels '
                        '(except nans).')
    parser.add_argument('-g', '--allowGap', metavar='NUMPIX',
                        type=int, default=3,
                        help='optional, only used with "-n remove". Removing '
                        'nan values can leave gaps of a few pixels in the '
                        'spectrum.  Gaps smaller than this will be ignored '
                        'when calculating the size of fragments to remove '
                        'from the spectrum.')
    parser.add_argument('-o', '--outputHeader', action='store_true',
                        help='save the header information from the fits file '
                        'to a text file (with a name ending in .out).')
    args = parser.parse_args()

    flist = args.observation
    ftype = args.type
    if isinstance(ftype, str): ftype = ftype.lower()
    nanTreatment = args.nan.lower()
    removeSegmentSize = args.removeSegment
    allowGapSize = args.allowGap
    saveFitsHeader = args.outputHeader

    if flist == []:
        print('No files given')
    else:
        #Run the conversion function
        pol.converters.spirou(flist, ftype=ftype, nanTreatment=nanTreatment,
                              removeSegmentSize=removeSegmentSize,
                              allowGapSize=allowGapSize,
                              saveFitsHeader=saveFitsHeader)
    return

################

#For running this Python script as a terminal program
if __name__ == "__main__":
    spirou_cli()
