#!/usr/bin/python3
#
# Convert ESPaDOnS .fits files from the UPENA format to a text .s format.

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


def espadons_cli():
    """
    Main function for running the ESPaDOnS converter as a terminal program.
    
    Takes no arguments, but uses command line parameters instead.
    """
    import argparse
    parser = argparse.ArgumentParser(description='Convert FITS file spectra '
            'from the ESPaDOnS CADC archive format to text .s files. Output '
            'files as [filename]n.s for the pipeline normalized spectra, '
            '[filename]u.s for unnormalized spectra, and [filename].out for '
            'header information. Supports i.fits (intensity spectrum) '
            'and p.fits (polarization & intensity spectrum) files.')
    parser.add_argument("observation", nargs='*',
                        help='a list of FITS files to process.')
    args = parser.parse_args()

    flist=args.observation

    if flist == []:
        print('No files given')
    else:
        #Run the conversion funciton
        pol.converters.espadons(flist)
    return

################

#For running this scirpt as a terminal program
if __name__ == "__main__":
    espadons_cli()
