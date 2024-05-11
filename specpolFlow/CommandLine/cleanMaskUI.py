#!/usr/bin/python3
#
# Interactive tool for cleaning poor lines from an LSD line mask
# also can fit line depths in the mask.

import argparse

try:
    import specpolFlow as pol
except ModuleNotFoundError:
    #If specpolFlow is not installed or in the python path
    #try guessing the location and adding it temporarily
    import sys, os
    loc0 = os.path.dirname(__file__)
    sys.path.insert(0, os.path.join(loc0, '..', '..'))
    import specpolFlow as pol

def clean_mask_cli():
    """Main function for running cleanMaskUI as a terminal program.
    """
    parser = argparse.ArgumentParser(description='Interactively remove poor lines from an LSD line mask, and optionally modify line depths. This plots the mask, a comparison observation, and the model LSD spectrum (the convolution of the mask and LSD profile). For the model LSD spectrum calculation, some more useful parameters can be set within the graphical UI, and there is an option for displaying the resulting profile.')
    parser.add_argument('mask', help='The line mask to clean')
    parser.add_argument('observation', help='An observed spectrum for comparison')
    parser.add_argument('outName', nargs='?', default=None, help='Optional name for the output cleaned line mask, defaults to [input_file_name].clean')
    parser.add_argument('-e', '--exclude', default='excludeRanges.dat', help='Optional, a file of wavelength ranges to exclude from the line mask.  If the file does not exist default initial values will be used.')
    parser.add_argument('-b', '--batch', action='store_true', help='Run in batch mode, just apply the exclude regions from file, skip the interactive UI and LSD calculations.')
    args = parser.parse_args()

    maskName = args.mask
    obsName = args.observation
    outMaskName = args.outName
    excludeFileName = args.exclude
    batchMode = args.batch
    
    cleanMask, excludeRanges = pol.cleanMaskUI(maskName, obsName,
                                               outMaskName=outMaskName,
                                               inExcludeName=excludeFileName,
                                               outExcludeName=excludeFileName,
                                               batchMode=batchMode)
    return

#For running this Python script as a terminal program
if __name__ == "__main__":
    clean_mask_cli()
