#!/usr/bin/python3
#
# Generate an LSD line mask from a VALD line list.

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

def make_mask_cli():
    """
    Main function for running make_mask as a terminal program.

    Takes no arguments, but uses command line parameters instead.
    """
    #Take input and output file names as command line arguments,
    #with some additional optional control parameters.
    parser = argparse.ArgumentParser(description='Generate an LSD line mask from a VALD line list. When Lande factors are not in the line list they will be estimated if possible.')
    parser.add_argument("lineList", help="Line list from VALD. This file should be an 'extract stellar' line list in the 'long' format.")
    parser.add_argument("mask", nargs='?', default='mask.dat', help='Save the line mask to this file.')
    parser.add_argument("-d", "--depth", default='0', help='Optional number, only include lines deeper than this depth cutoff in the mask (defaults to 0, all lines included)')
    parser.add_argument("-w1", default='0', help='Optional number, starting wavelength for lines included in the mask (in nm)')
    parser.add_argument("-w2", default='1e10', help='Optional number, ending wavelength for lines included in the mask (in nm)')
    parser.add_argument("-g1", default='-1e10', help='Optional number, lower cutoff in effective Lande factor for lines included in the mask')
    parser.add_argument("-g2", default='1e10', help='Optional number, upper cutoff in effective Lande factor')
    parser.add_argument("-e", "--elements", default='', help="Optional, elements to include in the mask, as list symbols. For multiple elements enclose in quotes, e.g. 'C O Fe Ni'. By default all elements are used.")
    parser.add_argument("-ex", "--excludeEl", default='', help="Optional, elements to exclude from the mask, as a list symbols. For multiple elements enclose in quotes, e.g. 'C O Fe Ni'. By default all elements are used.")
    args = parser.parse_args()
    
    llistName = args.lineList
    maskName = args.mask
    depthCutoff = float(args.depth)
    wl1 = float(args.w1)
    wl2 = float(args.w2)
    lande1 = float(args.g1)
    lande2 = float(args.g2)
    elementsUsed = args.elements.split()
    elementsExclude = args.excludeEl.split()
    
    #Run the line mask generating code
    print('reading line list from', llistName)
    mask = pol.make_mask(llistName, outMaskName = maskName, 
                         depthCutoff = depthCutoff, wlStart = wl1, wlEnd = wl2,
                         landeStart = lande1, landeEnd = lande2,
                         elementsUsed = elementsUsed,
                         elementsExclude = elementsExclude,
                         includeNoLande = False)
    print('mask saved to', maskName)
    return

#For running this script as a terminal program
if __name__ == "__main__":
    make_mask_cli()
