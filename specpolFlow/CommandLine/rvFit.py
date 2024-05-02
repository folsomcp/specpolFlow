#!/usr/bin/python3
#
# A simple tool for measuring radial velocities from LSD profiles.

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

def rv_fit_cli():
    """Main function for running LSD.fit_gaussian_rv() as a terminal program

    Takes no arguments, but uses command line parameters instead.
    """
    #Take input file names as command line arguments,
    #with some additional optional control parameters.
    parser = argparse.ArgumentParser(description='Fit a Gaussian to an LSD line profile, to determine radial velocity.  Prints the radial velocity in the same velocity units as the LSD profile.')
    parser.add_argument("fileList", nargs='*',
                        help='LSD profile files, to determine RVs for.  Can be more than one file.')
    parser.add_argument("-v", "--velRange", nargs=2, type=float, metavar=('VEL1', 'VEL2'),
                        help='Optional starting and ending velocity for the range used in the fit')
    parser.add_argument("-p", "--plotFit", action='store_true',
                        help='Optional plot the fit to the LSD profile')
    args = parser.parse_args()
    #Process the command line parameters
    fileList = []
    for fileName in args.fileList:
        fileList += [fileName]
    velRange = args.velRange
    plotFit = args.plotFit

    print('fileName                               RV    error')
    #Loop over the files provided, and fit a radial velocity for each
    for fileName in fileList:
        lsd = pol.read_lsd(fileName)
        
        #Run the radial velocity fitting routine
        fitVel, fitVelErr = lsd.fit_gaussian_rv(velRange, plotFit=plotFit)

        print('{:30} {:10.4f} {:8.4f}'.format(fileName, fitVel, fitVelErr))
    return

#For running as a terminal program
if __name__ == "__main__":
    rv_fit_cli()
