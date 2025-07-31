#!/usr/bin/python3
#
#A wrapper for quickly plotting spectra from .s files and VALD line lists

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

def plot_spec_lines_cli():
    """
    Function for quickly reading and plotting spectra from .s files
    and VALD long format line lists, ran from a terminal.
    """
    #Take command line arguments
    parser = argparse.ArgumentParser(description='Plot a set of spectra and line lists.')
    parser.add_argument("fileList", nargs='*',
                        help="Spectra .s file(s) to plot.")
    parser.add_argument("-l", "--lineList", nargs='*',
                        help="Line list file(s) to plot, in a VALD "
                        "'extract stellar' long format.")
    parser.add_argument("-v", "--velSpec", nargs='*', type=float, default=[0.0],
                        help="Optionally, Doppler shift the spectra by this "
                        "velocity (in km/s). This can be a single velocity for "
                        "all spectra, or a set of velocities one for each "
                        "spectrum.")
    parser.add_argument("-vl", "--velLines", nargs='*', type=float, default=[0.0],
                        help="Optionally, Doppler shift the line list by this "
                        "velocity (in km/s). This can be a single velocity for "
                        "all line lists, or a set of velocities one for each "
                        "line list.")
    parser.add_argument("-s", "--stokes", default='I',
                        choices=['I', 'V', 'N1', 'N2', 'IV'],
                        help="Which Stokes parameter from the spectrum to plot."
                        " can be: I, V, N1, N2, or IV. The 'IV' option "
                        "plots both Stokes I and Stokes V.")    
    parser.add_argument("-e", "--errors", action='store_true',
                        help="Optionally, plot errorbars for the spectra.")
    parser.add_argument("-d", "--depthCut", type=float, default=0.0,
                        help="Only plot lines with VALD depth values greater "
                        "than this cutoff depth.")
    parser.add_argument("-m", "--maxLabels", type=int, default=100,
                        help="Only draw labels for this number of the deepest "
                        "lines. Restricting the number of labels drawn can "
                        "improve performance.")
    parser.add_argument("-g", "--legend", action='store_true',
                        help="Optionally, show a legend of spectra file names.")
    ##parser.add_argument("-s", "--save", default=None,
    ##                    help="Optionally, save the plot to this file")
    args = parser.parse_args()
    #Process the command line parameters
    specList = []
    for fileName in args.fileList:
        specList += [fileName]
    llList = []
    if args.lineList is not None:
        for fileName in args.lineList:
            llList += [fileName]
    velSpec = args.velSpec
    if len(velSpec) == 1: velSpec = velSpec[0]
    velLines = args.velLines
    if len(velLines) == 1: velLines = velLines[0]
    depthCut = args.depthCut
    maxLabels = args.maxLabels
    stokes = args.stokes
    showErr = args.errors
    showLegend = args.legend

    if len(specList) == 0 and len(llList) == 0:
        print('No spectrum or line list given to plot. '
              'Use -h for some help information')
        return
    if len(args.velSpec) > 1 and len(specList) != len(args.velSpec):
        print('Input velocity must either be one single value, '
              'or must have one value for each spectrum')
        return
    if len(args.velLines) > 1 and len(llList) != len(args.velLines):
        print('Input line list velocity must either be one single value, '
              'or must have one value for each line list')
        return
    
    pol.plot_obs_lines(specList, llList, 
                       depthCut=depthCut, maxLabels=maxLabels,
                       velSpec=velSpec, velLines=velLines, stokes=stokes,
                       showErr=showErr, showLegend=showLegend)
    plt.show()
    return

#For running this Python script as a terminal program
if __name__ == "__main__":
    plot_lsd_cli()
