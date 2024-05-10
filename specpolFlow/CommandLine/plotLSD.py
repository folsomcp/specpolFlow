#!/usr/bin/python3
#
#A wrapper for quickly plotting a set of LSD profiles from the command line

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

def plot_lsd_cli():
    """
    Function for quickly reading and plotting LSD profiles from a terminal.
    """
    #Take input file names as command line arguments,
    #with some additional optional control parameters.
    parser = argparse.ArgumentParser(description='Plot a set of LSD profiles.')
    parser.add_argument("fileList", nargs='*',
                        help="LSD profile files to plot, can be more than one file.")
    parser.add_argument("-l", "--legend", action='store_true',
                        help="Optionally, show a legend of file names")
    parser.add_argument("-s", "--save", default=None,
                        help="Optionally, save the plot to this file")
    args = parser.parse_args()
    #Process the command line parameters
    fileList = []
    for fileName in args.fileList:
        fileList += [fileName]
    showLegend = args.legend
    saveName = args.save

    fig = None
    axs = None
    for fileName in fileList:
        lsd = pol.read_lsd(fileName)
        fig, axs = lsd.plot(fig=fig, ls='-', label=fileName)

    if showLegend:
        axs[-1].legend(loc='lower left')
    if saveName is not None:
        print('saving to', saveName)
        fig.savefig(saveName)
    
    plt.show()
    return

#For running this Python script as a terminal program
if __name__ == "__main__":
    plot_lsd_cli()
