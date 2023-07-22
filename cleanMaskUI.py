#!/usr/bin/python3
## @module cleanMaskUI.py
#
#Interactive tool for cleaning poor lines from an LSD line mask

import numpy as np
from matplotlib.figure import Figure
try:
    import specpolFlow.iolsd as iolsd
except ModuleNotFoundError:
    import iolsd
try:
    import specpolFlow.cleanMaskSubroutines as utils
except ModuleNotFoundError:
    import cleanMaskSubroutines as utils


def cleanMaskUI(maskName, obsName, outMaskName=None,
                excludeFileName='excludeRanges.dat', batchMode=False):
    """
    Run an interactive tool for LSD line mask cleaning.

    This opens a window with a plot of the line mask, a reference observation,
    and a model LSD spectrum (convolution of a line mask and LSD profile).
    The parameters for the LSD model spectrum are read from a file inlsd.dat

    :param maskName: File name for the input line mask to clean.
    :param obsName: File name for the reference observation to compare with.
    :param outMaskName: File name for the output cleaned mask, defaults to [maskName].clean
    :param excludeFileName: File name for a set of regions to be excluded from the line mask (read from and write to). If the file doesn't exist a set of default values will be used.
    :rtype: mask
    """
    
    if outMaskName is None:
        outMaskName = maskName+'.clean'
        
    #Read the input mask
    mask = iolsd.read_mask(maskName)
    
    #Read a list of wavelength regions to exclude from the mask 
    #(if the file exists), or generate one from default values.
    try:
        regions = np.loadtxt(excludeFileName)
        #convert the numpy array to the nested lists used here
        excludeRanges = regions.tolist()
    except OSError:
        excludeRanges = []
        try:
            import specpolFlow.cleanMask as cleanMask
        except ModuleNotFoundError:
            import cleanMask
        regionsDic = cleanMask.default_exclude_regions(1000.)
        #convert the dictionary to the nested lists used here
        for i in range(len(regionsDic['WLStart'])):
            excludeRanges += [[float(regionsDic['WLStart'][i]),
                               float(regionsDic['WLFinish'][i])]]
    #Update the mask with these initial exclude regions
    for region in excludeRanges:
        mask.iuse[(mask.wl >= region[0]) & (mask.wl <= region[1])] = 0
    
    if batchMode == False:
        #Read the comparison observation
        obs = iolsd.read_spectrum(obsName)
        
        #Make an initial LSD model spectrum
        print('Generating test LSD profile...')
        #first save the current updated mask, so we can run LSDpy on it
        mask.save(outMaskName)
        #then run LSDpy
        lsdProf, modelSpec = iolsd.run_lsdpy(obs=obsName, mask=outMaskName,
                                             fLSDPlotImg=0)
        
        #Make a plot of the reference observation, mask, and LSD model spectrum
        fig =  Figure()
        ax1 = fig.add_subplot(111)
        pltObsI = ax1.plot(obs.wl, obs.specI, 'k')
        maskUsed = mask[mask.iuse == 1]
        maskNot = mask[mask.iuse == 0]
        pltMaskN = ax1.vlines(maskNot.wl, ymin=1.0-maskNot.depth, ymax=1.0,
                              colors='r', lw=1)
        pltMaskU = ax1.vlines(maskUsed.wl, ymin=1.0-maskUsed.depth, ymax=1.0,
                              colors='b', lw=1)
        pltMaskF = ax1.vlines([], ymin=[], ymax=1.0,
                              colors='aqua', lw=1) #tab:green
        pltModelI = ax1.plot(modelSpec.wl,modelSpec.specI, 'm', lw=1)
        ax1.set_xlabel('Wavelength')
        ax1.set_ylabel('Flux')

        fitDepthFlags = np.zeros_like(mask.wl)
        
        #Run the main interactive program
        utils.makeWin(fig, ax1, mask, obs, obsName, outMaskName, pltMaskU, 
                      pltMaskN, pltMaskF, pltModelI, excludeRanges, 
                      excludeFileName, fitDepthFlags)
    
    #Save the final result
    mask.save(outMaskName)
    return mask


###############################################################
#For running cleanMaskUI as a terminal program
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Interactively remove poor lines from an LSD line mask.  This plots the mask, a comparison observation, and the model LSD spectrum (the convolution of the mask and LSD profile). For the model LSD spectrum calculation, LSDpy requires an inlsd.dat file for most input parameters.')
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
    
    cleanMask = cleanMaskUI(maskName, obsName, outMaskName=outMaskName,
                            excludeFileName=excludeFileName,
                            batchMode=batchMode)
    
