"""
An interactive program for cleaning problem lines from a LSD line mask
and tweaking (fitting) depths of some lines.  This opens its main
program in another window.  The LSDpy package must be installed
(in your Python path) to use this program.  
"""

import numpy as np
from matplotlib.figure import Figure

from . import cleanMaskUISubroutines as cleanSub
from . import profileLSD
from . import mask as maskTools
from . import obsSpec

###################################
###################################

def cleanMaskUI(maskName, obsName, outMaskName=None, inExcludeName=None,
                outExcludeName='excludeRanges.dat', batchMode=False):
    """
    Run an interactive tool for LSD line mask cleaning.

    This opens a window with a plot of the line mask, a reference observation,
    and a model LSD spectrum (convolution of a line mask and LSD profile).
    The parameters for the LSD calculation are taken from the defaults in
    the cleanMaskUISubroutines lsdParams() function.

    This saves a copy of the cleaned mask to a file, and also returns
    a copy as a Mask object and the regions excluded from the mask
    as an ExcludeMaskRegions object.

    :param maskName: File name, or a Mask object, for the input line mask
                     to clean.
    :param obsName: File name for the reference observation to compare with.
    :param outMaskName: File name for the output cleaned mask,
                        defaults to [maskName].clean
    :param inExcludeName: File name for reading a set of regions to be excluded
                          from the line mask.
    :param outExcludeName: File name for saving the set of regions excluded 
                           from the mask. (may be the same as inExcludeName)
    :param batchMode: Flag for skipping the GUI. False = open the GUI (the default).
                      True = skip the GUI, just read inExcludeName and apply 
                      those exclude regions.
    :return: a Mask object with the cleaned mask,
             and an ExcludeMaskRegions with the selected exclude regions.
    """
    
    if outMaskName is None:
        if isinstance(maskName, maskTools.Mask):
            outMaskName = 'mask.clean'
        else:
            outMaskName = maskName+'.clean'
    
    #Read the input mask
    if isinstance(maskName, maskTools.Mask):
        mask = maskName
    else:
        mask = maskTools.read_mask(maskName)
    
    #Read a list of wavelength regions to exclude from the mask 
    #(if the file exists), or generate one from default values.
    try:
        if inExcludeName is None: raise OSError #fall back to except block
        regions = maskTools.read_exclude_mask_regions(inExcludeName)
    except OSError:
        regionsBalmer = maskTools.get_Balmer_regions_default(1000.)
        regionsTelluric = maskTools.get_telluric_regions_default()
        regions = regionsBalmer + regionsTelluric
    #Update the mask with these initial exclude regions
    mask = mask.clean(regions)
    #Make a list of lists from the exclude regions, used internally.
    #The inner list should contain only 2 entries: a start and end wavelength.
    #(It could be cleaner to use ExcludeMaskRegions internally, but one would
    #have to be careful to modify the object, not replace it with a copy)
    excludeRanges = []
    for region in regions:
        excludeRanges += [[region.start, region.stop]]
    #Sort the exclude regions by starting wavelength (to be safe)
    excludeRanges = sorted(excludeRanges, key=lambda ran: ran[0])
    
    if batchMode == False:
        #Read the comparison observation
        obs = obsSpec.read_spectrum(obsName)
        
        #Make an initial LSD model spectrum
        print('Generating test LSD profile...')
        lsdp = cleanSub.lsdParams(obsName, outMaskName)
        #first save the current updated mask, so we can run LSDpy on it
        mask.save(outMaskName)
        #then run LSDpy, with specific 'default' parameters
        lsdProf, modelSpec = profileLSD.run_lsdpy(obs=lsdp.obs, mask=lsdp.mask,
            outLSDName=lsdp.outName, velStart=lsdp.velStart, velEnd=lsdp.velEnd,
            velPixel=lsdp.velPixel, normDepth=lsdp.normDepth,
            normLande=lsdp.normLande, normWave=lsdp.normWave,
            removeContPol=lsdp.removeContPol, trimMask=lsdp.trimMask,
            sigmaClipIter=lsdp.sigmaClipIter, sigmaClip=lsdp.sigmaClip,
            outModelName=lsdp.outModelName,
            plotLSD=lsdp.plotLSD, outPlotLSDName=lsdp.outPlotLSDName)
        
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
                              colors='aqua', lw=1)
        pltModelI = ax1.plot(modelSpec.wl,modelSpec.specI, 'm', lw=1)
        ax1.set_xlabel('Wavelength')
        ax1.set_ylabel('Flux')

        #initialize an array of flags for fitting LSD line depths
        fitDepthFlags = np.zeros_like(mask.wl, dtype=int)
        
        #Run the main interactive program
        cleanSub.makeWin(fig, ax1, mask, obs, lsdp, lsdProf, pltMaskU,
                         pltMaskN, pltMaskF, pltModelI, excludeRanges,
                         outExcludeName, fitDepthFlags)

    #make an ExcludeMaskRegions object to return
    wlStarts = np.array([ran[0] for ran in excludeRanges])
    wlEnds = np.array([ran[1] for ran in excludeRanges])
    labels = np.array(['cleanMaskUI']*len(excludeRanges),dtype=object)
    regions = maskTools.ExcludeMaskRegions(wlStarts, wlEnds, labels)

    #Save the final result
    mask.save(outMaskName)
    return mask, regions
