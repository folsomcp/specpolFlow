## @module cleanMaskUI.py
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

def cleanMaskUI(maskName, obsName, outMaskName=None,
                excludeFileName='excludeRanges.dat', batchMode=False):
    """
    Run an interactive tool for LSD line mask cleaning.

    This opens a window with a plot of the line mask, a reference observation,
    and a model LSD spectrum (convolution of a line mask and LSD profile).
    The parameters for the LSD calculation are taken from the defaults in
    the cleanMaskUISubroutines lsdParams() function.

    :param maskName: File name, or a Mask object, for the input line mask
                     to clean.
    :param obsName: File name for the reference observation to compare with.
    :param outMaskName: File name for the output cleaned mask,
                        defaults to [maskName].clean
    :param excludeFileName: File name for a set of regions to be excluded from
                            the line mask (read from and write to). If the file
                            doesn't exist a set of default values will be used.
    :rtype: Mask
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
        regions = maskTools.read_exclude_mask_regions(excludeFileName)
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
            outName=lsdp.outName, velStart=lsdp.velStart, velEnd=lsdp.velEnd,
            velPixel=lsdp.velPixel, normDepth=lsdp.normDepth,
            normLande=lsdp.normLande, normWave=lsdp.normWave,
            removeContPol=lsdp.removeContPol, trimMask=lsdp.trimMask,
            sigmaClipIter=lsdp.sigmaClipIter, sigmaClip=lsdp.sigmaClip,
            interpMode=lsdp.interpMode, outModelName=lsdp.outModelName,
            fLSDPlotImg=lsdp.fLSDPlotImg, fSavePlotImg=lsdp.fSavePlotImg,
            outPlotImgName=lsdp.outPlotImgName)
        
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
        cleanSub.makeWin(fig, ax1, mask, obs, lsdp, pltMaskU, 
                         pltMaskN, pltMaskF, pltModelI, excludeRanges, 
                         excludeFileName, fitDepthFlags)
    
    #Save the final result
    mask.save(outMaskName)
    return mask
