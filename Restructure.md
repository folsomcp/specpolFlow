Some convention:

filename : fname
reading function read_objectName()
use save in classes for the functions that save


[ ] Rename ioLSD to something better, and arrange the init for direct access to these functions?

Name for the classes ?
[ ] Don't forget to change name of class in _setitem_
* LSD 
* Mask 
* Spectrum 
* LineList


Class: lsd_obj -> 
    * _init_
    * def save(self, fname)
    * def __len__(self)
    * def __getitem__(self, key):
    * def __setitem__(self, key, newval):
    
    * def norm(self, normValue):
    * def vshift(self, velShift):  ??? name consistency in parameters?

    * def set_weights(self, wint_old, wpol_old, wint_new, wpol_new):
    * def scale(self, scale_int, scale_pol):

    * def plot(self, figsize=(10,10), sameYRange=True, plotZeroLevel=True, **kwargs):
        [ ] add fig=None, ax=None so that another call can overplot another LSD profile. 
        [ ] maybe add an option for scatter with errors and just lines. 

-> read_lsd()


def run_lsdpy(obs=None, mask=None, outName='prof.dat',
         velStart=None, velEnd=None, velPixel=None, 
         normDepth=None, normLande=None, normWave=None,
         removeContPol=None, trimMask=None, sigmaClipIter=None, sigmaClip=None, 
         interpMode=None, outModelName='',
         fLSDPlotImg=None, fSavePlotImg=None, outPlotImgName=None):



Class: mask_obj
    * _init_
    * def __getitem__(self, key):
    * def __setitem__(self, key, newval):
    * def __len__(self):
    * def __len__(self):
    * def prune(self):
    * def get_weights(self, normDepth, normWave, normLande):
    * def save(self, fname): 

    *[ ] def clean_mask(self,**exclude_region**):
            NOTE: change to return the clean mask instead of doing it in place
    *[ ] wrapper that takes a region file something that takes a file 
    *[ ] wrapper function input mask file and region file. 


-> read_mask()

[ ] Class: ExcludeMaskRegion
        * set item
        * write item
        * write
        * read 
        * concatenate multiple regions objects


-> read_region()
-> default_region()
    * with velocity range for around the Balmer lines
    * just the telluric















Class: spectrum_obj

    * [ ] Add a split order class function (see cleanMask.py split_order or Colin's normplot) 



-> read_spectrum()

Class: line_list_obj

-> read_VALD()


[ ] Put the rvfit into LSD?
[ ] Add a rvfit to spectrum (putting the line into a lsd object and calling rvfit?)


[ ] Clean Mask needs some more structure


## Standalone functions

[maybe] plot spectrum, mask, and regions and save (still version of the clean mask gui)
    see: Excluded_Regions_Visual(input_file, input_mask, WLRegions): in cleanMask.py

# Command line tools:


## GUI only tools

- cleanMaskUI
[ ] Can we run this from an open python shell
[ ] Make the dymanic use of LSDpy better ? (it wants to write file...)


## CommandLineScripts/ 
- Bz
- rvfit
- makeMask


## Structure of the repository at the end
[ ] Check first how to make pip work for installation. 

repo / specpolFlow / all of the source files
     / CommandLineScripts
     / Tutorials

[ ] arrange the __init__ file the way we want it. 
