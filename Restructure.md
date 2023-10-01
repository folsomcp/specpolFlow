TODO in order

Note from Vero: I tried to limit the use of external packages as much as possible, but I added the copy package to make a deepcopy of the Mask object in the cleaning routine. 
I could get away with it if you'd like, but it it will be less elegant. 
(Cannot use the constructor directly)

Vero: 
[x] finish the Regions class
    [x] defaults
    [x] read
    [x] test the save, slice, read, and defaults.
    [x] overload the + operator
[x] check the mask class for consistency (and also .clean instead of .clean_mask is class function)
[x] in mask class, use the region class to return mask
[ ] Bz into the LSD class. 
    [x] The cog calculations
    [x] The FAP
    [x] The helper integration functions
    [x] The helper plotting function
    [x] The Bz function itself
    [x] Test the code in the Bz standalone tutorial

Colin:
[x] tidy up the class names linelist and spectrum
[ ] move rvfits to the LSD class
[ ] moving the __main__ to scripts in the CommandLineScript
    wait until V has moved the Bz to the LSD class for that one..
[ ] check whether we can split off some part of ioLSD
* def plot(self, figsize=(10,10), sameYRange=True, plotZeroLevel=True, **kwargs):
    [x] add fig=None, ax=None so that another call can overplot another LSD profile. 
    [x] add a command line wrapper around LSD.plot()
    [ ] maybe add an option for scatter with errors and just lines. 



More TODOs that are for filling in the docs:
[ ] The ExcludeMaskRegions tutorial has all the important code, but needs some actual text. 


*** Colin will check the [:] in the set items for a
lsd[2] = lsd[3] type of command and check if it is needed
-> DONE!

Name changes to be aware of for other scripts:
* lsd_prof => LSD
* mask => Mask
* line_list => LineList
* observation => Spectrum
* mask.set_weights() => Mask.get_weights()
* mask.clean_mask() => Mask.clean()
* observation.write_s() => Spectrum.save()


Some convention:

filename : fname
reading function read_objectName()
use save in classes for the functions that save

Spectrum are in nm. 
Masks are in nm.
ExcludeMaskRegions in nm. 


[ ] Rename ioLSD to something better, and arrange the init for direct access to these functions?

Name for the classes ?
[ ] Don't forget to change name of class in _setitem_
* LSD 
* Mask 
* Spectrum 
* LineList


Class: LSD -> 
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
        [x] add fig=None, ax=None so that another call can overplot another LSD profile. 
        [ ] maybe add an option for scatter with errors and just lines. 

-> read_lsd()


def run_lsdpy(obs=None, mask=None, outName='prof.dat',
         velStart=None, velEnd=None, velPixel=None, 
         normDepth=None, normLande=None, normWave=None,
         removeContPol=None, trimMask=None, sigmaClipIter=None, sigmaClip=None, 
         interpMode=None, outModelName='',
         fLSDPlotImg=None, fSavePlotImg=None, outPlotImgName=None):



Class: Mask
    * _init_
    * def __getitem__(self, key):
    * def __setitem__(self, key, newval):
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
        * get item
        * save
        * [ ] read 
        * [ ] concatenate multiple regions objects


-> read_region()
-> default_region()
    * with velocity range for around the Balmer lines
    * just the telluric








Class: Spectrum

    * def __init__(self,wl, specI, specV, specN1, specN2, specSig, header=None)
    
    * def __getitem__(self, key): 
    * def __setitem__(self, key, newval):
    * def __len__(self):
    
    * def save(self, fname, noHeader=False):

    * [later] Add a split order class function (see cleanMask.py split_order or Colin's normplot) 
    * [later] Add a little wrapper that takes a wavelength (and a range?), cut a line, put in an lsd object and calls lsd.rvfit. 
    
 [x] Colin need to check the \n in the reading of the header to see whether a \n is needed or not when writing the header back in a .s file

 [ ] Would be better is the converters reading other format file, create a spectrum object and write it back with the class.save method. 




-> read_spectrum()
    [x] Right now is only works with the full polarization parameters. Colin will add some functionality for a 2 column version (just wave and flux)



Class: LineList
    *  def __init__(self, ion, wl, loggf, Elo, Jlo, Eup, Jup, landeLo, landeUp,
                 landeEff, rad, stark, waals, depth, configLo, configUp, refs):
    * def __getitem__(self, key):
    * def __setitem__(self, key, newval):

-> read_VALD()



[ ] Put the rvfit into LSD?


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


## To add later:

* ExcludeMaskRegions default for the Pashen and Bracken series, for spirou usage. 
* ExcludeMaskRegions default for the telluric regions in the IR? 