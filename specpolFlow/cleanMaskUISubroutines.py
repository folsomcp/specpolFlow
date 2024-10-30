#Additional subroutines used inside cleanMaskUI.py
#Not really meant to be stand alone utilities.

import numpy as np
import scipy
import tkinter as tk
import tkinter.ttk as ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#Try to import the default matplotlib toolbar
try: #for matplotlib 3.x
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
    _mplToolbar=3
except ImportError:  #If that fails this probably isn't matplotlib 3.x
    try: #for matplotlib 2.x
        from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
        _mplToolbar=2
    except ImportError: #if that fails, just don't use this extra toolbar
        _mplToolbar=0

from . import profileLSD
from . import mask as maskTools
from . import obsSpec

def makeWin(fig, ax, mask, obs, lsdp, lsdProf, pltMaskU, pltMaskN,
            pltMaskF, pltModelI, excludeRanges, outExcludeName, fitDepthFlags):
    #Build GUI with tkinter
    root = tk.Tk(className='Clean Masks')
    #the className seems to set an icon title (at least in Ubuntu)
    root.title("Clean Masks")

    ##For an icon image (used by window managers)
    #try:
    #    imgIcon = tk.Image("photo", file="./iconNorm.png")
    #    root.iconphoto(True, imgIcon)
    #except:
    #    pass

    #Generate a tk.DrawingArea using matplotlib's Tk backend,
    #from the given figure
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw() #this draw call may not be strictly necessary? 
    canvasWidget = canvas.get_tk_widget()
    canvasWidget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    canvas.figure.tight_layout()  #try to use matplotlib's tighter boarders

    #Add a matplotlib's standard toolbar.
    if _mplToolbar > 0:
        if _mplToolbar == 2:
            toolbarMpl = NavigationToolbar2TkAgg(canvas, root)
        if _mplToolbar == 3:
            toolbarMpl = NavigationToolbar2Tk(canvas, root)
        #Note, making the toolbar calls the parent container's .pack()
        toolbarMpl.update() #not strictly necessary
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    #Add custom UI controls for the plot
    
    #Create a little event manager object for processing
    #callbacks for view changing buttons
    viewM = viewFuncs(ax, canvas)
    #Make a frame to hold our tool buttons,
    #and allow us to use use the tk grid geometry manager for them
    tools = ttk.Frame(root, padding="3 3 3 3")
    tools.pack(side=tk.BOTTOM, fill=tk.X)
    #and put view changing buttons in their own frame
    viewtools = ttk.Frame(tools)
    viewtools.grid(row=3, column=0, columnspan=8, sticky=tk.W, padx=2)
    #
    butAuto = ttk.Button(master=viewtools, text="auto",
                         command=viewM.autoScale, width=7)
    ToolTip(butAuto, "Auto-scale axes ranges [a]")
    butAuto.grid(row=0, column=0, sticky=(tk.E,tk.W))
    butAutoY = ttk.Button(master=viewtools, text="auto-y",
                          command=viewM.autoScaleY, width=7)
    ToolTip(butAutoY, "Auto-scale y-axis range [A]")
    butAutoY.grid(row=0, column=1, sticky=(tk.E,tk.W))
    butZin = ttk.Button(master=viewtools, text="zoom",
                        command=viewM.zoomRec, width=5)
    ToolTip(butZin, "Zoom to selected region [z]")
    butZin.grid(row=0, column=2, sticky=(tk.E,tk.W))
    butZin = ttk.Button(master=viewtools, text="Zin",
                        command=viewM.zoomIn, width=5)
    ToolTip(butZin, "Zoom in [i]")
    butZin.grid(row=0, column=3, sticky=(tk.E,tk.W))
    butZout = ttk.Button(master=viewtools, text="Zout",
                         command=viewM.zoomOut, width=5)
    ToolTip(butZout, "Zoom out [o]")
    butZout.grid(row=0, column=4, sticky=(tk.E,tk.W))
    butLeft = ttk.Button(master=viewtools, text="\u2190",
                         command=viewM.panLeft, width=3)
    ToolTip(butLeft, "Pan left [arrow key]")
    butLeft.grid(row=0, column=5, sticky=(tk.E,tk.W))
    butRight = ttk.Button(master=viewtools, text="\u2192",
                          command=viewM.panRight, width=3)
    ToolTip(butRight, "Pan right [arrow key]")
    butRight.grid(row=0, column=6, sticky=(tk.E,tk.W))
    butUp = ttk.Button(master=viewtools, text="\u2191",
                       command=viewM.panUp, width=3)
    ToolTip(butUp, "Pan up [arrow key]")
    butUp.grid(row=0, column=7, sticky=(tk.E,tk.W))
    butDown = ttk.Button(master=viewtools, text="\u2193",
                         command=viewM.panDown, width=3)
    ToolTip(butDown, "Pan down [arrow key]")
    butDown.grid(row=0, column=8, sticky=(tk.E,tk.W))
    #Allow .grid() managed buttons to expand/contract if necessary
    viewtools.rowconfigure(0, weight=1)
    viewtools.columnconfigure(0, weight=1)
    viewtools.columnconfigure(1, weight=1)
    viewtools.columnconfigure(2, weight=1)
    viewtools.columnconfigure(3, weight=1)
    viewtools.columnconfigure(4, weight=1)
    viewtools.columnconfigure(5, weight=1)
    viewtools.columnconfigure(6, weight=1)
    viewtools.columnconfigure(7, weight=1)
    viewtools.columnconfigure(8, weight=1)
    
    #Set some general key bindings
    root.bind('<Key-Left>', viewM.panLeft)
    root.bind('<Key-Right>', viewM.panRight)
    root.bind('<Key-Up>', viewM.panUp)
    root.bind('<Key-Down>', viewM.panDown)
    root.bind('a', viewM.autoScale)
    root.bind('A', viewM.autoScaleY)
    root.bind('z', viewM.zoomRec)
    root.bind('i', viewM.zoomIn)
    root.bind('o', viewM.zoomOut)


    #let the .grid() managed cells expand
    tools.columnconfigure(6, weight=10)
    
    #Set regions to be used for adjusting line depths
    selectFitDepth = uiSelectFitDepth(canvas, mask, pltMaskU, pltMaskN, pltMaskF,
                                      excludeRanges, fitDepthFlags)
    butSelRange = ttk.Button(master=tools, text='select lines\nto fit depth',
                             command=selectFitDepth.runSpanSelect)
    ToolTip(butSelRange, 'Select lines to fit for depth in the mask')
    butSelRange.grid(row=0, column=0, sticky=tk.W, padx=2, pady=2)
    #Set regions to not be used for adjusting line depths
    unselectFitDepth = uiUnselectFitDepth(canvas, mask, pltMaskU, pltMaskN, pltMaskF,
                                          excludeRanges, fitDepthFlags)
    butUnselRange = ttk.Button(master=tools, text='unselect lines\nto fit depth',
                             command=unselectFitDepth.runSpanSelect)
    ToolTip(butUnselRange, 'Unselect lines to fit for depth in the mask')
    butUnselRange.grid(row=0, column=1, sticky=tk.W, padx=2, pady=2)
    #Link the select and unselect buttons so they can turn each other off,
    #done after creating the include & exclude buttons so they all can be linked
    #Apply the depth fitting
    fitDepthsM = fitDepths(mask, obs, lsdp, lsdProf, fitDepthFlags, root,
                           canvas, pltMaskU, pltMaskN, pltMaskF)
    butFitDepths = ttk.Button(master=tools, text='fit\ndepths',
                               command=fitDepthsM.runFit)
    ToolTip(butFitDepths, 'Fit the selected line depths, using the current LSD profile and observation')
    butFitDepths.grid(row=0, column=3, sticky=tk.W, padx=2, pady=2)
    #Undo the depth fitting
    undoDepthsM = undoDepths(mask, fitDepthFlags, canvas,
                             pltMaskU, pltMaskN, pltMaskF)
    butUndoDepths = ttk.Button(master=tools, text='undo\nfit',
                               command=undoDepthsM.undo)
    ToolTip(butUndoDepths, 'Undo all line depth fitting of currently selected lines, reverting them to the initial mask depths')
    butUndoDepths.grid(row=0, column=4, sticky=tk.W, padx=2, pady=2)
    
    speratorFitting = ttk.Separator(tools, orient=tk.VERTICAL)
    speratorFitting.grid(row=0, column=5, sticky='NS', padx=2, pady=2)

    #Set regions to be included in the mask
    incRangeM = uiIncludeRange(canvas, mask, pltMaskU, pltMaskN, pltMaskF,
                               excludeRanges, fitDepthFlags)
    butIncRange = ttk.Button(master=tools, text='include\nlines',
                             command=incRangeM.runSpanSelect)
    ToolTip(butIncRange, 'Include selected lines in the mask')
    butIncRange.grid(row=0, column=7, sticky=tk.E, padx=2, pady=2)
    #Set regions to be excluded in the mask
    excRangeM = uiExcludeRange(canvas, mask, pltMaskU, pltMaskN, pltMaskF,
                               excludeRanges, fitDepthFlags)
    butExcRange = ttk.Button(master=tools, text='exclude\nlines',
                             command=excRangeM.runSpanSelect)
    ToolTip(butExcRange, 'Exclude selected lines from the mask')
    butExcRange.grid(row=0, column=8, sticky=tk.E, padx=2, pady=2)
    #Link the include and exclude buttons so they can turn each other off,
    #done in the next block where all related buttons are linked
    #Save the exclude ranges to a text file
    saveRangesM = saveRanges(outExcludeName, excludeRanges)
    butSaveRanges = ttk.Button(master=tools, text='save\nranges',
                               command=saveRangesM.saveToFile)
    ToolTip(butSaveRanges, 'Save the excluded wavelength ranges to the {:} file'.format(outExcludeName))
    butSaveRanges.grid(row=0, column=9, sticky=tk.E, padx=2, pady=2)

    #Link the select, unselect, include, and exclude buttons
    #so they can all turn each other off, and have only one mode active
    selectFitDepth.linkButton(butSelRange, [unselectFitDepth, incRangeM, excRangeM])
    unselectFitDepth.linkButton(butUnselRange, [selectFitDepth, incRangeM, excRangeM])
    incRangeM.linkButton(butIncRange, [excRangeM, selectFitDepth, unselectFitDepth])
    excRangeM.linkButton(butExcRange, [incRangeM, selectFitDepth, unselectFitDepth])
    
    #Modify LSD control parameters
    modifyLSDparM = modifyLSDpar(root, lsdp)
    butModLSD = ttk.Button(master=tools, text='LSD param.',
                              command=modifyLSDparM.openWindow)
    ToolTip(butModLSD, 'Modify the parameters related to the LSD calculation')
    butModLSD.grid(row=3, column=8, sticky=tk.E, padx=2)
    
    #Update the LSD calculation and the plotted spectrum
    updateLSDM = updateLSD(canvas, root, mask, lsdp, lsdProf, pltModelI)
    butUpdateLSD = ttk.Button(master=tools, text='update LSD',
                              command=updateLSDM.rerunLSD)
    ToolTip(butUpdateLSD, 'Save the mask, run LSD, and update the model spectrum')
    butUpdateLSD.grid(row=3, column=9, sticky=(tk.N, tk.S, tk.E, tk.W), padx=2)
    
    #Run the main loop!
    root.mainloop()


def combineRanges(ranges1, ranges2):
    #ranges1 and ranges2 are lists of ranges, each containing 2 entry lists 
    #of start and end wavelengths for a range.
    #This combines the two sets of ranges and removes overlap by merging ranges.
    #combine ranges and sort by starting wavelength
    rangesT = sorted(ranges1 + ranges2, key=lambda ran: ran[0])
    #recursively merge adjacent entries
    update = True
    while update == True:
        update = False
        skip = False
        rangesO = []
        #Loop over all but the last entry (compare this entry with the next)
        for i, r1 in enumerate(rangesT[:-1]):
            #Skip this entry if it's been merged with the previous
            if skip:
                skip = False
                continue
            r2 = rangesT[i+1]
            #if there is overlap merge these windows
            if r1[1] >= r2[0]:
                update = True
                #if the next range ends after this range
                if r2[1] > r1[1]:
                    rangesO += [[r1[0], r2[1]]]
                    skip = True
                #if the next range ends before this range (is smaller)
                else:
                    rangesO += [[r1[0], r1[1]]]
                    skip = True
            #if there is no overlap
            else:
                rangesO += [r1]
        #Include the last entry
        if(rangesO[-1][1] < rangesT[-1][1]):
            rangesO += [rangesT[-1]]
        rangesT = rangesO
    return rangesO


def removeRange(ranges1, range2):
    #for a set of 'exclude wavelength ranges' in ranges1,
    #remove region in range2, so that the region in range2 is 'included'
    rangesT = sorted(ranges1, key=lambda ran: ran[0])
    rangesO = []
    for i, r1 in enumerate(rangesT):
        #if this 'exclude' range starts in the new 'include' range
        if (range2[0] <= r1[0]) and (r1[0] < range2[1]):
            #and this exclude range ends after the new include range
            if range2[1] < r1[1]:
                rangesO += [[range2[1], r1[1]]]
            #Otherwise the exclude region starts and ends in the include range
            #so skip it
        #if this 'exclude' range ends in the new 'include' range
        elif (range2[0] <= r1[1]) and (r1[1] < range2[1]):
            #and this exclude range starts before the new include range
            if r1[0] < range2[0]:
                rangesO += [[r1[0], range2[0]]]
        #If this 'exclude' ranges spans the new include range
        elif (r1[0] < range2[0]) and (range2[1] < r1[1]):
            if range2[0] < range2[1]: #sanity check
                rangesO += [[r1[0], range2[0]]]
                rangesO += [[range2[1], r1[1]]]
        #otherwise this exclude range should be outside the include range
        else:
            rangesO += [r1]
    return rangesO


#A fairly simple, fairly general tooltip
class ToolTip(object):
    """
    Create a tooltip for a given widget, using Tkinter.
    This is a fairly simple class, but it wraps text at a specified length,
    puts the tooltip below the widget, and moves the tooltip to keep it 
    within the screen limits.  
    Based on the recipe from:
    http://www.voidspace.org.uk/python/weblog/arch_d7_2006_07_01.shtml#e387
    and some discussion from: 
    https://stackoverflow.com/questions/3221956/how-do-i-display-tooltips-in-tkinter
    """
    def __init__(self, widget, text='widget info',
                 waittime = 750, wraplength = 180):
        self.waittime = waittime       #miliseconds
        self.wraplength = wraplength   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.waitid = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.waitid = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        waitid = self.waitid
        self.waitid = None
        if waitid:
            self.widget.after_cancel(waitid)

    def showtip(self, event=None):
        #Create a window showing the tooltip
        #shouldn't be necessary, but just in case remove any existing window
        self.hidetip() 
        #Make the tooltip appear outside the mouse position (below the button)
        #helps prevent accidental loops of showing and hiding the tooltip.  
        x = self.widget.winfo_pointerx() + 1
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 1
        #Creates a toplevel window to hold the tooltip
        #and remove all the standard window parts.
        self.tw = tk.Toplevel(self.widget)
        self.tw.wm_overrideredirect(True)
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background="#ffffff", relief='solid', borderwidth=1,
                       wraplength = self.wraplength)
        #extra 1 pixel spacing around the text
        label.pack(ipadx=1) 
        #Check that the tooltip will fit on the screen,
        #and if it won't move the window containing it.
        if (x + label.winfo_reqwidth() + 2) > self.widget.winfo_screenwidth():
            x = self.widget.winfo_screenwidth() - (label.winfo_reqwidth() + 2)
        if (y + label.winfo_reqheight()) > self.widget.winfo_screenheight():
            y = self.widget.winfo_rooty() - label.winfo_reqheight() - 2
        #The new window seems to not draw until update_idletasks is called
        #or when this function exits and control returns to the main loop.
        #So we can change the after everything else.
        self.tw.wm_geometry("+%d+%d" % (x, y))

    def hidetip(self):
        #set self.tw to None first, in case the destroy method throws an exception
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()


#Functions for panning and Zooming
#A little event manager thing for dealing with callbacks
class viewFuncs:
    def __init__(self, ax, canvas):
        self.ax = ax
        #Get the (first) Line2D object from this Axes object,
        #which should be the plot() of the observation.
        self.plObs = ax.get_lines()[0]
        self.canvas = canvas
        self.zoomRecActive = False
        self.lastrect = None
        return
    def autoScale(self, *event):
        fracSpace = 0.05
        xmin = np.min(self.plObs.get_xdata())
        xmax = np.max(self.plObs.get_xdata())
        ymin = np.min(self.plObs.get_ydata())
        ymax = np.max(self.plObs.get_ydata())
        xspan = xmax-xmin
        yspan = ymax-ymin
        self.ax.set_xlim(xmin-xspan*fracSpace, xmax+xspan*fracSpace)
        self.ax.set_ylim(ymin-yspan*fracSpace, ymax+yspan*fracSpace)
        self.canvas.draw()
        return
    def autoScaleY(self, *event):
        fracSpace = 0.05
        xmin, xmax = self.ax.get_xlim()
        iInViewX = (self.plObs.get_xdata() > xmin) & (self.plObs.get_xdata() < xmax)
        ymin = np.min(self.plObs.get_ydata()[iInViewX])
        ymax = np.max(self.plObs.get_ydata()[iInViewX])
        yspan = ymax-ymin
        self.ax.set_ylim(ymin-yspan*fracSpace, ymax+yspan*fracSpace)        
        self.canvas.draw()
        return
    def zoomIn(self, *event):
        fracZoom = 0.1
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        xspan = xmax-xmin
        yspan = ymax-ymin
        modFracZoom = fracZoom/(1.+2.*fracZoom) #modify so zoom in reverses zoom out exactly
        self.ax.set_xlim(xmin+xspan*modFracZoom, xmax-xspan*modFracZoom)
        self.ax.set_ylim(ymin+yspan*modFracZoom, ymax-yspan*modFracZoom)
        self.canvas.draw()
        return
    def zoomOut(self, *event):
        fracZoom = 0.1
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        xspan = xmax-xmin
        yspan = ymax-ymin
        self.ax.set_xlim(xmin-xspan*fracZoom, xmax+xspan*fracZoom)
        self.ax.set_ylim(ymin-yspan*fracZoom, ymax+yspan*fracZoom)
        self.canvas.draw()
        return
    def panLeft(self, *event):
        fracPan = 0.2
        xmin, xmax = self.ax.get_xlim()
        xspan = xmax-xmin
        self.ax.set_xlim(xmin-xspan*fracPan, xmax-xspan*fracPan)
        self.canvas.draw()
        return
    def panRight(self, *event):
        fracPan = 0.2
        xmin, xmax = self.ax.get_xlim()
        xspan = xmax-xmin
        self.ax.set_xlim(xmin+xspan*fracPan, xmax+xspan*fracPan)
        self.canvas.draw()
        return
    def panUp(self, *event):
        fracPan = 0.2
        ymin, ymax = self.ax.get_ylim()
        yspan = ymax-ymin
        self.ax.set_ylim(ymin+yspan*fracPan, ymax+yspan*fracPan)
        self.canvas.draw()
        return
    def panDown(self, *event):
        fracPan = 0.2
        ymin, ymax = self.ax.get_ylim()
        yspan = ymax-ymin
        self.ax.set_ylim(ymin-yspan*fracPan, ymax-yspan*fracPan)
        self.canvas.draw()
        return
    
    def zoomRec(self, *event):
        #Zoom to a rectangle selected by the user with the mouse
        if not self.zoomRecActive:
            self.zoomRecActive = True
            self.canvasWidget = self.canvas.get_tk_widget()
            self.numClick = 0
            self.x0rec, self.y0rec, self.x1rec, self.y1rec = 0., 0., 0., 0.
            self.oldBindClick = self.canvasWidget.bind('<Button-1>')
            self.canvasWidget.bind('<Button-1>', self.onClickRec)
            self.oldBindClick3 = self.canvasWidget.bind('<Button-3>')
            self.canvasWidget.bind('<Button-3>', self.onClickCancel)
            self.oldCursor = self.canvasWidget.cget('cursor')
            self.canvasWidget.config(cursor='crosshair')
        return
    def zoomRecDeactivate(self):
        self.canvasWidget.bind('<Button-1>', self.oldBindClick)
        self.canvasWidget.bind('<Button-3>', self.oldBindClick3)
        self.canvasWidget.config(cursor=self.oldCursor)
        self.zoomRecActive = False
        #when deactivating also set the plot view to the selected range
        if self.x1rec > 0:
            #get the canvas position in window pixel coordinates
            x0full, x1full = self.canvas.figure.axes[0].bbox.intervalx
            y0inv, y1inv = self.canvas.figure.axes[0].bbox.intervaly
            #and correct for matplotlib using y in the opposite direction from tk
            height = self.canvas.figure.bbox.height
            y0full = height - y1inv
            y1full = height - y0inv
            axesxrange = x1full - x0full
            axesyrange = y1full - y0full
            #get the current view range in data coordinates
            axdatax0, axdatax1, axdatay0, axdatay1 = self.canvas.figure.axes[0].axis()
            dataxrange = axdatax1 - axdatax0
            datayrange = axdatay1 - axdatay0
            #convert the selected window pixel range into data coordinates
            datax0 = (self.x0rec - x0full)/(axesxrange)*(dataxrange) + axdatax0
            datax1 = (self.x1rec - x0full)/(axesxrange)*(dataxrange) + axdatax0
            datay0 = (self.y0rec - y0full)/(axesyrange)*(-datayrange) + axdatay1
            datay1 = (self.y1rec - y0full)/(axesyrange)*(-datayrange) + axdatay1
            #set the plotted range into the selected range in data coordinates
            if abs(datax1 - datax0) > 0. and abs(datay1 - datay0) > 0.:
                self.ax.set_xlim(min(datax0, datax1), max(datax0, datax1))
                self.ax.set_ylim(min(datay0, datay1), max(datay0, datay1))
                self.canvas.draw()
        return
    def onClickRec(self, event):
        #On the first click start drawing a selection rectangle
        if self.numClick == 0:
            self.numClick += 1
            self.oldBindMove = self.canvasWidget.bind('<Motion>')
            self.canvasWidget.bind('<Motion>', self.onMoveRec)
            self.x0rec = event.x
            self.y0rec = event.y
            self.lastrect = self.canvasWidget.create_rectangle(
                self.x0rec, self.y0rec, self.x0rec+1, self.y0rec+1, width=1)
        #On the second click deactivate and finish the zoom
        else:
            self.x1rec = event.x
            self.y1rec = event.y
            if self.lastrect is not None:
                self.canvasWidget.delete(self.lastrect)
            self.canvasWidget.bind('<Motion>', self.oldBindMove)
            self.zoomRecDeactivate()
        return
    def onClickCancel(self, event):
        #cancel the zoom to rectangle, usually in right click
        self.canvasWidget.delete(self.lastrect)
        if self.numClick > 0:
            self.canvasWidget.bind('<Motion>', self.oldBindMove)
        self.zoomRecDeactivate()
        return
    def onMoveRec(self, event):
        #delete and redraw the rectangle when the mouse moves
        x1 = event.x
        y1 = event.y
        self.canvasWidget.delete(self.lastrect)
        self.lastrect = self.canvasWidget.create_rectangle(
            self.x0rec, self.y0rec, x1, y1, width=1)
        return


#A base class for uiExcludeRange and uiIncludeRange
class uiRangeBase:
    def __init__(self, canvas, mask, pltMaskU, pltMaskN, pltMaskF,
                 excludeRanges, fitDepthFlags):
        self.canvas = canvas
        self.canvasWidget = canvas.get_tk_widget()
        self.mask = mask
        self.pltMaskU = pltMaskU
        self.pltMaskN = pltMaskN
        self.pltMaskF = pltMaskF
        self.excludeRanges = excludeRanges
        self.fitDepthFlags = fitDepthFlags
        self.active = False
        self.rangeSelect = rangeSelect(canvas, self)

        #Set up a style for the 'actively selecting a range' button state
        #Seems like I need to set the base style changes first,
        #then use map to set state specific differences,
        #since map doesn't provide easy access to being in
        #not any of the special states (only not a specific state)
        self.style = ttk.Style()
        self.style.configure('ActRange.TButton', relief='sunken',
                        background='#80ccff')
        self.style.map('ActRange.TButton',
                  relief=[('active','sunken')],
                  background=[('pressed','#eeeeee'),('active','#80ccff')])
        #note: order matters for map here, 
        #when the button is 'active' (mouse over) and also being pressed
        #The scope of the style names seems to be pretty wide...

    def linkButton(self, button, otherRangeList):
        self.button = button
        self.otherRangeList = otherRangeList
    def runSpanSelect(self):
        #First turn off the other linked buttons (only one mode active!)
        for otherRange in self.otherRangeList:
            if (otherRange.active == True):
                otherRange.deactivate()
        self.active = not self.active
        if self.active:
            self.button.configure(style='ActRange.TButton')
            #get a string of tk info that we can use later to restore the current bindings
            self.oldBindClick = self.canvasWidget.bind('<Button-1>')
            #override the current binding
            self.fidClick = self.canvasWidget.bind('<Button-1>',
                                                   self.rangeSelect.onClick)
            self.oldBindClick3 = self.canvasWidget.bind('<Button-3>')
            self.canvasWidget.bind('<Button-3>', self.rangeSelect.cancel)
        else:
            self.deactivate()
    def deactivate(self):
        #stop span selector...
        self.rangeSelect.deactivate()
        self.button.configure(style='TButton')
        #reset to the previous bindings once we are done
        #(Safer to bind to something different, since unbind seems to have poorly defined behaviour)
        self.canvasWidget.bind('<Button-1>', self.oldBindClick)
        self.canvasWidget.bind('<Button-3>', self.oldBindClick3)
        self.active = False


class uiIncludeRange(uiRangeBase):
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update the mask and 'exclude regions'
        self.mask.iuse[(self.mask.wl >= minval) & (self.mask.wl <= maxval)] = 1
        self.excludeRanges[:] = removeRange(self.excludeRanges, [minval,maxval])
        
        updateLinePlots(self.mask, self.fitDepthFlags,
                        self.pltMaskU, self.pltMaskN, self.pltMaskF)
        self.canvas.draw()

        
class uiExcludeRange(uiRangeBase):
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update the mask and 'exclude regions'
        self.mask.iuse[(self.mask.wl >= minval) & (self.mask.wl <= maxval)] = 0
        self.excludeRanges[:] = combineRanges(self.excludeRanges,
                                              [[minval,maxval]])
        
        updateLinePlots(self.mask, self.fitDepthFlags,
                        self.pltMaskU, self.pltMaskN, self.pltMaskF)
        self.canvas.draw()


def updateLinePlots(mask, fitDepthFlags, pltMaskU, pltMaskN, pltMaskF):
    #Update the plot of lines from the mask,
    #with plots for lines used, lines excluded,
    #and lines with depths to fit.
    maskUsed = mask[mask.iuse == 1]
    maskNot = mask[mask.iuse == 0]
    #update the data in the plot of lines
    segments = []
    for miniM in maskUsed:
        segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                              [miniM.wl, 1.0]])]
    pltMaskU.set_segments(segments)
    segments = []
    for miniM in maskNot:
        segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                              [miniM.wl, 1.0]])]
    pltMaskN.set_segments(segments)
    #including the plot of lines highlighted for depth fitting
    maskFit = mask[(fitDepthFlags == 1) & (mask.iuse == 1)]
    segments = []
    for miniM in maskFit:
        segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                               [miniM.wl, 1.0]])]
    pltMaskF.set_segments(segments)
    return


class uiSelectFitDepth(uiRangeBase):
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update selected lines
        #update the data in the plot of lines
        fitDF = self.fitDepthFlags
        fitDF[(self.mask.wl >= minval) & (self.mask.wl <= maxval)] = 1
        maskFit = self.mask[(fitDF == 1) & (self.mask.iuse == 1)]
        segments = []
        for miniM in maskFit:
            segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                                   [miniM.wl, 1.0]])]
        self.pltMaskF.set_segments(segments)
        self.canvas.draw()


class uiUnselectFitDepth(uiRangeBase):
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update selected lines
        #update the data in the plot of lines
        fitDepthFlags = self.fitDepthFlags
        fitDepthFlags[(self.mask.wl >= minval) & (self.mask.wl <= maxval)] = 0
        maskFit = self.mask[(fitDepthFlags == 1) & (self.mask.iuse == 1)]
        segments = []
        for miniM in maskFit:
            segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                                   [miniM.wl, 1.0]])]
        self.pltMaskF.set_segments(segments)
        self.canvas.draw()


#Handle drawing a range selection on the screen, and getting the range selected
class rangeSelect:
    def __init__(self, canvas, parentUIRange):
        self.canvas = canvas
        self.canvasWidget = canvas.get_tk_widget()
        self.parentUIRange = parentUIRange
        self.bStartSelect = True
        self.lastrect = None
        
    def onClick(self, event):
        #Tkinter uses units from the top left, while matplotlib provides units from the bottom left.
        #So we need to flip matplotlilb's provided values (even in screen space)
        x0full, x1full = self.canvas.figure.axes[0].bbox.intervalx
        y0full, y1full = self.canvas.figure.axes[0].bbox.intervaly
        height = self.canvas.figure.bbox.height
        y0mod = height - y1full
        y1mod = height - y0full
        
        if self.bStartSelect:  #For the 1st click
            self.x0 = x0 = event.x
            self.y0 = y0 = event.y
            
            #drawn then update a rectangle as the mouse moves
            self.bStartSelect = False
            self.lastrect = self.canvasWidget.create_rectangle(x0, y0mod, x0, y1mod)
            self.oldBindMove = self.canvasWidget.bind('<Motion>')
            self.funcid = self.canvasWidget.bind('<Motion>', self.onMovePointer)
        else:
            self.x1 = event.x
            self.y1 = event.y

            #Get the cursor position in plotted data coordinates
            #(Testing this against matplotlib's toolbar it seems accurate to a pixel level.)
            axesdatax0, axesdatax1, axesdatay0, axesdatay1 = self.canvas.figure.axes[0].axis()
            axesxrange = x1full - x0full
            dataxrange = axesdatax1 - axesdatax0
            self.datax0 = (self.x0 - x0full)/(axesxrange)*(dataxrange) + axesdatax0
            self.datax1 = (self.x1 - x0full)/(axesxrange)*(dataxrange) + axesdatax0
            
            #Call the creating object's function for processing this range.
            #This is a bit recursive, so there may be a small possibility
            #to leak memory or odd behaviour?
            self.parentUIRange.selectedWl(min(self.datax0, self.datax1),
                                          max(self.datax0, self.datax1))
            self.deactivate()

    def cancel(self, *event):
            self.deactivate()
            
    def deactivate(self):
        #if there is an active selection region (not starting to select) undo it
        if not self.bStartSelect:
            self.bStartSelect = True
            self.canvasWidget.delete(self.lastrect)
            #Bind can take a string of tk binding info,
            #and bind returns such a string if called with just an event name,
            #so we can use that to restore a previous binding.
            self.funcid = self.canvasWidget.bind('<Motion>', self.oldBindMove)

    def onMovePointer(self, event):
        #Delete and redraw the rectangle when the mouse moves
        x1 = event.x
        y1 = event.y
        self.canvasWidget.delete(self.lastrect)
        #Get the full vertical range of the plot for the rectangle
        y0full, y1full = self.canvas.figure.axes[0].bbox.intervaly
        #Change from matplotlib to tk's y-axis direction
        height = self.canvas.figure.bbox.height
        y0mod = height - y1full
        y1mod = height - y0full
        self.lastrect = self.canvasWidget.create_rectangle(self.x0, y0mod, x1, y1mod, dash=[5,5], width=2)
        #self.lastrect = self.canvasWidget.create_rectangle(self.x0, y0mod, x1, y1mod, dash=[3,3])


class saveRanges:
    #mini manager for the button that saves the set of exclude ranges to file
    def __init__(self, fname, ranges):
        self.fname = fname
        self.ranges = ranges
    def saveToFile(self):
        #make an ExcludeMaskRegions object
        wlStarts = np.array([ran[0] for ran in self.ranges])
        wlEnds = np.array([ran[1] for ran in self.ranges])
        labels = np.array(['cleanMaskUI']*len(self.ranges),dtype=object)
        regions = maskTools.ExcludeMaskRegions(wlStarts, wlEnds, labels)
        #Save the set of exclude ranges to a file
        regions.save(self.fname)
        return


class fitDepths:
    #mini manager to run the fitting of line depths
    def __init__(self, mask, obs, lsdp, lsdProf, fitDepthFlags, root,
                 canvas, pltMaskU, pltMaskN, pltMaskF):
        self.mask = mask
        self.obs = obs
        self.lsdp = lsdp
        self.lsdProf = lsdProf
        self.fitDepthFlags = fitDepthFlags
        self.root = root
        self.canvas = canvas
        self.pltMaskU = pltMaskU
        self.pltMaskN = pltMaskN
        self.pltMaskF = pltMaskF
        
    def runFit(self):
        #Fit line depths using a linear least squares method
        #and then update the plots of lines from the mask.
        #This automatically removes degenerate or near-degenerate
        #lines (i.e. very near in wavelength) before fitting depths.

        #First set a 'wait' cursor
        oldCursor = self.root.cget('cursor')
        self.root.config(cursor='watch')
        self.root.update()

        #Get the lines that are flagged to be fit (and used)
        indUse1 = np.nonzero((self.fitDepthFlags == 1) & (self.mask.iuse == 1))
        useMask1 = self.mask[indUse1]
        if len(useMask1) > 0:
            #Get the reference LSD profile
            prof = self.lsdProf
            #Remove degenerate (or nearly) lines
            pixVel = prof.vel[1]-prof.vel[0]
            removePoorLines(useMask1, pixVel, fracPix = 3.0, sumDepths=False)
            self.mask.iuse[indUse1] = useMask1.iuse
            #and save the good lines in the mask for fitting
            indUse2 = np.nonzero(useMask1.iuse == 1)
            useMask2 = useMask1[indUse2]
            
            normDepth = self.lsdp.normDepth
            #Fit in two steps, and save the result back to the mask
            MD = buildMD(self.obs, useMask2, prof, normDepth)
            weights = linlsqDepths(MD, useMask2, self.obs)
            self.mask.depth[indUse1[0][indUse2]] = weights
            #update the plot
            updateLinePlots(self.mask, self.fitDepthFlags,
                            self.pltMaskU, self.pltMaskN, self.pltMaskF)
            self.canvas.draw()

        #Return the cursor to normal
        self.root.config(cursor=oldCursor)
        return

class undoDepths:
    #mini manager to undo the fitting of line depths
    def __init__(self, mask, fitDepthFlags, canvas,
                 pltMaskU, pltMaskN, pltMaskF):
        self.mask = mask
        self.fitDepthFlags = fitDepthFlags
        self.canvas = canvas
        self.pltMaskU = pltMaskU
        self.pltMaskN = pltMaskN
        self.pltMaskF = pltMaskF
        self.oldDepths = np.zeros_like(mask.depth)
        self.oldDepths[:] = mask.depth[:] #make a copy instead of a reference

    def undo(self):
        fUse = (self.fitDepthFlags == 1) & (self.mask.iuse == 1)
        self.mask.depth[fUse] = self.oldDepths[fUse]
        updateLinePlots(self.mask, self.fitDepthFlags,
                        self.pltMaskU, self.pltMaskN, self.pltMaskF)
        self.canvas.draw()
        return


#Make a window for modifying some LSD calculation parameters,
#includes associated functions for handling the UI
class modifyLSDpar:
    def __init__(self, parent, lsdp):
        self.parent = parent
        self.lsdp = lsdp
        self.active = False
        self.win = None

        self.varSaveMod = tk.IntVar()
        if self.lsdp.outModelName is None:
            self.varSaveMod.set(0)
        else:
            self.varSaveMod.set(1)
    
    def closeWindow(self, *event):
        #Update any input information
        self.setVelStart()
        self.setVelEnd()
        self.setVelPix()
        self.setNDepth()
        self.setNWav()
        self.setNLande()
        self.setModName()
        self.setProfName()
        
        self.active = False
        try:
            self.win.destroy()
        except:
            pass
        return
    
    def openWindow(self, *event):
        #If the window is already open do nothing
        if self.active:
            return
        self.active = True
        
        self.win = tk.Toplevel(self.parent)
        self.win.title("Set LSD parameters")
        #Associate this window with its parent,
        #and keep this window out of the window manager.
        self.win.transient(self.parent)
        #Send keyboard and mouse events only to this window,
        #mostly deactivating other windows.
        self.win.grab_set()
        #Set what happens when the window manager asks to close the window.
        self.win.protocol("WM_DELETE_WINDOW", self.closeWindow)

        wframe = ttk.Labelframe(self.win, text='Set parameters for LSD profile',
                                padding="10 10 10 10")
        wframe.pack(fill=tk.BOTH, expand=1, padx=4, pady=4)

        #Modify the velocity grid for the LSD profile
        labVelStart = ttk.Label(wframe, text='start vel', justify='right')
        self.txtVelStart = tk.StringVar()
        self.txtVelStart.set('{:.1f}'.format(self.lsdp.velStart))
        entVelStart = ttk.Entry(master=wframe, textvariable=self.txtVelStart,
                                width=8, validate='focusout',
                                validatecommand=self.setVelStart)
        entVelStart.bind('<Key-Return>', self.setVelStart)
        ToolTip(entVelStart, 'Set the starting velocity for the LSD profile (km/s)')
        labVelEnd = ttk.Label(wframe, text=' end vel', justify='right')
        self.txtVelEnd = tk.StringVar()
        self.txtVelEnd.set('{:.1f}'.format(self.lsdp.velEnd))
        entVelEnd = ttk.Entry(master=wframe, textvariable=self.txtVelEnd,
                              width=8, validate='focusout',
                              validatecommand=self.setVelEnd)
        entVelEnd.bind('<Key-Return>', self.setVelEnd)
        ToolTip(entVelEnd, 'Set the ending velocity for the LSD profile (km/s)')
        labVelPix = ttk.Label(wframe, text=' pixel size', justify='right')
        self.txtVelPix = tk.StringVar()
        self.txtVelPix.set('{:.2f}'.format(self.lsdp.velPixel))
        entVelPix = ttk.Entry(master=wframe, textvariable=self.txtVelPix,
                              width=8, validate='focusout',
                              validatecommand=self.setVelPix)
        entVelPix.bind('<Key-Return>', self.setVelPix)
        ToolTip(entVelPix, 'Set the pixel size of the LSD profile in km/s')

        #Modify LSD normalization values
        labNDepth = ttk.Label(wframe, text='norm. depth', justify='right')
        self.txtNDepth = tk.StringVar()
        self.txtNDepth.set('{:.2f}'.format(self.lsdp.normDepth))
        entNDepth = ttk.Entry(master=wframe, textvariable=self.txtNDepth,
                              width=8, validate='focusout',
                            validatecommand=self.setNDepth)
        entNDepth.bind('<Key-Return>', self.setNDepth)
        ToolTip(entNDepth, 'Set the normalizing depth for the LSD mask/profile')
        labNWav = ttk.Label(wframe, text=' norm. wavelength',
                            justify='right')
        self.txtNWav = tk.StringVar()
        self.txtNWav.set('{:.1f}'.format(self.lsdp.normWave))
        entNWav = ttk.Entry(master=wframe, textvariable=self.txtNWav,
                            width=8, validate='focusout',
                            validatecommand=self.setNWav)
        entNWav.bind('<Key-Return>', self.setNWav)
        ToolTip(entNWav, 'Set the normalizing wavelength for the LSD mask/profile (usually nm)')
        labNLande = ttk.Label(wframe, text=' norm. Lande',
                              justify='right')
        self.txtNLande = tk.StringVar()
        self.txtNLande.set('{:.2f}'.format(self.lsdp.normLande))
        entNLande = ttk.Entry(master=wframe, textvariable=self.txtNLande,
                              width=8, validate='focusout',
                            validatecommand=self.setNLande)
        entNLande.bind('<Key-Return>', self.setNLande)
        ToolTip(entNLande, 'Set the normalizing Landé factor for the LSD mask/profile')

        #remove closely spaced lines in the LSD calculation?
        self.varRemoveClose = tk.IntVar()
        self.varRemoveClose.set(int(self.lsdp.trimMask))
        checkRemoveClose = ttk.Checkbutton(wframe,
                                           text='remove closely\nspaced lines',
                                           variable=self.varRemoveClose,
                                           command=self.setRemoveClose)
        ToolTip(checkRemoveClose, "Remove very closely spaced lines in the mask when calculating the LSD profile")

        #Optionally save the model LSD spectrum
        checkSaveMod = ttk.Checkbutton(wframe,
                                       text='save model\nspectrum',
                                       variable=self.varSaveMod,
                                       command=self.setSaveMod)
        ToolTip(checkSaveMod, "Save a copy of the model LSD spectrum to a file?")
        labModName = ttk.Label(wframe, text='model\nfile name',
                               justify='right')
        self.txtModName = tk.StringVar()
        if self.lsdp.outModelName is None:
            self.txtModName.set('')
        else:
            self.txtModName.set('{:}'.format(self.lsdp.outModelName))
        entModName = ttk.Entry(master=wframe,
                               textvariable=self.txtModName, width=25,
                               validate='focusout',
                               validatecommand=self.setModName)
        entModName.bind('<Key-Return>', self.setModName)
        ToolTip(entModName, 'File name for the output model LSD spectrum')
        
        #Set the output LSD profile file name, and optionally plot it.
        labProfName = ttk.Label(wframe, text='profile\nfile name',
                                justify='right')
        self.txtProfName = tk.StringVar()
        self.txtProfName.set('{:}'.format(self.lsdp.outName))
        entProfName = ttk.Entry(master=wframe,
                                textvariable=self.txtProfName, width=25, 
                                validate='focusout',
                                validatecommand=self.setProfName)
        entProfName.bind('<Key-Return>', self.setProfName)
        ToolTip(entProfName, 'File name for the output LSD profile')

        self.varPltProf = tk.IntVar()
        self.varPltProf.set(int(self.lsdp.plotLSD))
        checkPltProf = ttk.Checkbutton(wframe, text='plot\nprofile',
                                       variable=self.varPltProf,
                                       command=self.setPlotProf)
        ToolTip(checkPltProf, "Plot the LSD profile in a new window when updated?")

        #Lay out the widgits in the window
        labVelStart.grid(row=0, column=0, sticky=tk.E, pady=4, padx=4)
        entVelStart.grid(row=0, column=1, sticky=tk.W, pady=4)
        labVelEnd.grid(row=1, column=0, sticky=tk.E, pady=4, padx=4)
        entVelEnd.grid(row=1, column=1, sticky=tk.W, pady=4)
        labVelPix.grid(row=2, column=0, sticky=tk.E, pady=4, padx=4)
        entVelPix.grid(row=2, column=1, sticky=tk.W, pady=4)

        labNDepth.grid(row=0, column=2, sticky=tk.E, pady=4, padx=4)
        entNDepth.grid(row=0, column=3, sticky=tk.W, pady=4)
        labNWav.grid(row=1, column=2, sticky=tk.E, pady=4, padx=4)
        entNWav.grid(row=1, column=3, sticky=tk.W, pady=4)
        labNLande.grid(row=2, column=2, sticky=tk.E, pady=4, padx=4)
        entNLande.grid(row=2, column=3, sticky=tk.W, pady=4)

        checkRemoveClose.grid(row=3, column=0, columnspan=4, padx=4, pady=4)

        checkSaveMod.grid(row=4, column=0, columnspan=1, pady=4)
        labModName.grid(row=4, column=1, sticky=tk.E, pady=4, padx=4)
        entModName.grid(row=4, column=2, columnspan=2, sticky=tk.W, pady=4)

        labProfName.grid(row=5, column=1, sticky=tk.E, pady=4, padx=4)
        entProfName.grid(row=5, column=2, columnspan=3, sticky=tk.W, pady=4)
        checkPltProf.grid(row=6, column=0, columnspan=4, pady=4, padx=4)

        wframe.columnconfigure(0, weight=1)
        wframe.columnconfigure(1, weight=1)
        wframe.columnconfigure(2, weight=1)
        wframe.columnconfigure(3, weight=1)
        #wframe.columnconfigure(4, weight=1)
        wframe.rowconfigure(0, weight=1)
        wframe.rowconfigure(1, weight=1)
        wframe.rowconfigure(2, weight=1)
        wframe.rowconfigure(3, weight=1)
        wframe.rowconfigure(4, weight=1)
        wframe.rowconfigure(5, weight=1)
        wframe.rowconfigure(6, weight=1)
        
        #Enter a local event loop, returning only when the given window
        #is destroyed.
        self.win.wait_window(self.win)

    def setVelStart(self, *event):
        try:  #Make sure this text is a number
            val = float(self.txtVelStart.get())
        except ValueError:
            self.txtVelStart.set('{:.1f}'.format(self.lsdp.velStart))
            return False
        if val <= -1000.:  #and make sure the value is sensible
            val = -999.0
            self.txtVelStart.set('{:.1f}'.format(val))
        elif (val > self.lsdp.velEnd - 10.):
            val = self.lsdp.velEnd - 10.
            self.txtVelStart.set('{:.1f}'.format(val))
        self.lsdp.velStart = val
        return True
    def setVelEnd(self, *event):
        try:  #Make sure this text is a number
            val = float(self.txtVelEnd.get())
        except ValueError:
            self.txtVelEnd.set('{:.1f}'.format(self.lsdp.velEnd))
            return False
        if val >= 1000.:  #and make sure the value is sensible
            val = 999.0
            self.txtVelEnd.set('{:.1f}'.format(val))
        elif (val < self.lsdp.velStart + 10.):
            val = self.lsdp.velStart + 10.
            self.txtVelEnd.set('{:.1f}'.format(val))
        self.lsdp.velEnd = val
        return True
    def setVelPix(self, *event):
        try:  #Make sure this text is a number
            val = float(self.txtVelPix.get())
        except ValueError:
            self.txtVelPix.set('{:.2f}'.format(self.lsdp.velPixel))
            return False
        if val <= 0:
            self.txtVelPix.set('{:.2f}'.format(self.lsdp.velPixel))
            return False
        elif val < 0.1:
            val = 0.1
            self.txtVelPix.set('{:.2f}'.format(val))
        elif val > (self.lsdp.velEnd - self.lsdp.velStart)/8.0:
            val = (self.lsdp.velEnd - self.lsdp.velStart)/8.0
            self.txtVelPix.set('{:.2f}'.format(val))
        self.lsdp.velPixel = val
        return True
    def setNDepth(self, *event):
        try:
            val = float(self.txtNDepth.get())
        except ValueError:
            self.txtNDepth.set('{:.2f}'.format(self.lsdp.normDepth))
            return False
        if val <= 0:
            self.txtNDepth.set('{:.2f}'.format(self.lsdp.normDepth))
            return False
        elif val > 1.0:
            val = 1.0
            self.txtNDepth.set('{:.2f}'.format(val))
        self.lsdp.normDepth = val
        return True
    def setNWav(self, *event):
        try:
            val = float(self.txtNWav.get())
        except ValueError:
            self.txtNWav.set('{:.1f}'.format(self.lsdp.normWave))
            return False
        if val <= 0:
            self.txtNWav.set('{:.1f}'.format(self.lsdp.normWave))
            return False
        elif val > 100000.0:
            val = 100000.0
            self.txtNWav.set('{:.1f}'.format(val))
        self.lsdp.normWave = val
        return True
    def setNLande(self, *event):
        try:
            val = float(self.txtNLande.get())
        except ValueError:
            self.txtNLande.set('{:.2f}'.format(self.lsdp.normLande))
            return False
        if val <= 0:
            self.txtNLande.set('{:.2f}'.format(self.lsdp.normLande))
            return False
        elif val > 10.0:
            val = 10.0
            self.txtNLande.set('{:.2f}'.format(val))
        self.lsdp.normLande = val
        return True

    def setRemoveClose(self, *event):
        val = self.varRemoveClose.get()
        self.lsdp.trimMask = bool(val)
        return
    
    def setSaveMod(self, *event):
        val = self.varSaveMod.get()
        if val == 0:
            self.lsdp.outModelName = None
        else:
            self.lsdp.outModelName = self.txtModName.get()
        return
    def setModName(self, *event):
        flag = self.varSaveMod.get()
        if flag != 0:
            self.lsdp.outModelName = self.txtModName.get()
        else:
            self.lsdp.outModelName = None
        return
    def setProfName(self, *event):
        val = self.txtProfName.get()
        if val != '':
            self.lsdp.outName = val
        return
    def setPlotProf(self, *event):
        val = self.varPltProf.get()
        if val == 0:
            self.lsdp.plotLSD = False
        else:
            self.lsdp.plotLSD = True
        return


class updateLSD:
    def __init__(self, canvas, root, mask, lsdp, lsdProf, pltModelI):
        self.canvas = canvas
        self.root = root
        self.mask = mask
        self.lsdp = lsdp
        self.lsdProf = lsdProf
        self.pltModelI = pltModelI
    def rerunLSD(self):
        #Recalculate the LSD profile and update the plot
        #First set a 'wait' cursor
        oldCursor = self.root.cget('cursor')
        self.root.config(cursor='watch')
        self.root.update()
        #Run LSD
        lsdp = self.lsdp
        self.mask.save(lsdp.mask)
        lsdProf, modelSpec = profileLSD.run_lsdpy(obs=lsdp.obs, mask=lsdp.mask,
            outLSDName=lsdp.outName, velStart=lsdp.velStart, velEnd=lsdp.velEnd,
            velPixel=lsdp.velPixel, normDepth=lsdp.normDepth,
            normLande=lsdp.normLande, normWave=lsdp.normWave,
            removeContPol=lsdp.removeContPol, trimMask=lsdp.trimMask,
            sigmaClipIter=lsdp.sigmaClipIter, sigmaClip=lsdp.sigmaClip,
            outModelName=lsdp.outModelName,
            plotLSD=lsdp.plotLSD, outPlotLSDName=lsdp.outPlotLSDName)
        #update the current LSD object, rather than replacing it
        #(since fitDepths needs to retain a reference to the same object)
        #self.lsdProf[:] = lsdProf[:]
        self.lsdProf.vel = lsdProf.vel
        self.lsdProf.specI = lsdProf.specI
        self.lsdProf.specSigI = lsdProf.specSigI
        self.lsdProf.specV = lsdProf.specV
        self.lsdProf.specSigV = lsdProf.specSigV
        self.lsdProf.specN1 = lsdProf.specN1
        self.lsdProf.specSigN1 = lsdProf.specSigN1
        self.lsdProf.specN2 = lsdProf.specN2
        self.lsdProf.specSigN2 = lsdProf.specSigN2
        self.lsdProf.npix = lsdProf.npix
        self.lsdProf.numParam = lsdProf.numParam
        self.lsdProf.header = lsdProf.header
        #Return the cursor to normal
        self.root.config(cursor=oldCursor)
        
        #re-draw points used in the fit
        for dline in self.pltModelI:
            dline.set_data(modelSpec.wl, modelSpec.specI)
        self.canvas.draw()
        return


def buildMD(obs, mask, prof, normDepth):
    #Set up an matrix that that can be multiplied with a vector of line weights
    #to produce a mode LSD spectrum, array dimensions are nObsPts x nMaskLines
    #This has columns for each line in the mask that each essentially 
    #contain the LSD profile, with the profile shifted to be centred 
    #on the line's wavelength, and interpolated onto the wavelengths of 
    #the observation (the row pixels).  The the dot product of this
    #matrix with the vector of line mask weights then scales each copy of 
    #the LSD profile by that line's weight, and then sums over all lines.
    #Thus this dot product the LSD model spectrum.
    cvel = 2.99792458e5 #c in km/s

    #calculate wavelengths for the profile at each line in the mask here, since it is reusable
    wlProfA = np.outer(prof.vel/cvel, mask.wl) + np.tile(mask.wl, (prof.npix,1))  #(prof.npix, mask.wl.shape)

    #We can use a sparse "List of Lists" matrix to efficiently build the desired
    #matrix, but we have to start with the transpose of what we want for
    #efficient access when adding (non-zero) entries to the lil_matrix.
    Dlil = scipy.sparse.lil_matrix((mask.wl.shape[0], obs.wl.shape[0]))

    for l in range(mask.wl.shape[0]):
        #Get observation points in range of this line in the mask
        iObsRange = np.where( (wlProfA[0,l] < obs.wl[:]) & (wlProfA[-1,l] > obs.wl[:]) )
        #Interpolate the LSD profile for this line onto the observation
        interpProfI = np.interp(obs.wl[iObsRange], wlProfA[:,l], 1.-prof.specI)

        Dlil[l, iObsRange] = interpProfI / normDepth

    #Convert the matrix into a more efficient sparse form (csr_matrix) for matrix math
    Dout = Dlil.tocsr().T
    
    ##Output the model spectrum, for testing purposes only.
    #weightI = mask.depth/normDepth
    #specT = 1.-Dout.dot(weightI)
    #fOutPre = open('outSpecBeforeTweak.dat', 'w')
    #for i in range(obs.wl.shape[0]):
    #    fOutPre.write('{:} {:}\n'.format(obs.wl[i], specT[i]))
    #fOutPre.close()    
    return Dout

def linlsqDepths(maskMD, mask, obs):
    #Perform the linear least squares fit to the observation
    #to find the optimal line depths for the LSD mask.
    #This follow the same logic as the LSD profile calculation in LSDpy,
    #or most other linear least squares calculations. Here the model is 
    #Y = D.W, for the vector of line depths W, and 2D matrix M.  For that 
    #model with observations Yo, with S a diagonal matrix of 1/sigma
    # chi^2 = (Yo - D.W)^T.(S^2).(Yo - D.W)
    #Taking the derivative of chi^2 wrt W and setting it equal to 0 
    #to find the minimum in chi^2 gives
    # 0 = (-D)^T.(S^2).(Yo - D.W)*2
    # W = (D^T.(S^2).D)^(-1).D^T.(S^2).Yo
    #so we mostly just need to invert the matrix (D^T.(S^2).D)
    #and also build D^T.(S^2).Yo

    #We are working in shifted units with the continuum at 0 and lines +ve
    obsI = 1. - obs.specI
    
    #Use some sparse matrices here to be efficient
    #The dot products for alpha and beta can be slower than matrix inversion
    #Use similar choices to LSDpy for sparse routines
    sparseS2 = scipy.sparse.diags(obs.specSig**(-2), offsets=0)
    MD = scipy.sparse.csr_matrix(maskMD) #probably already a sparse csr
    
    beta = MD.T.dot(sparseS2.dot(obsI))
    alpha = MD.T.dot(sparseS2.dot(MD))

    covar = np.linalg.inv(alpha.toarray())
    lineWeights = np.dot(covar, beta)
    #lineWeightErrs = np.sqrt(np.diag(covar)) #if errors are needed
    return lineWeights


def removePoorLines(mask, pixVel, fracPix = 1.0, sumDepths=True):
    #Remove nearly degenerate lines from the mask.
    #Reject lines separated by less than fracPix of an LSD (velocity) pixel.
    #(Based on LSDpy)
    cvel = 2.99792458e5 #c in km/s
    depthCutoff = 0.6
    nTrimmed = 0
    for l in range(1,mask.wl.shape[0]):
        #This loop is relatively inefficient but handles unsorted line masks
        #and lines with multiple bad blends.
        if mask.iuse[l] == 1:
            deltas = np.abs(mask.wl[l] - mask.wl)/mask.wl[l]*cvel
            iClose = np.nonzero((deltas < pixVel*fracPix) & (mask.iuse == 1))[0]
            if iClose.shape[0] > 1:
                #If other lines are too close to the current line
                mask.iuse[iClose] = 0
                deepestLine = np.argmax(mask.depth[iClose])
                mask.iuse[iClose[deepestLine]] = 1
                
                #If we want to sum line depths, limit the maximum depth lines 
                #can sum to, as a very rough approximation for saturation.
                if sumDepths:
                    summedDepth = np.sum(mask.depth[iClose])
                    if summedDepth < depthCutoff:
                        mask.depth[iClose[deepestLine]] = summedDepth
                    else:
                        mask.depth[iClose[deepestLine]] = max(depthCutoff,
                                                mask.depth[iClose[deepestLine]])
                nTrimmed += np.count_nonzero(iClose) - 1
    if nTrimmed > 0:
        print('Modified line mask, removed {:d} too closely spaced lines'.format(nTrimmed))

    ##If one wanted a shorter mask
    #mask = mask.prune()
    return


class lsdParams:
    def __init__(self, obs, mask, outName='',
                 velStart=-200., velEnd=200., velPixel=None,
                 normDepth=0.2, normLande=1.2, normWave=500.,
                 removeContPol=True, trimMask=False, 
                 sigmaClipIter=0, sigmaClip=500.,
                 outModelName=None,
                 plotLSD=False, outPlotLSDName=None):
        """
        A simple data structure to hold input parameters for LSD calculations,
        and to set some sensible defaults.
        """
        self.obs = obs
        self.mask = mask
        self.outName = outName
        self.velStart = velStart
        self.velEnd = velEnd
        self.normDepth = normDepth
        self.normLande = normLande
        self.normWave = normWave
        self.removeContPol = removeContPol
        self.trimMask = trimMask
        self.sigmaClipIter = sigmaClipIter
        self.sigmaClip = sigmaClip
        self.outModelName = outModelName
        self.plotLSD = plotLSD
        self.outPlotLSDName = outPlotLSDName
        
        cvel = 2.99792458e5 #c in km/s
        if velPixel == None:
            #Get average observed wavelength spacing
            obsSp = obsSpec.read_spectrum(self.obs)
            wlStep = obsSp.wl[1:] - obsSp.wl[:-1]
            medianVelPix = np.median(wlStep/obsSp.wl[:-1]*cvel)
            self.velPixel = round(medianVelPix, ndigits=2)
        else:
            self.velPixel = velPixel
