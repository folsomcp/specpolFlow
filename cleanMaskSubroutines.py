#Additional subroutines used inside cleanMaskUI.py
#Not realy meant to be stand alone utilities.

import numpy as np
import tkinter as tk
import tkinter.ttk as ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg)
#Try to import the default matplotlib toolbar
try: #for matplotlib 3.x
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
    mplToolbar=3
except ImportError:  #If that fails this probably isn't matplotlib 3.x
    try: #for matplotlib 2.x
        from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
        mplToolbar=2
    except ImportError: #not sure what failed, try just not using this extra toolbar
        mplToolbar=0
try:
    import specpolFlow.iolsd as iolsd
except ModuleNotFoundError:
    import iolsd


def makeWin(fig, ax, mask, obsName, maskName, pltMaskU, pltMaskN,
            pltModelI, excludeRanges, excludeFileName):
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
    if mplToolbar > 0:
        if mplToolbar == 2:
            toolbarMpl = NavigationToolbar2TkAgg(canvas, root)
        if mplToolbar == 3:
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
    #Set regions to be included in the mask
    incRangeM = uiIncludeRange(canvas, mask, pltMaskU, pltMaskN, excludeRanges)
    butIncRange = ttk.Button(master=tools, text='include\nlines',
                             command=incRangeM.runSpanSelect)
    ToolTip(butIncRange, 'Include selected lines in the mask')
    butIncRange.grid(row=0, column=7, sticky=tk.E, padx=2, pady=2)
    #Set regions to be excluded in the mask
    excRangeM = uiExcludeRange(canvas, mask, pltMaskU, pltMaskN, excludeRanges)
    butExcRange = ttk.Button(master=tools, text='exclude\nlines',
                             command=excRangeM.runSpanSelect)
    ToolTip(butExcRange, 'Exclude selected lines from the mask')
    butExcRange.grid(row=0, column=8, sticky=tk.E, padx=2, pady=2)
    #Link the include and exclude buttons so they can turn eachother off
    incRangeM.linkButton(butIncRange, excRangeM)
    excRangeM.linkButton(butExcRange, incRangeM)
    #Save the exclude ranges to a text file
    saveRangesM = saveRanges(excludeFileName, excludeRanges)
    butSaveRanges = ttk.Button(master=tools, text='save\nranges',
                               command=saveRangesM.saveToFile)
    ToolTip(butSaveRanges, 'Save the excluded wavelength ranges to the {:} file'.format(excludeFileName))
    butSaveRanges.grid(row=0, column=9, sticky=tk.E, padx=2, pady=2)
    
    #Update the LSD calculation and the plotted spectrum
    updateLSDM = updateLSD(canvas, root, mask, obsName, maskName, pltModelI)
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


class saveRanges:
    #mini manager for the button that saves the set of exclude ranges to file
    def __init__(self, fname, ranges):
        self.fname = fname
        self.ranges = ranges
    def saveToFile(self):
        #Save the set of exclude ranges to a file
        fout = open(self.fname,'w')
        fout.write('#exclude wavelength ranges\n')
        for ran in self.ranges:
            fout.write('{:9.3f} {:9.3f}\n'.format(ran[0], ran[1]))
        fout.close()
        return


#A fairly simple, fairly general tooltip
class ToolTip(object):
    """
    Create a tooltip for a given widget, using Tkinter
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
            #and correct for matplotlib using y in the opposide direction from tk
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
            #set the plotted range into the selected range in data coordiantes
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


class uiIncludeRange:
    def __init__(self, canvas, mask, pltMaskU, pltMaskN, excludeRanges):
        self.canvas = canvas
        self.canvasWidget = canvas.get_tk_widget()
        self.mask = mask
        self.pltMaskU = pltMaskU
        self.pltMaskN = pltMaskN
        self.excludeRanges = excludeRanges
        self.active = False
        self.rangeSelect = rangeSelect(canvas, self)

        #Set up a style for the 'actively selecting a range' button state
        #Seems like I need to set the base style changes first,
        #then use map to set state specific differences,
        #since map doesn't provide easy access to being in
        #not any of the special states (only not a specific state)
        self.style = ttk.Style()
        self.style.configure('ActRange.TButton', relief='sunken',
                        background='#fbfbfb')
        self.style.map('ActRange.TButton',
                  relief=[('active','sunken')],
                  background=[('pressed','#eeeeee'),('active','#fbfbfb')])
        #note: order matters for map here, 
        #when the button is 'active' (mouse over) and also being pressed
        #The scope of the style names seems to be pretty wide...

    def linkButton(self, button, otherRange):
        self.button = button
        self.otherRange = otherRange
    def runSpanSelect(self):
        #First turn off the other button (only include or exclude not both!)
        bOtherActive = self.otherRange.active
        if (bOtherActive == True): 
            self.otherRange.deactivate()
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
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update the mask and 'exclude regions'
        self.mask.iuse[(self.mask.wl >= minval) & (self.mask.wl <= maxval)] = 1
        maskUsed = self.mask[self.mask.iuse == 1]
        maskNot = self.mask[self.mask.iuse == 0]
        
        self.excludeRanges[:] = removeRange(self.excludeRanges, [minval,maxval])
        
        #update the data in the plot for each spectal order
        #segments should be a list of 2x2 arrays
        segments = []
        for miniM in maskUsed:
            segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                                  [miniM.wl, 1.0]])]
        self.pltMaskU.set_segments(segments)
        segments = []
        for miniM in maskNot:
            segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                                  [miniM.wl, 1.0]])]
        self.pltMaskN.set_segments(segments)
        self.canvas.draw()

        
class uiExcludeRange:
    def __init__(self, canvas, mask, pltMaskU, pltMaskN, excludeRanges):
        self.canvas = canvas
        self.canvasWidget = canvas.get_tk_widget()
        self.mask = mask
        self.pltMaskU = pltMaskU
        self.pltMaskN = pltMaskN
        self.excludeRanges = excludeRanges
        self.active = False
        self.rangeSelect = rangeSelect(canvas, self)
        
    def linkButton(self, button, otherRange):
        self.button = button
        self.otherRange = otherRange
    def runSpanSelect(self):
        #First turn off the other button (only include or exclude not both!)
        bOtherActive = self.otherRange.active
        if (bOtherActive == True):
            self.otherRange.deactivate()
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
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update the mask and 'exclude regions'
        self.mask.iuse[(self.mask.wl >= minval) & (self.mask.wl <= maxval)] = 0
        maskUsed = self.mask[self.mask.iuse == 1]
        maskNot = self.mask[self.mask.iuse == 0]

        self.excludeRanges[:] = combineRanges(self.excludeRanges,
                                              [[minval,maxval]])
        
        #update the data in the plot for each spectal order
        #segments should be a list of 2x2 arrays
        segments = []
        for miniM in maskUsed:
            segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                                  [miniM.wl, 1.0]])]
        self.pltMaskU.set_segments(segments)
        segments = []
        for miniM in maskNot:
            segments += [np.array([[miniM.wl, 1.0-miniM.depth],
                                  [miniM.wl, 1.0]])]
        self.pltMaskN.set_segments(segments)
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
            #to leak memory or odd behavior?
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
            #and bind retuns such a string if called with just an event name,
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


class updateLSD:
    def __init__(self, canvas, root, mask, obsName, maskName, pltModelI):
        self.canvas = canvas
        self.root = root
        self.mask = mask
        self.obsName = obsName
        self.maskName = maskName
        self.pltModelI = pltModelI
    def rerunLSD(self):
        #Recalculate the LSD profile and update the plot
        #First set a 'wait' cursor
        oldCursor = self.root.cget('cursor')
        self.root.config(cursor='watch')
        self.root.update()
        #Run LSD
        self.mask.save(self.maskName)
        lsdProf, modelSpec = iolsd.run_lsdpy(obs=self.obsName, mask=self.maskName, fLSDPlotImg=0)
        #Return the cursor to normal
        self.root.config(cursor=oldCursor)
        
        #re-draw points used in the fit
        for dline in self.pltModelI:
            dline.set_data(modelSpec.wl, modelSpec.specI)
        self.canvas.draw()
        return
