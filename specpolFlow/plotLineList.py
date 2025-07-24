import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def plot_LineList(llist, depthCut=0.0, maxLabels=None,
                  scaleDepths=1.0, cont=1.01, rise=0.05,
                  nrows=2, padding=4.0, vpadding=6.0,
                  romanIon=False, avoidOverlaps=True, dynamicUpdate=True,
                  linewidth=1.0, linecolor='grey', linestyle='-',
                  fontsize=8, rotation='horizontal',
                  bindKeys=True, ax=None, **kwargs):
    '''
    Plot a line list, with tick marks at the positions of lines and text labels
    for the species of line.

    This will either generate new matplotlib figure and axes objects with
    the plot, or if an existing axes is passed to the function (with the ax
    keyword) then the line list will be plotted in that axes.  
    
    By default, this function will try to adjust the positions of the labels
    so that they don't overlap.

    Labels plotted can be limited to only being deeper than a threshold depth
    with depthCut.  The maximum number of labels drawn can be limited with
    maxLabels, in this case tick marks are drawn for all lines, only the
    number text labels is limited, which can be useful for efficiency.
    
    :param llist: the line list to plot, as a LineList object
    :param depthCut: only lines with a depth value bigger than this will be plotted.
    :param maxLabels: the maximum number of line labels that will be drawn at
                      one time. Only labels for the maxLabels deepest lines
                      depths will be drawn, although tick marks for all lines
                      will be drawn.
                      Drawing labels is relatively slow, and drawing a large
                      number can make this routine run slowly (as well as being 
                      hard to read!).  Defaults to plotting all labels.
    :param scaleDepths: scale the depths of the tick marks by this value. Tick 
                        depths are proportional to the line depth parameter.
                        (Tick marks extend down to 1.0 - llist.depth*scaleDepths)
    :param cont: draw tick lines from the depth up to this value
                 (typically set this at or just above the continuum level).
    :param rise: place the labels this much above the cont level.
    :param nrows: the number of rows of labels to be drawn
    :param padding: spacing to leave between labels (in pixels)
    :param vpadding: vertical spacing to leave between rows of labels (in pixels),
                     only used when avoidOverlaps is True.
    :param romanIon: if True, convert the numbers in the ion strings to use
                     roman numerals. Roman numerals for ionization stage 
                     are typically better for publication quality figures.
    :param avoidOverlaps: if True shift, the positions of labels to avoid 
                          having the text overlap.  If too many labels are
                          plotted, overlaps may still occur (try changing
                          maxLabels or depthCut).  Defaults to True.
    :param dynamicUpdate: update the label positions of every time the figure changes.
                          This keeps labels from overlapping or being cut off
                          if the x-axis changes (e.g. when using ax.set_xlim() or
                          through panning or zooming in an interactive window).
                          (This connects to the Axes 'xlim_changed' event and the
                          Figure 'resize_event', and stores some additional data
                          in the Figure object as fig.dynamicLineList.)
                          Defaults to True.
    :param linewidth: matplotlib line width for the tick marks
    :param linecolor: matplotlib color for the tick marks
    :param linestyle: matplotlib line style for the tick marks
    :param fontsize: font size for the labels
    :param rotation: orientation of the line labels, can be 'horizontal',
                     'vertical', or a float with an angle in degrees.
    :param bindKeys: bind some keys (using matplotlib) to interactively
                     adjust some plot properties.  Set this False to prevent
                     this behaviour.  Only active if dynamicUpdate is True.
                     Keys used are:
                     arrow keys -- pan left, right, up, and down;
                     i -- zoom in; o -- zoom out;
                     z -- activate matplotlib's zoom tool;
                     a -- autoscale the zoom to show all data;
                     A -- autoscale the y-axis only;
    :param ax: the matplotlib axes object to plot the lines in. If this is 
               None, then a new figure and axes will be generated.
    :param **kwargs: any remaining keyword arguments are passed to the matplotlib 
                     ax.text() function when creating the labels for lines.
    :return: a matplotlib figure object, and an axes object,
             containing the plot.
    '''

    llist_all = llist[np.argsort(llist.wl)]

    #optionally convert the ion strings to use roman numerals
    if romanIon:
        fion = _fancy_ion_string(llist_all.ion)
        llist_all.ion = fion

    if maxLabels is None:
        maxLabels = llist.nLines
    
    riser = cont + rise #place labels at this y level
    rowSpacing = rise #spacing between rows of labels in data coordinates
    #(for fast mode, when the size of the labels in display coordinates isn't known)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12,6)) #layout='constrained'
    else:
        fig = ax.get_figure()

    # If this is a new empty Axes object, set the the x-range,
    # this is useful for placing non-overlapping labels later.
    if ax.has_data():
        wlMin, wlMax = ax.get_xlim()
    else:
        wlMin = np.min(llist.wl)
        wlMax = np.max(llist.wl)
        wlRange = wlMax - wlMin
        ax.set_xlim(wlMin - 0.05*wlRange, wlMax + 0.05*wlRange, auto=True)
        ax.set_ylim(0.0, 1.1, auto=True)

    # trim the shown lines to be deep enough
    show_lines = llist_all.depth > depthCut

    # trim the labels shown to be in the wavelength range,
    # and limit the total number shown
    show_labels, ind_labels_in_lines = _get_visible_labels(
        wlMin, wlMax, show_lines, llist_all, maxLabels)

    #generate a list of line labels first (hide them all now, show some later)
    llabels = [] # list of text objects
    for i, line in enumerate(llist_all):
        if dynamicUpdate or show_labels[i]:
            txt = ax.text(line.wl, riser, line.ion, 
                          ha='center', va='bottom', rotation=rotation,
                          fontsize=fontsize,
                          visible=False, clip_on=True, **kwargs)
            llabels += [txt]

    # flag each n-th label that is shown as being in the n-th row
    nrowText = np.arange(np.sum(show_labels)) % nrows
    labelFlx = riser + nrowText*rowSpacing
    if dynamicUpdate:  #leave hidden labels in row 0
        _nrowText = np.zeros(len(llabels), dtype=int)
        _nrowText[show_labels] = nrowText
        nrowText = _nrowText
    labelWls = llist_all.wl[show_labels]
    # shift those labels into their rows, and set as visable
    llabels2d = [[] for nrow in range(nrows)] # list of rows with lists of text
    j = 0
    for i, label in enumerate(llabels):
        if show_labels[i] or not dynamicUpdate:
            label.set_visible(True)
            label.set_y(labelFlx[j])
            llabels2d[nrowText[i]] += [label]
            j += 1

    #If not in 'fast' mode, set the line label positions to avoid overlapping
    if avoidOverlaps:
        labelWls, labelFlx = _set_label_pos(fig, ax, llabels2d,
                                    llist_all[show_labels], padding, vpadding)

    # Generate tick marks and connectors to labels,
    # use 3 points (2 line segments) to make a vertical tick & riser connector
    # the final array of points should have shape (nLines, 3, 2)
    ##wlpts = np.tile(llist_show.wl, (3,1))
    ##flxpts = np.tile(1.0 - llist_show.depth*scaleDepths, (3,1))
    wlpts = np.tile(llist_all.wl[show_lines], (3,1))
    flxpts = np.tile(1.0 - llist_all.depth[show_lines]*scaleDepths, (3,1))
    linePts = np.stack([wlpts.T,flxpts.T], axis=2)
    # adjust the points for the riser connectors, where there are labels
    linePts[:, 0:2, 1] = cont
    linePts[ind_labels_in_lines, 0, 1] = labelFlx
    linePts[ind_labels_in_lines, 0, 0] = labelWls

    # plot the tick marks as a matplotlib LineCollection
    linecollection = LineCollection(linePts, colors=linecolor,
                                    linewidths=linewidth, linestyles=linestyle)
    lcPlot = ax.add_collection(linecollection)

    if dynamicUpdate:
        ulp = _updateLinesPlot(llist_all, avoidOverlaps, depthCut, maxLabels, 
                              scaleDepths, cont, rise, nrows, padding, vpadding,
                              lcPlot, llabels, show_lines, show_labels)
        ulp.connectCallbacks(fig, ax)
        
        #Keyboard shortcuts
        if bindKeys:
            # numpy defaults to mapping arrow keys to previous and next view
            # (also it maps zoom to o)
            # remap keys for these functions so we can use them to pan
            plt.rcParams['keymap.back'] = [u',', u'<', u'backspace']
            plt.rcParams['keymap.forward'] = [u'.', u'>']
            plt.rcParams['keymap.zoom'] = [u'z']
            # catch key press events and call our function
            fig.canvas.mpl_connect('key_press_event', ulp.key_press_handler)

        # Save an object with callback functions and references to data for
        # dynamic updating of the plot.  Save it to the Figure object.
        # (to preserve a reference to the object after this function ends)
        if hasattr(fig, 'dynamicLineList'):
            fig.dynamicLineList += [ulp]
        else:
            fig.dynamicLineList = [ulp]
    
    return fig, ax


class _updateLinesPlot:
    '''
    Save data and provide functions for dynamically updating the line list plot.
    '''
    def __init__(self, llist_all, avoidOverlaps, depthCut, maxLabels, 
                 scaleDepths, cont, rise, nrows, padding, vpadding,
                 linecollection, llabels, show_lines, show_labels):
        self.llist_all = llist_all
        # format parameters
        self.avoidOverlaps = avoidOverlaps
        self.deptCut = depthCut
        self.maxLabels = maxLabels
        self.scaleDepths = scaleDepths
        self.cont = cont
        self.rise = rise
        self.nrows = nrows
        self.padding = padding
        self.vpadding = vpadding
        # derived format parameters
        self.riser =  cont + rise
        self.rowSpacing = rise
        # data
        self.linecollection = linecollection
        self.llabels = llabels
        self.show_lines = show_lines
        self.show_labels = show_labels
    def connectCallbacks(self, fig, ax):
        self.fig = fig
        self.ax = ax
        self.ax.callbacks.connect('xlim_changed', self.on_lims_change)
        ##self.ax.callbacks.connect('ylim_changed', self.on_lims_change)
        self.fig.canvas.mpl_connect('resize_event', self.on_fig_resize)
        #self.fig.canvas.mpl_connect('draw_event', self.on_fig_draw)
        
    def on_lims_change(self, ax):
        self.update_labels(redraw=False)
        return
    def on_fig_resize(self, event):
        self.update_labels(redraw=True)
        return
    def on_fig_draw(self, event):
        #self.update_labels(redraw=False)
        return

    def update_labels(self, redraw):
        '''
        Update which labels are visable, and update their positions,
        for the new window and x axis ranges.
        '''
        wlMin, wlMax = self.ax.get_xlim()

        # first reset all labels to not visible
        for i in np.nonzero(self.show_labels)[0]:
            self.llabels[i].set_visible(False)
            self.llabels[i].set_x(self.llist_all.wl[i])

        # get the labels that should be visible
        show_labels, ind_labels_in_lines = _get_visible_labels(
            wlMin, wlMax, self.show_lines, self.llist_all, self.maxLabels)

        # flag each n-th label that is shown as being in the n-th row
        nrowText = np.zeros(len(self.llabels), dtype=int)
        nrowText[show_labels] = np.arange(np.sum(show_labels)) % self.nrows
        labelFlx = self.riser + nrowText[show_labels]*self.rowSpacing
        labelWls = self.llist_all.wl[show_labels]
        # and shift those labels into their rows
        llabels2d = [[] for nrow in range(self.nrows)] # list of rows
        j = 0
        for i, label in enumerate(self.llabels):
            if show_labels[i]:
                label.set_visible(True)
                label.set_y(labelFlx[j])
                llabels2d[nrowText[i]] += [label]
                j += 1

        if self.avoidOverlaps:
            labelWls, labelFlx = _set_label_pos(self.fig, self.ax, llabels2d,
                                        self.llist_all[show_labels], self.padding,
                                        self.vpadding, redraw=redraw)
        # save the labels being shown for later reference
        self.show_labels = show_labels

        # Generate tick marks and connectors to labels
        wlpts = np.tile(self.llist_all.wl[self.show_lines], (3,1))
        flxpts = np.tile(1.0 - self.llist_all.depth[self.show_lines]*self.scaleDepths, (3,1))
        linePts = np.stack([wlpts.T,flxpts.T], axis=2)
        # adjust the points for the riser connectors, where there are labels
        linePts[:, 0:2, 1] = self.cont
        linePts[ind_labels_in_lines, 0, 1] = labelFlx
        linePts[ind_labels_in_lines, 0, 0] = labelWls

        # update the LineCollection object data
        self.linecollection.set_segments(linePts)
        return

    def key_press_handler(self, event):
        '''
        receive a key press event and call the appropriate function for that key
        '''
        if event.key == u'left':
            self.pan_left(event)
        if event.key == u'right':
            self.pan_right(event)
        if event.key == u'up':
            self.pan_up(event)
        if event.key == u'down':
            self.pan_down(event)
        if event.key == u'i':
            self.zoom_in(event)
        if event.key == u'o':
            self.zoom_out(event)
        if event.key == u'a':
            self.autoScale(event)
        if event.key == u'A':
            self.autoScaleY(event)
        return
    
    def pan_left(self, event):
        fracPan = 0.2
        xmin, xmax = self.ax.get_xlim()
        xspan = xmax - xmin
        self.ax.set_xlim(xmin - xspan*fracPan, xmax - xspan*fracPan)
        self.fig.canvas.draw()
        return
    def pan_right(self, event):
        fracPan = 0.2
        xmin, xmax = self.ax.get_xlim()
        xspan = xmax - xmin
        self.ax.set_xlim(xmin + xspan*fracPan, xmax + xspan*fracPan)
        self.fig.canvas.draw()
        return
    def pan_up(self, event):
        fracPan = 0.2
        ymin, ymax = self.ax.get_ylim()
        yspan = ymax - ymin
        self.ax.set_ylim(ymin + yspan*fracPan, ymax + yspan*fracPan)
        self.fig.canvas.draw()
        return
    def pan_down(self, event):
        fracPan = 0.2
        ymin, ymax = self.ax.get_ylim()
        yspan = ymax - ymin
        self.ax.set_ylim(ymin - yspan*fracPan, ymax - yspan*fracPan)
        self.fig.canvas.draw()
        return
    def zoom_in(self, event):
        fracZoom = 0.1
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        xspan = xmax-xmin
        yspan = ymax-ymin
        modFracZoom = fracZoom/(1.+2.*fracZoom) #modify so zoom in reverses zoom out exactly
        self.ax.set_xlim(xmin+xspan*modFracZoom, xmax-xspan*modFracZoom)
        self.ax.set_ylim(ymin+yspan*modFracZoom, ymax-yspan*modFracZoom)
        self.fig.canvas.draw()
        return
    def zoom_out(self, event):
        fracZoom = 0.1
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        xspan = xmax-xmin
        yspan = ymax-ymin
        self.ax.set_xlim(xmin-xspan*fracZoom, xmax+xspan*fracZoom)
        self.ax.set_ylim(ymin-yspan*fracZoom, ymax+yspan*fracZoom)
        self.fig.canvas.draw()
    def autoScale(self, event):
        fracSpace = 0.05
        self.ax.autoscale(enable=True, axis='both')
        self.fig.canvas.draw()
        return
    def autoScaleY(self, event):
        #self.ax.autoscale(enable=True, axis='y')
        fracSpace = 0.05
        xmin, xmax = self.ax.get_xlim()
        axLines = self.ax.get_lines()
        yminList = []
        ymaxList = []
        for lineData in axLines:
            lineDataX = lineData.get_xdata()
            lineDataY = lineData.get_ydata()
            lineDataY = lineDataY[(lineDataX > xmin) & (lineDataX < xmax)]
            if lineDataY.size > 0:
                yminList += [np.min(lineDataY)]
                ymaxList += [np.max(lineDataY)]

        if len(yminList) > 0 and len(ymaxList) > 0:
            ymin = min(yminList)
            ymax = max(ymaxList+[self.cont + self.rise*self.nrows])
        else:
            ymin = 0.0
            ymax = 1.1
        yspan = ymax - ymin
        self.ax.set_ylim(ymin - yspan*fracSpace, ymax + yspan*fracSpace)
        # either run set_xlim or self.update_labels(redraw=False) to update labels
        self.update_labels(redraw=False)
        self.fig.canvas.draw()
        return


def _fancy_ion_string(ion):
    '''
    Convert the ion string from VALD's format to a fancier format
    using roman numerals for the ionization stage.
    (Currently only supports ionization up to IX)
    '''
    fion = copy.deepcopy(ion)
    trans = str.maketrans({'1':'I', '2':'II', '3':'III', '4':'IV', '5':'V',
                           '6':'VI', '7':'VII', '8':'VIII', '9':'IX'})
    for i in range(len(fion)):
        text = fion[i].translate(trans)
        fion[i] = text
    return fion


def _get_visible_labels(wlMin, wlMax, show_lines, llist_all, maxLabels):
    '''
    Get an array of flags for line labels to be shown.
    Limit labels plotted to lines being plotted, in the right wavelength range,
    and only show the maxLabels deepest lines.
    Return an array of flags for visible labels,
    and indices for the shown labels in the set of shown lines.
    '''
    show_labels = (llist_all.wl > wlMin) & (llist_all.wl < wlMax) & show_lines
    
    # Limit the number of text labels plotted, only plotting the deepest lines
    # set flags for the sub-set of the line list used for drawing text labels
    # and get indices for their positions in the full line list being plotted
    if np.sum(show_labels) > maxLabels:
        # get only the maxLabels deepest lines
        indSortDepth = np.argsort(llist_all.depth[show_labels])
        indUse = indSortDepth[-maxLabels:]
        # sort the array of indices to keep back into wavelength order
        ind_show_labels_use = indUse[np.argsort(llist_all.wl[show_labels][indUse])]
        # set the show flag (in an array of all labels) for only those labels to keep
        ind_update = np.nonzero(show_labels)[0]
        show_labels[:] = False
        show_labels[ind_update[ind_show_labels_use]] = True
        
        # get indices for the labels shown in the array of lines shown (plotted)
        ind_labels_in_lines = np.nonzero(show_labels[show_lines])[0]
    else:
        ind_labels_in_lines = np.nonzero(show_labels[show_lines])[0]

    return show_labels, ind_labels_in_lines


def _label_overlap_regions(llabels2d, padding, closeRegions=None, textBBcorn=None):
    '''
    Check if line labels are overlapping in the plot, using screen space
    (pixels). Flag lines that all overlap as one overlapping
    region group, and give them the same value in the closeRegions array.
    Lines can be iteratively added to the closeRegions array,
    with subsequent function calls.
    
    * llabels2d - list of lists, containing rows of matplotlib text objects
    * padding - leave this many extra pixels around labels
    * closeRegions - an existing list of arrays of close regions, to be updated.
                     If None then a new list of arrays is created.
    * textBBcorn - list of 3d arrays with bounding box corners, used to check
                   for overlapping labels.  If None then new bounding boxes 
                   are gotten from text objects in llabels2d (this is slower 
                   by 10x, but ensures llabels2d and textBBcorn are consistent).
    
    returns
    * numOverlaps - total number of overlapping labels found
    * closeRegions - list of arrays of close/overlapping regions
    * textBBcorn - list of arrays of text label bounding box corners
    '''
    # (allow a margin since data - pixel conversions may not be perfect)
    marginForError = 1.0 #allowable overlap in pixels
    
    nrows = len(llabels2d)
    numOverlaps = 0
    # initialize close regions to 0 if we don't have them
    if closeRegions is None:
        closeRegions = []
        for m in range(nrows):
            closeRegions += [np.zeros(len(llabels2d[m]), dtype=int)]

    # get the bounding box corners for llabels2d if we don't have them
    if textBBcorn is None:
        textBBcorn = []
        for m in range(nrows):
            _textBBcorn = np.zeros((len(llabels2d[m]), 2, 2))
            for i, txt in enumerate(llabels2d[m]):
                # get_window_extent() is a slow function call...
                bbox = txt.get_window_extent()
                _textBBcorn[i,:,:] = bbox.get_points()
            textBBcorn += [_textBBcorn]

    # update close regions using the bounding box corners
    for m in range(nrows):
        _textBBcorn = textBBcorn[m]
        _closeRegions = closeRegions[m]
        if len(llabels2d[m]) > 0: #protect against empty rows/lists
            _numClose = np.max(_closeRegions)
            xlast = -100000.0
            for i, txt in enumerate(llabels2d[m]):
                x = _textBBcorn[i,0,0]
                
                if x + marginForError < xlast + padding:
                    numOverlaps += 1
                    if _closeRegions[i-1] == 0 and _closeRegions[i] == 0:
                        # if we need to define a new close region
                        _numClose += 1
                        _closeRegions[i-1] = _numClose
                        _closeRegions[i] = _numClose
                    elif _closeRegions[i-1] == 0:
                        # add the left line to this region 
                        _closeRegions[i-1] = _closeRegions[i]
                    elif _closeRegions[i] == 0:
                        # add the right line to this region
                        _closeRegions[i] = _closeRegions[i-1]                        
                    else:
                        # merge the right region into the left region
                        _closeRegions[_closeRegions == _closeRegions[i]] = _closeRegions[i-1]
                
                xlast = _textBBcorn[i,1,0] #set up for next iteration
            
    return numOverlaps, closeRegions, textBBcorn


def _set_label_pos(fig, ax, llabels2d, llist_labels, padding, vpadding, redraw=True):
    """
    Modify the positions of labels on the plot, to avoid having overlapping
    text labels.

    Takes:
    * fig - the figure object containing the plot
    * ax - the axes object containing the plot
    * llabels2d - a list of lists, with label text objects in lists for each row
    * padding - extra space to leave around the labels, in pixels
    * vpadding - extra vertical space to leave between rows of labels, in pixels

    Returns: labelWls, labelFlx - arrays with the x and y positions of the labels, 
                                  in data coordinates (wavelength and flux)
    """

    nrows = len(llabels2d)
    # Draw the figure so that overlapping labels can be identified
    if redraw: fig.draw_without_rendering() #this drawing call is very slow...
    ##if redraw: fig.canvas.draw()  #... but drawing and rendering is even slower
    # get a function for converting screen space coordinates back into data coordinates
    invertCoords = ax.transData.inverted()

    numOverlaps, closeRegions, textBBcorn = _label_overlap_regions(llabels2d, padding)

    # If the line labels wont all fit without overlap, don't try adjusting positions
    bbox_ax = ax.get_window_extent()
    win_x, win_y, win_w, win_h = bbox_ax.bounds
    if numOverlaps > 0:
        for m in range(nrows):
            tot_width = np.sum(textBBcorn[m][:,1,0] - textBBcorn[m][:,0,0])
            if win_w < tot_width + padding*len(llabels2d[m]):
                numOverlaps = -1

    # If labels are off the edge of the window, reposition them
    if numOverlaps > -1:
        for m in range(nrows):
            if textBBcorn[m].size > 0: #protect against empty rows/lists
                #check that the first label fits on the plot
                if textBBcorn[m][0, 0, 0] < win_x:
                    xwidth = textBBcorn[m][0, 1, 0] - textBBcorn[m][0, 0, 0]
                    xnew = win_x + 0.5*xwidth + padding
                    dataCoords = invertCoords.transform([xnew, textBBcorn[m][0, 0, 1]])
                    llabels2d[m][0].set_x(dataCoords[0])
                #check that the last label fits on the plot
                if textBBcorn[m][-1, 1, 0] > win_x + win_w:
                    xwidth = textBBcorn[m][-1, 1, 0] - textBBcorn[m][-1, 0, 0] 
                    xnew = win_x + win_w - 0.5*xwidth - padding
                    dataCoords = invertCoords.transform([xnew, textBBcorn[m][-1, 0, 1]])
                    llabels2d[m][-1].set_x(dataCoords[0])

    # adjust label vertical positions based on the current label size
    # in display coordinates (assume all labels have the same height)
    if textBBcorn[0].size > 0:  #protect against empty rows/lists
        lheigh = textBBcorn[0][0,1,1] - textBBcorn[0][0,0,1]
        for m in range(1, nrows):
            newy = m*(lheigh + vpadding) + textBBcorn[0][0,0,1]
            dataCoords = invertCoords.transform([textBBcorn[0][0,0,0], newy])
            datNewy = dataCoords[1]
            for txt in llabels2d[m]:
                txt.set_y(datNewy)
            #This leaves the y values in textBBcorn un-updated, but that's probably ok

    flagLEdge = [0] * nrows # list of flags, starting at 0
    flagREdge = [0] * nrows
    nIter = 0
    while(numOverlaps > 0 and nIter < 10):
        nIter += 1
        for m in range(nrows):
            
            # loop over groups of close labels and reposition all of them
            for n in np.unique(closeRegions[m][closeRegions[m] > 0]):
                # get the indices for labels in this overlap region
                iuse = np.nonzero(closeRegions[m] == n)[0]

                ## Get the position of the middle label (or midpoint of two labels)
                ## and the index for the label(s) at the middle of the region
                ## in the overlap (for an odd number of labels imidLo == imidUp)
                #imidLo = iuse[0] + (iuse[-1] - iuse[0])//2
                #imidUp = iuse[0] + (iuse[-1] - iuse[0] + 1)//2
                #xmid = 0.5*(textBBcorn[m][imidLo, 0, 0] + textBBcorn[m][imidUp, 1, 0])
                
                # Get the center of the labels in this overlap region
                xmid = 0.5*(textBBcorn[m][iuse[0], 0, 0] + textBBcorn[m][iuse[-1], 1, 0])
                
                #Position labels in this close group around the mid point
                tot_width = np.sum(textBBcorn[m][iuse,1,0] - textBBcorn[m][iuse,0,0])
                tot_width += padding*(iuse.size - 1)
                
                #set the left edge x position for the first label
                xnewL = xmid - 0.5*tot_width
                #check that it fits on the plot (if a group hits an edge stick to it)
                if (xnewL < win_x) or (flagLEdge[m] == 1 and n == 1):
                    xnewL = win_x + padding
                    flagLEdge[m] = 1
                #check that the last label fits on the plot (keep that group at that edge)
                elif (xnewL + tot_width > win_x + win_w) or (
                        flagREdge[m] == 1 and n == np.max(closeRegions[m])):
                    xnewL = win_x + win_w - tot_width - padding
                    flagREdge[m] = 1

                #adjust the positions for the labels in this group
                for i in iuse:
                    xwidth = textBBcorn[m][i, 1, 0] - textBBcorn[m][i, 0, 0] 
                    xnew = xnewL + 0.5*xwidth
                    dataCoords = invertCoords.transform([xnew, textBBcorn[m][i, 0, 1]])
                    # set the label position
                    llabels2d[m][i].set_x(dataCoords[0])
                    # update stored corner positions, for comparison with other labels 
                    textBBcorn[m][i, 0, 0] = xnewL
                    textBBcorn[m][i, 1, 0] = xnewL + xwidth
                    # set up for the next label
                    xnewL += xwidth + padding
        
        # done one full iteration, check if there are any remaining overlaps
        numOverlaps, closeRegions, textBBcorn = _label_overlap_regions(
            llabels2d, padding, closeRegions, textBBcorn)

        # if there are still overlaps, reset label positions for next iteration
        if numOverlaps > 0:
            screenx = ax.transData.transform(np.array([llist_labels.wl,
                                        np.ones_like(llist_labels.wl)]).T)
            for i in range(llist_labels.nLines):
                j = i%nrows
                k = i//nrows
                if(closeRegions[j][k] > 0):
                    llabels2d[j][k].set_x(llist_labels.wl[i])
                    screenw = textBBcorn[j][k,1,0] - textBBcorn[j][k,0,0]
                    textBBcorn[j][k,0,0] = screenx[i,0] - 0.5*screenw
                    textBBcorn[j][k,1,0] = screenx[i,0] + 0.5*screenw

    # Save the label positions for drawing ticks
    nPlotted = llist_labels.nLines
    labelWls = np.zeros(nPlotted)
    labelFlx = np.zeros(nPlotted)
    for i in range(nPlotted):
        j = i%nrows
        k = i//nrows
        labelPos = llabels2d[j][k].get_position()
        labelWls[i] = labelPos[0]
        labelFlx[i] = labelPos[1]
    
    return labelWls, labelFlx
