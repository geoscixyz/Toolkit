# -*- coding: utf-8 -*-
"""
Functions for creating and manipulating graphics, colormaps and plots.

@author: Joseph Barraud, Geophysics Labs.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from skimage import exposure

# import local modules
from . import colors # definition of new colormaps

#==============================================================================
# stats_boundaries
#==============================================================================
def stats_boundaries(data,nSigma=1,sigmaStep=1):
    '''
    Return a list of statistical values ordered in increasing order that can
    be used for ticks or boundaries.
    '''
    mu = np.nanmean(data)
    sigma = np.nanstd(data)
    newTicks = mu + sigma*np.arange(-nSigma,nSigma+sigmaStep,sigmaStep)

    return [np.nanmin(data)] + newTicks.tolist() + [np.nanmax(data)]

#===============================================================================
# makeColormap
#===============================================================================
def makeColormap(table,name='CustomMap'):
    """
    Return a LinearSegmentedColormap.
    table: a sequence of RGB tuples. Values need to be between 0 and 1.
    colormap is subsequently registered with plt.register_cmap(cmap=my_cmap)
    """
    if np.any(table > 1):
        table = table / 255.
    cdict = {'red': [], 'green': [], 'blue': []}
    N = float(len(table))-1
    for i,rgb in enumerate(table):
        r1, g1, b1 = rgb
        cdict['red'].append([i/N, r1, r1])
        cdict['green'].append([i/N, g1, g1])
        cdict['blue'].append([i/N, b1, b1])

    new_cmap = mcolors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=new_cmap)
    return new_cmap

#===============================================================================
# cmap_to_array
#===============================================================================
def cmap_to_array(cmap,N=256):
    """
    Return a Nx3 array of RGB values generated from a colormap.
    """
    return cmap(np.linspace(0, 1, N))[:,:3] # remove alpha column

#===============================================================================
# load_cmap
#===============================================================================
def load_cmap(name='geosoft'):
    """
    Create and register a color map. The name must be selected from the
    colormaps available in the local colors module.
    """
    try:
        cmList = colors.datad[name]
        new_cm = mcolors.LinearSegmentedColormap.from_list(name, cmList)
        plt.register_cmap(cmap=new_cm)
        return new_cm
    except:
        raise ValueError('Colormap {} has not been recognised'.format(name))

#==============================================================================
# autolevels
#==============================================================================
def autolevels(image,minPercent=2,maxPercent=98,funcName='mean',perChannel=False):
    '''
    Rescale intensity of an image. For RGB images, the new limits are calculated
    per channel and then mean or median of these limits are applied to the whole
    image (if perChannel option is False).
    '''
    # dictionary of functions
    funcs = {'mean':np.mean,'median':np.median,'min':np.min,'max':np.max}

    # calculate percentiles (returns 3 values for RGB pictures or vectors, 1 for grayscale images)
    if image.shape[1] == 3:
        pMin,pMax = np.percentile(image,(minPercent, maxPercent),axis=0)
    else:
        pMin,pMax = np.percentile(image,(minPercent, maxPercent),axis=(0,1))

    # Apply normalisation
    if not perChannel: # finds new min and max using selected function applied to all channels
        newMin = funcs[funcName](pMin)
        newMax = funcs[funcName](pMax)
        auto = exposure.rescale_intensity(image,in_range=(newMin,newMax))

    else: # applies a rescale on each channel separately
        r_channel = exposure.rescale_intensity(image[:,:,0], in_range=(pMin[0],pMax[0]))
        g_channel = exposure.rescale_intensity(image[:,:,1], in_range=(pMin[1],pMax[1]))
        b_channel = exposure.rescale_intensity(image[:,:,2], in_range=(pMin[2],pMax[2]))
        auto = np.stack((r_channel,g_channel,b_channel),axis=2)

    return auto

#===============================================================================
# equalizeColormap
#===============================================================================
def equalizeColormap(cmap,bins,cdf,name='EqualizedMap'):
    '''
    Re-map a colormap according to a cumulative distribution. This is used to
    perform histogram equalization of an image by changing the colormap
    instead of the image. *This is not strickly speaking the equalization of the
    colormap itself*.
    The cdf and bins should be calculated from an input image, as if carrying out
    the histogram equalization of that image. In effect, the cdf becomes integrated
    to the colormap as a mapping function by redistributing the indices of the
    input colormap.

    Parameters
    ----------
    cmap : string or colormap object
        Input colormap to remap.
    bins : array
        Centers of bins.
    cdf : array
        Values of cumulative distribution function.
    '''

    # first retrieve the color table (lists of RGB values) behind the input colormap
    if cmap in colors.datad: # one of the additional colormaps in colors module
        cmList = colors.datad[cmap]
    elif cmap in cm.cmap_d: # matplotlib colormaps + plus the new ones (viridis, inferno, etc.)
        cmList = cmap_to_array(cm.cmap_d[cmap])
    else:
        try:
            # in case cmap is a colormap object
            cmList = cmap_to_array(cmap)
        except:
            raise ValueError('Colormap {} has not been recognised'.format(cmap))

    # normalize the input bins to interval (0,1)
    bins_norm = (bins - bins.min())/np.float(bins.max() - bins.min())

    # calculate new indices by applying the cdf as a function on the old indices
    # which are initially regularly spaced.
    old_indices = np.linspace(0,1,len(cmList))
    new_indices = np.interp(old_indices,cdf,bins_norm)

    # make sure indices start with 0 and end with 1
    new_indices[0] = 0.0
    new_indices[-1] = 1.0

    # remap the color table
    cdict = {'red': [], 'green': [], 'blue': []}
    for i,n in enumerate(new_indices):
        r1, g1, b1 = cmList[i]
        cdict['red'].append([n, r1, r1])
        cdict['green'].append([n, g1, g1])
        cdict['blue'].append([n, b1, b1])

    return mcolors.LinearSegmentedColormap(name, cdict)

#===============================================================================
# normalizeColormap
#===============================================================================
def normalizeColormap(cmapName,norm='autolevels',**kwargs):
    '''
    Apply a normalising function to a colormap. Only "autolevels" is implemented
    for the moment.

    **kwargs are passed to the normalising function.
    '''
    try:
        cmap = cm.get_cmap(cmapName) # works even if cmapName is already a colormap
    except:
        # colormap is one of the extra ones added by the colors module
        cmap = load_cmap(cmapName)

    # convert cmap to array for normalisation
    cmList = cmap_to_array(cmap)

    # normalise
    if norm == 'autolevels':
        cmList_norm = autolevels(cmList,**kwargs)
    else:
        cmList_norm = cmList

    # create new colormap
    new_cm = mcolors.LinearSegmentedColormap.from_list(cmap.name + '_n', cmList_norm)

    return new_cm

#===============================================================================
# alpha_blend
#===============================================================================
def alpha_blend(rgb, intensity, alpha = 0.7):
    """
    Combines an RGB image with an intensity map using "alpha" transparent blending.
    https://en.wikipedia.org/wiki/Alpha_compositing

    Parameters
    ----------
    rgb : ndarray
        An MxNx3 RGB array of floats ranging from 0 to 1 (color image).
    intensity : ndarray
        An MxNx1 array of floats ranging from 0 to 1 (grayscale image).
    alpha : float
        This controls the transparency of the rgb image. 1.0 is fully opaque
        while 0.0 is fully transparent.

    Returns
    -------
    rgb : ndarray
        An MxNx3 RGB array representing the combined images.
    """

    return alpha*rgb + (1 - alpha)*intensity

#===============================================================================
# imshow_hs
#===============================================================================
def imshow_hs(data,ax=None,cmap='geosoft',cmap_norm='equalize',hs=True,
              zf=10,azdeg=45,altdeg=45,dx=1,dy=1,fraction=1.5,blend_mode='alpha',
              alpha=0.7,contours=False,levels=32,colorbar=True,cb_contours=False,
              cb_ticks='linear',nSigma=1,**kwargs):
    '''
    Display an array with optional hillshading and contours. The colormap can be
    normalised by equalisation or by clipping extremes (autolevels).

    Parameters
    ----------
    data : 2D array
        Grid to plot. Arrays with NaNs and masked arrays are supported.
    ax : matplotlib axes instance
        This indicates where to draw the figure. Create new figure if absent.
    cmap : string
        Name of the colormap to use to display the array. The default 'geosoft' is
        the blue to pink clra colormap from Geosoft Oasis Montaj.
    cmap_norm : string
        Type of normalisation of the colormap.
        Possible values are:
            'equalize' (or 'equalization')
                Increases contrast by distributing intensities across all the
                possible colours. With this option, it is not the data that is normalised
                but the colormap, based on the data.
            'auto' (or 'autolevels')
                Stretches the histogram of the colormap so that dark colours become
                darker and the bright colours become brighter. Two extra parameters control
                the amount of clipping at the extremes: minPercent (default to 10%) and
                maxPercent (default to 90%)
    hs : boolean
        If True, the array is displayed in colours over a grey hillshaded version
        of the data.
    zf : number
        Vertical exaggeration (Z factor) for hillshading.
    azdeg : number
        The azimuth (0-360, degrees clockwise from North) of the light source.
    altdeg : number
        The altitude (0-90, degrees up from horizontal) of the light source.
    dx : number, optional
        cell size in the x direction
    dy : number, optional
        cell size in the y direction
    fraction : number
        Increases or decreases the contrast of the hillshade.
    blend_mode :  {'alpha', 'hsv', 'overlay', 'soft'}
        The type of blending used to combine the colormapped data values with the
        illumination intensity. Default is 'alpha' and the effect is controlled
        by the alpha parameter.
    alpha : float
        Controls the transparency of the data overlaid over the hillshade.
        1.0 is fully opaque while 0.0 is fully transparent.
    contours : Boolean
        If True, adds contours to the map. The number of calculated contours is
        defined by:
            levels : integer
                Number of contour levels.
    colorbar : Boolean
        If True, draw a colorbar on the right-hand side of the map. The colorbar
        shows the distribution of colors, as modified by the normalization algorithm.
    cb_ticks : string
        If left as default ('linear') the ticks and labels on the colorbar are
        spaced linearly in the standard way. Otherwise (any other keyword, for example
        'stats'), the mean and two ticks at + and - nSigma*(standard deviation)
        are shown instead.
            nSigma : integer (default is 1)
                Size of the interval to show between ticks on the colorbar.
    cb_contours : Boolean
        Add lines corresponding to contours on the colorbar.
    kwargs : other optional arguments
        Can be used to pass other arguments to imshow, such as 'origin' and 'extent'.

    Notes
    -----
    This function exploits the hillshading capabilities implemented in
    matplotlib.colors.LightSource. It adds additional blending mode (alpha compositing,
    see https://en.wikipedia.org/wiki/Alpha_compositing) and normalising functions
    for the data (equalization).

    '''
    # modify colormap if required
    if cmap_norm in ['equalize','equalization']:
        # histogram equalization
        cdf, bins = exposure.cumulative_distribution(data[~np.isnan(data)].flatten(),nbins=256)
        my_cmap = equalizeColormap(cmap,bins,cdf)
    elif cmap_norm in ['auto','autolevels']:
        # autolevels
        minP = kwargs.pop('minPercent',10) # also removes the key from the dictionary
        maxP = kwargs.pop('maxPercent',90)
        my_cmap = normalizeColormap(cmap,norm='autolevels',minPercent=minP,maxPercent=maxP)
    elif cmap in plt.colormaps():
        # colormap defined as string (recognised name)
        my_cmap = plt.get_cmap(cmap)
    else:
        # colormap is one of the extra ones added by the colors module
        my_cmap = load_cmap(cmap)  # raises error if not recognised

    # create figure or retrieve the one already there
    if ax:
        fig = ax.get_figure()
    else:
        fig,ax = plt.subplots(figsize=(8,8))

    # convert input data to masked array
    data = np.ma.masked_array(data, np.isnan(data))

    # add array to figure with hillshade or not
    if hs:
        # flip azimuth upside down if grid is also flipped
        if 'origin' in kwargs:
            if kwargs['origin'] == 'lower':
                azdeg = 180 - azdeg

        # extract keyword arguments that can be passed to ls.shade
        kwargs_norm = {}
        kwargs_norm['norm'] = kwargs.get('norm')
        kwargs_norm['vmin'] = kwargs.get('vmin')
        kwargs_norm['vmax'] = kwargs.get('vmax')

        # create light source
        ls = mcolors.LightSource(azdeg, altdeg)

        # calculate hillshade and combine the colormapped data with the intensity
        if alpha == 0:
            # special case when only the shaded relief is needed without blending
            rgb = ls.hillshade(data,vert_exag=zf,dx=dx,dy=dy,fraction=fraction)
            kwargs['cmap'] = 'gray'

        elif blend_mode == 'alpha':
            # transparency blending
            rgb = ls.shade(data,cmap=my_cmap,blend_mode=alpha_blend,vert_exag=zf,dx=dx,dy=dy,
                           fraction=fraction,alpha=alpha,**kwargs_norm)

        else:
            # other blending modes from matplotlib function
            rgb = ls.shade(data,cmap=my_cmap,blend_mode=blend_mode,vert_exag=zf,dx=dx,dy=dy,
                           fraction=fraction,**kwargs_norm)

        ax.imshow(rgb,**kwargs)

    else:
        # display data without hillshading
        im = ax.imshow(data,cmap=my_cmap,**kwargs)

    # add contours
    if contours:
        ct = plt.contour(data,levels,linewidths=0.5,colors='k',linestyles='solid',**kwargs)

    # add colorbar
    if colorbar and alpha != 0:
        if hs:
            # Use a proxy artist for the colorbar
            im = ax.imshow(data,cmap=my_cmap,**kwargs)
            im.remove()

        if cb_ticks=='linear': # normal equidistant ticks on a linear scale
            cb1 = fig.colorbar(im,shrink=0.8)
        else: # show ticks at min, max, mean and standard deviation interval
            newTicks = stats_boundaries(data,nSigma,nSigma)
            cb1 = fig.colorbar(im,shrink=0.8,ticks=newTicks)

        # add optional contour lines on colorbar
        if contours and cb_contours:
            cb1.add_lines(ct)

        cb1.update_normal(im)

    # final show
    plt.show()

#==============================================================================
# saveMap
#==============================================================================
def saveMap(outfile,fig=None,orig_size=None,dpi=100):
    '''
    Save a Matplotlib figure as an image without borders or frames.
    Parameters
    ----------
    outfile (string): Path to output file.

    fig (Matplotlib figure instance): figure you want to save as the image

    orig_size (tuple): width, height of the original image used to maintain
    aspect ratio.

    dpi (integer): image resolution.
    '''
    if fig==None:
        fig = plt.gcf()
    ax = fig.gca()
    ax.set_axis_off()
    ax.set_position([0,0,1,1])
    ax.set_aspect('auto')
    fig.set_frameon(False)

    if orig_size: # Aspect ratio scaling if required
        w,h = orig_size
        fig.set_size_inches(w/float(dpi),h/float(dpi),forward=False)

    fig.savefig(outfile,dpi=dpi)

