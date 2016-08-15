# -----------------------------------------------------------------------------
#
from numpy import arange, fromstring, array

import os.path
from os.path import isfile, isdir, join, split, splitext
from os.path import getsize, isabs, exists, abspath

import zipfile
import gzip


def addext(filename, extension):
    """Returns *filename*, with *extension* if it does not have one."""

    return filename + ('' if splitext(filename)[1] else extension)

def gzip_open(filename, *args, **kwargs):
    if args and "t" in args[0]:
        args = (args[0].replace("t", ""), ) + args[1:]
    if isinstance(filename, str):
        return gzip.open(filename, *args, **kwargs)
    else:
        return gzip.GzipFile(filename, *args, **kwargs)

OPEN = {
    '.gz': gzip_open,
    '.zip': zipfile.ZipFile,
}

def openFile(filename, *args, **kwargs):
    """Open *filename* for reading, writing, or appending.  First argument in
    *args* is treated as the mode.  Opening :file:`.gz` and :file:`.zip` files
    for reading and writing is handled automatically.

    :arg backup: backup existing file using :func:`.backupFile` when opening
        in append or write modes, default is obtained from package settings
    :type backup: bool

    :arg backup_ext: extension for backup file, default is :file:`.BAK`
    :type backup_ext: str"""

    try:
        exists = isfile(filename)
    except Exception as err:
        raise TypeError('filename must be a string ({0})'.format(str(err)))

    folder = kwargs.pop('folder', None)
    if folder:
        filename = join(folder, filename)

    backup = kwargs.pop('backup', None)
    if backup is not None and backup and args and args[0][0] in ('a', 'w'):
        backupFile(filename, backup=backup,
                   backup_ext=kwargs.pop('backup_ext', None))

    ext = splitext(filename)[1]
    return OPEN.get(ext.lower(), open)(filename, *args, **kwargs)

def startswith(this, that):
    """Returns **True** if *this* or *that* starts with the other."""

    if len(this) < len(that):
        return that.startswith(this)
    else:
        return this.startswith(that)

def intorfloat(x):
    """Returns ``int(x)``, or ``float(x)`` upon :exc:`ValueError`."""

    try:
        return int(x)
    except ValueError:
        return float(x)

HMTYPES = {
    'max': float,
    'min': float,
    'numbering': lambda line: line.split(':'),
    'title': str,
    'xlabel': str,
    'ylabel': str,
    'xorigin': intorfloat,
    'xstep': intorfloat,
}

def showHeatmap(heatmap, *args, **kwargs):
    """Show *heatmap*, which can be an two dimensional array or a Heat Mapper
    :file:`.hm` file.

    Heatmap is plotted using :func:`~matplotlib.pyplot.imshow` function.
    Default values passed to this function are ``interpolation='nearest'``,
    ``aspect='auto'``, and ``origin='lower'``."""

    try:
        ndim, shape = heatmap.ndim, heatmap.shape
    except AttributeError:
        heatmap, headers = parseHeatmap(heatmap)
        ndim, shape = heatmap.ndim, heatmap.shape

        xorigin = headers.pop('xorigin', 0)
        xextent = headers.pop('xstep', 1) * shape[0]

        ylabel = kwargs.get('ylabel', '').lower()
        indices = None
        if ylabel:
            for key in headers.get('numbering', []):
                if startswith(ylabel, key.lower()):
                    indices = headers.get(key)
        if indices is not None:
            extent = [indices[0] - .5, indices[0] + len(indices) - .5,
                      xorigin - .5, xextent - .5]
        else:
            extent = [-.5, shape[1] * 2 - .5, xorigin - .5, xextent - .5]
        kwargs.setdefault('extent', extent)

        for key in ['numbering', 'min', 'max'] + headers.get('numbering', []):
            headers.pop(key, None)
        headers.update(kwargs)
        kwargs = headers


    if ndim != 2:
        raise ValueError('mutinfo must be a 2D matrix')

    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault('origin', 'lower')
    kwargs.setdefault('aspect', 'auto')

    xlabel = kwargs.pop('xlabel', None)
    ylabel = kwargs.pop('ylabel', None)
    title = kwargs.pop('title', None)

    import matplotlib.pyplot as plt
    show = plt.imshow(heatmap, *args, **kwargs), plt.colorbar()

    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    return show


def parseHeatmap(heatmap, **kwargs):
    """Returns a two dimensional array and a dictionary with information parsed
    from *heatmap*, which may be an input stream or an :file:`.hm` file in VMD
    plugin Heat Mapper format."""

    try:
        readline, close = heatmap.readline, lambda: None
    except AttributeError:
        heatmap = openFile(heatmap)
        readline, close = heatmap.readline, heatmap.close

    meta = {}
    arrs = []

    line = readline()
    while line:
        if line.startswith('-'):
            label, data = line[1:].split(None, 1)
            data = data.strip()
            if data[0] == data[-1] == '"':
                data = data[1:-1]
            label = label.strip()
            try:
                meta[label] = HMTYPES[label](data)
            except KeyError:
                LOGGER.warn('Unrecognized label encountered: {0}'
                            .format(repr(label)))
                meta[label] = HMTYPES[label](data)
            except TypeError:
                LOGGER.warn('Could not parse data with label {0}.'
                            .format(repr(label)))
        else:
            arrs.append(line.rstrip())
        line = readline()
    close()
    nnums = len(meta.get('numbering', ''))
    heatmap = []
    numbers = []

    for arr in arrs:
        if nnums:
            items = arr.split(':', nnums + 1)
            numbers.append(items[:nnums])
        else:
            items = [arr]
        heatmap.append(fromstring(items[-1], float, sep=';'))

    heatmap = array(heatmap)
    if nnums:
        numbering = meta['numbering']
        try:
            numbers = array(numbers, int)
        except ValueError:
            try:
                numbers = array(numbers, float)
            except ValueError:
                LOGGER.warn('Numbering for y-axis could not be parsed.')
                numbering = []
        for i, label in enumerate(numbering):
            meta[label] = numbers[:, i].copy()

    return heatmap, meta

def helpTextWizard():
    text = ("ProDy Interface\n" 
      "===============\n\n"
      "\nActive Mode\n"
      "--------------\n\n"
      "Select the active mode for which you want to draw arrows or make an animation. "
      "Direction of arrows depicting the normal mode can be changed using +/- button. "
      "Arrows can be drawn along both directions by changing the options Mode Graphics Options panel. "
      "The selected color effects both arrow graphics and square fluctuation plots."
      "\n\n**RMSD**\n\n"
      "The RMSD corresponding to the displacement described by the arrows is displayed. User can change the RMSD value to rescale the arrows. "
      "The scaling factor that produces the specified RMSD is printed to the VMD console (along with the magnitude of the mode provided in NMD file). "
      "\n\n**Selection**\n\n"
      "Selection entry allows the user to display arrows for a subset of atoms.\n\n"
      "*TIP*: If the arrow graphics are too crowded or the display is slow, draw arrows for an evenly spaced subset of residues, e.g try 'name CA and residue % 4 == 0', which will draw an arrow for every fourth residue."
      "\n\n\n"
      "Mode Graphics\n"
      "--------------\n\n"
      "Id of the molecule that contains the arrow graphics of the active mode is shown in parentheses.\n\n"
      "Buttons:\n\n"
      " * Draw: draw/redraw arrow graphics for the active mode\n"
      " * Clean: remove most recently drawn arrow graphics\n"
      " * Hide/Show: hide/show most recently drawn arrow graphics\n"
      " * Options: show/hide arrow graphics option panel\n"
      "\nOptions:\n\n"
      "User can change arrow graphics properties and how NMWiz behave upon such changes in this panel.\n"
      "\nBy default:\n\n"
      " * arrow graphics are set to automatically change when graphics properties are changed by the user\n"
      " * current graphics are hidden the active mode is changed\n"
      "\nOptionally:\n\n"
      " * arrows can be drawn in both directions to look like a double headed arrow\n"
      " * arrows shorter than a length (A) can be hidden\n"
      "\nAdditionally, user can change:\n\n"
      " * width of the arrow cylinder\n"
      " * width/height of the arrow head code\n"
      " * graphics material and resolution"
      "\n\n\n"
      "Animations\n"
      "----------\n\n"
      "Id of the molecule that contains the most recently generated animation is shown in parentheses.\n\n"
      "Buttons:\n\n"
      " * Draw: animate fluctuations along the active mode\n"
      " * Play : play/pause the animation\n"
      " * Hide : hide/show the animation\n"
      " * Options: show/hide animation option panel\n"
      "\nOptions:\n\n"
      "User can elect automatic generation and continuous play of animations when the active mode changes. User can also select the number of frames in the animation."
      "\n\n\n"
      "Plotting\n"
      "--------\n\n"
      "Id of the molecule for displaying selected residues is shown in parentheses.\n\n"
      "Buttons:\n\n"
      " * Plot: plot squared-fluctuations along the active mode\n"
      " * Clear: clear all selections and selected atom labels\n"
      " * Hide/Show: hide/show the selected residues\n"
      " * Options: change plotting options\n"
      "\n\n\n"
      "Molecule Representations\n"
      "------------------------\n\n"
      "Id of the molecule that contains the structure is shown in parentheses.\n\n"
      "Buttons:\n\n"
      " * Update: update molecule representation\n"
      " * Focus: reset view to focus on the structure\n"
      " * Hide/Show: hide/show strudture\n"
      " * Options: change molecular system representation\n"
      "\nOptions:\n\n"
      "User can select the representation and coloring scheme. User can change the molecule representation settings manually, by setting 'Show structure as' to 'Custom'.\n\n"
      "Structure can be colored based on the `Mobility` of the residues in the active mode, based on 'Bfactors' that came in NMD file, or based on residue/atom 'Index'.\n\n"
      "In addition to the standard representations (e.g. Tube/Trace/Licorice), structure can be represented as an elastic network. Color scale method and midpoint can be used to adjust mobility and Bfactors based coloring."
      "User can set the cutoff distance, width of dynamic bonds, and node spheres. Note that changing the cutoff distance distance only affects representation, not the precalculated normal mode data.\n\n"
      "*TIP*: When visualizing a large system, display molecule at lower resolutions and/or try displaying fewer atoms if all atoms are displayed.")
    return text