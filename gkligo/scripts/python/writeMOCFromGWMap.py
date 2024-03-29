#!/usr/bin/env python
"""
Write MOCS from GW Maps - fixes FITS header too.
https://emfollow.docs.ligo.org/userguide/tutorial/multiorder_skymaps.html

This code written by Ken Smith based on code by Leo Singer and Roy Williams 

Requires:
numpy
ligo.skymap
astropy
mocpy

Usage:
  %s <skymap> [--directory=<directory>] [--contours=<contours>] [--logfile=<logfile>] [--organise]
  %s (-h | --help)
  %s --version

where action is start|stop|restart|listen

If action is listen, start in non-daemon mode. Otherwise use daemon mode.

Options:
  -h --help                         Show this screen.
  --version                         Show version.
  --directory=<directory>           Directory to where the maps and MOCs will be written [default: /tmp].
  --contours=<contours>             Which MOC contours do you want? Multiple contours should be separated by commas, with no spaces [default: 90]
  --logfile=<logfile>               PID file [default: /tmp/ligo.log]
  --organise                        Instead of writing a long unique filename, write a directory structure.

E.g.:
  %s /home/ligo/GWmap.fits --directory=/home/ligo/maps_gkligo --writemeta --contours=90,50,10

"""
import sys
__doc__ = __doc__ % (sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
from docopt import docopt

import json
import base64
from gkutils.commonutils import Struct, cleanOptions
from io import BytesIO
import os
import time
import logging
from copy import deepcopy


def getContourArea(inputFilePointer, contour, logger):

    from astropy.table import Table
    from astropy import units as u
    import numpy as np
    import math
    from ligo.skymap.moc import uniq2pixarea

    # Read and verify the input
    skymap = Table.read(inputFilePointer, format='fits')

    # Sort by prob density of pixel
    skymap.sort('PROBDENSITY', reverse=True)

    # Get area*probdensity for each pixel
    pixel_area = uniq2pixarea(skymap['UNIQ'])

    # Probability per pixel
    prob = pixel_area*skymap['PROBDENSITY']
    cumprob = np.cumsum(prob)

    # Should be 1.0. But need not be.
    sumprob = np.sum(prob)

    # Find the index where contour of prob is inside
    i = cumprob.searchsorted(contour*sumprob)
    area = float(pixel_area[:i].sum() * (180/math.pi)**2)

    return area



def writeMOC(skymap, outputMOCName, contour, logger):
    """writeMOC.

    Args:
        inputFilePointer:
        outputMOCName:
        contour:
        logger:
    """
    from astropy.table import Table
    from astropy.io import fits
    from astropy import units as u
    import numpy as np
    import math
    from ligo.skymap.moc import uniq2pixarea

    # Read and verify the input
    skymap = Table.read(open(skymap, 'rb'), format='fits')
    #print('Input multi-order skymap:')
    logger.info(skymap.info)

    # Sort by prob density of pixel
    skymap.sort('PROBDENSITY', reverse=True)

    # Get area*probdensity for each pixel
    pixel_area = uniq2pixarea(skymap['UNIQ'])
    #print('Total area = %.1f\n' % (np.sum(pixel_area) * (180/math.pi)**2))

    # Probability per pixel
    prob = pixel_area*skymap['PROBDENSITY']
    cumprob = np.cumsum(prob)

    # Should be 1.0. But need not be.
    sumprob = np.sum(prob)
    logger.info('Sum probability = %.3f\n' % sumprob)

    # Find the index where contour of prob is inside
    i = cumprob.searchsorted(contour*sumprob)
    area_wanted = pixel_area[:i].sum()
    logger.info('Area of %.2f contour is %.2f sq deg' % \
        (contour, area_wanted * (180/math.pi)**2))

    # A MOC is just an astropy Table with one column of healpix indexes
    skymap = skymap[:i]
    skymap = skymap['UNIQ',]
    logger.info(skymap.info)
    skymap.write(outputMOCName, format='fits', overwrite=True)

    # 2023-09-21 KWS There's a bug in MOCpy that requires format to be '1K'
    #                rather than the default 'K' (which is perfectly valid).
    h = fits.open(outputMOCName)
    header = h[1].header
    header['TFORM1'] = '1K'
    h.writeto(outputMOCName, overwrite=True)
    h.close()

    logger.info('MOC file %s written' % outputMOCName)


 
def writeFiles(options, logger):
    for contour in options.contours.split(','):
        skymap = options.skymap
        try:
            c = float(contour)/100.0
            if options.organise:
                os.makedirs(options.directory, exist_ok = True)
                writeMOC(options.skymap, options.directory + '/' + contour + '.moc', c, logger)
            else:
                writeMOC(options.skymap, options.directory + '/' + contour + '.moc', c, logger)
        except ValueError as e:
            logger.error("Contour %s is not a float" % contour)
            print(e)


def main():
    """main.
    """
    opts = docopt(__doc__, version='0.0.10')
    opts = cleanOptions(opts)

    # Use utils.Struct to convert the dict into an object for compatibility with old optparse code.
    options = Struct(**opts)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    writeFiles(options, logger)





if __name__ == '__main__':
    main()
