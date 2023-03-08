#!/usr/bin/env python
"""
Download GW Alerts (multi-order skymap) and convert to MOC file if required.
https://emfollow.docs.ligo.org/userguide/tutorial/multiorder_skymaps.html

This code written by Ken Smith based on code by Leo Singer and Roy Williams 

Requires:
numpy
ligo.skymap
astropy
mocpy

Usage:
  %s <configFile> [--writeMap] [--writeMOC] [--directory=<directory>] [--contours=<contours>]
  %s (-h | --help)
  %s --version

Options:
  -h --help                    Show this screen.
  --version                    Show version.
  --writeMap                   Write the Map file
  --writeMOC                   Write the MOC file (selected contour)
  --directory=<directory>      Directory to where the maps and MOCs will be written [default: /tmp].
  --contours=<contours>        Which MOC contours do you want? Multiple contours should be separated by commas, with no spaces [default: 90]

E.g.:
  %s config.yaml --directory=/home/atls/ligo --writeMap
  %s config.yaml --writeMap --writeMOC --directory=/tmp --contours=90
  %s config.yaml --writeMOC --directory=/tmp --contours=90,50,10

"""
import sys
__doc__ = __doc__ % (sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
from docopt import docopt

from gcn_kafka import Consumer
import json
import base64
from gkutils.commonutils import Struct, cleanOptions
from io import BytesIO

def writeMOC(inputFilePointer, outputMOCName, contour):
    from astropy.table import Table
    from astropy import units as u
    import astropy_healpix as ah
    import numpy as np
    import math
    from ligo.skymap.moc import uniq2pixarea

    # Read and verify the input
    skymap = Table.read(inputFilePointer, format='fits')
    #print('Input multi-order skymap:')
    print(skymap.info)

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
    print('Sum probability = %.3f\n' % sumprob)

    # Find the index where contour of prob is inside
    i = cumprob.searchsorted(contour*sumprob)
    area_wanted = pixel_area[:i].sum()
    print('Area of %.2f contour is %.2f sq deg' % \
        (contour, area_wanted * (180/math.pi)**2))

    # A MOC is just an astropy Table with one column of healpix indexes
    skymap = skymap[:i]
    skymap = skymap['UNIQ',]
    print(skymap.info)
    skymap.write(outputMOCName, format='fits', overwrite=True)
    print('MOC file %s written' % outputMOCName)


def main():
    opts = docopt(__doc__, version='0.1')
    opts = cleanOptions(opts)

    # Use utils.Struct to convert the dict into an object for compatibility with old optparse code.
    options = Struct(**opts)

    import yaml
    with open(options.configFile) as yaml_file:
        config = yaml.safe_load(yaml_file)

    client_id = config['client_id']
    client_secret = config['client_secret']
    topic = config['topics']

    # Connect as a consumer.
    # Warning: don't share the client secret with others.
    consumer = Consumer(client_id=client_id,
                        client_secret=client_secret)

    # Subscribe to topics and receive alerts
    consumer.subscribe(['igwn.gwalert'])

    while True:
        for message in consumer.consume():
            dataDict = json.loads(message.value().decode('utf-8'))
            #for k, v in dataDict.items():
            #    print(k, v)
            try:
                superEventId = dataDict['superevent_id']
                alertTimeStamp = dataDict['time_created'].replace(' ','T')
                alertType = dataDict['alert_type']
                alertName = superEventId + '_' + alertType + '_' + alertTimeStamp
                print("Alert Received: %s" % alertName) # Future version of this script will log the output.
                if dataDict['event'] is not None and options.writeMap:
                    skymap = dataDict['event']['skymap']
                    with open(options.directory + '/' + alertName + '.fits', 'wb') as fitsFile:
                        fitsFile.write(base64.b64decode(skymap))

                if dataDict['event'] is not None and options.writeMOC:
                    for contour in options.contours.split(','):
                        skymap = dataDict['event']['skymap']
                        try:
                            c = float(contour)/100.0
                            writeMOC(BytesIO(base64.b64decode(skymap)), options.directory + '/' + alertName + '_' + contour + '.moc', c)
                        except ValueError as e:
                            print("Contour %s is not a float" % contour)


            except KeyError as e:
                print(e)
                pass
            print()

if __name__ == '__main__':
    main()
