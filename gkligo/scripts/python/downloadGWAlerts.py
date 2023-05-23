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
  %s <configFile> <action> [--writeMap] [--writeMOC] [--directory=<directory>] [--contours=<contours>] [--pidfile=<pidfile>] [--logfile=<logfile>] [--daemonErrFile=<deamonErrFile>] [--daemonOutFile=<deamonOutFile>] [--terminal] [--writemeta] [--superevents]
  %s (-h | --help)
  %s --version

where action is start|stop|restart|listen

If action is listen, start in non-daemon mode. Otherwise use daemon mode.

Options:
  -h --help                         Show this screen.
  --version                         Show version.
  --writeMap                        Write the Map file
  --writeMOC                        Write the MOC file (selected contour)
  --directory=<directory>           Directory to where the maps and MOCs will be written [default: /tmp].
  --contours=<contours>             Which MOC contours do you want? Multiple contours should be separated by commas, with no spaces [default: 90]
  --pidfile=<pidfile>               PID file [default: /tmp/ligo.pid]
  --logfile=<logfile>               PID file [default: /tmp/ligo.log]
  --daemonErrFile=<daemonErrFile>   Daemon Error File - for recording unexpected errors [default: /tmp/ligodaemonerr.log].
  --daemonOutFile=<daemonOutFile>   Daemon Out File - for recording unexpected output [default: /tmp/ligodaemonout.log].
  --terminal                        Override the Daemon error and out files and write to the terminal. 
  --writemeta                       Attempt to write the event meta into a file.
  --superevents                     Only deal with superevents. Ignore Mock and Test events.

E.g.:
  %s config.yaml start --directory=/home/atls/ligo --writeMap
  %s config.yaml start --writeMap --writeMOC --directory=/tmp --contours=90
  %s config.yaml start --writeMOC --directory=/tmp --contours=90,50,10
  %s /home/ligo/maps_gkligo/config_downloads.yaml start --directory=/home/ligo/maps_gkligo --writeMap --writeMOC --writemeta --contours=90,50,10 --superevents

"""
import sys
__doc__ = __doc__ % (sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
from docopt import docopt

from gcn_kafka import Consumer
import yaml
import json
import base64
from gkutils.commonutils import Struct, cleanOptions
from io import BytesIO
import daemon
from daemon import pidfile
import signal
import os
import time
import logging

def shutdown(signum, frame):  # signum and frame are mandatory
    """shutdown.

    Args:
        signum:
        frame:
    """
    sys.stdout.write("\nOuch...\n")
    sys.exit(0)


def writeMOC(inputFilePointer, outputMOCName, contour, logger):
    """writeMOC.

    Args:
        inputFilePointer:
        outputMOCName:
        contour:
        logger:
    """
    from astropy.table import Table
    from astropy import units as u
    import astropy_healpix as ah
    import numpy as np
    import math
    from ligo.skymap.moc import uniq2pixarea

    # Read and verify the input
    skymap = Table.read(inputFilePointer, format='fits')
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
    logger.info('MOC file %s written' % outputMOCName)


def writeMeta(options, dataDict, logger):
    #import MySQLdb
    from astropy.io import fits

    eventid = dataDict['superevent_id']
    url = dataDict['urls']['gracedb']
    far = None
    instruments = None

    mjd = None
    distance = None
    distanceStd = None
    classification = None
    properties = None
    significant = None
    time = None
    group = None
    pipeline = None
    alertType = dataDict['alert_type']
    timeCreated = dataDict['time_created']

    try:
        time = dataDict['event']['time']
    except KeyError as e:
        logger.error("Can't find the time info.")

    try:
        group = dataDict['event']['group']
    except KeyError as e:
        logger.error("Can't find the group info.")

    try:
        pipeline = dataDict['event']['pipeline']
    except KeyError as e:
        logger.error("Can't find the pipeline info.")

    try:
        classification = dataDict['event']['classification']
    except KeyError as e:
        logger.error("Can't find the classification info.")

    try:
        significant = dataDict['event']['significant']
    except KeyError as e:
        logger.error("Can't find the significance info.")

    try:
        properties = dataDict['event']['properties']
    except KeyError as e:
        logger.error("Can't find the properties info.")

    eventMeta = {}

    try:
        far = float(dataDict['event']['far'])
        # years = 1/(far*(24*3600*365.2425))
    except ValueError as e:
        logger.error("Can't find the FAR info.")
    except KeyError as e:
        logger.error("Can't find the FAR info.")

    try:
        instruments = dataDict['event']['instruments']
    except KeyError as e:
        logger.error("Can't find the instuments info.")


    # Some info (e.g. distance) only in the FITS file
    h = fits.open(BytesIO(base64.b64decode(dataDict['event']['skymap'])))
    header = h[1].header

    try:
        mjd = header['MJD-OBS']
    except KeyError as e:
        logger.error("The MJD-OBS variable is missing.")

    try:
        distance = header['DISTMEAN']
    except KeyError as e:
        logger.error("The DISTMEAN variable is missing.")

    try:
        distanceStd = header['DISTSTD']
    except KeyError as e:
        logger.error("The DISTSTD variable is missing.")

    # Do NOT write to the database, which could potentially be locked and crash the
    # daemon. Write another script (e.g the reports script) that does that. Just
    # record the metadata for loading.

    eventMeta = {'ALERT': {'event': {'far': far,
                                     'instruments': instruments,
                                     'significant': significant,
                                     'classification': classification,
                                     'properties': properties,
                                     'time': time,
                                     'group': group,
                                     'pipeline': pipeline},
                           'urls': {'gracedb': url},
                           'superevent_id': eventid},
                           'alert_type': alertType,
                           'time_created': timeCreated,
                 'HEADER': {'MJD-OBS': mjd,
                            'DISTMEAN': distance,
                            'DISTSTD': distanceStd},
                }

    return eventMeta

 
def listen(options):
    """listen.

    Args:
        options:
    """
    pid = os.getpid()
    sys.stdout.write("\nListening... PID is %s\n" % pid)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(options.logfile)
    fh.setLevel(logging.INFO)

    formatstr = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(formatstr)

    fh.setFormatter(formatter)

    logger.addHandler(fh)

    with open(options.configFile) as yaml_file:
        config = yaml.safe_load(yaml_file)

    client_id = config['client_id']
    client_secret = config['client_secret']
    topic = config['topics']

    consumer = Consumer(client_id=client_id,
                        client_secret=client_secret)

    consumer.subscribe(['igwn.gwalert'])

    while True:
        # Non zero timeout is required, otherwise consume blocks termination signals.
        for message in consumer.consume(timeout=5):
            dataDict = json.loads(message.value().decode('utf-8'))
            try:
                superEventId = dataDict['superevent_id']
                eventType = superEventId[0] # M, T or S

                alertTimeStamp = dataDict['time_created'].replace(' ','T')
                alertType = dataDict['alert_type']
                alertName = superEventId + '_' + alertType + '_' + alertTimeStamp
                logger.info("Alert Received: %s" % alertName) # Future version of this script will log the output.

                # Act on all events unless we only want superevents. Need to think about the logic!
                if eventType != 'S' and options.superevents:
                    logger.info("Skipping over this event.")
                    # Skip to the next event
                    continue

                # Write the event meta into the databases
                if dataDict['event'] is not None and options.writemeta:
                    meta = writeMeta(options, dataDict, logger)
                    if meta:
                        for k,v in meta.items():
                            logger.info("%s = %s" % (k, str(v)))
                        # Overwrite the superevent info every time a new update arrives.
                        with open(options.directory + '/' + superEventId + '.yaml', 'w') as yamlFile:
                            yamlFile.write(yaml.dump(meta))


                if dataDict['event'] is not None and options.writeMap:
                    skymap = dataDict['event']['skymap']
                    with open(options.directory + '/' + alertName + '.fits', 'wb') as fitsFile:
                        fitsFile.write(base64.b64decode(skymap))

                if dataDict['event'] is not None and options.writeMOC:
                    for contour in options.contours.split(','):
                        skymap = dataDict['event']['skymap']
                        try:
                            c = float(contour)/100.0
                            writeMOC(BytesIO(base64.b64decode(skymap)), options.directory + '/' + alertName + '_' + contour + '.moc', c, logger)
                        except ValueError as e:
                            logger.error("Contour %s is not a float" % contour)


            except KeyError as e:
                logger.error(e)
                pass
            logger.info("")


def startDaemon(options):
    """startDaemon.

    Args:
        options:
    """

    if options.terminal:
        err = sys.stderr
        out = sys.stdout
    else:
        out = open(options.daemonOutFile, 'w')
        err = open(options.daemonErrFile, 'w')

    with daemon.DaemonContext(
        working_directory='/tmp',
        umask=0o002,
        pidfile=pidfile.TimeoutPIDLockFile(options.pidfile),
        stdout=out,
        stderr=err,
        signal_map={ signal.SIGTERM: shutdown }
        ) as context:
        listen(options)


def main():
    """main.
    """
    opts = docopt(__doc__, version='0.0.10')
    opts = cleanOptions(opts)

    # Use utils.Struct to convert the dict into an object for compatibility with old optparse code.
    options = Struct(**opts)


    if options.action not in ['start', 'stop', 'restart', 'listen']:
        sys.stderr.write('Valid options for action are start|stop|restart|listen\n')
        sys.exit(1)

    # First check that the listen option has been chosen. If so, start listening to the Kafka
    # stream immediately, without daemonisation.

    if options.action == 'listen':
        if os.path.exists(options.pidfile):
            with open(options.pidfile, mode='r') as f:
                pid = f.read().strip()
                sys.stderr.write("\nDaemon is already running (PID = %s). Kill the existing daemon first. E.g. use the stop option.\n" % pid)
        else:
            listen(options)
    
    if options.action == 'start':
        if os.path.exists(options.pidfile):
            with open(options.pidfile, mode='r') as f:
                pid = f.read().strip()
                sys.stderr.write("\nDaemon is already running (PID = %s). Stop and restart if you want to restart it.\n" % pid)
        else:
            startDaemon(options)

    if options.action == 'stop':
        if os.path.exists(options.pidfile):
            with open(options.pidfile, mode='r') as f:
                pid = f.read().strip()
            print("Stopping daemon (PID = %s)." % pid)
            os.kill(int(pid), signal.SIGTERM)
        else:
            sys.stderr.write("\nDaemon is not running.\n")


    if options.action == 'restart':
        if os.path.exists(options.pidfile):
            with open(options.pidfile, mode='r') as f:
                pid = f.read().strip()
            print("Stopping daemon (PID = %s)." % pid)
            os.kill(int(pid), signal.SIGTERM)
            time.sleep(5)
            startDaemon(options)
        else:
            startDaemon(options)



if __name__ == '__main__':
    main()
