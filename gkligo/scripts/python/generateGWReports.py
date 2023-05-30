#!/usr/bin/env python
"""Generate GW reports for the various ATLAS and Pan-STARRS databases.

Usage:
  %s <configfile> [--ndays=<ndays>] [--outputlocation=<outputlocation>] [--loadalerts] [--metalocation=<metalocation>] [--recursive] [--significant]
  %s (-h | --help)
  %s --version

Options:
  -h --help                           Show this screen.
  --version                           Show version.
  --ndays=<ndays>                     Number of days [default: 21].
  --outputlocation=<outputlocation>   Location to store the results [default: /tmp].
  --metalocation=<metalocation>       Metafile location [default: /tmp].
  --loadalerts                        Load the Alerts into the database (needs read/write access).
  --recursive                         Check the directory recursively for alerts. 
  --significant                       Only load significant events.

E.g.:
  %s config.yaml --ndays=21
  %s config.yaml --loadalerts --metalocation=/home/ligo/maps_gkligo
"""

import sys
__doc__ = __doc__ % (sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
from docopt import docopt
import os, shutil, re, csv, subprocess
from gkutils.commonutils import Struct, cleanOptions, dbConnect
import MySQLdb
import glob
import yaml

def getFiles(regex, directory, recursive = False):

    directoryRecurse = '/'
    if recursive:
        directoryRecurse = '/**/'
    fileList = glob.glob(directory + directoryRecurse + regex, recursive=recursive)

    fileList.sort()

    return fileList

def parseAlerts(options):

    alerts = []
    if options.recursive:
        yamlFiles = getFiles('meta.yaml', options.metalocation, recursive = True)
    else:
        yamlFiles = getFiles('S*.yaml', options.metalocation, recursive = False)

    for y in yamlFiles:
        with open(y) as yaml_file:
            alert = yaml.safe_load(yaml_file)
            alerts.append(alert)

    return alerts


def flattenAlert(alert):
    flatAlert = {'superevent_id': None,
                 'significant': None,
                 'alert_type': None,
                 'alert_time': None,
                 'mjd_obs': None,
                 'far': None,
                 'distmean': None,
                 'diststd': None,
                 'class_bbh': None,
                 'class_bns': None,
                 'class_nsbh': None,
                 'class_terrestrial': None,
                 'prop_hasns': None,
                 'prop_hasremnant': None,
                 'prop_hasmassgap': None,
                 'area10': None,
                 'area50': None,
                 'area90': None,
                 'creator': None,
                 'group': None,
                 'pipeline': None,
                 'map_iteration': None,
                 'map': None}

    try:
        flatAlert['superevent_id'] = alert['ALERT']['superevent_id']
    except KeyError as e:
        pass

    try:
        flatAlert['alert_type'] = alert['ALERT']['alert_type']
    except KeyError as e:
        pass

    try:
        flatAlert['alert_time'] = alert['ALERT']['time_created']
        if flatAlert['alert_time'] is not None:
            flatAlert['alert_time'] = flatAlert['alert_time'].replace('T', ' ').replace('Z', '')
    except KeyError as e:
        pass

    try:
        flatAlert['map_iteration'] = alert['ALERT']['time_created'].replace('-','').replace(':','').replace('Z','') + '_' + alert['ALERT']['alert_type'].lower()
    except KeyError as e:
        pass

    # Is this a retraction? If so, exit now.

    if flatAlert['alert_type'] == 'RETRACTION':
        return flatAlert

    try:
        flatAlert['significant'] = alert['ALERT']['event']['significant']
    except KeyError as e:
        pass

    try:
        flatAlert['mjd_obs'] = alert['HEADER']['MJD-OBS']
    except KeyError as e:
        pass

    try:
        flatAlert['far'] = alert['ALERT']['event']['far']
    except KeyError as e:
        pass

    try:
        flatAlert['distmean'] = alert['HEADER']['DISTMEAN']
    except KeyError as e:
        pass

    try:
        flatAlert['diststd'] = alert['HEADER']['DISTSTD']
    except KeyError as e:
        pass

    try:
        flatAlert['class_bbh'] = alert['ALERT']['event']['classification']['BBH']
    except KeyError as e:
        pass

    try:
        flatAlert['class_bns'] = alert['ALERT']['event']['classification']['BNS']
    except KeyError as e:
        pass

    try:
        flatAlert['class_nsbh'] = alert['ALERT']['event']['classification']['NSBH']
    except KeyError as e:
        pass

    try:
        flatAlert['class_terrestrial'] = alert['ALERT']['event']['classification']['Terrestrial']
    except KeyError as e:
        pass

    try:
        flatAlert['prop_hasns'] = alert['ALERT']['event']['properties']['HasNS']
    except KeyError as e:
        pass

    try:
        flatAlert['prop_hasremnant'] = alert['ALERT']['event']['properties']['HasRemnant']
    except KeyError as e:
        pass

    try:
        flatAlert['prop_hasmassgap'] = alert['ALERT']['event']['properties']['HasMassGap']
    except KeyError as e:
        pass

    try:
        flatAlert['area10'] = alert['EXTRA']['area10']
    except KeyError as e:
        pass

    try:
        flatAlert['area50'] = alert['EXTRA']['area50']
    except KeyError as e:
        pass

    try:
        flatAlert['area90'] = alert['EXTRA']['area90']
    except KeyError as e:
        pass

    try:
        flatAlert['creator'] = alert['HEADER']['CREATOR']
    except KeyError as e:
        pass

    try:
        flatAlert['group'] = alert['ALERT']['event']['group']
    except KeyError as e:
        pass

    try:
        flatAlert['pipeline'] = alert['ALERT']['event']['pipeline']
    except KeyError as e:
        pass

    try:
        flatAlert['map'] = alert['ALERT']['map']
    except KeyError as e:
        pass


    return flatAlert

def insertAlert(conn, flatAlert):
    import MySQLdb

    try:
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)
        cursor.execute ("""
             insert into tcs_gravity_alerts (`superevent_id`,
                                             `significant`,
                                             `alert_type`,
                                             `alert_time`,
                                             `mjd_obs`,
                                             `far`,
                                             `distmean`,
                                             `diststd`,
                                             `class_bbh`,
                                             `class_bns`,
                                             `class_nsbh`,
                                             `class_terrestrial`,
                                             `prop_hasns`,
                                             `prop_hasremnant`,
                                             `prop_hasmassgap`,
                                             `area10`,
                                             `area50`,
                                             `area90`,
                                             `creator`,
                                             `group`,
                                             `pipeline`,
                                             `map_iteration`,
                                             `map`)
             values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
             """, (flatAlert['superevent_id'],
                   flatAlert['significant'],
                   flatAlert['alert_type'],
                   flatAlert['alert_time'],
                   flatAlert['mjd_obs'],
                   flatAlert['far'],
                   flatAlert['distmean'],
                   flatAlert['diststd'],
                   flatAlert['class_bbh'],
                   flatAlert['class_bns'],
                   flatAlert['class_nsbh'],
                   flatAlert['class_terrestrial'],
                   flatAlert['prop_hasns'],
                   flatAlert['prop_hasremnant'],
                   flatAlert['prop_hasmassgap'],
                   flatAlert['area10'],
                   flatAlert['area50'],
                   flatAlert['area90'],
                   flatAlert['creator'],
                   flatAlert['group'],
                   flatAlert['pipeline'],
                   flatAlert['map_iteration'],
                   flatAlert['map']))


    except MySQLdb.Error as e:
        sys.stderr.write("Error %d: %s\n" % (e.args[0], e.args[1]))

    cursor.close ()
    return conn.insert_id()

# Flatten the alert into key, value pairs so we can load up into the database.
def ingestAlerts(options):
    configFile = options.configfile

    with open(configFile) as yaml_file:
        config = yaml.safe_load(yaml_file)

    usernameatlas = config['databases']['atlas']['username']
    passwordatlas = config['databases']['atlas']['password']
    databaseatlas = config['databases']['atlas']['database']
    hostnameatlas = config['databases']['atlas']['hostname']

    usernameps = config['databases']['ps']['username']
    passwordps = config['databases']['ps']['password']
    databaseps = config['databases']['ps']['database']
    hostnameps = config['databases']['ps']['hostname']

    usernamepso4 = config['databases']['pso4']['username']
    passwordpso4 = config['databases']['pso4']['password']
    databasepso4 = config['databases']['pso4']['database']
    hostnamepso4 = config['databases']['pso4']['hostname']

    connatlas = dbConnect(hostnameatlas, usernameatlas, passwordatlas, databaseatlas)
    connps = dbConnect(hostnameps, usernameps, passwordps, databaseps)
    connpso4 = dbConnect(hostnamepso4, usernamepso4, passwordpso4, databasepso4)

    connatlas.autocommit(True)
    connps.autocommit(True)
    connpso4.autocommit(True)

    alerts = parseAlerts(options)
    for a in alerts:
        flatAlert = flattenAlert(a)
        insertId = insertAlert(connatlas, flatAlert)
        insertId = insertAlert(connps, flatAlert)
        insertId = insertAlert(connpso4, flatAlert)

    connatlas.close()
    connps.close()
    connpso4.close()

def getATLASExposures(conn, options):

    try:
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)

        cursor.execute ("""
            SELECT 
                `dec` AS `decDeg`,
                `texp` AS `exp_time`,
                `filt` AS `filter`,
                `mjd`,
                `ra` AS `raDeg`,
                mag5sig AS `limiting_magnitude`,
                `obs` AS `expname`,
                `obj`
            FROM
                atlas_metadataddc
            WHERE
                mjd > mjdnow() - %s
            ORDER BY mjd DESC
        """, (float(options.ndays),))
        resultSet = cursor.fetchall ()

        cursor.close ()

    except MySQLdb.Error as e:
        print("Error %d: %s" % (e.args[0], e.args[1]))
        sys.exit (1)


    return resultSet


def getPanSTARRSExposures(conn, options, warpstack = 'W'):
    try:
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)

        cursor.execute ("""
            SELECT
                imageid,
                skycell,
                m.exptime exp_time,
                TRUNCATE(mjd_obs, 8) mjd,
                LEFT(fpa_filter, 1) AS filter,
                IF(deteff_counts < 200,
                    m.zero_pt + m.deteff_magref+2.5*log(10,exptime),
                    m.zero_pt + m.deteff_magref + m.deteff_calculated_offset+2.5*log(10,exptime)) AS limiting_mag
            FROM
                tcs_cmf_metadata m
            WHERE
                filename LIKE concat('%%.', %s, 'S.%%')
            AND
                mjd_obs > mjdnow() - %s
            ORDER BY mjd_obs DESC
        """, (warpstack, float(options.ndays),))
        resultSet = cursor.fetchall ()

        cursor.close ()

    except MySQLdb.Error as e:
        print("Error %d: %s" % (e.args[0], e.args[1]))
        sys.exit (1)


    return resultSet


def getCoverage(options):
    configFile = options.configfile

    with open(configFile) as yaml_file:
        config = yaml.safe_load(yaml_file)

    usernameatlas = config['databases']['atlas']['username']
    passwordatlas = config['databases']['atlas']['password']
    databaseatlas = config['databases']['atlas']['database']
    hostnameatlas = config['databases']['atlas']['hostname']

    usernameps = config['databases']['ps']['username']
    passwordps = config['databases']['ps']['password']
    databaseps = config['databases']['ps']['database']
    hostnameps = config['databases']['ps']['hostname']

    usernamepso4 = config['databases']['pso4']['username']
    passwordpso4 = config['databases']['pso4']['password']
    databasepso4 = config['databases']['pso4']['database']
    hostnamepso4 = config['databases']['pso4']['hostname']

    connatlas = dbConnect(hostnameatlas, usernameatlas, passwordatlas, databaseatlas)
    connps = dbConnect(hostnameps, usernameps, passwordps, databaseps)
    connpso4 = dbConnect(hostnamepso4, usernamepso4, passwordpso4, databasepso4)

    atlasExps = getATLASExposures(connatlas, options)
    panstarrsExpsWS = getPanSTARRSExposures(connps, options, warpstack = 'W')
    panstarrsExpsSS = getPanSTARRSExposures(connps, options, warpstack = 'S')
    pso4ExpsWS = getPanSTARRSExposures(connpso4, options, warpstack = 'W')
    pso4ExpsSS = getPanSTARRSExposures(connpso4, options, warpstack = 'S')

    # Write files even if they are empty.
    with open(options.outputlocation + '/' + databaseatlas + 'Exps.csv', 'w') as f:
        if len(atlasExps) > 0:
            w = csv.DictWriter(f, atlasExps[0].keys(), delimiter = ',')
            w.writeheader()
            for row in atlasExps:
                w.writerow(row)

    with open(options.outputlocation + '/' + databaseps + 'WSExps.csv', 'w') as f:
        if len(panstarrsExpsWS) > 0:
            w = csv.DictWriter(f, panstarrsExpsWS[0].keys(), delimiter = ',')
            w.writeheader()
            for row in panstarrsExpsWS:
                w.writerow(row)

    with open(options.outputlocation + '/' + databaseps + 'SSExps.csv', 'w') as f:
        if len(panstarrsExpsSS) > 0:
            w = csv.DictWriter(f, panstarrsExpsSS[0].keys(), delimiter = ',')
            w.writeheader()
            for row in panstarrsExpsSS:
                w.writerow(row)

    with open(options.outputlocation + '/' + databasepso4 + 'WSExps.csv', 'w') as f:
        if len(pso4ExpsWS) > 0:
            w = csv.DictWriter(f, pso4ExpsWS[0].keys(), delimiter = ',')
            w.writeheader()
            for row in pso4ExpsWS:
                w.writerow(row)

    with open(options.outputlocation + '/' + databasepso4 + 'SSExps.csv', 'w') as f:
        if len(pso4ExpsSS) > 0:
            w = csv.DictWriter(f, pso4ExpsSS[0].keys(), delimiter = ',')
            w.writeheader()
            for row in pso4ExpsSS:
                w.writerow(row)

    connatlas.close()
    connps.close()
    connpso4.close()



def main():
    opts = docopt(__doc__, version='0.1')
    opts = cleanOptions(opts)
    options = Struct(**opts)

    if options.loadalerts:
        ingestAlerts(options)
    else:
        getCoverage(options)


if __name__ == '__main__':
    main()
