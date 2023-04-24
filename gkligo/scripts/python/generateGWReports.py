#!/usr/bin/env python
"""Generate GW reports for the various ATLAS and Pan-STARRS databases.

Usage:
  %s <configfile> <numberOfDays> [--outputlocation=<outputlocation>]
  %s (-h | --help)
  %s --version

Options:
  -h --help                           Show this screen.
  --version                           Show version.
  --outputlocation=<outputlocation>   Location to store the results [default: /tmp].

E.g.:
  %s config.yaml 21
"""

import sys
__doc__ = __doc__ % (sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
from docopt import docopt
import os, shutil, re, csv, subprocess
from gkutils.commonutils import Struct, cleanOptions, dbConnect
import MySQLdb


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
                `obs` AS `expname`
            FROM
                atlas_metadataddc
            WHERE
                mjd > mjdnow() - %s
            ORDER BY mjd DESC
        """, (float(options.numberOfDays),))
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
        """, (warpstack, float(options.numberOfDays),))
        resultSet = cursor.fetchall ()

        cursor.close ()

    except MySQLdb.Error as e:
        print("Error %d: %s" % (e.args[0], e.args[1]))
        sys.exit (1)


    return resultSet





def main():
    opts = docopt(__doc__, version='0.1')
    opts = cleanOptions(opts)
    options = Struct(**opts)

    configFile = options.configfile

    import yaml
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


if __name__ == '__main__':
    main()
