# gkligo
Code to download GW Alert skymaps from Kafka, and optionally write them
out and/or convert them to MOC files at specified contours.

The code is based on code written by Roy Williams and Leo Singer.

It requires the following python packages:

* gcn-kafka
* healpy
* gkutils
* docopt
* pyYAML
* numpy
* ligo.skymap
* astropy
* mocpy

The command line utilities are:
* downloadGWAlert (Download the alert from Kafka)

The config_example.yaml file is an example config file. Copy it and rename to config.yaml.
Get your credentials by following the first part of the [Kafka Notices via GCN](https://emfollow.docs.ligo.org/userguide/tutorial/receiving/gcn.html) tutorial. Leave the topics configuration unchanged.

TO DO:
* Convert into a Daemon
* Create a log file

