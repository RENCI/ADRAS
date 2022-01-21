#!/usr/bin/env python

#############################################################
#
# Retain from original utilities, only file IO, URL processing, and logging methods
# Add timing calls
# That can be used by any of the ADCIRC support tools
#
# RENCI 2020
#############################################################

import datetime as dt
import pandas as pd
import sys,os
import yaml
import logging
import json
#from argparse import ArgumentParser

#LOGGER = None


class Utilities:
    """

    """
    def __init__(self):
        """
        Initialize the Utilities class, set up logging
        """
#        global LOGGER
        self.config = self.load_config()

#        if LOGGER is None and self.config["DEFAULT"]["LOGGING"]:
#            log = self.initialize_logging()
#            LOGGER = log
#        self.log = LOGGER

#############################################################
# Logging

    def initialize_logging(self):
        """
        Initialize project logging
        """
        # logger = logging.getLogger(__name__)
        logger = logging.getLogger("adras_services")
        log_level = self.config["DEFAULT"].get('LOGLEVEL', 'DEBUG')
        # log_level = getattr(logging, self.config["DEFAULT"].get('LOGLEVEL', 'DEBUG'))
        logger.setLevel(log_level)

        # LogFile = self.config['LOG_FILE']
        # LogFile = '{}.{}.log'.format(thisDomain, currentdatecycle.cdc)
        LogFile = 'log.hazus'
        formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(funcName)s : %(module)s : %(name)s : %(message)s ')
        dirname = os.path.dirname(LogFile)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)
        file_handler = logging.FileHandler(LogFile, mode='a')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # logging stream
        # formatter = logging.Formatter('%(asctime)s - %(process)d - %(name)s - %(module)s:%(lineno)d - %(levelname)s - %(message)s')
        # stream_handler = logging.StreamHandler()
        # stream_handler.setFormatter(formatter)
        # logger.addHandler(stream_handler)

        return logger


#############################################################
# YAML
    def load_config(self, yaml_file=os.path.join(os.path.dirname(__file__), '..', 'config', 'main.yml')):
        #yaml_file = os.path.join(os.path.dirname(__file__), "../config/", "main.yml")
        if not os.path.exists(yaml_file):
            raise IOError("Failed to load yaml config file {}".format(yaml_file))
        with open(yaml_file, 'r') as stream:
            config = yaml.safe_load(stream)
        # check to see if user has an aws cred file in $HOME
        #awsfile=os.path.join(os.path.expanduser("~"),'aws_adcirc_credentials.csv')
        #if os.path.exists(awsfile):
        #    df = pd.read_csv(awsfile, usecols=['User name', 'Password', 'Access key ID', 'Secret access key', 'Console login link'])
        #    temp = {}
        #    temp['s3'] = df.to_dict(orient='records')
        #    config.update(temp)

        self.config = config
        return config

    def print_dict(self, t, s):
        if not isinstance(t, dict) and not isinstance(t, list):
            if not isinstance(t, float):
                print("\t" * s + str(t))
        else:
            for key in t:
                if not isinstance(t, float) and not isinstance(t, list) \
                        and not isinstance(t, int) and not isinstance(t, unicode):
                    print("\t" * s + str(key))
                if not isinstance(t, list):
                    self.print_dict(t[key], s + 1)

    def serializeMe(self, o):
        if isinstance(o, dt.datetime):
            return o.__str__()

    def readConfigYml(self, yamlfilename):
        if not os.path.exists(yamlfilename):
            raise IOError("Failed to find config file %s" % yamlfilename)
        with open(yamlfilename, 'r') as stream:
            config_file = yaml.safe_load(stream)
        return config_file

#############################################################
# IO uses the base YAML config to do its work

    def fetchBasedir(self, inconfig, basedirExtra=None):
        try:
            rundir = os.environ[inconfig.replace('$', '')]  # Yaml call to be subsequently removed
        except:
            rundir = os.getcwd()
        if basedirExtra is not None:
            rundir = rundir+'/'+basedirExtra
            if not os.path.exists(rundir):
                try:
                    os.makedirs(rundir)
                except OSError:
                    sys.exit("Creation of the high level run directory %s failed" % rundir)
        return rundir

    def setBasedir(self, indir, basedirExtra=None):
        if basedirExtra is not None:
            indir = indir+'/'+basedirExtra
        if not os.path.exists(indir):
            try:
                os.makedirs(indir)
            except OSError:
                sys.exit("Creation of the high level run directory %s failed" % indir)
        return indir

    def getSubdirectoryFileName(self, basedir, subdir, fname ):
        """
        Check and existence of and construct filenames for
        storing the image data. basedir/subdir/filename 
        subdir is created as needed.
        """
        if not os.path.exists(basedir):
            try:
                os.makedirs(basedir)
            except OSError:
                sys.exit("Creation of the basedir %s failed" % basedir)
        fulldir = os.path.join(basedir, subdir)
        if not os.path.exists(fulldir):
            try:
                os.makedirs(fulldir)
            except OSError:
                sys.exit("Creation of the directory %s failed" % fulldir)
            #else:
            #    print("Successfully created the directory %s " % fulldir)
        return os.path.join(fulldir, fname)

    def writePickle(self, df, filename=None):
        """ 
        """
        try:
            df.to_pickle(filename)
        except IOError:
            raise IOError("Failed to write PKL file %s" % (filename))
        return 0

    def writeDictToJson(self, dictdata, rootdir='.',subdir='errorfile',fileroot='filename',iometadata='Nometadata'):
        """
        Write out current self.merged_dict as a Json. Must not use a datetime  as keys
        """
        newfilename=None
        try:
            mdir = rootdir
            newfilename = self.getSubdirectoryFileName(mdir, subdir, fileroot+iometadata+'.json')
            with open(newfilename, 'w') as fp:
                json.dump(dictdata, fp)
            utilities.log.info('Wrote JSON file {}'.format(newfilename))
        except IOError:
            raise IOError("Failed to write file %s" % (newfilename))
        return newfilename

    def read_json_file(self, filepath):
        # Read data from JSON file specified by full path
        data = {}
        try:
            with open(filepath, 'r') as fp:
                data = json.load(fp)
        except FileNotFoundError:
            raise FileNotFoundError("Failed to read file %s" % (filepath))
        return data

    def write_json_file(self, data, filepath):
        # write data from JSON file specified by full path
        try:
            with open(filepath, 'w') as fp:
                json.dump(data, fp)
        except IOError:
            raise IOError("Failed to write JSON file %s" % (filepath))


def validate_url(url):  # TODO
    pass
    # status = True  # Innocent until proven guilty
    # # check for validity
    # if not status:
    #     utilities.log.error("Invalid URL:".format(url))
    # return status


#############################################################
# instance the utilities class on import so that logging is
# immediately available.
utilities = Utilities()
