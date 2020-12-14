##
## Methods to manage the config yaml and env 
## eg., startup up logging
##
## APSVIZ, April 2020.
##

import datetime as dt
import pandas as pd
import sys
import os
import yaml
import json
import logging

LOGGER = None


class Utilities:
    """
    Class to manage access to the requisite YAML file (default name=main.yml)
    and manage the environment
    """

    def __init__(self):
        """
        Initialize the Utilities class, set up logging
        """
        global LOGGER
        self.config = self.load_config()

        if LOGGER is None and self.config["DEFAULT"]["LOGGING"]:
            log = self.initialize_logging()
            LOGGER = log
        self.log = LOGGER

    def initialize_logging(self):
        """
        Initialize logging
        """
        # logger = logging.getLogger(__name__)
        logger = logging.getLogger("ADCIRC_support_services")
        log_level = getattr(logging, self.config["DEFAULT"].get('LOGLEVEL', 'DEBUG'))
        logger.setLevel(log_level)

        # LogFile = self.config['LOG_FILE']
        LogFile = 'log'
        # print('Use a log filename of '+LogFile)
        formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(funcName)s : %(module)s : %(name)s : %(message)s ')

        dirname = os.path.dirname(LogFile)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)
        file_handler = logging.FileHandler(LogFile, mode='w')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # logging stream
        # formatter = logging.Formatter('%(asctime)s - %(process)d - %(name)s - %(module)s:%(lineno)d - %(levelname)s - %(message)s')
        # stream_handler = logging.StreamHandler()
        # stream_handler.setFormatter(formatter)
        # logger.addHandler(stream_handler)

        return logger

    def load_config(self):
        yaml_file = os.path.join(os.path.dirname(__file__), "../config/", "main.yml")
        if not os.path.exists(yaml_file):
            raise IOError("Failed to load yaml config file")
        with open(yaml_file, 'r') as stream:
            config = yaml.safe_load(stream)
        self.config = config
        return config

    def print_dict(self, t, s):
        if not isinstance(t, dict) and not isinstance(t, list):
            # pass
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

    def fetchBasedir(self, inconfig, basedirExtra='None'):
        try:
            rundir = os.environ[inconfig.replace('$', '')]  # Yaml call to be subsequently removed
        except:
            print('Chosen basedir invalid: '+str(inconfig['DEFAULT']['RDIR']))
            print('reset to CWD')
            rundir = os.getcwd()
        if basedirExtra is not None:
            rundir = rundir+'/'+basedirExtra
            if not os.path.exists(rundir):
                #print("Create high level Cycle dir space at "+rundir)
                try:
                    #os.mkdir(rundir)
                    os.makedirs(rundir)
                except OSError:
                    sys.error("Creation of the high level run directory %s failed" % rundir)
        return rundir

    def getSubdirectoryFileName(self, basedir, subdir, fname):
        """Check and existence of and construct filenames for
        storing the image data. basedir/subdir/filename 
        subdir is created as needed.
        """
        # print(basedir)
        if not os.path.exists(basedir):
            sys.error("Basepath for dumping files non existent "+str(basedir))
        fulldir = os.path.join(basedir, subdir)
        if not os.path.exists(fulldir):
            # print("Create datastation dir space at "+fulldir)
            try:
                os.mkdir(fulldir)
            except OSError:
                sys.error("Creation of the directory %s failed" % fulldir)
            #else:
            #    print("Successfully created the directory %s " % fulldir)
        return os.path.join(fulldir, fname)

    def read_json_file(self, filepath):
        # Read data from JSON file specified by full path
        data = {}
        with open(filepath, 'r') as fp:
            data = json.load(fp)
        return data

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')


utilities = Utilities()
