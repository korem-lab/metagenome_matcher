import errno
import subprocess
from os.path import expandvars, basename
import logging
from subprocess import check_output as sbcheck_output
from gzip import open as opengz
from os.path import splitext, exists
from os import makedirs
import os
import sys
from itertools import zip_longest
import numpy as np
from time import sleep

GZIP_EXT = '.gz'
log_ = logging.getLogger(__name__)

def get_logger(log_level=logging.DEBUG, abs_path=False):
    """Sets up a logger with name set to the calling file.
    :abs_path: If true, will use the absolute path for the caller name."""
    import inspect
    caller_name = inspect.stack()[1].filename
    caller_name = caller_name if abs_path else os.path.splitext(os.path.basename(caller_name))[0]
    log = logging.getLogger(caller_name)
    log.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    log.addHandler(ch)
    return log

def shell_command(com, verbose=False, loglevel=None, no_exceptions = False, expand_vars=True):
    if expand_vars == True:
        com = expandvars(com).replace('$', '\$')
    if loglevel is not None:
        if verbose:
            log_.log(loglevel, com)
        errnout = ''
        try:
            errnout = sbcheck_output(com, shell=True, stderr=subprocess.STDOUT)
        except:
            if errnout != '':
                log_.log(loglevel, errnout)
            if not no_exceptions:
                raise
        else:
            if errnout != '':
                log_.log(loglevel, errnout)
        return errnout.decode(sys.stdout.encoding)
    else:
        log_.debug(com)
        try:
            return sbcheck_output(com, shell=True, stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        except:
            if not no_exceptions:
                raise