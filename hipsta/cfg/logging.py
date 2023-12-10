"""

"""

import logging
import os
import sys
import tempfile
import time

from .version import get_version

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS


def setup_logging(args):
    """
    Set up logging
    """

    # check if output directory exists or can be created
    if os.path.isdir(args.outputdir):
        print("Found output directory " + args.outputdir)
    else:
        try:
            os.mkdir(args.outputdir)
        except OSError as e:
            print("Cannot create output directory " + args.outputdir)
            print("Reason: " + str(e))
            raise

    # check if logfile can be written in output directory
    try:
        testfile = tempfile.TemporaryFile(dir=args.outputdir)
        testfile.close()
    except OSError as e:
        print(args.outputdir + " not writeable")
        print("Reason: " + str(e))
        raise

    # setup logging
    logfile =  os.path.join(args.outputdir, 'logfile.txt')
    logfile_format = "[%(levelname)s: %(filename)s] %(message)s"
    logfile_handlers_stream = logging.StreamHandler(sys.stdout)
    logfile_handlers_stream.setFormatter(logging.Formatter(logfile_format))
    logfile_handlers_stream.setLevel(logging.INFO)
    logfile_handlers_file = logging.FileHandler(filename=logfile, mode="w")
    logfile_handlers_file.setFormatter(logging.Formatter(logfile_format))
    logfile_handlers_file.setLevel(logging.INFO)

    logging.basicConfig(level=logging.INFO, format=logfile_format, handlers=[logfile_handlers_file, logfile_handlers_stream])

    # intial messages (do not use logging earlier, since this will do the
    # default setup, and subsequent custom setup will have no effects)
    LOGGER.info("Starting logging ...")
    LOGGER.info("Logfile: %s", logfile)
    LOGGER.info("Version: %s", get_version())
    LOGGER.info("Date: %s", time.strftime('%d/%m/%Y %H:%M:%S'))

    # log args
    LOGGER.info("Command: " + " ".join(sys.argv))