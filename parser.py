import pyteomics.mzid as py_mzid
# import pyteomics.mzml as py_mzml
import json
import sys
import os
import shutil
import logging
import xiSPEC_mzid as mzidParser
import ntpath


try:
    # set working directory
    try:
        abspath = os.path.abspath(__file__)
        dname = os.path.dirname(abspath)
        os.chdir(dname)
    except NameError:
        dname = ''

    # logging
    try:
        dev = False
        logFile = dname + "/log/" + sys.argv[3] + ".log"

    except IndexError:
        dev = True
        logFile = "log/parser.log"

    if not os.path.isfile(logFile):
        os.fdopen(os.open(logFile, os.O_WRONLY | os.O_CREAT, 0o777), 'w').close()

    # create logger
    logging.basicConfig(filename=logFile, level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)

except Exception as e:
    print e
    sys.exit(1)


# paths and file names
try:
    import xiSPEC_sqlite as db

    unimodPath = 'obo/unimod.obo'

    # development testfiles
    if dev:
        baseDir = "/home/lars/work/xiSPEC/"
        identifications_file = baseDir + "DSSO_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD_CID-only.mzid"
        # mzidFile = baseDir + 'OpenxQuest_example_added_annotations.mzid'
        peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"
        # peakList_file = baseDir + "B170918_12_Lumos_LK_IN_90_HSA-DSSO-HCD_Rep1.mgf"
        dbName = 'test.db'

    else:
        identifications_file = sys.argv[1]
        peakList_file = sys.argv[2]
        upload_folder = "../uploads/" + sys.argv[3]
        dbfolder = "dbs/tmp/"
        try:
            os.stat(dbfolder)
        except:
            os.mkdir(dbfolder)
        dbName = dbfolder + sys.argv[3] + '.db'
except Exception as e:
    logger.error(e)
    print(e)
    sys.exit(1)

# connect to DB
try:
    con = db.connect(dbName)
    cur = con.cursor()

except db.DBException as e:
    logger.error(e)
    print(e)
    sys.exit(1)

# create Database tables
try:
    db.create_tables(cur, con)
except db.DBException as e:
    logger.error(e)
    print(e)
    sys.exit(1)


returnJSON = {
    "response": "",
    "errors": []
}

# parsing
try:
    logger.info('reading mzid - start')
    # schema: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0.xsd
    mzidReader = py_mzid.MzIdentML(identifications_file)
    logger.info('reading mzid - done')

    mzidParser.parse(mzidReader, peakList_file, unimodPath, cur,  con, logger)

    # delete uploaded files after they have been parsed
    if not dev:
        logger.info('deleting uploaded files')
        shutil.rmtree(upload_folder)

    if len(returnJSON["errors"]) > 0:
        returnJSON['response'] = "Warning: %i error(s) occured!" % len(
            returnJSON['errors'])
        for e in returnJSON['errors']:
            logger.error(e)
    else:
        returnJSON['response'] = "No errors, smooth sailing!"

    print(json.dumps(returnJSON))
    if con:
        con.close()
        logger.info('all done!')

except Exception as e:
    logger.exception(e)
    returnJSON['errors'].append(
        {"type": "Error", "message": e})
