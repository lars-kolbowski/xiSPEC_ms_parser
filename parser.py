import json
import sys
import os
import shutil
import logging
import ntpath
from zipfile import BadZipfile
from time import time

try:
    # set working directory
    try:
        abspath = os.path.abspath(__file__)
        dname = os.path.dirname(abspath)
        os.chdir(dname)
    except NameError:
        dname = ''

    # import local files
    import xiSPEC_mzid as mzidParser
    import xiSPEC_csv as csvParser
    from xiSPEC_peakList import unzip_peak_lists

    # logging
    try:
        dev = False
        logFile = dname + "/log/%s_%s.log" % (sys.argv[3], int(time()))

    except IndexError:
        dev = True
        logFile = "log/parser.log"

    try:
        os.remove(logFile)
    except OSError:
        pass
    os.fdopen(os.open(logFile, os.O_WRONLY | os.O_CREAT, 0o777), 'w').close()

    # create logger
    logging.basicConfig(filename=logFile, level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)

except Exception as e:
    print (e)
    sys.exit(1)

try:
    if sys.argv[4] == "pg":
        import xiUI_pg as db
    else:
        import xiSPEC_sqlite as db
except IndexError:
    import xiSPEC_sqlite as db
    
# paths and file names
try:

    unimodPath = 'obo/unimod.obo'

    # development testfiles
    if dev:
        baseDir = "/media/data/xiSPEC_test_files/xiSPEC/"
        # identifications_file = baseDir + 'OpenxQuest_example_added_annotations.mzid'
        # peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"
        peakList_file = baseDir + "B170918_12_Lumos_LK_IN_90_HSA-DSSO-HCD_Rep1.mgf"

        # small mzid dataset
        identifications_file = baseDir + "DSSO_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD_CID-only.mzid"
        # peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"

        # large mzid dataset
        # identifications_file = baseDir + "test/Tmuris_exosomes1.mzid"
        # peakList_file = baseDir + "test/20171027_DDA_JC1.zip"

        #csv file
        # identifications_file = baseDir + "example.csv"

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
    logger.error(e.args[0])
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
    "modifications": [],
    "errors": [],
    "warnings": []
}

# parsing
startTime = time()
try:

    # check for peak list zip file
    peakList_fileName = ntpath.basename(peakList_file)
    if peakList_fileName.lower().endswith('.zip'):
        try:
            unzipStartTime = time()
            logger.info('unzipping start')
            peakList_fileList = unzip_peak_lists(peakList_file)
            logger.info('unzipping done. Time: ' + str(round(time() - unzipStartTime, 2)) + " sec")
        except IOError as e:
            logger.error(e.args[0])
            returnJSON['errors'].append({
                "type": "zipParseError",
                "message": e.args[0],
            })
            print(json.dumps(returnJSON))
            sys.exit(1)
        except BadZipfile as e:
            logger.error(e.args[0])
            returnJSON['errors'].append({
                "type": "zipParseError",
                "message": "Looks something went wrong with the upload! Try uploading again.\n",
            })
            print(json.dumps(returnJSON))
            sys.exit(1)

    else:
        peakList_fileList = [peakList_file]

    # Identification File
    identifications_fileName = ntpath.basename(identifications_file)
    if identifications_fileName.lower().endswith('.mzid'):
        logger.info('parsing mzid start')
        identifications_fileType = 'mzid'
        returnJSON = mzidParser.parse(identifications_file, peakList_fileList, unimodPath, cur,  con, logger)

    elif identifications_fileName.endswith('.csv'):
        logger.info('parsing csv start')
        identifications_fileType = 'csv'
        returnJSON = csvParser.parse(identifications_file, peakList_fileList, cur, con, logger)
        # mgfReader = py_mgf.read(peak_list_file)
        # peakListArr = [pl for pl in mgfReader]

    # delete uploaded files after they have been parsed
    if not dev:
        logger.info('deleting uploaded files')
        shutil.rmtree(upload_folder)

except Exception as e:
    logger.exception(e)
    returnJSON['errors'].append(
        {"type": "Error", "message": e.args[0]})


if len(returnJSON["errors"]) > 0 or len(returnJSON["errors"]) > 0:
    returnJSON['response'] = "%i warning(s) and %i error(s) occurred!" % (len(returnJSON['warnings']), len(returnJSON['errors']))
    for warn in returnJSON['warnings']:
        logger.error(warn)
    for err in returnJSON['errors']:
        logger.error(err)
    returnJSON["log"] = logFile.split('/')[-1]

else:
    returnJSON['response'] = "No errors, smooth sailing!"

if len(returnJSON["errors"]) > 100:
    returnJSON["errors"] = returnJSON["errors"][:100]

print(json.dumps(returnJSON))

if con:
    con.close()
logger.info('all done! Total time: ' + str(round(time() - startTime, 2)) + " sec")
