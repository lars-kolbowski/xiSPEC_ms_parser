import json
import sys
import os
import shutil
import logging
import ntpath
from zipfile import BadZipfile
from time import time
import re
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], "f", [])
except getopt.GetoptError:
    # print 'test.py -i <inputfile> -o <outputfile>'
    sys.exit(2)

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
        logFile = dname + "/log/%s_%s.log" % (args[2], int(time()))

    except IndexError:
        dev = True
        logFile = "log/parser_%s.log" % int(time())

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
    if args[3] == "pg":
        import xiUI_pg as db
    else:
        import xiSPEC_sqlite as db
except IndexError:
    import xiSPEC_sqlite as db


returnJSON = {
    "response": "",
    "modifications": [],
    "errors": [],
    "warnings": [],
    "log": logFile.split('/')[-1]
}


# paths and file names
try:

    unimodPath = 'obo/unimod.obo'

    # development testfiles
    if dev:
        baseDir = "/media/data/work/xiSPEC_test_files/"
        # identifications_file = baseDir + 'OpenxQuest_example_added_annotations.mzid'
        # peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"
        # peakList_file = baseDir + "B170918_12_Lumos_LK_IN_90_HSA-DSSO-HCD_Rep1.mgf"

        # identifications_file = "/media/data/work/xiSPEC_test_files/SL/02_wogroups.mzid"
        # peakList_file = "/media/data/work/xiSPEC_test_files/SL/mscon_PF_20_100_0_B160803_02_new.mgf"

        # # mzid has duplicate ids!!! - fixed now with non-flat index
        # identifications_file = "/media/data/work/xiSPEC_test_files/PXD006767/MTases_Trypsin_ETD_search.mzid"
        # peakList_file = "/media/data/work/xiSPEC_test_files/PXD006767/PXD006767.zip"

        # small mzid dataset
        identifications_file = baseDir + "DSSO_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD_CID-only.mzid"
        peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"

        # # large mzid dataset
        # identifications_file = baseDir + "Tmuris_exo/Tmuris_exosomes1.mzid"
        # peakList_file = baseDir + "Tmuris_exo/20171027_DDA_JC1.zip"

        #csv file
        # identifications_file = baseDir + "example.csv"

        dbName = 'test.db'

    else:
        if '-f' in [o[0] for o in opts]:
            import ftplib

            upload_folder = "../uploads/%s/" % int(time())
            try:
                os.stat(upload_folder)
            except:
                os.mkdir(upload_folder)

            id_file_path = "/".join(args[0].split("/")[3:-1])
            id_file_path = "/%s/" % id_file_path
            id_file_name = args[0].split("/")[-1]
            identifications_file = upload_folder + id_file_name

            pl_file_path = "/".join(args[1].split("/")[3:-1])
            pl_file_path = "/%s/" % pl_file_path
            pl_file_name = args[1].split("/")[-1]
            peakList_file = upload_folder + pl_file_name

            ftp = ftplib.FTP('ftp.pride.ebi.ac.uk')
            ftp.login()

            try:
                ftp.cwd(id_file_path)
                ftp.retrbinary("RETR " + id_file_name, open(identifications_file, 'wb').write)
            except ftplib.error_perm as e:
                error_msg = "%s: %s" % (id_file_name, e.args[0])
                logger.error(error_msg)
                returnJSON['errors'].append({
                    "type": "ftpError",
                    "message": error_msg,
                })
                print(json.dumps(returnJSON))
                sys.exit(1)

            try:
                ftp.cwd(pl_file_path)
                ftp.retrbinary("RETR " + pl_file_name, open(peakList_file, 'wb').write)
                ftp.quit()
            except ftplib.error_perm as e:
                error_msg = "%s: %s" % (pl_file_name, e.args[0])
                logger.error(error_msg)
                returnJSON['errors'].append({
                    "type": "ftpError",
                    "message": error_msg,
                })
                print(json.dumps(returnJSON))
                sys.exit(1)
        else:

            identifications_file = args[0]
            peakList_file = args[1]
            upload_folder = "../uploads/" + args[2]

        dbfolder = "dbs/tmp/"
        try:
            os.stat(dbfolder)
        except:
            os.mkdir(dbfolder)
        dbName = dbfolder + args[2] + '.db'
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

# parsing
startTime = time()
try:
    # check for peak list zip file
    peakList_fileName = ntpath.basename(peakList_file)
    if re.search(".*\.(gz|zip)$", peakList_fileName):
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
    if re.match(".*\.mzid(\.gz)?$", identifications_fileName):
        logger.info('parsing mzid start')
        identifications_fileType = 'mzid'
        id_returnJSON = mzidParser.parse(identifications_file, peakList_fileList, unimodPath, cur,  con, logger)
        returnJSON.update(id_returnJSON)

    elif identifications_fileName.endswith('.csv'):
        logger.info('parsing csv start')
        identifications_fileType = 'csv'
        id_returnJSON = csvParser.parse(identifications_file, peakList_fileList, cur, con, logger)
        returnJSON.update(id_returnJSON)
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


if len(returnJSON["errors"]) > 0 or len(returnJSON["warnings"]) > 0:
    returnJSON['response'] = "%i warning(s) and %i error(s) occurred!" % (len(returnJSON['warnings']), len(returnJSON['errors']))
    for warn in returnJSON['warnings']:
        logger.error(warn)
    for err in returnJSON['errors']:
        logger.error(err)

else:
    returnJSON['response'] = "No errors, smooth sailing!"

if len(returnJSON["errors"]) > 100:
    returnJSON["errors"] = returnJSON["errors"][:100]

print(json.dumps(returnJSON))

if con:
    con.close()
logger.info('all done! Total time: ' + str(round(time() - startTime, 2)) + " sec")
