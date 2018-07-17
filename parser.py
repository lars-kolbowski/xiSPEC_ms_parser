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


dev = False
use_ftp, use_postgreSQL, user_id = False, False, False
identifications_file, peakList_file, identifier = False, False, False

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:p:s:u:", ["ftp", "postgresql"])
except getopt.GetoptError:
    print('parser.py (-f -pg) -i <identifications file> -p <peak list file> -s <session identifier>')
    sys.exit(2)

for o, a in opts:
    if o in ("-f", "--ftp"):
        import ftplib
        use_ftp = True

    if o == "-i":
        identifications_file = a

    if o == "-p":
        peakList_file = a

    if o == "-s":
        identifier = a

    if o == '--postgresql':
        use_postgreSQL = True

    if o == '-u': # user_id
        user_id = a


if identifications_file is False or identifier is False:
    dev = True
    print ("dev test mode...")

if use_postgreSQL:
    import PostgreSQL as db
else:
    import SQLite as db

if use_ftp:
    import ftplib

try:
    # set working directory
    try:
        abspath = os.path.abspath(__file__)
        dname = os.path.dirname(abspath)
        os.chdir(dname)
    except NameError:
        dname = ''

    # import local files
    import MzIdParser
    import CsvParser
    import PeakListParser

    # logging
    logFile = dname + "/log/%s_%s.log" % (identifier, int(time()))
    if not dev:
        try:
            os.remove(logFile)
        except OSError:
            pass
        os.fdopen(os.open(logFile, os.O_WRONLY | os.O_CREAT, 0o777), 'w').close()

        # create logger
        logging.basicConfig(filename=logFile, level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(name)s %(message)s')
    else:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(name)s %(message)s')


    logger = logging.getLogger(__name__)



except Exception as e:
    print (e)
    sys.exit(1)

logger.info('argv:' + " ".join(sys.argv))

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

    if dev:
        # development testfiles
        baseDir = "/media/data/work/xiSPEC_test_files/"

        # identifications_file = baseDir + 'OpenxQuest_example_added_annotations.mzid'
        # peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"
        # peakList_file = baseDir + "B170918_12_Lumos_LK_IN_90_HSA-DSSO-HCD_Rep1.mgf"

        # identifications_file = "/media/data/work/xiSPEC_test_files/SL/02_wogroups.mzid"
        # peakList_file = "/media/data/work/xiSPEC_test_files/SL/mscon_PF_20_100_0_B160803_02_new.mgf"

        # # mzid has duplicate ids!!! - fixed now with non-flat index
        # identifications_file = "/media/data/work/xiSPEC_test_files/PXD006767/MTases_Trypsin_ETD_search.mzid"
        # peakList_file = "/media/data/work/xiSPEC_test_files/PXD006767/PXD006767.zip"
        # peakList_file = "/media/data/work/xiSPEC_test_files/PXD006767/as.zip"

        # # HSA-BS3 dataset
        # identifications_file = baseDir + "/cross-link/xiFDR/E171207_15_Lumos_AB_DE_160_VI186_B1_xiFDR_1.0.23.48/E171207_15_Lumos_AB_DE_160_VI186_B1.mzid"
        # peakList_file = baseDir + "/cross-link/xiFDR/E171207_15_Lumos_AB_DE_160_VI186_B1_xiFDR_1.0.23.48/E171207_15_Lumos_AB_DE_160_VI186_B1.mzML"

        # small mzid dataset
        # identifications_file = baseDir + "DSSO_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD_CID-only.mzid"
        # peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"

        # # large mzid dataset
        # identifications_file = baseDir + "Tmuris_exo/Tmuris_exosomes1.mzid"
        # peakList_file = baseDir + "Tmuris_exo/20171027_DDA_JC1.zip"

        # SL
        # identifications_file = baseDir + "/cross-link/xiFDR/pc_revision/B160803_02_26_57.mzid"
        # peakList_file = baseDir + "/cross-link/xiFDR/pc_revision/B160803_02_Lumos_LK_IN_190_PC_BS3_HCD_DT_1.mzML"

        # PXD006574
        # identifications_file = baseDir + "PXD006574/monomerResults.mzid.gz"
        # peakList_file = baseDir + "PXD006574/monomerResults-specId.pride.mgf.gz"
        # identifications_file = baseDir + "PXD006574/dimerResultsToPRIDE.mzid.gz"
        # peakList_file = baseDir + "PXD006574/dimerResultsToPRIDE-specId.pride.mgf.gz"

        # # PXD001677 - ms2 peak list file
        # identifications_file = baseDir + "cross-link/PXD001677/result_DynamicDBReduction_Plus_Report_Ions-specId.pride.mgf.gz"
        # peakList_file = baseDir + "cross-link/PXD001677/result_DynamicDBReduction-specId.ms2"

        # PXD007836 - mzid 1.1.0
        # identifications_file = baseDir + "PXD007836/data.mzid"
        # peakList_file = baseDir + "PXD007836/c.zip"

        # csv file
        identifications_file = "/home/col/test2/result_new_XiVersion1.6.739_PSM_xiFDR1.1.27.csv"
        peakList_file = "/home/col/test2/Rappsilber_CLMS_PolII_mgfs.zip"

        database = 'test.db'
        upload_folder = "/".join(identifications_file.split("/")[:-1]) + "/"

    else:

        database = "dbs/tmp/%s.db" % identifier
        upload_folder = "../uploads/" + identifier

        if use_ftp:

            upload_folder = "../uploads/%s/" % int(time())
            try:
                os.stat(upload_folder)
            except:
                os.mkdir(upload_folder)

            id_file_path = "/".join(identifications_file.split("/")[3:-1])
            id_file_path = "/%s/" % id_file_path
            id_file_name = identifications_file.split("/")[-1]
            identifications_file = upload_folder + id_file_name

            pl_file_path = "/".join(peakList_file.split("/")[3:-1])
            pl_file_path = "/%s/" % pl_file_path
            pl_file_name = peakList_file.split("/")[-1]
            peakList_file = upload_folder + pl_file_name

            try:
                ftp = ftplib.FTP('ftp.pride.ebi.ac.uk')
                ftp.login()
            # ToDO: more specific except clause
            except:
                error_msg = "general ftp connection error! Please try again later."
                logger.error(error_msg)
                returnJSON['errors'].append({
                    "type": "ftpError",
                    "message": error_msg,
                })
                print(json.dumps(returnJSON))
                sys.exit(1)

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


except Exception as e:
    logger.error(e.args[0])
    print(e)
    sys.exit(1)


# parsing
startTime = time()
try:
    peak_list_folder = upload_folder
    if peakList_file.endswith('.zip'):
        try:
            unzipStartTime = time()
            logger.info('unzipping start')
            # peakList_fileList = peakListParser.PeakListParser.unzip_peak_lists(peakList_file)
            peak_list_folder = PeakListParser.PeakListParser.unzip_peak_lists(peakList_file)
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

    identifications_fileName = ntpath.basename(identifications_file)
    if re.match(".*\.mzid(\.gz)?$", identifications_fileName):
        logger.info('parsing mzid start')
        identifications_fileType = 'mzid'
        if use_postgreSQL:
            id_parser = MzIdParser.MzIdParser(identifications_file, upload_folder, peak_list_folder, db, logger,
                                              user_id=user_id)
        else:
            id_parser = MzIdParser.xiSPEC_MzIdParser(identifications_file, upload_folder, peak_list_folder, db, logger,
                                              db_name=database)

    elif identifications_fileName.endswith('.csv'):
        logger.info('parsing csv start')
        identifications_fileType = 'csv'
        if use_postgreSQL:
            id_parser = CsvParser.CsvParser(identifications_file, upload_folder, peak_list_folder, db, logger,
                                            user_id=user_id)
        else:
            id_parser = CsvParser.xiSPEC_CsvParser(identifications_file, upload_folder, peak_list_folder, db, logger,
                                            db_name=database)


    else:
        raise Exception('Unknown identifications file format!')


    # create Database tables
    if not use_postgreSQL:
        try:
            db.create_tables(id_parser.cur, id_parser.con)
        except db.DBException as e:
            logger.error(e)
            print(e)
            sys.exit(1)

    id_parser.parse()

    returnJSON['identifier'] = str(id_parser.upload_id) + "-" + str(id_parser.random_id)
    returnJSON['modifications'] = id_parser.unknown_mods
    returnJSON["warnings"] = id_parser.warnings

    # delete uploaded files after they have been parsed
    if not dev:
        logger.info('deleting uploaded files')
        shutil.rmtree(upload_folder)

except Exception as e:
    # print(e)
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

print(json.dumps(returnJSON, indent=4))

# if con:
#     con.close()
logger.info('all done! Total time: ' + str(round(time() - startTime, 2)) + " sec")
