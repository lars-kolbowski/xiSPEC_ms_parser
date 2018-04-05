import ftplib
import sys
import json
import logging
import psycopg2
import os
import gc
import shutil

from MzIdParser import MzIdParser
from MzIdParser import NumpyEncoder
import PostgreSQL as db


class TestLoop:

    def __init__(self):

        self.exclusion_list = [
            '2016/04/PXD003564',
            '2016/04/PXD003565',
            '2016/04/PXD003566',
            '2016/04/PXD003567',
            '2016/04/PXD003568',
            '2016/05/PXD002905',
            '2016/10/PXD003935',
            '2016/10/PXD004572',
            '2017/05/PXD005403',
            '2017/06/PXD001767'  # big zip
        ]
        # logging
        # try:
        #     dev = False
        #     logFile = dname + "/log/%s_%s.log" % (args[2], int(time()))
        #
        # except IndexError:
        #     dev = True
        #     logFile = "log/parser_%s.log" % int(time())
        #
        # try:
        #     os.remove(logFile)
        # except OSError:
        #     pass
        # os.fdopen(os.open(logFile, os.O_WRONLY | os.O_CREAT, 0o777), 'w').close()

        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(name)s %(message)s')
        self.logger = logging.getLogger(__name__)

        self.ip = "193.62.192.9"
        self.base = "pride/data/archive"
        self.mzId_count = 0
        self.unimod_path = 'obo/unimod.obo'

        self.temp_dir = os.path.expanduser('~') + "/parser_temp"

        # connect to DB
        # try:
        #     con = db.connect('')
        #     cur = con.cursor()
        #
        # except db.DBException as e:
        #     self.logger.error(e)
        #     print(e)
        #     sys.exit(1)
        #
        # # create Database tables
        # try:
        #     db.create_tables(cur, con)
        # except db.DBException as e:
        #     self.logger.error(e)
        #     print(e)
        #     sys.exit(1)
        #
        # con.close

    def all_years(self):
        files = self.get_file_list(self.base)
        for f in files:
            self.year(f)

    def year(self, y):
        target_dir = self.base + '/' + y
        files = self.get_file_list(target_dir)
        for f in files:
            self.month(y + '/' + f)

    def month(self, ym):
        target_dir = self.base + '/' + ym
        files = self.get_file_list(target_dir)
        for f in files:
            ymp = ym + '/' + f
            if ymp not in self.exclusion_list:
                self.project(ymp)
            else:
                print('skipping ' + ymp)

    def project(self, ymp):
        target_dir = self.base + '/' + ymp
        files = self.get_file_list(target_dir)
        print ('>> ' + ymp)
        try:
            os.mkdir(self.temp_dir)
        except OSError:
            pass

        for f in files:
            if f.lower().endswith('mzid') or f.lower().endswith('mzid.gz'):
                print(f)
                path = self.temp_dir + '/' + f
                ftp = ftplib.FTP(self.ip)
                ftp.login()  # Username: anonymous password: anonymous@

                try:
                    ftp.cwd(target_dir)
                    ftp.retrbinary("RETR " + f, open(path, 'wb').write)
                except ftplib.error_perm as e:
                    ftp.quit()
                    error_msg = "%s: %s" % (f, e.args[0])
                    self.logger.error(error_msg)
                    # returnJSON['errors'].append({
                    #     "type": "ftpError",
                    #     "message": error_msg,
                    # })
                    # print(json.dumps(returnJSON))
                    sys.exit(1)
                ftp.quit()

                mzId_parser = MzIdParser(path, self.temp_dir, ymp, db, self.ip, self.base, self.logger)
                try:
                    mzId_parser.parse()
                except Exception as mzId_error:
                    self.logger.exception(mzId_error)

                    error = json.dumps(mzId_error.args, cls=NumpyEncoder)

                    spec_id_formats = ','.join(mzId_parser.spectrum_id_formats)
                    file_formats = ','.join(mzId_parser.file_formats)

                    warnings = json.dumps(mzId_parser.warnings, cls=NumpyEncoder)

                    peak_list_files = json.dumps(mzId_parser.peak_list_file_names, cls=NumpyEncoder)

                    con = db.connect('')
                    cur = con.cursor()
                    try:
                        cur.execute("""
                    UPDATE uploads SET
                        error_type=%s,
                        upload_error=%s,
                        spectrum_id_format=%s,
                        file_format=%s,
                        upload_warnings=%s,
                        peak_list_file_names=%s
                    WHERE id = %s""", [type(mzId_error).__name__, error, spec_id_formats, file_formats, warnings, peak_list_files, mzId_parser.upload_id])
                        con.commit()

                    except psycopg2.Error as e:
                        raise db.DBException(e.message)
                    con.close()

                try:
                    shutil.rmtree(self.temp_dir)
                except OSError:
                    pass
                self.mzId_count = self.mzId_count + 1
                mzId_parser = None
                gc.collect()
                break

    def get_file_list (self, dir):
        ftp = ftplib.FTP(self.ip)
        ftp.login()
        try:
            ftp.cwd(dir)
        except ftplib.error_perm as e:
            error_msg = "%s: %s" % (dir, e.args[0])
            print error_msg
            ftp.quit()
            return []

        files = []

        try:
            files = ftp.nlst()
        except ftplib.error_perm, resp:
            if str(resp) == "550 No files found":
                print "No files in this directory"
            else:
                error_msg = "%s: %s" % (dir, ftplib.error_perm.args[0])
                print error_msg

        ftp.quit()
        return files

test_loop = TestLoop()

# test_loop.project("2014/07/PXD000710")  # MS:1000776: scan number only nativeID format -  -  ["requested scanID 38581 not found in peakList file"]
# test_loop.project("2014/09/PXD000966")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["ParseError", ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=154"]
# test_loop.project("2014/09/PXD001000")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["ParseError", ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=726"]
# test_loop.project("2014/09/PXD001006")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["ParseError", ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=250"]
# test_loop.project("2014/10/PXD001034")  # MS:1000774: multiple peak list nativeID format - MS:1001062: Mascot MGF file - ["invalid spectrum ID format!"]
# test_loop.project("2015/02/PXD001213")  # MS:1000774: multiple peak list nativeID format - MS:1001062: Mascot MGF file - ["invalid spectrum ID format!"]
# test_loop.project("2015/03/PXD000719")  # MS:1000774: multiple peak list nativeID format - MS:1001062: Mascot MGF file - ["requested scanID 280934 not found in peakList file"]

# missing FragmentTolerance
# test_loop.project("2015/05/PXD002161")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=8819"]

test_loop.project("2016/11/PXD004785")  # MS:1001530:  mzML unique identifier - MS:1000584: mzML file

# test_loop.project("2015/06/PXD002041")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=456"]
# test_loop.project("2015/06/PXD002042")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=262"]
# test_loop.project("2015/06/PXD002043")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=356"]
# test_loop.project("2015/06/PXD002044")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=221"]
# test_loop.project("2015/06/PXD002045")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=672"]
# test_loop.project("2015/06/PXD002046")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=364"]
# test_loop.project("2015/06/PXD002047")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=480"]
# test_loop.project("2015/06/PXD002048")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=388"]
# test_loop.project("2015/06/PXD002049")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=302"]
# test_loop.project("2015/06/PXD002050")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=17"]
# test_loop.project("2015/07/PXD002080")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=443"]
# test_loop.project("2015/07/PXD002081")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=297"]
# test_loop.project("2015/07/PXD002082")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=537"]
# test_loop.project("2015/07/PXD002083")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=303"]
# test_loop.project("2015/07/PXD002084")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=318"]
# test_loop.project("2015/07/PXD002085")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=384"]
# test_loop.project("2015/07/PXD002086")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=773"]
# test_loop.project("2015/07/PXD002087")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=234"]
# test_loop.project("2015/07/PXD002088")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=320"]
# test_loop.project("2015/07/PXD002089")  # MS:1000768: Thermo nativeID format - MS:1000584: mzML file - ["failed to parse spectrumID from controllerType=0 controllerNumber=1 scan=214"]
# test_loop.project("2016/02/PXD001376")  # MS:1000774: multiple peak list nativeID format - MS:1001062: Mascot MGF file - ["invalid spectrum ID format!"]
# test_loop.project("2016/02/PXD002963")  # MS:1000774: multiple peak list nativeID format - MS:1001062: Mascot MGF file - ["invalid spectrum ID format!"]
# test_loop.project("2016/04/PXD003232")  # MS:1000774: multiple peak list nativeID format - MS:1001062: Mascot MGF file - ["invalid spectrum ID format!"]
# test_loop.project("2016/05/PXD001953")  # MS:1000774: multiple peak list nativeID format - MS:1001062: Mascot MGF format - ["requested scanID 207443 not found in peakList file"]
# test_loop.project("2016/05/PXD002967")  # MS:1000774: multiple peak list nativeID format - MS:1000584: mzML format - ["invalid spectrum ID format!"]


# test_loop.allYears()  # no point, starts 2012/12

# test_loop.month('2012/12')
# test_loop.year('2013')
# test_loop.year('2014')
# test_loop.year('2015')
# test_loop.year('2016')

# test_loop.month('2017/08')
# test_loop.month('2017/09')
# test_loop.month('2017/10')

# >> 2017/10/PXD005168 crashed here

# test_loop.month('2017/11')
# test_loop.month('2017/12')
# # test_loop.year('2017')
# test_loop.year('2018')




# test_loop.project('2016/05/PXD002967')  # missing version

# normal
test_loop.project('2015/04/PXD001885')

# pyteomics lib had prob reading SpectraData?
# test_loop.project('2014/09/PXD001006')
# test_loop.project('2014/09/PXD001000')
# test_loop.project('2014/09/PXD000966')
# test_loop.project('2014/07/PXD000710')
# test_loop.project('2016/02/PXD001376')
# test_loop.project('2012/12/PXD000112')  # empty mod name
# test_loop.project('2017/06/PXD006757')
# test_loop.project('2018/01/PXD005859')  # empty mod name

# test_loop.project('2017/10/PXD004883')
# test_loop.project('2017/04/PXD004748')  # no id for DataCollection
# test_loop.project('2014/09/PXD001054')  # contains bib ref
# test_loop.project('2012/12/PXD000039')  # 1.0.0
# test_loop.project('2017/09/PXD005119')  # key error:  PeptideEvidence
# test_loop.project('2017/08/PXD004706')  # raw files
# test_loop.project('2017/06/PXD001683')  # windows file paths
# test_loop.project('2017/09/PXD007267')  # xiUI_pg.DBException: integer out of range, should be fixed
# header stuff that can mess up mgf reader
# test_loop.project('2016/02/PXD001997')
# test_loop.project('2016/03/PXD002078')
# test_loop.project('2016/01/PXD002855')
# test_loop.project('2015/04/PXD001877')
# test_loop.project('2015/02/PXD001357')

# index=null
# test_loop.project('2016/02/PXD001376')
# test_loop.project('2015/02/PXD001213')
# test_loop.project('2016/02/PXD002963')
# test_loop.project('2014/10/PXD001034')
# test_loop.project('2015/03/PXD000719')
# test_loop.project('2016/04/PXD003232')

# missing scan
# test_loop.project('2014/10/PXD001034')
# test_loop.project('2015/02/PXD001213')
# test_loop.project('2015/03/PXD000719')
# test_loop.project('2016/02/PXD001376')

# missing files
# test_loop.project('2012/12/PXD000112')
# test_loop.project('2013/09/PXD000443')
# test_loop.project('2013/12/PXD000623')
# test_loop.project('2014/01/PXD000198')
# test_loop.project('2014/01/PXD000456')
# test_loop.project('2014/04/PXD000521')
# test_loop.project('2014/04/PXD000565')
# test_loop.project('2014/04/PXD000566')
# test_loop.project('2014/04/PXD000567')
# test_loop.project('2014/04/PXD000579')
# test_loop.project('2014/05/PXD000223')
# test_loop.project('2014/05/PXD000568')
# test_loop.project('2014/07/PXD000662')
# test_loop.project('2015/05/PXD000941')
# test_loop.project('2015/05/PXD000942')
#
# # were missing file errors
# test_loop.project('2017/01/PXD004764')
# test_loop.project('2017/01/PXD004778')
# test_loop.project('2017/01/PXD004788')
# test_loop.project('2017/01/PXD004796')
# test_loop.project('2016/01/PXD003445')
# test_loop.project('2016/03/PXD002759')
# test_loop.project('2016/03/PXD003132')

# out of memory (or struggles)
# test_loop.project('2015/05/PXD002117')
# test_loop.project('2014/11/PXD001422')
# test_loop.project('2016/04/PXD002417')
# test_loop.project('2018/01/PXD006308')

# wiff peak data
# test_loop.project('2016/01/PXD003445')


# >> 2014/10/PXD001390 # loads of big mgf (fails coz missing scan)

print("mzId count:" + str(test_loop.mzId_count))
