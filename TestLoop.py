import ftplib
import sys
import json
import logging
import psycopg2
import os
import gc
import shutil

from xiUI_mzid import MzIdParser
from xiUI_mzid import NumpyEncoder
import xiUI_pg as db


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

        # # connect to DB
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
            if f != '2011':
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
            if f.endswith('mzid') or f.endswith('mzid.gz'):
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
                    formats = []
                    for fmat in mzId_parser.spectrum_id_formats:
                        formats.append(fmat)
                    formats = json.dumps(formats, cls=NumpyEncoder)
                    con = db.connect('')
                    cur = con.cursor()
                    try:
                        cur.execute("""
                    UPDATE uploads SET
                        error_type=%s,
                        upload_error=%s,
                        spectrum_id_format=%s
                    WHERE id = %s""", [type(mzId_error).__name__, json.dumps(mzId_error.args), formats, mzId_parser.upload_id])
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
                raise
        ftp.quit()
        return files

test_loop = TestLoop()
# test_loop.allYears()  # no point, starts 2012/12

# # test_loop.month('2012/12')
# # test_loop.year('2013')
# # test_loop.year('2014')
# test_loop.year('2015')
# test_loop.year('2016')
# test_loop.year('2017')
# test_loop.year('2018')

test_loop.project('2015/04/PXD001877')
# test_loop.project('2015/02/PXD001357') - missing scan - how did it pass before

# test_loop.project('2014/09/PXD001054') # - strange one, contains bib ref
# 2014/09/PXD001311
# 2014/10/PXD001403


# 2014/04/PXD000579
# 2014/07/PXD000662
# 2014/07/PXD000710
# 2012/12/PXD000112
# 2013/09/PXD000443
# 2013/12/PXD000623
# 2014/01/PXD000198
# 2014/01/PXD000456
# 2014/04/PXD000521
# 2014/04/PXD000565
# 2014/04/PXD000566
# 2014/04/PXD000567
# 2014/05/PXD000223
# 2014/05/PXD000568
# 2014/09/PXD000966
# 2014/09/PXD001000
# 2014/09/PXD001006
# 2012/12/PXD000039

# test_loop.project('2017/10/PXD004883')
# test_loop.project('2017/04/PXD004748') # no id for DataCollection
# test_loop.project('2012/12/PXD000039') # 1.0.0
# test_loop.project('2017/09/PXD005119') # key error:  PeptideEvidence
# test_loop.project('2017/08/PXD004706') # raw files
# test_loop.project('2017/06/PXD001683') # windows file paths

# 2017/09/PXD007267 xiUI_pg.DBException: integer out of range, should be fixed

# >> 2014/10/PXD001390 # loads of big mgf (fails coz missing scan)

# successful
# 2013/07/PXD000225
# 2014/01/PXD000647
# 2014/02/PXD000202
# 2014/02/PXD000682
# 2014/05/PXD000237
# 2014/05/PXD000807
# 2014/06/PXD000783
# 2014/06/PXD001077
# 2014/07/PXD000923
# 2014/07/PXD001072
# 2014/07/PXD001073
# 2014/10/PXD000961
# 2014/10/PXD001402
# 2014/11/PXD001045
# 2014/11/PXD001089
# 2014/11/PXD001090
# 2014/11/PXD001267
# 2014/11/PXD001386


print("mzId count:" + str(test_loop.mzId_count))
