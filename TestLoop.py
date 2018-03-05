import ftplib
import sys
import json
import logging
import psycopg2
import os
import gc

import xiUI_mzid as mzidParser
import xiUI_pg as db


class TestLoop:

    def __init__(self):

        #exclusion list

        # create logger
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(name)s %(message)s')
        self.logger = logging.getLogger(__name__)

        self.ip = "193.62.192.9"
        self.base = "pride/data/archive"
        self.mzId_count = 0
        self.unimod_path = 'obo/unimod.obo'

        self.temp_dir = "/home/col/parser_temp"

        # connect to DB
        try:
            self.con = db.connect('')
            self.cur = self.con.cursor()

        except db.DBException as e:
            self.logger.error(e)
            print(e)
            sys.exit(1)

        # create Database tables
        try:
            db.create_tables(self.cur, self.con)
        except db.DBException as e:
            self.logger.error(e)
            print(e)
            sys.exit(1)

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
            self.project(ym + '/' + f)

    def project(self, ymp):
        target_dir = self.base + '/' + ymp
        files = self.get_file_list(target_dir)
        print ('>> ' + ymp)
        for f in files:
            # print(f)
            # regex = r'(.*)mzid'
            # if re.search(regex, f):
            if f.endswith('mzid') or f.endswith('mzid.gz'):
                print(f)
                path = self.temp_dir + '/' + f
                ftp = ftplib.FTP(self.ip)
                ftp.login()  # Username: anonymous password: anonymous@

                try:
                    ftp.cwd(target_dir)
                    ftp.retrbinary("RETR " + f, open(path, 'wb').write)
                except ftplib.error_perm as e:
                    error_msg = "%s: %s" % (f, e.args[0])
                    self.logger.error(error_msg)
                    # returnJSON['errors'].append({
                    #     "type": "ftpError",
                    #     "message": error_msg,
                    # })
                    # print(json.dumps(returnJSON))
                    sys.exit(1)

                returned_json = {}
                try:
                    returned_json = mzidParser.parse(path, target_dir, self.unimod_path, self.cur, self.con, self.logger)
                except Exception as mzId_error:
                    self.logger.exception(mzId_error)
                    try:
                        self.cur.execute("""
                    INSERT INTO uploads (
                        base_dir,
                        filename,
                        error_type,                        
                        upload_error
                    )
                    VALUES (%s, %s, %s, %s)""", [ymp + '/' + f, f, type(mzId_error).__name__, json.dumps(mzId_error.args)])
                        self.con.commit()

                    except psycopg2.Error as e:
                        raise db.DBException(e.message)

                print(json.dumps(returned_json, indent=4))
                try:
                    os.remove(path)
                    if (path.endswith('.gz')):
                        os.remove(path[0:len(path) - 3])
                except OSError:
                    pass
                self.mzId_count = self.mzId_count + 1
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

test_loop.month('2012/12')
test_loop.year('2013')
test_loop.year('2014')
test_loop.year('2015')
test_loop.year('2016')
test_loop.year('2017')
test_loop.year('2018')

# test_loop.month('2017/04')
# test_loop.project('2017/04/PXD004748')
# test_loop.project('2012/12/PXD000039')
# test_loop.project('2013/09/PRD000647')
# 2016/04/PXD003564 # biggie
# 2016/04/PXD003565 # big 2016/04/PXD003565 2016/04/PXD003566, "67 "68
#  2016/05/PXD002905 ?
# >> 2016/10/PXD003935 ?
# 2016/10/PXD004572
# 2017/05/PXD005403 6gb
print("mzId count:" + str(test_loop.mzId_count))
