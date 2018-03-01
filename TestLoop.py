import ftplib
import sys
import json
import logging
import psycopg2

import xiUI_mzid as mzidParser
import xiUI_pg as db


class TestLoop:

    def __init__(self):

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
        ftp = ftplib.FTP(self.ip)
        ftp.login()
        target_dir = self.base
        try:
            ftp.cwd(target_dir)
        except ftplib.error_perm as e:
            error_msg = "%s: %s" % (target_dir, e.args[0])
            print error_msg
            ftp.quit()
            sys.exit(1)

        files = []

        try:
            files = ftp.nlst()
        except ftplib.error_perm, resp:
            if str(resp) == "550 No files found":
                print "No files in this directory"
            else:
                raise
        ftp.quit()

        for f in files:
            self.year(f)

    def year(self, y):
        ftp = ftplib.FTP(self.ip)
        ftp.login()
        target_dir = self.base + '/' + y
        try:
            ftp.cwd(target_dir)
        except ftplib.error_perm as e:
            error_msg = "%s: %s" % (target_dir, e.args[0])
            print error_msg
            ftp.quit()
            sys.exit(1)

        files = []

        try:
            files = ftp.nlst()
        except ftplib.error_perm, resp:
            if str(resp) == "550 No files found":
                print "No files in this directory"
            else:
                raise
        ftp.quit()

        for f in files:
            if f != '2011':
                self.month(y + '/' + f)

    def month(self, ym):
        ftp = ftplib.FTP(self.ip)
        ftp.login()
        target_dir = self.base + '/' + ym
        try:
            ftp.cwd(target_dir)
        except ftplib.error_perm as e:
            error_msg = "%s: %s" % (target_dir, e.args[0])
            print error_msg
            ftp.quit()
            sys.exit(1)

        files = []

        try:
            files = ftp.nlst()
        except ftplib.error_perm, resp:
            if str(resp) == "550 No files found":
                print "No files in this directory"
            else:
                raise
        ftp.quit()

        for f in files:
            self.project(ym + '/' + f)

    def project(self, ymp):
        ftp = ftplib.FTP(self.ip)
        ftp.login()
        target_dir = self.base + '/' + ymp
        try:
            ftp.cwd(target_dir)
        except ftplib.error_perm as e:
            error_msg = "%s: %s" % (target_dir, e.args[0])
            print error_msg
            ftp.quit()
            sys.exit(1)

        files = []

        try:
            files = ftp.nlst()
        except ftplib.error_perm, resp:
            if str(resp) == "550 No files found":
                print "No files in this directory"
            else:
                raise
        ftp.quit()
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
                except mzidParser.MzIdParseException as mzId_error:
                    self.logger.exception(mzId_error.message)
                    try:
                        self.cur.execute("""
                    INSERT INTO uploads (
                        base_dir,
                        filename,
                        upload_error,
                        error_type                        
                    )
                    VALUES (%s, %s, %s, %s)""",
                                    [ymp, f, 'MzIdParseException', mzId_error.message])
                        self.con.commit()

                    except psycopg2.Error as e:
                        raise db.DBException(e.message)

                print(json.dumps(returned_json, indent=4))
                # try:
                #     os.remove(filename)
                # except OSError:
                #     pass
                self.mzId_count = self.mzId_count + 1


test_loop = TestLoop()
# test_loop.allYears()  # no point, starts 2012/12

# test_loop.month('2012/12')
# test_loop.year('2013')
# test_loop.year('2014')
# test_loop.year('2015')
# test_loop.year('2016')
# test_loop.year('2017')
# test_loop.year('2018')

# test_loop.month('2017/03')
test_loop.project('2017/03/PXD004025')
print("mzId count:" + str(test_loop.mzId_count))
