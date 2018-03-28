import ftplib
import sys
import json
import logging
import psycopg2
import os
import urllib
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

        self.temp_dir = os.path.expanduser('~') + "/parser_temp/"

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
        files = self.get_ftp_file_list(self.base)
        for f in files:
            self.year(f)

    def year(self, y):
        target_dir = self.base + '/' + y
        files = self.get_ftp_file_list(target_dir)
        for f in files:
            self.month(y + '/' + f)

    def month(self, ym):
        target_dir = self.base + '/' + ym
        files = self.get_ftp_file_list(target_dir)
        for f in files:
            ymp = ym + '/' + f
            if ymp not in self.exclusion_list:
                self.project(ymp)
            else:
                print('skipping ' + ymp)

    def project(self, ymp):
        pxd = ymp.split('/')[-1]
        prideAPI = urllib.urlopen('https://www.ebi.ac.uk:443/pride/ws/archive/project/' + pxd).read()
        pride = json.loads(prideAPI)

        if pride['submissionType'] == 'COMPLETE':
            target_dir = self.base + '/' + ymp
            files = self.get_ftp_file_list(target_dir)
            print ('>> ' + ymp)
            try:
                os.mkdir(self.temp_dir)
            except OSError:
                pass

            for f in files:
                if f.lower().endswith('mzid') or f.lower().endswith('mzid.gz'):
                    print(f)
                    self.file(ymp, f)

    def file(self, ymp, file_name):
        path = self.temp_dir + '/' + file_name
       #         ftp = ftplib.FTP(self.ip)
        #         ftp.login()  # Username: anonymous password: anonymous@
        #
        #         try:
        #             ftp.cwd(target_dir)
        #             ftp.retrbinary("RETR " + f, open(path, 'wb').write)
        #         except ftplib.error_perm as e:
        #             ftp.quit()
        #             error_msg = "%s: %s" % (f, e.args[0])
        #             self.logger.error(error_msg)
        #             # returnJSON['errors'].append({
        #             #     "type": "ftpError",
        #             #     "message": error_msg,
        #             # })
        #             # print(json.dumps(returnJSON))
        #             sys.exit(1)
        #         ftp.quit()
        base = 'ftp://' + self.ip + '/' + self.base + '/' + ymp + '/'
        urllib.urlretrieve(base + file_name, path)
        mzId_parser = MzIdParser(path, self.temp_dir, db, self.logger)
        peak_files = mzId_parser.get_peak_list_file_names()
        for peak_file in peak_files:
            # urllib.urlretrieve(base + peak_file, self.temp_dir + '/' + peak_file)

            ftp = ftplib.FTP(self.ip)
            ftp.login()  # Uses password: anonymous@
            try:
                target_dir = '/' + self.base + '/' + ymp
                ftp.cwd(target_dir)
                self.logger.info('getting ' + peak_file)
                ftp.retrbinary("RETR " + peak_file,
                               open(self.temp_dir + '/' + peak_file, 'wb').write)
            except ftplib.error_perm as e:
                raise e
                # #  check for gzipped
                # try:
                #     self.logger.info('getting ' + peak_list_file_name + '.gz')
                #     ftp.retrbinary("RETR " + peak_list_file_name + '.gz',
                #                    open(self.temp_dir + '/' + peak_list_file_name + '.gz', 'wb').write)
                # except ftplib.error_perm as e:
                #     ftp.close()
                #     raise e
                # ftp.close()
                # peak_list_file_name = ntpath.basename(
                #     peakListParser.unzip_peak_lists(self.temp_dir + '/' + peak_list_file_name + '.gz')[0])
            ftp.close()

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


    def get_ftp_file_list (self, dir):
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

test_loop.project("2017/05/PXD006574")
test_loop.project("2015/02/PXD001677")

print("mzId count:" + str(test_loop.mzId_count))
