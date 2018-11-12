import ftplib
import sys
import json
import logging
import psycopg2
import os
import gc
import shutil
import time

from MzIdParser import MzIdParser
from NumpyEncoder import NumpyEncoder
import PostgreSQL as db


class TestLoop:

    def __init__(self):
        # count parsed id files
        self.mzId_count = 0
        # logging
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s %(name)s %(message)s')
        self.logger = logging.getLogger(__name__)
        # config
        self.ip = "193.62.192.9"
        self.base = "pride/data/archive"
        self.unimod_path = 'obo/unimod.obo'
        self.temp_dir = os.path.expanduser('~') + "/parser_temp/"
        # connect to DB
        # try:
        #     con = db.connect('')
        # except db.DBException as e:
        #     self.logger.error(e)
        #     print(e)
        #     sys.exit(1)

    def all_years(self):
        files = self.get_ftp_file_list(self.base)
        for f in files:
            self.year(f)

    def year(self, year):
        target_dir = self.base + '/' + year
        files = self.get_ftp_file_list(target_dir)
        for f in files:
            self.month(year + '/' + f)

    def month(self, year_month):
        target_dir = self.base + '/' + year_month
        files = self.get_ftp_file_list(target_dir)
        for f in files:
            ymp = year_month + '/' + f
            self.project(ymp)

    def project(self, year_month_project):
        target_dir = self.base + '/' + year_month_project
        files = self.get_ftp_file_list(target_dir)
        print ('>> ' + year_month_project)
        for f in files:
            if f.lower().endswith('mzid') or f.lower().endswith('mzid.gz'):
                print(f)
                self.file(year_month_project, f)
                break

    def file(self, ymp, file_name):
        #  make temp dir, it is entirely removed again at end of this function
        try:
            os.mkdir(self.temp_dir)
        except OSError:
            pass

        path = self.temp_dir + file_name
        target_dir = '/' + self.base + '/' + ymp
        ftp = self.get_ftp_login()

        # fetch mzId file from pride
        try:
            ftp.cwd(target_dir)
            ftp.retrbinary("RETR " + file_name, open(path, 'wb').write)
        except ftplib.error_perm as e:
            ftp.quit()
            error_msg = "%s: %s" % (file_name, e.args[0])
            self.logger.error(error_msg)
            raise e
        ftp.quit()

        # init parser
        try:
            # MzIdParser.MzIdParser(identifications_file, upload_folder, peak_list_folder, db, logger,
            #                       user_id=user_id)
            mzId_parser = MzIdParser(path, self.temp_dir,  self.temp_dir, db, self.logger, 0, origin=ymp)
        except Exception as mzId_error:
            error = json.dumps(mzId_error.args, cls=NumpyEncoder)

            con = db.connect('')
            cur = con.cursor()
            try:
                cur.execute("""
                        INSERT INTO uploads (
                            user_id,
                            origin,
                            filename,
                            error_type,
                            upload_error)
                        VALUES (%s, %s, %s, %s, %s)""",
                            [5, ymp, file_name, type(mzId_error).__name__, error])
                con.commit()
            except psycopg2.Error as e:
                raise db.DBException(e.message)
            con.close()
            return


        # fetch peak list files from pride
        peak_files = mzId_parser.get_peak_list_file_names()
        for peak_file in peak_files:
            # don't download raw files, neater to download everything else even if not supported peak list format
            if not peak_file.endswith('raw'):
                ftp = self.get_ftp_login()
                try:
                    ftp.cwd(target_dir)
                    print('getting ' + peak_file)
                    ftp.retrbinary("RETR " + peak_file,
                                   open(self.temp_dir + peak_file, 'wb').write)
                except ftplib.error_perm as e:
                    print('missing file: ' + peak_file + " (checking for .gz)")
                    #  check for gzipped
                    try:
                        os.remove(self.temp_dir + peak_file)
                        print('getting ' + peak_file + '.gz')
                        # ftp.cwd(target_dir + '/generated/')
                        ftp.retrbinary("RETR " + peak_file + '.gz',
                                       open(self.temp_dir + '/' + peak_file + '.gz', 'wb').write)
                    except ftplib.error_perm as e:
                        ftp.close()
                        print('missing file: ' + peak_file + '.gz')

                        warnings = json.dumps(mzId_parser.warnings, cls=NumpyEncoder)

                        con = db.connect('')
                        cur = con.cursor()
                        try:
                            cur.execute("""
                            UPDATE uploads SET
                                error_type=%s,
                                upload_error=%s,
                                upload_warnings=%s
                            WHERE id = %s""", ["Missing file?", peak_file, warnings, mzId_parser.upload_id])
                            con.commit()
                        except psycopg2.Error as e:
                            raise db.DBException(e.message)
                        con.close()
                        return
                    ftp.close()
                ftp.close()

        # actually parse
        try:
            mzId_parser.parse()
        except Exception as mzId_error:
            self.logger.exception(mzId_error)

            error = json.dumps(mzId_error.args, cls=NumpyEncoder)

            con = db.connect('')
            cur = con.cursor()
            db.write_error(mzId_parser.upload_id, type(mzId_error).__name__, error, cur, con)

        try:
            shutil.rmtree(self.temp_dir)
        except OSError:
            pass
        self.mzId_count = self.mzId_count + 1
        mzId_parser = None
        gc.collect()

    def get_ftp_login(self):
        while True:
            try:
                ftp = ftplib.FTP(self.ip)
                ftp.login()  # Uses password: anonymous@
                return ftp
            except:
                print('FTP fail at '+time.strftime("%c")+'... waiting an hour')
                time.sleep(60 * 60)

    def get_ftp_file_list (self, dir):
        ftp = self.get_ftp_login()
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
        files.reverse()
        return files



test_loop = TestLoop()


# test_loop.year('2018')
# test_loop.year('2017')
# test_loop.year('2016')
# test_loop.year('2015')
# test_loop.year('2014')
# test_loop.year('2013')
# test_loop.month('2012/12')

test_loop.project("2018/11/PXD009010")




# mzML
# test_loop.project("2017/11/PXD007748")
# test_loop.project("2016/11/PXD004785")
# test_loop.project("2016/05/PXD002967")
# test_loop.project("2016/09/PXD004499")
# test_loop.project("2015/06/PXD002045")
# test_loop.project("2017/08/PXD007149")
# test_loop.project("2015/06/PXD002048")
# test_loop.project("2015/06/PXD002047")
# test_loop.project("2014/11/PXD001267")

# 2015/06/PXD002046
# 2014/09/PXD001006
# 2014/09/PXD001000
# 2016/09/PXD002317
# 2014/09/PXD000966
# 2015/06/PXD002044
# 2015/06/PXD002043
# 2015/06/PXD002042
# 2015/06/PXD002041
# 2016/06/PXD004163
# 2015/05/PXD002161
# 2018/01/PXD007913
# 2017/11/PXD006204
# 2015/07/PXD002089
# 2015/07/PXD002088
# 2015/07/PXD002087
# 2015/07/PXD002086
# 2017/07/PXD002901
# 2015/07/PXD002085
# 2017/11/PXD007689
# 2015/07/PXD002084
# 2015/05/PXD002161
# 2015/05/PXD002161
# 2015/07/PXD002083
# 2015/07/PXD002082
# 2015/07/PXD002081
# 2015/07/PXD002080
# 2015/06/PXD002050
# 2015/06/PXD002049

#sim-xl
# test_loop.project("2017/05/PXD006574")
# test_loop.project("2015/02/PXD001677")

#missing file
# test_loop.project("2013/09/PXD000443")

#prob
# test_loop.project("2014/04/PXD000579")

print("mzId count:" + str(test_loop.mzId_count))

# @staticmethod
# def get_pride_info (pxd):
#     time.sleep(1)
#     try:
#         prideAPI = urllib.urlopen('https://www.ebi.ac.uk:443/pride/ws/archive/project/' + pxd).read()
#         pride = json.loads(prideAPI)
#         return pride
#     except Exception:
#         print ("failed to get " + pxd + "from pride api. Will try again in 5 secs.")
#         time.sleep(5)
#         return TestLoop.get_pride_info(pxd)
