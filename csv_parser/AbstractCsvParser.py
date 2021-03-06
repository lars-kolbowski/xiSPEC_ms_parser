import sys
import numpy as np
from time import time
import pandas as pd
from PeakListParser import PeakListParser
import os
#import pyteomics.fasta as py_fasta
import SimpleFASTA


class CsvParseException(Exception):
    pass


class MissingFileException(Exception):
    pass


class AbstractCsvParser:
    """

    """

    default_values = {
        'rank': 1,
        'pepseq1': '',
        'pepseq2': '',
        'linkpos1': -1,
        'linkpos2': -1,
        'crosslinkermodmass': 0,
        'passthreshold': True,
        'fragmenttolerance': '10 ppm',
        'iontypes': 'peptide;b;y',
        'score': 0,
        'decoy1': -1,
        'decoy2': -1,
        'protein2': '',
        'peppos2': -1,
        'expmz': -1,  # ToDo: required in mzid - also make required col?
        'calcmz': -1
    }

    def __init__(self, csv_path, temp_dir, peak_list_dir, db, logger, db_name='', user_id=0):
        """

        :param csv_path: path to csv file
        :param temp_dir: absolute path to temp dir for unzipping/storing files
        :param db: database python module to use (xiUI_pg or xiSPEC_sqlite)
        :param logger: logger to use
        """

        self.csv_path = csv_path
        self.upload_id = 0
        self.peak_list_readers = {}  # peak list readers indexed by spectraData_ref

        self.temp_dir = temp_dir
        if not self.temp_dir.endswith('/'):
            self.temp_dir += '/'
        self.peak_list_dir = peak_list_dir
        if peak_list_dir and not peak_list_dir.endswith('/'):
            self.peak_list_dir += '/'

        self.user_id = user_id

        self.db = db
        self.logger = logger

        # self.spectra_data_protocol_map = {}
        # ToDo: Might change to pyteomics unimod obo module
        # ToDo: check self.modlist against unimod?
        self.unimod_path = 'obo/unimod.obo'
        self.modlist = []
        self.unknown_mods = []

        self.contains_crosslinks = False
        self.fasta = False
        self.random_id = 0

        self.warnings = []

        # connect to DB
        try:
            self.con = db.connect(db_name)
            self.cur = self.con.cursor()

        except db.DBException as e:
            self.logger.error(e)
            print(e)
            sys.exit(1)

        self.logger.info('reading csv - start')
        self.start_time = time()
        # schema: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0.xsd
        self.csv_reader = pd.read_csv(self.csv_path)

        # check for duplicate columns
        col_list = self.csv_reader.columns.tolist()
        duplicate_cols = set([x for x in col_list if col_list.count(x) > 1])
        if len(duplicate_cols) > 0:
            raise CsvParseException("duplicate column(s): %s" % '; '.join(duplicate_cols))

        self.csv_reader.columns = [x.lower().replace(" ", "") for x in self.csv_reader.columns]
        self.meta_columns = [col for col in self.csv_reader.columns if col.startswith('meta')][:3]

        # remove unused columns
        for col in self.csv_reader.columns:
            if col not in self.required_cols + self.optional_cols + self.meta_columns:
                try:
                    del self.csv_reader[col]
                except KeyError:
                    pass



        # check required cols
        # for required_col in self.required_cols:
        #     if required_col not in self.csv_reader.columns:
        #         raise CsvParseException("Required csv column %s missing" % required_col)

        # create missing non-required cols and fill with NaN (will then be fill with default values)
        for optional_col in self.optional_cols:
            if optional_col not in self.csv_reader.columns:
                self.csv_reader[optional_col] = np.nan

        self.csv_reader.fillna(value=self.default_values, inplace=True)

        # self.csv_reader.fillna('Null', inplace=True)

    def check_required_columns(self):
        for required_col in self.required_cols:
            if required_col not in self.csv_reader.columns:
                raise CsvParseException("Required csv column %s missing" % required_col)
        return True

    def get_missing_required_columns(self):
        missing_cols = []
        for required_col in self.required_cols:
            if required_col not in self.csv_reader.columns:
                missing_cols.append(required_col)
        return missing_cols


    # ToDo: not used atm - can be used for checking if all files are present in temp dir
    def get_peak_list_file_names(self):
        """
        :return: list of all used peak list file names
        """
        return self.csv_reader.peaklistfilename.unique()

    def get_sequenceDB_file_names(self):
        fasta_files = []
        for file in os.listdir(self.temp_dir):
            if file.endswith(".fasta") or file.endswith(".FASTA"):
                fasta_files.append(self.temp_dir + "/" + file)
        return fasta_files

    def set_peak_list_readers(self):
        """
        sets self.peak_list_readers
        dictionary:
            key: peak list file name
            value: associated peak list reader
        """

        peak_list_readers = {}
        for peak_list_file_name in self.csv_reader.peaklistfilename.unique():

            # ToDo: what about .ms2?
            if peak_list_file_name.lower().endswith('.mgf'):
                file_format_accession = 'MS:1001062'        # MGF
                spectrum_id_format_accesion = 'MS:1000774'  # MS:1000774 multiple peak list nativeID format - zero based

            elif peak_list_file_name.lower().endswith('.mzml'):
                file_format_accession = 'MS:1000584'        # mzML
                spectrum_id_format_accesion = 'MS:1001530'  # mzML unique identifier
            else:
                raise CsvParseException("Unsupported peak list file type for: %s" % peak_list_file_name)

            peak_list_file_path = self.peak_list_dir + peak_list_file_name

            try:
                peak_list_reader = PeakListParser(
                    peak_list_file_path,
                    file_format_accession,
                    spectrum_id_format_accesion
                )
            except IOError:
                # try gz version
                try:
                    peak_list_reader = PeakListParser(
                        PeakListParser.extract_gz(peak_list_file_path + '.gz'),
                        file_format_accession,
                        spectrum_id_format_accesion
                    )
                except IOError:
                    # ToDo: output all missing files not just first encountered. Use get_peak_list_file_names()?
                    raise CsvParseException('Missing peak list file: %s' % peak_list_file_name)

            peak_list_readers[peak_list_file_name] = peak_list_reader

        self.peak_list_readers = peak_list_readers

    def parse(self):

        start_time = time()

        # ToDo: more gracefully handle missing files
        if self.peak_list_dir:
            self.set_peak_list_readers()

        self.upload_info() # overridden (empty function) in xiSPEC subclass
        self.parse_db_sequences() # overridden (empty function) in xiSPEC subclass
        self.main_loop()

        meta_col_names = [col.replace("meta_", "") for col in self.meta_columns]
        while len(meta_col_names) < 3:
            meta_col_names.append(-1)
        meta_data = [self.upload_id] + meta_col_names + [self.contains_crosslinks]
        self.db.write_meta_data(meta_data, self.cur, self.con)

        self.logger.info('all done! Total time: ' + str(round(time() - start_time, 2)) + " sec")



    # @staticmethod
    # def get_unimod_masses(unimod_path):
    #     masses = {}
    #     mod_id = -1
    #
    #     with open(unimod_path) as f:
    #         for line in f:
    #             if line.startswith('id: '):
    #                 mod_id = ''.join(line.replace('id: ', '').split())
    #
    #             elif line.startswith('xref: delta_mono_mass ') and not mod_id == -1:
    #                 mass = float(line.replace('xref: delta_mono_mass ', '').replace('"', ''))
    #                 masses[mod_id] = mass
    #
    #     return masses

    def parse_db_sequences(self):
        self.logger.info('reading fasta - start')
        self.start_time = time()
        self.fasta = SimpleFASTA.get_db_sequence_dict(self.get_sequenceDB_file_names())
        self.logger.info('reading fasta - done. Time: ' + str(round(time() - self.start_time, 2)) + " sec")

    def upload_info(self):
        self.logger.info('new csv upload')
        # ident_file_size = os.path.getsize(self.csv_path)
        # peak_list_file_names = json.dumps(self.get_peak_list_file_names(), cls=NumpyEncoder)
        self.upload_id = self.db.new_upload([self.user_id, os.path.basename(self.csv_path), "-"],
                       self.cur, self.con,
                       )
        self.random_id = self.db.get_random_id(self.upload_id, self.cur, self.con)


# class NumpyEncoder(json.JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, np.ndarray):
#             return obj.tolist()
#         return json.JSONEncoder.default(self, obj)
