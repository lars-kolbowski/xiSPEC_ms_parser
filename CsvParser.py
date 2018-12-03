import re
import json
import sys
import numpy as np
from time import time
import pandas as pd
from PeakListParser import PeakListParser
import zipfile
import gzip
import os
#import pyteomics.fasta as py_fasta
import SimpleFASTA

class CsvParseException(Exception):
    pass


class MissingFileException(Exception):
    pass


class CsvParser:
    """

    """
    # ToDo: adjust to xiUI needs
    required_cols = [
        'scanid',
        'charge',
        'pepseq1',
        'protein1',
        'peppos1',
        'peaklistfilename',
        # 'expMZ'
    ]

    optional_cols = [
        # 'spectrum_id' $ ToDo: get rid of this? select alternatives by scanid and peaklistfilename?
        'rank',
        'fragmenttolerance',
        'iontypes',
        'pepseq2',
        'linkpos1',
        'linkpos2',
        'crosslinkermodmass',
        'passthreshold',
        'score',
        'decoy1',
        'decoy2',
        'protein2',
        'peppos2',
        'expmz',    # ToDo: required in mzid - also make required col?
        'calcmz'
    ]

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

        self.upload_id = 0
        if csv_path.endswith('.gz') or csv_path.endswith('.zip'):
            self.csv_path = CsvParser.extract_csv(csv_path)
        else:
            self.csv_path = csv_path
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
        for required_col in self.required_cols:
            if required_col not in self.csv_reader.columns:
                raise CsvParseException("Required csv column %s missing" % required_col)

        # create missing non-required cols and fill with NaN (will then be fill with default values)
        for optional_col in self.optional_cols:
            if optional_col not in self.csv_reader.columns:
                self.csv_reader[optional_col] = np.nan

        self.csv_reader.fillna(value=self.default_values, inplace=True)

        # self.csv_reader.fillna('Null', inplace=True)

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

    # ToDo: split into two functions?
    @staticmethod
    def extract_csv(archive):
        if archive.endswith('zip'):
            zip_ref = zipfile.ZipFile(archive, 'r')
            unzip_path = archive + '_unzip/'
            zip_ref.extractall(unzip_path)
            zip_ref.close()

            return_file_list = []

            for root, dir_names, file_names in os.walk(unzip_path):
                file_names = [f for f in file_names if not f[0] == '.']
                dir_names[:] = [d for d in dir_names if not d[0] == '.']
                for file_name in file_names:
                    os.path.join(root, file_name)
                    if file_name.lower().endswith('.mzid'):
                        return_file_list.append(root+'/'+file_name)
                    else:
                        raise IOError('unsupported file type: %s' % file_name)

            if len(return_file_list) > 1:
                raise StandardError("more than one csv file found!")

            return return_file_list[0]

        elif archive.endswith('gz'):
            in_f = gzip.open(archive, 'rb')
            archive = archive.replace(".gz", "")
            out_f = open(archive, 'wb')
            out_f.write(in_f.read())
            in_f.close()
            out_f.close()

            return archive

        else:
            raise StandardError('unsupported file type: %s' % archive)

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
        # pyteomics.fasta is too strict in what headers it accepts (must be one of https://www.uniprot.org/help/fasta-headers)
        # for file in self.get_sequenceDB_file_names():
           # fasta_iterator = py_fasta.read(self.temp_dir + "/" + file)
           # for (a, b) in fasta_iterator:
           #     # self.logger.info("" + a  b)
           #     header = py_fasta.parse(a)
           #     self.fasta[header['id']] = b
        self.fasta = SimpleFASTA.get_db_sequence_dict(self.get_sequenceDB_file_names())
        self.logger.info('reading fasta - done. Time: ' + str(round(time() - self.start_time, 2)) + " sec")

    def main_loop(self):
        main_loop_start_time = time()
        self.logger.info('main loop - start')

        peptide_evidences = []
        spectrum_identifications = []
        spectra = []
        peptides = []
        proteins = set()
        # list of spectra that were already seen - index in list is spectrum_id
        # combination of peaklistfilename and scanid is a unique identifier
        seen_spectra = []

        # list of peptides that were already seen - index in list is peptide_id
        # pep sequence including cross-link site and cross-link mass is unique identifier
        seen_peptides = []

        cross_linker_pair_count = 0

        # # ID VALIDITY CHECK - unique ids
        # if len(self.csv_reader['id'].unique()) < len(self.csv_reader):
        #     duplicate_ids = self.csv_reader[self.csv_reader.duplicated('id', keep=False)].id.tolist()
        #     duplicate_ids = [str(i) for i in duplicate_ids]
        #     raise CsvParseException('Duplicate ids found: %s' % "; ".join(duplicate_ids))

        for identification_id, id_item in self.csv_reader.iterrows():  # identification_id, id_item = id_df.iterrows().next()

            # 1 based row number
            row_number = identification_id + 1

            #
            # VALIDITY CHECKS & TYPE CONVERSIONS - ToDo: move type checks/conversions to col level in parse()?
            #
            # rank - ToDo: more elaborate checks?
            try:
                rank = int(id_item['rank'])
            except KeyError:
                rank = 1    # ToDo: create default values in parse()
            except ValueError:
                raise CsvParseException('Invalid rank: %s for row: %s' % (id_item['rank'], row_number))

            # pepSeq

            # ToDo: reorder peptides by length and alphabetical?
            # add cross-linker always to first peptide?
            # From mzIdentML schema 1.2.0:
            # the cross-link donor SHOULD contain the complete mass delta introduced by the cross-linking reagent,
            # and that the cross-link acceptor reports a mass shift
            # delta of zero. It is RECOMMENDED that the 'donor' peptide SHOULD be the longer peptide, followed by
            # alphabetical order for equal length peptides.

            invalid_char_pattern_pepseq = '([^GALMFWKQESPVICYHRNDTXa-z:0-9(.)\-]+)'
            # pepSeq - 1
            if id_item['pepseq1'] == '':
                raise CsvParseException('Missing PepSeq1 for row: %s' % row_number)

            invalid_char_match = re.match(invalid_char_pattern_pepseq, id_item['pepseq1'])
            if invalid_char_match:
                invalid_chars = "; ".join(invalid_char_match.groups())
                raise CsvParseException(
                    'Invalid character(s) found in PepSeq1: %s for row: %s' % (invalid_chars, row_number)
                )
            pepseq1 = id_item['pepseq1']
            # pepSeq - 2
            if id_item['pepseq2'] == '':
                cross_linked_id_item = False
            else:
                self.contains_crosslinks = True
                cross_linked_id_item = True
                invalid_char_match = re.match(invalid_char_pattern_pepseq, id_item['pepseq2'])
                if invalid_char_match:
                    invalid_chars = "; ".join(invalid_char_match.groups())
                    raise CsvParseException(
                        'Invalid character(s) found in PepSeq2: %s for row: %s' % (invalid_chars, row_number)
                    )
            pepseq2 = id_item['pepseq2']

            # LinkPos
            # LinkPos - 1
            try:
                linkpos1 = int(id_item['linkpos1'])
            except ValueError:
                raise CsvParseException('Invalid LinkPos1: %s for row: %s' % (id_item['linkpos1'], row_number))

            # LinkPos - 2
            try:
                linkpos2 = int(id_item['linkpos2'])
            except ValueError:
                raise CsvParseException('Invalid LinkPos2: %s for row: %s' % (id_item['linkpos2'], row_number))

            if (linkpos1 == -1 and not linkpos2 == -1) or (linkpos1 == -1 and not linkpos2 == -1):
                raise CsvParseException('Incomplete cross-link site information for row: %s' % row_number)

            # CrossLinkerModMass
            try:
                cross_link_mod_mass = float(id_item['crosslinkermodmass'])
            except ValueError:
                raise CsvParseException('Invalid CrossLinkerModMass: %s for row: %s' % (id_item['crosslinkermodmass'], row_number))

            # charge
            try:
                charge = int(id_item['charge'])
            except ValueError:
                raise CsvParseException('Invalid charge state: %s for row: %s' % (id_item['charge'], row_number))

            # passthreshold
            if isinstance(id_item['passthreshold'], bool):
                pass_threshold = id_item['passthreshold']
            else:
                raise CsvParseException('Invalid passThreshold value: %s for row: %s' % (id_item['passthreshold'], row_number))

            # fragmenttolerance
            if not re.match('^([0-9.]+) (ppm|Da)$', str(id_item['fragmenttolerance'])):
                raise CsvParseException(
                    'Invalid FragmentTolerance value: %s in row: %s' % (id_item['fragmenttolerance'], row_number))
            else:
                fragment_tolerance = id_item['fragmenttolerance']

            # iontypes
            ions = id_item['iontypes'].split(';')
            valid_ions = [
                'peptide',
                'a',
                'b',
                'c',
                'x',
                'y',
                'z',
                '' # split will add an empty sell if string ends with ';'
            ]
            if any([True for ion in ions if ion not in valid_ions]):
                raise CsvParseException(
                    'Unsupported IonType in: %s in row %s! Supported ions are: peptide;a;b;c;x;y;z.'
                    % (id_item['iontypes'], row_number)
                )
            ion_types = id_item['iontypes']

            # score
            try:
                score = float(id_item['score'])
            except ValueError:
                raise CsvParseException('Invalid score: %s in row %s' % (id_item['score'], row_number))

            # protein1
            protein_list1 = id_item['protein1'].split(";")
            protein_list1 = [s.strip() for s in protein_list1]
            for p in protein_list1:
                proteins.add(p)

            # decoy1 - if decoy1 is not set fill list with default value (0)
            if id_item['decoy1'] == -1:
                is_decoy_list1 = [False] * len(protein_list1)
            else:
                is_decoy_list1 = []
                for decoy in str(id_item['decoy1']).split(";"):
                    if decoy.lower().strip() == 'true':
                        is_decoy_list1.append(True)
                    elif decoy.lower().strip() == 'false':
                        is_decoy_list1.append(False)
                    else:
                        raise CsvParseException(
                            'Invalid value in Decoy 1: %s in row %s. Allowed values: True, False.'
                            % (id_item['decoy1'], row_number)
                        )

            # pepPos1 - if pepPos1 is not set fill list with default value (-1)
            # ToDo: might need changing for xiUI where pepPos is not optional
            if id_item['peppos1'] == -1:
                pep_pos_list1 = [-1] * len(protein_list1)
            else:
                pep_pos_list1 = str(id_item['peppos1']).split(";")
                pep_pos_list1 = [s.strip() for s in pep_pos_list1]

            # protein - decoy - pepPos sensibility check
            if not len(protein_list1) == len(is_decoy_list1):
                raise CsvParseException(
                    'Inconsistent number of protein to decoy values for Protein1 and Decoy1 in row %s!' % row_number)
            if not len(protein_list1) == len(pep_pos_list1):
                raise CsvParseException(
                    'Inconsistent number of protein to pepPos values for Protein1 and PepPos1 in row %s!' % row_number)

            # protein2
            protein_list2 = id_item['protein2'].split(";")
            protein_list2 = [s.strip() for s in protein_list2]
            for p in protein_list2:
                proteins.add(p)

            # decoy2 - if decoy2 is not set fill list with default value (0)
            if id_item['decoy2'] == -1:
                is_decoy_list2 = [False] * len(protein_list2)
            else:
                is_decoy_list2 = []
                for decoy in str(id_item['decoy2']).split(";"):
                    if decoy.lower().strip() == 'true':
                        is_decoy_list2.append(True)
                    elif decoy.lower().strip() == 'false':
                        is_decoy_list2.append(False)
                    else:
                        raise CsvParseException(
                            'Invalid value in Decoy 2: %s in row %s. Allowed values: True, False.'
                            % (id_item['decoy2'], row_number)
                        )

            # pepPos2 - if pepPos2 is not set fill list with default value (-1)
            # ToDo: might need changing for xiUI where pepPos is not optional
            if id_item['peppos2'] == -1:
                pep_pos_list2 = [-1] * len(protein_list2)
            else:
                pep_pos_list2 = str(id_item['peppos2']).split(";")
                pep_pos_list2 = [s.strip() for s in pep_pos_list2]

            # protein - decoy - pepPos sensibility check
            if not len(protein_list2) == len(is_decoy_list2):
                raise CsvParseException(
                    'Inconsistent number of protein to decoy values for Protein2 and Decoy2 in row %s!' % row_number)
            if not len(protein_list2) == len(pep_pos_list2):
                raise CsvParseException(
                    'Inconsistent number of protein to pepPos values for Protein2 and PepPos2! in row %s!' % row_number)

            # scanId
            try:
                scan_id = int(id_item['scanid'])
            except ValueError:
                raise CsvParseException('Invalid scanid: %s in row %s' % (id_item['scanid'], row_number))

            # peakListFilename

            # expMZ
            try:
                exp_mz = float(id_item['expmz'])
            except ValueError:
                raise CsvParseException('Invalid expMZ: %s in row %s' % (id_item['exmpmz'], row_number))
            # calcMZ
            try:
                calc_mz = float(id_item['calcmz'])
            except ValueError:
                raise CsvParseException('Invalid calcMZ: %s in row %s' % (id_item['calcmz'], row_number))

            #
            # -----Start actual parsing------
            #
            # SPECTRA
            peak_list_file_name = id_item['peaklistfilename']

            unique_spec_identifier = "%s-%s" % (peak_list_file_name, scan_id)

            if unique_spec_identifier not in seen_spectra:
                seen_spectra.append(unique_spec_identifier)
                spectrum_id = len(seen_spectra) - 1
                peak_list = None
                precursor_mz = None
                precursor_charge = None
                if self.peak_list_dir:
                    # get peak list
                    try:
                        peak_list_reader = self.peak_list_readers[peak_list_file_name]
                    except KeyError:
                        raise CsvParseException('Missing peak list file: %s' % peak_list_file_name)

                    scan = peak_list_reader.get_scan(scan_id)
                    peak_list = scan['peaks']
                    precursor_mz = scan['precursor']['mz']
                    precursor_charge = scan['precursor']['charge']

                spectrum = [
                    spectrum_id,                    # 'id',
                    peak_list,                      # 'peak_list',
                    peak_list_file_name,            # 'peak_list_file_name',
                    scan_id,                        # 'scan_id',
                    fragment_tolerance,             # 'frag_tol',
                    self.upload_id,                 # 'upload_id',
                    'Spec_%s' % spectrum_id,        # 'spectrum_ref'
                    precursor_mz,                   # 'precursor_mz',
                    precursor_charge,               # 'precursor_charge'
                ]
                spectra.append(spectrum)
            else:
                spectrum_id = seen_spectra.index(unique_spec_identifier)

            #
            # PEPTIDES
            if cross_linked_id_item:
                cross_linker_pair_id = cross_linker_pair_count
                cross_linker_pair_count += 1
            else:
                cross_linker_pair_id = -1  # linear ToDo: -1 or None?

            # peptide - 1
            unique_pep_identifier1 = "%s-%s" % (pepseq1, cross_linker_pair_id)

            if unique_pep_identifier1 not in seen_peptides:
                seen_peptides.append(unique_pep_identifier1)
                pep1_id = len(seen_peptides) - 1

                peptide1 = [
                    pep1_id,                        # id,
                    pepseq1,                        # seq_mods,
                    linkpos1,                       # link_site,
                    cross_link_mod_mass,            # crosslinker_modmass, declare peptide 1 as cl donor: full mass
                    self.upload_id,                 # upload_id,
                    cross_linker_pair_id            # crosslinker_pair_id
                ]
                peptides.append(peptide1)
            else:
                pep1_id = seen_peptides.index(unique_pep_identifier1)

            if cross_linked_id_item:
                # peptide - 2
                unique_pep_identifier2 = "%s-%s" % (pepseq2, cross_linker_pair_id)

                if unique_pep_identifier2 not in seen_peptides:
                    seen_peptides.append(unique_pep_identifier2)
                    pep2_id = len(seen_peptides) - 1
                    peptide2 = [
                        pep2_id,                        # id,
                        pepseq2,                        # seq_mods,
                        linkpos2,                       # link_site,
                        0,                              # crosslinker_modmass, declare peptide 2 as cl acceptor: 0 mass
                        self.upload_id,                 # upload_id,
                        cross_linker_pair_id            # crosslinker_pair_id
                    ]
                    peptides.append(peptide2)
                else:
                    pep2_id = seen_peptides.index(unique_pep_identifier2)
            else:
                pep2_id = None

            #
            # PEPTIDE EVIDENCES
            # peptide evidence - 1
            for i in range(len(protein_list1)):

                pep_evidence1 = [
                    pep1_id,                # peptide_ref
                    protein_list1[i],       # dbsequence_ref - ToDo: might change to numerical id
                    protein_list1[i],       # protein_accession
                    pep_pos_list1[i],       # pep_start
                    is_decoy_list1[i],      # is_decoy
                    self.upload_id          # upload_id
                ]

                peptide_evidences.append(pep_evidence1)

            if cross_linked_id_item:
                # peptide evidence - 2

                if pep2_id is None:
                    raise StandardError('Fatal! peptide id error!')

                for i in range(len(protein_list2)):

                    pep_evidence2 = [
                        pep2_id,                # peptide_ref
                        protein_list2[i],       # dbsequence_ref - ToDo: might change to numerical id
                        protein_list2[i],       # protein_accession
                        pep_pos_list2[i],       # pep_start
                        is_decoy_list2[i],      # is_decoy
                        self.upload_id          # upload_id
                    ]

                    peptide_evidences.append(pep_evidence2)

            #
            # SPECTRUM IDENTIFICATIONS
            # ToDo: experimental_mass_to_charge, calculated_mass_to_charge
            scores = json.dumps({'score': score})

            try:
                meta1 = id_item[self.meta_columns[0]]
            except IndexError:
                meta1 = ""
            try:
                meta2 = id_item[self.meta_columns[1]]
            except IndexError:
                meta2 = ""
            try:
                meta3 = id_item[self.meta_columns[2]]
            except IndexError:
                meta3 = ""

            spectrum_identification = [
                identification_id,          # 'id',
                self.upload_id,             # 'upload_id',
                spectrum_id,                # 'spectrum_id',
                pep1_id,                    # 'pep1_id',
                pep2_id,                    # 'pep2_id',
                charge,                     # 'charge_state',
                rank,                       # 'rank',
                pass_threshold,             # 'pass_threshold',
                ion_types,                  # 'ions',
                scores,                     # 'scores',
                exp_mz,                     # 'experimental_mass_to_charge',
                calc_mz,                    # 'calculated_mass_to_charge'
                meta1,
                meta2,
                meta3
            ]
            spectrum_identifications.append(spectrum_identification)

            #
            # MODIFICATIONS
            # ToDo: check against unimod?

            try:
                modifications = re.findall('[^A-Z]+', ''.join([pepseq1, pepseq2]))
            except AttributeError:
                modifications = []

            for mod in modifications:
                if mod not in self.unknown_mods:
                    self.unknown_mods.append(mod)


        # DBSEQUENCES
        if self.fasta:
            db_sequences = []
            for prot in proteins:
               try:
                   #data = [prot] + self.fasta[prot] + [self.upload_id]
                   temp = self.fasta[prot]
                   data = [prot, temp[0], temp[1], temp[2], temp[3], self.upload_id] # surely there's a better way
               except Exception as ke:
                   seq = "NO SEQUENCE"
                   data = [prot, prot, prot, "", "NO SEQUENCE", self.upload_id]

               # is_decoy - not there
               # data.append("false")

               # data.append(self.upload_id)

               db_sequences.append(data)


        # end main loop
        self.logger.info('main loop - done. Time: ' + str(round(time() - main_loop_start_time, 2)) + " sec")

        # once loop is done write data to DB
        db_wrap_up_start_time = time()
        self.logger.info('write spectra to DB - start')
        try:

            self.db.write_peptide_evidences(peptide_evidences, self.cur, self.con)
            self.db.write_peptides(peptides, self.cur, self.con)
            self.db.write_spectra(spectra, self.cur, self.con)
            self.db.write_spectrum_identifications(spectrum_identifications, self.cur, self.con)
            if self.fasta:
                self.db.write_db_sequences(db_sequences, self.cur, self.con)
            self.con.commit()
        except Exception as e:
            raise e

        self.logger.info('write spectra to DB - start - done. Time: '
                         + str(round(time() - db_wrap_up_start_time, 2)) + " sec")

    def upload_info(self):
       self.logger.info('write upload info')
       # peak_list_file_names = json.dumps(self.get_peak_list_file_names(), cls=NumpyEncoder)
       self.upload_id = self.db.write_upload([self.user_id, os.path.basename(self.csv_path), "{}", "{}",
                        "{}", "{}", "{}", "{}", "{}", "{}", "{}", "{}", self.warnings],
                       self.cur, self.con,
                       )
       self.random_id = self.db.get_random_id(self.upload_id, self.cur, self.con)



class xiSPEC_CsvParser(CsvParser):
    required_cols = [
        'scanid',
        'charge',
        'pepseq1',
        'protein1',
        'peaklistfilename',
        # 'expMZ'
    ]

    optional_cols = [
        'rank',
        'fragmenttolerance',
        'iontypes',
        'pepseq2',
        'linkpos1',
        'linkpos2',
        'crosslinkermodmass',
        'passthreshold',
        'score',
        'decoy1',
        'decoy2',
        'protein2',
        'peppos1',
        'peppos2',
        'expmz',    # ToDo: required in mzid - also make required col?
        'calcmz'
    ]

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
        'peppos1': -1,
        'peppos2': -1,
        'expmz': -1,  # ToDo: required in mzid - also make required col?
        'calcmz': -1
    }

    def upload_info(self):
        pass

    def parse_db_sequences(self):
        pass


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
