import re
import ntpath
import json
import sys
import numpy as np
from time import time
import pandas as pd
from xiSPEC_peakList import PeakListReader
import zipfile
import gzip
import os
import psycopg2


class CsvParseException(Exception):
    pass


class MissingFileException(Exception):
    pass


class CsvParser:
    """

    """
    def __init__(self, csv_path, temp_dir, db, logger, db_name=''):
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
        self.db = db
        self.logger = logger

        # self.spectra_data_protocol_map = {}
        # ToDo: Might change to pyteomics unimod obo module
        self.unimod_path = 'obo/unimod.obo'

        self.modlist = []

        self.contains_crosslinks = False

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
        try:
            self.csv_reader = pd.read_csv(self.csv_path)
            self.csv_reader.columns = [x.lower().replace(" ", "") for x in self.csv_reader.columns]

            # ToDo: check required cols
            required_cols = ['id', 'scannumber', 'charge', 'pepseq1', 'protein1', 'runname']

            for required_col in required_cols:
                if required_col not in self.csv_reader.columns:
                    raise CsvParseException("Required csv column %s missing" % required_col)

        except Exception as e:
            raise CsvParseException(type(e).__name__, e.args)

        self.logger.info('reading csv - done. Time: ' + str(round(time() - self.start_time, 2)) + " sec")

    # ToDo: not used atm - can be used for checking if all files are present in temp dir
    def get_peak_list_file_names(self):
        """
        :return: list of all used peak list file names
        """
        return self.csv_reader.runname.unique()

    def get_sequenceDB_file_names(self):
        pass

    def set_peak_list_readers(self):
        """
        sets self.peak_list_readers
        dictionary:
            key: peak list file name
            value: associated peak list reader
        """

        peak_list_readers = {}
        for peak_list_file_name in self.csv_reader.runname.unique():

            if peak_list_file_name.lower().endswith('.mgf'):
                file_format_accession = 'MS:1001062'        # MGF
                spectrum_id_format_accesion = 'MS:1000774'  # MS:1000774 multiple peak list nativeID format - zero based

            elif peak_list_file_name.lower().endswith('.mzml'):
                file_format_accession = 'MS:1000584'        # mzML
                spectrum_id_format_accesion = 'MS:1001530'  # mzML unique identifier
            else:
                raise CsvParseException("Unsupported peak list file type for: %s" % peak_list_file_name)

            peak_list_file_path = self.temp_dir + peak_list_file_name

            try:
                peak_list_reader = PeakListReader(
                    peak_list_file_path,
                    file_format_accession,
                    spectrum_id_format_accesion
                )
            except IOError:
                # try gz version
                try:
                    peak_list_reader = PeakListReader(
                        PeakListReader.extract_gz(peak_list_file_path + '.gz'),
                        file_format_accession,
                        spectrum_id_format_accesion
                    )
                except IOError:
                    # ToDo: output all missing files not just first encountered. Use get_peak_list_file_names()
                    raise CsvParseException('Missing peak list file: %s' % ntpath.basename(peak_list_file_path))

            peak_list_readers[peak_list_file_name] = peak_list_reader

        self.peak_list_readers = peak_list_readers

    def parse(self):

        start_time = time()

        # ToDo: more gracefully handle missing files
        self.set_peak_list_readers()

        #
        # upload info
        #
        # upload_info_start_time = time()
        # self.logger.info('getting upload info (provider, etc) - start')
        # self.parse_upload_info(csv_reader)
        # self.logger.info(
        #     'getting upload info - done. Time: ' + str(round(time() - upload_info_start_time, 2)) + " sec")

        #
        # Sequences, Peptides, Peptide Evidences (inc. peptide positions), Modifications

        # self.parse_db_sequences(self.csv_reader)
        # self.parse_peptides(self.csv_reader)
        # self.parse_peptide_evidences(self.csv_reader)
        # self.map_spectra_data_to_protocol()     # ToDo: does not use self.csv_reader
        self.main_loop()

        #
        # Fill missing scores with
        # score_fill_start_time = time()
        # self.logger.info('fill in missing scores - start')
        # self.db.fill_in_missing_scores(self.cur, self.con)
        # self.logger.info('fill in missing scores - done. Time: ' + str(round(time() - score_fill_start_time, 2)) + " sec")

        self.logger.info('all done! Total time: ' + str(round(time() - start_time, 2)) + " sec")

    # def get_ion_types_mzid(self, sid_item):
    #     try:
    #         ion_names_list = [i['name'] for i in sid_item['IonType']]
    #         ion_names_list = list(set(ion_names_list))
    #     except KeyError:
    #         return []
    #
    #     # ion_types = ["P"]
    #     ion_types = []
    #     for ion_name in ion_names_list:
    #         try:
    #             ion = re.search('frag: ([a-z]) ion', ion_name).groups()[0]
    #             ion_types.append(ion)
    #         except (IndexError, AttributeError) as e:
    #             self.logger.info(e, ion_name)
    #             continue
    #
    #     return ion_types

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
                raise StandardError("more than one mzid file found!")

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

    # # @staticmethod
    # def map_spectra_data_to_protocol(self):
    #     """
    #     extract and map spectrumIdentificationProtocol which includes annotation data like fragment tolerance
    #     only fragment tolerance is extracted for now
    #     # ToDo: improve error handling
    #     #       extract modifications, cl mod mass, ...
    #
    #     Parameters:
    #     ------------------------
    #     csv_reader: pyteomics csv_reader
    #     """
    #
    #     self.logger.info('generating spectra data protocol map - start')
    #     start_time = time()
    #
    #     spectra_data_protocol_map = {}
    #
    #     sid_protocols = []
    #
    #     analysis_collection = self.csv_reader.iterfind('AnalysisCollection').next()
    #     for spectrumIdentification in analysis_collection['SpectrumIdentification']:
    #         sid_protocol_ref = spectrumIdentification['spectrumIdentificationProtocol_ref']
    #         sid_protocol = self.csv_reader.get_by_id(sid_protocol_ref, tag_id='SpectrumIdentificationProtocol', detailed=True)
    #         sid_protocols.append(sid_protocol)
    #         try:
    #             frag_tol = sid_protocol['FragmentTolerance']
    #             frag_tol_plus = frag_tol['search tolerance plus value']
    #             frag_tol_value = re.sub('[^0-9,.]', '', str(frag_tol_plus['value']))
    #             if frag_tol_plus['unit'].lower() == 'parts per million':
    #                 frag_tol_unit = 'ppm'
    #             elif frag_tol_plus['unit'].lower() == 'dalton':
    #                 frag_tol_unit = 'Da'
    #             else:
    #                 frag_tol_unit = frag_tol_plus['unit']
    #
    #             if not all([
    #                 frag_tol['search tolerance plus value']['value'] == frag_tol['search tolerance minus value']['value'],
    #                 frag_tol['search tolerance plus value']['unit'] == frag_tol['search tolerance minus value']['unit']
    #             ]):
    #                 self.warnings.append(
    #                     {"type": "mzidParseError",
    #                      "message": "search tolerance plus value doesn't match minus value. Using plus value!"})
    #
    #         except KeyError:
    #             self.warnings.append({
    #                 "type": "mzidParseError",
    #                 "message": "could not parse ms2tolerance. Falling back to default values.",
    #                 # 'id': id_string
    #             })
    #             frag_tol_value = '10'
    #             frag_tol_unit = 'ppm'
    #             # spectra_data_protocol_map['errors'].append(
    #             #     {"type": "mzidParseError",
    #             #      "message": "could not parse ms2tolerance. Falling back to default values."})
    #
    #         for inputSpectra in spectrumIdentification['InputSpectra']:
    #             spectra_data_ref = inputSpectra['spectraData_ref']
    #
    #             spectra_data_protocol_map[spectra_data_ref] = {
    #                 'protocol_ref': sid_protocol_ref,
    #                 'fragmentTolerance': ' '.join([frag_tol_value, frag_tol_unit])
    #             }
    #
    #     self.csv_reader.reset()
    #     self.spectra_data_protocol_map = spectra_data_protocol_map
    #     self.logger.info('generating spectraData_ProtocolMap - done. Time: ' + str(round(time() - start_time, 2)) + " sec")
    #     # self.db.write.protocols()

    # def add_to_modlist(self, mod):
    #     # modlist_test = [
    #     #     {'name': "bs3nh2", 'monoisotopicMassDelta': 123, 'residues': ["A"]}
    #     # ]
    #     # mod1_test = {'name': "unknown_modification", 'monoisotopicMassDelta': 123, 'residues': ["B"]}
    #     # mod2_test = {'name': "bs3nh2", 'monoisotopicMassDelta': 1232, 'residues': ["B"]}
    #     # mod3_test = {'name': "bs3nh2", 'monoisotopicMassDelta': 12323, 'residues': ["A"]}
    #     # print add_to_modlist(mod1_test, modlist_test)
    #     # print modlist_test
    #     # print add_to_modlist(mod2_test, modlist_test)
    #     # print modlist_test
    #     # print add_to_modlist(mod3_test, modlist_test)
    #     # print modlist_test
    #
    #     if mod['name'] == "unknown_modification":
    #         mod['name'] = "({0:.2f})".format(mod['monoisotopicMassDelta'])
    #
    #     mod['monoisotopicMassDelta'] = float(mod['monoisotopicMassDelta'])
    #
    #     mod['residues'] = [aa for aa in mod['residues']]
    #
    #     if mod['name'] in [m['name'] for m in self.modlist]:
    #         old_mod = self.modlist[[m['name'] for m in self.modlist].index(mod['name'])]
    #         # check if modname with different mass exists already
    #         if mod['monoisotopicMassDelta'] != old_mod['monoisotopicMassDelta']:
    #             mod['name'] += "*"
    #             self.add_to_modlist(mod)
    #         else:
    #             for res in mod['residues']:
    #                 if res not in old_mod['residues']:
    #                     old_mod['residues'].append(res)
    #     else:
    #         self.modlist.append(mod)
    #
    #     return mod['name']

    # def parse_db_sequences(self, csv_reader):

        # self.logger.info('parse db sequences - start')
        # start_time = time()
        # # DBSEQUENCES
        # inj_list = []
        # # for db_sequence in sequence_collection['DBSequence']:
        # for db_id in csv_reader._offset_index["DBSequence"].keys():
        #     db_sequence = csv_reader.get_by_id(db_id, tag_id='DBSequence', detailed=True)
        #
        #     data = [db_sequence["id"], db_sequence["accession"]]
        #
        #     # name, optional elem att
        #     if "name" in db_sequence :
        #         data.append(db_sequence["name"])
        #     else :
        #         data.append(db_sequence["accession"])
        #
        #     # description, officially not there?
        #     if "protein description"  in db_sequence:
        #         data.append(json.dumps(db_sequence["protein description"], cls=NumpyEncoder))
        #     else:
        #         data.append(None)
        #
        #     #searchDatabase_ref
        #
        #     # sequence
        #     if "Seq" in db_sequence and isinstance(db_sequence["Seq"], basestring):  # Seq is optional child elem of DBSequence
        #         seq = db_sequence["Seq"]
        #         data.append(seq)
        #     else:
        #         # todo: get sequence
        #         data.append("no sequence")
        #
        #     # is_decoy - not there
        #     #data.append("false")
        #
        #     data.append(self.upload_id)
        #
        #     inj_list.append(data)
        #
        # self.db.write_db_sequences(inj_list, self.cur, self.con)
        #
        # self.logger.info('parse db sequences - done. Time: ' + str(round(time() - start_time, 2)) + " sec")

    # def parse_peptides(self, csv_reader):
    #     # start_time = time()
    #     # self.logger.info('parse peptides, modifications - start')
    #     #
    #     # self.peptide_id_lookup = {}
    #     #
    #     # # ToDo: might be stuff in pyteomics lib for this?
    #     # unimod_masses = self.get_unimod_masses(self.unimod_path)
    #     # mod_aliases = {
    #     #     # "amidated_bs3": "bs3nh2",
    #     #     # "carbamidomethyl": "cm",
    #     #     # "hydrolyzed_bs3": "bs3oh",
    #     #     # "oxidation": "ox"
    #     # }
    #
    #     #
    #     # # PEPTIDES
    #     # peptide_index = 0
    #     # peptide_inj_list = []
    #     # for pep_id in csv_reader._offset_index["Peptide"].keys():
    #     #     peptide = csv_reader.get_by_id(pep_id, tag_id='Peptide', detailed=True)
    #     #     # peptide2 = csv_reader.get_by_id(pep_id, tag_id='Peptide', detailed=True, accession_key=True)
    #     #     pep_seq_dict = []
    #     #     for aa in peptide['PeptideSequence']:
    #     #         pep_seq_dict.append({"Modification": "", "aminoAcid": aa})
    #     #
    #     #     link_site = -1
    #     #     crosslinker_modmass = -1
    #     #
    #     #     value = -1
    #     #
    #     #     # MODIFICATIONS
    #     #     # add in modifications
    #     #     if 'Modification' in peptide.keys():
    #     #         for mod in peptide['Modification']:
    #     #
    #     #             if 'monoisotopicMassDelta' not in mod.keys():
    #     #                 try:
    #     #                     mod['monoisotopicMassDelta'] = unimod_masses[mod['accession']]
    #     #
    #     #                 except KeyError:
    #     #                     #seq_ref_prot_map['errors'].append({
    #     #                     #     "type": "mzidParseError",
    #     #                     #     "message": "could not get modification mass for modification {}" % mod,
    #     #                     #     "id": mod["id"]
    #     #                     # })
    #     #                     continue
    #     #
    #     #             # link_index = 0  # TODO: multilink support
    #     #
    #     #             if mod['location'] == 0:
    #     #                 mod_location = 0
    #     #                 n_terminal_mod = True
    #     #             elif mod['location'] == len(peptide['PeptideSequence']) + 1:
    #     #                 mod_location = mod['location'] - 2
    #     #                 c_terminal_mod = True
    #     #             else:
    #     #                 mod_location = mod['location'] - 1
    #     #                 n_terminal_mod = False
    #     #                 c_terminal_mod = False
    #     #             if 'residues' not in mod:
    #     #                 mod['residues'] = peptide['PeptideSequence'][mod_location]
    #     #
    #     #             if 'name' in mod.keys():
    #     #                 # fix mod names
    #     #                 if isinstance(mod['name'], list):  # todo: have a look at this  - cc
    #     #                     mod['name'] = ','.join(mod['name'])
    #     #                 mod['name'] = mod['name'].lower()
    #     #                 mod['name'] = mod['name'].replace(" ", "_")
    #     #                 if mod['name'] in mod_aliases.keys():
    #     #                     mod['name'] = mod_aliases[mod['name']]
    #     #                 if 'cross-link donor' not in mod.keys() and 'cross-link acceptor' not in mod.keys():
    #     #                     cur_mod = pep_seq_dict[mod_location]
    #     #                     # join modifications into one for multiple modifications on the same aa
    #     #                     if not cur_mod['Modification'] == '':
    #     #                         mod['name'] = '_'.join(sorted([cur_mod['Modification'], mod['name']], key=str.lower))
    #     #                         cur_mod_mass = [x['monoisotopicMassDelta'] for x in self.modlist if x['name'] == cur_mod['Modification']][0]
    #     #                         mod['monoisotopicMassDelta'] += cur_mod_mass
    #     #
    #     #                     mod['name'] = self.add_to_modlist(mod)  # save to all mods list and get back new_name
    #     #                     cur_mod['Modification'] = mod['name']
    #     #
    #     #             # error handling for mod without name
    #     #             else:
    #     #                 # cross-link acceptor doesn't have a name
    #     #                 if 'cross-link acceptor' not in mod.keys():
    #     #                     pass
    #     #                     # logger.error('modification without name!')
    #     #                     # logger.error(mod)
    #     #
    #     #             # add CL locations
    #     #             if 'cross-link donor' in mod.keys() or 'cross-link acceptor' in mod.keys() or 'cross-link receiver' in mod.keys():
    #     #                 link_site = mod_location
    #     #                 # return_dict['linkSites'].append(
    #     #                 #     {"id": link_index, "peptideId": pep_index, "linkSite": mod_location - 1})
    #     #             if 'cross-link acceptor' in mod.keys():
    #     #                 crosslinker_modmass = mod['monoisotopicMassDelta']
    #     #                 value = mod['cross-link acceptor']['value']
    #     #             if 'cross-link donor' in mod.keys():
    #     #                 crosslinker_modmass = mod['monoisotopicMassDelta']
    #     #                 value = mod['cross-link donor']['value']
    #     #
    #     #     # we should consider swapping these over because modX format has modification before AA
    #     #     peptide_seq_with_mods = ''.join([''.join([x['aminoAcid'], x['Modification']]) for x in pep_seq_dict])
    #     #
    #     #     # data.append(peptide["PeptideSequence"])  # PeptideSequence, required child elem
    #     #     data = [peptide_index, peptide_seq_with_mods, link_site, crosslinker_modmass, self.upload_id, value]
    #     #
    #     #     peptide_inj_list.append(data)
    #     #     self.peptide_id_lookup[peptide['id']] = peptide_index
    #     #     peptide_index += 1
    #
    #     self.db.write_peptides(peptide_inj_list, self.cur, self.con)
    #
    #     #
    #     mod_index = 0
    #     modifications_inj_list = []
    #     for mod in self.modlist:
    #         try:
    #             mod_accession = mod['accession']
    #         except KeyError:
    #             mod_accession = ''
    #         modifications_inj_list.append([
    #             mod_index,
    #             self.upload_id,
    #             mod['name'],
    #             mod['monoisotopicMassDelta'],
    #             ''.join(mod['residues']),
    #             mod_accession
    #         ])
    #         mod_index += 1
    #     self.db.write_modifications(modifications_inj_list, self.cur, self.con)
    #
    #     self.logger.info('parse peptides, modifications - done. Time: ' + str(round(time() - start_time, 2)) + " sec")

    # def parse_peptide_evidences(self, csv_reader):
    #
    #     db_seq_ref_prot_map = {}
    #
    #     sequence_collection = csv_reader.iterfind('SequenceCollection').next()
    #     for sequence in sequence_collection['DBSequence']:
    #         db_seq_ref_prot_map[sequence["id"]] = sequence["accession"]
    #     csv_reader.reset()
    #
    #     start_time = time()
    #     self.logger.info('parse peptide evidences - start')
    #     #PEPTIDE EVIDENCES
    #     inj_list = []
    #     # for peptide_evidence in sequence_collection['PeptideEvidence']:
    #     for pep_ev_id in csv_reader._offset_index["PeptideEvidence"].keys():
    #         peptide_evidence = csv_reader.get_by_id(pep_ev_id, tag_id='PeptideEvidence', detailed=True)
    #
    #         data = []  # peptide_ref, dBSequence_ref, protein_accession, start, upload_id
    #         data.append(self.peptide_id_lookup[peptide_evidence["peptide_ref"]])  # peptide_ref att, required
    #         data.append(peptide_evidence["dBSequence_ref"])  # DBSequence_ref att, required
    #         data.append(db_seq_ref_prot_map[peptide_evidence["dBSequence_ref"]])
    #
    #         if "start" in peptide_evidence:
    #             data.append(peptide_evidence["start"])  # start att, optional
    #         else:
    #             data.append(-1)
    #
    #         if "isDecoy" in peptide_evidence:
    #             data.append(peptide_evidence["isDecoy"])  # isDecoy att, optional
    #         else:
    #             data.append(None)
    #
    #         data.append(self.upload_id)
    #
    #         inj_list.append(data)
    #
    #     try:
    #         self.db.write_peptide_evidences(inj_list, self.cur, self.con)
    #     except AttributeError:
    #         pass
    #
    #     self.con.commit()
    #     csv_reader.reset()
    #
    #     self.logger.info('parse peptide evidences - done. Time: ' + str(round(time() - start_time, 2)) + " sec")

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

    def main_loop(self):
        # main loop
        main_loop_start_time = time()
        self.logger.info('main loop - start')

        inj_list_peptide_evidences = []
        inj_list_spectrum_identifications = []
        inj_list_spectra = []
        inj_list_modifications = []
        inj_list_peptides = []

        peptide_id = 0
        cross_linker_pair_id = 0

        for identification_id, id_item in self.csv_reader.iterrows():  # identification_id, id_item = id_df.iterrows().next()

            #
            # PEPTIDE EVIDENCES
            # ToDo: pep_start, dbsequence_ref and isDecoy1-2 split
            pep_evidence1 = [
                self.upload_id,     # upload_id
                id_item.id,         # peptide_ref
                '',                 # dbsequence_ref
                id_item.protein1,   # protein_accession
                '',                 # pep_start
                id_item.isdecoy     # is_decoy
            ]

            inj_list_peptide_evidences.append(pep_evidence1)

            if id_item.pepseq2:
                pep_evidence2 = [
                    self.upload_id,  # upload_id
                    id_item.id,  # peptide_ref
                    '',  # dbsequence_ref
                    id_item.protein2,  # protein_accession
                    '',  # pep_start
                    id_item.isdecoy  # is_decoy
                ]

                inj_list_peptide_evidences.append(pep_evidence2)

            #
            # PEPTIDES
            if id_item.pepseq2:
                my_cross_linker_pair_id = cross_linker_pair_id
                cross_linker_pair_id += 1
            else:
                my_cross_linker_pair_id = -1

            peptide1 = [
                peptide_id,                     # id,
                id_item.pepseq1,                # seq_mods,
                id_item.linkpos1,               # link_site,
                id_item.crosslinkermodmass,     # crosslinker_modmass,
                self.upload_id,                 # upload_id,
                my_cross_linker_pair_id         # crosslinker_pair_id
            ]
            inj_list_peptides.append(peptide1)
            pep1_id = peptide_id    # save cur pep_id to pep1_id for spectrum_identification
            peptide_id += 1

            if id_item.pepseq2:
                peptide2 = [
                    peptide_id,  # id,
                    id_item.pepseq2,  # seq_mods,
                    id_item.linkpos2,  # link_site,
                    id_item.crosslinkermodmass,  # crosslinker_modmass,
                    self.upload_id,  # upload_id,
                    my_cross_linker_pair_id  # crosslinker_pair_id
                ]
                inj_list_peptides.append(peptide2)
                pep2_id = peptide_id  # save cur pep_id to pep1_id for spectrum_identification
                peptide_id += 1
            else:
                pep2_id = None

            #
            # SPECTRUM IDENTIFCATIONS
            # ToDo: experimental_mass_to_charge, calculated_mass_to_charge, spectrum_id
            scores = json.dumps({'score': id_item['score']})

            spectrum_identification = [
                id_item.id,                 # 'id',     #ToDo: check for uniqueness or use separate id?
                self.upload_id,             # 'upload_id',
                spectrum_id,                # 'spectrum_id',
                pep1_id,                    # 'pep1_id',
                pep2_id,                    # 'pep2_id',
                id_item.charge,             # 'charge_state',
                id_item.rank,               # 'rank',
                id_item.passthreshold,      # 'pass_threshold',
                id_item.iontypes,           # 'ions',
                scores,   # 'scores',
                '',                         # 'experimental_mass_to_charge',
                ''                          # 'calculated_mass_to_charge'
            ]
            inj_list_spectrum_identifications.append(spectrum_identification)

            #
            # SPECTRA
            # ToDo: everything
            spectrum = [
                # 'id',
                # 'peak_list',
                # 'peak_list_file_name',
                # 'scan_id',
                # 'frag_tol',
                # 'upload_id',
                # 'spectrum_ref'
            ]
            inj_list_spectra.append(spectrum)

            #
            # MODIFICATIONS
            # ToDo: everything - need global modList again

            # spec_id_set = set()

            # extract scan_id
            # try:
            scan_id = int(id_item['scannumber'])


        #     # peakList
        #     try:
        #         peak_list_reader = peakListParser.get_reader(peak_list_readers, raw_file_name)
        #     except peakListParser.PeakListParseError as e:
        #         return_json['errors'].append({
        #             "type": "peakListParseError",
        #             "message": e.args[0],
        #             'id': id_item['id']
        #         })
        #         continue
        #     try:
        #         scan = peakListParser.get_scan(peak_list_reader, scan_id)
        #     except peakListParser.PeakListParseError:
        #         try:
        #             scan_not_found_error[raw_file_name].append(scan_id)
        #         except KeyError:
        #             scan_not_found_error[raw_file_name] = [scan_id]
        #         continue
        #
        #     peak_list = peakListParser.get_peak_list(scan, peak_list_reader['fileType'])
        #
        #     multiple_inj_list_peak_lists.append([identification_id, peak_list])
        #
        #     # identification id
        #     identification_id = id_item['id']
        #
        #     # rank
        #     try:
        #         rank = id_item['rank']
        #     except KeyError:
        #         rank = 1
        #
        #     # peptides and link positions
        #     pep1 = id_item['pepseq1']
        #     try:
        #         # ToDo: improve error handling for cl peptides
        #         pep2 = str(id_item['pepseq2'])
        #         if pep2 == 'nan':
        #             pep2 = ''
        #         link_pos1 = id_item['linkpos1']
        #         link_pos2 = id_item['linkpos2']
        #         cl_mod_mass = id_item['crosslinkermodmass']
        #     except KeyError:
        #         # linear
        #         pep2 = ''
        #
        #     if pep2 == '':
        #         link_pos1 = -1
        #         link_pos2 = -1
        #         cl_mod_mass = 0
        #
        #     # charge
        #     charge = id_item['charge']
        #
        #     # passThreshold
        #     try:
        #         pass_threshold = 1 if id_item['passthreshold'] else 0
        #     except KeyError:
        #         pass_threshold = 1
        #
        #     # fragment tolerance
        #     frag_tol = id_item['fragmenttolerance'].split(' ', 1)
        #     if frag_tol[1].lower() == 'parts per million':
        #         frag_tol[1] = 'ppm'
        #     elif frag_tol[1].lower() == 'dalton':
        #         frag_tol[1] = 'Da'
        #     if frag_tol[1] not in ['ppm', 'Da']:
        #         return_json['errors'].append({
        #             "type": "fragTolParseError",
        #             "message": "unknown fragment tolerance unit: %s\nSupported values are: ppm, Da" % frag_tol[1],
        #             'id': id_item['id']
        #         })
        #         continue
        #     ms2_tol = ' '.join(frag_tol)
        #
        #     # ion types
        #     ion_types = id_item['iontypes'].lower()
        #     unknown_ions = [ion for ion in ion_types.split(';') if ion not in ['peptide', 'a', 'b', 'c', 'x', 'y', 'z']]
        #     if len(unknown_ions) > 0:
        #         return_json['errors'].append({
        #             "type": "ionTypeParseError",
        #             "message": "unknown ion(s) : %s\nSupported values are: peptide, a, b, c, x, y, z" % ';'.join(
        #                 unknown_ions),
        #             'id': id_item['id']
        #         })
        #     # ToDo: could check against mzml fragmentation type and display warning if ions don't match
        #
        #     # score
        #     # score = id_item['score']
        #     all_scores = json.dumps({'score': id_item['score']})
        #
        #     # isDecoy
        #     try:
        #         is_decoy = 1 if id_item['isdecoy'] else 0
        #     except KeyError:
        #         is_decoy = 0
        #
        #     # protein
        #     protein1 = id_item['protein1']
        #     try:
        #         protein2 = id_item['protein2']
        #     except KeyError:
        #         protein2 = ''
        #
        #     # create entry
        #     multiple_inj_list_identifications.append(
        #         [identification_id,
        #          identification_id,
        #          pep1,
        #          pep2,
        #          link_pos1,
        #          link_pos2,
        #          charge,
        #          pass_threshold,
        #          ms2_tol,
        #          ion_types,
        #          cl_mod_mass,
        #          rank,
        #          # score,
        #          all_scores,
        #          is_decoy,
        #          protein1,
        #          protein2,
        #          raw_file_name,
        #          scan_id,
        #          identification_id]
        #     )
        #
        #     try:
        #         modifications = re.search('([^A-Z]+)', ''.join([pep1, pep2])).groups()
        #     except AttributeError:
        #         modifications = []
        #
        #     for mod in modifications:
        #         if mod not in return_json['modifications']:
        #             return_json['modifications'].append(mod)
        #
        #     #  write to DB
        #     if identification_id % 1000 == 0:
        #         logger.info('writing 1000 entries to DB')
        #         try:
        #             db.write_identifications(multiple_inj_list_identifications, cur, con)
        #             multiple_inj_list_identifications = []
        #
        #             db.write_peaklists(multiple_inj_list_peak_lists, cur, con)
        #             multiple_inj_list_peak_lists = []
        #
        #         except db.DBException as e:
        #             return_json['errors'].append(
        #                 {"type": "dbError",
        #                  "message": e.args[0],
        #                  'id': identification_id
        #                  })
        #             return return_json
        #
        #         # commit changes
        #         con.commit()
        #
        # # end main loop
        # self.logger.info('main loop - done. Time: ' + str(round(time() - main_loop_start_time, 2)) + " sec")
        #
        # # once loop is done write remaining data to DB
        # db_wrap_up_start_time = time()
        # self.logger.info('write spectra to DB - start')
        # try:
        #     self.db.write_spectra(spectra, self.cur, self.con)
        #     self.db.write_spectrum_identifications(spectrum_identifications, self.cur, self.con)
        #     self.con.commit()
        # except Exception as e:
        #     raise e
        #
        # self.logger.info('write spectra to DB - start - done. Time: '
        #             + str(round(time() - db_wrap_up_start_time, 2)) + " sec")
        #
        # # warnings
        # if len(fragment_parsing_error_scans) > 0:
        #     if len(fragment_parsing_error_scans) > 50:
        #         id_string = '; '.join(fragment_parsing_error_scans[:50]) + ' ...'
        #     else:
        #         id_string = '; '.join(fragment_parsing_error_scans)
        #         self.warnings.append({
        #             "type": "IonParsing",
        #             "message": "mzidentML file does not specify fragment ions.",
        #             'id': id_string
        #         })

    # def upload_info(self, csv_reader):
    #     # AnalysisSoftwareList - optional element
    #     # see https://groups.google.com/forum/#!topic/pyteomics/Mw4eUHmicyU
    #     csv_reader.schema_info['lists'].add("AnalysisSoftware")
    #     try:
    #         analysis_software = json.dumps(csv_reader.iterfind('AnalysisSoftwareList').next()['AnalysisSoftware'])
    #     except StopIteration:
    #         analysis_software = '{}'
    #     csv_reader.reset()
    #
    #     # Provider - optional element
    #     provider = '{}'
    #     # try:
    #     #     provider = json.dumps(csv_reader.iterfind('Provider').next())
    #     # except StopIteration:
    #     #     pass
    #     # csv_reader.reset()
    #
    #     # AuditCollection - optional element
    #     audits = '{}'
    #     # try:
    #     #     audits = json.dumps(csv_reader.iterfind('AuditCollection').next())
    #     # except StopIteration:
    #     #     audits = '{}'
    #     # csv_reader.reset()
    #
    #     # AnalysisSampleCollection - optional element
    #     samples = '{}'
    #     # try:
    #     #     samples = json.dumps(csv_reader.iterfind('AnalysisSampleCollection').next()['Sample'])
    #     # except StopIteration:
    #     #     samples = '{}'
    #     # csv_reader.reset()
    #
    #     # AnalysisCollection - required element
    #     analyses = '{}'
    #     # analyses = json.dumps(csv_reader.iterfind('AnalysisCollection').next()['SpectrumIdentification'])
    #     # csv_reader.reset()
    #
    #     # AnalysisProtocolCollection - required element
    #     protocols ='{}'
    #     # protocol_collection = csv_reader.iterfind('AnalysisProtocolCollection').next()
    #     # protocols = protocol_collection['SpectrumIdentificationProtocol']
    #     protocols = json.dumps(protocols, cls=NumpyEncoder)
    #     # csv_reader.reset()
    #
    #     # BibliographicReference - optional element
    #     bibRefs = []
    #     # for bib in csv_reader.iterfind('BibliographicReference'):
    #     #     bibRefs.append(bib)
    #     # bibRefs = json.dumps(bibRefs)
    #     # csv_reader.reset()
    #
    #
    #     # # AnalysisSoftwareList - optional element
    #     # # see https://groups.google.com/forum/#!topic/pyteomics/Mw4eUHmicyU
    #     # csv_reader.schema_info['lists'].add("AnalysisSoftware")
    #     # analysis_software = {}
    #     # for analysis_software_list_id in csv_reader._offset_index["AnalysisSoftwareList"].keys():
    #     #     analysis_software = csv_reader.get_by_id(analysis_software_list_id, tag_id='AnalysisSoftwareList',
    #     #                                               detailed=True)['AnalysisSoftware']
    #     # analysis_software = json.dumps(analysis_software)
    #     #
    #     # # try:
    #     # #     analysis_software = json.dumps(csv_reader.iterfind('AnalysisSoftwareList').next()['AnalysisSoftware'])
    #     # # except StopIteration:
    #     # #     analysis_software = '{}'
    #     # # csv_reader.reset()
    #     #
    #     # # Provider - optional element
    #     # try:
    #     #     provider = json.dumps(csv_reader.iterfind('Provider').next())
    #     # except StopIteration:
    #     #     provider = '{}'
    #     # csv_reader.reset()
    #     # analysis_software = {}
    #     # for analysis_software_list_id in csv_reader._offset_index["AnalysisSoftwareList"].keys():
    #     #     analysis_software = csv_reader.get_by_id(analysis_software_list_id, tag_id='AnalysisSoftwareList',
    #     #                                               detailed=True)['AnalysisSoftware']
    #     # analysis_software = json.dumps(analysis_software)
    #     #
    #     # # AuditCollection - optional element
    #     # # try:
    #     # #     audits = json.dumps(csv_reader.iterfind('AuditCollection').next())
    #     # # except StopIteration:
    #     # #     audits = '{}'
    #     # # csv_reader.reset()
    #     #
    #     # audits = {}
    #     # for audit_collection_id in csv_reader._offset_index["AuditCollection"].keys():
    #     #     audits = csv_reader.get_by_id(audit_collection_id, tag_id='AuditCollection',
    #     #                                               detailed=True)
    #     # audits = json.dumps(audits)
    #     #
    #     # upload_id = self.db.write_upload([self.user_id, filename, json.dumps(peak_list_file_names),
    #     #                   analysis_software, provider, audits, samples, analyses, protocols, bibRefs, 'country'],
    #     #                  self.cur, self.con,
    #     #                  )
    #
    # def fetch_supporting_files_from_PRIDE(self, origin, ip, base):
    #     """
    #     e.g.
    #     ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/09/PXD001054/
    #
    #     origin: 2014/09/PXD001054/
    #     ip: ftp.pride.ebi.ac.uk
    #     base: pride/data/archive/
    #
    #     :param origin: year month project folder on ftp
    #     :param ip: ip of the ftp server
    #     :param base: ftp base path
    #     """
    #     #  ignoring source files
    #     # todo - self.search_databases
    #
    #     # get spectra data
    #     spectra_data = {}
    #     # for spectra_data_id in self.csv_reader._offset_index["SpectraData"].keys():
    #     #     sp_datum = self.csv_reader.get_by_id(spectra_data_id, tag_id='SpectraData', detailed=True)
    #     #     sd_id = sp_datum['id']
    #     #     # peak list file name
    #     #     # if 'location' in spectra_data.keys():
    #     #     #  location is required attribute
    #     #     peak_list_file_name = ntpath.basename(sp_datum['location'])
    #     #     self.peak_list_file_names.append(peak_list_file_name)
    #     #     # elif 'name' in spectra_data.keys():
    #     #     #     peak_list_file_name = spectra_data['name']
    #     #
    #     #     file_format_readable = ''
    #     #     if 'FileFormat' in sp_datum:
    #     #         if 'accession' in sp_datum['FileFormat']:
    #     #             file_format_readable = sp_datum['FileFormat']['accession'] + ': ' + sp_datum['FileFormat']['name']
    #     #             # missing .mgf file extension?
    #     #             if sp_datum['FileFormat']['accession'] == 'MS:1001062' and not peak_list_file_name.lower().endswith(
    #     #                         '.mgf'):
    #     #                 self.warnings.append('location of mgf file is missing .mgf extension: ' + peak_list_file_name)
    #     #                 peak_list_file_name = peak_list_file_name + '.mgf'
    #     #             # also might be some fails coz .MGF is being converted to .mgf in ftp archive
    #     #         else:
    #     #             self.warnings.append('SpectraData>FileFormat is missing accession.')
    #     #             file_format_readable = sp_datum['FileFormat']
    #     #
    #     #     else:
    #     #         self.warnings.append('SpectraData missing required element FileFormat.')
    #     #     self.file_formats.add(file_format_readable)
    #     #
    #     #     spec_id_format_readable = ''
    #     #     if 'SpectrumIDFormat' in sp_datum:
    #     #         if 'accession' in sp_datum['SpectrumIDFormat']:
    #     #             spec_id_format_readable = sp_datum['SpectrumIDFormat']['accession'] + ': ' + \
    #     #                                       sp_datum['SpectrumIDFormat']['name']
    #     #         else:
    #     #             self.warnings.append('SpectraData>SpectrumIDFormat is missing accession.')
    #     #             spec_id_format_readable = sp_datum['SpectrumIDFormat']
    #     #     else:
    #     #         self.warnings.append('SpectraData missing required element SpectrumIDElement.')
    #     #     self.spectrum_id_formats.add(spec_id_format_readable)
    #     #
    #     #     ftp = ftplib.FTP(ip)
    #     #     ftp.login()  # Uses password: anonymous@
    #     #     try:
    #     #         target_dir = '/' + base + '/' + origin
    #     #         ftp.cwd(target_dir)
    #     #         self.logger.info('getting ' + peak_list_file_name)
    #     #         ftp.retrbinary("RETR " + peak_list_file_name,
    #     #                        open(self.temp_dir + '/' + peak_list_file_name, 'wb').write)
    #     #     except ftplib.error_perm as e:
    #     #         #  check for gzipped
    #     #         try:
    #     #             self.logger.info('getting ' + peak_list_file_name + '.gz')
    #     #             ftp.retrbinary("RETR " + peak_list_file_name + '.gz',
    #     #                            open(self.temp_dir + '/' + peak_list_file_name + '.gz', 'wb').write)
    #     #         except ftplib.error_perm as e:
    #     #             ftp.close()
    #     #             raise MissingFileException(type(e).__name__, peak_list_file_name, e.args)
    #     #         ftp.close()
    #     #         peak_list_file_name = ntpath.basename(
    #     #             peakListParser.unzip_peak_lists(self.temp_dir + '/' + peak_list_file_name + '.gz')[0])
    #     #     ftp.close()
    #     #
    #     #     sp_datum['peak_list_reader'] = peakListParser.get_peak_list_reader(self.temp_dir + '/' + peak_list_file_name)
    #     #     spectra_data[sd_id] = sp_datum
    #
    #     return spectra_data
    #
    #
    # def log_preliminary_info(self):
    #
    #     # info used for debugging info in db
    #     spectrum_id_formats = set([])
    #     file_formats = set([])
    #     peak_list_file_names = []
    #
    #     # get some preliminary info before actual parsing starts for error logging to db
    #     mzId_file_size = os.path.getsize(self.csv_path)
    #     mzId_stream = open(self.csv_path, 'r')
    #
    #     file_start = mzId_stream.read(10000)  # read first 10000 chars in. Should include all necessary info
    #     mzId_stream.close()
    #     version_match = re.search('mzIdentML.*?version="(.*?)"', file_start,  flags=re.IGNORECASE)
    #     if version_match is not None:
    #         xml_version = version_match.group(1)
    #     else:
    #         xml_version = ''
    #         self.warnings.append("Missing mzid version info.")
    #
    #     software_iter = re.finditer('<SoftwareName>.*?<cvParam.*?name="(.*?)".*?</SoftwareName>', file_start, re.DOTALL)
    #     analysis_software = []
    #     for i in software_iter:
    #         analysis_software.append(i.group(1))
    #     analysis_software = json.dumps(analysis_software, cls=NumpyEncoder)
    #
    #     try:
    #         self.cur.execute("""
    #             INSERT INTO uploads (
    #                 origin,
    #                 filename,
    #                 xml_version,
    #                 file_size,
    #                 analysis_software
    #             )
    #             VALUES (%s, %s, %s, %s, %s) RETURNING id AS upload_id""",
    #                          [self.origin, ntpath.basename(self.csv_path), xml_version, mzId_file_size, analysis_software])
    #         self.con.commit()
    #
    #     except psycopg2.Error as e:
    #         raise self.db.DBException(e.message)
    #
    #     rows = self.cur.fetchall()
    #     self.upload_id = rows[0][0]
    #
    # def log_post_parse_info(self):
    #     pass
    #     # parse_time = str(round(time() - mzid_start_time, 2))
    #     # store some info on how things went
    #     # try:
    #     #     self.cur.execute("""
    #     #         UPDATE uploads SET
    #     #             peak_list_file_names=%s,
    #     #             upload_warnings=%s,
    #     #             spectrum_id_format=%s,
    #     #             file_format=%s,
    #     #             parse_time=%s,
    #     #             contains_crosslinks=%s
    #     #          WHERE id = %s""",
    #     #                      [
    #     #                          json.dumps(self.peak_list_file_names, cls=NumpyEncoder),
    #     #                          json.dumps(self.warnings),
    #     #                          ','.join(self.spectrum_id_formats),
    #     #                          ','.join(self.file_formats),
    #     #
    #     #                          parse_time,
    #     #                          self.contains_crosslinks,
    #     #                          self.upload_id
    #     #                      ])
    #     #     self.con.commit()
    #     # except psycopg2.Error as e:
    #     #     raise e


class xiSPEC_CsvParser(CsvParser):
    def parse_db_sequences(self, csv_reader):
        pass


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
