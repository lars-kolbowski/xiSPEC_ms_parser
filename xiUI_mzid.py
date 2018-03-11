import pyteomics.mzid as py_mzid
import re
import ntpath
import json
import sys
import numpy as np
from time import time
import xiSPEC_peakList as peakListParser
import zipfile
import gzip
import os
import psycopg2
import ftplib

class MzIdParseException(Exception):
    pass

class MzIdMissingFileException(Exception):
    pass

class MzIdScanNotFoundException(Exception):
    pass

class MzIdParser:

    def __init__(self, mzId_path, temp_dir, origin, db, ip, base, logger):

        self.temp_dir = temp_dir
        self.origin = origin
        self.db = db
        self.ip = ip
        self.base = base
        self.logger = logger

        self.unimod_path = 'obo/unimod.obo'
        self.modlist = []


        self.peak_list_readers = {}  # peakListParser.create_peak_list_readers(peak_list_file_list)
        self.spectrum_id_formats = set([])
        self.contains_crosslinks = False

        self.warnings = []

        # connect to DB
        try:
            self.con = db.connect('')
            self.cur = self.con.cursor()

        except db.DBException as e:
            self.logger.error(e)
            print(e)
            sys.exit(1)

        if mzId_path.endswith('.gz') or mzId_path.endswith('.zip'):
            mzId_path = MzIdParser.extract_mzid(mzId_path)

        self.mzId_path = mzId_path
        # get some preliminary info
        self.mzId_file_size = os.path.getsize(mzId_path)
        mzId_stream = open(mzId_path, 'r')

        file_start = mzId_stream.read(5000)  # read by character
        mzId_stream.close()
        self.xml_version = re.search('mzIdentML.*?version="(.*?)"', file_start,  flags=re.IGNORECASE).group(1)

        software_iter = re.finditer('<SoftwareName>.*?<cvParam.*?name="(.*?)".*?</SoftwareName>', file_start, re.DOTALL)
        analysis_software = []
        for i in software_iter:
            analysis_software.append(i.group(1))
        analysis_software = json.dumps(analysis_software, cls=NumpyEncoder)

        try:
            self.cur.execute("""
                INSERT INTO uploads (
                    filename,
                    origin,
                    xml_version,
                    file_size,
                    analysis_software
                )
                VALUES (%s, %s, %s, %s, %s) RETURNING id AS upload_id""",
                             [mzId_path, origin, self.xml_version, self.mzId_file_size, analysis_software])
            self.con.commit()

        except psycopg2.Error as e:
            raise db.DBException(e.message)

        rows = self.cur.fetchall()
        self.upload_id = rows[0][0]

    def get_ion_types_mzid(self, sid_item):
        try:
            ion_names_list = [i['name'] for i in sid_item['IonType']]
            ion_names_list = list(set(ion_names_list))
        except KeyError:
            return []

        # ion_types = ["P"]
        ion_types = []
        for ion_name in ion_names_list:
            try:
                ion = re.search('frag: ([a-z]) ion', ion_name).groups()[0]
                ion_types.append(ion)
            except (IndexError, AttributeError) as e:
                self.logger.info(e, ion_name)
                continue

        return ion_types

    # split into two functions
    @staticmethod
    def extract_mzid(archive):
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

    # ToDo: clear confusion about 0 & 1 based formats
    @staticmethod
    def get_scan(spec_id, spec_id_format):
        # mzml 0 based

        #
        # if (fileIdFormat == Constants.SpecIdFormat.MASCOT_QUERY_NUM) {
        #     String rValueStr = spectrumID.replaceAll("query=", "");
        #     String id = null;
        #     if(rValueStr.matches(Constants.INTEGER)){
        #         id = Integer.toString(Integer.parseInt(rValueStr) + 1);
        #     }
        #     return id;
        # } else if (fileIdFormat == Constants.SpecIdFormat.MULTI_PEAK_LIST_NATIVE_ID) {
        #     String rValueStr = spectrumID.replaceAll("index=", "");
        #     String id;
        #     if(rValueStr.matches(Constants.INTEGER)){
        #         id = Integer.toString(Integer.parseInt(rValueStr) + 1);
        #         return id;
        #     }
        #     return spectrumID;
        # } else if (fileIdFormat == Constants.SpecIdFormat.SINGLE_PEAK_LIST_NATIVE_ID) {
        #     return spectrumID.replaceAll("file=", "");
        # } else if (fileIdFormat == Constants.SpecIdFormat.MZML_ID) {
        #     return spectrumID.replaceAll("mzMLid=", "");
        # } else if (fileIdFormat == Constants.SpecIdFormat.SCAN_NUMBER_NATIVE_ID) {
        #     return spectrumID.replaceAll("scan=", "");
        # } else {
        #     return spectrumID;
        # }


        # e.g.: MS:1000768(Thermo        nativeID        format)
        # e.g.: MS:1000769(Waters        nativeID        format)
        # e.g.: MS:1000770(WIFF        nativeID        format)
        # e.g.: MS:1000771(Bruker / Agilent        YEP        nativeID        format)
        # e.g.: MS:1000772(Bruker        BAF        nativeID        format)
        # e.g.: MS:1000773(Bruker        FID        nativeID        format)
        # e.g.: MS:1000774(multiple       peak        list        nativeID        format)
        # e.g.: MS:1000775(single        peak        list        nativeID        format)
        # e.g.: MS:1000776(scan        number        only        nativeID        format)
        # e.g.: MS:1000777(spectrum        identifier        nativeID        format)

        if spec_id_format['accession'] == 'MS:1000774':  # (multiple peak list nativeID format - zero based)

            matches = re.findall("(?:query|index)=([0-9]+)", spec_id)
            if len(matches) == 1:
                return int(matches[0])

    @staticmethod
    def map_spectra_data_to_protocol(mzid_reader):
        """
        extract and map spectrumIdentificationProtocol which includes annotation data like fragment tolerance
        only fragment tolerance is extracted for now
        # ToDo: improve error handling
        #       extract modifications, cl mod mass, ...

        Parameters:
        ------------------------
        mzid_reader: pyteomics mzid_reader
        """

        spectra_data_protocol_map = {
            'errors': [],
        }

        analysis_collection = mzid_reader.iterfind('AnalysisCollection').next()
        for spectrumIdentification in analysis_collection['SpectrumIdentification']:
            sid_protocol_ref = spectrumIdentification['spectrumIdentificationProtocol_ref']
            sid_protocol = mzid_reader.get_by_id(sid_protocol_ref, tag_id='SpectrumIdentificationProtocol', detailed=True)

            try:
                frag_tol = sid_protocol['FragmentTolerance']
                frag_tol_plus = frag_tol['search tolerance plus value']
                frag_tol_value = re.sub('[^0-9,.]', '', str(frag_tol_plus['value']))
                if frag_tol_plus['unit'].lower() == 'parts per million':
                    frag_tol_unit = 'ppm'
                elif frag_tol_plus['unit'].lower() == 'dalton':
                    frag_tol_unit = 'Da'
                else:
                    frag_tol_unit = frag_tol_plus['unit']

                if not all([
                    frag_tol['search tolerance plus value']['value'] == frag_tol['search tolerance minus value']['value'],
                    frag_tol['search tolerance plus value']['unit'] == frag_tol['search tolerance minus value']['unit']
                ]):
                    spectra_data_protocol_map['errors'].append(
                        {"type": "mzidParseError",
                         "message": "search tolerance plus value doesn't match minus value. Using plus value!"})

            except KeyError:
                frag_tol_value = '10'
                frag_tol_unit = 'ppm'
                spectra_data_protocol_map['errors'].append(
                    {"type": "mzidParseError",
                     "message": "could not parse ms2tolerance. Falling back to default values."})

            for inputSpectra in spectrumIdentification['InputSpectra']:
                spectra_data_ref = inputSpectra['spectraData_ref']

                spectra_data_protocol_map[spectra_data_ref] = {
                    'protocol_ref': sid_protocol_ref,
                    'fragmentTolerance': ' '.join([frag_tol_value, frag_tol_unit])
                }

        mzid_reader.reset()

        return spectra_data_protocol_map

    @staticmethod
    def get_cross_link_identifier(sid_item):
        # For reporting the evidence associated with the identification, within a given <SpectrumIdentificationResult>,
        # a pair of cross-linked peptides MUST be reported as two instances of <SpectrumIdentificationItem> through having
        # a shared local unique identifier as the value for the CV term "cross-link spectrum identification item" MS:1002511
        # The two instances of <SpectrumIdentificationItem> MUST also share the same value for the rank attribute.
        #
        # If a cross-linked pair of peptides has been identified, there MUST be two
        # <SpectrumIdentificationItem> elements with the same rank value. Both MUST have the
        # "cross-link spectrum identification item" cvParam, and the value acts as a local identifier
        # within the <SpectrumIdentificationResult> to group these two elements together. The
        # experimentalMassToCharge, calculateMassToCharge and chargeState MUST be identical
        # over both SII elements, indicating the overall values for the pair.

        # cl_id_item = str(int(sid_item['cross-link spectrum identification item']))
        # rank = str(sid_item['rank'])

        return sid_item['cross-link spectrum identification item']
        # return '_'.join([cl_id_item, rank])

    def add_to_modlist(self, mod):
        # modlist_test = [
        #     {'name': "bs3nh2", 'monoisotopicMassDelta': 123, 'residues': ["A"]}
        # ]
        # mod1_test = {'name': "unknown_modification", 'monoisotopicMassDelta': 123, 'residues': ["B"]}
        # mod2_test = {'name': "bs3nh2", 'monoisotopicMassDelta': 1232, 'residues': ["B"]}
        # mod3_test = {'name': "bs3nh2", 'monoisotopicMassDelta': 12323, 'residues': ["A"]}
        # print add_to_modlist(mod1_test, modlist_test)
        # print modlist_test
        # print add_to_modlist(mod2_test, modlist_test)
        # print modlist_test
        # print add_to_modlist(mod3_test, modlist_test)
        # print modlist_test

        if mod['name'] == "unknown_modification":
            mod['name'] = "({0:.2f})".format(mod['monoisotopicMassDelta'])

        mod['monoisotopicMassDelta'] = float(mod['monoisotopicMassDelta'])

        mod['residues'] = [aa for aa in mod['residues']]

        if mod['name'] in [m['name'] for m in self.modlist]:
            old_mod = self.modlist[[m['name'] for m in self.modlist].index(mod['name'])]
            # check if modname with different mass exists already
            if mod['monoisotopicMassDelta'] != old_mod['monoisotopicMassDelta']:
                mod['name'] += "*"
                self.add_to_modlist(mod)
            else:
                for res in mod['residues']:
                    if res not in old_mod['residues']:
                        old_mod['residues'].append(res)
        else:
            self.modlist.append(mod)

        return mod['name']

    def parse_sequence_collection (self, mzid_reader):

        # ToDo: might be stuff in pyteomics lib for this?
        unimod_masses = self.get_unimod_masses(self.unimod_path)
        mod_aliases = {
            # "amidated_bs3": "bs3nh2",
            # "carbamidomethyl": "cm",
            # "hydrolyzed_bs3": "bs3oh",
            # "oxidation": "ox"
        }

        # first get search databases?
        # search_databases = {}
        #
        # try:
        #     for search_database in mzid_reader.iterfind("Inputs").next()["SearchDatabase"]:
        #         print(json.dumps(search_database))
        # except StopIteration:
        #     pass
        # mzid_reader.reset()

        sequence_collection = mzid_reader.iterfind('SequenceCollection').next()

        #DBSEQUENCES
        inj_list = []
        for db_sequence in sequence_collection['DBSequence']:

            data = [db_sequence["id"], db_sequence["accession"]];

            # name, optional elem att
            if "name" in db_sequence :
                data.append(db_sequence["name"])
            else :
                data.append(db_sequence["accession"])

            # description, officially not there?
            if "protein description"  in db_sequence:
                data.append(db_sequence["protein description"])
            else:
                data.append("protein description")

            #searchDatabase_ref

            # sequence
            if "Seq" in db_sequence and isinstance(db_sequence["Seq"], basestring):  # Seq is optional child elem of DBSequence
                seq = db_sequence["Seq"]
                data.append(seq)
            else:
                # todo: get sequence
                data.append("no sequence")

            # is_decoy - not there
            #data.append("false")

            data.append(self.upload_id)

            inj_list.append(data)

        self.db.write_db_sequences(inj_list, self.cur, self.con)

        #PEPTIDES
        inj_list = []
        for pep_id in mzid_reader._offset_index["Peptide"].keys():
            # peptide = mzid_reader.get_by_id('1712126853_1712127183_5_13_p1', tag_id='Peptide', detailed=True)
            peptide = mzid_reader.get_by_id(pep_id, tag_id='Peptide', detailed=True)

            pep_seq_dict = []
            for aa in peptide['PeptideSequence']:
                pep_seq_dict.append({"Modification": "", "aminoAcid": aa})

            link_site = -1
            crosslinker_modmass = -1

            value = -1

            #MODIFICATIONS
            # add in modifications
            if 'Modification' in peptide.keys():
                for mod in peptide['Modification']:

                    if 'monoisotopicMassDelta' not in mod.keys():
                        try:
                            mod['monoisotopicMassDelta'] = unimod_masses[mod['accession']]

                        except KeyError:
                            #seq_ref_prot_map['errors'].append({
                            #     "type": "mzidParseError",
                            #     "message": "could not get modification mass for modification {}" % mod,
                            #     "id": mod["id"]
                            # })
                            continue

                    # link_index = 0  # TODO: multilink support

                    if mod['location'] == 0:
                        mod_location = 0
                        n_terminal_mod = True
                    elif mod['location'] == len(peptide['PeptideSequence']) + 1:
                        mod_location = mod['location'] - 2
                        c_terminal_mod = True
                    else:
                        mod_location = mod['location'] - 1
                        n_terminal_mod = False
                        c_terminal_mod = False
                    if 'residues' not in mod:
                        mod['residues'] = peptide['PeptideSequence'][mod_location]

                    if 'name' in mod.keys():
                        # fix mod names
                        mod['name'] = mod['name'].lower()
                        mod['name'] = mod['name'].replace(" ", "_")
                        if mod['name'] in mod_aliases.keys():
                            mod['name'] = mod_aliases[mod['name']]
                        if 'cross-link donor' not in mod.keys() and 'cross-link acceptor' not in mod.keys():
                            cur_mod = pep_seq_dict[mod_location]
                            # join modifications into one for multiple modifications on the same aa
                            if not cur_mod['Modification'] == '':
                                mod['name'] = '_'.join(sorted([cur_mod['Modification'], mod['name']], key=str.lower))
                                cur_mod_mass = [x['monoisotopicMassDelta'] for x in self.modlist if x['name'] == cur_mod['Modification']][0]
                                mod['monoisotopicMassDelta'] += cur_mod_mass

                            mod['name'] = self.add_to_modlist(mod)  # save to all mods list and get back new_name
                            cur_mod['Modification'] = mod['name']

                    # error handling for mod without name
                    else:
                        # cross-link acceptor doesn't have a name
                        if 'cross-link acceptor' not in mod.keys():
                            pass
                            # logger.error('modification without name!')
                            # logger.error(mod)

                    # add CL locations
                    if 'cross-link donor' in mod.keys() or 'cross-link acceptor' in mod.keys() or 'cross-link receiver' in mod.keys():
                        link_site = mod_location
                        # return_dict['linkSites'].append(
                        #     {"id": link_index, "peptideId": pep_index, "linkSite": mod_location - 1})
                    if 'cross-link acceptor' in mod.keys():
                        crosslinker_modmass = mod['monoisotopicMassDelta']
                        value = mod['cross-link acceptor']['value']
                    if 'cross-link donor' in mod.keys():
                        crosslinker_modmass = mod['monoisotopicMassDelta']
                        value = mod['cross-link donor']['value']

            # we should consider swapping these over because modX format has modification before AA
            peptide_seq_with_mods = ''.join([''.join([x['aminoAcid'], x['Modification']]) for x in pep_seq_dict])

            # data.append(peptide["PeptideSequence"])  # PeptideSequence, required child elem
            data = [peptide["id"], peptide_seq_with_mods, link_site, crosslinker_modmass, self.upload_id, value]

            inj_list.append(data)

        self.db.write_peptides(inj_list, self.cur, self.con)

        # modifications - ToDo
        # modifications = []
        # for mod in all_mods:
        #     try:
        #         mod_accession = mod['accession']
        #     except KeyError:
        #         mod_accession = ''
        #         modifications.append({
        #         'aminoAcids': mod['residues'],
        #         'id': mod['name'],
        #         'mass': mod['monoisotopicMassDelta'],
        #         'accession': mod_accession
        #     })

        # add mods to global modList
        # for mod in pep_info['annotation']['modifications']:
        #     if mod['id'] not in [m['id'] for m in modifications]:
        #         modifications.append(mod)
        #     else:
        #         old_mod = modifications[[m['id'] for m in modifications].index(mod['id'])]
        #         # check if modname with different mass exists already
        #         for res in mod['aminoAcids']:
        #             if res not in old_mod['aminoAcids']:
        #                 old_mod['aminoAcids'].append(res)
        #
        # mod_index = 0
        # multiple_inj_list_modifications = []
        # for mod in all_mods:
        #     multiple_inj_list_modifications.append([
        #         mod_index,
        #         mod['id'],
        #         mod['mass'],
        #         ''.join(mod['aminoAcids']),
        #         mod['accession']
        #     ])
        #     mod_index += 1
        #
        # db.write_modifications(multiple_inj_list_modifications, cur, con)

        #PEPTIDE EVIDENCES
        inj_list = []
        for peptide_evidence in sequence_collection['PeptideEvidence']:
            data = [];  # peptide_ref, dBSequence_ref, start, upload_id
            data.append(peptide_evidence["peptide_ref"])  # peptide_ref att, required
            data.append(peptide_evidence["dBSequence_ref"])  # DBSequence_ref att, required
            if "start" in peptide_evidence:
                data.append(peptide_evidence["start"])  # start att, optional
            else:
                data.append(-1)

            if "isDecoy" in peptide_evidence:
                data.append(peptide_evidence["isDecoy"])  # isDecoy att, optional
            else:
                data.append("false")  # hmm, not right ToDo: fix , there's comments in Lars' code about it

            data.append(self.upload_id)

            inj_list.append(data)

        self.db.write_peptide_evidences(inj_list, self.cur, self.con)

        self.con.commit()
        mzid_reader.reset()

    @staticmethod
    def get_unimod_masses(unimod_path):
        masses = {}
        mod_id = -1

        with open(unimod_path) as f:
            for line in f:
                if line.startswith('id: '):
                    mod_id = ''.join(line.replace('id: ', '').split())

                elif line.startswith('xref: delta_mono_mass ') and not mod_id == -1:
                    mass = float(line.replace('xref: delta_mono_mass ', '').replace('"', ''))
                    masses[mod_id] = mass

        return masses


    def parse(self):
        self.logger.info('reading mzid - start')
        mzid_start_time = time()
        # schema: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0.xsd
        try:
            mzid_reader = py_mzid.MzIdentML(self.mzId_path)
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)

        self.logger.info('reading mzid - done. Time: ' + str(round(time() - mzid_start_time, 2)) + " sec")

        #
        # upload info
        #
        # upload_info_start_time = time()
        # self.logger.info('getting upload info (provider, etc) - start')
        # self.parse_upload_info(mzid_reader)
        # self.logger.info(
        #     'getting upload info - done. Time: ' + str(round(time() - upload_info_start_time, 2)) + " sec")

        #
        # Sequences, Peptides, Peptide Evidences (inc. peptide positions), Modifications
        #
        sequences_start_time = time()
        self.logger.info('getting sequences, peptides, peptide evidences, modification - start')
        self.parse_sequence_collection(mzid_reader)
        self.logger.info('getting sequences, etc - done. Time: ' + str(round(time() - sequences_start_time, 2)) + " sec")

        #
        # init spectra to protocol lookup
        #
        spectra_map_start_time = time()
        self.logger.info('generating spectra data protocol map - start')
        spectra_data_protocol_map = self.map_spectra_data_to_protocol(mzid_reader)
        # return_json['errors'] += spectra_data_protocol_map['errors']
        # ToDo: save FragmentTolerance to annotationsTable
        self.logger.info('generating spectraData_ProtocolMap - done. Time: ' + str(round(time() - spectra_map_start_time, 2)) + " sec")

        mzid_item_index = 0
        spec_id_item_index = 0
        spectra =[]
        spectrum_identifications = []

        fragment_parsing_error_scans = []

        #
        # main loop
        #
        main_loop_start_time = time()
        self.logger.info('main loop - start')

        for sid_result in mzid_reader:

            # get spectra data
            try:
                spectra_data = mzid_reader.get_by_id(sid_result['spectraData_ref'], tag_id='SpectraData', detailed=True)
                # raw_file_name = id_item['spectraData_ref'].split('/')[-1]
                # raw_file_name = re.sub('\.(mgf|mzml)', '', raw_file_name, flags=re.IGNORECASE)

            except KeyError:
                raise MzIdParseException("no spectraData_ref specified: " + sid_result['id'])

            # peak list file name
            if 'location' in spectra_data.keys():
                peak_list_file_name = ntpath.basename(spectra_data['location'])
            elif 'name' in spectra_data.keys():
                    peak_list_file_name = spectra_data['name']
            else:
                peak_list_file_name = ntpath.basename(['spectraData_ref'])

            if peak_list_file_name.endswith('raw'):
                raise MzIdParseException('.raw files not supported')

            # todo: move this stuff into init
            # missing .mgf file extension?
            if 'FileFormat' in spectra_data:
                fstring = str(json.dumps([spectra_data['SpectrumIDFormat'], spectra_data['FileFormat']])) # todo: tidy / remove
                if 'accession' in spectra_data['FileFormat']:
                    if spectra_data['FileFormat']['accession'] == 'MS:1001062' and not peak_list_file_name.lower().endswith('.mgf'):
                        self.warnings.append('location of mgf file is missing .mgf extension: ' + peak_list_file_name)
                        peak_list_file_name = peak_list_file_name + '.mgf'
                else:
                    self.warnings.append('SpectraData>FileFormat cvParam is missing accession? (Not really looks like its pytoeomics gone wrong.)')
            else:
                fstring = str(json.dumps(spectra_data['SpectrumIDFormat']))  # todo: tidy / remove
                self.warnings.append('SpectraData missing required element FileFormat.')

            self.spectrum_id_formats.add(fstring)  # todo: tidy / remove

            if self.peak_list_readers.has_key(peak_list_file_name):
                peak_list_reader = self.peak_list_readers[peak_list_file_name]
            else:
                ftp = ftplib.FTP("193.62.192.9")
                ftp.login()  # Username: anonymous password: anonymous@
                try:
                    target_dir = '/' + self.base + '/' + self.origin
                    ftp.cwd(target_dir)
                    self.logger.info('getting ' + peak_list_file_name)
                    ftp.retrbinary("RETR " + peak_list_file_name, open(self.temp_dir + '/' + peak_list_file_name, 'wb').write)
                except ftplib.error_perm as e:
                    #  check for gzipped
                    try:
                        self.logger.info('getting ' + peak_list_file_name + '.gz')
                        ftp.retrbinary("RETR " + peak_list_file_name + '.gz',
                                       open(self.temp_dir + '/' + peak_list_file_name + '.gz', 'wb').write)
                    except ftplib.error_perm as e:
                        ftp.close()
                        raise MzIdMissingFileException(type(e).__name__, peak_list_file_name, e.args)
                    ftp.close()
                    peak_list_file_name = ntpath.basename(peakListParser.unzip_peak_lists(self.temp_dir + '/' + peak_list_file_name + '.gz')[0])
                ftp.close()

                peak_list_reader = peakListParser.get_peak_list_reader(self.temp_dir + '/' + peak_list_file_name)
                self.peak_list_readers[peak_list_file_name] = peak_list_reader


            try:
                scan = peakListParser.get_scan(peak_list_reader, sid_result["spectrumID"], spectra_data['SpectrumIDFormat'])
            except Exception as e:
                raise MzIdScanNotFoundException(type(e).__name__, peak_list_file_name, sid_result["spectrumID"], e.args)

            peak_list = peakListParser.get_peak_list(scan, peak_list_reader['fileType'])
            # print(sid_result["spectrumID"])
            protocol = spectra_data_protocol_map[sid_result['spectraData_ref']]

            spectra.append([mzid_item_index, peak_list, peak_list_file_name, sid_result["spectrumID"],
                            protocol['fragmentTolerance'], self.upload_id, sid_result['id']])

            spectrum_ident_dict = dict()
            linear_index = -1  # negative index values for linear peptides

            for spec_id_item in sid_result['SpectrumIdentificationItem']:
                # get suitable id
                if 'cross-link spectrum identification item' in spec_id_item.keys():
                    self.contains_crosslinks = True
                    id = MzIdParser.get_cross_link_identifier(spec_id_item)
                else:  # assuming linear
                    # misusing 'cross-link spectrum identification item' for linear peptides with negative index
                    #specIdItem['cross-link spectrum identification item'] = linear_index
                    #spec_id_set.add(get_cross_link_identifier(specIdItem))

                    id = linear_index
                    linear_index -= 1

                # check if seen it before
                if id in spectrum_ident_dict.keys():
                    # do crosslink specific stuff
                    ident_data = spectrum_ident_dict.get(id)
                    ident_data[4] = spec_id_item['peptide_ref']
                else:
                    # do stuff common to linears and crosslinks
                    charge_state = spec_id_item['chargeState']
                    pass_threshold = spec_id_item['passThreshold']
                    scores = {
                        k: v for k, v in spec_id_item.iteritems()
                        if 'score' in k.lower() or
                           'pvalue' in k.lower() or
                           'evalue' in k.lower() or
                           'sequest' in k.lower() or
                           'scaffold' in k.lower()
                    }
                    #
                    # # fragmentation ions ToDO
                    ions = self.get_ion_types_mzid(spec_id_item)
                    # if no ion types are specified in the id file check the mzML file
                    # if len(ions) == 0 and peak_list_reader['fileType'] == 'mzml':
                    #     ions = peakListParser.get_ion_types_mzml(scan)

                    ions = list(set(ions))

                    if len(ions) == 0:
                        ions = ['peptide', 'b', 'y']
                        # ToDo: better error handling for general errors - bundling together of same type errors
                        fragment_parsing_error_scans.append(sid_result['id'])

                    ions = ';'.join(ions)

                    # extract other useful info to display
                    rank = spec_id_item['rank']

                    ident_data = [spec_id_item_index,
                                  self.upload_id,
                                  mzid_item_index,
                                  spec_id_item['peptide_ref'],
                                  '',  # pep2
                                  charge_state,
                                  rank,
                                  pass_threshold,
                                  json.dumps(ions),
                                  json.dumps(scores)]

                    spectrum_ident_dict[id] = ident_data

                    spec_id_item_index += 1

            #ToDO: better way to concat arrays, numpy?
            for spec_ident in spectrum_ident_dict.values():
                spectrum_identifications.append(spec_ident)

            mzid_item_index += 1

        # end main loop
        self.logger.info('main loop - done. Time: ' + str(round(time() - main_loop_start_time, 2)) + " sec")

        # once loop is done write remaining data to DB
        db_wrap_up_start_time = time()
        self.logger.info('write remaining entries and modifications to DB - start')
        try:
            self.db.write_spectra(spectra, self.cur, self.con)
            self.db.write_spectrum_identifications(spectrum_identifications, self.cur, self.con)
            self.con.commit()
        except psycopg2.Error as e:
            raise e

        self.logger.info('write remaining entries and modifications to DB - done. Time: '
                    + str(round(time() - db_wrap_up_start_time, 2)) + " sec")

        # fill in missing score information
        # score_fill_start_time = time()
        # logger.info('fill in missing scores - start')
        # db.fill_in_missing_scores(cur, con)
        # logger.info('fill in missing scores - done. Time: ' + str(round(time() - score_fill_start_time, 2)) + " sec")

        # warnings
        if len(fragment_parsing_error_scans) > 0:
            if len(fragment_parsing_error_scans) > 50:
                id_string = '; '.join(fragment_parsing_error_scans[:50]) + ' ...'
            else:
                id_string = '; '.join(fragment_parsing_error_scans)
                self.warnings.append({
                    "type": "IonParsing",
                    "message": "mzidentML file does not specify fragment ions.",
                    'id': id_string
                })

        parse_time = str(round(time() - mzid_start_time, 2))
        formats = []
        for fmat in self.spectrum_id_formats:
            formats.append(fmat)
        formats = json.dumps(formats, cls=NumpyEncoder)
        # store some info on how things went
        try:
            self.cur.execute("""
                UPDATE uploads SET
                    peak_list_file_names=%s, 
                    upload_warnings=%s,
                    spectrum_id_format=%s,
                    parse_time=%s,
                    contains_crosslinks=%s
                 WHERE id = %s""",
                             [
                                 json.dumps(self.peak_list_readers.keys(), cls=NumpyEncoder),
                                 json.dumps(self.warnings),
                                 formats,
                                 parse_time,
                                 self.contains_crosslinks,
                                 self.upload_id
                             ])
            self.con.commit()
        except psycopg2.Error as e:
            raise e

        self.logger.info('all done! Total time: ' + str(round(time() - mzid_start_time, 2)) + " sec")


    def upload_info(self, mzid_reader):
        # AnalysisSoftwareList - optional element
        # see https://groups.google.com/forum/#!topic/pyteomics/Mw4eUHmicyU
        mzid_reader.schema_info['lists'].add("AnalysisSoftware")
        try:
            analysis_software = json.dumps(mzid_reader.iterfind('AnalysisSoftwareList').next()['AnalysisSoftware'])
        except StopIteration:
            analysis_software = '{}'
        mzid_reader.reset()

        # Provider - optional element
        provider = '{}'
        # try:
        #     provider = json.dumps(mzid_reader.iterfind('Provider').next())
        # except StopIteration:
        #     pass
        # mzid_reader.reset()

        # AuditCollection - optional element
        audits = '{}'
        # try:
        #     audits = json.dumps(mzid_reader.iterfind('AuditCollection').next())
        # except StopIteration:
        #     audits = '{}'
        # mzid_reader.reset()

        # AnalysisSampleCollection - optional element
        samples = '{}'
        # try:
        #     samples = json.dumps(mzid_reader.iterfind('AnalysisSampleCollection').next()['Sample'])
        # except StopIteration:
        #     samples = '{}'
        # mzid_reader.reset()

        # AnalysisCollection - required element
        analyses = '{}'
        # analyses = json.dumps(mzid_reader.iterfind('AnalysisCollection').next()['SpectrumIdentification'])
        # mzid_reader.reset()

        # AnalysisProtocolCollection - required element
        protocols ='{}'
        # protocol_collection = mzid_reader.iterfind('AnalysisProtocolCollection').next()
        # protocols = protocol_collection['SpectrumIdentificationProtocol']
        protocols = json.dumps(protocols, cls=NumpyEncoder)
        # mzid_reader.reset()

        # BibliographicReference - optional element
        bibRefs = []
        # for bib in mzid_reader.iterfind('BibliographicReference'):
        #     bibRefs.append(bib)
        # bibRefs = json.dumps(bibRefs)
        # mzid_reader.reset()


        # # AnalysisSoftwareList - optional element
        # # see https://groups.google.com/forum/#!topic/pyteomics/Mw4eUHmicyU
        # mzid_reader.schema_info['lists'].add("AnalysisSoftware")
        # analysis_software = {}
        # for analysis_software_list_id in mzid_reader._offset_index["AnalysisSoftwareList"].keys():
        #     analysis_software = mzid_reader.get_by_id(analysis_software_list_id, tag_id='AnalysisSoftwareList',
        #                                               detailed=True)['AnalysisSoftware']
        # analysis_software = json.dumps(analysis_software)
        #
        # # try:
        # #     analysis_software = json.dumps(mzid_reader.iterfind('AnalysisSoftwareList').next()['AnalysisSoftware'])
        # # except StopIteration:
        # #     analysis_software = '{}'
        # # mzid_reader.reset()
        #
        # # Provider - optional element
        # try:
        #     provider = json.dumps(mzid_reader.iterfind('Provider').next())
        # except StopIteration:
        #     provider = '{}'
        # mzid_reader.reset()
        # analysis_software = {}
        # for analysis_software_list_id in mzid_reader._offset_index["AnalysisSoftwareList"].keys():
        #     analysis_software = mzid_reader.get_by_id(analysis_software_list_id, tag_id='AnalysisSoftwareList',
        #                                               detailed=True)['AnalysisSoftware']
        # analysis_software = json.dumps(analysis_software)
        #
        # # AuditCollection - optional element
        # # try:
        # #     audits = json.dumps(mzid_reader.iterfind('AuditCollection').next())
        # # except StopIteration:
        # #     audits = '{}'
        # # mzid_reader.reset()
        #
        # audits = {}
        # for audit_collection_id in mzid_reader._offset_index["AuditCollection"].keys():
        #     audits = mzid_reader.get_by_id(audit_collection_id, tag_id='AuditCollection',
        #                                               detailed=True)
        # audits = json.dumps(audits)
        #

       # upload_id = self.db.write_upload([self.user_id, filename, json.dumps(peak_list_file_names),
       #                   analysis_software, provider, audits, samples, analyses, protocols, bibRefs, 'country'],
       #                  self.cur, self.con,
       #                  )

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)