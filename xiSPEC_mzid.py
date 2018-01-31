import pyteomics.mzid as py_mzid
import re
import ntpath
import json
import sys
from time import time
import xiSPEC_peakList as peakListParser

try:
    if sys.argv[4] == "pg":
        import xiUI_pg as db
    else:
        import xiSPEC_sqlite as db
except IndexError:
    import xiSPEC_sqlite as db


def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def get_ion_types_mzid(sid_item, logger):
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
            logger.info(e, ion_name)
            continue

    return ion_types


# ToDo: clear confusion about 0 & 1 based formats
def get_scan_id(spec_id, spec_id_format):
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


def map_seq_ref_to_protein(mzid_reader, cur, con):
    """
    extract and map - which includes data like -
     ToDo: improve error handling

    Parameters:
    ------------------------
    mzid_reader: pyteomics mzid_reader
    """

    seq_ref_prot_map = {
        'errors': [],
    }


    sequence_collection = mzid_reader.iterfind('SequenceCollection').next()

    inj_list = [];
    for db_sequence in sequence_collection['DBSequence']:
        seq_ref_prot_map[db_sequence["id"]] = db_sequence["accession"]
        inj_list.append(db_sequence)
    db.write_db_sequences(inj_list, cur, con)

    inj_list = [];
    for peptide in sequence_collection['Peptide']:
        inj_list.append(peptide)
    db.write_peptides(inj_list, cur, con)

    inj_list = [];
    for peptide_evidence in sequence_collection['PeptideEvidence']:
        inj_list.append(peptide_evidence)
    db.write_peptide_evidences(inj_list, cur, con)

    mzid_reader.reset()

    # takes more processing time but uses less memory ~15 Mb difference from memory snapshots on large Tmuris example
    # this might not be accurate way to measure memory...
    # for db_sequence_id in mzid_reader._offset_index["DBSequence"].keys():
    #     db_sequence = mzid_reader.get_by_id(db_sequence_id)
    #     try:
    #         seq_ref_prot_map[db_sequence_id] = db_sequence['accession']
    #     except KeyError:
    #         seq_ref_prot_map[db_sequence_id] = db_sequence['name']

    return seq_ref_prot_map


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


def add_to_modlist(mod, modlist):
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

    mod['monoisotopicMassDelta'] = round(float(mod['monoisotopicMassDelta']), 6)

    mod['residues'] = [aa for aa in mod['residues']]

    if mod['name'] in [m['name'] for m in modlist]:
        old_mod = modlist[[m['name'] for m in modlist].index(mod['name'])]
        # check if modname with different mass exists already
        if mod['monoisotopicMassDelta'] != old_mod['monoisotopicMassDelta']:
            mod['name'] += "*"
            add_to_modlist(mod, modlist)
        else:
            for res in mod['residues']:
                if res not in old_mod['residues']:
                    old_mod['residues'].append(res)
    else:
        modlist.append(mod)

    return mod['name']


def get_peptide_info(sid_items, mzid_reader, unimod_masses, seq_ref_protein_map, logger):

    return_dict = {
        'peptides': [],
        'linkSites': [],
        'annotation': {'modifications': []},
        'cross-linker modMass': 0,
        'errors': []
    }

    all_mods = []  # Modifications list
    mod_aliases = {
        # "amidated_bs3": "bs3nh2",
        # "carbamidomethyl": "cm",
        # "hydrolyzed_bs3": "bs3oh",
        # "oxidation": "ox"
    }

    proteins = []

    pep_index = 0
    target_decoy = []
    for sid_item in sid_items:  # len = 1 for linear

        # Target-Decoy
        peptide_evidences = [mzid_reader.get_by_id(s['peptideEvidence_ref'], tag_id='PeptideEvidence') for s in sid_item['PeptideEvidenceRef']]
        # ToDo: isDecoy might not be defined. How to handle? (could make use of pyteomics.mzid.is_decoy())
        try:
            decoy = peptide_evidences[0]['isDecoy']
        except KeyError:
            decoy = None
        target_decoy.append({"peptideId": pep_index, 'isDecoy': decoy})  # TODO: multiple PeptideEvidenceRefs TD?

        # proteins
        proteins.append(seq_ref_protein_map[peptide_evidences[0]['dBSequence_ref']])

        # convert pepsequence to dict
        pepId = sid_item['peptide_ref']
        peptide = mzid_reader.get_by_id(pepId, tag_id='Peptide', detailed=True)

        pep_seq_dict = []
        for aa in peptide['PeptideSequence']:
            pep_seq_dict.append({"Modification": "", "aminoAcid": aa})

        # add in modifications
        if 'Modification' in peptide.keys():
            for mod in peptide['Modification']:

                if 'monoisotopicMassDelta' not in mod.keys():
                    try:
                        mod['monoisotopicMassDelta'] = unimod_masses[mod['accession']]

                    except KeyError:
                        return_dict['errors'].append({
                            "type": "mzidParseError",
                            "message": "could not get modification mass for modification {}" % mod,
                            "id": sid_item
                        })
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
                            cur_mod_mass = [x['monoisotopicMassDelta'] for x in all_mods if x['name'] == cur_mod['Modification']][0]
                            mod['monoisotopicMassDelta'] += cur_mod_mass

                        mod['name'] = add_to_modlist(mod, all_mods)  # save to all mods list and get back new_name
                        cur_mod['Modification'] = mod['name']

                # error handling for mod without name
                else:
                    # cross-link acceptor doesn't have a name
                    if 'cross-link acceptor' not in mod.keys():
                        logger.error('modification without name!')
                        logger.error(mod)

                # add CL locations
                if 'cross-link donor' in mod.keys() or 'cross-link acceptor' in mod.keys():
                    return_dict['linkSites'].append(mod_location - 1)
                    # return_dict['linkSites'].append(
                    #     {"id": link_index, "peptideId": pep_index, "linkSite": mod_location - 1})
                if 'cross-link donor' in mod.keys():
                    return_dict["cross-linker modMass"] = round(mod['monoisotopicMassDelta'], 6)

        pep_index += 1

        peptide_seq_with_mods = ''.join([''.join([x['aminoAcid'], x['Modification']]) for x in pep_seq_dict])
        return_dict['peptides'].append(peptide_seq_with_mods)

        # ToDo: use searchDatabase_ref
        try:
            protein_accessions = proteins
            #protein_accessions = [p['accession'] for p in proteins]
        except KeyError:
            return_dict['errors'].append({
                "type": "mzidParseError",
                "message": "could not get protein",
                "id": sid_item
            })
            protein_accessions = ['']

        # other parameters - should be the same for each paired sid
        return_dict['precursorCharge'] = sid_item['chargeState']
        return_dict['isDecoy'] = target_decoy
        return_dict['protein1'] = protein_accessions[0]
        try:
            return_dict['protein2'] = protein_accessions[1]
        except IndexError:
            return_dict['protein2'] = ''
        return_dict['passThreshold'] = sid_item['passThreshold']
        return_dict['scores'] = {
            k: v for k, v in sid_item.iteritems()
            if 'score' in k.lower() or
             'pvalue' in k.lower() or
             'evalue' in k.lower() or
             'sequest' in k.lower() or
             'scaffold' in k.lower()
        }
        for mod in all_mods:
            try:
                mod_accession = mod['accession']
            except KeyError:
                mod_accession = ''
            return_dict['annotation']['modifications'].append({
                'aminoAcids': mod['residues'],
                'id': mod['name'],
                'mass': mod['monoisotopicMassDelta'],
                'accession': mod_accession
            })

    return return_dict


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


def get_analysis_software(mzid_file):
    # hack to get analysis software name
    mzid_stream = open(mzid_file, 'r')
    file_start = mzid_stream.read(5000)  # read by character
    mzid_stream.close()
    r = re.finditer('<SoftwareName>.*?<cvParam.*?name="(.*?)"', file_start, re.DOTALL)
    analysis_software = []
    for i in r:
        analysis_software.append(i.group(1))

    return analysis_software


def parse(mzid_file, peak_list_file_list, unimod_path, cur, con, logger):
    logger.info('reading mzid - start')
    mzid_start_time = time()

    analysis_software = get_analysis_software(mzid_file)

    return_json = {
        "response": "",
        "modifications": [],
        "errors": [],
        "warnings": [],
        "analysisSoftware": analysis_software
    }

    # schema: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0.xsd
    try:
        mzid_reader = py_mzid.MzIdentML(mzid_file)
    except Exception as e:
        return_json['errors'].append({
            "type": "mzidParseError",
            "message": e
        })
        return return_json

    logger.info('reading mzid - done. Time: ' + str(round(time() - mzid_start_time, 2)) + " sec")

    unimod_masses = get_unimod_masses(unimod_path)

    spectra_map_start_time = time()
    logger.info('generating spectra data protocol map - start')
    spectra_data_protocol_map = map_spectra_data_to_protocol(mzid_reader)
    return_json['errors'] += spectra_data_protocol_map['errors']
    # ToDo: save FragmentTolerance to annotationsTable
    logger.info('generating spectraData_ProtocolMap - done. Time: ' + str(round(time() - spectra_map_start_time, 2)) + " sec")

    protein_map_start_time = time()
    logger.info('generating dBSequence to protein map - start')
    seq_ref_protein_map = map_seq_ref_to_protein(mzid_reader, cur, con)
    return_json['errors'] += seq_ref_protein_map['errors']
    # ToDo: save FragmentTolerance to annotationsTable
    logger.info('generating dBSequence to protein map - done. Time: ' + str(round(time() - protein_map_start_time, 2)) + " sec")


    mzid_item_index = 0
    spec_id_item_index = 0
    multiple_inj_list_identifications = []
    multiple_inj_list_peak_lists = []
    modifications = []

    # peakList readers
    peak_list_start_time = time()
    logger.info('reading peakList files - start')
    peak_list_readers = peakListParser.create_peak_list_readers(peak_list_file_list)
    logger.info('reading peakList files - done. Time: ' + str(round(time() - peak_list_start_time, 2)) + " sec")

    # ToDo: better error handling for general errors - bundling errors of same type errors together
    fragment_parsing_error_scans = []
    scan_not_found_error = {}

    # main loop
    # main_loop_start_time = time()
    # logger.info('main loop - start')
    #
    # for id_item in mzid_reader:  # mzid_item = mzid_reader.next()
    #
    #     # make_spec_id_pairs(mzid_item['SpectrumIdentificationItem'])
    #     spec_id_set = set()
    #     linear_index = -1  # negative index values for linear peptides
    #
    #     for specIdItem in id_item['SpectrumIdentificationItem']:
    #         if 'cross-link spectrum identification item' in specIdItem.keys():
    #             spec_id_set.add(get_cross_link_identifier(specIdItem))
    #         else:  # assuming linear
    #             # misusing 'cross-link spectrum identification item' for linear peptides with negative index
    #             specIdItem['cross-link spectrum identification item'] = linear_index
    #             spec_id_set.add(get_cross_link_identifier(specIdItem))
    #             linear_index -= 1
    #
    #     # get spectra data
    #     try:
    #         spectra_data = mzid_reader.get_by_id(id_item['spectraData_ref'], tag_id='SpectraData', detailed=True)
    #         # raw_file_name = id_item['spectraData_ref'].split('/')[-1]
    #         # raw_file_name = re.sub('\.(mgf|mzml)', '', raw_file_name, flags=re.IGNORECASE)
    #
    #     except KeyError:
    #         return_json['errors'].append({
    #             "type": "mzidParseError",
    #             "message": "no spectraData_ref specified",
    #             'id': id_item['id']
    #         })
    #
    #     # get scan id ToDo: clear up 1/0-based confusion
    #     # scan_id = get_scan_id(id_item["spectrumID"], spectra_data['SpectrumIDFormat'])
    #     try:
    #         scan_id = int(id_item['peak list scans'])
    #     except KeyError:
    #         matches = re.findall("([0-9]+)", id_item["spectrumID"])
    #         if len(matches) > 1:
    #             # ToDo: this might not work for all mzids. Check more file formats. 0 vs 1 based mess
    #             matches = re.findall("(?:scan|index|query|mzMLid)?=?([0-9]+)", id_item["spectrumID"])
    #         if len(matches) > 0:
    #             # ToDo: handle multiple scans? Is this standard compliant?
    #             # found in https://github.com/HUPO-PSI/mzIdentML/blob/master/examples/1_2examples/crosslinking/OpenxQuest_example_added_annotations.mzid
    #             scan_ids = [int(m) for m in matches]
    #             if len(scan_ids) > 1:
    #                 return_json['errors'].append(
    #                     {"type": "mzidParseError",
    #                      "message": "More than one scan found for SpectrumIdentificationItem: %s"
    #                                 % id_item["spectrumID"],
    #                      'id': id_item['id']
    #                      })
    #                 continue
    #             else:
    #                 scan_id = scan_ids[0]
    #         else:
    #             return_json['errors'].append({
    #                 "type": "mzidParseError",
    #                 "message": "Error parsing scanID from mzidentml: %s" % id_item["spectrumID"],
    #                 "id": id_item['id']
    #             })
    #             continue
    #
    #     # raw file name
    #     if 'name' in spectra_data.keys():
    #         raw_file_name = spectra_data['name']
    #     elif 'location' in spectra_data.keys():
    #         raw_file_name = spectra_data['location'].split('/')[-1]
    #     else:
    #         raw_file_name = id_item['spectraData_ref'].split('/')[-1]
    #
    #     raw_file_name = re.sub('\.(mgf|mzml)', '', raw_file_name, flags=re.IGNORECASE)
    #
    #     # peak list
    #     try:
    #         peak_list_reader = peakListParser.get_reader(peak_list_readers, raw_file_name)
    #     except peakListParser.ParseError as e:
    #         return_json['errors'].append({
    #             "type": "peakListParseError",
    #             "message": e.args[0],
    #             'id': id_item['id']
    #         })
    #         continue
    #     try:
    #         scan = peakListParser.get_scan(peak_list_reader, scan_id)
    #     except peakListParser.ParseError:
    #         try:
    #             scan_not_found_error[raw_file_name].append(scan_id)
    #         except KeyError:
    #             scan_not_found_error[raw_file_name] = [scan_id]
    #         continue
    #
    #     peak_list = peakListParser.get_peak_list(scan, peak_list_reader['fileType'])
    #
    #     multiple_inj_list_peak_lists.append([mzid_item_index, peak_list])
    #
    #     # ms2 tolerance
    #     ms2_tol = spectra_data_protocol_map[id_item['spectraData_ref']]['fragmentTolerance']
    #
    #     # alternatives = []
    #     for SpecId in spec_id_set:
    #         paired_spec_id_items = [sid_item for sid_item in id_item['SpectrumIdentificationItem'] if
    #                                 sid_item['cross-link spectrum identification item'] == SpecId]
    #
    #         # if len(paired_specIdItems) > 2:
    #         #     returnJSON['errors'].append({
    #         #         "type": "PeptideParseError",
    #         #         "message": "more than 2 peptides with the same cross-link id found",
    #         #         'id': mzidItem['id']
    #         #     })
    #         #     continue
    #         pep_info = {
    #             'peptides': [],
    #             'linkSites': [],
    #             'annotation': {'modifications': []},
    #             'cross-linker modMass': 0,
    #             'errors': []
    #         }
    #         pep_info['precursorCharge'] = 0#sid_item['chargeState']
    #         pep_info['isDecoy'] = []#target_decoy
    #         pep_info['protein1'] = "que?"  # protein_accessions[0]
    #         pep_info['protein2'] = "que?2"  # protein_accessions[0]
    #         pep_info['passThreshold'] = 0  # protein_accessions[0]
    #         pep_info['scores'] = []  # protein_accessions[0]
    #
    #         # fragmentation ions
    #         pep_info['ions'] = get_ion_types_mzid(paired_spec_id_items[0], logger)
    #         # if no ion types are specified in the id file check the mzML file
    #         if len(pep_info['ions']) == 0 and peak_list_reader['fileType'] == 'mzml':
    #             pep_info['ions'] = peakListParser.get_ion_types_mzml(scan)
    #
    #         pep_info['ions'] = list(set(pep_info['ions']))
    #
    #         if len(pep_info['ions']) == 0:
    #             pep_info['ions'] = ['peptide', 'b', 'y']
    #             # ToDo: better error handling for general errors - bundling together of same type errors
    #             fragment_parsing_error_scans.append(id_item['id'])
    #
    #         pep_info['ions'] = ';'.join(pep_info['ions'])
    #
    #         # extract other useful info to display
    #         rank = paired_spec_id_items[0]['rank']
    #
    #         # ToDo: handling for mzid that don't include isDecoy
    #         is_decoy = any([pep['isDecoy'] for pep in pep_info['isDecoy']])
    #         # accessions = ";".join(pep_info['proteins'])
    #         protein1 = pep_info['protein1']
    #         protein2 = pep_info['protein2']
    #
    #         # passThreshold
    #         if pep_info['passThreshold']:
    #             pass_threshold = 1
    #         else:
    #             pass_threshold = 0
    #
    #         # peptides and linker position
    #         pep1 = "nopep"#pep_info['peptides'][0]
    #
    #         if len(pep_info['peptides']) > 1:
    #             link_pos1 = pep_info['linkSites'][0]
    #             pep2 = pep_info['peptides'][1]
    #             link_pos2 = pep_info['linkSites'][1]
    #
    #         else:
    #             pep2 = ""
    #             link_pos1 = -1
    #             link_pos2 = -1
    #
    #         multiple_inj_list_identifications.append(
    #             [spec_id_item_index,
    #              id_item['id'],
    #              pep1,
    #              pep2,
    #              link_pos1,
    #              link_pos2,
    #              pep_info['precursorCharge'],
    #              pass_threshold,
    #              ms2_tol,
    #              pep_info['ions'],
    #              pep_info["cross-linker modMass"],
    #              rank,
    #              # score,
    #              json.dumps(pep_info['scores']),
    #              is_decoy,
    #              protein1,
    #              protein2,
    #              raw_file_name,
    #              scan_id,
    #              mzid_item_index]
    #         )
    #
    #         # add mods to global modList
    #         for mod in pep_info['annotation']['modifications']:
    #             if mod['id'] not in [m['id'] for m in modifications]:
    #                 modifications.append(mod)
    #             else:
    #                 old_mod = modifications[[m['id'] for m in modifications].index(mod['id'])]
    #                 # check if modname with different mass exists already
    #                 for res in mod['aminoAcids']:
    #                     if res not in old_mod['aminoAcids']:
    #                         old_mod['aminoAcids'].append(res)
    #
    #         spec_id_item_index += 1
    #
    #     mzid_item_index += 1
    #
    #     if spec_id_item_index % 1000 == 0:
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
    #                  "message": e.message,
    #                  'id': spec_id_item_index
    #                  })
    #             return return_json
    #
    #         # commit changes
    #         con.commit()
    #
    # # end main loop
    # logger.info('main loop - done. Time: ' + str(round(time() - main_loop_start_time, 2)) + " sec")

    # once loop is done write remaining data to DB
    db_wrap_up_start_time = time()
    logger.info('write remaining entries and modifications to DB - start')
    try:
        db.write_identifications(multiple_inj_list_identifications, cur, con)
        db.write_peaklists(multiple_inj_list_peak_lists, cur, con)

        # modifications
        mod_index = 0
        multiple_inj_list_modifications = []
        for mod in modifications:
            multiple_inj_list_modifications.append([
                mod_index,
                mod['id'],
                mod['mass'],
                ''.join(mod['aminoAcids']),
                mod['accession']
            ])
            mod_index += 1
        db.write_modifications(multiple_inj_list_modifications, cur, con)

    except db.DBException as e:
        return_json['errors'].append(
            {"type": "dbError",
             "message": e.message,
             "id": spec_id_item_index
             })
        return return_json

    logger.info('write remaining entries and modifications to DB - done. Time: '
                + str(round(time() - db_wrap_up_start_time, 2)) + " sec")

    # fill in missing score information
    score_fill_start_time = time()
    logger.info('fill in missing scores - start')
    db.fill_in_missing_scores(cur, con)
    logger.info('fill in missing scores - done. Time: ' + str(round(time() - score_fill_start_time, 2)) + " sec")

    # multi error handler
    if len(fragment_parsing_error_scans) > 0:
        if len(fragment_parsing_error_scans) > 50:
            id_string = '; '.join(fragment_parsing_error_scans[:50]) + ' ...'
        else:
            id_string = '; '.join(fragment_parsing_error_scans)
        return_json['warnings'].append({
            "type": "IonParsing",
            "message": "could not parse fragment ions assuming precursor-, b- and y-ion",
            'id': id_string
        })

    for pl_file, scan_id_list in scan_not_found_error.iteritems():
        return_json['errors'].append({
            "type": "",
            "message": "requested scanID(s) not found in peakList file %s" % pl_file,
            'id': '; '.join([str(scan_id) for scan_id in scan_id_list])
        })

    return return_json
