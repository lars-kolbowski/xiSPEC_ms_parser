import pyteomics.mzid as py_mzid
import pyteomics.mgf as py_mgf
# import pyteomics.mzml as py_mzml
import pymzml
import re
import json
import sys
import sqlite3
import os
import shutil
import logging
import ntpath

try:
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)


    dev = False
    unimodPath = dname+'/obo/unimod.obo'
    try:
        mzidFile = sys.argv[1]
        peakList_file = sys.argv[2]
        upload_folder = "../uploads/" + sys.argv[3]
        # unimodPath = 'obo/unimod.obo'

    except IndexError:
        dev = True
        # unimodPath = 'unimod.obo'
        baseDir = "/home/lars/work/xiSPEC/"
        mzidFile = baseDir + "DSSO_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD_CID-only.mzid"
        # mzidFile = baseDir + 'OpenxQuest_example_added_annotations.mzid'

        peakList_file = baseDir + "centroid_B170808_08_Lumos_LK_IN_90_HSA-DSSO-Sample_Xlink-CID-EThcD.mzML"
        # peakList_file = baseDir + "B170918_12_Lumos_LK_IN_90_HSA-DSSO-HCD_Rep1.mgf"

    if dev:
        logFile = "/home/lars/Xi/xiSPEC/log/parser.log"
    else:
        logFile = dname+"/log/" + sys.argv[3] + ".log"
        if not os.path.exists(logFile):
            open(logFile, 'w').close()

    if not os.path.isfile(logFile):
        os.fdopen(os.open(logFile, os.O_WRONLY | os.O_CREAT, 0o777), 'w').close()
        # open(logFile, 'a').close()

    logging.basicConfig(filename=logFile, level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)

    unimod_masses = {}

    with open(unimodPath) as f:
        for line in f:
            if line.startswith('id: '):
                id = ''.join(line.replace('id: ', '').split())

            elif line.startswith('xref: delta_mono_mass '):
                mass = float(line.replace('xref: delta_mono_mass ', '').replace('"', ''))
                unimod_masses[id] = mass
except Exception as e:
    print e

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def write_to_db(inj_list, cur):
    try:
        cur.executemany("""
INSERT INTO identifications (
    'id',
    'mzid',
    'pep1',
    'pep2',
    'linkpos1',
    'linkpos2',
    'charge',
    'passThreshold',
    'fragTolerance',
    'ionTypes',
    'crosslinker_modMass',
    'rank',
    'scores',
    'isDecoy',
    'protein',
    'file',
    'scanID',
    'peakList_id'
)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", inj_list)
    except sqlite3.Error as err:
        return err
    return True


def add_to_modlist(mod, modlist):
    if mod['name'] == "unknown_modification":
        mod['name'] = "({0:.2f})".format(mod['monoisotopicMassDelta'])

    mod['monoisotopicMassDelta'] = round(float(mod['monoisotopicMassDelta']), 6)

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


# def get_peaklist_from_mzml(scan):
#     """
#     Function to extract peaklist in xiAnnotator JSON format dict from mzml
#
#     Parameters:
#     ------------------------
#     scan, pymzml reader spectrum
#     id: scanID
#     """
#
#     if "profile spectrum" in scan.keys():
#         peaks = scan.centroidedPeaks
#     else:
#         peaks = scan.centroidedPeaks
#
#     return [{"mz": peak[0], "intensity": peak[1]} for peak in peaks if peak[1] > 0]


def get_peptide_info(sid_items, mzid_reader):

    return_dict = {
        'peptides': [],
        'linkSites': [],
        'annotation': {'modifications': []},
        'cross-linker modMass': 0
    }

    all_mods = []  # Modifications list
    mod_aliases = {
        "amidated_bs3": "bs3nh2",
        "carbamidomethyl": "cm",
        "hydrolyzed_bs3": "bs3oh",
        "oxidation": "ox"
    }

    pep_index = 0
    target_decoy = []
    for sid_item in sid_items:  # len = 1 for linear

        # Target-Decoy
        peptide_evidences = [mzid_reader.get_by_id(s['peptideEvidence_ref']) for s in sid_item['PeptideEvidenceRef']]
        # ToDo: isDecoy might not be defined. How to handle? (could make use of pyteomics.mzid.is_decoy())
        try:
            decoy = peptide_evidences[0]['isDecoy']
        except KeyError:
            decoy = None
        target_decoy.append({"peptideId": pep_index, 'isDecoy': decoy})  # TODO: multiple PeptideEvidenceRefs TD?

        # proteins
        proteins = [mzid_reader.get_by_id(p['dBSequence_ref']) for p in peptide_evidences]
        protein_accessions = [p['accession'] for p in proteins]

        # convert pepsequence to dict
        pepId = sid_item['peptide_ref']
        peptide = mzid_reader.get_by_id(pepId, detailed=True)

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
                        returnJSON['errors'].append({"type": "mzidParseError",
                                                     "message": "could not get modification mass for modification {}" % mod})
                        continue

                link_index = 0  # TODO: multilink support

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
                        mod['name'] = add_to_modlist(mod, all_mods)  # save to all mods list and get back new_name

                        if pep_seq_dict[mod_location]['Modification'] == '':
                            pep_seq_dict[mod_location]['Modification'] = mod['name']
                        else:
                            logger.error('double modification on aa')
                            logger.error(mod)
                            logger.error(pep_seq_dict[mod_location])

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

        # other parameters - should be the same for each paired sid
        return_dict['precursorCharge'] = sid_item['chargeState']
        return_dict['isDecoy'] = target_decoy
        return_dict['proteins'] = protein_accessions
        return_dict['passThreshold'] = sid_item['passThreshold']
        return_dict['scores'] = {k: v for k, v in sid_item.iteritems()
                   if 'score' in k.lower() or 'pvalue' in k.lower() or 'evalue' in k.lower()}
        for mod in all_mods:
            return_dict['annotation']['modifications'].append({
                'aminoAcids': mod['residues'],
                'id': mod['name'],
                'mass': mod['monoisotopicMassDelta']
            })

    return return_dict


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

    for analysisCollection in mzid_reader.iterfind('AnalysisCollection'):
        for spectrumIdentification in analysisCollection['SpectrumIdentification']:
            sid_protocol_ref = spectrumIdentification['spectrumIdentificationProtocol_ref']
            sid_protocol = mzidReader.get_by_id(sid_protocol_ref, detailed=True)
            frag_tol = sid_protocol['FragmentTolerance']
            try:
                frag_tol_plus = frag_tol['search tolerance plus value']
                frag_tol_value = re.sub('[^0-9,.]', '', str(frag_tol_plus['value']))
                if frag_tol_plus['unit'].lower() == 'parts per million':
                    frag_tol_unit = 'ppm'
                elif frag_tol_plus['unit'].lower() == 'dalton':
                    frag_tol_unit = 'Da'
                else:
                    frag_tol_unit = frag_tol_plus['unit']

            except KeyError:
                frag_tol_value = '10'
                frag_tol_unit = 'ppm'
                spectra_data_protocol_map['errors'].append(
                    {"type": "mzidParseError",
                     "message": "could not parse ms2tolerance. Falling back to default values."})

            if not all([
                frag_tol['search tolerance plus value']['value'] == frag_tol['search tolerance minus value']['value'],
                frag_tol['search tolerance plus value']['unit'] == frag_tol['search tolerance minus value']['unit']
            ]):
                spectra_data_protocol_map['errors'].append(
                    {"type": "mzidParseError",
                     "message": "search tolerance plus value doesn't match minus value. Using plus value!"})

            for inputSpectra in spectrumIdentification['InputSpectra']:
                spectra_data_ref = inputSpectra['spectraData_ref']

                spectra_data_protocol_map[spectra_data_ref] = {
                    'protocol_ref': sid_protocol_ref,
                    'fragmentTolerance': ' '.join([frag_tol_value, frag_tol_unit])
                }

    mzid_reader.reset()

    return spectra_data_protocol_map


def get_ion_types_mzid(sid_item):
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


def get_ion_types_mzml(mzml_scan):
    frag_methods = {
        'beam-type collision-induced dissociation': ["b", "y"],
        'collision-induced dissociation': ["b", "y"],
        'electron transfer dissociation': ["c", "z"],
    }
    # get fragMethod and translate that to Ion Types
    ion_types = []
    for key in mzml_scan.keys():
        if key in frag_methods.keys():
            ion_types += frag_methods[key]
    return ion_types


def get_cross_link_identifier(sid_item):
    # For reporting the evidence associated with the identification, within a given <SpectrumIdentificationResult>,
    # a pair of cross-linked peptides MUST be reported as two instances of <SpectrumIdentificationItem> through having
    # a shared local unique identifier as the value for the CV term "cross-link spectrum identification item" (MS:1002511).
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


if not dev:
    dbfolder = "dbs/tmp/"
    try:
        os.stat(dbfolder)
    except:
        os.mkdir(dbfolder)

try:
    if dev:
        con = sqlite3.connect('test.db')
    else:
        con = sqlite3.connect(dbfolder + sys.argv[3] + '.db')
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS identifications")
    cur.execute(
        "CREATE TABLE identifications("
        "id INT PRIMARY KEY, "
        "mzid TEXT, "
        "pep1 TEXT, "
        "pep2 TEXT, "
        "linkpos1 INT, "
        "linkpos2 INT, "
        "charge INT, "
        "passThreshold INT, "
        "fragTolerance TEXT, "
        "ionTypes TEXT, "
        "crosslinker_modMass FLOAT, "
        "rank INT, "
        "scores TEXT, "
        "isDecoy INT, "
        "protein TEXT, "
        "file TEXT, "
        "scanID INT, "
        "peakList_id INT)"
    )
    cur.execute("DROP TABLE IF EXISTS modifications")
    cur.execute(
        "CREATE TABLE modifications("
        "id INT PRIMARY KEY, "
        "name TEXT, "
        "mass FLOAT, "
        "residues TEXT)"
    )
    cur.execute("DROP TABLE IF EXISTS peakLists")
    cur.execute(
        "CREATE TABLE peakLists("
        "id INT PRIMARY KEY, "
        "peakList TEXT)"
    )


except sqlite3.Error as e:
    logger.error(e)
    print(json.dumps({"error": e.args[0]}))
    sys.exit(1)

returnJSON = {
    "response": "",
    "errors": []
}

try:
    # Identification File
    logger.info('reading mzid - start')
    # schema: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0.xsd
    mzidReader = py_mzid.MzIdentML(mzidFile)
    logger.info('reading mzid - done')

    logger.info('generating spectraData_ProtocolMap - start')
    spectraData_ProtocolMap = map_spectra_data_to_protocol(mzidReader)
    returnJSON['errors'] += spectraData_ProtocolMap['errors']
    # ToDo: save FragmentTolerance to annotationsTable
    logger.info('generating spectraData_ProtocolMap - done')

    # peakList file
    logger.info('reading peakList file - start')
    peakList_fileName = ntpath.basename(peakList_file).lower()
    if peakList_fileName.endswith('.mzml'):
        peakList_fileType = 'mzml'
        # premzml = mzml.PreIndexedMzML(mzml_file)
        # mzmlReader = py_mzml.read(peakList_file)
        pymzmlReader = pymzml.run.Reader(peakList_file)

    elif peakList_fileName.endswith('.mgf'):
        peakList_fileType = 'mgf'
        mgfReader = py_mgf.read(peakList_file)
        peakListArr = [pl for pl in mgfReader]
    logger.info('reading peakList file - done')

    mzidItem_index = 0
    specIdItem_index = 0
    multipleInjList_identifications = []
    multipleInjList_peakLists = []
    modifications = []


    # main loop
    logger.info('entering main loop')

    for mzidItem in mzidReader:  # mzidItem = mzid_reader.next()

        # make_spec_id_pairs(mzid_item['SpectrumIdentificationItem'])
        SpecIdSet = set()
        linear_index = -1  # negative index values for linear peptides

        for specIdItem in mzidItem['SpectrumIdentificationItem']:
            if 'cross-link spectrum identification item' in specIdItem.keys():
                SpecIdSet.add(get_cross_link_identifier(specIdItem))
            else:  # assuming linear
                # misusing 'cross-link spectrum identification item' for linear peptides with negative index
                specIdItem['cross-link spectrum identification item'] = linear_index
                SpecIdSet.add(get_cross_link_identifier(specIdItem))
                linear_index -= 1

        # extract scanID
        try:
            scanID = int(mzidItem['peak list scans'])
        except KeyError:
            # ToDo: this might not work for all mzids. ProteomeDiscoverer 2.2 format 'scan=xx file=xx'
            matches = re.findall("scan=([0-9]+)", mzidItem["spectrumID"])
            if len(matches) > 0:

                # ToDo: handle multiple scans? Is this standard compliant?
                # found in https://github.com/HUPO-PSI/mzIdentML/blob/master/examples/1_2examples/crosslinking/OpenxQuest_example_added_annotations.mzid
                scanIDs = [int(m) for m in matches]
                if len(scanIDs) > 1:
                    returnJSON['errors'].append(
                        {"type": "mzidParseError",
                         "message": "More than one scan found for SpectrumIdentificationItem: %s" % mzidItem["spectrumID"],
                         'id': mzidItem['id']
                        })
                    continue
                else:
                    scanID = scanIDs[0]
            else:
                returnJSON['errors'].append(
                {"type": "mzidParseError", "message": "Error parsing scanID from mzidentml!"})
                continue

        # peakList
        if peakList_fileType == 'mzml':
            try:
                scan = pymzmlReader[scanID]
            except KeyError:
                returnJSON['errors'].append(
                    {"type": "mzmlParseError",
                     "message": "requested scanID %i not found in peakList file" % scanID,
                     'id': mzidItem['id']
                    })
                continue
            if scan['ms level'] == 1:
                returnJSON['errors'].append(
                    {"type": "mzmlParseError",
                     "message": "requested scanID %i is not a MSn scan" % scanID,
                     'id': mzidItem['id']
                    })
                continue

            peakList = "\n".join(["%s %s" % (mz, i) for mz, i in scan.peaks if i > 0])
            # peaklist = get_peaklist_from_mzml(scan)

        elif peakList_fileType == 'mgf':
            try:
                scan = peakListArr[scanID]
            except IndexError:
                returnJSON['errors'].append(
                    {"type": "mzmlParseError",
                     "message": "requested scanID %i not found in peakList file" % scanID,
                     'id': mzidItem['id']
                    })
                continue
            peaks = zip(scan['m/z array'], scan['intensity array'])
            peakList = "\n".join(["%s %s" % (mz, i) for mz, i in peaks if i > 0])

        multipleInjList_peakLists.append([mzidItem_index, peakList])

        ms2Tol = spectraData_ProtocolMap[mzidItem['spectraData_ref']]['fragmentTolerance']

        # alternatives = []
        for SpecId in SpecIdSet:
            paired_specIdItems = [sid_item for sid_item in mzidItem['SpectrumIdentificationItem'] if
                                  sid_item['cross-link spectrum identification item'] == SpecId]

            # if len(paired_specIdItems) > 2:
            #     returnJSON['errors'].append({
            #         "type": "PeptideParseError",
            #         "message": "more than 2 peptides with the same cross-link id found",
            #         'id': mzidItem['id']
            #     })
            #     continue

            pep_info = get_peptide_info(paired_specIdItems, mzidReader)


            # fragmentation ions
            pep_info['ions'] = get_ion_types_mzid(paired_specIdItems[0])
            # if no ion types are specified in the id file check the mzML file
            if len(pep_info['ions']) == 0 and peakList_fileType == 'mzml':
                pep_info['ions'] = get_ion_types_mzml(scan)

            pep_info['ions'] = list(set(pep_info['ions']))

            if len(pep_info['ions']) == 0:
                pep_info['ions'] = ['peptide', 'b', 'y']
                returnJSON['errors'].append(
                    {"type": "IonParsing", "message": "could not parse fragment ions assuming precursor-, b- and y-ion", 'id': mzidItem['id']})

            pep_info['ions'] = ';'.join(pep_info['ions'])

                # extract other useful info to display
            rank = paired_specIdItems[0]['rank']

            # ToDo: handling for mzid that don't include isDecoy
            isDecoy = any([pep['isDecoy'] for pep in pep_info['isDecoy']])
            accessions = ";".join(pep_info['proteins'])

            # raw file name
            try:
                rawFileName = mzidItem['spectraData_ref']
            except KeyError:
                returnJSON['errors'].append(
                    {"type": "mzidParseError", "message": "no spectraData_ref specified", 'id': mzidItem['id']})
                rawFileName = ""

            try:
                rawFileName = path_leaf(mzidReader.get_by_id(rawFileName)['location'])
            except KeyError:
                pass

            # passThreshold
            if pep_info['passThreshold']:
                passThreshold = 1
            else:
                passThreshold = 0

            # peps and linkpos
            pep1 = pep_info['peptides'][0]

            if len(pep_info['peptides']) > 1:
                linkpos1 = pep_info['linkSites'][0]
                pep2 = pep_info['peptides'][1]
                linkpos2 = pep_info['linkSites'][1]

            else:
                pep2 = ""
                linkpos1 = -1
                linkpos2 = -1

            multipleInjList_identifications.append(
                [specIdItem_index,
                 mzidItem['id'],
                 pep1,
                 pep2,
                 linkpos1,
                 linkpos2,
                 pep_info['precursorCharge'],
                 passThreshold,
                 ms2Tol,
                 pep_info['ions'],
                 pep_info["cross-linker modMass"],
                 rank,
                 json.dumps(pep_info['scores']),
                 isDecoy,
                 accessions,
                 rawFileName,
                 scanID,
                 mzidItem_index]
            )

            # add mods to global modList
            for mod in pep_info['annotation']['modifications']:
                if mod['id'] not in [m['id'] for m in modifications]:
                    modifications.append(mod)
                else:
                    old_mod = modifications[[m['id'] for m in modifications].index(mod['id'])]
                    # check if modname with different mass exists already
                    for res in mod['aminoAcids']:
                        if res not in old_mod['aminoAcids']:
                            old_mod['aminoAcids'].append(res)

            specIdItem_index += 1

        mzidItem_index += 1

        if specIdItem_index % 500 == 0:
            logger.info('writing 500 entries to DB')
            write_to_db(multipleInjList_identifications, cur)
            multipleInjList_identifications = []

            cur.executemany("""INSERT INTO peakLists ('id', 'peaklist') VALUES (?, ?)""", multipleInjList_peakLists)
            multipleInjList_peakLists = []

            con.commit()
            # if dev:
            #     break

    # once its done submit the last reqs to DB
    logger.info('writing remaining entries to DB')
    if len(multipleInjList_identifications) > 0:
        write_to_db(multipleInjList_identifications, cur)
        multipleInjList_identifications = []

        cur.executemany("""INSERT INTO peakLists ('id', 'peaklist') VALUES (?, ?)""", multipleInjList_peakLists)
        multipleInjList_peakLists = []

        # modifications
        mod_index = 0
        multipleInjList_modifications = []
        for mod in modifications:
            multipleInjList_modifications.append([mod_index, mod['id'], mod['mass'], ''.join(mod['aminoAcids'])])
            mod_index += 1
        cur.executemany("""INSERT INTO modifications ('id', 'name', 'mass', 'residues') VALUES (?, ?, ?, ?)""", multipleInjList_modifications)

        con.commit()

    # delete uploaded files after they have been parsed
    if not dev:
        logger.info('deleting uploaded files')
        shutil.rmtree(upload_folder)

    if len(returnJSON["errors"]) > 0:
        returnJSON['response'] = "Warning: %i error(s) occured!" % len(
            returnJSON['errors'])
        for e in returnJSON['errors']:
            logger.error(e)
    else:
        returnJSON['response'] = "No errors, smooth sailing!"

    print(json.dumps(returnJSON))
    if con:
        con.close()
        logger.info('all done!')

except Exception as e:
    logger.exception(e)
    returnJSON['errors'].append(
        {"type": "Error", "message": e})
