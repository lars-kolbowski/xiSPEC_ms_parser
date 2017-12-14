import pyteomics.mgf as py_mgf
import pymzml
import pandas as pd
import ntpath
import re
import xiSPEC_sqlite as db
import json

def parse(csv_file, peak_list_file_list, cur, con, logger):

    return_json = {
        "response": "",
        "modifications": [],
        "errors": []
    }

    logger.info('reading csv - start')
    # schema: https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0.xsd
    id_df = pd.read_csv(csv_file)
    id_df.columns = [x.lower() for x in id_df.columns]

    required_cols = ['id', 'scannumber', 'charge', 'pepseq 1', 'protein 1']

    for header in required_cols:
        if header not in id_df.columns:
            return_json['errors'].append({
                "type": "csvParseError",
                "message": "Required csv column %s missing" % header,
            })
            return return_json

    logger.info('reading csv - done')

    # unimod_masses = get_unimod_masses(unimod_path)

    multiple_inj_list_identifications = []
    multiple_inj_list_peak_lists = []
    # modifications = []

    # peakList file
    peak_list_readers = {}
    for peak_list_file in peak_list_file_list:
        logger.info('reading peakList file - start')
        peak_list_file_name = ntpath.basename(peak_list_file)
        if peak_list_file_name.lower().endswith('.mzml'):
            peak_list_file_type = 'mzml'
            peak_list_readers[peak_list_file_name] = pymzml.run.Reader(peak_list_file)

        elif peak_list_file_name.lower().endswith('.mgf'):
            peak_list_file_type = 'mgf'
            mgf_reader = py_mgf.read(peak_list_file)
            scans = [s for s in mgf_reader]
            peak_list_readers[peak_list_file_name] = scans

        logger.info('reading peakList file - done')

    # main loop
    logger.info('entering main loop')

    for id_item_index, id_item in id_df.iterrows():  # id_item_index, id_item = id_df.iterrows().next()

        # spec_id_set = set()

        # extract scan_id
        # try:
        scan_id = int(id_item['scannumber'])

        # raw file name
        try:
            raw_file_name = id_item['runname']
        except KeyError:
            raw_file_name = ''

        # peakList ToDo: duplicate code blocks - needs refactoring
        if peak_list_file_type == 'mzml':
            try:
                pl_reader = peak_list_readers[raw_file_name]
            except KeyError:
                if len(peak_list_readers.keys()) == 1:
                    pl_reader = peak_list_readers[peak_list_readers.keys()[0]]
                else:
                    return_json['errors'].append({
                        "type": "mzidParseError",
                        "message": "peak list filename %s from csv does not match any of your peaklist files" % raw_file_name,
                        'id': id_item['id']
                    })
                    continue
            try:
                scan = pl_reader[scan_id]
            except IndexError:
                return_json['errors'].append({
                    "type": "mzmlParseError",
                    "message": "requested scanID %i not found in peakList file" % scan_id,
                    'id': id_item['id']
                })
                continue
            if scan['ms level'] == 1:
                return_json['errors'].append({
                    "type": "mzmlParseError",
                    "message": "requested scanID %i is not a MSn scan" % scan_id,
                    'id': id_item['id']
                })
                continue

            peak_list = "\n".join(["%s %s" % (mz, i) for mz, i in scan.peaks if i > 0])
            # peak_list = get_peaklist_from_mzml(scan)

        elif peak_list_file_type == 'mgf':
            try:
                pl_reader = peak_list_readers[raw_file_name]
            except KeyError:
                if len(peak_list_readers.keys()) == 1:
                    pl_reader = peak_list_readers[peak_list_readers.keys()[0]]
                else:
                    return_json['errors'].append({
                        "type": "mzidParseError",
                        "message": "peak list filename %s from csv does not match any of your peaklist files" % raw_file_name,
                        'id': id_item['id']
                    })
                    continue
            try:
                scan = pl_reader[scan_id]
            except IndexError:
                return_json['errors'].append({
                    "type": "mzmlParseError",
                    "message": "requested scanID %i not found in peakList file" % scan_id,
                    'id': id_item['id']
                })
                continue

            peaks = zip(scan['m/z array'], scan['intensity array'])
            peak_list = "\n".join(["%s %s" % (mz, i) for mz, i in peaks if i > 0])

        multiple_inj_list_peak_lists.append([id_item_index, peak_list])

        # identification id
        identification_id = id_item['id']

        # rank
        try:
            rank = id_item['rank']
        except KeyError:
            rank = 1

        # peptides and link positions
        pep1 = id_item['pepseq 1']
        try:
            # ToDo: improve error handling for cl peptides
            pep2 = id_item['pepseq 2']
            linkpos1 = id_item['linkpos 1'] - 1
            linkpos2 = id_item['linkpos 2'] - 1
            cl_mod_mass = id_item['crosslinkermodmass']
        except KeyError:
            # linear
            pep2 = ""
            linkpos1 = -1
            linkpos2 = -1
            cl_mod_mass = 0

        # charge
        charge = id_item['charge']

        # passThreshold
        try:
            pass_threshold = 1 if id_item['passthreshold'] else 0
        except KeyError:
            pass_threshold = 1

        # fragment tolerance
        frag_tol = id_item['fragmenttolerance'].split(' ', 1)
        if frag_tol[1].lower() == 'parts per million':
            frag_tol[1] = 'ppm'
        elif frag_tol[1].lower() == 'dalton':
            frag_tol[1] = 'Da'
        if frag_tol[1] not in ['ppm', 'Da']:
            return_json['errors'].append({
                "type": "fragTolParseError",
                "message": "unknown fragment tolerance unit: %s\nSupported values are: ppm, Da" % frag_tol[1],
                'id': id_item['id']
            })
            continue
        ms2_tol = ' '.join(frag_tol)

        # ion types
        ion_types = id_item['iontypes'].lower()
        unknown_ions = [ion for ion in ion_types.split(';') if ion not in ['peptide', 'a', 'b', 'c', 'x', 'y', 'z']]
        if len(unknown_ions) > 0:
            return_json['errors'].append({
                "type": "ionTypeParseError",
                "message": "unknown ion(s) : %s\nSupported values are: peptide, a, b, c, x, y, z" % ';'.join(unknown_ions),
                'id': id_item['id']
            })
        # ToDo: could check against mzml fragmentation type and display warning if ions don't match

        # score
        score = json.dumps({'score': id_item['score']})

        # isDecoy
        try:
            is_decoy = 0 if id_item['isdecoy'] else 1
        except KeyError:
            is_decoy = 0

        # protein
        protein1 = id_item['protein 1']
        try:
            protein2 = id_item['protein 2']
        except KeyError:
            protein2 = ''

        # create entry
        multiple_inj_list_identifications.append(
            [id_item_index,
             identification_id,
             pep1,
             pep2,
             linkpos1,
             linkpos2,
             charge,
             pass_threshold,
             ms2_tol,
             ion_types,
             cl_mod_mass,
             rank,
             score,
             is_decoy,
             protein1,
             protein2,
             raw_file_name,
             scan_id,
             id_item_index]
        )

        modifications = re.search('([^A-Z]+)', ''.join([pep1, pep2])).groups()
        for mod in modifications:
            if mod not in return_json['modifications']:
                return_json['modifications'].append(mod)

        #  write to DB
        if id_item_index % 500 == 0:
            logger.info('writing 500 entries to DB')
            try:
                multiple_inj_list_identifications = db.write_identifications(multiple_inj_list_identifications, cur, con)
                multiple_inj_list_peak_lists = db.write_peaklists(multiple_inj_list_peak_lists, cur, con)

            except db.DBException as e:
                return_json['errors'].append(
                    {"type": "dbError",
                     "message": e.args[0],
                     'id': id_item_index
                     })
                return return_json

            # commit changes
            con.commit()

    # once loop is done write remaining data to DB
    logger.info('writing remaining entries to DB')
    try:
        db.write_identifications(multiple_inj_list_identifications, cur, con)
        db.write_peaklists(multiple_inj_list_peak_lists, cur, con)

    except db.DBException as e:
        return_json['errors'].append(
            {"type": "dbError",
             "message": e.args[0],
             'id': id_item_index
             })
        return return_json

    return return_json
