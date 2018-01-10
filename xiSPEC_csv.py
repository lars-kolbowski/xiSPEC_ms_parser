import pandas as pd
import re
import json
import sys
import xiSPEC_peakList as peakListParser
try:
    if sys.argv[4] == "pg":
        import xiUI_pg as db
    else:
        import xiSPEC_sqlite as db
except IndexError:
    import xiSPEC_sqlite as db


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

    # peakList readers
    logger.info('reading peakList files - start')
    peak_list_readers = peakListParser.create_peak_list_readers(peak_list_file_list)
    logger.info('reading peakList files - done')

    scan_not_found_error = {}

    # main loop
    logger.info('entering main loop')

    for id_item_index, id_item in id_df.iterrows():  # id_item_index, id_item = id_df.iterrows().next()

        # spec_id_set = set()

        # extract scan_id
        # try:
        scan_id = int(id_item['scannumber'])

        # raw file name
        try:
            raw_file_name = id_item['runname'].split('/')[-1]
            raw_file_name = re.sub('\.(mgf|mzml)', '', raw_file_name, flags=re.IGNORECASE)

        except KeyError:
            raw_file_name = ''

        # peakList
        try:
            peak_list_reader = peakListParser.get_reader(peak_list_readers, raw_file_name)
        except peakListParser.ParseError as e:
            return_json['errors'].append({
                "type": "peakListParseError",
                "message": e.args[0],
                'id': id_item['id']
            })
            continue
        try:
            scan = peakListParser.get_scan(peak_list_reader, scan_id)
        except peakListParser.ParseError:
            try:
                scan_not_found_error[raw_file_name].append(scan_id)
            except KeyError:
                scan_not_found_error[raw_file_name] = [scan_id]
            continue

        peak_list = peakListParser.get_peak_list(scan, peak_list_reader['fileType'])

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
            pep2 = str(id_item['pepseq 2'])
            if pep2 == 'nan':
                pep2 = ""
            linkpos1 = id_item['linkpos 1'] - 1
            linkpos2 = id_item['linkpos 2'] - 1
            cl_mod_mass = id_item['crosslinkermodmass']
        except KeyError:
            # linear
            pep2 = ""

        if pep2 == '':
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
        score = id_item['score']
        all_scores = json.dumps({'score': id_item['score']})

        # isDecoy
        try:
            is_decoy = 1 if id_item['isdecoy'] else 0
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
             all_scores,
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

    # multi error handler

    for pl_file, scan_id_list in scan_not_found_error.iteritems():
        return_json['errors'].append({
            "type": "",
            "message": "requested scanID(s) not found in peakList file %s" % pl_file,
            'id': ';'.join([str(scan_id) for scan_id in scan_id_list])
        })

    return return_json
