import psycopg2
import json


class DBException(Exception):
    pass


def connect(dbname):
    import credentials
    try:
        con = psycopg2.connect(host=credentials.hostname, user=credentials.username, password=credentials.password,
                               dbname=credentials.database)
    except psycopg2.Error as e:
        raise DBException(e.message)

    return con


def create_tables(cur, con):
    try:
        cur.execute("DROP TABLE IF EXISTS uploads")
        cur.execute(
            "CREATE TABLE uploads("
            "id SERIAL PRIMARY KEY, "
            "user_id INT,"
            "filename TEXT, "
            "peak_list_file_names TEXT, "
            "analysis_software JSON,"
            "provider JSON,"
            "audits JSON,"
            "samples JSON,"
            "analyses JSON,"
            "protocol JSON,"
            "bib JSON,"
            "upload_time DATE, "
            "upload_loc TEXT,"          
            "default_pdb TEXT,"
            "contains_crosslink BOOLEAN,"
            "upload_errors JSON)"
        )
        cur.execute("DROP TABLE IF EXISTS protocols")
        cur.execute(
            "CREATE TABLE protocols("
            "id text PRIMARY KEY, "
            "upload_id INT,"
            "protocol JSON,"
            "ms2_tol FLOAT)"
        )
        cur.execute("DROP TABLE IF EXISTS db_sequences")
        cur.execute(
            "CREATE TABLE db_sequences("
            "id text, "
            "upload_id INT,"
            "accession TEXT, "
            "name TEXT, "
            "description TEXT, "
            "sequence TEXT, "
            "is_decoy BOOLEAN)"
        )
        cur.execute("DROP TABLE IF EXISTS peptides")
        cur.execute(
            "CREATE TABLE peptides("
            "id text, "
            "upload_id INT,"
            #"sequence TEXT,"
            "seq_mods TEXT,"
            "link_site INT,"
            "crosslinker_modmass FLOAT,"
            "value INT)"

        )
        cur.execute("DROP TABLE IF EXISTS modifications")
        cur.execute(
            "CREATE TABLE modifications("
            "id INT PRIMARY KEY, "
            "upload_id INT,"
            "name TEXT, "
            "mass FLOAT, "
            "residues TEXT, "
            "accession TEXT)"
        )
        cur.execute("DROP TABLE IF EXISTS peptide_evidences")
        cur.execute(
            "CREATE TABLE peptide_evidences("
            "upload_id INT,"
            "peptide_ref text, "
            "dbsequence_ref text, "
            "start int, "
            "is_decoy BOOLEAN)"
        )
        cur.execute("DROP TABLE IF EXISTS spectra")
        cur.execute(
            "CREATE TABLE spectra("
            "id INT, "
            "upload_id INT,"
            "peak_list text, "
            "peak_list_file_name text, "
            "scan_id TEXT, "
            "frag_tol TEXT,"
            "spectrum_id TEXT)"
        )
        cur.execute("DROP TABLE IF EXISTS spectrum_identifications")
        cur.execute(
            "CREATE TABLE spectrum_identifications("
            "id INT, "
            "upload_id INT,"
            "spectrum_id INT, "
            "pep1_id text, "
            "pep2_id text, "
            "charge_state int, "
            "pass_threshold text, "
            "rank int,"
            "ions JSON, "
            "scores JSON)"
        )
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    return True


def write_upload(inj_list, cur, con):
    try:
        cur.execute("""
    INSERT INTO uploads (
        user_id,
        filename,
        peak_list_file_names,
        analysis_software,
        provider,
        audits,
        samples,
        analyses,
        protocol,
        bib,
        upload_time,
        upload_loc
    )
    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, CURRENT_TIMESTAMP, %s) RETURNING id AS upload_id""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    rows = cur.fetchall()
    return rows[0]


def write_protocols(inj_list, cur, con):
    return True

def write_db_sequences(inj_list, cur, con):

    try:
        cur.executemany("""
        INSERT INTO db_sequences (
            id,
            accession,
            name,
            description,
            sequence,
            upload_id
        )
        VALUES (%s, %s, %s, %s, %s, %s) """, inj_list)
        #     con.commit()
        #
    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_peptides(inj_list, cur, con):
    # try:
    cur.executemany("""
    INSERT INTO peptides (
        id,
        /*sequence,*/
        seq_mods,
        link_site,
        crosslinker_modmass,
        upload_id,
        value
    )
    VALUES (%s, %s, %s, %s, %s, %s)""", inj_list)
    #     con.commit()
    #
    # except psycopg2.Error as e:
    #     raise DBException(e.message)

    return True


def write_modifications(inj_list, cur, con):
    # try:
    cur.executemany("""INSERT INTO modifications (id, name, mass, residues, accession) VALUES (%s, %s, %s, %s, %s)""",
                        inj_list)
    #     con.commit()
    # except psycopg2.Error as e:
    #     raise DBException(e.message)

    return True


def write_peptide_evidences(inj_list, cur, con):
    # try:
    cur.executemany("""
    INSERT INTO peptide_evidences (
        peptide_ref,
        dbsequence_ref,
        start,
        is_decoy,
        upload_id
    )
    VALUES (%s, %s, %s, %s, %s)""", inj_list)
    #     con.commit()
    #
    # except psycopg2.Error as e:
    #     raise DBException(e.message)

    return True


def write_spectra(inj_list, cur, con):
    # try:
    cur.executemany("""INSERT INTO spectra (id, peak_list, peak_list_file_name, scan_id, frag_tol, upload_id, spectrum_id)
                        VALUES (%s, %s, %s, %s, %s, %s, %s)""", inj_list)
    #     con.commit()
    #
    # except psycopg2.Error as e:
    #     raise DBException(e.message)

    return True


def write_spectrum_identifications(inj_list, cur, con):
    # try:
    cur.executemany("""INSERT INTO spectrum_identifications (id, upload_id, spectrum_id, pep1_id, pep2_id,
                            charge_state, rank, pass_threshold, ions, scores) VALUES (%s, %s, %s, %s, %s, %s, %s, %s , %s, %s)""", inj_list)
    #     con.commit()
    #
    # except psycopg2.Error as e:
    #     raise DBException(e.message)

    return True


# con = connect(/home/lars/Xi/xiSPEC_ms_parser/dbs/saved/Tmuris_exosomes1.db)
# cur = con.cursor()


def fill_in_missing_scores(cur, con):
    try:
        cur.execute("""
      SELECT DISTINCT scoresJSON.key as scoreKey
      FROM spectrum_identifications, json_each(spectrum_identifications.scores) AS scoresJSON""")

        all_scores = cur.fetchall()
        all_scores = set([str(x[0]) for x in all_scores])

        multiple_inj_list = []

        cur.execute('SELECT id, scores FROM spectrum_identifications')
        res = cur.fetchall()

        for row in res:
            row_scores = json.loads(row[1])
            missing = all_scores - set(row_scores.keys())
            missing_dict = {key: -1 for key in missing}

            if len(missing) > 0:
                row_scores.update(missing_dict)
                multiple_inj_list.append([json.dumps(row_scores), row[0]])
                # cur.execute('UPDATE identifications SET allScores=%s WHERE id = row[0]', json.dumps(row_scores))

        cur.executemany("""
        UPDATE spectrum_identifications
        SET `allScores` = %s
        WHERE `id` = %s""", multiple_inj_list)

        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    pass