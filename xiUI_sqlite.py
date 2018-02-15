import sqlite3
import json

class DBException(Exception):
    pass


def connect(dbname):
    try:
        con = sqlite3.connect(dbname)
    except sqlite3.Error as e:
        raise DBException(e.message)

    return con


def create_tables(cur, con):
    try:
        cur.execute("DROP TABLE IF EXISTS uploads")
        cur.execute(
            "CREATE TABLE uploads("
            "id INT PRIMARY KEY, "
            "user_id INT,"
            "filename TEXT, "
            "raw_file_names TEXT, "
            "analysis_software JSON,"
            "provider JSON,"
            "audits JSON,"
            "samples JSON,"
            "analyses JSON,"
            "protocol JSON,"
            "bib JSON,"
            "upload_time DATE, "
            "upload_loc TEXT,"          
            "default_pdb TEXT)"
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
            "id text PRIMARY KEY, "
            "upload_id INT,"
            "accession VARCHAR(10), "
            "name TEXT, "
            "description TEXT, "
            "sequence TEXT, "
            "is_decoy BOOLEAN)"
        )
        cur.execute("DROP TABLE IF EXISTS peptides")
        cur.execute(
            "CREATE TABLE peptides("
            "id text PRIMARY KEY, "
            "upload_id INT,"
            "sequence TEXT,"
            "seq_mods TEXT,"
            "link_site int,"
            "crosslinker_modmass FLOAT)"

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
        cur.execute("DROP TABLE IF EXISTS spectrum")
        cur.execute(
            "CREATE TABLE spectrum("
            "id INT, "
            "upload_id INT,"
            "peak_list text, "
            "raw_file_name text, "
            "scan_id INT, "
            "frag_tol)"
        )
        cur.execute("DROP TABLE IF EXISTS spectrum_identifications")
        cur.execute(
            "CREATE TABLE spectrum_identifications("
            "id INT, "
            "upload_id INT,"
            "spectrum_id INT, "
            "pep1_id text, "
            "pep2_id text, "
            "pass_threshold text, "
            "rank int,"
            "ions JSON, "
            "scores JSON)"
        )
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)
    return True


def write_upload(inj_list, cur, con):
    try:
        cur.executemany("""
    INSERT INTO uploads (
        'user_id',
        'filename',
        'analysis_software',
        'provider',
        'audits',
        'samples',
        'analyses',
        'protocol',
        'bib',
        'upload_time',
        'upload_loc'
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_protocols(inj_list, cur, con):

    # try:
    #     cur.executemany("""
    # INSERT INTO db_sequences (
    #     'id',
    #     'accession',
    #     'name',
    #     'description',
    #     'sequence',
    #     'is_decoy'
    # )
    # VALUES (?, ?, ?, ?, ?, ?)""", inj_list)
    #     con.commit()
    #
    # except sqlite3.Error as e:
    #     raise DBException(e.message)

    return True

def write_db_sequences(inj_list, cur, con):

    try:
        cur.executemany("""
    INSERT INTO db_sequences (
        'id',
        'accession',
        'name',
        'description',
        'sequence',
        'is_decoy'
    )
    VALUES (?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_peptides(inj_list, cur, con):
    try:
        cur.executemany("""
    INSERT INTO peptides (
        'id',
        'sequence',
        'seq_mods',
        'link_site',
        'crosslinker_modmass'
    )
    VALUES (?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_modifications(inj_list, cur, con):
    try:
        cur.executemany("""INSERT INTO modifications ('id', 'name', 'mass', 'residues', 'accession') VALUES (?, ?, ?, ?, ?)""",
                        inj_list)
        con.commit()
    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_peptide_evidences(inj_list, cur, con):
    try:
        cur.executemany("""
    INSERT INTO peptide_evidences (
        'peptide_ref',
        'dbsequence_ref',
        'start',
        'is_decoy'
    )
    VALUES (?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_spectra(inj_list, cur, con):
    try:
        cur.executemany("""INSERT INTO spectrum ('id', 'upload_id', 'peak_list', 'raw_file_name', 'scan_id', 'frag_tol')
                        VALUES (?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_spectrum_identifications(inj_list, cur, con):
    try:
        cur.executemany("""INSERT INTO spectrum_identifications ('id', 'upload_id', 'spectrum_id', 'pep1_id', 'pep2_id',
                            rank, 'pass_threshold', 'ions', 'scores') VALUES (?, ?, ?, ?, ?, ?, ? , ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


# con = connect('/home/lars/Xi/xiSPEC_ms_parser/dbs/saved/Tmuris_exosomes1.db')
# cur = con.cursor()


def fill_in_missing_scores(cur, con):
    # try:
    #     cur.execute("""
    #   SELECT DISTINCT scoresJSON.key as scoreKey
    #   FROM identifications, json_each(identifications.allScores) AS scoresJSON""")
    #
    #     all_scores = cur.fetchall()
    #     all_scores = set([str(x[0]) for x in all_scores])
    #
    #     multiple_inj_list = []
    #
    #     cur.execute('SELECT id, allScores FROM identifications')
    #     res = cur.fetchall()
    #
    #     for row in res:
    #         row_scores = json.loads(row[1])
    #         missing = all_scores - set(row_scores.keys())
    #         missing_dict = {key: -1 for key in missing}
    #
    #         if len(missing) > 0:
    #             row_scores.update(missing_dict)
    #             multiple_inj_list.append([json.dumps(row_scores), row[0]])
    #             # cur.execute('UPDATE identifications SET allScores=? WHERE id = row[0]', json.dumps(row_scores))
    #
    #     cur.executemany("""
    #     UPDATE identifications
    #     SET `allScores` = ?
    #     WHERE `id` = ?""", multiple_inj_list)
    #
    #     con.commit()
    #
    # except sqlite3.Error as e:
    #     raise DBException(e.message)
    pass




