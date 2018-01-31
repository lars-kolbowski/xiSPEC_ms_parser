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
            "allScores TEXT,"
            "isDecoy INT, "
            "protein1 TEXT, "
            "protein2 TEXT, "
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
            "residues TEXT, "
            "accession TEXT)"
        )
        cur.execute("DROP TABLE IF EXISTS peakLists")
        cur.execute(
            "CREATE TABLE peakLists("
            "id INT PRIMARY KEY, "
            "peakList TEXT)"
        )

        cur.execute("DROP TABLE IF EXISTS db_sequences")
        cur.execute(
            "CREATE TABLE db_sequences("
            "id text PRIMARY KEY, "
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
            "sequence TEXT)"
        )
        cur.execute("DROP TABLE IF EXISTS peptide_evidences")
        cur.execute(
            "CREATE TABLE peptide_evidences("
            "peptide_ref text, "
            "dbsequence_ref text, "
            "start int, "
            "is_decoy BOOLEAN)"
        )

        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)
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
    INSERT INTO db_sequences (
        'id',
        'sequence'
    )
    VALUES (?, ?)""", inj_list)
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

def write_identifications(inj_list, cur, con):

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
        'allScores',
        'isDecoy',
        'protein1',
        'protein2',
        'file',
        'scanID',
        'peakList_id'
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_peaklists(inj_list, cur, con):
    try:
        cur.executemany("""INSERT INTO peakLists ('id', 'peaklist') VALUES (?, ?)""",
                        inj_list)
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

# con = connect('/home/lars/Xi/xiSPEC_ms_parser/dbs/saved/Tmuris_exosomes1.db')
# cur = con.cursor()


def fill_in_missing_scores(cur, con):
    try:
        cur.execute("""
      SELECT DISTINCT scoresJSON.key as scoreKey 
      FROM identifications, json_each(identifications.allScores) AS scoresJSON""")

        all_scores = cur.fetchall()
        all_scores = set([str(x[0]) for x in all_scores])

        multiple_inj_list = []

        cur.execute('SELECT id, allScores FROM identifications')
        res = cur.fetchall()

        for row in res:
            row_scores = json.loads(row[1])
            missing = all_scores - set(row_scores.keys())
            missing_dict = {key: -1 for key in missing}

            if len(missing) > 0:
                row_scores.update(missing_dict)
                multiple_inj_list.append([json.dumps(row_scores), row[0]])
                # cur.execute('UPDATE identifications SET allScores=? WHERE id = row[0]', json.dumps(row_scores))

        cur.executemany("""
        UPDATE identifications 
        SET `allScores` = ?
        WHERE `id` = ?""", multiple_inj_list)

        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)





