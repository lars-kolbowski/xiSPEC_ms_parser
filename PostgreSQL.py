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
    # don't create tables here
    # use file postgreSQL_schema.sql to init db
    #
    # you will need to search and replace 'username' in the sql file,
    # replacing it with the role name you use to access the database
    return True


def write_upload(inj_list, cur, con):
    try:
        cur.execute("""
    INSERT INTO uploads (
        user_id,
        filename,
        peak_list_file_names,
        spectra_formats,
        analysis_software,
        provider,
        audits,
        samples,
        analyses,
        protocol,
        bib,
        upload_time,
        origin,
        upload_warnings
    )
    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, CURRENT_TIMESTAMP, %s, %s) RETURNING id AS upload_id""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    rows = cur.fetchall()
    return rows[0][0]


def get_random_id(upload_id, cur, con):
    try:
        cur.execute("SELECT random_id FROM uploads WHERE id = " + str(upload_id) + ";")
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    rows = cur.fetchall()
    return rows[0][0]


# def write_protocols(inj_list, cur, con):
#     return True

def write_db_sequences(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO db_sequences (
            id,
            accession,
            protein_name,
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

def write_meta_data(values, cur, con):
    pass
    # try:
    #     cur.execute("""
    #       INSERT INTO meta_data (
    #         'upload_id',
    #         'sid_meta1_name',
    #         'sid_meta2_name',
    #         'sid_meta3_name',
    #         'contains_crosslink'
    #       )
    #       VALUES (?, ?, ?, ?, ?)""",  values)
    #     con.commit()
    #
    # except sqlite3.Error as e:
    #     raise DBException(e.message)
    #
    # return True

def write_peptides(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO peptides (
            id,
            seq_mods,
            link_site,
            crosslinker_modmass,
            upload_id,
            crosslinker_pair_id
        )
        VALUES (%s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_modifications(inj_list, cur, con):
    try:
        cur.executemany("""
          INSERT INTO modifications (
            id,
            upload_id,
            mod_name,
            mass,
            residues,
            accession
          )
          VALUES (%s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()
    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_peptide_evidences(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO peptide_evidences (
            peptide_ref,
            dbsequence_ref,
            protein_accession,
            pep_start,
            is_decoy,
            upload_id
        )
        VALUES (%s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_spectra(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO spectra (
        id, 
        peak_list, 
        peak_list_file_name, 
        scan_id, 
        frag_tol, 
        upload_id, 
        spectrum_ref,
        precursor_mz,
        precursor_charge
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_spectrum_identifications(inj_list, cur, con):
    try:
        cur.executemany("""
          INSERT INTO spectrum_identifications (
              id,
              upload_id,
              spectrum_id,
              pep1_id,
              pep2_id,
              charge_state,
              rank,
              pass_threshold,
              ions,
              scores,
              exp_mz,
              calc_mz,
              meta1,
              meta2,
              meta3
          ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s , %s, %s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True
