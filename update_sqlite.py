import SQLite


def update_database(con):

    cur = con.cursor()

    cur.execute(
        "CREATE TABLE IF NOT EXISTS meta_data("
        "upload_id INT,"
        "sid_meta1_name TEXT,"
        "sid_meta2_name TEXT,"
        "sid_meta3_name TEXT,"
        "contains_crosslink BOOLEAN)"
    )

    try:
        cur.execute('ALTER TABLE spectrum_identifications ADD COLUMN meta1 TEXT')
        cur.execute('ALTER TABLE spectrum_identifications ADD COLUMN meta2 TEXT')
        cur.execute('ALTER TABLE spectrum_identifications ADD COLUMN meta3 TEXT')
    except:
        pass  # columns already updated

    con.commit()

    return True


import glob

for db_name in glob.glob("./dbs/saved/*.db"):
    con = SQLite.connect(db_name)
    update_database(con)
