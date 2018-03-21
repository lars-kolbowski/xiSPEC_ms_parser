class DBException(Exception):
    pass


def connect(dbname):
    pass


def create_tables(cur, con):
    pass


def write_upload(inj_list, cur, con):
    pass


def write_db_sequences(inj_list, cur, con):
    pass


def write_peptides(inj_list, cur, con):
    pass


def write_modifications(inj_list, cur, con):
    pass


def write_peptide_evidences(inj_list, cur, con):
    pass


def write_spectra(inj_list, cur, con):
    pass


def write_spectrum_identifications(inj_list, cur, con):
    pass


def fill_in_missing_scores(cur, con):
    pass
