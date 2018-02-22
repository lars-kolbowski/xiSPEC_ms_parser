class DBException(Exception):
    pass

def connect(dbname):
    return Bogus_db();


def create_tables(cur, con):
    return True


def write_identifications(inj_list, cur, con):
    return True


def write_peaklists(inj_list, cur, con):
    return True


def write_modifications(inj_list, cur, con):
    return True

def fill_in_missing_scores(cur, con):
    pass

class Bogus_db(object):
    # The class "constructor" - It's actually an initializer
    def __init__(self):
        pass

    def cursor(self):
        pass

    def commit(self):
        pass

    def close(self):
        pass