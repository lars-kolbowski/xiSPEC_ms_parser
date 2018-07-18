import re

def get_db_sequence_dict(fasta_file_list):
    db_sequence_dict = {}
    identifier = None
    sequence = ""
    description = None
    for fasta_file in fasta_file_list:
        with open(fasta_file) as f:
            for line in f:
                # semi-colons indicate comments, ignore them
                if not line.startswith(";"):
                    if line.startswith(">"):
                        if identifier is not None:
                            add_entry(identifier, sequence, description, db_sequence_dict)

                            #clear sequence
                            sequence = ""

                        # get new identifier
                        identifier = line
                        if " " not in line:
                            identifier = line[1:].rstrip()
                        else:
                            iFirstSpace = line.index(" ")
                            identifier = line[1:iFirstSpace].rstrip()
                            description = line[iFirstSpace:].rstrip()
                    else:
                        sequence += line.rstrip()

    # add last entry
    add_entry(identifier, sequence, description, db_sequence_dict)

    return db_sequence_dict


def add_entry(identifier, sequence, description, seq_dict):
    m = re.search("..\|(.*)\|(.*)\s?", identifier)
    # id = identifier
    accession = identifier
    name = identifier
    if m:
        accession = m.groups()[0]
        name = m.groups()[1]

    data = [accession, name, description, sequence]
    seq_dict[identifier] = data
    seq_dict[accession] = data
