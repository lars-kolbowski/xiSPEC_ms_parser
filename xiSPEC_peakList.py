import ntpath
import zipfile
import xiSPEC_mgfReader as py_mgf
import pymzml
import re
import os


class ParseError(Exception):
    pass


def unzip_peak_lists(zip_file):
    zip_ref = zipfile.ZipFile(zip_file, 'r')
    unzip_path = zip_file + '_unzip/'
    zip_ref.extractall(unzip_path)
    zip_ref.close()

    return_file_list = []

    for root, dir_names, file_names in os.walk(unzip_path):
        file_names = [f for f in file_names if not f[0] == '.']
        dir_names[:] = [d for d in dir_names if not d[0] == '.']
        for file_name in file_names:
            os.path.join(root, file_name)
            if file_name.lower().endswith('.mgf') or file_name.lower().endswith('.mzml'):
                return_file_list.append(root+'/'+file_name)
            else:
                raise IOError('unsupported file type: %s' % file_name)

    return return_file_list


def get_ion_types_mzml(scan):
    frag_methods = {
        'beam-type collision-induced dissociation': ["b", "y"],
        'collision-induced dissociation': ["b", "y"],
        'electron transfer dissociation': ["c", "z"],
    }
    # get fragMethod and translate that to Ion Types
    ion_types = []
    for key in scan.keys():
        if key in frag_methods.keys():
            ion_types += frag_methods[key]
    return ion_types


def get_peak_list(scan, pl_file_type):

    if pl_file_type == 'mzml':
        if scan['ms level'] == 1:
            raise ParseError("requested scanID %i is not a MSn scan" % scan['id'])

        peak_list = "\n".join(["%s %s" % (mz, i) for mz, i in scan.peaks if i > 0])

    elif pl_file_type == 'mgf':
        peak_list = scan['peaks']
        # peak_list = "\n".join(["%s %s" % (mz, i) for mz, i in scan['peaks'] if i > 0])

    else:
        raise ParseError("unsupported peak list file type: %s" % pl_file_type)

    return peak_list


def get_reader(readers, file_name):
    try:
        reader = readers[file_name]
    except KeyError:
        #add warning?
        if len(readers.keys()) == 1:
            reader = readers[readers.keys()[0]]
        else:
            raise ParseError("%s from identifications file does not match any of your peaklist files: %s" %
                             (file_name, '; '.join(readers.keys())))
    return reader


def get_scan(reader, scan_id):
    try:
        if reader['fileType'] == 'mgf':
            scan = reader['reader'][scan_id+1]
        elif reader['fileType'] == 'mzml':
            scan = reader['reader'][scan_id]
    except (IndexError, KeyError):
        raise ParseError("requested scanID %i not found in peakList file" % scan_id)

    return scan


def create_peak_list_readers(pl_file_list):
    return_dict = {}

    for pl_file in pl_file_list:

        pl_file_name = ntpath.basename(pl_file)

        if pl_file_name.lower().endswith('.mzml'):
            peak_list_file_type = 'mzml'
            reader = pymzml.run.Reader(pl_file)

        elif pl_file_name.lower().endswith('.mgf'):
            peak_list_file_type = 'mgf'
            reader = py_mgf.Reader(pl_file)

        else:
            raise ParseError("unsupported peak list file type for: %s" % pl_file_name)

        peak_list_reader_index = re.sub("\.(mzml|mgf)\Z", "", pl_file_name, flags=re.I)
        return_dict[peak_list_reader_index] = {
            'reader': reader,
            'fileType': peak_list_file_type
        }

    return return_dict

def create_peak_list_reader(pl_file_name):
    pass
