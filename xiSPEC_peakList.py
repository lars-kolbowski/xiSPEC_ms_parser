import ntpath
import zipfile
import glob


class ParseError(Exception):
    pass


def unzip_peak_lists(zip_file):
    zip_ref = zipfile.ZipFile(zip_file, 'r')
    unzip_path = zip_file + '_unzip/'
    zip_ref.extractall(unzip_path)
    zip_ref.close()

    return_file_list = []

    for file_path in glob.glob(unzip_path + '/*'):
        file_name = ntpath.basename(file_path)
        if file_name.lower().endswith('.mgf') or file_name.lower().endswith('.mzml'):
            return_file_list.append(file_path)
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
        peaks = zip(scan['m/z array'], scan['intensity array'])
        peak_list = "\n".join(["%s %s" % (mz, i) for mz, i in peaks if i > 0])

    else:
        raise ParseError("unsupported peak list file type: %s" % pl_file_type)

    return peak_list


def get_scan(readers, file_name, scan_id):
    try:
        reader = readers[file_name]
    except KeyError:
        if len(readers.keys()) == 1:
            reader = readers[readers.keys()[0]]
        else:
            raise ParseError("%s from identifications file does not match any of your peaklist files: %s" %
                             (file_name, ';'.join(readers.keys())))
    try:
        scan = reader[scan_id]
    except (IndexError, KeyError):
        raise ParseError("requested scanID %i not found in peakList file" % scan_id)

    return scan



