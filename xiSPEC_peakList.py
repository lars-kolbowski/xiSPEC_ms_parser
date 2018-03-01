import ntpath
import zipfile
import xiSPEC_mgfReader as py_mgf
import pymzml
import re
import os
import gzip


class ParseError(Exception):
    pass


def unzip_peak_lists(zip_file):

    if zip_file.endswith(".zip"):
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

    elif zip_file.endswith('.gz'):
        in_f = gzip.open(zip_file, 'rb')
        zip_file = zip_file.replace(".gz", "")
        out_f = open(zip_file, 'wb')
        out_f.write(in_f.read())
        in_f.close()
        out_f.close()

        return [zip_file]
    else:
        raise StandardError("unsupported file extension for: %s" % zip_file)


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


def get_scan(reader, spec_id, file_id_format=None):


    # file_id_format_accession = file_id_format['accession']
    # #
    # # if (fileIdFormat == Constants.SpecIdFormat.MASCOT_QUERY_NUM) {
    # #     String rValueStr = spectrumID.replaceAll("query=", "");
    # #     String id = null;
    # #     if(rValueStr.matches(Constants.INTEGER)){
    # #         id = Integer.toString(Integer.parseInt(rValueStr) + 1);
    # #     }
    # #     return id;
    # # } else if (fileIdFormat == Constants.SpecIdFormat.MULTI_PEAK_LIST_NATIVE_ID) {
    # #     String rValueStr = spectrumID.replaceAll("index=", "");
    # #     String id;
    # #     if(rValueStr.matches(Constants.INTEGER)){
    # #         id = Integer.toString(Integer.parseInt(rValueStr) + 1);
    # #         return id;
    # #     }
    # #     return spectrumID;
    # # } else if (fileIdFormat == Constants.SpecIdFormat.SINGLE_PEAK_LIST_NATIVE_ID) {
    # #     return spectrumID.replaceAll("file=", "");
    # # } else if (fileIdFormat == Constants.SpecIdFormat.MZML_ID) {
    # #     return spectrumID.replaceAll("mzMLid=", "");
    # # } else if (fileIdFormat == Constants.SpecIdFormat.SCAN_NUMBER_NATIVE_ID) {
    # #     return spectrumID.replaceAll("scan=", "");
    # # } else {
    # #     return spectrumID;
    # # }
    #
    # if file_id_format == :
    #     rValueStr = scan_id.replaceAll("query=", "")
    #     id = null
    #     if(rValueStr.matches(Constants.INTEGER))
    #         id = Integer.toString(Integer.parseInt(rValueStr) + 1)
    #     return id
    # elif file_id_format == Constants.SpecIdFormat.MULTI_PEAK_LIST_NATIVE_ID:
    #     rValueStr = scan_id.replaceAll("index=", "")
    # if rValueStr.matches(Constants.INTEGER):
    #         id = Integer.toString(Integer.parseInt(rValueStr) + 1)
    #         return id
    #     return scan_id
    # elif file_id_format == Constants.SpecIdFormat.SINGLE_PEAK_LIST_NATIVE_ID:
    #     return spectrumID.replaceAll("file=", "")
    # elif file_id_format == Constants.SpecIdFormat.MZML_ID:
    #     return spectrumID.replaceAll("mzMLid=", "")
    # elif file_id_format == Constants.SpecIdFormat.SCAN_NUMBER_NATIVE_ID:
    #     return spectrumID.replaceAll("scan=", "")
    # else:
    #     return spectrumID;
    #
    #
    # # e.g.: MS:1000768(Thermo        nativeID        format)
    # # e.g.: MS:1000769(Waters        nativeID        format)
    # # e.g.: MS:1000770(WIFF        nativeID        format)
    # # e.g.: MS:1000771(Bruker / Agilent        YEP        nativeID        format)
    # # e.g.: MS:1000772(Bruker        BAF        nativeID        format)
    # # e.g.: MS:1000773(Bruker        FID        nativeID        format)
    # # e.g.: MS:1000774(multiple       peak        list        nativeID        format)
    # # e.g.: MS:1000775(single        peak        list        nativeID        format)
    # # e.g.: MS:1000776(scan        number        only        nativeID        format)
    # # e.g.: MS:1000777(spectrum        identifier        nativeID        format)

    if file_id_format is not None:

        if file_id_format['accession'] == 'MS:1000774':  # (multiple peak list nativeID format - zero based)
            matches = re.findall("index=([0-9]+)", spec_id)
            spec_id = int(matches[0])

        # MS:1000775
        # The nativeID must be the same as the source file ID.
        # Used for referencing peak list files with one spectrum per file,
        # typically in a folder of PKL or DTAs, where each sourceFileRef is different.
        elif file_id_format['accession'] == 'MS:1000775':
            spec_id = 0

        # MS:1000776
        # Used for referencing mzXML, or a DTA folder where native scan numbers can be derived.
        elif file_id_format['accession'] == 'MS:1000776':
            matches = re.findall("scan=([0-9]+)", spec_id)
            spec_id = int(matches[0])

        else:
            matches = re.findall("([0-9]+)", spec_id)

            if len(matches) == 1:
                spec_id = int(matches[0])
            else:
                raise ParseError("failed to parse spectrumID from %s" % spec_id)

        if reader['fileType'] == 'mgf':
            try:
                return reader['reader'].get_by_id(spec_id, ignore_dict_index=True)
            except (IndexError, KeyError):
                raise ParseError("requested scanID %s not found in peakList file" % spec_id)

    try:
        return reader['reader'][spec_id]
    except (IndexError, KeyError):
        raise ParseError("requested scanID %s not found in peakList file" % spec_id)


    # try:
    #     if reader['fileType'] == 'mgf':
    #         scan = reader['reader'][scan_id+1]  # ToDo: clear up 0/1 confusion
    #     elif reader['fileType'] == 'mzml':
    #         scan = reader['reader'][scan_id]
    # except (IndexError, KeyError):
    #     raise ParseError("requested scanID %i not found in peakList file" % scan_id)
    #
    # return scan


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


def add_peak_list_reader(pl_file, return_dict, key):

    pl_file_name = ntpath.basename(pl_file)

    if pl_file_name.lower().endswith('.mzml'):
        peak_list_file_type = 'mzml'
        reader = pymzml.run.Reader(pl_file)

    elif pl_file_name.lower().endswith('.mgf'):
        peak_list_file_type = 'mgf'
        reader = py_mgf.Reader(pl_file)

    else:
        raise ParseError("unsupported peak list file type for: %s" % pl_file_name)

    peak_list_reader_index = key  # re.sub("\.(mzml|mgf)\Z", "", pl_file_name, flags=re.I)
    return_dict[peak_list_reader_index] = {
        'reader': reader,
        'fileType': peak_list_file_type
    }

    return return_dict


# def create_peak_list_reader(pl_file_name):
#     pass
