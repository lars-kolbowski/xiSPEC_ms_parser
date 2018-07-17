import ntpath
import zipfile
import Ms2Reader as py_msn
import MGF as py_mgf
import pymzml
import re
import gzip


class PeakListParseError(Exception):
    pass


class ScanNotFoundException(Exception):
    pass


class PeakListParser:
    def __init__(self, pl_path, file_format_accession, spectrum_id_format_accession):
        # self.spectra_data = spectra_data
        self.file_format_accession = file_format_accession
        self.spectrum_id_format_accession = spectrum_id_format_accession
        self.peak_list_path = pl_path

        if self.is_mzML():
            self.reader = pymzml.run.Reader(pl_path)
        elif self.is_mgf():
            self.reader = py_mgf.Reader(pl_path)
        elif self.is_ms2():
            self.reader = py_msn.Reader(pl_path)
        else:
            raise PeakListParseError("unsupported peak list file type for: %s" % ntpath.basename(pl_path))

    def is_mgf(self):
        return self.file_format_accession == 'MS:1001062'

    def is_mzML(self):
        return self.file_format_accession == 'MS:1000584'

    def is_ms2(self):
        return self.file_format_accession == 'MS:1001466'


    @staticmethod
    def extract_gz(in_file):
        if in_file.endswith('.gz'):
            in_f = gzip.open(in_file, 'rb')
            in_file = in_file.replace(".gz", "")
            out_f = open(in_file, 'wb')
            out_f.write(in_f.read())
            in_f.close()
            out_f.close()

            return in_file

        else:
            raise StandardError("unsupported file extension for: %s" % in_file)

    @staticmethod
    def unzip_peak_lists(zip_file):
        """
        unzips and returns resulting folder
        :param zip_file: path to archive to unzip
        :return: resulting folder
        """

        if zip_file.endswith(".zip"):
            zip_ref = zipfile.ZipFile(zip_file, 'r')
            unzip_path = zip_file + '_unzip/'
            zip_ref.extractall(unzip_path)
            zip_ref.close()

            return unzip_path

        else:
            raise StandardError("unsupported file extension for: %s" % zip_file)

    @staticmethod
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

    def get_peak_list(self, scan_id):

        try:
            scan = self.reader[scan_id]
        except Exception as e:
            # raise ScanNotFoundException(type(e).__name__,
            #                             ntpath.basename(self.peak_list_path), e.args)
            raise ScanNotFoundException("%s - for file: %s - scanId: %s" % (e.args[0], ntpath.basename(self.peak_list_path), scan_id))

        if self.is_mzML():
            # if scan['ms level'] == 1:
            #     raise ParseError("requested scanID %i is not a MSn scan" % scan['id'])

            peak_list = "\n".join(["%s %s" % (mz, i) for mz, i in scan.peaks if i > 0])

        elif self.is_mgf():
            peak_list = scan['peaks']
            # peak_list = "\n".join(["%s %s" % (mz, i) for mz, i in scan['peaks'] if i > 0])

        elif self.is_ms2():
            peak_list = scan['peaks']

        else:   # this should never happen is it would have raise error in constructor
            raise PeakListParseError("unsupported peak list file type")

        return peak_list

    def parse_scan_id(self, spec_id):

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

        # ignore_dict_index = False
        identified_spec_id_format = False

        # if spec_id_format is not None and 'accession' in spec_id_format: # not needed, checked in constructor

        # MS:1000774 multiple peak list nativeID format - zero based
        if self.spectrum_id_format_accession == 'MS:1000774':
            identified_spec_id_format = True
            # ignore_dict_index = True
            try:
                matches = re.match("index=([0-9]+)", spec_id).groups()
                spec_id = int(matches[0])

            # try to cast spec_id to int if re doesn't match -> PXD006767 has this format
            # ToDo: do we want to be stricter?
            except (AttributeError, IndexError):
                try:
                    spec_id = int(spec_id)
                except ValueError:
                    raise PeakListParseError("invalid spectrum ID format!")

        # MS:1000775 single peak list nativeID format
        # The nativeID must be the same as the source file ID.
        # Used for referencing peak list files with one spectrum per file,
        # typically in a folder of PKL or DTAs, where each sourceFileRef is different.
        elif self.spectrum_id_format_accession == 'MS:1000775':
            identified_spec_id_format = True
            # ignore_dict_index = True
            spec_id = 0

        # MS:1000776 scan number only nativeID format
        # Used for referencing mzXML, or a DTA folder where native scan numbers can be derived.
        elif self.spectrum_id_format_accession == 'MS:1000776':
            identified_spec_id_format = True
            try:
                matches = re.match("scan=([0-9]+)", spec_id).groups()
                spec_id = int(matches[0])
            except (IndexError, AttributeError):
                raise PeakListParseError("invalid spectrum ID format!")

        # MS:1000768 Thermo nativeID format:
        # controllerType=xsd:nonNegativeIntege controllerNumber=xsd:positiveInteger scan=xsd:positiveInteger
        elif self.spectrum_id_format_accession == 'MS:1000768':
            identified_spec_id_format = True
            try:
                matches = re.search("scan=([0-9]+)", spec_id).groups()
                spec_id = int(matches[0])
            except (IndexError, AttributeError):
                raise PeakListParseError("invalid spectrum ID format!")

        # MS:1001530 mzML unique identifier:
        # Used for referencing mzML. The value of the spectrum ID attribute is referenced directly.
        elif self.spectrum_id_format_accession == 'MS:1001530':
            # ToDo: pymzML uses scan number as index and not the spectrumID attribute therefore
            # ToDo: we resort to using the default way of parsing the spec_id, i.e. using the last number in the string
            # ToDo: This might change if we are going to use a different (more standards compliant approach)
            identified_spec_id_format = False
            # matches = re.search("scan=([0-9]+)", spec_id).groups()
            # try:
            #     spec_id = int(matches[0])
            #     identified_spec_id_format = True
            # except IndexError:
            #     pass

        if not identified_spec_id_format:
            # ToDo: display warning or throw error? depending on strict mode or not?
            matches = re.findall("([0-9]+)", spec_id)
            # match = re.match("(.*)([0-9]+)", spec_id)
            try:
                spec_id = int(matches[-1])
                # spec_id = match.group(2)
            except IndexError:
                raise PeakListParseError("failed to parse spectrumID from %s" % spec_id)

            #
            # spec_id = match.group(2)

        return spec_id

