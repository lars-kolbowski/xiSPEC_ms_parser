# xiSPEC mgf reader - Lars Kolbowski
#
#    adopted from pymzml.run.reader
#   (Copyright (C) 2010-2014 T. Bald, J. Barth, A. Niehues, M. Specht, H. Roest, C. Fufezan)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import re
import os
import bisect
import codecs

from collections import defaultdict as ddict


class RegexPatterns(object):
    params_pattern = re.compile('([A-Z]+)=(.*)')
    peak_list_pattern = re.compile('(^(?:[0-9.]+\s[0-9.]+\s+)+)', re.M)


class ParseError(Exception):
    pass


class Reader(object):
    """

    Initializes an indexed mgf reader.

    :param path: path to mgf file.
    :type path: string

    :param file_object: file object or any other iterable stream, this will make
                        path obsolete, seeking is disabled
    :type file_object: File_object like

    Example:

    """

    def __init__(
            self,
            path=None,
            file_object=None,
    ):

        # self.info contains information extracted from the mgf file
        self.info = dict()

        self.info['offsetList'] = []

        # self.info['spectra_count'] = 0

        # self.info['encoding'] = None

        assert path is not None or file_object is not None, \
            'Must provide either a path or a file object to parse'

        self.info['fileObject'], self.info['seekable'] = self.__open_file(
            path,
            file_object
        )
        self.info['filename'] = path

        self.seeker = self._build_index()

        self.spectrum = {}

        return

    def _open_file(self, path, given_file_object=None):
        return self.__open_file(path, given_file_object=given_file_object)

    def __open_file(self, path, given_file_object=None):
        # Arbitrary supplied file objects are not seekable
        file_object = given_file_object
        seekable = False
        if file_object is None:
            import codecs
            if path.endswith('.gz'):
                # Gzipped files are not seekable
                import gzip
                file_object = codecs.getreader("utf-8")(
                    gzip.open(path)
                )
            else:
                file_object = codecs.open(
                    path,
                    mode='r'
                )
                seekable = True

        return file_object, seekable

    def _build_index(self):
        """
        .. method:: _build_index()

        Builds an index: a list of offsets to which a file pointer can seek
        directly to access a particular spectrum without parsing the entire file.

        :returns: A file-like object used to access the indexed content by
                  seeking to a particular offset for the file.
        """

        # Declare the seeker
        seeker = open(self.info['filename'], 'rb')

        self.info['offsets'] = None
        seeker.seek(0, 2) #  what's this for? - cc

        self._build_index_from_scratch(seeker)

        seeker.close()
        seeker = codecs.open(
            self.info['filename'],
            mode='rb'
        )

        return seeker

    def _build_index_from_scratch(self, seeker):
        """Build an index of spectra data with offsets by parsing the file."""

        def get_data_indices(fh):
            """Get a list with binary file indices of spectra in mgf file."""
            spec_positions = []

            # go to start of file
            fh.seek(0)
            pos = 0
            scan_start_pos = None
            for line in fh:
                if line[0] == "S":
                    if scan_start_pos:
                        scan_end_pos = pos
                        spec_positions.append((scan_start_pos, scan_end_pos))
                    scan_start_pos = pos
                pos = pos + len(line)

            # last one
            spec_positions.append((scan_start_pos, pos))

            return spec_positions



            # go to start of file
            # fh.seek(0)
            # pos = 0
            # peak_list_start_pos = None
            # for line in fh:
            #     if line[0].isdigit():
            #         if peak_list_start_pos is None or peak_list_start_pos == -1:
            #             peak_list_start_pos = pos
            #
            #     elif peak_list_start_pos is not None and peak_list_start_pos != -1:
            #             spec_positions.append((peak_list_start_pos, pos))
            #             peak_list_start_pos = -1
            #
            #     pos = pos + len(line)
            #
            # # last one
            # if peak_list_start_pos != -1:
            #     spec_positions.append((peak_list_start_pos, pos))
            # return spec_positions

        indices = get_data_indices(seeker)
        if indices is None:
            raise ParseError()
        self.info['offsetList'] = indices
        self.info['seekable'] = True

        return

    def get_by_id(self, scan_id):
        """"
         Random access to spectrum peak list in mgf by scanId

         """

        position = self.info['offsetList'][scan_id]
        start_pos = position[0]
        end_pos = position[1]

        if start_pos == -1:  # empty scan
            self.spectrum['peaks'] = ''
            # self.spectrum['params'] = params
            return self.spectrum

        self.seeker.seek(start_pos, 0)
        scan = self.seeker.read(end_pos - start_pos)

        if scan is None:
            raise KeyError("MS2 file does not contain a spectrum with index {0}.".format(scan_id))
        else:
            self.spectrum['peaks'] = self.parse_peak_list(scan)
            self.spectrum['precursor'] = self.parse_precursor(scan)
            return self.spectrum

    def __getitem__(self, scan_id):
        """"
        Random access to spectrum peak list in mgf by scanId

        """
        return self.get_by_id(scan_id)

    @staticmethod
    def parse_peak_list(raw_scan):
        lines = raw_scan.splitlines()
        peaks = []

        for line in lines:
            if re.match('[0-9\.]+\s[0-9\.]+', line):
                peaks.append(line)
            # if not line.startswith('#') and len(line.split('=')) == 1:
            #     peaks.append(line)

        return '\n'.join(peaks)

    @staticmethod
    def parse_precursor(raw_scan):
        lines = raw_scan.splitlines()
        precursor = {
            'mz': None,
            'charge': None
        }

        proton_mass = 1.007277

        for line in lines:
            if line.startswith('Z'):
                precursor_match = re.match('Z\s+([0-9]+)\s+([0-9\.]+)', line)
                if precursor_match:
                    precursor['charge'] = int(precursor_match.groups()[0])
                    precursor_mass = float(precursor_match.groups()[1])
                    precursor['mz'] = precursor_mass / precursor['charge'] + proton_mass
                else:
                    raise Exception("Error parsing precursor from scan")

        return precursor

