from FullCsvParser import FullCsvParser
# import MinimalCsvParser
# import SuperMinimalCsvParser

import zipfile
import gzip
import os

# import unittest
#
#
# class MyTestCase(unittest.TestCase):
#     def test_something(self):
#         self.assertEqual(True, False)
#
#
# if __name__ == '__main__':
#     unittest.main()

def makeCsvParser(csv_path, temp_dir, peak_list_dir, db, logger, db_name='', user_id=0):
    if csv_path.endswith('.gz') or csv_path.endswith('.zip'):
        csv_path = extract_csv(csv_path)

    full_parser = FullCsvParser(csv_path, temp_dir, peak_list_dir, db, logger, db_name, user_id)
    if full_parser.check_required_columns():
        return full_parser

    # minimal_parser = MinimalCsvParser(csv_path, temp_dir, peak_list_dir, db, logger, db_name, user_id)
    # if minimal_parser.check_required_columns():
    #     return minimal_parser
    #
    # super_minimal_parser = SuperMinimalCsvParser(csv_path, temp_dir, peak_list_dir, db, logger, db_name, user_id)
    # if super_minimal_parser.check_required_columns():
    #     return super_minimal_parser

# ToDo: split into two functions?
@staticmethod
def extract_csv(archive):
    if archive.endswith('zip'):
        zip_ref = zipfile.ZipFile(archive, 'r')
        unzip_path = archive + '_unzip/'
        zip_ref.extractall(unzip_path)
        zip_ref.close()

        return_file_list = []

        for root, dir_names, file_names in os.walk(unzip_path):
            file_names = [f for f in file_names if not f[0] == '.']
            dir_names[:] = [d for d in dir_names if not d[0] == '.']
            for file_name in file_names:
                os.path.join(root, file_name)
                if file_name.lower().endswith('.mzid'):
                    return_file_list.append(root+'/'+file_name)
                else:
                    raise IOError('unsupported file type: %s' % file_name)

        if len(return_file_list) > 1:
            raise StandardError("more than one mzid file found!")

        return return_file_list[0]

    elif archive.endswith('gz'):
        in_f = gzip.open(archive, 'rb')
        archive = archive.replace(".gz", "")
        out_f = open(archive, 'wb')
        out_f.write(in_f.read())
        in_f.close()
        out_f.close()

        return archive

    else:
        raise StandardError('unsupported file type: %s' % archive)