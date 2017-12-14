import ntpath
import zipfile
import glob


def unzip_peak_lists(zip_file):
    zip_ref = zipfile.ZipFile(zip_file, 'r')
    unzip_path = '/'.join(zip_file.split('/')[:-1]) + '/unzip'
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

