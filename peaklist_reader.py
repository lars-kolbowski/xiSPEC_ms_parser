

# peakList file
logger.info('reading peakList file - start')
peak_list_file_name = ntpath.basename(peak_list_file).lower()
if peak_list_file_name.endswith('.mzml'):
    peak_list_file_type = 'mzml'
    pymzml_reader = pymzml.run.Reader(peak_list_file)

elif peak_list_file_name.endswith('.mgf'):
    peak_list_file_type = 'mgf'
    mgf_reader = py_mgf.read(peak_list_file)
    peak_list_arr = [pl for pl in mgf_reader]
logger.info('reading peakList file - done')