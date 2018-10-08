from FullCsvParser import FullCsvParser
# from AbstractCsvParser import CsvParseException

class MimimalCsvParser(FullCsvParser):
    required_cols = [
        'pepseq1',
        'peppos1',
        'linkpos1',
        'protein1',
        'pepseq2',
        'peppos2',
        'linkpos2',
        'protein2',
        # 'expMZ'
    ]

    optional_cols = [
        'scanid',
        'charge',
        'peaklistfilename',
        'rank',
        'fragmenttolerance',
        'iontypes',
        'crosslinkermodmass',
        'passthreshold',
        'score',
        'decoy1',
        'decoy2',
        'expmz',  # ToDo: required in mzid - also make required col?
        'calcmz'
    ]
