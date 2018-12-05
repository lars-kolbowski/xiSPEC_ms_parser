from csv_parser.FullCsvParser import FullCsvParser
# from AbstractCsvParser import CsvParseException

class NoPeakListsCsvParser(FullCsvParser):
    required_cols = [
        'pepseq1',
        'peppos1',
        'linkpos1',
        'protein1',
        'pepseq2',
        'peppos2',
        'linkpos2',
        'protein2',
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
        'expmz',
        'calcmz'
    ]
