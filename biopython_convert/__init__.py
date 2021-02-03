"""
Biopython Convert
Convert between any formats that Biopython supports or gffutils.
Provides a means of querying/filtering documents using JMESPath query language.
"""
import sys
import pathlib
import itertools
from collections import defaultdict

import getopt

from Bio import SeqIO
import gffutils
from gffutils import biopython_integration

from . import JMESPathGen

gff_types = ['gff', 'gff3']
extended_types = ['txt', 'json', 'yaml', 'yml']
SeqIO_types = ['abi', 'abi-trim', 'ace', 'cif-atom', 'cif-seqres', 'clustal', 'embl', 'fasta', 'fasta-2line',
               'fastq-sanger', 'fastq', 'fastq-solexa', 'fastq-illumina', 'genbank', 'gb', 'ig', 'imgt', 'nexus',
               'pdb-seqres', 'pdb-atom', 'phd', 'phylip', 'pir', 'seqxml', 'sff', 'sff-trim', 'stockholm', 'swiss',
               'tab', 'qual', 'uniprot-xml']
stat_annotations = ['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi',
                    'keywords', 'source', 'organism']

usage = """\
Use: biopython.convert [-s] [-v] [-i] [-q JMESPath] input_file input_type output_file output_type
\t-s Split records into seperate files
\t-q JMESPath to select records. Must return list of SeqIO records. Root is list of input SeqIO records.
\t-i Print out details of records during conversion
\t-v Print version and exit
""" + "\nValid types: " + ', '.join(SeqIO_types + gff_types) + "\n"


def get_args(sysargs: list):
    """
    Parse command line arguments
    :param sysargs: list of command line arguments (sys.argv[1:])
    :return: (input_path, input_type, output_path, output_type, split, jmespath, stats)
    """
    split = False
    jpath = None
    stats = None
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sysargs, 'vsiq:')
        for opt, val in opts:
            if opt == '-v':
                from . import __version
                print(__version.__version__)
                exit(0)
            elif opt == '-s':
                split = True
            elif opt == '-q':
                jpath = val
            elif opt == '-i':
                stats = sys.stdout

    except getopt.GetoptError as err:
        print("Argument error(" + str(err.opt) + "): " + err.msg, file=sys.stderr)
        args = []

    # Check for minimum number of arguments
    if len(args) < 4:
        print(usage, file=sys.stderr)
        exit(1)

    input_path = pathlib.Path(args[0])
    input_type = args[1]
    output_path = pathlib.Path(args[2])
    output_type = args[3]

    return input_path, input_type, output_path, output_type, split, jpath, stats


def to_stats(record: SeqIO.SeqRecord) -> str:
    """
    Build GFF record representing summary of SeqRecord
    :param record: SeqIO.SeqRecord to represent
    :return: string containing GFF record
    """
    # 'if' statements can be removed after https://github.com/daler/gffutils/pull/144
    if record.name:
        attributes = {'Name': [record.name]}
    else:
        attributes = {}
    for k, v in record.annotations.items():
        if k in stat_annotations:
            if isinstance(v, list):
                v = [str(a) for a in v]
            else:
                v = [str(v)]
            if v:
                attributes[k] = v

    # Count features of each type
    feat_count = defaultdict(int)
    for f in record.features:
        feat_count[f.type] += 1
    attributes['features'] = [f"{k}:{v}" for k, v in feat_count.items()]

    if record.description:
        attributes['desc'] = [record.description]
    return str(gffutils.Feature(record.id, "biopython.convert", "sequence", start=1, end=len(record), attributes=attributes))


def get_records(input_handle, input_type: str, jpath: str = ''):
    """
    Read in records and apply optional jmespath
    :param input_handle: File handle to read data from
    :param input_type: one of abi,abi-trim,ace,cif-atom,cif-seqres,clustal,embl,fasta,fasta-2line,fastq-sanger,fastq,
        fastq-solexa,fastq-illumina,genbank,gb,ig,imgt,nexus,pdb-seqres,pdb-atom,phd,phylip,pir,seqxml,sff,sff-trim,
        stockholm,swiss,tab,qual,uniprot-xml,gff3
    :param jpath: JMESPath selecting records to keep. The root is the list of records. The path must return a list of records.
    :return: iterable of resulting records
    """
    def gentype(x):
        # Newer versions of biopython return iterator that is not of type Generator, wrap in generator type so jmespath can detect
        for a in x:
            yield a

    if input_type in gff_types:
        # If input is GFF load with gffutils library
        db = gffutils.create_db(input_handle, ":memory:", merge_strategy="create_unique")
        # Wrap features in generator that converts to BioPython SeqRecords
        input_records = map(
            lambda x: SeqIO.SeqRecord("", features=list(x)),
            itertools.groupby(
                map(biopython_integration.to_seqfeature, db.all_features(order_by="seqid")),
                lambda x: x.id
            )
        )
    else:
        input_records = SeqIO.parse(input_handle, input_type)

    # Wrap input in JMESPath selector if provided
    if jpath:
        input_records = JMESPathGen.search(jpath, gentype(input_records))

    if isinstance(input_records, dict):
        # Support generating a new record in JMESPath
        input_records = SeqIO.SeqRecord(**input_records)

    if isinstance(input_records, SeqIO.SeqRecord):
        # Support returning single record from JMESPath
        input_records = (input_records,)

    return input_records


def _generate_suffixes(path: pathlib.Path) -> pathlib.Path:
    """
    Helper to generate a new file path from a base path on each iteration
    :param path: base path to add suffix to
    :return: new path
    """
    i = 0
    while True:
        # Open output file with file name suffix if splitting
        yield path.with_suffix(f".${i}${path.suffix}")
        i += 1


def gff_writer(records: [SeqIO.SeqRecord], handle, output_type: str):
    """
    Convert SeqRecord to gffutils GFF3 record and output to handle
    :param handle: file handle to write to
    :param records: iterable of SeqRecord instances
    :param output_type: output format, ignored
    :return: None
    """
    for record in records:
        # TODO extend gffutils SeqFeature support
        for feature in record.features:
            feature = biopython_integration.from_seqfeature(feature)
            feature.seqid = record.id
            feature.source = 'biopython.convert'
            print(feature, file=handle)


def _to_SeqRecord(record):
    """
    Helper to convert all output records to SeqRecords
    :param record: dict or SeqRecord
    :return: SeqRecord
    """
    if isinstance(record, dict):
        # Support generating new records in JMESPath
        record = SeqIO.SeqRecord(**record)
    return record


def _print_stats(record, stats):
    """
    Helper to print stats of record
    :param record: SeqRecord to print stats of
    :param stats: IO handle or None
    :return: record, unaltered
    """
    if stats:
        print(to_stats(record), file=stats)
    return record


def convert(input_path, input_type, output_path, output_type, split=None, jpath='', stats=None):
    with input_path.open("r") as handle:
        if output_type == 'txt':
            raise NotImplemented
        elif output_type == 'json':
            raise NotImplemented
        elif output_type in ('yml', 'yaml'):
            raise NotImplemented
        elif output_type in gff_types:
            writer = gff_writer
        else:
            writer = SeqIO.write

        if stats:
            print("##gff-version 3", file=stats)

        seq_records = map(_to_SeqRecord, get_records(handle, input_type, jpath))
        if split:
            for record, path in zip(seq_records, _generate_suffixes(output_path)):
                _print_stats(record, stats)
                with path.open('w') as output_handle:
                    writer((record,), output_handle, output_type)
        else:
            with output_path.open('w') as output_handle:
                writer(
                    map(
                        lambda r: _print_stats(r, stats),
                        seq_records
                    ),
                    output_handle,
                    output_type
                )

