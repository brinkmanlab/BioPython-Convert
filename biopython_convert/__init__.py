"""
Biopython Convert
Convert between any formats that Biopython supports or gffutils.
Provides a means of querying/filtering documents using JMESPath query language.
"""
import sys
import pathlib
import itertools
import types
from collections import defaultdict

import getopt
from typing import Callable, Generator, OrderedDict

from Bio import SeqIO, StreamModeError, SeqFeature
import gffutils
from gffutils import biopython_integration

from . import JMESPathGen

gff_types = ['gff', 'gff3']
extended_types = ['text', 'json', 'yaml', 'yml']
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
""" + "\nValid types: " + ', '.join(SeqIO_types + gff_types + extended_types) + "\n"


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

    feat_count = defaultdict(int)
    for f in record.features:
        if f.type == 'source':
            # Include source coordinate
            attributes[f"source__location"] = [str(part) for part in f.location.parts]
            # Include source qualifiers
            for k, v in f.qualifiers.items():
                attr = attributes.get(f"source_{k}", [])
                if not isinstance(attr, list):
                    attr = [attr]
                attr.append(v)
                attributes[f"source_{k}"] = attr
        # Count features of each type
        feat_count[f.type] += 1
    attributes['features'] = [f"{k}:{v}" for k, v in feat_count.items()]

    if record.description:
        attributes['desc'] = [record.description]
    return str(gffutils.Feature(record.id, "biopython.convert", "sequence", start=1, end=len(record), attributes=attributes))


def _allow_single(records):
    """
    Helper to allow returing single record from JMESPath
    :param records: SeqIO.SeqRecord instance
    :return: tuple containing SeqIO.SeqRecord
    """
    if isinstance(records, SeqIO.SeqRecord):
        # Support returning single record from JMESPath
        return (records,)
    return records


def _to_SeqRecord(records):
    """
    Helper to convert all output records to SeqRecords
    :param records: dict or SeqRecord
    :return: SeqRecord
    """
    if isinstance(records, dict):
        # Support generating a single new record in JMESPath
        records = SeqIO.SeqRecord(**records)

    records = _allow_single(records)

    return map(lambda r: SeqIO.SeqRecord(**r) if isinstance(records, dict) else r, records)


def get_records(input_handle, input_type: str, jpath: str = '', xform: Callable = _to_SeqRecord):
    """
    Read in records and apply optional jmespath
    :param input_handle: File handle to read data from
    :param input_type: one of abi,abi-trim,ace,cif-atom,cif-seqres,clustal,embl,fasta,fasta-2line,fastq-sanger,fastq,
        fastq-solexa,fastq-illumina,genbank,gb,ig,imgt,nexus,pdb-seqres,pdb-atom,phd,phylip,pir,seqxml,sff,sff-trim,
        stockholm,swiss,tab,qual,uniprot-xml,gff3
    :param jpath: JMESPath selecting records to keep. The root is the list of records. The path must return a list of records.
    :param xform: Callable that takes the result of the jmespath and does anything necessary to convert to a iterable of output records
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

    # Apply xform to both entire return value
    input_records = xform(input_records)

    return input_records


def _generate_suffixes(path: pathlib.Path) -> Generator[pathlib.Path, None, None]:
    """
    Helper to generate a new file path from a base path on each iteration
    :param path: base path to add suffix to
    :return: new path
    """
    i = 0
    while True:
        # Open output file with file name suffix if splitting
        yield path.with_suffix(f".{i}{path.suffix}")
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


def _print_stats(record, stats):
    """
    Helper to print stats of record
    :param record: SeqRecord to print stats of
    :param stats: IO handle or None
    :return: record, unaltered
    """
    if stats and isinstance(record, SeqIO.SeqRecord):
        print(to_stats(record), file=stats)
    return record


def to_strings(v):
    """
    Helper to recursively convert Generators to lists, stringifing all else
    :param v: Parent object/list
    :return: list/dict with all children converted to the same or a string
    """
    if isinstance(v, str):
        return v

    if isinstance(v, (types.GeneratorType, map, filter, tuple)):
        v = list(v)

    if hasattr(v, 'keys'):
        keys = v.keys()
    elif hasattr(v, '__getitem__'):
        keys = range(len(v))
    else:
        return str(v)

    for i in keys:
        v[i] = to_strings(v[i])
    return v


def to_dicts(v):
    """
    Helper to recursively convert Objects and Generators to dicts and lists
    :param v: Parent object/list
    :return: list/dict with all children converted to the same
    """
    if isinstance(v, str):
        try:
            return int(v)
        except ValueError:
            pass
        return v

    if isinstance(v, (types.GeneratorType, map, filter, tuple)):
        v = list(v)

    if isinstance(v, SeqIO.SeqRecord):
        v = {
            **v.__dict__,
            'seq': str(v.seq)
        }
        del v['_seq']
    elif isinstance(v, SeqFeature.FeatureLocation):
        v = {
            **v.__dict__,
            'start': v.start,
            'end': v.end,
            'strand': v.strand,
        }
        del v['_start']
        del v['_end']
        del v['_strand']
    elif isinstance(v, SeqFeature.AbstractPosition):
        return to_dicts(str(v))
    elif isinstance(v, OrderedDict):
        v = dict(v)
    elif hasattr(v, '__dict__'):
        v = v.__dict__

    if hasattr(v, 'keys'):
        keys = v.keys()
    elif hasattr(v, '__getitem__'):
        keys = range(len(v))
    else:
        return v

    for i in keys:
        v[i] = to_dicts(v[i])
    return v


def convert(input_path: pathlib.Path, input_type: str, output_path: pathlib.Path, output_type: str, split: bool = False, jpath: str = '', stats=None):
    """
    Convert document from one format to another, optionally querying via JMESPath or splitting into separate outputs
    :param input_path: Path to input dataset
    :param input_type: Format of input dataset
    :param output_path: Path to output dataset
    :param output_type: Format of output dataset
    :param split: Split each record into a different output dataset. Adds index suffix to output path.
    :param jpath: JMESPath query to apply to input dataset before outputting
    :param stats: File handle to output GFF3 summary of output records
    :return: None
    """
    xform = _to_SeqRecord
    with input_path.open("r") as handle:
        if output_type == 'text':
            writer = lambda records, fh, t: fh.write("\n".join(map(str, to_strings(records))) + "\n")
            xform = lambda x: x
        elif output_type == 'json':
            import json
            writer = lambda records, fh, t: json.dump(to_dicts(records), fh, skipkeys=True, indent=True)
            xform = _allow_single
        elif output_type in ('yml', 'yaml'):
            from ruamel.yaml import YAML
            yml = YAML(typ='unsafe')
            writer = lambda records, fh, t: yml.dump(to_dicts(records), fh)
            xform = _allow_single
        elif output_type in gff_types:
            writer = gff_writer
        else:
            writer = SeqIO.write

        if stats:
            print("##gff-version 3", file=stats)

        seq_records = get_records(handle, input_type, jpath, xform)
        binary = ''
        if split:
            for record, path in zip(seq_records, _generate_suffixes(output_path)):
                _print_stats(record, stats)
                while True:
                    try:
                        with path.open('w' + binary) as output_handle:
                            writer((record,), output_handle, output_type)
                    except StreamModeError:
                        if binary == 'b':
                            raise
                        binary = 'b'
                        continue
                    break
        else:
            while True:
                try:
                    with output_path.open('w' + binary) as output_handle:
                        writer(
                            map(
                                lambda r: _print_stats(r, stats),
                                seq_records
                            ),
                            output_handle,
                            output_type
                        )
                except StreamModeError:
                    if binary == 'b':
                        raise
                    binary = 'b'
                    continue
                break

