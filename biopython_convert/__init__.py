#!/usr/bin/env python
"""
Biopython Convert
Convert between any formats that Biopython supports or gffutils.
Provides a means of querying/filtering documents using JMESPath query language.
"""
from Bio import SeqIO
from .__version import __version__, __versionstr__
import itertools
import getopt
import gffutils
from gffutils import biopython_integration
from . import JMESPathGen
import sys

usage = "Use: biopython.convert [-s] [-v] [-i] [-q JMESPath] input_file input_type output_file output_type\n" \
        "\t-s Split records into seperate files\n" \
        "\t-q JMESPath to select records. Must return list of SeqIO records. Root is list of input SeqIO records.\n" \
        "\t-i Print out details of records during conversion\n" \
        "\t-v Print version and exit\n"

gff_types = ["gff", "gff3"]
stat_annotations = ['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism']


def append_filename(path: str, s: str):
    """
    Append a string to a file name before the extension
    :param path: file path
    :param s: string to append
    :return: appended file name
    """
    seg = path.rsplit(".", 1)
    if len(seg) > 1:
        return seg[0] + s + "." + seg[1]
    else:
        return path + s


def get_args(sysargs: list):
    """
    Parse command line arguments
    :param sysargs: list of command line arguments (sys.argv[1:])
    :return: (input_path, input_type, output_path, output_type, jmespath, split, stats)
    """
    split = False
    jpath = None
    stats = False
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sysargs, 'vsiq:')
        for opt, val in opts:
            if opt == '-v':
                print(__versionstr__)
                exit(0)
            elif opt == '-s':
                split = True
            elif opt == '-q':
                jpath = val
            elif opt == '-i':
                stats = True

    except getopt.GetoptError as err:
        print("Argument error(" + str(err.opt) + "): " + err.msg, file=sys.stderr)
        args = []

    # Check for minimum number of arguments
    if len(args) < 4:
        print(usage, file=sys.stderr)
        exit(1)

    input_path = args[0]
    input_type = args[1]
    output_path = args[2]
    output_type = args[3]

    return input_path, input_type, output_path, output_type, jpath, split, stats


def convert(input_handle, input_type: str, output_path: str, output_type: str, jpath: str, split: bool, stats: bool):
    """
    Convert the input data to a file of specified format
    :param input_handle: File handle to read data from
    :param input_type: one of abi,abi-trim,ace,cif-atom,cif-seqres,clustal,embl,fasta,fasta-2line,fastq-sanger,fastq,
        fastq-solexa,fastq-illumina,genbank,gb,ig,imgt,nexus,pdb-seqres,pdb-atom,phd,phylip,pir,seqxml,sff,sff-trim,
        stockholm,swiss,tab,qual,uniprot-xml,gff3
    :param output_path: Path to output to
    :param output_type: Format to output as
    :param jpath: JMESPath selecting records to keep. The root is the list of records. The path must return a list of records.
    :param split: True to split into separate files, False otherwise
    :param stats: True to output record info to stdout, False otherwise
    :return: None
    """
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
        input_records = JMESPathGen.search(jpath, input_records)

    if isinstance(input_records, dict):
        # Support generating a new record in JMESPath
        input_records = SeqIO.SeqRecord(**input_records)

    if isinstance(input_records, SeqIO.SeqRecord):
        # Support returning single record from JMESPath
        input_records = (input_records,)

    if stats:
        print("##gff-version 3")

    output = None

    for i, record in enumerate(input_records):  # type: int, SeqIO.SeqRecord
        # TODO allow objects other than SeqRecord, transform to SeqRecord or handle special output (like if output format == txt|json, pretty print object)
        if isinstance(record, dict):
            # Support generating new records in JMESPath
            record = SeqIO.SeqRecord(**record)

        if stats:
            attributes = {'Name': [record.name]}
            for k, v in record.annotations.items():
                if k in stat_annotations:
                    if isinstance(v, list):
                        v = [str(a) for a in v]
                    else:
                        v = [str(v)]
                    attributes[k] = v
            attributes['desc'] = [record.description]
            print(gffutils.Feature(record.id, "biopython-convert", "sequence", start=1, end=len(record), attributes=attributes))

        if not output:
            # Open output file with file name suffix if splitting
            output = open(append_filename(output_path, "." + str(i)) if split else output_path, "w")

        if output_type in gff_types:
            # If output type is GFF, use gffutils library
            for feature in record.features:
                print(biopython_integration.from_seqfeature(feature), file=output)
        else:
            SeqIO.write(record, output, output_type)

        if split:
            # If splitting, open next file
            output.close()
            output = None


if __name__ == '__main__':
    in_path, *remaining_args = get_args(sys.argv[1:])

    with open(in_path, "r") as handle:
        convert(handle, *remaining_args)
