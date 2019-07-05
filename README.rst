==================
BioPython-Convert
==================
Interconvert various file formats supported by BioPython. Supports querying records with JMESPath.
--------------------------------------------------------------------------------------------------

Installation
------------
::

    pip install biopython.convert

or::

    conda install biopython.convert

Use
---
::

    biopython.convert [-s] [-v] [-i] [-q JMESPath] input_file input_type output_file output_type
        -s Split records into seperate files
        -q JMESPath to select records. Must return list of SeqIO records. Root is list of input SeqIO records.
        -i Print out details of records during conversion
        -v Print version and exit

Supported formats
    abi, abi-trim, ace, cif-atom, cif-seqres, clustal, embl, fasta, fasta-2line, fastq-sanger, fastq,
    fastq-solexa, fastq-illumina, genbank, gb, ig, imgt, nexus, pdb-seqres, pdb-atom, phd, phylip, pir, seqxml,
    sff, sff-trim, stockholm, swiss, tab, qual, uniprot-xml, gff3

JMESPath examples:
    Filter out any plasmids::

        [?!(features[?type==`source`].qualifiers.plasmid)]

    Keep only the first record::

        [0]

