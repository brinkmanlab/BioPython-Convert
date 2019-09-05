.. image:: https://zenodo.org/badge/195302632.svg
    :target: https://zenodo.org/badge/latestdoi/195302632

==================
BioPython-Convert
==================

Interconvert various file formats supported by BioPython.

Supports querying records with JMESPath.

Installation
------------
::

    pip install biopython-convert

or::

    conda install biopython-convert

or::

    git clone https://github.com/brinkmanlab/BioPython-Convert.git
    cd BioPython-Convert
    ./setup.py install

Use
---
::

    biopython.convert [-s] [-v] [-i] [-q JMESPath] input_file input_type output_file output_type
        -s Split records into seperate files
        -q JMESPath to select records. Must return list of SeqIO records or mappings. Root is list of input SeqIO records.
        -i Print out details of records during conversion
        -v Print version and exit

Supported formats
    abi, abi-trim, ace, cif-atom, cif-seqres, clustal, embl, fasta, fasta-2line, fastq-sanger, fastq,
    fastq-solexa, fastq-illumina, genbank, gb, ig, imgt, nexus, pdb-seqres, pdb-atom, phd, phylip, pir, seqxml,
    sff, sff-trim, stockholm, swiss, tab, qual, uniprot-xml, gff3

JMESPath_
---------
The root node for a query is a list of SeqRecord_ objects. The query can return a list with a subset of these or
a mapping, keying to the `constructor parameters`_ of a SeqRecord object.


Examples:
    Append a new record::

        [@, [{`seq`: `AAAA`, `name`: `my_new_record`}]] | []

    Filter out any plasmids::

        [?!(features[?type==`source`].qualifiers.plasmid)]

    Keep only the first record::

        [0]


See CONTRIBUTING.rst_ for information on contributing to this repo.

.. _CONTRIBUTING.rst: CONTRIBUTING.rst
.. _JMESPath: http://jmespath.org/
.. _SeqRecord: https://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html
.. _constructor parameters: https://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html#__init__
