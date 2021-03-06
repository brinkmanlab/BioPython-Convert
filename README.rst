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
    sff, sff-trim, stockholm, swiss, tab, qual, uniprot-xml, gff3, txt, json, yaml

JMESPath_
---------
The root node for a query is a list of SeqRecord_ objects. The query can return a list with a subset of these or
a mapping, keying to the `constructor parameters`_ of a SeqRecord object.

If the formats are txt, json, or yaml, then the JMESPath resulting object will simply be dumped in those formats.

A web based tool is available to experiment with constructing queries in real time on your data. Simply convert your
dataset to JSON and load it into the `JMESPath playground`_ to begin composing your query. It supports loading JSON files
directly rather than trying to copy/paste the data.

`split()`_ and `let()`_ functions are available in addition to the JMESPath standard functions

`extract(Seq, SeqFeature)` is also made available to allow access to the `SeqFeature.extract()`_ function within the query

Examples:
    Append a new record::

        [@, [{'seq': 'AAAA', 'name': 'my_new_record'}]] | []

    Filter out any plasmids::

        [?!(features[?type=='source'].qualifiers.plasmid)]

    Keep only the first record::

        [0]

    Output taxonomy of each record (txt output)::

        [*].annotations.taxonomy

    Output json object containing id and molecule type::

        [*].{id: id, type: annotations.molecule_type}

    Convert dataset to PTT format using text output::

        [0].[join(' - 1..', [description, to_string(length(seq))]), join(' ', [to_string(length(features[?type=='CDS' && qualifiers.translation])), 'proteins']), join(`"\t"`, ['Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym', 'Code', 'COG', 'Product']), (features[?type=='CDS' && qualifiers.translation].[join('..', [to_string(sum([location.start, `1`])), to_string(location.end)]), [location.strand][?@==`1`] && '+' || '-', length(qualifiers.translation[0]), (qualifiers.db_xref[?starts_with(@, 'GI')].split(':', @)[1])[0] || '-', qualifiers.gene[0] || '-', qualifiers.locus_tag[0] || '-', '-', '-', qualifiers.product[0] ] | [*].join(`"\t"`, [*].to_string(@)) )] | []

		Convert dataset to faa format using fasta output::

				[0].let({org: (annotations.organism || annotations.source)}, &(features[?type=='CDS' && qualifiers.translation].{id:
				join('|', [
					(qualifiers.db_xref[?starts_with(@, 'GI')].['gi', split(':', @)[1]]),
					(qualifiers.protein_id[*].['ref', @]),
					(qualifiers.locus_tag[*].['locus', @]),
					join('', [':', [location][?strand==`-1`] && 'c' || '', to_string(sum([location.start, `1`])), '..', to_string(location.end)])
				][][]),
				seq: qualifiers.translation[0],
				description: (org && join('', [qualifiers.product[0], ' [', org, ']']) || qualifiers.product[0])}))

See CONTRIBUTING.rst_ for information on contributing to this repo.

.. _CONTRIBUTING.rst: CONTRIBUTING.rst
.. _JMESPath: http://jmespath.org/
.. _SeqRecord: https://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html
.. _constructor parameters: https://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html#__init__
.. _JMESPath playground: https://glenveegee.github.io/jmespath-edit/
.. _split(): https://github.com/jmespath/jmespath.py/issues/159
.. _let(): https://github.com/jmespath/jmespath.site/pull/6
.. _SeqFeature.extract(): https://biopython.org/docs/latest/api/Bio.SeqFeature.html#Bio.SeqFeature.SeqFeature.extract