import difflib
import io
from unittest import TestCase
from hashlib import sha256
from tempfile import TemporaryDirectory
from pathlib import Path

from biopython_convert import convert


class TestConvert(TestCase):
    input_path = Path('test-data/has_plasmids.gbff')
    noseq_path = Path('test-data/no_seq.gbff')
    output_path = Path('test-data/outputs/')
    input_type = 'genbank'
    convert_type = 'embl'
    filter = "[?!(features[?type=='source'].qualifiers.plasmid)]"

    def setUp(self) -> None:
        self.workdir = TemporaryDirectory()

    def hash_file(self, path: Path):
        with path.open('rb') as output_handle:
            return sha256(output_handle.read()).hexdigest()

    def compare_files(self, a: Path, b: Path):
        with open(a) as a_stream:
            with open(b) as b_stream:
                diffs = list(difflib.unified_diff(a_stream.readlines(), b_stream.readlines(), fromfile=str(a), tofile=str(b)))
                self.assertEqual(0, len(diffs), ''.join(diffs))

    def tearDown(self) -> None:
        self.workdir.cleanup()

    def test_basic(self):
        output_path = Path(self.workdir.name, 'basic')
        convert(self.input_path, self.input_type, output_path, self.input_type)
        self.compare_files(Path.joinpath(self.output_path, 'basic'), output_path)

    def test_info(self):
        output_path = Path(self.workdir.name, 'basic')
        stats = io.StringIO()
        convert(self.input_path, self.input_type, output_path, self.input_type, stats=stats)
        self.compare_files(Path.joinpath(self.output_path, 'basic'), output_path)
        a = Path.joinpath(self.output_path, 'info')
        with open(a) as a_stream:
            diffs = list(difflib.unified_diff(a_stream.readlines(), stats.getvalue().splitlines(keepends=True), fromfile=str(a), tofile='stdout'))
            self.assertEqual(0, len(diffs), ''.join(diffs))

    def test_convert(self):
        output_path = Path(self.workdir.name, 'convert')
        convert(self.input_path, self.input_type, output_path, self.convert_type)
        self.compare_files(Path.joinpath(self.output_path, 'convert'), output_path)

    def test_noseq(self):
        output_path = Path(self.workdir.name, 'no_seq')
        convert(self.noseq_path, self.input_type, output_path, 'fasta')
        self.compare_files(Path.joinpath(self.output_path, 'no_seq'), output_path)

    def test_filter(self):
        output_path = Path(self.workdir.name, 'filter')
        convert(self.input_path, self.input_type, output_path, self.input_type, jpath=self.filter)
        self.compare_files(Path.joinpath(self.output_path, 'filter'), output_path)

    def test_split(self):
        output_path = Path(self.workdir.name, 'record')
        convert(self.input_path, self.input_type, output_path, self.input_type, split=True)
        files = list(x for x in Path(self.workdir.name).glob('*') if x.is_file())
        truth_files = list(Path.joinpath(self.output_path, 'split').glob('*'))
        files.sort()
        truth_files.sort()
        self.assertListEqual([f.name for f in truth_files], [f.name for f in files])
        for a, b in zip(truth_files, files):
            self.compare_files(a, b)

    def test_gff(self):
        output_path = Path(self.workdir.name, 'gff')
        convert(self.input_path, self.input_type, output_path, 'gff3')
        self.compare_files(Path.joinpath(self.output_path, 'gff'), output_path)

    def test_txt(self):
        output_path = Path(self.workdir.name, 'txt')
        convert(self.input_path, self.input_type, output_path, 'text', jpath='[*].annotations.taxonomy')
        self.compare_files(Path.joinpath(self.output_path, 'txt'), output_path)

    def test_ptt(self):
        output_path = Path(self.workdir.name, 'ptt')
        convert(self.input_path, self.input_type, output_path, 'text', jpath="[0].[join(' - 1..', [description, to_string(length(seq))]), join(' ', [to_string(length(features[?type=='CDS' && qualifiers.translation])), 'proteins']), join(`\"\\t\"`, ['Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym', 'Code', 'COG', 'Product']), (features[?type=='CDS' && qualifiers.translation].[join('..', [to_string(sum([location.start, `1`])), to_string(location.end)]), [location.strand][?@==`1`] && '+' || '-', length(qualifiers.translation[0]), (qualifiers.db_xref[?starts_with(@, 'GI')].split(':', @)[1])[0] || '-', qualifiers.gene[0] || '-', qualifiers.locus_tag[0] || '-', '-', '-', qualifiers.product[0] ] | [*].join(`\"\\t\"`, [*].to_string(@)) )] | []")
        self.compare_files(Path.joinpath(self.output_path, 'ptt'), output_path)

    def test_txt_stats(self):
        output_path = Path(self.workdir.name, 'txt')
        stats = io.StringIO()
        convert(self.input_path, self.input_type, output_path, 'text', jpath='[*].annotations.taxonomy', stats=stats)
        self.compare_files(Path.joinpath(self.output_path, 'txt'), output_path)
        self.assertEqual(1, len(stats.getvalue().splitlines(keepends=True)))

    def test_json(self):
        output_path = Path(self.workdir.name, 'json')
        convert(self.input_path, self.input_type, output_path, 'json')
        self.compare_files(Path.joinpath(self.output_path, 'json'), output_path)

    def test_yaml(self):
        output_path = Path(self.workdir.name, 'yaml')
        convert(self.input_path, self.input_type, output_path, 'yaml')
        self.compare_files(Path.joinpath(self.output_path, 'yaml'), output_path)

    def test_json_jpath(self):
        output_path = Path(self.workdir.name, 'json_jpath')
        convert(self.input_path, self.input_type, output_path, 'json', jpath='[*].{id: id, type: annotations.molecule_type}')
        self.compare_files(Path.joinpath(self.output_path, 'json_jpath'), output_path)

    def test_json_jpath_split(self):
        output_path = Path(self.workdir.name, 'json_jpath_split')
        convert(self.input_path, self.input_type, output_path, 'json', jpath='[*].{id: split(\'.\', id), type: annotations.molecule_type}')
        self.compare_files(Path.joinpath(self.output_path, 'json_jpath_split'), output_path)

    def test_json_jpath_extract(self):
        output_path = Path(self.workdir.name, 'json_jpath_extract')
        convert(self.input_path, self.input_type, output_path, 'json', jpath="[0].let({seq:seq}, &features[?type=='gene'].extract(seq, @))") # "[0].let({seq: seq}, &features[?type=='gene']|[0].{id: id, description: 'extracted sequence', seq: seq[::location]})")
        self.compare_files(Path.joinpath(self.output_path, 'json_jpath_extract'), output_path)

    def test_json_jpath_let(self):
        output_path = Path(self.workdir.name, 'json_jpath_let')
        convert(self.input_path, self.input_type, output_path, 'json', jpath='[*].let({type: annotations.molecule_type}, &{id: id, type: type})')
        self.compare_files(Path.joinpath(self.output_path, 'json_jpath_let'), output_path)

    def test_yaml_jpath(self):
        output_path = Path(self.workdir.name, 'yaml_jpath')
        convert(self.input_path, self.input_type, output_path, 'yaml', jpath='[*].{id: id, annotations: annotations}')
        self.compare_files(Path.joinpath(self.output_path, 'yaml_jpath'), output_path)

    def test_gentype(self):
        """
        Test handling the generator type within JMESPath functions
        """
        output_path = Path(self.workdir.name, 'gentype')
        convert(self.input_path, self.input_type, output_path, 'text', jpath="[*].[(annotations.organism || annotations.source), 'foo'] | [*].join('\t', @)")
        self.compare_files(Path.joinpath(self.output_path, 'gentype'), output_path)

    def test_creation(self):
        """
        Test handling JMESPath creating new SeqRecords
        """
        output_path = Path(self.workdir.name, 'faa')
        convert(self.input_path, self.input_type, output_path, 'fasta', jpath="""
[0].let({organism: (annotations.organism || annotations.source)}, &features[?type=='CDS' && qualifiers.translation].{id:
join('|', [
  (qualifiers.db_xref[?starts_with(@, 'GI')].['gi', split(':', @)[1]]),
  (qualifiers.protein_id[*].['ref', @]),
  (qualifiers.locus_tag[*].['locus', @]),
  join('', [':', [location][?strand==`-1`] && 'c' || '', to_string(sum([location.start, `1`])), '..', to_string(location.end)])
][][]),
seq: qualifiers.translation[0],
description: (organism && join('', [qualifiers.product[0] || qualifiers.protein_id[0], ' [', organism, ']']) || qualifiers.product[0])})
        """)
        self.compare_files(Path.joinpath(self.output_path, 'faa'), output_path)

    def test_creation2(self):
        """
        Test handling JMESPath creating new SeqRecords
        """
        output_path = Path(self.workdir.name, 'ffn')
        convert(self.input_path, self.input_type, output_path, 'fasta', jpath="""
[0].let({desc: description, seq: seq}, &features[?type=='gene'].{id:
join('|', [
  (qualifiers.db_xref[?starts_with(@, 'GI')].['gi', split(':', @)[1]]),
  (qualifiers.protein_id[*].['ref', @]),
  (qualifiers.locus_tag[*].['locus', @]),
  join('', [':', [location][?strand==`-1`] && 'c' || '', to_string(sum([location.start, `1`])), '..', to_string(location.end)])
][][]),
seq: extract(seq, @),
description: desc})
        """)
        self.compare_files(Path.joinpath(self.output_path, 'ffn'), output_path)