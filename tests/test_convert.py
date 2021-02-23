import difflib
import io
from unittest import TestCase
from hashlib import sha256
from tempfile import TemporaryDirectory
from pathlib import Path

from biopython_convert import convert


class TestConvert(TestCase):
    input_path = Path('test-data/has_plasmids.gbff')
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