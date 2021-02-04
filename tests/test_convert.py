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
    filter = '[?!(features[?type==`source`].qualifiers.plasmid)]'

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
        convert(Path(self.input_path), self.input_type, output_path, self.input_type)
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
        convert(self.input_path, self.input_type, output_path, self.input_type)
        self.compare_files(Path.joinpath(self.output_path, 'filter'), output_path)

    def test_split(self):
        output_path = Path(self.workdir.name, 'record')
        convert(self.input_path, self.input_type, output_path, self.input_type, split=True)
        files = tuple(x for x in Path(self.workdir.name).glob('*') if x.is_file())
        truth_files = list(Path.joinpath(self.output_path, 'split').glob('*'))
        self.assertEqual(len(truth_files), len(files))
        for a, b in zip(truth_files, files):
            self.compare_files(a, b)

    def test_gff(self):
        output_path = Path(self.workdir.name, 'gff')
        convert(self.input_path, self.input_type, output_path, 'gff3')
        self.compare_files(Path.joinpath(self.output_path, 'gff'), output_path)