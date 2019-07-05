import io
from unittest import TestCase
from hashlib import sha256
from tempfile import TemporaryDirectory
from pathlib import Path

from biopython_convert import convert

class TestConvert(TestCase):
    input_type = 'genbank'
    basic_hash = '2808187bb8e2231545e4d2d7a27dc802df4d1f7c0e953a8399300b2df6b0c737'
    convert_type = 'embl'
    convert_hash = '5598cb679f5f6c31349968ddde3646fe97296da42ee528ed3f46dec3f5490cbd'
    info_hash = 'a611656c5a7e7f719c3d64f6b348b67c1abcb8ed56fa82f51fc90cbe2125e5f0'
    filter = '[?!(features[?type==`source`].qualifiers.plasmid)]'
    record_hash = ('8d02b2087c4cea42da7c5f0a69b7a40d544d953c1a9d611b97bd116cc1f8cd7f',
                   'e37ecc4288ae8b2c3bea25484326a69ced9679fa791162ed593064fdf535944d',
                   'e142d7e1fbd103c96e3b728e3b75f7af6955c97cdbddb87c3202f2c1e2f133d4')

    def setUp(self) -> None:
        self.input_handle = open('../test-data/has_plasmids.gbff')
        self.workdir = TemporaryDirectory()

    def hash_file(self, path: Path):
        with path.open('rb') as output_handle:
            return sha256(output_handle.read()).hexdigest()

    def tearDown(self) -> None:
        self.input_handle.close()
        self.workdir.cleanup()

    def test_basic(self):
        output_path = Path(self.workdir.name, 'basic')
        convert(self.input_handle, self.input_type, output_path, self.input_type)
        self.assertEqual(
            self.basic_hash,
            self.hash_file(output_path),
        )

    def test_info(self):
        output_path = Path(self.workdir.name, 'basic')
        stats = io.StringIO()
        convert(self.input_handle, self.input_type, output_path, self.input_type, stats=stats)
        self.assertEqual(
            self.basic_hash,
            self.hash_file(output_path),
        )
        self.assertEqual(
            self.info_hash,
            sha256(stats.getvalue().encode()).hexdigest(),
        )

    def test_convert(self):
        output_path = Path(self.workdir.name, 'convert')
        convert(self.input_handle, self.input_type, output_path, self.convert_type)
        self.assertEqual(
            self.convert_hash,
            self.hash_file(output_path),
        )

    def test_filter(self):
        output_path = Path(self.workdir.name, 'filter')
        convert(self.input_handle, self.input_type, output_path, self.input_type)
        self.assertEqual(
            self.basic_hash,
            self.hash_file(output_path),
        )

    def test_split(self):
        output_path = Path(self.workdir.name, 'record')
        convert(self.input_handle, self.input_type, output_path, self.input_type, split=True)
        files = tuple(x for x in Path(self.workdir.name).glob('record.*') if x.is_file())
        self.assertEqual(len(self.record_hash), len(files))
        for record_path, record_hash in zip(files, self.record_hash):
            self.assertEqual(
                record_hash,
                self.hash_file(record_path),
            )