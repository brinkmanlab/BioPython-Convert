import io
from unittest import TestCase
from hashlib import sha256
from tempfile import TemporaryDirectory
from pathlib import Path

from biopython_convert import convert

class TestConvert(TestCase):
    input_path = Path('test-data/has_plasmids.gbff')
    input_type = 'genbank'
    basic_hash = 'c62e16b6233f24c7da0356d11bf37691d17aaf6ea19becaae0c7a264a010a282'
    convert_type = 'embl'
    convert_hash = 'b3331bfb8952db5dd043f2eaa216ca8369b14977f7cd46cd36bfa7d4626135b1'
    info_hash = '16d51442f49c12184c51b47241d2c3b52252f4e920e89538226452b9a9e3c548'
    filter = '[?!(features[?type==`source`].qualifiers.plasmid)]'
    record_hash = ('566e86a245e21c8d0610b26b86cbe9d0aa9a44d03994acd998700b37155ba5b5',
                   '23ac84c9b0367ff791eba3dc24ce580cafc105a97e60d26771c5ee71651b05fa',
                   '2dc99ebdd320b7a1ce5fa1e2d6ca15c63d16f92f8bec629755cb7de2bb22021a')
    gff_hash = '3f92c6fd87bba681866e17e210a96ec0e4745891f02813dea78d64af20c1ac8d'

    def setUp(self) -> None:
        self.workdir = TemporaryDirectory()

    def hash_file(self, path: Path):
        with path.open('rb') as output_handle:
            return sha256(output_handle.read()).hexdigest()

    def tearDown(self) -> None:
        self.workdir.cleanup()

    def test_basic(self):
        output_path = Path(self.workdir.name, 'basic')
        convert(Path(self.input_path), self.input_type, output_path, self.input_type)
        self.assertEqual(
            self.basic_hash,
            self.hash_file(output_path),
        )

    def test_info(self):
        output_path = Path(self.workdir.name, 'basic')
        stats = io.StringIO()
        convert(self.input_path, self.input_type, output_path, self.input_type, stats=stats)
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
        convert(self.input_path, self.input_type, output_path, self.convert_type)
        self.assertEqual(
            self.convert_hash,
            self.hash_file(output_path),
        )

    def test_filter(self):
        output_path = Path(self.workdir.name, 'filter')
        convert(self.input_path, self.input_type, output_path, self.input_type)
        self.assertEqual(
            self.basic_hash,
            self.hash_file(output_path),
        )

    def test_split(self):
        output_path = Path(self.workdir.name, 'record')
        convert(self.input_path, self.input_type, output_path, self.input_type, split=True)
        files = tuple(x for x in Path(self.workdir.name).glob('record.*') if x.is_file())
        self.assertEqual(len(self.record_hash), len(files))
        for record_path, record_hash in zip(files, self.record_hash):
            self.assertEqual(
                record_hash,
                self.hash_file(record_path),
            )

    def test_gff(self):
        output_path = Path(self.workdir.name, 'gff')
        convert(self.input_path, self.input_type, output_path, 'gff3')
        self.assertEqual(
            self.gff_hash,
            self.hash_file(output_path),
        )