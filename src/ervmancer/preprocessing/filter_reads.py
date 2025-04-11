import logging
import os
import gzip
from typing import Tuple, Optional


class FastqValidator:
    @staticmethod
    def is_gzipped(filepath: str) -> bool:
        with open(filepath, 'rb') as f:
            return f.read(2).startswith(b'\x1f\x8b')

    @staticmethod
    def is_fastq(filepath: str) -> bool:
        try:
            opener = gzip.open if FastqValidator.is_gzipped(filepath) else open
            with opener(filepath, 'rt') as f:
                header = f.readline().strip()
                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()

                return (header.startswith('@') and
                        len(sequence) > 0 and
                        plus.startswith('+') and
                        len(quality) == len(sequence))
        except Exception as e:
            logging.error(f"Error validating FASTQ file {filepath}: {str(e)}")
            return False

    @staticmethod
    def validate_input_file(input_path: str) -> str:
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file not found: {input_path}")

        if not input_path.lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
            raise ValueError(
                f"Invalid file extension for {input_path}. Must be .fastq, .fq, .fastq.gz, or .fq.gz")

        if not FastqValidator.is_fastq(input_path):
            raise ValueError(f"File {input_path} is not in valid FASTQ format")

        return input_path


class ReadFilter:
    def __init__(self, output_dir: str, r1_path: Optional[str] = None,
                 r2_path: Optional[str] = None, s1_path: Optional[str] = None):
        self.output_dir = output_dir
        self.validator = FastqValidator()

        if r1_path and r2_path:
            self.paired = True
            self.r1_path = self.validator.validate_input_file(r1_path)
            self.r2_path = self.validator.validate_input_file(r2_path)
            self.base_name = self._get_base_name(r1_path)
        elif s1_path:
            self.paired = False
            self.s1_path = self.validator.validate_input_file(s1_path)
            self.base_name = self._get_base_name(s1_path)

    def _get_base_name(self, filepath: str) -> str:
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        return base_name[:-6] if base_name.endswith('.fastq') else base_name

    def get_path(self, subdir: str, filename: str, ext: str) -> str:
        return os.path.join(self.output_dir, subdir, f"{filename}.{ext}")

    def validate_inputs(self) -> Tuple[str, bool]:
        return self.base_name, self.paired
