import logging
import os
import gzip
from typing import Tuple, Optional


class FastqValidator:
    @staticmethod
    def is_gzipped(filepath: str) -> bool:
        """Determines whether a filepath contains a gunzipped file/folder using magic numbers

        Args:
            filepath (str): Absolute filepath.

        Returns:
            bool: Whether or not the destination filepath is gunzipped.
        """
        with open(filepath, 'rb') as f:
            return f.read(2).startswith(b'\x1f\x8b')

    @staticmethod
    def is_fastq(filepath: str) -> bool:
        """Determines whether the given filepath is a fastq file with valid format.

        Args:
            filepath (str): Absolute filepath.

        Returns:
            bool: Whether or not the filepath contains a fastq file.
        """
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
        """Helper function to validate if the fastq file is in the correct format.

        Args:
            input_path (str): Absolute filepath.

        Raises:
            FileNotFoundError: Path does not resolve to a valid filepath.
            ValueError: Invalid file extension was given, must be the four kwown fastq extensions
            ValueError: Input path is invalid.

        Returns:
            str: Returns input path for further use.
        """
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
                 r2_path: Optional[str] = None, s1_path: Optional[str] = None,
                 csv_path: Optional[str] = None):
        self.output_dir = output_dir
        self.validator = FastqValidator()
        self.paired = False
        self.base_name = None

        # Handle different input combinations
        if r1_path and r2_path:
            # paired
            self.paired = True
            self.r1_path = self.validator.validate_input_file(r1_path)
            self.r2_path = self.validator.validate_input_file(r2_path)
            self.base_name = self._get_base_name(r1_path)
        elif s1_path:
            # any type of single strand
            self.paired = False
            self.s1_path = self.validator.validate_input_file(s1_path)
            self.base_name = self._get_base_name(s1_path)
        elif csv_path:
            # entrypoint 3 - csv/other method
            self.paired = False
            self.base_name = self._get_base_name(csv_path)

    def _get_base_name(self, filepath: str) -> str:
        """Gets the base name without the extension of the filepath.

        Args:
            filepath (str): Absolute filepath.

        Returns:
            str: base name of designated file.
        """
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        if base_name.endswith('.fastq'):
            return base_name[:-6]
        elif base_name.endswith('.fastq.gz'):
            return base_name[:-9]
        elif base_name.endswith('.fq'):
            return base_name[:-3]
        elif base_name.endswith('.fq.gz'):
            return base_name[:-6]
        elif base_name.endswith('.csv'):
            return base_name[:-4]
        return base_name

    def get_path(self, subdir: str, filename: str, ext: str) -> str:
        """Grabs a path using Python's os library to create an abspath

        Args:
            subdir (str): subdirectory to join
            filename (str): file name to use for output
            ext (str): extension to use for output

        Returns:
            str: _description_
        """
        return os.path.join(self.output_dir, subdir, f"{filename}.{ext}")

    def validate_inputs(self) -> Tuple[str, bool]:
        """Validates the inputs by ensuring the minimum required attributes for a ReadFilter object are instantiated.

        Returns:
            Tuple[str, bool]: Base name of designated file path and whether or not this is paired or not (R1 + R2 vs S1)
        """
        return self.base_name, self.paired
