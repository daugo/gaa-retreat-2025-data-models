import logging
from datetime import datetime
from pathlib import Path

from pydantic import (
    BaseModel,
    FilePath,
    DirectoryPath,
    ValidationError,
)


class CDSCountsRow(BaseModel):
    seqid: str
    start: int
    end: int
    name: str
    score: str
    strand: str


COL_NAMES = ["seqid", "start", "end", "name", "score", "strand"]


def _validate_tss_rows(
    file_path: FilePath,
) -> (int, list):
    errors = []

    with open(file_path, "r") as f:
        for line_idx, line in enumerate(f, start=1):
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) != len(COL_NAMES):
                errors.append(
                    f"Line {line_idx}: Incorrect number of "
                    f"columns. Expected {len(COL_NAMES)}."
                )
                continue

            col_vals = dict(zip(COL_NAMES, cols))

            try:
                CDSCountsRow(
                    **col_vals,
                )
            except ValidationError as e:
                errors.append(f"Line {line_idx}: {e}")

    return 0, errors


def _valid_extension(file_path: FilePath) -> bool:
    allowed_extensions = [".bed"]
    path = Path(file_path)
    return path.suffixes[-1] in allowed_extensions


def validate(file_path: FilePath, out_dir: DirectoryPath) -> int:
    if _valid_extension(file_path):
        try:
            _, errors = _validate_tss_rows(file_path)
        except OSError as e:
            logging.error(f"Error reading file: {e}")
            return
    else:
        logging.error("File with .bed extension expected.")
        return

    if errors:
        time_info = datetime.now().strftime("%Y%m%d-%H%M%S")
        error_file = Path(out_dir) / f"validation_errors_{time_info}.txt"

        with open(error_file, "w") as f:
            f.write("\n".join(errors))

        logging.error(f"Found {len(errors)} validation errors.")
        logging.error(f"Validation errors written to {error_file}")
        return 1
    else:
        logging.info("No validation errors found.")
        return 0
