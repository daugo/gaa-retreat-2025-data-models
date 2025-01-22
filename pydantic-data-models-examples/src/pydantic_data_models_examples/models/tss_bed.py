import logging
from datetime import datetime
from pathlib import Path
from typing import (
    Literal,
    Annotated,
)

from pydantic import (
    BaseModel,
    Field,
    field_validator,
    model_validator,
    FilePath,
    DirectoryPath,
    ValidationError,
)
from typing_extensions import Self


class BedRange(BaseModel):
    start: Annotated[int, Field(ge=0)]
    end: Annotated[int, Field(ge=1)]

    @model_validator(mode="after")
    def validate_coordinates(self) -> Self:
        if self.start >= self.end:
            raise ValueError(
                f"Start coordinate ({self.start}) is greater or equal than "
                f"end coordinate ({self.end})."
            )
        return self


class TssBedRange(BedRange):
    @model_validator(mode="after")
    def validate_tss_coordinates(self) -> Self:
        if (self.end - self.start) != 1:
            raise ValueError(
                f"TSSs should be a one base feature. (end ({self.end}) -"
                f" start({self.start})) != 1."
                f"{self.end})."
            )
        return self


class TssRow(BaseModel):
    seqid: str
    location: TssBedRange
    name: Annotated[str, Field(pattern="^ENSG\\d{11}\\(\\+|\\-\\)$")]
    score: Literal["*"]
    strand: Literal["+", "-"]

    @field_validator("seqid", mode="after")
    @classmethod
    def chr_prefix_not_present(cls, value: str) -> str:
        if value.lower().startswith("chr"):
            raise ValueError(
                f"Seqid value ('{value}') from an "
                f"Ensembl Generated GFF3 not expected to start "
                f"with'chr'"
            )
        return value


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
            col_vals["location"] = TssBedRange(
                start=cols[1],
                end=cols[2],
            )

            try:
                TssRow(
                    **col_vals,
                )
            except ValidationError as e:
                errors.append(f"Line {line_idx}: {e}")

    return 0, errors


def _valid_extension(file_path: FilePath) -> bool:
    allowed_extensions = [".bed"]
    path = Path(file_path)
    return path.suffixes[-1] in allowed_extensions


def validate(file_path: FilePath, out_dir: DirectoryPath) -> None:
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
    else:
        logging.info("No validation errors found.")
