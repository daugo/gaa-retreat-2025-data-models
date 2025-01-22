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

from pydantic_data_models_examples.models.ensembl_gff3 import GenomicRange


class BedRange(BaseModel):
    start: Annotated[int, Field(ge=0)]
    end: Annotated[int, Field(ge=1)]

    @model_validator(mode="after")
    def validate_coordinates(self) -> Self:
        if self.start >= self.end:
            raise ValueError(
                f"Start coordinate ({self.start}) is greater than end coordinate ({self.end})."
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


def validate(file: FilePath, out_dir: DirectoryPath) -> None:
    errors = []

    col_names = ["seqid", "start", "end", "name", "score", "strand"]

    with open(file, "r") as f:
        for line_idx, line in enumerate(f, start=1):
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) != len(col_names):
                errors.append(
                    f"Line {line_idx}: Incorrect number of "
                    f"columns. Expected {len(col_names)}."
                )
                continue

            col_vals = dict(zip(col_names, cols))
            col_vals["location"] = GenomicRange(
                start=cols[1],
                end=cols[2],
            )

            try:
                TssRow(
                    **col_vals,
                )
            except ValidationError as e:
                errors.append(f"Line {line_idx}: {e}")

    if errors:
        time_info = datetime.now().strftime("%Y%m%d-%H%M%S")
        error_file = Path(out_dir) / f"validation_errors_{time_info}.txt"

        with open(error_file, "w") as f:
            f.write("\n".join(errors))

        logging.error(f"Found {len(errors)} validation errors.")
        logging.error(f"Validation errors written to {error_file}")
    else:
        logging.info("No validation errors found.")
