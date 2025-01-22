import logging
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import (
    Annotated,
    Self,
    Literal,
)

import polars as pl
from pydantic import (
    FilePath,
    DirectoryPath,
    ValidationError,
    Field,
    BaseModel,
    PositiveInt,
    model_validator,
)


class TSSParquetByRow(BaseModel):
    reference_name: str = Field(validation_alias="referenceName")
    start: Annotated[PositiveInt, Field(validation_alias="tssStart")]
    end: Annotated[PositiveInt, Field(validation_alias="tssEnd")]
    gene_id: str = Field(validation_alias="geneId")
    strand: Literal["FORWARD", "REVERSE"]

    @model_validator(mode="after")
    def validate_coordinates(self) -> Self:
        if self.start > self.end:
            raise ValueError(
                f"Start coordinate ({self.start}) "
                f"is greater than end coordinate ({self.end})."
            )
        return self


class TSSParquetSchemaByCol(BaseModel):
    reference_name: Sequence[str] = Field(validation_alias="referenceName")
    start: Sequence[PositiveInt] = Field(validation_alias="tssStart")
    end: Sequence[PositiveInt] = Field(validation_alias="tssEnd")
    gene_id: Sequence[str] = Field(validation_alias="geneId")
    strand: Sequence[Literal["FORWARD", "REVERSE"]]


def write_validation_errors(errors: list[str], out_dir: Path) -> None:
    time_info = datetime.now().strftime("%Y%m%d-%H%M%S")
    error_file = out_dir / f"validation_errors_{time_info}.txt"

    with open(error_file, "w") as f:
        f.write("\n".join(errors))

    logging.error(f"Found {len(errors)} validation errors.")
    logging.error(f"Validation errors written to {error_file}")


def validate(file_path: FilePath, out_dir: DirectoryPath) -> None:
    errors = []

    df = pl.read_parquet(file_path)

    for row_idx, row in enumerate(df.to_dicts(), start=1):
        try:
            TSSParquetByRow.model_validate(row)
        except ValidationError as e:
            errors.append(f"Line {row_idx}: {e}")

    try:
        TSSParquetSchemaByCol.model_validate(df.to_dict(as_series=False))
    except ValidationError as e:
        errors.append(f"Column-related error: {e}")

    if errors:
        write_validation_errors(errors, out_dir)
