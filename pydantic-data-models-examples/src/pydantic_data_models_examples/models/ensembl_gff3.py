import gzip
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import (
    Annotated,
    Literal,
    Any,
    Self,
)

from pydantic import (
    BaseModel,
    Field,
    field_validator,
    model_validator,
    AliasChoices,
    ValidationError,
    FilePath,
    DirectoryPath,
    TypeAdapter,
)


# More info:
# http://www.ensembl.org/Help/Faq?id=468#:~:text=The%20Ensembl%20automatic%20annotation%20system,incorporate%20manual%20annotation%20from%20Havana.


EnsemblBiotype = Annotated[
    Literal[
        "IG_C_gene",
        "IG_C_pseudogene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_J_pseudogene",
        "IG_V_gene",
        "IG_V_pseudogene",
        "IG_pseudogene",
        "Mt_rRNA",
        "Mt_tRNA",
        "TEC",
        "TR_C_gene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_J_pseudogene",
        "TR_V_gene",
        "TR_V_pseudogene",
        "artifact",
        "lncRNA",
        "miRNA",
        "misc_RNA",
        "nonsense_mediated_decay",
        "processed_pseudogene",
        "processed_transcript",
        "protein_coding",
        "protein_coding_LoF",
        "pseudogene",
        "rRNA",
        "rRNA_pseudogene",
        "retained_intron",
        "ribozyme",
        "sRNA",
        "scRNA",
        "scaRNA",
        "snRNA",
        "snoRNA",
        "transcribed_processed_pseudogene",
        "transcribed_unitary_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "translated_processed_pseudogene",
        "unitary_pseudogene",
        "unprocessed_pseudogene",
        "vault_RNA",
    ],
    Field(),
]

TranscriptTag = Annotated[
    Literal[
        "gencode_basic",
        "Ensembl_canonical",
        "gencode_primary",
        "MANE_Select",
        "MANE_Plus_Clinical",
    ],
    Field(),
]


def na_to_none(value: Any) -> Any:
    if value.lower() == "na":
        return None
    return value


class GenomicRange(BaseModel):
    start: Annotated[int, Field(ge=1)]
    end: Annotated[int, Field(ge=1)]

    @model_validator(mode="after")
    def validate_coordinates(self) -> Self:
        if self.start > self.end:
            raise ValueError(
                f"Start coordinate ({self.start}) is greater than end coordinate ({self.end})."
            )
        return self


class Row(BaseModel):
    seqid: Annotated[str, Field(pattern="^[a-zA-Z0-9.:^*$@!+_?-|]+$")]
    source: str
    type: str
    location: GenomicRange
    score: float | None
    strand: Literal["+", "-", ".", "?"]
    phase: Literal["0", "1", "2"] | None
    attributes: dict[str, list[str]]

    @field_validator("score", "phase", mode="before")
    @classmethod
    def handle_dot_as_none(cls, value: Any) -> Any:
        if value == ".":
            return None

        return value

    @field_validator("seqid", mode="after")
    @classmethod
    def chr_prefix_not_present(cls, value: str) -> str:
        if value.lower().startswith("chr"):
            raise ValueError(
                f"Seqid (col 1) value ('{value}') from an "
                f"Ensembl Generated GFF3 not expected to start "
                f"with'chr'"
            )

    @field_validator("attributes", mode="before")
    @classmethod
    def parse_attributes(cls, value: Any) -> dict[str, list[str]]:
        attributes = {}
        for attribute in value.split(";"):
            key, value = attribute.split("=")
            attributes[key] = value.split(",")
        return attributes


# ID=gene:ENSG00000239945;biotype=lncRNA;description=novel transcript;gene_id=ENSG00000239945;logic_name=havana_homo_sapiens;version=1
class GeneLikeAttributes(BaseModel):
    id: Annotated[str, Field(validation_alias=AliasChoices("ID"), pattern="^gene:")]
    biotype: EnsemblBiotype
    description: str | None = None
    gene_id: str
    logic_name: str
    version: str


class GeneLikeRow(Row):
    attributes: GeneLikeAttributes

    @field_validator("attributes", mode="before")
    @classmethod
    def parse_attributes(cls, value: Any) -> GeneLikeAttributes:
        attributes = {}
        for attribute in value.split(";"):
            key, value = attribute.split("=")
            vals = value.split(",")

            try:
                assert len(vals) == 1, (
                    f"Expected single value for {key} in attributes column, got {vals}"
                )
            except AssertionError as e:
                raise ValueError(e)

            attributes[key] = vals[0]

        return GeneLikeAttributes(**attributes)


# ID=transcript:ENST00000623180;Parent=gene:ENSG00000280279;Name=LINC02887-201;biotype=lncRNA;tag=gencode_basic,Ensembl_canonical;transcript_id=ENST00000623180;transcript_support_level=5;version=1
class GencodeBasicTranscriptAttributes(BaseModel):
    id: Annotated[
        str, Field(validation_alias=AliasChoices("ID"), pattern="^transcript:")
    ]
    parent: str = Field(validation_alias=AliasChoices("Parent"), pattern="^gene:")
    name: str | None = Field(validation_alias=AliasChoices("Name"), default=None)
    biotype: EnsemblBiotype
    tags: Annotated[list[TranscriptTag], Field(alias="tag")]
    transcript_id: str
    transcript_support_level: str | None = None
    version: int


class GencodeBasicTranscriptRow(Row):
    attributes: GencodeBasicTranscriptAttributes

    @field_validator("attributes", mode="before")
    @classmethod
    def parse_attributes(cls, value: Any) -> GencodeBasicTranscriptAttributes:
        attributes = {}
        for attribute in value.split(";"):
            key, value = attribute.split("=")
            vals = value.split(",")
            if key == "tag":
                attributes[key] = vals
            else:
                try:
                    assert len(vals) == 1, (
                        f"Expected single value for {key} "
                        f"in attributes column,"
                        f" got {vals}"
                    )
                except AssertionError as e:
                    raise ValueError(e)

                attributes[key] = vals[0]

        return GencodeBasicTranscriptAttributes(**attributes)


GFF3FilePath = Annotated[
    FilePath,
    Field(
        description="Ensembl GFF3 file path.",
        pattern=".*\\.gff3\\.gz$",
    ),
]


gff3_file_adapter = TypeAdapter(GFF3FilePath)


def _validate_gff3_rows(f, errors: list) -> None:
    for line_idx, line in enumerate(f, start=1):
        if line.startswith("#"):
            continue
        columns = line.strip().split("\t")
        if len(columns) != 9:
            errors.append(
                f"Line {line_idx}: Incorrect number of "
                f"columns. Expected {len(COL_NAMES)}."
            )
            continue

        cols_vals = dict(zip(COL_NAMES, columns))

        cols_vals["location"] = GenomicRange(start=columns[3], end=columns[4])

        if re.search(r"ID=gene:", columns[8]):
            try:
                GeneLikeRow(**cols_vals)
            except ValidationError as e:
                errors.append(f"Line {line_idx}: {str(e)}")
        elif re.search(r"tag=gencode_basic", columns[8]):
            try:
                GencodeBasicTranscriptRow(**cols_vals)

            except ValidationError as e:
                errors.append(f"Line {line_idx}: {str(e)}")
        else:
            try:
                Row(**cols_vals)
            except ValidationError as e:
                errors.append(f"Line {line_idx}: {str(e)}")


COL_NAMES = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]


def _is_valid_gff3(file_path: FilePath):
    allowed_compressions = [".gz"]

    path = Path(file_path)
    suffixes = path.suffixes[-2:]

    if suffixes[0] == ".gff3" and suffixes[1] in allowed_compressions:
        return True

    if len(suffixes) == 1 and suffixes[0] == ".gff3":
        return True

    return False


def _is_compressed_gff3(file_path: FilePath):
    path = Path(file_path)
    suffixes = path.suffixes[-2:]

    logging.info(f"suffixes: {suffixes}")

    if _is_valid_gff3(path) and len(suffixes) == 2 and suffixes[1] == ".gz":
        return True

    return False


def validate(ensembl_gff3_file: FilePath, out_dir: DirectoryPath) -> None:
    errors = []

    if _is_compressed_gff3(ensembl_gff3_file):
        with gzip.open(ensembl_gff3_file, "rt") as f:
            _validate_gff3_rows(f, errors)

    elif _is_valid_gff3(ensembl_gff3_file):
        with open(ensembl_gff3_file, "rt") as f:
            _validate_gff3_rows(f, errors)
    else:
        logging.error("File must be a GFF3 file with a .gff3 or .gff3.gz extension.")
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
