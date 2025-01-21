import logging.config
import sys
from typing import (
    Annotated,
    Literal,
)

from pydantic import FilePath, Field, ValidationError, AliasChoices, DirectoryPath
from pydantic_settings import CliApp, CliPositionalArg, BaseSettings

from pydantic_data_models_examples.models import (
    ensembl_gff3,
    tss_bed,
    tss_parquet,
)

from rich.logging import RichHandler

MD5SumHash = Annotated[
    str, Field(description="MD5 checksum hash.", pattern="^[a-fA-F0-9]{32}$")
]

InputFilePath = Annotated[
    FilePath,
    Field(description="Input file path to validate"),
]


InputFileTypes = Annotated[
    Literal[
        "ensembl-genome-gff3",
        "tss-parquet",
        "tss-bed",
        "merged-exons-parquet",
        "merged-exons-bed",
    ],
    Field(description="Input file type"),
]


class EnsAnnotationValidate(
    BaseSettings,
    cli_prog_name="ens-annot-validate",
    cli_kebab_case=True,
    cli_use_class_docs_for_groups=True,
):
    """
    Example: ens-annot-validate -t ensembl-genome-gff3 <FILE-PATH>

    Validate:
        - Ensembl genomic annotation files (Ensembl genome GFF3s)
        - Computed annotation features files used for Ensembl Regulation annotation Db
            - TSSs
            - Merged Exons
            - CDSs counts
            - Representative exons and CDSs
    """

    file_path: CliPositionalArg[InputFilePath]
    type: Annotated[
        Literal[
            "ensembl-genome-gff3",
            "tss-parquet",
            "tss-bed",
            # "merged-exons-parquet",
            # "merged-exons-bed",
        ],
        Field(
            validation_alias=AliasChoices("t"),
            default="ensembl-genome-gff3",
            description="Input file type",
        ),
    ]
    output_dir: DirectoryPath = "./"


def parse_cli() -> EnsAnnotationValidate:
    return CliApp.run(
        EnsAnnotationValidate,
        cli_exit_on_error=True,
        cli_args=sys.argv[1:],
    )


def _setup_logging():
    logging.basicConfig(
        format="[ %(asctime)s ] — [ %(funcName)s:%(lineno)d ] — %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            RichHandler(
                rich_tracebacks=True,
                tracebacks_show_locals=False,
                show_time=False,
            )
        ],
    )


def cli() -> int:
    _setup_logging()

    if len(sys.argv) == 1:
        sys.argv.append("--help")

    try:
        args = parse_cli()
    except ValidationError as e:
        msgs = [
            f"{[e_details['loc'][0]]} {e_details['msg']}" for e_details in e.errors()
        ]
        print("\nError:")
        print("\n".join(msgs))
        print("-" * 79)

        sys.argv = [sys.argv[0], "-h"]
        parse_cli()
        return 1

    except SystemExit:
        return 1

    if args.type == "ensembl-genome-gff3":
        logging.info(f"Validating Ensembl genome GFF3 file: {args.file_path}")
        ensembl_gff3.validate(args.file_path, args.output_dir)

    elif args.type == "tss-bed":
        logging.info(f"Validating TSS BED file: {args.file_path}")
        tss_bed.validate(args.file_path, args.output_dir)

    elif args.type == "tss-parquet":
        logging.info(f"Validating TSS Parquet file: {args.file_path}")
        tss_parquet.validate(args.file_path)

    return 0


if __name__ == "__main__":
    sys.exit(cli())
