import logging
import sys
from typing import (
    Annotated,
    Literal,
)
from rich import print

from pydantic import (
    FilePath,
    Field,
    ValidationError,
)
from pydantic_settings import CliApp, CliPositionalArg, BaseSettings

from models import ensembl_gff3

log = logging.getLogger(__name__)

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


# class Validate(BaseSettings, cli_parse_args=True):
#     """Custom validation of ensembl annotation files. XX"""
#
#     file_path: CliPositionalArg[InputFilePath]
#     file_type: InputFileTypes
#     md5sum: MD5SumHash | None = None
#
#     def cli_cmd(self) -> None:
#         pass
#

#
# class EnsAnnotation(
#     BaseSettings,
#     cli_parse_args=True,
#     cli_prog_name="ens-annotation",
#     cli_kebab_case=True,
#     cli_use_class_docs_for_groups=True,
# ):
#     """Validate ensembl annotation files."""
#
#     check: CliSubCommand[Validate]
#     version: str
#
#     def cli_cmd(self) -> None:
#         CliApp.run_subcommand(self)


class EnsAnnotationValidate(
    BaseSettings,
    cli_prog_name="ens-annotation-validate",
    cli_kebab_case=True,
    cli_use_class_docs_for_groups=True,
):
    """
    Example: ens-annotation-validate ensembl-genome-gff3 <file_path>

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
            "merged-exons-parquet",
            "merged-exons-bed",
        ],
        Field(description="Input file type"),
    ]


def main() -> int:
    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) == 1:
        sys.argv.append("--help")

    try:
        cmd = CliApp.run(
            EnsAnnotationValidate,
            cli_exit_on_error=True,
            cli_args=sys.argv[1:],
        )
    except ValidationError as e:
        print("Args validation errors:\n")
        msgs = [
            f"{[e_details['loc'][0]]} {e_details['msg']}" for e_details in e.errors()
        ]
        print("\n".join(msgs))

        return 1
    except SystemExit:
        return 1

    print(cmd.file_path)

    if cmd.type == "ensembl-genome-gff3":
        print(f"Validating Ensembl genome GFF3 file: {cmd.file_path}")
        ensembl_gff3.validate(cmd.file_path)

    # try:
    #     cmd = CliApp.run(
    #         EnsAnnotation,
    #         cli_exit_on_error=True,
    #     )
    #
    # except (SystemExit, ValidationError):
    #     # print(e)
    #     return 1
    # args: EnsAnnotationValidate = cappa.parse(
    #     EnsAnnotationValidate,
    #     argv=argv,
    #     version="0.0.1",
    # )
    #
    # log.info(args)

    # cappa.invoke(
    #     EnsAnnotationValidate,
    #     argv=argv,
    #     version="0.0.1",
    # )

    # app = cappa.parse(EnsAnnotationValidate)
    #
    # print(app)
    # try:
    #     args = parse(EnsAnnotationValidate)
    # except ValidationError as e:
    #     print(e)
    #     return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
