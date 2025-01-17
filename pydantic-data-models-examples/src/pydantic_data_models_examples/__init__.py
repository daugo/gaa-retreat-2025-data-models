import gzip
import re
from pathlib import Path
from pprint import pprint

from pydantic_core._pydantic_core import ValidationError

from models import (
    ensembl_gff3,
)


def validate() -> str:
    file_path = Path("/Users/davidu/Repos/region-of-interest-data/data/human/Homo_sapiens.GRCh38.113.gff3.gz")

    errors = []
    gencode_basic_biotypes= set()
    gencode_basic_tags = set()

    with (gzip.open(file_path, "rt") as f):
        for line_idx, line in enumerate(f, start=1):
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            if len(columns) != 9:
                errors.append(f"Line {line_idx}: Incorrect number of columns.")
                continue

            cols_vals = {
                "seqid": columns[0],
                "source": columns[1],
                "type": columns[2],
                "location": ensembl_gff3.GenomicRange(
                    start=columns[3],
                    end=columns[4]
                ),
                "score": columns[5],
                "strand": columns[6],
                "phase": columns[7],
                "attributes": columns[8],
            }

            if re.search(r'tag=gencode_basic', columns[8]):
                try:
                    gencode_basic_transcript = (
                        ensembl_gff3.GencodeBasicTranscriptRow(**cols_vals)
                    )

                    # gencode_basic_biotypes.add(
                    #     gencode_basic_transcript.attributes.biotype)

                    [gencode_basic_tags.add(t) for t in
                        gencode_basic_transcript.attributes.tags]


                except ValidationError as e:
                    errors.append(f"Line {line_idx}: {str(e)}")
            else:
                try:
                    ensembl_gff3.Row(**cols_vals)
                except ValidationError as e:
                    errors.append(f"Line {line_idx}: {str(e)}")

    # pprint(gencode_basic_tags)

    with open("validation_errors.txt", "w") as f:
        f.write("\n".join(errors))



if __name__ == "__main__":
    print(validate())
