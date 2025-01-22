from unittest import TextTestResult

# Custom data validation with Pydantic

## Introduction to Pydantic

Main purposes:

- Serialization/Derialization
- Parsing/Coercion
- Validation (byproduct, different from jsonschema; strict mode available)

We'll try to cover:

- Pydantic most common features
- Application to data validation in a semi-realistic scenario

Won't cover:

- Pydantic advanced features
  Pydantic implementation details and their impact on performance.
    - If interested: Samuel's PyCON US 2023
      talk: [How Pydantic V2 leverages Rust's Superpowers](https://www.youtube.com/watch?app=desktop&v=pWZw7hYoRVU)
- Pydantic integration with other libraries (e.g. FastAPI, pandera, etc.)
    - If curious: [Awesome Pydantic](https://github.com/Kludex/awesome-pydantic)
    - Interesting integration [`datamodel-code-generator`](https://github.com/koxudaxi/datamodel-code-generator/)

Useful links:

- [Official documentation](https://docs.pydantic.dev/latest/)
  Recommendations:
    - Make sure you are reading the docs for V2 ([migration guide](https://docs.pydantic.dev/dev/migration/)). Major
      changes to library API, additionally to performance gains by
      implementing pydantic core logic in Rust ([pydantic-core](https://github.com/pydantic/pydantic-core)).
    - The [concepts](https://docs.pydantic.dev/latest/concepts/models/) section is an excellent place
      to get a good understanding of all the features available.

- [Pydantic GitHub](https://github.com/pydantic/pydantic)
- [Type Hinting cheat sheet](https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html)

Relevant PEPs:
[pep-484](https://www.python.org/dev/peps/pep-0484/)
[pep-526](https://www.python.org/dev/peps/pep-0526/)
[pep-563](https://www.python.org/dev/peps/pep-0563/)
[pep-585](https://www.python.org/dev/peps/pep-0585/)
[pep-604](https://www.python.org/dev/peps/pep-0604/)

## Main concepts

### Models

- Primary way of defining schemas.
- As close as you can get in Python to **structs** in C like languages.
- Model definition using type hinting, similar to dataclasses (class typically containing mainly data)
- Classes that inherit from `pydantic.BaseModel`

Let's look at a simple example- validating a BED file containing transcription start sites (TSSs), one of the
intermediate outputs of the example workflow:

The file looks as follows:

```
1	11425	11426	ENSG00000290825(+)	*	+
1	28588	28589	ENSG00000243485(+)	*	+
1	29342	29343	ENSG00000310526(-)	*	-
1	29344	29345	ENSG00000310526(-)	*	-
1	29345	29346	ENSG00000310526(-)	*	-
1	29346	29347	ENSG00000310526(-)	*	-
1	29354	29355	ENSG00000310526(-)	*	-
1	29355	29356	ENSG00000310526(-)	*	-
1	29357	29358	ENSG00000310526(-)	*	-
1	29358	29359	ENSG00000310526(-)	*	-
...
```

_As a reminder, the BED file format is a tab-delimited text file format used to store genomic regions as coordinates and
associated annotations. The first three columns are the **chromosome (seq_id)**, **start**, and **end** positions of the
feature,
respectively. The fourth column is the **name** of the feature (displayed by genome browsers). The fifth column is the *
*score**
of the feature, and the sixth column is the **strand** of the
feature [(More information)](https://www.ensembl.org/info/website/upload/bed.html)._

We can use Pydantic to describe the schema of this file:

Let's focus on describing a row of the file:

An initial vague description could be:

```python
from pydantic import BaseModel


class TssRow(BaseModel):
    seqid: str
    start: int
    end: int
    name: str
    score: str
    strand: str

```

If we were to validate a row of the file, we could use the TssRow model as follows:

```python

external_data = {
    "seqid": "1",
    "start": 11425,
    "end": 11426,
    "name": "ENSG00000290825(+)",
    "score": "*",
    "strand": "+"
}

try:
    TssRow(**external_data)
except ValidationError as e:
    print(e.errors())

```

By default, Pydantic will try to coerce the data to the expected type. If coercion fails, a `ValidationError` will be
raised.

A string "453058" as the input to the `start`, an `int` field, will be converted/coerced to an `453058`.

Depending on the case this could be desirable or not. This behaviour can be disabled by setting the `strict` parameter
to
`BaseModel` validation methods.

```python
# This example will succeed
external_data = {
    "seqid": "1",
    "start": 11425,
    "end": 11426,
    "name": "ENSG00000290825(+)",
    "score": "*",
    "strand": "+"
}

try:
    TssRow.model_validate(external_data, strict=True)
except ValidationError as e:
    print(e.errors())

# This example will not
external_data = {
    "seqid": "1",
    "start": "11425",
    "end": 11426,
    "name": "ENSG00000290825(+)",
    "score": "*",
    "strand": "+"
}

try:
    TssRow.model_validate(external_data, strict=True)
except ValidationError as e:
    print(e.errors())


```

In the example of TssRow above, the model is too permissive. We can be more specific about the data types and
constraints of each field.

Also, the data model doesn't need to be a 1:1 mapping of the file format, in this case tabular.

```python

class BedRange(BaseModel):
    start: Annotated[int, Field(ge=0)]
    end: Annotated[int, Field(ge=1)]


class TssRow(BaseModel):
    seqid: str
    location: BedRange,
    name: Annotated[str, Field(pattern="^ENSG\\d{11}\\(\\+|\\-\\)$")]
    score: Literal["*"]
    strand: Literal["+", "-"]

```

### Validators

- Pydantic provides a couple of ways to define custom validators for fields.

docs: [Validators](https://docs.pydantic.dev/latest/concepts/validators/#field-validators)

```python
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

```

## Setting up this project

### 1. Clone this repository

```

git clone https://github.com/daugo/gaa-retreat-2025-data-models.git

```

### 2. Install uv (recommended to have more flexibility when interacting with the examples)

Using [uv](https://docs.astral.sh/uv/). Python package and project manager:

```

curl -LsSf https://astral.sh/uv/install.sh | sh

```

or

```

wget -qO- https://astral.sh/uv/install.sh | sh

```

More specific instructions for uv installation can be found [here](https://docs.astral.sh/uv/installation/).

If using codon, we can point you to an existent venv.

uv will be installed in `$HOME/.local/bin`. To add `$HOME/.local/bin` to your PATH, either restart your shell or run:

```

> source $HOME/.local/bin/env (sh, bash, zsh)
> source $HOME/.local/bin/env.fish (fish)

```

### 3. Build and run example project

From the root of the project (repository you clonned in step 1), run:

```

cd pydantic-data-models-examples
uv run ens-annot-validate

```

If everything is set up correctly, you should see the following help menu of the command line interface that we will use
as the entrypoint to interact with the data validation examples.

```
usage: ens-annot-validate [-h] [-t {ensembl-genome-gff3,tss-parquet,tss-bed}] [--output-dir Path] FILE-PATH

Example: ens-annot-validate -t ensembl-genome-gff3 <FILE-PATH>

Validate:
    - Ensembl genomic annotation files (Ensembl genome GFF3s)
    - Computed annotation features files used for Ensembl Regulation annotation Db
        - TSSs
        - Merged Exons
        - CDSs counts
        - Representative exons and CDSs

positional arguments:
  FILE-PATH             Input file path to validate

options:
  -h, --help            show this help message and exit
  -t {ensembl-genome-gff3,tss-parquet,tss-bed}
                        Input file type (default: ensembl-genome-gff3)
  --output-dir Path     (default: ./)

```

Project structure

```
.
├── pyproject.toml
├── README.md
├── src
│   └── pydantic_data_models_examples
│       ├── cli.py
│       ├── __init__.py
│       ├── models
│       │   ├── ensembl_gff3.py
│       │   ├── __init__.py
│       │   ├── tss_bed.py
│       │   └── tss_parquet.py
│       └── py.typed
└── uv.lock
```















