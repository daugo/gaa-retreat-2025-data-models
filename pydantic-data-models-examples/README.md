# Custom data validation with Pydantic

## Setting up this project

### 1. Clone this repository

```
git clone https://github.com/daugo/gaa-retreat-2025-data-models.git
```

### 2. Install uv (Recommended for have more flexibility when interacting with the examples)

Using [uv](https://docs.astral.sh/uv/). Python package and project manager:

```
curl -LsSf https://astral.sh/uv/install.sh | sh
```

or

```
wget -qO- https://astral.sh/uv/install.sh | sh
```

More specific instructions for uv installation can be found [here](https://docs.astral.sh/uv/installation/).

If using codon, we can point you to a existent venv.

uv will be installed in `$HOME/.local/bin`. To add `$HOME/.local/bin` to your PATH, either restart your shell or run:

```
 > source $HOME/.local/bin/env (sh, bash, zsh)
 > source $HOME/.local/bin/env.fish (fish)
```

### 3. Build and run example project

Being in the root of the project (repository you clonned in step 1), run:

```
cd pydantic-data-models-examples
uv run ens-annot-validate
```

If everything is set up correctly, you should see the following the help menu of the command line interface we will use
as entrypoint for the data validation examples.

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

### Project structure

```
.
├── pyproject.toml
├── README.md
├── src
│   └── pydantic_data_models_examples
│       ├── cli.py => [entrypoint to the examples]
│       ├── __init__.py
│       ├── models => [Pydantic data models]
│       │   ├── ensembl_gff3.py
│       │   ├── __init__.py
│       │   ├── tss_bed.py
│       │   └── tss_parquet.py
│       └── py.typed
└── uv.lock
```












