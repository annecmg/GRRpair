# GRRpair

This program calculates the weighted gene repertoire relatedness (wGRR) between a pair of genomes. The measurement is described in [Kupczok et al. (2022)](https://doi.org/10.1101/2022.06.30.498228) and is based on [de Sousa et al. (2021)](https://doi.org/10.1093/molbev/msab044).

## Requirements
- BLAST+ v2.10 or higher, available [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) or in [bioconda](https://anaconda.org/bioconda/blast)
- powerneedle, see installation instructions in ``powerneedle`` directory

## Usage
```
python3 grr_for_pair.py genome1.faa genome2.faa output_folder
```
- ``genome1.faa`` and ``genome2.faa`` are multifasta file containg the protein sequences of each genome
- ``output_folder``: here the results will be stored, it must be created before running the script

## Results

A ``.grr`` file is generated in the ``output_folder`` containing 9 columns:
- Name of first genome
- Name of second genome
- GRR
- unnormalized GRR (the GRR values before dividing by the number of proteins)
- number of proteins in first genome
- number of proteins in second genome
- number of best bidirectional blast hits
- number of protein pairs that have an identity of at least 35%
- number of protein pairs that have an identity of at least 80%

To compare many pairs of genomes, it is possible to run ``grr_for_pair.py`` in a loop with the same output folder and concatenate the ``.grr`` files into a tab-delimited output file.

Additional files with intermediate results are generated:
- ``.blast`` - the blast results in output format 7, named ``subject_query.blast``
- ``.bbh`` - a list of best bidirectional blast hits
- ``.needle`` - the output of powerneedle containing the pairwise identities between all homologous pairs
