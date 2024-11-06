# Marcus Viscardi, Nov 5th 2024
## Plasmidsaurus Application - Genome Assembly Assignment
***
### Requirements:
- Singularity/Docker
- Unix-based system (Linux, MacOS, WSL)
- Internet connection

### Running:
1. For the simplest implementation, run the Snakefile with the following command:
```bash
snakemake --use-singularity
```
### Notes:
1. Assembly and annotation files will be saved in the `output` directory.

Note: All build assemblies and annotations are retained in the `output` directory.
The `contig_1` appears to be the biologically relevant assembly, as it is the only one with any rRNA annotated.

2. The SRR30810013 (or other short read archive samples defined in the Snakemake file) directory and the identically
named .fastq file(s) can be deleted after the run is complete.
They are retained here for reproducibility and step checking.