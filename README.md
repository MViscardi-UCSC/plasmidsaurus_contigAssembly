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
1. The `SRR30810013.fastq` file will be loaded in the first rule of the Snakefile using fastq-dump from SRA Toolkit.
To avoid this, just have the `SRR30810013.fastq` file in the same directory as the Snakefile.
2. Assembly and annotation files will be saved in the `output` directory.

Note: All build assemblies and annotations are retained in the `output` directory.
The `contig_1` appears to be the biologically relevant assembly, as it is the only one with any rRNA annotated.

3. The SRR30810013 (or other short read archive samples defined in the Snakemake file) directory and the identically
named .fastq file(s) can be deleted after the run is complete.
They are retained here for reproducibility and step checking.

### General Process Overview:
1. Download the SRR30810013 sample from the SRA database using fastq-dump.
2. Remove short reads and low-quality sequences using Filtlong (R. Wick).
3. Assemble initial contigs using Flye (M. Kolmogorov).
4. Polish the contigs with Medaka (ONT).
5. Count rough coverage using Minimap2 (H. Li) and Samtools. *Here we drop any contigs with low coverage.*
6. Annotate the contigs using Prokka (T. Seemann).