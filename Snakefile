SAMPLES = ["SRR30810013",]
MIN_READ_LENGTH = 1000  # This will help drop some short reads that are likely to be noise and/or hard to assemble.
MIN_READS_PER_KB_PER_CONTIG = 10  # This will be highly dependent on the N50 of the input ONT reads.

rule all:
    input:
        final_assembly = expand("output/{sample}_assembly.fasta", sample=SAMPLES),
        assembly_mapped_reads = expand("output/{sample}_assembly.counts.tsv", sample=SAMPLES),
        gbk = expand("output/{sample}_annotations.gbk", sample=SAMPLES),
        gff = expand("output/{sample}_annotations.gff", sample=SAMPLES)


rule get_sample_fastqs:
    input:
        sample_name = expand("sample", sample=SAMPLES)
    output:
        sample_fastq = expand("{sample}.fastq", sample=SAMPLES)
    singularity:
        "docker://ncbi/sra-tools:latest"
    shell:
        """
        fastq-dump {sample_name}
        """

# Filtlong: https://github.com/rrwick/Filtlong
rule filter_reads_by_length:
    input:
        sample_fastq = expand("{sample}.fastq", sample=SAMPLES)
    output:
        filtered_reads = expand("{sample}/{sample}.filtered.fastq", sample=SAMPLES)
    singularity:
        "docker://staphb/filtlong:latest"
    shell:
        """
        filtlong --min_length {MIN_READ_LENGTH} {input.sample_fastq} > {output.filtered_reads}
        """


# Flye Assembler: https://github.com/mikolmogorov/Flye
rule raw_assembly_with_flye:
    input:
        filtered_reads = expand("{sample}/{sample}.filtered.fastq", sample=SAMPLES)
    output:
        flye_assembly_dir = directory(expand("{sample}/raw_assembly", sample=SAMPLES))
    singularity:
        "docker://staphb/flye:latest"
    shell:
        """
        flye --nano-hq {input.filtered_reads} --out-dir {output.flye_assembly_dir}
        """


rule reorganize_flye_results:
    input:
        flye_assembly_dir = expand("{sample}/raw_assembly", sample=SAMPLES)
    output:
        raw_assembly_file = expand("{sample}/flye_assembly.fasta", sample=SAMPLES)
    shell:
        """
        mv {input.flye_assembly_dir}/assembly.fasta {output.raw_assembly_file} &&
        rm -fr {input.flye_assembly_dir}
        """


# Medaka Assembly Polisher: https://github.com/nanoporetech/medaka
rule assembly_polishing_with_medaka:
    input:
        raw_assembly_file = expand("{sample}/flye_assembly.fasta", sample=SAMPLES),
        filtered_reads = expand("{sample}/{sample}.filtered.fastq", sample=SAMPLES)
    output:
        medaka_assembly_dir = directory(expand("{sample}/final_assembly", sample=SAMPLES))
    singularity:
        "docker://staphb/medaka:latest"
    shell:
        """
        medaka_consensus -i {input.filtered_reads} -d {input.raw_assembly_file} -o {output.medaka_assembly_dir} --bacteria
        """


rule map_reads_to_final_assembly_for_scoring:
    input:
        medaka_assembly_dir = expand("{sample}/final_assembly", sample=SAMPLES),
        filtered_reads = expand("{sample}/{sample}.filtered.fastq", sample=SAMPLES)
    output:
        mapped_to_assembly_bam = expand("{sample}/mapped_to_medaka.bam", sample=SAMPLES)
    singularity:
        "docker://staphb/minimap2:latest"
    shell:  # unsure about the --secondary=yes/no usage. Allowing might help collapse similar contigs, preventing would be more "sepcific"
        """
        minimap2 -x map-ont --secondary=no -a {input.medaka_assembly_dir}/consensus.fasta {input.filtered_reads} -o {output.mapped_to_assembly_bam}
        """


rule count_contig_mapping_with_samtools:
    input:
        mapped_to_assembly_bam = expand("{sample}/mapped_to_medaka.bam", sample=SAMPLES)
    output:
        sorted_bam = expand("{sample}/mapped_to_medaka.sorted.bam",sample=SAMPLES),
        contig_count_file = expand("{sample}/contig_counts.txt", sample=SAMPLES),
        mapped_to_assembly_sam = expand("{sample}/mapped_to_medaka.sam", sample=SAMPLES),
        mapped_to_assembly_index = expand("{sample}/mapped_to_medaka.sorted.bam.bai",sample=SAMPLES),
    singularity:
        # "docker://chrishah/samtools:latest"  # This docker is unexpectedly large and slow. Could make my own.
        # "docker://biopython/biopython:latest"  # This has samtools and is much faster?! But older and out of date.
        "docker://mgibio/samtools:1.15.1-buster"  # This is a smaller image, but still has samtools.
    shell:
        """
        samtools sort {input.mapped_to_assembly_bam} > {output.sorted_bam}
        samtools view {output.sorted_bam} > {output.mapped_to_assembly_sam}
        samtools index {output.sorted_bam}
        samtools idxstats {output.sorted_bam} > {output.contig_count_file}
        """


rule final_selection:
    input:
        contig_count_file = expand("{sample}/contig_counts.txt", sample=SAMPLES),
        sorted_bam= expand("{sample}/mapped_to_medaka.sorted.bam",sample=SAMPLES),
        medaka_assembly_dir = expand("{sample}/final_assembly", sample=SAMPLES)
    output:
        final_assembly = expand("output/{sample}_assembly.fasta", sample=SAMPLES),
        assembly_mapped_reads = expand("output/{sample}_assembly.counts.tsv", sample=SAMPLES)
    singularity:
        "docker://pandas/pandas:pip-all"
    shell:
        """
        python3 final_assembly_assessment.py \
        {input.contig_count_file} \
        {input.medaka_assembly_dir}/consensus.fasta \
        {input.sorted_bam} \
        {output.final_assembly} \
        {output.assembly_mapped_reads} \
        --min_reads_per_kb {MIN_READS_PER_KB_PER_CONTIG}
        """


# Prokka Annotation: https://github.com/tseemann/prokka
rule genome_annotation:
    input:
        final_assembly = expand("output/{sample}_assembly.fasta",sample=SAMPLES),
    params:
        sample = expand("{sample}", sample=SAMPLES)
    output:
        gbk = expand("output/{sample}_annotations.gbk", sample=SAMPLES),
        gff = expand("output/{sample}_annotations.gff", sample=SAMPLES)
    singularity:
        "docker://staphb/prokka:latest"
    shell:
        """
        prokka --outdir {params.sample}/annotations --prefix annotation {input.final_assembly} &&
        cp {params.sample}/annotations/annotation.gbk {output.gbk} &&
        cp {params.sample}/annotations/annotation.gff {output.gff}
        """