import os


rule filter_contigs:
    input:
        condition=lambda wildcards: "iso"
        if os.path.exists(
            "results/{wildcards.accession}/02_{wildcards.accession}_spades_isolate_contigs.fasta"
        )
        and os.path.getsize(
            "results/{wildcards.accession}/02_{wildcards.accession}_spades_isolate_contigs.fasta"
        )
        > 0
        else "mega",
    output:
        add_accession="results/{accession}/04_{accession}_contigs_add_accession.fasta",
        filtered="results/{accession}/04_{accession}_contigs_filtered.fasta",
    log:
        logO="logs/filter_contigs/{accession}.log",
        logE="logs/filter_contigs/{accession}.err.log",
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/add_accession_to_fasta.py {input} {output.add_accession} {wildcards.accession} > {log.logO} 2> {log.logE}
        python workflow/scripts/filter_fasta_by_len.py {output.add_accession} {output.filtered} 600 8000 >> {log.logO} 2>> {log.logE}
        """
