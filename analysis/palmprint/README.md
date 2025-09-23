# download getorf (emboss)
 conda install bioconda::emboss

# build docker (gcc:7.2)
 docker build -t palmscan .

# docker help
 docker run -it --rm -v "$PWD":/usr/src/palmscan/data -w /usr/src/palmscan/data palmscan -help

# run docker
 docker run -it --rm -v "$PWD":/usr/src/palmscan/data -w /usr/src/palmscan/data palmscan -search_pssms seqs.fasta -tsv hits.tsv -report_pssms report -fasta pp.fa
    #seqs.fasta = aa zaporedje
    #report = seqA, seqB, seqC 
    #pp.fa = whole palmprint