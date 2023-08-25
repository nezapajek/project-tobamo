#!/bin/bash
s=$1

source /etc/profile.d/sra-tools.sh

if [ -f "${s}" ] && [ -f "${s}_1.fastq.gz" ] && [ -f "${s}_2.fastq.gz" ]; then
    echo "$s already downloaded... paired-end"
elif [ -f "${s}" ] && [ -f "${s}.fastq.gz" ]; then
    echo "$s already downloaded... single-end"
else
    echo "Downloading $s"
    fasterq-dump -f -e 12 -p --split-files $(basename "$s") -O $(dirname "$s")
    if [ $? -ne 0 ]; then
        exit 1
    fi

    # check for wrong paired-end annotation of single-end SRRs
    # if only one file retrieved from SRA, make sure that properly named
    shopt -s nullglob
    generated_files=(${s}?*)
    numf=${#generated_files[@]}
    shopt -u nullglob

    if [ "$numf" -eq "1" ]; then
        # rename if needed
        if [ "${generated_files[0]}" != "${s}.fastq" ]; then
            echo "renaming to proper filename"
            mv -fv "${generated_files[0]}" "${s}.fastq"
        fi
    fi

    # compress in background
    if [ -f "${s}_1.fastq" ]; then
        # paired-end
        echo "Compressing paired-end"
        gzip -f ${s}_1.fastq
        if [ $? -ne 0 ]; then
            exit 2
        fi
        gzip -f ${s}_2.fastq
        if [ $? -ne 0 ]; then
            exit 3
        fi
    else
        # single-end
        echo "Compressing single-end"
        gzip -f ${s}.fastq
        if [ $? -ne 0 ]; then
            exit 4
        fi
    fi
fi

touch $s
echo "done"