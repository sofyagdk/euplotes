#!/bin/bash
#$ -cwd
#$ -pe smp 5
# python animal_collector.py
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_crassus/stephen/Euplotes/Transcriptomes/E.crassus/EU6_CAGATC_L005_R1.fastq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/cras_1_R1.fastq
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_crassus/stephen/Euplotes/Transcriptomes/E.crassus/EU5_GCCAAT_L005_R1.fastq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/cras_2_R1.fastq
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_focardii/stephen/Euplotes/E.focardii/MMETSP0205-Euplotes-focardii-TN1.1.fastq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/foca_R1.fastq
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_focardii/stephen/Euplotes/E.focardii/MMETSP0205-Euplotes-focardii-TN1.2.fastq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/foca_R2.fastq
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_harpa/stephen/Euplotes/Transcriptomes/E.harpa/MMETSP0213-Euplotes-harpa-FSP1_4.1.fastq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/harpa_R1.fastq
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_harpa/stephen/Euplotes/Transcriptomes/E.harpa/MMETSP0213-Euplotes-harpa-FSP1_4.2.fastq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/harpa_R2.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_minuta/stephen/Euplotes/Transcriptomes/minuta/*_1.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/minu_R1.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_minuta/stephen/Euplotes/Transcriptomes/minuta/*_2.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/minu_R2.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_octocarinatus/stephen/Euplotes/Transcriptomes/E.octocarinatus/octo/*_1.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/octo_R1.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_octocarinatus/stephen/Euplotes/Transcriptomes/E.octocarinatus/octo/*_2.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/octo_R2.fastq
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_raikovi/stephen/Euplotes/Transcriptomes/E.raikovi/E.raikovi_I_L1_1.fq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/raik_R1.fastq
# cp /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_raikovi/stephen/Euplotes/Transcriptomes/E.raikovi/E.raikovi_I_L1_2.fq /mnt/gamma/user/sofya/data/RNA_seq/concatenated/raik_R2.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_rariseta/stephen/Euplotes/Transcriptomes/rariseta/*_1.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/rari_R1.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_rariseta/stephen/Euplotes/Transcriptomes/rariseta/*_2.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/rari_R2.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_euryhalinus/*_1.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/eury_R1.fastq
# cat /mnt/gamma/user/sofya/data/RNA_seq/extracted/E_euryhalinus/*_2.fq > /mnt/gamma/user/sofya/data/RNA_seq/concatenated/eury_R2.fastq
# bowtie2-build ./animal_fastas/cras.fasta ./animal_fastas/cras_index 
# bowtie2-build ./animal_fastas/minu.fasta ./animal_fastas/minu_index 
# bowtie2-build ./animal_fastas/foca.fasta ./animal_fastas/foca_index 
# bowtie2-build ./animal_fastas/octo.fasta ./animal_fastas/octo_index 
# bowtie2-build ./animal_fastas/rari.fasta ./animal_fastas/rari_index 
# bowtie2-build ./animal_fastas/raik.fasta ./animal_fastas/raik_index 
# bowtie2-build ./animal_fastas/harp.fasta ./animal_fastas/harp_index 
# bowtie2-build ./animal_fastas/petz.fasta ./animal_fastas/petz_index 
# bowtie2-build ./animal_fastas/eury.fasta ./animal_fastas/eury_index 
# bowtie2 --sensitive-local --threads=5 --local --no-unal -x cras_index  -U /mnt/gamma/user/sofya/data/RNA_seq/concatenated/cras_1_R1.fastq,/mnt/gamma/user/sofya/data/RNA_seq/concatenated/cras_2_R1.fastq -S cras.sam  
# bowtie2 --sensitive-local --threads=5 --local --no-unal -x minu_index  -1 /mnt/gamma/user/sofya/data/RNA_seq/concatenated/minu_R1.fastq -2 /mnt/gamma/user/sofya/data/RNA_seq/concatenated/minu_R2.fastq -S minu.sam 
#kallisto quant -i minu.idx -o minu ../../data/Euplotes/RNA_seq/concatenated/minu_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/minu_R2.fastq


####done on local computer
# kallisto quant -i foca.idx -o foca ../../data/Euplotes/RNA_seq/concatenated/foca_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/foca_R2.fastq

# kallisto quant -i rari.idx -o rari ../../data/Euplotes/RNA_seq/concatenated/rari_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/rari_R2.fastq

# kallisto quant -i raik.idx -o raik ../../data/Euplotes/RNA_seq/concatenated/raik_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/raik_R2.fastq

# kallisto quant -i octo.idx -o octo ../../data/Euplotes/RNA_seq/concatenated/octo_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/octo_R2.fastq

# kallisto quant -i harp.idx -o harp ../../data/Euplotes/RNA_seq/concatenated/harpa_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/harpa_R2.fastq

# kallisto quant -i eury.idx -o eury ../../data/Euplotes/RNA_seq/concatenated/eury_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/eury_R2.fastq

# kallisto quant -i cras.idx -o cras --single -l 50 -s 1  ../../data/Euplotes/RNA_seq/concatenated/cras_1_R1.fastq ../../data/Euplotes/RNA_seq/concatenated/cras_2_R1.fastq



## abundance_full.tsv is for the whole transcriptome whereas abundance_all.tsv is only for orthogroups 