#!/usr/bin/perl
perl /home/mmoldovan//pangenomes_2/proteinortho_v5.15/proteinortho5.pl -p=blastn -e=1e-25 -identity=70 -project=filtered_e25_id70 -cpus=8  ./data/transcriptomes/filtered/*.fasta
