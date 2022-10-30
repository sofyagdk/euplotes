
python create_crassus_translated.py
makeblastdb -in ./crassus_transcripts.fasta -parse_seqids -out ./transcripts_db -dbtype 'nucl' -hash_index

python make_ribobits.py
blastn -db ./transcripts_db -query ./ribobits.fasta -outfmt 6 -out ./ribototrans.tsv -perc_identity 95
python kart.py > mapcrassus.txt

python architect.py > cras_arch.txt


