# plaston needs fasta suffix. so remove files when done
for file in gtdbtk_fna/*.fasta ;
do 
    platon --verbose \
    --db /data/san/data0/databases/platon/db/ \
    --output platon_out \
    --mode accuracy \
    --threads 8 $file
done