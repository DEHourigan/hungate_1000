
### output for bacteria.
  phylophlan \
   -i /data/users/david/hungate_genomes \
   -d phylophlan \
   -f 02_tol.cfg \
   --diversity high \
   --fast \
   -o hun410 \
   --nproc 5 \
   --verbose 2>&1 | tee logs/phylophlan.log