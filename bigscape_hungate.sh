python /data/users/david/programmes/BiG-SCAPE/bigscape.py \
 --inputdir region_files \
 -o out_local \
 --pfam_dir pfam \
 -c 30 \
 --mix \
 --exclude_gbk_str PROKKA --cutoffs 0.05 0.1 0.3 0.7 0.9 --clan_cutoff 0.05 0.5  \
 --banned_classes PKSI PKSother NRPS Saccharides Terpene PKS-NRP_Hybrids \
 --mode glocal
