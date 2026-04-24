# Run from the svforge project root 

mkdir -p data_local/gnomad data_local/blacklist 

## gnomAD SV v4.1
wget -P data_local/gnomad/ https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz 

wget -P data_local/gnomad/ https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi

## ENCODE blacklist v2
wget -P data_local/blacklist/ https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz 

