# Find_gRNA

### requirments
- crisprposfinder.py
- bedtools
- genome fasta file
- genome annotation gff3(gff) file
  

### find whole genome gRNA (Take NGG, for example)
`python crisprposfinder.py -s genome_file  -l 20 -p NGG -m 5 --offpath cas-offinder --threads 16 -o NGG_sgRNA.txt`

### convert to bed file 
`cut -f 1-3 NGG_sgRNA.txt > NGG_sgRNA.bed`

### Extract the CDS region 
`grep 'CDS' annotation_gff_file  > CDS.gff`

### Obtain gRNAs that overlap with the CDS region
`bedtools intersect -a NGG_sgRNA.bed -b CDS.gff | sort -V | uniq > NGG_sgRNA_intersect_CDS.bed`


### Obtain CDS region gRNAs number
`less NGG_sgRNA_intersect_CDS.bed | wc -l`

### Obtain whole genome gRNAs number
`less NGG_sgRNA.bed | wc -l`


### Obtain non-coding region gRNAs number
`non-coding region gRNA number = whole genome gRNAs number - CDS region gRNAs number`


