############# Methylation site analysis #############

# Trim
trim_galore -q 25 --phred33 --stringency 3 --length 36  --paired R1.fq.gz R2.fq.gz --gzip -o

# Identify methylation site
bsmap -a R1_val_1.fq.gz -b R2_val_2.fq.gz -d Homo_sapiens.GRCh38.dna_chromosome.fa -p 10 -v 0.05 2>report.txt | samtools view -b -o sample.bam
sambamba sort \
--tmpdir tmp \
-t 50 \
-o sample_sorted.bam sample.bam

sambamba markdup \
--overflow-list-size 1000000 \
--tmpdir tmp \
-t 50 \
sample_sorted.bam \
sample_dedup.bam \
2> sample_markdup_report.txt

MethylDackel extract \
Homo_sapiens.GRCh38.dna_chromosome.fa \
sample_dedup.bam \
--opref sample

ethratio.py -r -z -p -u -d Homo_sapiens.GRCh38.dna_chromosome.fa -m 1 sample_dedup.bam -o sample_meth.BSmap.txt -w sample.BSmap.wig

## Filter out CpG sites with coverage no less than 10-fold 
awk '{{if($8>=10){{printf("%s\t%d\t%s\t%s\t%d\t%d\t%.4f\n",$1,$2,$3,$4,$7,$8,$5)}}}}' sample_meth.BSmap.txt | sed '/chr/d'|sed '1ichr\tpos\tstrand\ttype\tnum_C\tdepth\tmeth_level'> sample_meth.BSmap_dep10.txt
awk '{ if (NR > 1) {print $1"\t"$2"\t"$2+1"\t"$7}}' sample_meth.BSmap_dep10.txt > sample_meth.BSmap_dep10.bedGraph

# Distribution of CpG sites
awk '{if ($8>=10) {printf("%s %d %d %.3f\n",$1,$2-1,$2,$5)}}' sample_meth.BSmap.txt | sed '/chr/d' > sample_meth.bsmap.bedgraph
awk '{print $1"\t"$2"\t"$3"\t"$4}' ${id}.bsmap.bedgraph2 > ${id}.bsmap.bed
intersectBed -a ${id}.bsmap.bed -b Homo_sapiens.GRCh38.dna_chromosome.all.exon.sorted.merged.bed -wa > sample.bsmap_exon.bed
intersectBed -a ${id}.bsmap.bed -b Homo_sapiens.GRCh38.dna_chromosome.all.intron.bed -wa > sample.bsmap.intron.bed
intersectBed -a ${id}.bsmap.bed -b Homo_sapiens.GRCh38.dna_chromosome.all.five_prime_utr.sorted.merged.bed -wa > sample.bsmap.five_prime_utr.bed
intersectBed -a ${id}.bsmap.bed -b Homo_sapiens.GRCh38.dna_chromosome.all.three_prime_utr.sorted.merged.bed -wa > sample.bsmap.three_prime_utr.bed
intersectBed -a ${id}.bsmap.bed -b Homo_sapiens.GRCh38.dna_chromosome.all.intergenic.bed -wa > sample.bsmap.intergenic.bed
intersectBed -a ${id}.bsmap.bed -b Homo_sapiens.GRCh38.dna_chromosome.all.promotor.sorted.merged.bed -wa > sample.bsmap_promotor.bed

# Identify DMRs or DhMRs
metilene_input.pl --in1 group1_1_meth.BSmap_dep10.bedGraph,group1_2_meth.BSmap_dep10.bedGraph,group1_3_meth.BSmap_dep10.bedGraph,group1_4_meth.BSmap_dep10.bedGraph \
                  --in2 group2_1_meth.BSmap_dep10.bedGraph,group2_2_meth.BSmap_dep10.bedGraph,group2_3_meth.BSmap_dep10.bedGraph,group2_4_meth.BSmap_dep10.bedGraph \
                                  --h1 group1 --h2 group2 --out metilene_group1_group2.input
metilene -a group1 -b group2 -d 0.2 -t 50 metilene_group1_group2.input > metilene_group1_group2.output