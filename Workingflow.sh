# This is GATK workflow for resequence data
# Map to reference genome using BWA 
for id in $(cat id.txt)
do 
  bwa mem -t 24 -M -P -R '@RG\tID:$id\tSM:$id\tLB:$id\tPL:Illumina'  Reference.genome.fasta ${id}_1.clean.fq.gz ${id}_2.clean.fq.gz | samtools sort -@ 12 -o bwa_${id}.sort.bam
done
# flagstat to get mapping results
for id in $(cat id.txt) 
  do samtools flagstat bwa_${id}.sort.bam > bwa_${id}.flagstat.tax
done
# Call depth and coverage
for id in $(cat id.txt) 
do samtools coverage bwa_${id}.sort.bam  > bwa_${id}.sort.dep_cov
done
# MarkDuplicates
for id in $(cat id.txt)
  do java -Xmx4g -jar gatk-package-4.5.0.0-local.jar -I bwa_${id}.sort.bam -O bwa_${id}.sort.markdup.bam  -M bwa_${id}.sort.markdup_metrics.txt
done
# HaplotypeCaller with different Chrs
for id in $(cat id.txt)
do
  samtools index  bwa_${id}.sort.markdup.bam  
  java -Xmx4g -jar gatk-package-4.5.0.0-local.jar HaplotypeCaller -R Reference.genome.fasta -I bwa_${id}.sort.markdup.bam   -O ${id}.chr1.g.vcf.gz  -ERC GVCF -L chr1 &
  java -Xmx4g -jar gatk-package-4.5.0.0-local.jar HaplotypeCaller -R Reference.genome.fasta -I bwa_${id}.sort.markdup.bam   -O ${id}.chr2.g.vcf.gz  -ERC GVCF -L chr2 &
 done
# CombineGVCFs chrs to individuals to population
for id in $(cat id.txt)
do
java -Xmx4g -jar gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CombineGVCFs -R  Reference.genome.fasta --variant ${id}.chr1.g.vcf.gz  --variant ${id}.chr2.g.vcf.gz  -O ${id}.g.vcf.gz 
done
java -Xmx4g -jar gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CombineGVCFs -R  Reference.genome.fasta --variant individuals1.g.vcf.gz --variant individuals2.g.vcf.gz -O  POP1.g.vcf.gz
# JointGenotyping
Java -Xmx4g -jar gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R Reference.genome.fasta -V POP1.g.vcf.gz -O POP1.genotype.vcf
# SelectVariants
java -Xmx4g -jar gatk-package-4.5.0.0-local.jar SelectVariants -R  Reference.genome.fasta -V POP1.genotype.vcf --select-type SNP -O POP1.SNP.raw.vcf 
# HardFiltered
java -Xmx4g -jar gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar VariantFiltration -R Reference.genome.fasta -V POP1.genotype.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP < 8.0 " --filter-name 'SNP_filter' -O POP1.SNP.raw.filtered.vcf 
#Get filtered SNPs
java -Xmx4g -jar gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants -R Reference.genome.fasta -V POP1.SNP.raw.filtered.vcf --exclude-filtered -O POP1.SNP.raw.filtered.pass.vcf 
# VCFtools to filter
vcftools --vcf POP1.SNP.raw.filtered.pass.vcf --maf 0.05 --max-missing 0.8 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out POP1
