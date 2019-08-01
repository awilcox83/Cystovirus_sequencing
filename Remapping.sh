#Mapping reads back to assembled genomes to see average coverage.

for file in *_R1_001.fastq
do
	prefix=${file%_R1_001.fastq}
	trimmomatic PE -phred33  ${prefix}_R1_001.fastq ${prefix}_R2_001.fastq ${prefix}_R1.trimmed.fastq ${prefix}_R1.unpaired.trimmed.fastq \
	${prefix}_R2.trimmed.fastq ${prefix}_R2.unpaired.trimmed.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
done

for file in *_R1.trimmed.fastq
do
	prefix=${file%_R1.trimmed.fastq}
	sickle pe -f ${prefix}_R1.trimmed.fastq -r ${prefix}_R2.trimmed.fastq -t sanger -o ${prefix}_R1.qc.fastq -p ${prefix}_R2.qc.fastq -s ${prefix}.qcsingles.fastq -n > sickle_${prefix}.txt
done

for file in *.unpaired.trimmed.fastq
do
	prefix=${file%.unpaired.trimmed.fastq}
	sickle se -f ${prefix}.unpaired.trimmed.fastq -t sanger -o ${prefix}.unpaired.qc.fastq  -n > sickle_${prefix}.txt
done

# build references


for file in *.fasta
do
	prefix=${file%.fasta}
	bowtie2-build ${prefix}.fasta ${prefix}_ref
done


for file in *R1.qc.fastq
do
	prefix=${file%R1.qc.fastq}
	IFS='_' read -r sample string <<< "$prefix"
	echo ${sample}
	bowtie2 -x ${sample}_ref -U ${prefix}R1.unpaired.qc.fastq -U ${prefix}R2.unpaired.qc.fastq -1 ${prefix}R1.qc.fastq -2 ${prefix}R2.qc.fastq -S ${prefix}.mapped.SAM
done	

for file in *.mapped.SAM 
do
	prefix=${file%.mapped.SAM}
	samtools view -bS ${prefix}.mapped.SAM > ${prefix}.mapped.BAM
	samtools sort ${prefix}.mapped.BAM -o ${prefix}.mapped.sorted.BAM
	samtools index ${prefix}.mapped.sorted.BAM
	bam-readcount -b 20 -w 1 ${prefix}.mapped.sorted.BAM > ${prefix}.tab
	samtools depth ${prefix}.mapped.sorted.BAM | awk '{sum+=$3} END { print "Average = ",sum/NR}'
done


###  average coverage
samtools depth 28i.mapped.sorted.BAM | awk '{sum+=$3} END { print "Average = ",sum/NR}'
