#This script takes all fastq files that did not assemble and maps them to all assemblies to find the closes match.


#Index all fastas for mapping


for file in *.fasta
do
	prefix=${file%.fasta}
	bowtie2-build ${prefix}.fasta ${prefix}_ref
done

#Trimming and quality for reads

cd unassembled_fastqs
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

#For each pair, try and map to all references

for file in *R1.qc.fastq
do
    prefix=${file%R1.qc.fastq}
    #We are appending all bowtie output to a single file, so delete if it exists
    rm bowtie_${prefix}.txt
    touch bowtie_${prefix}.txt
    for xxx in *.fasta
    do
        y=${xxx%.fasta}
        echo $y >> bowtie_${prefix}.txt
        bowtie2 -x ${y}_ref -U ${prefix}R1.unpaired.qc.fastq -U ${prefix}R2.unpaired.qc.fastq -1 ${prefix}R1.qc.fastq -2 ${prefix}R2.qc.fastq -S ${prefix}.mapped.SAM 2>> bowtie_${prefix}.txt
        printf "\n\n\n" >> bowtie_${prefix}.txt
    done
done
