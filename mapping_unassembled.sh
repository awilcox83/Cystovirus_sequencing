#This script takes all fastq files that did not assemble and maps them to all assemblies to find the closes match.


#Index all fastas for mapping


for file in *.fasta
do
	prefix=${file%.fasta}
	bowtie2-build ${prefix}.fasta ${prefix}_ref
done
