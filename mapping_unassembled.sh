#This script takes all fastq files that did not assemble and maps them to all assemblies to find the closest match.


#Index all fastas for mapping


for file in *.fasta
do
	prefix=${file%.fasta}
	bowtie2-build ${prefix}.fasta ${prefix}_ref
done

#Trimming and quality for reads

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


rm phage_mapping.tsv
touch phage_mapping.tsv
printf "sample\t105b\t105c\t105d\t120a\t120c\t28d\t28i\t42c\t42d\t42e\t42f\t42g\t42k\t42m\t64505\t64506\t64507\t64511\t85a\t85c\t85d_phi2954\t85d_phi6\t90a\t90b\tca51\tca58\tca61\tca62\tca64\tca64d\tca65\tca65c_1\tca65d\tca66\tca66d\tca68c\tca69\tca73\tca76\tkri289\tkri301\tphi10\tphi14\tphi6\tphi9\tpt106\tpt110\tv1-1\tv1-3\tv1-7\tv1_4_phi2954\tv1_4_phi6\tv2\tv3\n" > phage_mapping.tsv

# Create a table showing how well each sample maps to each reference genome/assembly 

for file in bowtie*.txt
do
   #Get sample name from filename
	A="$(cut -d'_' -f2 <<<$file)"
	printf $A >> phage_mapping.tsv
	printf "\t" >> phage_mapping.tsv
   #get % from all lines containing "overall" 
    cat $file | while read line
    do
        if [[ $line == *"overall"* ]]
        then
            printf  $line | cut -d' ' -f1 >> phage_mapping.tsv
            #errors because printf treats % as a special character, but it doesn't affect the output
            printf "\t" >> phage_mapping.tsv
        fi
    done
    printf "\n" >> phage_mapping.tsv
done
