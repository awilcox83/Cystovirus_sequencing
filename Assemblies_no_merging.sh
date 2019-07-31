#redo spades without merging

#cutadapt + sickle

for file in *_R1_001.fastq
do
	prefix=${file%_R1_001.fastq}
	cutadapt ${prefix}_R1_001.fastq ${prefix}_R2_001.fastq -m 25 -n 2 -o ${prefix}_R1.trimmed.fastq -p ${prefix}_R2.trimmed.fastq \
	-b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
	-a CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG -a CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG \
	-a CAAGCAGAAGACGGCATACGAGATCTAGTACGGTCTCGTGGGCTCGG -a CAAGCAGAAGACGGCATACGAGATTTCTGCCTGTCTCGTGGGCTCGG \
	-a CAAGCAGAAGACGGCATACGAGATGCTCAGGAGTCTCGTGGGCTCGG -a CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTCTCGTGGGCTCGG \
	-a CAAGCAGAAGACGGCATACGAGATCATGCCTAGTCTCGTGGGCTCGG -a AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC \
	-a AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC > ${prefix}.txt
done

for file in *_R1.trimmed.fastq.gz
do
	prefix=${file%_R1.trimmed.fastq.gz}
	sickle pe -f ${prefix}_R1.trimmed.fastq.gz -r ${prefix}_R2.trimmed.fastq.gz -t sanger -o ${prefix}_R1.qc.fastq.gz -p ${prefix}_R2.qc.fastq.gz -s ${prefix}.qcsingles.fastq.gz -n > sickle_${prefix}.txt
done
