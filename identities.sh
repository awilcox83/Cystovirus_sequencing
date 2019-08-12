BlastN all sequences against each other





for query in *.fasta
	for subject in *.fasta
		blastn -query query -subject subject > output.txt

