#BlastX all sequences against protein fasta references


for file in *.fasta
do
    querypath=${file}
    queryname=${file%.fasta}
    echo ${queryname}
    for file2 in ./proteins/*.fasta
    do
        subjectpath=${file2}
        subjectname=$(basename "$file2")
        subjectname=${subjectname%.fasta}
        echo ${subjectname}
        blastx -query ${querypath} -subject ${subjectpath} -outfmt 6 > ./palignments/${queryname}_vs_${subjectname}.txt
    done
done
