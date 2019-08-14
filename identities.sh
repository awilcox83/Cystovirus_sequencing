#BlastN all sequences against each other


for file in *.fasta
do
    querypath=${file}
    queryname=${file%.fasta}
    echo ${queryname}
    for file2 in ./complete_genomes/*.fasta
    do
        subjectpath=${file2}
        subjectname=$(basename "$file2")
        subjectname=${subjectname%.fasta}
        echo ${subjectname}
        tblastx -query ${querypath} -subject ${subjectpath} -outfmt "6 qseqid length pident" > ./palignments/${queryname}_vs_${subjectname}.txt
    done
done




rm phi6_identity.tsv
touch phi6_identity.tsv
printf "sample\tL_length\tL_match_percentage\tM_length\tM_match_percentage\tS_length\tS_match_percentage\n" > phi6_identity.tsv
for file in *_vs_phi6.txt
do
    sample=${file%_vs_phi6.txt}
    printf ${sample} >> phi6_identity.tsv
    printf "\t" >> phi6_identity.tsv
    firstword=""
    count=0
    cat $file | while read line
    do
        A=$(cut -d$'\t' -f 1 <<<"$line")
        length=$(cut -d$'_' -f 4 <<<"$A")        
        matches=$(cut -d$'\t' -f 2 <<<"$line")
        match_perc=$(cut -d$'\t' -f 3 <<<"$line")
        if [ "$A" != "$firstword" ]
        then
            printf $matches >> phi6_identity.tsv
            printf "\t" >> phi6_identity.tsv
            printf $match_perc >> phi6_identity.tsv
            printf "\t" >> phi6_identity.tsv
        fi
        firstword=$A
        echo $match_perc
    done
    printf "\n" >> phi6_identity.tsv
done

rm phi2954_identity.tsv
touch phi2954_identity.tsv
printf "sample\tL_length\tL_match_percentage\tM_length\tM_match_percentage\tS_length\tS_match_percentage\n" > phi2954_identity.tsv
for file in *_vs_phi2954.txt
do
    sample=${file%_vs_phi2954.txt}
    printf ${sample} >> phi2954_identity.tsv
    printf "\t" >> phi2954_identity.tsv
    firstword=""
    count=0
    cat $file | while read line
    do
        A=$(cut -d$'\t' -f 1 <<<"$line")
        length=$(cut -d$'_' -f 4 <<<"$A")        
        matches=$(cut -d$'\t' -f 2 <<<"$line")
        match_perc=$(cut -d$'\t' -f 3 <<<"$line")
        if [ "$A" != "$firstword" ]
        then
            printf $matches >> phi2954_identity.tsv
            printf "\t" >> phi2954_identity.tsv
            printf $match_perc >> phi2954_identity.tsv
            printf "\t" >> phi2954_identity.tsv
        fi
        firstword=$A
        echo $match_perc
    done
    printf "\n" >> phi2954_identity.tsv
done

rm phi13_identity.tsv
touch phi13_identity.tsv
printf "sample\tL_length\tL_match_percentage\tM_length\tM_match_percentage\tS_length\tS_match_percentage\n" > phi13_identity.tsv
for file in *_vs_phi13.txt
do
    sample=${file%_vs_phi13.txt}
    printf ${sample} >> phi13_identity.tsv
    printf "\t" >> phi13_identity.tsv
    firstword=""
    count=0
    cat $file | while read line
    do
        A=$(cut -d$'\t' -f 1 <<<"$line")
        length=$(cut -d$'_' -f 4 <<<"$A")        
        matches=$(cut -d$'\t' -f 2 <<<"$line")
        match_perc=$(cut -d$'\t' -f 3 <<<"$line")
        if [ "$A" != "$firstword" ]
        then
            printf $matches >> phi13_identity.tsv
            printf "\t" >> phi13_identity.tsv
            printf $match_perc >> phi13_identity.tsv
            printf "\t" >> phi13_identity.tsv
        fi
        firstword=$A
        echo $match_perc
    done
    printf "\n" >> phi13_identity.tsv
done


rm phiNN_identity.tsv
touch phiNN_identity.tsv
printf "sample\tL_length\tL_match_percentage\tM_length\tM_match_percentage\tS_length\tS_match_percentage\n" > phiNN_identity.tsv
for file in *_vs_phiNN.txt
do
    sample=${file%_vs_phiNN.txt}
    printf ${sample} >> phiNN_identity.tsv
    printf "\t" >> phiNN_identity.tsv
    firstword=""
    count=0
    cat $file | while read line
    do
        A=$(cut -d$'\t' -f 1 <<<"$line")
        length=$(cut -d$'_' -f 4 <<<"$A")        
        matches=$(cut -d$'\t' -f 2 <<<"$line")
        match_perc=$(cut -d$'\t' -f 3 <<<"$line")
        if [ "$A" != "$firstword" ]
        then
            printf $matches >> phiNN_identity.tsv
            printf "\t" >> phiNN_identity.tsv
            printf $match_perc >> phiNN_identity.tsv
            printf "\t" >> phiNN_identity.tsv
        fi
        firstword=$A
        echo $match_perc
    done
    printf "\n" >> phiNN_identity.tsv
done
