#!/bin/bash -l

. /usr/share/Modules/init/bash  
module load Python/2.7
module load prodigal/2.60
module load ncbi-blast+/2.2.29
module load BEDTools/2.17.0
module load perl/5.12.3
#export PATH=/scicomp/home/ycm6/bjm_bin/:$PATH              ###Remove This Line###
#export PATH=/scicomp/home/ycm6/cutadapt-1.4.2/bin/:$PATH   ###Remove This Line###

###This program will extract a set of loci defined in a query mutli-fasta file given in the '-q' option from a set of paired-end fastq files given by the '-1' and '-2' options.###
###The query file doesn't need to contain the complete gene sequence. If a fragment of the gene seq is provided than the script will, by default, extract the whole gene.###
###However, if the '-p' flag is provided than the code will also extract just the matching fragment.###
###The percent coverage and percent identity thresholds for a true blast hit are given by the '-L' and '-i' arguments.  The default value is 70% identity and 25% coverage.###
###The genome size is needed for proper genome assembly. It is given by the '-s' argument with a default value of 2.1M.###

while getopts :1:2:q:o:n:L:i:s:p option
do
    case $option in
        1) fastq1_path=$OPTARG;;
        2) fastq2_path=$OPTARG;;
        q) query_file=$OPTARG;;
        o) output_dir=$OPTARG;;
        p) geneFrag="TRUE";;
	L) perLength=$OPTARG;;
	i) perIdent=$OPTARG;;
	s) genomeSize=$OPTARG;;
	n) outName=$OPTARG;;
    esac
done

if [[ -z "$fastq1_path" ]]
then
    echo "No paired end 1 fastq file path argument given."
    exit 1
fi

if [[ -z "$fastq2_path" ]]
then
    echo "No paired end 2 fastq file path argument given."
    exit 1
fi

if [[ -z "$query_file" ]]
then
    echo "The query fasta file path argument has not been given."
    exit 1
fi

if [[ -z "$output_dir" ]]
then
    echo "The files will be output into the current directory."
    out_dir="./"
elif [[ ! -d "$output_dir" ]]
then
    mkdir "$output_dir"
    out_dir="$output_dir"
    echo "The output directory has been created: $out_dir"
else
    out_dir="$output_dir"
fi

if [[ ! -z "$perLength" ]]
then
    fracLength=$(bc <<< "scale=2;$perLength/100")
    echo "The minimum query matching length is: $perLength%"
else
    fracLength=0.25
    echo "The minimum query matching length is the default value of 25%"
fi

if [[ ! -z "$perIdent" ]]
then
    echo "The minimum query matching identity is: $perIdent%"
else
    perIdent=70
    echo "The minimum query matching identity is the default value of 70%"
fi

if [[ ! -z "$genomeSize" ]]
then
    if [[ "$genomeSize" =~ [0-9]+[k|M|G]$ || "$genomeSize" =~ [0-9]+$ ]]
    then
	echo "The given genome size is: $genomeSize"
    else
	echo "The genome size is not in the proper format."
        echo "Please provide the genome size in bp (can use k/M/G suffix)."
	exit
    fi
else
    echo "The genome size argument has not been given."
    echo "The program will use a default value of 2.1M."
    genomeSize=2.1M
fi

if [[ -e "$fastq1_path" ]]
then
    fastq1="${fastq1_path}"
    echo "Paired-end Read-1 is: $fastq1"
else
    echo "This sequence directory is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/fastq_file)."
    exit 1
fi

if [[ -e "$fastq2_path" ]]
then
    fastq2="${fastq2_path}"
    echo "Paired-end Read-2 is: $fastq2"
else
    echo "This sequence directory is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/fastq_file)."
    exit 1
fi

if [[ -e "$query_file" ]]
then
    query_path="${query_file}"
    echo "The gene typing database file is in the following location: $query_path"
else
    echo "This gene typing database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/query_file)."
    exit 1
fi

if [[ ! -z "$outName" ]]
then
    echo "The output file name prefix: $outName"
else
    outName=$(echo "$fastq1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')
fi

eval fastq1=$fastq1
eval fastq2=$fastq2
eval query_path=$query_path

cd "$out_dir"
#fastq_extension=$(basename "$fastq1")                       ###Remove This Line###
#fastq_name=$(echo "$fastq1" | awk -F"/" '{print $(NF)}')    ###Remove This Line###
fastq1_trimd=cutadapt_"$outName"_S1_L001_R1_001.fastq
fastq2_trimd=cutadapt_"$outName"_S1_L001_R2_001.fastq

echo "Beginning cutadapt"
if [[ -e "$fastq1_trimd" ]]
then
    echo "Fastq files have already been preprocessed"
else
    cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output temp2.fastq -o temp1.fastq "$fastq1" "$fastq2"
    cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output "$fastq1_trimd" -o "$fastq2_trimd" temp2.fastq temp1.fastq
    rm temp1.fastq temp2.fastq
fi

velvetkVal=$(velvetk.pl --best --size "$genomeSize" "$fastq1_trimd" "$fastq2_trimd")
echo "Beginning Velvet"
if [[ -d ./velvet_output ]]
then
    echo "Velvet assembly has already been completed"
else
    VelvetOptimiser.pl -s "$velvetkVal" -e "$velvetkVal" -o "-scaffolding no" -f "-shortPaired -separate -fastq $fastq1_trimd $fastq2_trimd" -d velvet_output
fi

echo "Beginning Prodigal"
if [[ -e prodigal_gene_file ]]
then
    echo "Gene prediction has already been completed"
else
    prodigal -c -f gff -i ./velvet_output/contigs.fa -a PRE_"$outName".faa -o prodigal_"$outName".gff -d PRE_"$outName".fasta
    cat PRE_"$outName".faa | sed 's/ # .*//g' > prodigal_"$outName".faa
    cat PRE_"$outName".fasta | sed 's/ # .*//g' > prodigal_"$outName".fna
    rm PRE_"$outName".faa
    rm PRE_"$outName".fasta
fi

echo "Beginning Blast"
###Create a blast nucleotide database using the predicted genes obtained from Prodigal.###
if [[ -e extract_nucl_blast_db.nhr ]]
then
    echo "Blast database has already been created"
else
    makeblastdb -in prodigal_"$outName".fna -dbtype nucl -out extract_nucl_blast_db
fi

###Blast each sequence given in the query fasta file against the blast nucleotide database.###
cat "$query_path" | grep '>' | sed 's/>//g' > fasta_name.txt
while read geneLine
do
    extractFastaByID.pl $geneLine < "$query_path" > Extract_query.txt
    blastn -db extract_nucl_blast_db -query Extract_query.txt -outfmt 6 -word_size 7 -out nucl_"$geneLine"_blast.txt

    ###Get the best blast hit by sorting the blast output by bit score, then % ID, then alignment length and select the first hit as the best match.### 
    cat nucl_"$geneLine"_blast.txt | sort -r -n -k12,12 -k3,3 -k 4,4 | sed -n 1p > nucl_best_blast.txt
    #cat prot_"$geneLine"_blast.txt | sort -r -n -k12,12 | sed -n 1p > prot_best_blast.txt           ###Remove This Line###

    ###Caclulate the length threshold for a true query match. The script requires a potential match align at least 25% of the query length.###
    geneLen=$(cat Extract_query.txt | grep -v '>' | tr '\n' ' ' | sed 's/ //g' | wc -c)
    matchLen=$(echo "$fracLength * $geneLen" | bc)
    ###Obtain the identity, length and name of the best blast hit.###
    blastIden=$(cat nucl_best_blast.txt | awk -F"\t" '{print $3}')
    blastLen=$(cat nucl_best_blast.txt | awk -F"\t" '{print $4}')
    blastName=$(cat nucl_best_blast.txt | awk -F"\t" '{print $2}')

    echo "Query sequence name is $geneLine"
    echo "Query Gene Length is $geneLen"
    printf "\n"
    echo "Best blast hit name is $blastName"
    echo "Best blast hit identity is $blastIden"
    echo "Best blast hit length is $blastLen"
    printf "\n"
    echo "Min. Matched identity threshold is $perIdent"
    echo "Min. Matched length threshold is $matchLen"

    ###The script will call a match if the best blast hit aligns above a user provided identity (default: 70%) to the query 
    ###and across a user provided minimum percentage (default: 25%) of the query length.###
    if [[ $(echo "$blastIden >= $perIdent" | bc) -eq 1 && $(echo "$blastLen >= $matchLen" | bc) -eq 1 ]]
    then
        extractFastaByID.pl $blastName < prodigal_"$outName".fna > Final-Full_"$geneLine"_nucl.fna
        extractFastaByID.pl $blastName < prodigal_"$outName".faa > Final-Full_"$geneLine"_prot.faa
    else
        echo "For $geneLine:" | tee -a gene_extraction_ERROR_log.txt
        printf "The best blast hit ($blastName) for $geneLine didn't meet minimum criteria of length and identity to call a true match\n" | tee -a gene_extraction_ERROR_log.txt
	printf "\n" | tee -a gene_extraction_ERROR_log.txt
        continue
    fi

    ###If the '-p' flag is called, then the script will also extract just the section of the target sequence that corresponds to the query fragment.###
    ###If the best blast hit didn't include the entire query fragment, then the code will calculate the expected start/end coordinates of the complete fragment 
    ###and will attempt to extract the full fragment from the predicted gene sequence.###  
    ###The number of non-aligning bases at each end of the matching target sequence will recorded in the header name.###
    if [[ "$geneFrag" ]]
    then
        echo "Extracting gene fragment"
	blast_qStrt=$(cat nucl_best_blast.txt | awk -F"\t" '{print $7}')
	blast_qEnd=$(cat nucl_best_blast.txt | awk -F"\t" '{print $8}')
	blast_tStrt=$(cat nucl_best_blast.txt | awk -F"\t" '{print $9}')
	blast_tEnd=$(cat nucl_best_blast.txt | awk -F"\t" '{print $10}')
	blast_Strt_Diff=$(echo "$blast_qStrt - 1" | bc)
	blast_End_Diff=$(echo "$geneLen - $blast_qEnd" | bc)

        nucl_grab_Strt=$(echo "$blast_tStrt - $blast_qStrt" | bc)
        nucl_grab_End=$(echo "$blast_tEnd + $blast_End_Diff" | bc)

        prot_grab_Strt=$(echo "$nucl_grab_Strt / 3" | bc)
        prot_grab_End=$(echo "($nucl_grab_End / 3) + 1" | bc)

        echo "The target fragment sequence start is $nucl_grab_Strt"
        echo "The target fragment sequence end is $nucl_grab_End"
        #echo "Blast Start Diff is $blast_Strt_Diff &&&& Blast End Diff is $blast_End_Diff"  ###Remove This Line###

	###Create a fasta header name for the extracted sequence###
	headr_Append="$outName"-"$geneLine::S$blast_Strt_Diff|E$blast_End_Diff::"

        ###Use BEDtools with these coordinates to grab the sequence###
	###If BEDtools is unable to extract the sequence than the fragment sequence coordinates fell outside of the range of the gene sequence coordinates.###
	###In this case the script will throw an error.###
        printf "$blastName\t$nucl_grab_Strt\t$nucl_grab_End\n" > nucl_temp.bed
        bedtools getfasta -fi Final-Full_"$geneLine"_nucl.fna -bed nucl_temp.bed -fo PRE-Frag_"$geneLine"_nucl.txt
	if [[ -s PRE-Frag_"$geneLine"_nucl.txt ]]
	then
	    cat PRE-Frag_"$geneLine"_nucl.txt | sed 's/>.*/>'"$headr_Append"'/g' > Final-Frag_"$geneLine"_nucl.fna
        else
	    echo "For $geneLine:" | tee -a gene_extraction_ERROR_log.txt
	    echo "One or both genomic coordinates of the fragment fall outside the coordinates of the gene." | tee -a gene_extraction_ERROR_log.txt
	    echo "The extracted predicted gene nucleotide and protein sequence will be renamed to: Final-Full_$geneLine-extract-Anomaly_nucl.fna" | tee -a gene_extraction_ERROR_log.txt
	    echo "\n" | tee -a gene_extraction_ERROR_log.txt
	    mv Final-Full_"$geneLine"_nucl.fna Final-Full_"$geneLine"_extract-Anomaly_nucl.fna
	fi
        printf "$blastName\t$prot_grab_Strt\t$prot_grab_End\n" > prot_temp.bed
        bedtools getfasta -fi Final-Full_"$geneLine"_prot.faa -bed prot_temp.bed -fo PRE-Frag_"$geneLine"_prot.txt
	if [[ -s PRE-Frag_"$geneLine"_prot.txt ]]
        then
	    cat PRE-Frag_"$geneLine"_prot.txt | sed 's/>.*/>'"$headr_Append"'/g' > Final-Frag_"$geneLine"_prot.faa
	else
            #echo "One or both genomic coordinates of the fragment fall outside the coordinates of the gene."                      ###Remove This Line###
            #echo "This may be caused by an error in assembly or gene prediction so the full predicted gene will be renamed blah." ###Remove This Line###
            mv Final-Full_"$geneLine"_prot.faa Final-Full_"$geneLine"_extract-Anomaly_prot.faa
        fi
	rm PRE-Frag_"$geneLine"_nucl.txt
	rm PRE-Frag_"$geneLine"_prot.txt
	printf "\n"
    fi
done < fasta_name.txt

rm fasta_name.txt
rm nucl_temp.bed
rm prot_temp.bed
rm Extract_query.txt
rm nucl_best_blast.txt

module unload Python/2.7
module unload prodigal/2.60
module unload ncbi-blast+/2.2.29
module unload BEDTools/2.17.0
module unload perl/5.12.3
