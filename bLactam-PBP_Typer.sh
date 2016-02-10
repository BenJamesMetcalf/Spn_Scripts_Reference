#!/bin/bash -l

. /usr/share/Modules/init/bash
module load Python/2.7
module load ncbi-blast+/2.2.29

#This program takes in a 1) fasta amino acid sequence file

while getopts :a:b:x:r:o:n: option
do
    case $option in
	a) query_pbp1A=$OPTARG;;
	b) query_pbp2B=$OPTARG;;
	x) query_pbp2X=$OPTARG;;
        r) refSeq_path=$OPTARG;;
	o) output_dir=$OPTARG;;
	n) out_Name=$OPTARG;;
    esac
done

if [[ -z "$refSeq_path" ]]
then
    echo "No PBP reference database path argument given."
    exit 1
else
    refSeq="${refSeq_path}"
fi

Ref_1A=$(ls $refSeq/*_pbp1A.faa)
Ref_2B=$(ls $refSeq/*_pbp2B.faa)
Ref_2X=$(ls $refSeq/*_pbp2X.faa)

if [[ -e "$Ref_1A" && -e "$Ref_2B" && -e "$Ref_2X" ]]
then
    echo "The PBP 1A reference database is: $Ref_1A"
    echo "The PBP 2B reference database is: $Ref_2B"
    echo "The PBP 2X reference database is: $Ref_2X"
else
    echo "At least one of the PBP reference databases are not in the correct format or don't exist."
    echo "Make sure they are in the $refSeq directory and have the correct naming format." 
    echo "More specifically, each PBP sequence must end in either '_pbp1A.fasta', '_pbp2B.fasta' or '_pbp2X.fasta'"
    exit 1
fi

if [[ -z "$query_pbp1A" || -z "$query_pbp2B" || -z "$query_pbp2X" ]]
then
    echo "At least one of the PBP extracted sequence arguments have not been given"
    exit 1
else
    pbp1A="${query_pbp1A}"
    pbp2B="${query_pbp2B}"
    pbp2X="${query_pbp2X}"
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

if [[ -z "$out_Name" ]]
then 
    echo "Please provide an output prefix name using the -n argument."
    exit 1
else
    outName=$out_Name
fi  

eval refSeq=$refSeq
eval pbp1A=$pbp1A
eval pbp2B=$pbp2B
eval pbp2X=$pbp2X

cd "$out_dir"
#fastq_extension=$(basename "$fastq1")

#touch pbpID_output.txt
blastTyper () {
    pbpName=$(echo "$2" | sed 's/.*bLactam_\(.*\)\.faa/\1/g')
    SampleName=$(cat "$1" | grep '>' | sed 's/>\(.*\)-[1|2][A|B|X]-S2::.*|.*::/\1/g')

    pbpStop="no"
    if [[ ! -f "$1" ]]
    then
	printf "The following extracted PBP path argument isn't in the correct format or doesn't exist:\n"
	printf "$1\n"
	printf "$SampleName\t$pbpName\tNF\n" >> pbpID_output.txt
	pbpStop="yes"
    fi

    isBlastDB=$(ls "$refSeq"/"$pbpName"_prot_blast_db*)
    if [[ "$pbpStop" == "no" ]]
    then
	if [[ -z "$isBlastDB" ]]
	then
	    ### Make blast database
	    makeblastdb -in "$2" -dbtype prot -out "$refSeq"/"$pbpName"_prot_blast_db
	    blastp -db "$refSeq"/"$pbpName"_prot_blast_db -query "$1" -outfmt 6 -out "$SampleName"_blast-out_"$pbpName".txt
	else
	    #outName=$(echo "$1" | sed 's/\..*//g')
	    blastp -db "$refSeq"/"$pbpName"_prot_blast_db -query "$1" -outfmt 6 -out "$SampleName"_blast-out_"$pbpName".txt
	fi
    
	cat "$SampleName"_blast-out_"$pbpName".txt | sort -r -n -k12,12 -k3,3 -k4,4 | sed -n 1p > prot_best_blast.txt
	cat prot_best_blast.txt
	geneLen=$(cat "$1" | grep -v '>' | tr '\n' ' ' | sed 's/ //g' | wc -c)
	blastIden=$(cat prot_best_blast.txt | awk -F"\t" '{print $3}')
	blastLen=$(cat prot_best_blast.txt | awk -F"\t" '{print $4}')
	blastName=$(cat prot_best_blast.txt | awk -F"\t" '{print $2}')

	#echo "Best blast identity is: $blastIden"
	#echo "Best blast length is: $blastLen"
	#echo "PBP gene length is: $geneLen"

	if [[ $(echo "$blastIden == 100" | bc) -eq 1 && $(echo "$blastLen == $geneLen" | bc) -eq 1 ]]
	then
	    echo "Found a match"
	    pbpID=$(echo "$blastName" | sed 's/^\([0-9]\+\)||.*/\1/g')
	    printf "$outName\t$pbpName\t$pbpID\n" >> pbpID_output.txt
	    pbpOnlyID=$(echo "$2" | sed 's/.*bLactam_pbp\(.*\)\.faa/\1/g')
	    echo "The ID value is: $pbpOnlyID"
	    #rm Final-Full_"$pbpOnlyID"-S2*
	    #rm Final-Frag_"$pbpOnlyID"-S2*
	else
	    echo "Didn't find match.  The sequence needs to be added to the database using 'bLactam-PBP_Updater.sh'"
            fragPath=$(readlink -e "$1")
	    printf "$outName\t$fragPath\t$pbpName\n" >> Final_newRef_updater.txt
            printf "$outName\t$pbpName\tNA\n" >> pbpID_output.txt
	fi
    fi
}

blastTyper "$pbp1A" "$Ref_1A"
blastTyper "$pbp2B" "$Ref_2B"
blastTyper "$pbp2X" "$Ref_2X"

touch Final_pbpID_output.txt
printf "Sample_Name\tPBP_1A\tPBP_2B\tPBP_2X\tPBP_Code\n" > Final_pbpID_output.txt
#outName=$(cat pbpID_output.txt | head -n1 | awk -F"\t" '{print $1}')
Code1A=$(cat pbpID_output.txt | head -n1 | awk -F"\t" '{print $3}')
Code2B=$(cat pbpID_output.txt | head -n2 | tail -n1 | awk -F"\t" '{print $3}')
Code2X=$(cat pbpID_output.txt | head -n3 | tail -n1 | awk -F"\t" '{print $3}')
printf "$outName\t$Code1A\t$Code2B\t$Code2X\t$Code1A:$Code2B:$Code2X\n" >> Final_pbpID_output.txt

module unload Python/2.7
module unload ncbi-blast+/2.2.29
