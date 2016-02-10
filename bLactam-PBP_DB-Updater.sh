#!/bin/bash -l

. /usr/share/Modules/init/bash
module load Python/2.7
module load ncbi-blast+/2.2.29

#Comment blah...#

while getopts :r:u:t: option
do
    case $option in
        r) refSeq_path=$OPTARG;;
        u) new_PBPseq=$OPTARG;;
        t) typed_PBPs=$OPTARG;;
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

if [[ -z "$new_PBPseq" ]]
then
    echo "A file containing the new PBP sequences hasn't been given."
else
    if [[ -e "$new_PBPseq" ]]
    then
        newPBPs="${new_PBPseq}"
    else
        echo "The file containing the new PBP sequences doesn't exist."
    fi
fi

if [[ -z "$typed_PBPs" ]]
then
    echo "A file containing the typed PBP sequences hasn't been given."
else
    if [[ -e "$typed_PBPs" ]]
    then
        typedPBPs="${typed_PBPs}"
    else
        echo "The file containing the typed PBP sequences doesn't exist."
    fi
fi

eval refSeq=$refSeq
eval pbp1A=$pbp1A
eval pbp2B=$pbp2B
eval pbp2X=$pbp2X
sampl_out=$(echo "$typedPBPs" | sed 's/TABLE/SAMPL/g')

blastTyper () {
    pbpName=$(echo "$2" | sed 's/.*bLactam_\(.*\)\.faa/\1/g')
    SampleName=$(cat "$1" | grep '>' | sed 's/>\(.*\)-[1|2][A|B|X]-S2::.*|.*::/\1/g')
    typedPath=$(dirname "$5")
    pbpAlleleID=""
    #sampl_out=$(echo "$5" | sed 's/TABLE/SAMPL/g')

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
            blastp -db "$refSeq"/"$pbpName"_prot_blast_db -query "$1" -outfmt 6 -out "$SampleName"_blast-out_"$pbpName".txt
        fi

        cat "$SampleName"_blast-out_"$pbpName".txt | sort -r -n -k12,12 -k3,3 -k4,4 | sed -n 1p > prot_best_blast.txt
	printf "\n"
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
            pbpAlleleID=$(echo "$blastName" | sed 's/^\([0-9]\+\)||.*/\1/g')
        else
            echo "Didn't find match.  Will add sequence to database."
            ###Modify header (including adding number ID)###
            oldMaxID=$(cat "$2" | grep '>' | sed 's/>\([0-9]\+\)||.*/\1/g' | sort -n | tail -n1)
            pbpAlleleID=$(echo "$oldMaxID + 1" | bc)

            while read line
            do
                if [[ $line =~ ^\>.* ]]
                then
                    echo $line | sed "s/>\(.*\)/>$pbpAlleleID||\1/g" >> Single_newRef_Seq.faa
                else
                    echo $line >> Single_newRef_Seq.faa
                    #continue
                fi
            done < "$1"

	    printf "New PBP Sequence\n"
	    cat Single_newRef_Seq.faa
            ###Output the new allele to Final_newRef_Seq.faa and add it to the database###
            #cat Single_newRef_Seq.faa >> "$typedPath"/Final_newRef_Seq.faa
            cat Single_newRef_Seq.faa >> "$2"
            echo "Remaking blast database"
            rm "$refSeq"/"$pbpName"_prot_blast_db*
            makeblastdb -in "$2" -dbtype prot -out "$refSeq"/"$pbpName"_prot_blast_db
            rm Single_newRef_Seq.faa
        fi

       ###Add in the code to update the TABLE pbp output file with the pbp IDs for the new alleles###
       while read line
       do
           lineName=$(echo "$line" | awk -F"\t" '{print $1}')
           if [[ "$3" == "$lineName" ]]
           then
               if [[ "$4" == "pbp1A" ]]
               then
                   index=11
                   echo "$line" | awk -v OFS="\t" '{$'$index'='$pbpAlleleID'; print }' >> "$5"_PRE
               elif [[ "$4" == "pbp2B" ]]
               then
                    index=12
                    echo "$line" | awk -v OFS="\t" '{$'$index'='$pbpAlleleID'; print }' >> "$5"_PRE
               elif [[ "$4" == "pbp2X" ]]
               then
                    index=13
                    echo "$line" | awk -v OFS="\t" '{$'$index'='$pbpAlleleID'; print }' >> "$5"_PRE
               fi
           else
               echo "$line" >> "$5"_PRE
           fi
       done < "$5"
       rm "$5"
       mv "$5"_PRE "$5"
    fi

    ###Remove Files###
    rm "$SampleName"_blast-out_"$pbpName".txt
    rm prot_best_blast.txt
    pbpGene=$(echo "$pbpName" | sed 's/pbp//g')
    pathPrefx=$(dirname "$1")
    #rm "$pathPrefx"/Final-Full_"$pbpGene"-S2*
    #rm "$pathPrefx"/Final-Frag_"$pbpGene"-S2*
}



cp "$typedPBPs" "$typedPBPs"_not-updated
while read line
do
    sampleID=$(echo "$line" | awk -F"\t" '{print $1}')
    pbpPath=$(echo "$line" | awk -F"\t" '{print $2}')
    pbpAllele=$(echo "$line" | awk -F"\t" '{print $3}')
    if [[ "$pbpAllele" == "pbp1A" ]]
    then
        blastTyper "$pbpPath" "$Ref_1A" "$sampleID" "$pbpAllele" "$typedPBPs"
    elif [[ "$pbpAllele" == "pbp2B" ]]
    then
        blastTyper "$pbpPath" "$Ref_2B" "$sampleID" "$pbpAllele" "$typedPBPs"
    elif [[ "$pbpAllele" == "pbp2X" ]]
    then
        blastTyper "$pbpPath" "$Ref_2X" "$sampleID" "$pbpAllele" "$typedPBPs"
    fi
done < "$newPBPs"

###Use the updated TABLE typing output to update the SAMPL typing output###
cp "$sampl_out" "$sampl_out"_not-updated
while read line
do
    sampleID=$(echo "$line" | awk -F"\t" '{print $1}')
    pbpPath=$(echo "$line" | awk -F"\t" '{print $2}')
    pbpAllele=$(echo "$line" | awk -F"\t" '{print $3}')

    table_val=$(cat "$typedPBPs" | grep "^$sampleID" | cut -f11-13)
    #echo "TABLE PBPs: $table_val"

    touch "$sampl_out"_PRE
    is_right_sampl="no"
    pbpCount=0
    while IFS= read -r line
    do
	line_count=${#line}
	if [[ "$line" =~ "$sampleID" ]]
	then
	    is_right_sampl="yes"
	    echo "$line" >> "$sampl_out"_PRE
	    continue
	elif [[ "$is_right_sampl" == "no" ]]
	then
	    echo "$line" >> "$sampl_out"_PRE
	    continue
	fi

	if [[ "$is_right_sampl" == "yes" ]]
	then
	    if [[ "$line" =~ "PBP_ID Code:" || "$line" =~ "PBP_1A" ]]
	    then
		pbpCount=$[pbpCount + 1]
                echo "$line" >> "$sampl_out"_PRE
            elif [[ "$pbpCount" -eq 2 ]]
            then
		printf "\t\t$table_val\n" >> "$sampl_out"_PRE
		pbpCount=$[pbpCount + 1]
	    elif [[ "$pbpCount" -ne 2 ]]
	    then
		if [[ "$line" =~ "Misc. Resistance:" ]]
		then
		     is_right_sampl="no"
		     echo "$line" >> "$sampl_out"_PRE
		else
		    echo "$line" >> "$sampl_out"_PRE
		fi
	    fi
	fi
    done < "$sampl_out"
    rm "$sampl_out"
    mv "$sampl_out"_PRE "$sampl_out"
done < "$newPBPs"

module unload Python/2.7
module unload ncbi-blast+/2.2.29
