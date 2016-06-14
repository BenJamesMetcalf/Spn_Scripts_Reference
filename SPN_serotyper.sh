#!/bin/bash -l

#export PATH=/scicomp/home/ycm6/TEMP_GBS-Typing:$PATH
temp_path=$(pwd)
export PATH=$PATH:$temp_path

###Load Modules###
. /usr/share/Modules/init/bash
module load freebayes/0.9.21
module load srst2/0.1.7

###This script is used to predict SPN serotype using SRST2.###

while getopts :1:2:r:o:n: option
do
    case $option in
        1) fastq1_path=$OPTARG;;
        2) fastq2_path=$OPTARG;;
        r) SeroT_ref=$OPTARG;;
        o) output_dir=$OPTARG;;
        n) out_name=$OPTARG;;
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

if [[ -z "$SeroT_ref" ]]
then
    echo "The serotype database path argument has not been given."
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

if [[ -e "$fastq1_path" ]]
then
    readPair_1="${fastq1_path}"
    echo "Paired-end Read-1 is: $readPair_1"
else
    echo "This sequence directory is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/fastq_file)."
    exit 1
fi

if [[ -e "$fastq2_path" ]]
then
    readPair_2="${fastq2_path}"
    echo "Paired-end Read-2 is: $readPair_2"
else
    echo "This sequence directory is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/fastq_file)."
    exit 1
fi

if [[ -e "$SeroT_ref"  ]]
then
    echo "The serotype database file is in the following location: $SeroT_ref"
else
    echo "This serotype database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/query_file)."
    exit 1
fi

eval readPair_1=$readPair_1
eval readPair_2=$readPair_2
eval SeroT_ref=$SeroT_ref

if [[ ! -z "$out_name" ]]
then
    echo "The output file name prefix: $out_name"
    out_nameTYPE="TYPE_${out_name}"
else
    out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-2)"--"$(NF-1)}')
    #out_nameTYPE=TYPE_"$out_name"
fi




###Start Doing Stuff###
cd "$out_dir"
srst2 --samtools_args "\-A" --input_pe "$readPair_1" "$readPair_2" --output "$out_name" --save_scores --min_coverage 99.9 --max_divergence 5 --gene_db "$SeroT_ref"

###mpileup the '.*_TYPE__.*.sorted.bam and create the called variants file with freebayes. Don't use vcf2fq b/c it won't call indels###
seroT_bam=$(ls SERO_*sorted.bam)
seroT_pileup=$(ls SERO_*sorted.bam | sed 's/\.bam/\.pileup/g')
seroT_vcf=$(ls SERO_*sorted.bam | sed 's/\.bam/\.vcf/g')
seroT_bai=$(ls SERO_*sorted.bam | sed 's/\.bam/\.bai/g')
samtools mpileup -f "$SeroT_ref" "$seroT_bam" > "$seroT_pileup"
samtools index "$seroT_bam" "$seroT_bai"
freebayes -q 20 -p 1 -f "$SeroT_ref" "$seroT_bam" -v "$seroT_vcf"
###Create the variant-called consensus fasta###
bgzip "$seroT_vcf"
tabix -p vcf "$seroT_vcf".gz
cat "$SeroT_ref" | vcf-consensus "$seroT_vcf".gz | sed 's/>[0-9]\+__.*__\(.*\)__[0-9]\+/>\1/g' > seroT_target_consensus.fna
geneType=$(ls *__fullgenes__*__results.txt)
tail -n +2 "$geneType" > temp_geneType.txt

easyCall () {
    while read line
    do
	aVal=$(echo "$line" | sed -r 's/^[A-Z]+//g')
	gName=$(echo "$line" | sed 's/\([A-Z]\+\).*/\1/')

	if [[ "$1" =~ "$line"-[0-9]+ && ! "$1" =~ "*" ]]
	then
	    #echo "Allele is $aVal and Gene is $gName"
            printf "$1\tidentical\t$aVal\t$2\n" >> TEMP_SeroType_Results.txt
	elif [[ "$1" =~ "$line"-[0-9]+"*" ]]
	then
            #echo "Allele is $aVal and Gene is $gName"
            printf "$1\timperfect\t$aVal\t$2\n" >> TEMP_SeroType_Results.txt
	fi
    #done < "$easyPath"/easyCall_Alleles.txt
    done < "$easyPath"
}

#eval easyPath1=$easyPath1
#easyPath=$(echo "$easyPath1" | sed 's/\(.*\)\//\1::/g' | sed 's/::.*/\//g')
easyPath=$(dirname "$SeroT_ref" | sed 's/$/\/easyCall_Alleles.txt/g')
echo "easy path is: $easyPath"

touch TEMP_SeroType_Results.txt
printf "Matched_Allele\tMatch_Type\tSerotype\tAvgDepth\n" >> TEMP_SeroType_Results.txt

while read geneLine
do
    tempGenePRE=$(echo "$geneLine" | awk -F'\t' '{ print $4 }')
    geneDiff=$(echo "$geneLine" | awk -F'\t' '{ print $7 }')
    if [[ -n "$geneDiff" ]]
    then
	tempGene="$tempGenePRE*"
    else 
	tempGene=$tempGenePRE
    fi
    geneDepth=$(echo "$geneLine" | awk -F'\t' '{ print $6 }')

    if [[ -e "$easyPath" ]]
    then
	easyCall $tempGene $geneDepth
    else
	echo "The 'easyCall_Alleles.txt' file wasn't given or couldn't be found"
	echo "Make sure this file is in the same directory as the SRST2 gene database"
    fi

    if [[ "$tempGene" =~ WZY19A-[0-9]+ || "$tempGene" =~ WZY19AVAR-[0-9]+ && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\t19A\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    elif [[ "$tempGene" =~ WZY19A-[0-9]+"*" || "$tempGene" =~ WZY19AVAR-[0-9]+"*" ]]
    then
        printf "$tempGene\timperfect\t19A_(Extract_WZY19A/Avar)\t$geneDepth\n" >> TEMP_SeroType_Results.txt
	extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
    fi 

    if [[ "$tempGene" =~ WZY19F-[0-9]+ || "$tempGene" =~ WZY19FVAR-[0-9]+ && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\t19F\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    elif [[ "$tempGene" =~ WZY19F-[0-9]+"*" || "$tempGene" =~ WZY19FVAR-[0-9]+"*" ]]
    then
        printf "$tempGene\timperfect\t19F_(Extract_WZY19A/Avar)\t$geneDepth\n" >> TEMP_SeroType_Results.txt
	extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
    fi

    if [[ "$tempGene" =~ WZY15A-[0-9]+ ]]
    then
	isWCIZ15F="no"
	while read geneLine2
	do
	    tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
	    geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
	    if [[ -n "$geneDiff2" ]]
	    then
		tempGene2="$tempGenePRE2*"
	    else
		tempGene2=$tempGenePRE2
	    fi
	    geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

   	    if [[ "$tempGene2" =~ WCIZ15F-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
		printf "$tempGene/$tempGene2\tWCIZ15F=identical\t15F\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
		isWCIZ15F="yes"
            elif [[ "$tempGene2" =~ WCIZ15F-[0-9]+"*" ]]
	    then
		printf "$tempGene/$tempGene2\tWCIZ15F=imperfect\t15A\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
		extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
		isWCIZ15F="yes"
	    fi
        done < temp_geneType.txt
	if [[ "$isWCIZ15F" == "no" ]]
        then
	    printf "$tempGene\tWCIZ15F=not_present\t15A\t$geneDepth\n" >> TEMP_SeroType_Results.txt
	fi
    fi

    if [[ "$tempGene" =~ WZY15B-[0-9]+ ]]
    then
	isWCIZ15B="no"
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCIZ15B-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIZ15B=identical\t15B\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
		isWCIZ15B="yes"
            elif [[ "$tempGene2" =~ WCIZ15B-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIZ15B=imperfect\t15C\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
		isWCIZ15B="yes"
            fi
        done < temp_geneType.txt
        if [[ "$isWCIZ15B" == "no" ]]
        then
            printf "$tempGene\tWCIZ15B=not_present\t15C\t$geneDepth\n" >> TEMP_SeroType_Results.txt
        fi
    fi

    if [[ "$tempGene" =~ TTS-[0-9]+ && ! "$tempGene" =~ "*" ]]
    then
	printf "$tempGene\tidentical\t37\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    elif [[ "$tempGene" =~ TTS-[0-9]+"*"$ ]]
    then
	printf "$tempGene\timperfect\t37\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    fi

    if [[ "$tempGene" =~ RRGA-[0-9]+ && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\tPilus-1\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    elif [[ "$tempGene" =~ RRGA-[0-9]+"*"$ ]]
    then
        printf "$tempGene\timperfect\tPilus-1\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    fi

    if [[ "$tempGene" =~ PITB-[0-9]+ && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\tPilus-2\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    elif [[ "$tempGene" =~ PITB-[0-9]+"*"$ ]]
    then
        printf "$tempGene\timperfect\tPilus-2\t$geneDepth\n" >> TEMP_SeroType_Results.txt
    fi	

    if [[ "$tempGene" =~ WCIP6AC-[0-9]+ ]]
    then
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCIN6AB-[0-9]+ && ! "$tempGene" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6AC=identical\t6A\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
	    elif [[ "$tempGene2" =~ WCIN6CD-[0-9]+ && ! "$tempGene" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6AC=identical\t6C\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
	    fi	

	    if [[ "$tempGene2" =~ WCIN6AB-[0-9]+ && "$tempGene" =~ WCIP6AC-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6AC=imperfect\t6A_(Extract_WCIP6AC)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
	    elif [[ "$tempGene2" =~ WCIN6CD-[0-9]+ && "$tempGene" =~ WCIP6AC-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6AC=imperfect\t6C_(Extract_WCIP6AC)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
	    fi
        done < temp_geneType.txt
    fi

    if [[ "$tempGene" =~ WCIP6BD-[0-9]+ ]]
    then
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCIN6AB-[0-9]+ && ! "$tempGene" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6BD=identical\t6B\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
            elif [[ "$tempGene2" =~ WCIN6CD-[0-9]+ && ! "$tempGene" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6BD=identical\t6D\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
            fi

            if [[ "$tempGene2" =~ WCIN6AB-[0-9]+ && "$tempGene" =~ WCIP6BD-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6BD=imperfect\t6B_(Extract_WCIP6BD)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
            elif [[ "$tempGene2" =~ WCIN6CD-[0-9]+ && "$tempGene" =~ WCIP6BD-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIP6CD=imperfect\t6D_(Extract_WCIP6BD)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
		extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt                
            fi
        done < temp_geneType.txt
    fi

    if [[ "$tempGene" =~ WZY7F-[0-9]+ ]]
    then
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCWD7F-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCWD7F=identical\t7F\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
            elif [[ "$tempGene2" =~ WCWD7F-[0-9]+"*" ]]
            then
 		printf "$tempGene/$tempGene2\tWCWD7F=imperfect\t7A:7F_(Extract WCWD7F)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
		extractFastaByID.pl $tempGenePRE2 < seroT_target_consensus.fna >> Extract_results.txt
 	    fi
        done < temp_geneType.txt
    fi

    if [[ "$tempGene" =~ WZY7C-[0-9]+ ]]
    then
        isWCHF7C="no"
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCHF7C-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCHF7C=identical\t7C\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                isWCIZ15B="yes"
            elif [[ "$tempGene2" =~ WCHF7C-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCHF7C=imperfect\t7C\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                isWCIZ15B="yes"
            fi
        done < temp_geneType.txt
        if [[ "$isWCHF7C" == "no" ]]
        then
            printf "$tempGene\tWCHF7C=not_present\t7B\t$geneDepth\n" >> TEMP_SeroType_Results.txt
        fi
    fi

    if [[ "$tempGene" =~ WZY9N-[0-9]+ ]]
    then
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCJA9N-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCJA9N=identical\t9N\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
            elif [[ "$tempGene2" =~ WCJA9N-[0-9]+"*" ]]
            then
 		printf "$tempGene/$tempGene2\tWCJA9N=imperfect\t9N_(Extract_WCJA9N)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
		extractFastaByID.pl $tempGenePRE2 < seroT_target_consensus.fna >> Extract_results.txt
 	    fi

 	    if [[ "$tempGene2" =~ WCJA9L-[0-9]+ && ! "$tempGene2" =~ "*" ]]
 	    then
 		printf "$tempGene/$tempGene2\tWCJA9L=identical\t9L\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
            elif [[ "$tempGene2" =~ WCJA9L-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCJA9L=imperfect\t9L_(Extract_WCJA9L)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
		extractFastaByID.pl $tempGenePRE2 < seroT_target_consensus.fna >> Extract_results.txt
            fi
        done < temp_geneType.txt
    fi

    if [[ "$tempGene" =~ WZY9V-[0-9]+ ]]
    then
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCJE9V-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCJE9V=identical\t9V\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
            elif [[ "$tempGene2" =~ WCJE9V-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCJE9V=imperfect\t9A:9V_(Extract WCJE9V)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                extractFastaByID.pl $tempGenePRE2 < seroT_target_consensus.fna >> Extract_results.txt
            fi
        done < temp_geneType.txt
    fi

    if [[ "$tempGene" =~ WZY22F-[0-9]+ ]]
    then
	while read geneLine2
	do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

	    if [[ "$tempGene2" =~ WCWA22A-[0-9]+ ]]
	    then
		printf "$tempGene/$tempGene2\tWCWA22A=present\t22A\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
	    elif [[ "$tempGene2" =~ WCWA22F-[0-9]+ ]]
	    then
		printf "$tempGene/$tempGene2\tWCWA22F=present\t22F\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
	    fi
	done < temp_geneType.txt
    fi

    if [[ "$tempGene" =~ WZY18C-[0-9]+ && ! "$tempGene" =~ "*" ]]
    then
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCIX18C-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIX18C=identical\t18C\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
            elif [[ "$tempGene2" =~ WCIX18C-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCIX18C=imperfect\t18B:18C_(Extract_WCIX18C)\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
            fi
        done < temp_geneType.txt
    fi

    if [[ "$tempGene" =~ WZY33F-[0-9]+ ]]
    then
        isWCJE33A="no"
        while read geneLine2
        do
            tempGenePRE2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ -n "$geneDiff2" ]]
            then
                tempGene2="$tempGenePRE2*"
            else
                tempGene2=$tempGenePRE2
            fi
            geneDepth2=$(echo "$geneLine2" | awk -F'\t' '{ print $6 }')

            if [[ "$tempGene2" =~ WCJE33A-[0-9]+ && ! "$tempGene2" =~ "*" ]]
            then
                printf "$tempGene/$tempGene2\tWCJE33A=identical\t33A\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                isWCJE33A="yes"
            elif [[ "$tempGene2" =~ WCJE33A-[0-9]+"*" ]]
            then
                printf "$tempGene/$tempGene2\tWCJE33A=imperfect\t33F\t$geneDepth/$geneDepth2\n" >> TEMP_SeroType_Results.txt
                extractFastaByID.pl $tempGenePRE < seroT_target_consensus.fna >> Extract_results.txt
                isWCJE33A="yes"
            fi
        done < temp_geneType.txt
        if [[ "$isWCJE33A" == "no" ]]
        then
            printf "$tempGene\tWZY33F=not_present\t33F\t$geneDepth\n" >> TEMP_SeroType_Results.txt
        fi
    fi
done < temp_geneType.txt

if [[ -e Extract_results.txt ]]
then
    printf '#%.0s' {1..60} >> Check_Target_Sequence.txt
    printf "%s" "--Serotype Extracted Sequence--" >> Check_Target_Sequence.txt
    printf '#%.0s' {1..60} >> Check_Target_Sequence.txt
    printf "\n" >> Check_Target_Sequence.txt
    cat Extract_results.txt >> Check_Target_Sequence.txt
    printf '#%.0s' {1..151} >> Check_Target_Sequence.txt
    printf "\n\n" >> Check_Target_Sequence.txt
    rm Extract_results.txt
fi

###Delete temp files and unload modules###
rm temp_geneType.txt
module unload freebayes/0.9.21
module unload srst2/0.1.7

