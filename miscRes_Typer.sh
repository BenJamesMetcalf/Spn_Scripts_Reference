#!/bin/bash -l

#. /usr/share/Modules/init/bash
#module load samtools/0.1.18
#module load bowtie2/2.1.0
#module load Python/2.7
#module load tabix/0.2.6
#module load vcftools/0.1.11
#module load freebayes/9.9.2-6
module load freebayes/0.9.21
module load srst2/0.1.7

###This script is used to predict non-bLactam antibiotic resistance using both read mapping and assembly.###

while getopts :1:2:r:o:n:v: option
do
    case $option in
        1) fastq1_path=$OPTARG;;
        2) fastq2_path=$OPTARG;;
        r) miscDrug_ref=$OPTARG;;
        o) output_dir=$OPTARG;;
        n) out_name=$OPTARG;;
	v) vanc_ref=$OPTARG;;
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

if [[ -z "$miscDrug_ref" ]]
then
    echo "The misc. drug resistance database path argument has not been given."
    exit 1
fi

if [[ -z "$vanc_ref" ]]
then
    echo "The vancomycin drug resistance database path argument has not been given."
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

if [[ -e "$miscDrug_ref"  ]]
then
    echo "The gene typing database file is in the following location: $miscDrug_ref"
else
    echo "This gene typing database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/query_file)."
    exit 1
fi

if [[ -e "$vanc_ref"  ]]
then
    echo "The vancomycin resistance typing database file is in the following location: $vanc_ref"
else
    echo "This gene typing database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/query_file)."
    exit 1
fi

eval readPair_1=$readPair_1
eval readPair_2=$readPair_2
eval miscDrug_ref=$miscDrug_ref
eval vanc_ref=$vanc_ref

if [[ ! -z "$out_name" ]]
then
    echo "The output file name prefix: $out_name"
    out_nameMiscRes="MISC_${out_name}"
    out_nameVancRes="VANC_${out_name}"
else
    out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-2)"--"$(NF-1)}')
    out_nameMiscRes=MISC_"$out_name"
    out_nameVancRes=VANC_"$out_name"
fi

cd "$out_dir"
###Detect Spn serotype sequence###
#mod-srst2.py --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMiscRes" --log --save_scores --min_coverage 99 --max_divergence 5 --gene_db "$miscDrug_ref"
#mod-srst2.py --input_pe "$readPair_1" "$readPair_2" --output "$out_nameVancRes" --log --save_scores --min_coverage 99 --max_divergence 20 --gene_db "$vanc_ref"
srst2 --samtools_args "\-A" --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMiscRes" --log --save_scores --min_coverage 99.9 --max_divergence 5 --gene_db "$miscDrug_ref"
srst2 --samtools_args "\-A" --input_pe "$readPair_1" "$readPair_2" --output "$out_nameVancRes" --log --save_scores --min_coverage 99.9 --max_divergence 20 --gene_db "$vanc_ref"

###mpileup the 'MISC__.*.sorted.bam and create the called variants file with freebayes. Don't use vcf2fq b/c it won't call indels###
miscR_bam=$(ls MISC_*sorted.bam)
miscR_vcf=$(ls MISC_*sorted.bam | sed 's/\.bam/\.vcf/g')
miscR_bai=$(ls MISC_*sorted.bam | sed 's/\.bam/\.bai/g')
#samtools mpileup -f "$miscDrug_ref" srst2_miscRes_out_sorted.bam > srst2_miscRes_out_sorted.pileup
samtools index "$miscR_bam" "$miscR_bai"
freebayes -q 20 -p 1 -f "$miscDrug_ref" "$miscR_bam" -v "$miscR_vcf"
###Create the variant-called consensus fasta and convert to amino acid sequence###
bgzip "$miscR_vcf"
tabix -p vcf "$miscR_vcf".gz
cat "$miscDrug_ref" | vcf-consensus "$miscR_vcf".gz | sed 's/>[0-9]\+__.*__\(.*\)__[0-9]\+/>\1/g' > miscR_target_consensus.fna ###JanOw_"$out_nameMiscRes"__Extract.fna###
Translate_DNA-6Frame.pl -s miscR_target_consensus.fna -f 1 -l 3 > miscR_target_consensus.faa
rm *fna.TranslatedProtein.fasta

vanc_bam=$(ls VANC_*sorted.bam)
vanc_vcf=$(ls VANC_*sorted.bam |  sed 's/\.bam/\.vcf/g')
freebayes -q 20 -p 1 -f "$vanc_ref" "$vanc_bam" -v "$vanc_vcf"

###Loop thru the SRST2 'fullgenes' gene typing output file and look for the presence of non-bLactam resistance markers###
tail -n +2 MISC_*__fullgenes__* > temp_miscDrug_type.txt
printf "Target\tMatch_Type\tResistance\tCoverage\n" >> Misc_Resistance_results.txt
sxR_tmRarray=()
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
    #echo "tempGene is $tempGene"

    ###RPLD-2 :: MQD-2 Resistance###
    if [[ "$tempGene" =~ RPLD-2"*" ]]
    then
	tempSeq=$(extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa | grep -v '>')
	if [[ "$tempSeq" =~ DAVFGIEPNKSVVFDVI ]]
        then
            printf "$tempGene\timperfect\tPossible_MQD-2_(Extract Seq)\t$geneDepth\n" >> Misc_Resistance_results.txt
            extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa >> Misc_Resistance_Extraction.txt
        fi
    fi

    ###ERMB-1 :: MLSB-1 Resistance  && ! ermBupst-1 :: telithromycinR Resistant##
    if [[ "$tempGene" =~ ERMB-1 && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\tMLSB-1\t$geneDepth\n" >> Misc_Resistance_results.txt
	isERMBupst="no"
	while read geneLine2
	do
	    tempGene2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
	    geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
	    if [[ "$tempGene2" =~ ERMBUPST-1 ]]
	    then
		isERMBupst="yes"
		break
	    fi
	done < temp_miscDrug_type.txt

	if [[ "$isERMBupst" == "no" ]]
	then
	    printf "ERMBUPST-1\tno_match\ttelithromycinR\t$geneDepth\n" >> Misc_Resistance_results.txt
	fi
    elif [[ "$tempGene" =~ ERMB-1"*" ]]
    then
        printf "$tempGene\timperfect\tMLSB-1\t$geneDepth\n" >> Misc_Resistance_results.txt
        isERMBupst="no"
        while read geneLine2
        do
            isERMBupst="no"
            tempGene2=$(echo "$geneLine2" | awk -F'\t' '{ print $4 }')
            geneDiff2=$(echo "$geneLine2" | awk -F'\t' '{ print $7 }')
            if [[ "$tempGene2" =~ ERMBUPST-1 ]]
            then
                isERMBupst="yes"
                break
            fi
        done < temp_miscDrug_type.txt

        if [[ "$isERMBupst" == "no" ]]
        then
            printf "ERMBUPST-1\tno_match\tTelithromycinR\t$geneDepth\n" >> Misc_Resistance_results.txt
        fi
    fi

    ###ERMTR-1 :: MLSB-2 Resistance###
    if [[ "$tempGene" =~ ERMTR-1 && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\tMLSB-2\t$geneDepth\n" >> Misc_Resistance_results.txt
    elif [[ "$tempGene" =~ ERMTR-1"*" ]]
    then
        printf "$tempGene\timperfect\tMLSB-2\t$geneDepth\n" >> Misc_Resistance_results.txt
    fi

    ###ERMBS-1 :: ermB+/eryS Resistance###
    if [[ "$tempGene" =~ ERMBS-1 && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\termB+/eryS\t$geneDepth\n" >> Misc_Resistance_results.txt
    fi

    ###MEF-1 :: M-1 Resistance###
    if [[ "$tempGene" =~ MEF-1 && ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\tM-1\t$geneDepth\n" >> Misc_Resistance_results.txt
    elif [[ "$tempGene" =~ MEF-1"*" ]]
    then
        printf "$tempGene\timperfect\tM-1\t$geneDepth\n" >> Misc_Resistance_results.txt
    fi

    ###RPLV-2 :: MQD-3###
    if [[ "$tempGene" =~ RPLV-2"*" ]]
    then
        tempSeq=$(extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa | grep -v '>')
        if [[ "$tempSeq" != PTMKRFRPRA ]]
        then
            printf "$tempGene\timperfect\tPossible_MQD-3_(Extract Seq)\t$geneDepth\n" >> Misc_Resistance_results.txt
            extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa >> Misc_Resistance_Extraction.txt
	fi
    fi

    ###RPLV-1 :: MQD-4###
    if [[ "$tempGene" =~ RPLV-1"*" ]]
    then
        tempSeq=$(extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa | grep -v '>')
        if [[ "$tempSeq" != KRTAHITVA ]]
        then
            printf "$tempGene\timperfect\tPossible_MQD-4_(Extract Seq)\t$geneDepth\n" >> Misc_Resistance_results.txt
            extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa >> Misc_Resistance_Extraction.txt
	fi
    fi

    ###23S-WT1 :: MLS-R1/MLS-R2###
    if [[ "$tempGene" =~ 23SWT1-1 ]]
    then
	extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.fna > temp_MiscRes_Extract.txt
	geneSeq=$(cat temp_MiscRes_Extract.txt | grep -v '>' | tr '\n' ' ' | sed 's/ //g')
	extractFastaByID.pl 1__23SWT1__23SWT1-1__1 < "$miscDrug_ref" > temp_23S-Ref_Extract.txt
	ref23S=$(cat temp_23S-Ref_Extract.txt | grep -v '>' | tr '\n' ' ' | sed 's/ //g')
	A21G=$(echo ${geneSeq:20:1})
	A20C=$(echo ${geneSeq:19:1})

	if [[ "$A21G" == "G" ]]
	then
	    printf "$tempGene\tbp21=G;\tMLS-R1\t$geneDepth\n" >> Misc_Resistance_results.txt	    
	elif [[ "$A20C" == "C" ]]
	then
	    printf "$tempGene\tbp20=C;\tMLS-R2\t$geneDepth\n" >> Misc_Resistance_results.txt
	elif [[ "$geneSeq" != "$ref23S" ]]
	then
	    printf "$tempGene\timperfect\tUnknown_(Extract Seq)\t$geneDepth\n" >> Misc_Resistance_results.txt
	    extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa >> Misc_Resistance_Extraction.txt
	fi
	rm temp_MiscRes_Extract.txt
	rm temp_23S-Ref_Extract.txt
    fi

    ###23S-WT2 :: MLS-R3###
    if [[ "$tempGene" =~ 23SWT2-1"*" ]]
    then
        extractFastaByID.pl "$tempGene" < miscR_target_consensus.fna > temp_MiscRes_Extract.txt
        geneSeq=$(cat temp_MiscRes_Extract.txt | grep -v '>' | tr '\n' ' ' | sed 's/ //g')
	#echo "C14G is ${geneSeq:13:20}"
        C14G=$(echo ${geneSeq:13:1})
	
	if [[ "$C14G" == "G" ]]
	then
            printf "$tempGene\tbp14=G\tMLS-R3\t$geneDepth\n" >> Misc_Resistance_results.txt
	else
            printf "$tempGene\timperfect\tUnknown_(Extract Seq)\t$geneDepth\n" >> Misc_Resistance_results.txt
            extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa >> Misc_Resistance_Extraction.txt
        fi
    fi

    ###GYRA-1 :: F1###
    if [[ "$tempGene" =~ GYRA-1 ]]
    then
        extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa > temp_MiscRes_Extract.txt
        geneSeq=$(cat temp_MiscRes_Extract.txt | grep -v '>' | tr '\n' ' ' | sed 's/ //g')
        #echo "S11Y-F is ${geneSeq:10:20}"
        S11YF=$(echo ${geneSeq:10:1})

        if [[ "$S11YF" == "Y" || "$S11YF" == "F" ]]
        then
	    subArray+=("AA11=$S11YF;")
	    subList=$(echo "${subArray[*]}" | sed 's/ //g')
	    printf "$tempGene\t$subList\tF1\t$geneDepth\n" >> Misc_Resistance_results.txt
        fi
    subArray=()
    rm temp_MiscRes_Extract.txt
    fi

    ###PARC-12 :: F2###
    if [[ "$tempGene" =~ PARC-12 ]]
    then
        extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa > temp_MiscRes_Extract.txt
        geneSeq=$(cat temp_MiscRes_Extract.txt | grep -v '>' | tr '\n' ' ' | sed 's/ //g')
        S2FY=$(echo ${geneSeq:1:1})
        D6NY=$(echo ${geneSeq:5:1})
        N14D=$(echo ${geneSeq:13:1})
	if [[ "$S2FY" =~ [F|Y] ]]
	then
	    subArray+=("AA2=$S2FY;")
        fi
        if [[ "$D6NY" =~ [N|Y] ]]
        then
	    subArray+=("AA6=$D6NY;")
        fi
        if [[ "$N14D" == "D" ]]
        then
            subArray+=("AA14=D;")
        fi	

	#echo "subArray is: ${subArray[*]}"
        if [[ ${#subArray[@]} -gt 0 ]]
        then
	    subList=$(echo "${subArray[*]}" | sed 's/ //g')
	    printf "$tempGene\t$subList\tF2\t$geneDepth\n" >> Misc_Resistance_results.txt
        fi
    subArray=()
    fi

    ###TETM-1/TETO-1 :: tetR###
    if [[ "$tempGene" =~ TETM-1 || "$tempGene" =~ TET0-1 ]] && [[ ! "$tempGene" =~ "*" ]]
    then
        printf "$tempGene\tidentical\ttetR\t$geneDepth\n" >> Misc_Resistance_results.txt
    elif [[ "$tempGene" =~ TETM-1"*" || "$tempGene" =~ TETO-1"*" ]]
    then
	printf "$tempGene\tidentical\ttetR\t$geneDepth\n" >> Misc_Resistance_results.txt
    fi

    ###CAT-1 :: cmR###
    if [[ "$tempGene" =~ CAT-1 && ! "$tempGene" =~ "*" ]]
    then
	printf "$tempGene\tidentical\tcmR\t$geneDepth\n" >> Misc_Resistance_results.txt
    elif [[ "$tempGene" =~ CAT-1"*" ]]
    then
	printf "$tempGene\tidentical\ttetR\t$geneDepth\n" >> Misc_Resistance_results.txt
    fi

    ###RPOB-1 :: rifR###
    if [[ "$tempGene" =~ RPOB-1 ]]
    then
        extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa > temp_MiscRes_Extract.txt
        geneSeq=$(cat temp_MiscRes_Extract.txt | grep -v '>' | tr '\n' ' ' | sed 's/ //g')
        Q8L=$(echo ${geneSeq:7:1})
	M10I=$(echo ${geneSeq:9:1})
	D11V=$(echo ${geneSeq:10:1})
	S17F=$(echo ${geneSeq:16:1})
	H21FYN=$(echo ${geneSeq:20:1})
	S26F=$(echo ${geneSeq:25:1})
	if [[ "$Q8L" == "L" ]]
	then
	    subArray+=("AA8=L;")
	fi
        if [[ "$M10I" == "I" ]]
        then
            subArray+=("AA10=I;")
        fi
        if [[ "$D11V" == "V" ]]
        then
            subArray+=("AA11=V;")
        fi
        if [[ "$S17F" == "F" ]]
        then
            subArray+=("AA17=F;")
        fi
        if [[ "$H21FYN" =~ [F|Y|N] ]]
        then
            subArray+=("AA21=$H21FYN;")
        fi
        if [[ "$S26F" == "F" ]]
        then
            subArray+=("AA26=F;")
        fi

	if [[ ! "$Q8L" == "L" && ! "$M10I" == "I" && ! "$D11V" == "V" && ! "$S17F" == "F" && ! "$H21FYN" =~ [F|Y|N] && ! "$S26F" == "F" ]]
	then
	    tempSeq=$(extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa | grep -v '>')
	    if [[ "$tempSeq" != FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL ]]
	    then
		printf "$tempGene\timperfect\tPossible_rifR_(Extract Seq)\t$geneDepth\n" >> Misc_Resistance_results.txt
		extractFastaByID.pl "$tempGenePRE" < miscR_target_consensus.faa >> Misc_Resistance_Extraction.txt
	    fi
	else
	    subList=$(echo "${subArray[*]}" | sed 's/ //g')
	    printf "$tempGene\t$subList\trifR\t$geneDepth\n" >> Misc_Resistance_results.txt
	    subArray=()
	fi
    fi
done < temp_miscDrug_type.txt


###The FOLA-1 marker for tmR resistance may contain mosiac regions and must be extracted from the genome assembly###
###FOLA-1 :: tmR###
extractFastaByID.pl 7__FOLA__FOLA-1__7 < "$miscDrug_ref" | sed 's/>.*/>FOLA/g' > temp_MiscRes_Extract.txt
bash loTrac_gene.sh -1 "$readPair_1" -2 "$readPair_2" -q temp_MiscRes_Extract.txt -p -n "$out_name"
if [[ -e Final-Frag_FOLA_prot.faa ]]
then
    geneSeq=$(cat Final-Frag_FOLA_prot.faa | grep -v '>' | tr '\n' ' ' | sed 's/ //g')
    D12A=$(echo ${geneSeq:11:1})
    I20L=$(echo ${geneSeq:19:1})

    if [[ ! "$D12A" =~ [D|A] ]]
    then
    	subArray+=("AA12=$D12A;")
    fi

    if [[ "$I20L" == "L" ]]
    then
	subArray+=("AA20=L;")
    fi
    #echo "subArray is: ${subArray[*]}"
    if [[ ${#subArray[@]} -gt 0 ]]
    then
	sxR_tmRarray+=("yes")
        subList=$(echo "${subArray[*]}" | sed 's/ //g')
	printf "FOLA-1\t$subList\ttmR\tNA\n" >> Misc_Resistance_results.txt
    fi
else
    printf "FOLA-1\tFOLA_Extraction_Error\t---\tNA\n" >> Misc_Resistance_results.txt
subArray=()
fi

###The RPLD-1 marker for MQD-1 resistance may contain mosiac regions and must be extracted from the genome assembly###
###RPLD-1 :: MQD-1 Resistance###
extractFastaByID.pl 11__RPLD__RPLD-1__11 < "$miscDrug_ref" | sed 's/>.*/>RPLD1/g' > temp_MiscRes_Extract.txt
bash loTrac_gene.sh -1 "$readPair_1" -2 "$readPair_2" -q temp_MiscRes_Extract.txt -p -n "$out_name" -L 80 -i 80
if [[ -e Final-Frag_RPLD1_prot.faa ]]
then
    tempSeqPRE=$(cat Final-Frag_RPLD1_prot.faa | grep -v '>')
    tempSeq=${tempSeqPRE:0:12}
    if [[ "$tempSeq" != KPWRQKGTGRAR ]]
    then
	printf "RPLD-1\timperfect\tPossible_MQD-1_(Extract Seq)\t$geneDepth\n" >> Misc_Resistance_results.txt
	cat Final-Frag_RPLD1_prot.faa >> Misc_Resistance_Extraction.txt
    fi
else
    printf "RPLD-1\tRPLD1_Extraction_Error\t---\tNA\n" >> Misc_Resistance_results.txt
rm temp_MiscRes_Extract.txt
fi

###vanA-G :: vanR###
if [[ -e $(ls VANC_*__fullgenes__*) ]]
then
    tail -n +2 VANC_*__fullgenes__* > temp_vanc_type.txt
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

	if [[ "$tempGene" == VANA-1 || "$tempGene" == VANB-1 || "$tempGene" == VANC-1 || "$tempGene" == VAND-1 || "$tempGene" == VANE-1 || "$tempGene" == VANG-1 ]] && [[ ! "$tempGene" =~ "*" ]]
	then
	    printf "$tempGene\timperfect\tVanR\t$geneDepth\n" >> Misc_Resistance_results.txt
	elif [[ "$tempGene" == VANA-1 || "$tempGene" == VANB-1 || "$tempGene" == VANC-1 || "$tempGene" == VAND-1 || "$tempGene" == VANE-1 || "$tempGene" == VANG-1 ]]
	then
	    printf "$tempGene\tidentical\tVanR\t$geneDepth\n" >> Misc_Resistance_results.txt
	fi
    done < temp_vanc_type.txt
fi
    
###FOLP-1 :: sxR###
if [[ -e $(ls VANC*_Final.sorted.vcf) ]]
then
    cat VANC*_Final.sorted.vcf | grep "FOLP" > TEMP_vanc-FOLP.vcf
    if [[ -s TEMP_vanc-FOLP.vcf ]]
    then
	while read snpLine
	do
	    ref_allele=$(echo "$snpLine" | awk -F"\t" '{print $4}')
	    ref_length=$(echo ${#ref_allele})
	    alt_allele=$(echo "$snpLine" | awk -F"\t" '{print $5}')
	    alt_length=$(echo ${#alt_allele})
	    lociDP=$(echo "$snpLine" | awk -F"\t" '{print $8}' | sed 's/.*;DP=\([0-9]\+\);.*/\1/g') 
	    lociPos=$(echo "$snpLine" | awk -F"\t" '{print $2}')
	
	    if [[ "$ref_length" -ne "$alt_length" ]]
	    then
		sxR_tmRarray+=("yes")
		printf "FOLP\tIndel=nt$lociPos\tsxR\t$lociDP\n" >> Misc_Resistance_results.txt
	    fi
	done < TEMP_vanc-FOLP.vcf
    fi
fi

###tmR+scR :: cotR###
if [[ ${#sxR_tmRarray[@]} -eq 2 ]]
then
    printf "FOLA/FOLP\ttmR+/scR+\tcotR\tNA\n" >> Misc_Resistance_results.txt
fi

printf '#%.0s' {1..56} >> Check_Target_Sequence.txt
printf "%s" "--Misc. Resistance Extracted Sequence--" >> Check_Target_Sequence.txt
printf '#%.0s' {1..56} >> Check_Target_Sequence.txt
printf "\n" >> Check_Target_Sequence.txt
cat Misc_Resistance_Extraction.txt >> Check_Target_Sequence.txt
printf '#%.0s' {1..151} >> Check_Target_Sequence.txt
printf "\n\n" >> Check_Target_Sequence.txt

rm temp_miscDrug_type.txt
rm varCall_TEMP__Extract.faa
rm temp_MiscRes_Extract.txt
rm temp_vanDrug-scores.txt
#rm Misc_Resistance_Extraction.txt
#rm miscR_target_consensus.faa

module unload samtools/0.1.18
module unload bowtie2/2.1.0
module unload Python/2.7
#module unload tabix/0.2.6
#module unload vcftools/0.1.11
#module unload freebayes/9.9.2-6
module unload freebayes/0.9.21
