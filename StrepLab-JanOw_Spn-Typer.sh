#!/bin/bash -l

PATH=/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/Spn_JanOw_Typing_Scripts/scripts:$PATH

## -- begin embedded SGE options --
read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1/job-control.txt)
## -- end embedded SGE options --

###Load Modules###
#. /usr/share/Modules/init/bash
module load samtools/0.1.18
module load bowtie2/2.1.0
module load Python/2.7
#module load tabix/0.2.6
#module load vcftools/0.1.11
module load freebayes/0.9.21

###This script is called for each job in the qsub array. The purpose of this code is to read in and parse a line of the job-control.txt file 
###created by 'StrepLab-JanOw_Spn-wrapr.sh' and pass that information, as arguments, to other programs responsible for various parts of strain 
###characterization (MLST, serotype, pili and antibiotic drug resistance prediction).

readPair_1=${PARAM[0]}
readPair_2=${PARAM[1]}
SeroT_ref=${PARAM[2]}
mlst_ref=${PARAM[3]}
mlst_def=${PARAM[4]}
bLactam_ref=${PARAM[5]}
miscDrug_ref=${PARAM[6]}
vancDrug_ref=${PARAM[7]}
allDB_dir=${PARAM[8]}
batch_out=${PARAM[9]}
sampl_out=${PARAM[10]}

###Start Doing Stuff###
cd "$sampl_out"
batch_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)}')
out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')  ###Use This For Batches off the MiSeq###
#out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-1)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')   ###Otherwise Use This###
just_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')
out_nameMLST=MLST_"$out_name"
out_nameTYPE=TYPE_"$out_name"
sampleName=$(echo "$readPair_1" | awk -F"/" '{print $NF}' | sed 's/_S1_.*_001.fastq.gz//g')
#echo "The sampleName is: $sampleName"

###Call MLST###
mod-srst2.py --mlst_delimiter '_' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMLST" --save_scores --mlst_db "$mlst_ref" --mlst_definitions "$mlst_def"
###Check and extract new MLST alleles###
MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_pneumoniae__results.txt "$out_nameMLST"__*.Streptococcus_pneumoniae.sorted.bam "$mlst_ref"

###Detect Spn serotype sequence###
mod-srst2.py --input_pe "$readPair_1" "$readPair_2" --output "$out_nameTYPE" --save_scores --min_coverage 99 --max_divergence 5 --gene_db "$SeroT_ref"

###mpileup the '.*_TYPE__.*.sorted.bam and create the called variants file with freebayes. Don't use vcf2fq b/c it won't call indels###
seroT_bam=$(ls TYPE_*sorted.bam)
seroT_pileup=$(ls TYPE_*sorted.bam | sed 's/\.bam/\.pileup/g')
seroT_vcf=$(ls TYPE_*sorted.bam | sed 's/\.bam/\.vcf/g')
seroT_bai=$(ls TYPE_*sorted.bam | sed 's/\.bam/\.bai/g')
samtools mpileup -f "$SeroT_ref" "$seroT_bam" > "$seroT_pileup"
samtools index "$seroT_bam" "$seroT_bai"
freebayes -q 20 -p 1 -f "$SeroT_ref" "$seroT_bam" -v "$seroT_vcf"
###Create the variant-called consensus fasta###
bgzip "$seroT_vcf"
tabix -p vcf "$seroT_vcf".gz
cat "$SeroT_ref" | vcf-consensus "$seroT_vcf".gz | sed 's/>[0-9]\+__.*__\(.*\)__[0-9]\+/>\1/g' > seroT_target_consensus.fna

###Call serotypes using the mod-srst2.py output###
bash seroType_pred.sh -f *__fullgenes__*__results.txt -p "$SeroT_ref"

###Locate and extract the 3 PBP protein sequences###
bash loTrac_gene.sh -1 "$readPair_1" -2 "$readPair_2" -q "$bLactam_ref" -p -n "$out_name"
###Type the PBP sequences that were extracted from loTrac###
bash bLactam-PBP_Typer.sh -a ./Final-Frag_1A-S2_prot.faa -b ./Final-Frag_2B-S2_prot.faa -x ./Final-Frag_2X-S2_prot.faa -r "$allDB_dir" -n "$just_name"
###If 'Final_newRef_updater.txt' exists print content out to the parent batch directory###
if [[ -e Final_newRef_updater.txt ]]
then
    cat Final_newRef_updater.txt >> "$batch_out"/UPDATR_Spn_"$batch_name"_bLactam_updater.txt
fi

###Predict the other drug resistances###
miscRes_Typer.sh -1 "$readPair_1" -2 "$readPair_2" -r "$miscDrug_ref" -v "$vancDrug_ref" -n "$out_name"

#: <<'END'
###Output the serotype/MLST/drug resistance data for this sample to it's results output file###
sampl_out="Spn_Typing_Results.txt"
#table_out="Spn_Typing_Results.tbl"

printf "$just_name\n" >> "$sampl_out"
###MLST OUTPUT###
printf "\tMLST:\n" >> "$sampl_out"
count=0
while read -r line
do
    count=$(( $count + 1 ))
    if [[ "$count" -eq 1 ]]
    then
	printf "\t\t$line\n" >> "$sampl_out"
	#printf "$line\t" >> TEMP_table_title.txt
    else
	printf "\t\t$line\n" >> "$sampl_out"
	MLST_tabl=$(echo "$line" | cut -f1-9)
        printf "$MLST_tabl\t" >> TEMP_table_results.txt
    fi
done < "$out_nameMLST"__mlst__Streptococcus_pneumoniae__results.txt

###SEROTYPE OUTPUT###
printf "\tSerotype:\n" >> "$sampl_out"
lineNum=$(cat Serotype_results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    firstLine=$(head -n1 Serotype_results.txt)
    printf "\t\t$firstLine\n\t\tNo_Serotype\n" >> "$sampl_out"
    #printf "$firstLine\t" >> TEMP_table_title.txt
    printf "No_Serotype\t" >> TEMP_table_results.txt
else
    count=0
    while read -r line
    do	
	count=$(( $count + 1 ))
	if [[ "$count" -eq 1 ]]
	then
	    printf "\t\t$line\n" >> "$sampl_out"
	    #printf "$line\t" >> TEMP_table_title.txt
	else
	    printf "\t\t$line\n" >> "$sampl_out"
	    justTarget=$(echo "$line" | awk -F"\t" '{print $3}')
	    printf "$justTarget;" >> TEMP_table_results.txt
	fi
    done < Serotype_results.txt
fi
printf "\t" >> TEMP_table_results.txt

###PBP_ID OUTPUT###
printf "\tPBP_ID Code:\n" >> "$sampl_out"
count=0
while read -r line
do
    count=$(( $count + 1 ))
    justPBPs=$(echo "$line" | cut -f2-4)
    if [[ "$count" -eq 1 ]]
    then
        printf "\t\t$justPBPs\n" >> "$sampl_out"
        #printf "$justPBPs\t" >> TEMP_table_title.txt
    else
        printf "\t\t$justPBPs\n" >> "$sampl_out"
        printf "$justPBPs\t" >> TEMP_table_results.txt
    fi
done < Final_pbpID_output.txt

###MISC. RESISTANCE###
printf "\tMisc. Resistance:\n" >> "$sampl_out"
lineNum=$(cat Misc_Resistance_results.txt | wc -l)
if [[ "$lineNum" -eq 1 ]]
then
    firstLine=$(head -n1 Misc_Resistance_results.txt)
    printf "\t\t$firstLine\n\t\tNo_Resistance\n" >> "$sampl_out"
    #printf "$firstLine\t" >> TEMP_table_title.txt
    printf "No_Resistance" >> TEMP_table_results.txt
else
    count=0
    while read -r line
    do
	count=$(( $count + 1 ))
	if [[ "$count" -eq 1 ]]
	then
	    printf "\t\t$line\n" >> "$sampl_out"
	    #printf "$line\t" >> TEMP_table_title.txt
	else
	    printf "\t\t$line\n" >> "$sampl_out"
            justTarget=$(echo "$line" | awk -F"\t" '{print $3}')
            printf "$justTarget;" >> TEMP_table_results.txt
	fi
    done < Misc_Resistance_results.txt
fi
printf "\n" >> "$sampl_out"
printf "\n" >> TEMP_table_results.txt
cat TEMP_table_results.txt

cat "$sampl_out" >> "$batch_out"/SAMPL_Spn_"$batch_name"_Typing_Results.txt
#cat TEMP_table_results.txt >> TEMP_table_title.txt
cat TEMP_table_results.txt >> "$batch_out"/TABLE_Spn_"$batch_name"_Typing_Results.txt
#END


###Remove files here###
#rm TEMP_table_results
#rm TEMP_table_title.txt
#rm seroT_target_consensus.fna
#rm Spn_Typing_Results.txt
#rm Serotype_results.txt
#rm seroT_target_consensus.fna
#rm Final_pbpID_output.txt
#rm Misc_Resistance_results.txt
#rm miscR_target_consensus.fna
#rm miscR_target_consensus.faa
rm Final_newRef_updater.txt
rm cutadapt_*
rm extract_nucl_blast_db*
rm nucl_*_blast.txt
rm *.pileup*
rm *.vcf*
rm srst2_SeroT_out_sorted.vcf.gz
rm *.tbi
rm *blast-out_pbp*
rm *.log
rm *_Gene-DB_Final.sam.mod
rm *_Gene-DB_Final.unsorted.bam
rm temp*
rm varCall_TEMP_Extract.faa
rm prot_best_blast.txt
rm pbpID_output.txt
rm *.fai
rm Final-Frag_FOLA_*
##rm Final-Full_FOLA_*
##rm Final-Full_RPLD1*
rm Final-Frag_RPLD1*
##rm *__results.txt
rm ./velvet_output/Graph
rm ./velvet_output/Log
rm ./velvet_output/PreGraph
rm ./velvet_output/Sequences
rm ./velvet_output/stats.txt


###Unload Modules###
module unload samtools/0.1.18
module unload bowtie2/2.1.0
module unload Python/2.7
#module unload tabix/0.2.6
#module unload vcftools/0.1.11
#module unload freebayes/9.9.2-6
module unload freebayes/0.9.21
