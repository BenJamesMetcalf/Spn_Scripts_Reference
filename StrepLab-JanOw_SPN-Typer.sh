#!/bin/bash -l

temp_path=$(pwd)
export PATH=$PATH:$temp_path

## -- begin embedded SGE options --
read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1/job-control.txt)
## -- end embedded SGE options --

###Load Modules###
#. /usr/share/Modules/init/bash
module load perl/5.22.1
module load ncbi-blast+/2.2.29
module load BEDTools/2.17.0
module load freebayes/0.9.21
module load prodigal/2.60
module load cutadapt/1.8.3
module load srst2/0.1.7

###This script is called for each job in the qsub array. The purpose of this code is to read in and parse a line of the job-control.txt file
###created by 'StrepLab-JanOw_GAS-wrapr.sh' and pass that information, as arguments, to other programs responsible for various parts of strain
###characterization (MLST, emm type and antibiotic drug resistance prediction).

readPair_1=${PARAM[0]}
readPair_2=${PARAM[1]}
allDB_dir=${PARAM[2]}
batch_out=${PARAM[3]}
sampl_out=${PARAM[4]}


###Start Doing Stuff###
cd "$sampl_out"
batch_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)}')
out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-4)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')  ###Use This For Batches off the MiSeq###
#out_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF-1)"--"$(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')   ###Otherwise Use This###
just_name=$(echo "$readPair_1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\+_L[0-9]\+_R[0-9]\+.*//g')
out_nameMLST=MLST_"$just_name"

###Pre-Process Paired-end Reads###
fastq1_trimd=cutadapt_"$just_name"_S1_L001_R1_001.fastq
fastq2_trimd=cutadapt_"$just_name"_S1_L001_R2_001.fastq
cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 --paired-output temp2.fastq -o temp1.fastq $readPair_1 $readPair_2
cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 --paired-output $fastq1_trimd -o $fastq2_trimd temp2.fastq temp1.fastq
rm temp1.fastq
rm temp2.fastq
module load fastqc/0.11.5
mkdir "$just_name"_R1_cut
mkdir "$just_name"_R2_cut
fastqc "$fastq1_trimd" --outdir=./"$just_name"_R1_cut
fastqc "$fastq2_trimd" --outdir=./"$just_name"_R2_cut
module unload fastqc/0.11.5

###Call MLST###
srst2 --samtools_args '\\-A' --mlst_delimiter '_' --input_pe "$readPair_1" "$readPair_2" --output "$out_nameMLST" --save_scores --mlst_db "$allDB_dir/Streptococcus_pneumoniae.fasta" --mlst_definitions "$allDB_dir/spneumoniae.txt" --min_coverage 99.999
###Check and extract new MLST alleles###
MLST_allele_checkr.pl "$out_nameMLST"__mlst__Streptococcus_pneumoniae__results.txt "$out_nameMLST"__*.Streptococcus_pneumoniae.sorted.bam "$allDB_dir/Streptococcus_pneumoniae.fasta"

###Call GBS Serotype###
SPN_Serotyper.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/SPN_Sero_Gene-DB_Final.fasta" -n "$just_name"

###Call GBS bLactam Resistances###
module unload perl/5.22.1
module load perl/5.16.1-MT
PBP-Gene_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -r "$allDB_dir/MOD_bLactam_resistance.fasta" -n "$just_name" -s SPN -p 1A,2B,2X
module unload perl/5.16.1-MT
module load perl/5.22.1

###Predict bLactam MIC###
scr1="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/share/PBP_AA_to_MIC/scripts//PBP_AA_sampledir_to_MIC_20180710.sh"
bash "$scr1" "$sampl_out"

###Call GBS Misc. Resistances###
SPN_Res_Typer.pl -1 "$readPair_1" -2 "$readPair_2" -d "$allDB_dir" -r SPN_Res_Gene-DB_Final.fasta -n "$just_name"
SPN_Target2MIC.pl OUT_Res_Results.txt "$just_name"

###Output the emm type/MLST/drug resistance data for this sample to it's results output file###
tabl_out="TABLE_Isolate_Typing_results.txt"
bin_out="BIN_Isolate_Typing_results.txt"
printf "$just_name\t" >> "$tabl_out"
printf "$just_name," >> "$bin_out"

###Serotype Output###
sero_out="NF"
pili_out="neg"
while read -r line
do
    if [[ -n "$line" ]]
    then
        justTarget=$(echo "$line" | awk -F"\t" '{print $4}')
	if [[ "$justTarget" == "PI-1" ]]
	then
            if [[ "$pili_out" == "neg" ]]
            then
		pili_out="1"
            elif [[ "$pili_out" == "2" ]]
	    then
		pili_out="1:2"
            fi
	elif [[ "$justTarget" == "PI-2" ]]
	then
            if [[ "$pili_out" == "neg" ]]
            then
                pili_out="2"
            elif [[ "$pili_out" == "1" ]]
	    then
                pili_out="1:2"
            fi
        else
            if [[ "$sero_out" == "NF" ]]
            then
		sero_out="$justTarget"
            else
		sero_out="$sero_out;$justTarget"
            fi
	fi
    fi
done <<< "$(sed 1d OUT_SeroType_Results.txt)"
printf "$sero_out\t$pili_out\t" >> "$tabl_out"
printf "$sero_out,$pili_out\t" >> "$bin_out"

###MLST OUTPUT###
sed 1d "$out_nameMLST"__mlst__Streptococcus_pneumoniae__results.txt | while read -r line
do
    MLST_tabl=$(echo "$line" | cut -f2-9)
    echo "MLST line: $MLST_tabl\n";
    printf "$MLST_tabl\t" >> "$tabl_out"
    MLST_val=$(echo "$line" | awk -F" " '{print $2}')
    printf "$MLST_val," >> "$bin_out"
done

###PBP_ID Output###
justPBPs="NF"
sed 1d TEMP_pbpID_Results.txt | while read -r line
do
    if [[ -n "$line" ]]
    then
        justPBPs=$(echo "$line" | awk -F"\t" '{print $2}' | tr ':' '\t')
        justPBP_BIN=$(echo "$line" | awk -F"\t" '{print $2}' | tr ':' ',')
    fi
    printf "$justPBPs\t" >> "$tabl_out"
    printf "$justPBP_BIN," >> "$bin_out"
done

###bLactam Predictions###
#sed 1d "BLACTAM_MIC_RF.txt" | while read -r line
#sed 1d "BLACTAM_MIC_RF_with_SIR.txt" | while read -r line
#do
#    pbpID=$(tail -n1 "TEMP_pbpID_Results.txt" | awk -F"\t" '{print $2}')
#    if [[ ! "$pbpID" =~ .*NF.* ]] && [[ ! "$pbpID" =~ .*NEW.* ]]
#    then
#	echo "No NF or NEW outputs for PBP Type"
#	bLacTab=$(echo "$line" | tr ' ' '\t')
#	printf "$bLacTab\t" >> "$tabl_out"
#	bLacCom=$(echo "$line" | tr ' ' ',')
#	printf "$bLacCom," >> "$bin_out"
#    else
#	echo "One of the PBP types has an NF or NEW"
#	printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\t" >> "$tabl_out"
#	printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF," >> "$bin_out"
#    fi
#done

pbpID=$(tail -n1 "TEMP_pbpID_Results.txt" | awk -F"\t" '{print $2}')
if [[ ! "$pbpID" =~ .*NF.* ]] && [[ ! "$pbpID" =~ .*NEW.* ]]
then
    echo "No NF or NEW outputs for PBP Type"
    bLacTab=$(tail -n1 "BLACTAM_MIC_RF_with_SIR.txt" | tr ' ' '\t')
    printf "$bLacTab\t" >> "$tabl_out"
    #bLacCom=$(echo "$line" | tr ' ' ',')
    #printf "$bLacCom," >> "$bin_out"
else
    echo "One of the PBP types has an NF or NEW"
    printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\t" >> "$tabl_out"
    #printf "NF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF," >> "$bin_out"
fi


###Resistance Targets###
while read -r line
do
    #RES_targ=$(echo "$line" | cut -f2)
    #printf "$RES_targ\t" >> "$tabl_out"
    printf "$line\t" | tr ',' '\t' >> "$tabl_out"
done < RES-MIC_"$just_name"

if [[ -e $(echo ./velvet_output/*_Logfile.txt) ]]
then
    vel_metrics=$(echo ./velvet_output/*_Logfile.txt)
    print "velvet metrics file: $vel_metrics\n";
    velvetMetrics.pl -i "$vel_metrics";
    line=$(cat velvet_qual_metrics.txt | tr ',' '\t')
    printf "$line\t" >> "$tabl_out"

    printf "$readPair_1\t" >> "$tabl_out";
    pwd | xargs -I{} echo {}"/velvet_output/contigs.fa" >> "$tabl_out"
else
    printf "NA\tNA\tNA\tNA\t$readPair_1\tNA\n" >> "$tabl_out"
fi
#printf "\n" >> "$tabl_out"

#printf "\n" >> "$bin_out"
#cat BIN_Res_Results.txt | sed 's/$/,/g' >> "$bin_out"


###Remove Temporary Files###
#rm cutadapt*.fastq
#rm *.pileup
#rm *.bam
#rm *.sam
#rm TEMP*

###Unload Modules###
module unload perl/5.22.1
module unload ncbi-blast+/2.2.29
module unload BEDTools/2.17.0
module unload freebayes/0.9.21
module unload prodigal/2.60
module unload cutadapt/1.8.3
module unload srst2/0.1.7
