#!/bin/bash -l

#. /usr/share/Modules/init/bash
###This wrapper script validates the input arguments and creates the job-control.txt file which is needed to submit the qsub array job to the cluster.###

while getopts :s:t:l:d:b:m:v:r:o: option
do
    case $option in
        s) batch_dir=$OPTARG;;
        t) SeroT_ref=$OPTARG;;
        l) mlst_ref=$OPTARG;;
        d) mlst_def=$OPTARG;;
        b) bLactam_ref=$OPTARG;;
        m) miscDrug_ref=$OPTARG;;
        v) vancDrug_ref=$OPTARG;;
        r) allDB_dir=$OPTARG;;
        o) output_dir=$OPTARG;;
    esac
done

###Check if batch directory and reference database directories arguments were given and if they exist###
if [[ -z "$batch_dir" ]]
then
    echo "No sequence data directory path argument given."
    exit 1
fi

if [[ -z "$allDB_dir" ]]
then
    echo "No reference database directory path argument given."
    exit 1
fi

if [[ -d "$batch_dir" ]]
then
    batch_dir=$(echo "$batch_dir" | sed 's/\/$//g')
    echo "The sequence directory is in the following location: $batch_dir"
else
    echo "This sequence directory is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/sequence_directory)."
    exit 1
fi

if [[ -d "$allDB_dir" ]]
then
    allDB_dir=$(echo "$allDB_dir" | sed 's/\/$//g')
    echo "The references directory is in the following location: $allDB_dir"
else
    echo "This reference directory is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/reference_directory)."
    exit 1
fi


###Check if individual typing reference arguements were given and if they are found within the reference database directory###
###If any of the files are missing, the script will try to find the default reference files in the given reference database directory###
if [[ -z "$SeroT_ref" ]]
then
    echo "The serotype target reference file was not given as an argument. Will look for the file 'seroT_Gene-DB_Final.fasta' in the reference database directory"
    SeroT_ref="seroT_Gene-DB_Final.fasta"
fi
if [[ -z "$mlst_ref" ]]
then
    echo "The MLST target reference file was not given as an argument. Will look for the file 'Streptococcus_pneumoniae.fasta' in the reference database directory"
    mlst_ref="Streptococcus_pneumoniae.fasta"
fi
if [[ -z "$bLactam_ref" ]]
then
    echo "The beta Lactam represenatative sequence file was not given as an argument. Will look for the file 'MOD_bLactam_resistance.fasta' in the reference database directory"
    bLactam_ref="MOD_bLactam_resistance.fasta"
fi
if [[ -z "$miscDrug_ref" ]]
then
    echo "The misc. drug target reference file was not given as an argument. Will look for the file 'miscDrug_Gene-DB_Final.fasta' in the reference database directory"
    miscDrug_ref="miscDrug_Gene-DB_Final.fasta"
fi
if [[ -z "$vancDrug_ref" ]]
then
    echo "The Vancomycin drug target reference file was not given as an argument. Will look for the file 'vanDrug_Gene-DB_Final.fasta' in the reference database directory"
    vancDrug_ref="vanDrug_Gene-DB_Final.fasta"
fi
if [[ -z "$mlst_def" ]]
then
    echo "The MLST definitions file was not given as an argument. Will look for the file 'spneumoniae.txt' in the reference database directory"
    mlst_def="spneumoniae.txt"
fi

ST_db_path="$allDB_dir/$SeroT_ref"
mlst_db_path="$allDB_dir/$mlst_ref"
mlst_def_path="$allDB_dir/$mlst_def"
bLactam_db_path="$allDB_dir/$bLactam_ref"
misc_db_path="$allDB_dir/$miscDrug_ref"
vanc_db_path="$allDB_dir/$vancDrug_ref"

if [[ -e "$ST_db_path" ]]
then
    echo "The gene typing database file is in the following location: $ST_db_path"
else
    echo "This gene typing database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/file)."
    exit 1
fi

if [[ -e "$mlst_db_path" ]]
then
    echo "The MLST database file is in the following location: $mlst_db_path"
else
    echo "This MLST database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/file)."
    exit 1
fi

if [[ -e "$mlst_def_path" ]]
then
    echo "The MLST definition file is in the following location: $mlst_def_path"
elif [[ -n "$mlst_def_path" && ! -e "$mlst_def_path" ]]
then
    echo "This MLST definitions argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/file)."
    exit 1
fi

if [[ -e "$bLactam_db_path" ]]
then
    echo "The bLactam database file is in the following location: $bLactam_db_path"
else
    echo "This bLactam database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/file)."
    exit 1
fi

if [[ -e "$misc_db_path" ]]
then
    echo "The misc. drug resistances database file is in the following location: $misc_db_path"
else
    echo "This misc. drug resistances database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/file)."
    exit 1
fi

if [[ -e "$vanc_db_path" ]]
then
    echo "The vancomycin resistance database file is in the following location: $vanc_db_path"
else
    echo "This vancomycin resistance database argument is not in the correct format or doesn't exist."
    echo "Make sure you provide the full directory path (/root/path/file)."
    exit 1
fi

###Delete the 'job-control.txt' file if it already exists.###
if [[ -e $HOME/job-control.txt ]]
then
    echo "job-control.txt already exists. Will create a new job-control.txt file."
    rm $HOME/job-control.txt
fi

###Check if the output directory argument has been given. If yes, create the 'Analysis_Folder' and 'qsub_files' folders within the output dir###
###If no, output the results into a subdirectory of '~/Out_Spn_Sur'. The subdirectory name is extracted from the batch sequence full path###
if [[ -z "$output_dir" ]]
then
    echo "The files will be output into the default directory 'Spn_Typing_Analysis'."
    if [[ ! -d ~/Spn_Typing_Analysis ]]
    then
        mkdir ~/Spn_Typing_Analysis
        out_dir="~/Spn_Typing_Analysis"
        eval out_dir=$out_dir
        echo "The output directory has been created: $out_dir"
    else
        out_dir="~/Spn_Typing_Analysis"
        eval out_dir=$out_dir
    fi
    batch_name=$(echo "$batch_dir" | awk -F"/" '{print $(NF-3)}')
    out_analysis="${out_dir}"/"${batch_name}"/Spn_Typing_Output
    out_qsub="${out_dir}"/"${batch_name}"/qsub_files/
    out_jobCntrl="${out_dir}/${batch_name}/"
    eval out_analysis=$out_analysis
    eval out_qsub=$out_qsub
    eval out_jobCntrl=$out_jobCntrl
    mkdir -p "$out_analysis"
    mkdir -p "$out_qsub"
elif [[ ! -d "$output_dir" ]]
then
    output_dir=$(echo "$output_dir" | sed 's/\/$//g')
    mkdir "$output_dir"
    out_dir="$output_dir"
    eval out_dir=$out_dir
    echo "The output directory has been created: $out_dir"
    out_analysis="${out_dir}"/Spn_Typing_Output
    out_qsub="${out_dir}"/qsub_files/
    out_jobCntrl="${out_dir}/"
    eval out_analysis=$out_analysis
    eval out_qsub=$out_qsub
    eval out_jobCntrl=$out_jobCntrl
    mkdir -p "$out_analysis"
    mkdir -p "$out_qsub"
else
    output_dir=$(echo "$output_dir" | sed 's/\/$//g')
    out_dir="$output_dir"
    eval out_dir=$out_dir
    out_analysis="${out_dir}"/Spn_Typing_Output
    out_qsub="${out_dir}"/qsub_files/
    out_jobCntrl="${out_dir}/"
    eval out_analysis=$out_analysis
    eval out_qsub=$out_qsub
    eval out_jobCntrl=$out_jobCntrl
    mkdir -p "$out_analysis"
    mkdir -p "$out_qsub"
fi

eval ST_db_path=$ST_db_path
eval mlst_db_path=$mlst_db_path
eval mlst_def_path=$mlst_def_path
eval batch_dir=$batch_dir
eval bLactam_db_path=$bLactam_db_path
eval misc_db_path=$misc_db_path
eval vanc_db_path=$vanc_db_path

###Create the batch output files###
batch_name=$(echo "$batch_dir" | awk -F"/" '{print $(NF-3)}')
printf "Sample\tST\taroe\tgdh_\tgki_\trecP\tspi_\txpt_\tddl_\tSerotype\tPBP_1A\tPBP_2B\tPBP_2X\tMisc_Resistance\n" >> "$out_analysis"/TABLE_Spn_"$batch_name"_Typing_Results.txt

###Will search thru every file in the batch directory and check if it matches the following regexs: _L.*_R1_001.fastq and _L.*_R2_001.fastq###
###If both paired end fastq files are found then the full paths of each file will be written to the 'job-control.txt' file###
batch_dir_star="${batch_dir}/*"
for sample in $batch_dir_star
do
    sampl_name=$(echo "$sample" | sed 's/^.*\///g' | sed 's/_S[0-9]\+_.*_001.fastq.gz//g')
    #echo "Sample_name is: $sampl_name"
    #sampl_out="${batch_out="${out_dir}/Analysis_Folder"}/${sampl_name}"
    sampl_out="${out_analysis}"/"${sampl_name}"
    eval sampl_out=$sampl_out
    echo The sample file is: $sample
    if [[ $sampl_name =~ ^Undetermined ]]
    then
        echo "Skipping the 'Undetermined' fastq files"
        continue
    fi

    if [[ $sample =~ _L.*_R1_001.fastq && ! $sample =~ S[0-9]+ ]]
    then
        readPair_1=$(echo "$sample" | sed 's/_L\([0-9]\+\)_R1/_S1_L\1_R1/g')
        mv $sample $readPair_1
    elif [[ $sample =~ _L.*_R1_001.fastq && $sample =~ S[0-9]+ ]]
    then
        readPair_1=$sample
    fi

    if [[ $sample =~ _L.*_R2_001.fastq && ! $sample =~ S[0-9]+ ]]
    then
        readPair_2=$(echo "$sample" | sed 's/_L\([0-9]\+\)_R2/_S1_L\1_R2/g')
        mv $sample $readPair_2
    elif [[ $sample =~ _L.*_R2_001.fastq && $sample =~ S[0-9]+ ]]
    then
        readPair_2=$sample
    fi

    if [ -n "$readPair_1" -a -n "$readPair_2" ]
    then
	name1=$(echo $readPair_1 | sed 's/^.*\///g' | sed 's/_S[0-9]\+_.*_001.fastq.gz//g')
	name2=$(echo $readPair_2 | sed 's/^.*\///g' | sed 's/_S[0-9]\+_.*_001.fastq.gz//g')
	if [[ "$name1" == "$name2" ]]
	then
            if [[ ! -d "$sampl_out" ]]
            then
		mkdir "$sampl_out"
            fi
            echo "Both Forward and Reverse Read files exist."
            echo "Paired-end Read-1 is: $readPair_1"
            echo "Paired-end Read-2 is: $readPair_2"
            printf "\n"
            #qsubName="$sampl_name"_clust
            echo "$readPair_1 $readPair_2 $ST_db_path $mlst_db_path $mlst_def_path $bLactam_db_path $misc_db_path $vanc_db_path $allDB_dir $out_analysis $sampl_out" >> $out_jobCntrl/job-control.txt
            ###Prepare script for next sample###
            readPair_1=""
            readPair_2=""
	else 
	    echo "$readPair_1 and $readPair_2 did not match"
	    readPair_1=""
            readPair_2=""
	fi
    fi
done
qsub -t 1-$(cat $out_jobCntrl/job-control.txt | wc -l) -cwd -o "$out_qsub" -e "$out_qsub" ./StrepLab-JanOw_Spn-Typer.sh $out_jobCntrl
