#!/bin/bash

data_dir=/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/SPN_Typing_Output/\
18-B-429_Illumina-MiSeq-M03220/190621_M03220_0049_000000000-CHCYF_test/\
SPN_Typing_Output/20192226-S-ABC/ 

cd /scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/SPN_Scripts_Reference/bLactam_MIC_Rscripts

bash PBP_AA_sampledir_to_MIC_20180710.sh $data_dir
