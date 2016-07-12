#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;

###MODULE LOAD###
#module load samtools/0.1.18
#module load bowtie2/2.1.0
#module load Python/2.7
#module load freebayes/0.9.21

sub checkOptions {
    my %opts;
    getopts('h1:2:r:m:v:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $ref_dir, $misc_DB, $vanc_DB, $outDir, $outName);

    if($opts{h}) {
        $help = $opts{h};
        help();
    }

    if($opts{1}) {
        $fastq1 = $opts{1};
        if( -e $fastq1) {
            print "Paired-end Read 1 is: $fastq1\n";
        } else {
            print "The forward paired-end file name is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "No paired end 1 fastq file path argument given.\n";
        help();
    }

    if($opts{2}) {
        $fastq2 = $opts{2};
        if( -e $fastq2) {
            print "Paired-end Read 2 is: $fastq2\n";
        } else {
            print "The reverse paired-end file name is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "No paired end 2 fastq file path argument given.\n";
        help();
    }

    if($opts{r}) {
        $ref_dir = $opts{r};
        $ref_dir =~ s/\/$//g;
        if (-d $ref_dir) {
            print "Directory containing the misc. and vancomycin resistance reference sequences: $ref_dir\n";
        } else {
            print "The directory containing the misc. and vancomycin resistance references is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The misc. and vancomycin resistance reference sequence directory (including full path) has not been given.\n";
        help();
    }

    if($opts{m}) {
        $misc_DB = "$ref_dir/$opts{m}";
        if ($misc_DB) {
            print "The misc. resistance reference sequence file: $opts{m}\n";
        } else {
            print "The misc resistance reference sequence file is not in the correct format or doesn't exist.\n";
            #print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The misc. resistance reference sequence file has not been given.\n";
        help();
    }

    if($opts{v}) {
        $vanc_DB = "$ref_dir/$opts{v}";
        if ($vanc_DB) {
            print "The vancomycin resistance reference sequence file: $opts{v}\n";
        } else {
            print "The vancomycin resistance reference sequence file is not in the correct format or doesn't exist.\n";
            #print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The vancomycin resistance reference sequence file has not been given.\n";
        help();
    }

    $outDir = "./";
    if($opts{o}) {
        if (-d $opts{o}) {
            $outDir = $opts{o};
            print "The output directory is: $outDir\n";
        } else {
            $outDir = $opts{o};
            mkdir $outDir;
            print "The output directory has been created: $outDir\n";
        }
    } else {
        print "The files will be output into the current directory.\n";
    }

    if($opts{n}) {
        $outName = $opts{n};
        print "The output file name prefix: $outName\n";
    } else {
        $outName = `echo "$fastq1" | awk -F"/" '{print \$(NF)}' | sed 's/_S[0-9]\\+_L[0-9]\\+_R[0-9]\\+.*//g'`;
        print "The default output file name prefix is: $outName";
    }

    return ($help, $fastq1, $fastq2, $ref_dir, $misc_DB, $vanc_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
SPN_miscRes_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <ref directory: dir> -m <misc. resistance seq: file> -v <vanc. resistance seq: file> -o <output directory name: string> -n <output name prefix: string> [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   reference sequence directory (including full path)
    -m   misc. resistance reference sequence file
    -v   vancomycin resistance reference sequence file
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $ref_dir, $misc_DB, $vanc_DB, $outDir, $outName) = checkOptions( @ARGV );




###Subroutines###
sub sixFrame_Translate {
    my ($seq_input,$opt_f) = @_;

    sub codon2aa{
        my($codon)=@_;
        $codon=uc $codon;
        my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

        if(exists $g{$codon}){return $g{$codon};}
        elsif($codon=~/GC./i){return 'A';}
        elsif($codon=~/GG./i){return 'G';}
        elsif($codon=~/CC./i){return 'P';}
        elsif($codon=~/AC./i){return 'T';}
        elsif($codon=~/GT./i){return 'V';}
        elsif($codon=~/CG./i){return 'R';}
        elsif($codon=~/TC./i){return 'S';}
        else {
            return('x');
            print "Bad codon \"$codon\"!!\n";
        }
    }

    (my $DNAheader,my @DANseq)=split(/\n/,$seq_input);
    chomp $DNAheader;
    $DNAheader=~s/\s+$//g;
    my $DNAseq = join( '',@DANseq);
    $DNAseq =~ s/\s//g;
    $DNAheader=~s/>//g;
    $DNAseq=~s/>//g;
    my$DNA_length=length$DNAseq;
    #print "\nSeq:$DNAheader\t:$DNA_length nt\n\n";
    my $DNArevSeq = reverse($DNAseq);
    $DNArevSeq=~tr/ATGCatgc/TACGtacg/;
    #print "\nThe original DNA sequence is:\n$DNAseq \nThe reverse of DNA sequence is:\n$DNArevSeq\n";
    my @protein='';
    my @dna='';
    my $codon1;

    if ($opt_f == 1) {
        for(my $i=0;$i<(length($DNAseq)-2);$i+=3){
            $codon1=substr($DNAseq,$i,3);
            $protein[1].= codon2aa($codon1);
            #$dna[1].=codon2nt($codon1);
            }
        }
    if ($opt_f == 2) {
            my $codon2;
            for(my $i=1;$i<(length($DNAseq)-2);$i+=3){
                $codon2=substr($DNAseq,$i,3);
                $protein[2].= codon2aa($codon2);
                #$dna[2].=codon2nt($codon2);
                }
    }
    if ($opt_f == 3) {
            my $codon3;
            for(my $i=2;$i<(length($DNAseq)-2);$i+=3){
                $codon3=substr($DNAseq,$i,3);
                $protein[3].= codon2aa($codon3);
                #$dna[3].=codon2nt($codon3);
                }
    }
    if ($opt_f == 4) {
            my $codon4;
            for(my $i=0;$i<(length($DNArevSeq)-2);$i+=3){
                $codon4=substr($DNArevSeq,$i,3);
                $protein[4].= codon2aa($codon4);
                #$dna[4].=codon2nt($codon4);
                }
    }
    if ($opt_f == 5) {
            my $codon5;
            for(my $i=1;$i<(length($DNArevSeq)-2);$i+=3){
                $codon5=substr($DNArevSeq,$i,3);
                $protein[5].= codon2aa($codon5);
                #$dna[5].=codon2nt($codon5);
                }
    }
    if ($opt_f == 6) {
            my $codon6;
                for(my $i=2;$i<(length($DNArevSeq)-2);$i+=3){
                    $codon6=substr($DNArevSeq,$i,3);
                    $protein[6].= codon2aa($codon6);
                    #$dna[6].=codon2nt($codon6);
                    }
    }
    #print "translate result\n$protein[$opt_f]\n";
return $protein[$opt_f];
}

sub extractSeqByID {
    my ($lookup, $reference) = @_;
    open my $fh, "<", $reference or die $!;
    local $/ = "\n>";  # read by FASTA record

    my $output;
    while (my $seq = <$fh>) {
        chomp $seq;
        #print "while seq:\n$seq\n";
        my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        if ($id eq $lookup) {
            $seq =~ s/^>*.+\n//;  # remove FASTA header
            #$seq =~ s/\n//g;  # remove endlines
            #print ">$id\n";
            #print "$seq\n";
            #$output = ">$id\n$seq\n";
            $output = $seq;
            last;
        }
    }
return $output;
}

sub extractFastaByID {
    my ($lookup, $reference) = @_;
    open my $fh, "<", $reference or die $!;
    local $/ = "\n>";  # read by FASTA record

    my $output;
    while (my $seq = <$fh>) {
        chomp $seq;
        #print "while seq:\n$seq\n";
        my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        if ($id eq $lookup) {
            $seq =~ s/^>*.+\n//;  # remove FASTA header
            $seq =~ s/\n//g;  # remove endlines
            #print ">$id\n";
            #print "$seq\n";
            $output = ">$id\n$seq\n";
            #$output = $seq;
            last;
        }
    }
return $output;
}

sub freebayes_prior_fix {
    my ($bamFile, $refFile, $target) = @_;
    (my $samFile = $bamFile) =~ s/\.bam/\.sam/g;
    system("samtools view -h $bamFile > $samFile");
    system("cat $samFile | grep -E \"^\@HD|^\@SQ.*$target|^\@PG\" > CHECK_target_seq.sam");
    system("awk -F'\t' '\$3 == \"$target\" {print \$0}' $samFile >> CHECK_target_seq.sam");
    system("samtools view -bS CHECK_target_seq.sam > CHECK_target_seq.bam");
    system("samtools index CHECK_target_seq.bam CHECK_target_seq.bai");
    my $REF_seq = extractFastaByID("$target","$refFile");
    open(my $rf,'>',"CHECK_target_ref.fna");
    print $rf "$REF_seq\n";
    close $rf;
    system("freebayes -q 20 -p 1 -f CHECK_target_ref.fna CHECK_target_seq.bam -v CHECK_target_seq.vcf");
    system("bgzip CHECK_target_seq.vcf");
    system("tabix -p vcf CHECK_target_seq.vcf.gz");
    my $extractSeq = `echo "$REF_seq" | vcf-consensus CHECK_target_seq.vcf.gz`;
    chomp($extractSeq);
    #print "$target-----------------------------------\n";
    #system("cat CHECK_target_seq.sam");
    #system("zcat CHECK_target_seq.vcf.gz");
    #print "reference seq:\n$REF_seq\n";
    #print "extracted Seq:\n$extractSeq\n";
    #print "$target-----------------------------------\n";
    system("rm CHECK_target*");
return $extractSeq;
}





###Start Doing Stuff###
chdir "$outDir";
my $miscR_output = "TEMP_miscR_Results.txt";
open(my $fh,'>',$miscR_output) or die "Could not open file '$miscR_output' $!";
print $fh "Target\tMatch_Type\tResistance\tCoverage\n";

my $justName = `echo $outName | sed 's/MISC_//g'`;
my $outName_VANC = "VANC_".$justName;
chomp($outName_VANC);
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $outName --log --save_scores --min_coverage 99.9 --max_divergence 5 --gene_db $misc_DB");
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $outName_VANC --log --save_scores --min_coverage 95.0 --max_divergence 15 --gene_db $vanc_DB");

my @TEMP_MISC_bam = glob("MISC_*\.sorted\.bam");
my @TEMP_MISC_fullgene = glob("MISC_*__fullgenes__*__results\.txt");
my $MISC_bam = $TEMP_MISC_bam[0];
my $MISC_full_name = $TEMP_MISC_fullgene[0];
print "misc bam is: $MISC_bam || misc full gene $MISC_full_name\n";
(my $MISC_vcf = $MISC_bam) =~ s/\.bam/\.vcf/g;
(my $MISC_bai = $MISC_bam) =~ s/\.bam/\.bai/g;

my %drugRes_Targets = (
    "ERMBS-1" => "eryS",
    "ERMB-1" => "MLSB-1",
    "ERMTR-1" => "MLSB-2",
    "MEF-1" => "M-1",
    "CAT-1" => "cmR",
    "TETM-1" => "tetR",
    "TETO-1" => "tetR",
    "VANA-1" => "vanR",
    "VANB-1" => "vanR",
    "VANC-1" => "vanR",
    "VAND-1" => "vanR",
    "VANE-1" => "vanR",
    "VANG-1" => "vanR",
    );

my $isNotResistant = "yes";
###Type the non-vancomycin resistance targets###
my %miscR_extract;
open(MYINPUTFILE, "$MISC_full_name");
my %miscR_Type;
my @miscR_fullgene;
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    @miscR_fullgene = split('\t',$line);
    my $isMismatch = "no";
    if ($miscR_fullgene[6]) {
        $isMismatch = "yes";
    }
    my $value_out = "$miscR_fullgene[5]:$isMismatch";
    $miscR_Type{$miscR_fullgene[3]} = $value_out;
}
print Dumper(\%miscR_Type);

foreach my $key (keys(%miscR_Type)) {
    #print "$miscR_Type{$key}\n";
    if ($drugRes_Targets{$key}) {
        my @miscR_value = split(':',$miscR_Type{$key});
        my $status = "identical";
        if ($miscR_value[1] eq "yes") {
            $status = "imperfect";
        }
        print $fh "$key\t$status\t$drugRes_Targets{$key}\t$miscR_value[0]\n";
        $isNotResistant = "no";
    }
}


###Now Type the Non-Presence/Absence Targets###
my @sxR_tmR_array;

################################################################
###Type the ERMB-1 Telithromycin targets###
if ($miscR_Type{"ERMB-1"} && (! $miscR_Type{"ERMBUPST-1"})) {
    my @miscR_value = split(':',$miscR_Type{"ERMB-1"});
    print $fh "ERMB-1/ERMBUPST-1\tmatch/no_match\ttelR\t$miscR_value[0]\n";
}
################################################################


################################################################
###Type the RPOB-1 Rifampicin targets###
if ($miscR_Type{"RPOB-1"}) {
    my $RPOB_seq = freebayes_prior_fix($MISC_bam, $misc_DB, "15__RPOB__RPOB-1__15");
    my $RPOB_aaSeq = sixFrame_Translate($RPOB_seq,1);
    my $RPOB_aaRef = "FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL";
    print "RPOB ref: $RPOB_aaRef || RPOB seq: $RPOB_aaSeq\n";
    if ($RPOB_aaSeq ne "FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL") {
	my $mask = $RPOB_aaSeq ^ $RPOB_aaRef;
	my @seq_diffs;
	while ($mask =~ /[^\0]/g) {
	    print substr($RPOB_aaRef,$-[0],1), ' ', substr($RPOB_aaSeq,$-[0],1), ' ', $-[0], "\n";
	    my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
	    push(@seq_diffs,$diff_element);
	}
	print "RPOB seq: $RPOB_seq\n";
	my $diff_output = join(';',@seq_diffs);
	my @miscR_value = split(':',$miscR_Type{"RPOB-1"});
	print $fh "RPOB-1\t$diff_output\tPossible_Rif_(Extract Seq)\t$miscR_value[0]\n";
	$miscR_extract{"RPOB-1"} = $RPOB_seq;
    }
    ###Just in case I need these for the future###
    #my $RPOB_pos8 = substr($RPOB_aaSeq,7,1);
    #my $RPOB_pos10 = substr($RPOB_aaSeq,9,1);
    #my $RPOB_pos11 = substr($RPOB_aaSeq,10,1);
    #my $RPOB_pos17 = substr($RPOB_aaSeq,16,1);
    #my $RPOB_pos21 = substr($RPOB_aaSeq,20,1);
    #my $RPOB_pos26 = substr($RPOB_aaSeq,25,1);
}
#################################################################


#################################################################
###Type the FOLA-1 Trimethoprim targets###
###The FOLA-1 marker for tmR resistance may contain mosiac regions and must be extracted from the genome assembly###
my $REF_seq = extractFastaByID("7__FOLA__FOLA-1__7","$misc_DB");
`echo "$REF_seq" > TEMP_FOLA_Ref.fna`;
system("LoTrac_target.pl -1 $fastq1 -2 $fastq2 -q TEMP_FOLA_Ref.fna -S 2.2M -f -n $justName -o $outDir");
my $FOLA_file = glob("EXTRACT_*FOLA*.fasta");
my $FOLA_error = glob("ERROR_*FOLA*.fasta");
my @FOLA_output;
if ($FOLA_file && ! $FOLA_error) {
    my $FOLA_seq_hedr = `cat $FOLA_file | grep ">" | tail -n1 | sed 's/>//g'`;
    chomp($FOLA_seq_hedr);
    my $FOLA_seq = extractFastaByID("$FOLA_seq_hedr","$FOLA_file");
    chomp($FOLA_seq);
    #print "FOLA extract header: $FOLA_seq_hedr | FOLA extract file: $FOLA_file\n";
    my $FOLA_aaSeq = sixFrame_Translate($FOLA_seq,1);
    my $FOLA_pos12 = substr($FOLA_aaSeq,11,1);
    my $FOLA_pos20 = substr($FOLA_aaSeq,19,1);
    print "FOLA aaSeq: $FOLA_seq || pos12: $FOLA_pos12 || pos20: $FOLA_pos20\n";
    if ($FOLA_pos12 ne "D" && $FOLA_pos12 ne "A" && $FOLA_pos12 ne "V") {
        my $pos_output = "pos12=$FOLA_pos12";
        push(@FOLA_output,$pos_output);
    }
    if ($FOLA_pos20 eq "L") {
        my $pos_output = "pos20=$FOLA_pos20";
        push(@FOLA_output,$pos_output);	
    }

    if (@FOLA_output) {
	push(@sxR_tmR_array,"yes");
        my $FOLA_final = join(':',@FOLA_output);
        print $fh "FOLA-1\t$FOLA_final\ttmpR\tNA\n";
    }
} else {
    print $fh "FOLA-1\tFOLA_Extraction_Error\t---\tNA\n";
}
#################################################################


#################################################################
###Now Type the Fluoroquinolone Targets###
if ($miscR_Type{"PARC-12"}) {
    my @PARC_output;
    my @miscR_value = split(':',$miscR_Type{"PARC-12"});
    my $PARC_seq = freebayes_prior_fix($MISC_bam, $misc_DB,"10__PARC__PARC-12__10");
    my $PARC_aaSeq = sixFrame_Translate($PARC_seq,1);
    print "PARCgbs sequence: $PARC_seq || $PARC_aaSeq\n";
    my $PARC_pos2 = substr($PARC_aaSeq,1,1);
    my $PARC_pos6 = substr($PARC_aaSeq,5,1);
    my $PARC_pos14 = substr($PARC_aaSeq,13,1);
    if ($PARC_pos2 eq "F" || $PARC_pos2 eq "Y") {
        my $pos_output = "pos2=$PARC_pos2";
        push(@PARC_output,$pos_output);
        $isNotResistant = "no";
    }
    if ($PARC_pos6 eq "N" || $PARC_pos6 eq "Y") {
        my $pos_output = "pos6=$PARC_pos6";
        push(@PARC_output,$pos_output);
        $isNotResistant = "no";
    }
    if ($PARC_pos14 eq "D") {
        my $pos_output = "pos14=$PARC_pos14";
        push(@PARC_output,$pos_output);
        $isNotResistant = "no";
    }

    if (@PARC_output) {
        my @miscR_value = split(':',$miscR_Type{"PARC-12"});
        my $PARC_final = join(':',@PARC_output);
        print $fh "PARC-12\t$PARC_final\tFLQ\t$miscR_value[0]\n";
    } elsif ($PARC_aaSeq !~ /DSSIYDAMVRMSQNWKNREILVEMHGNNGSMDGDPPAAMRYTEARLSEIAGYLLQDIEK[K|N]TVPFAWNFDD/) {
        print $fh "PARC-12\timperfect\tPossible_FLQ_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"PARC-12"} = $PARC_seq;
    }
#    } elsif ($miscR_value[1] eq "yes") {
#        print $fh "PARC-12\timperfect\tPossible_FLQ_(Extract Seq)\t$miscR_value[0]\n";
#        $miscR_extract{"PARC-12"} = $PARC_seq;
#    }
}

if ($miscR_Type{"GYRA-1"}) {
    my @GYRA_output;
    my @miscR_value = split(':',$miscR_Type{"GYRA-1"});
    my $GYRA_seq = freebayes_prior_fix($MISC_bam, $misc_DB, "8__GYRA__GYRA-1__8");
    my $GYRA_aaSeq = sixFrame_Translate($GYRA_seq,1);
    print "GYRAgbs sequence: $GYRA_aaSeq\n";
    my $GYRA_pos11 = substr($GYRA_aaSeq,10,1);
    if ($GYRA_pos11 eq "F" || $GYRA_pos11 eq "Y") {
        my $pos_output = "pos11=$GYRA_pos11";
        push(@GYRA_output,$pos_output);
        $isNotResistant = "no";
    }
    if (@GYRA_output) {
        my @miscR_value = split(':',$miscR_Type{"GYRA-1"});
        my $GYRA_final = join(':',@GYRA_output);
        print $fh "GYRA-1\t$GYRA_final\tFLQ\t$miscR_value[0]\n";
    } elsif ($miscR_value[1] eq "yes") {
        print $fh "GYRA-1\timperfect\tPossible_FLQ_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"GYRA-1"} = $GYRA_seq;
    }
}
#################################################################


#################################################################
###Type the RPLD and RPLV Macrolide and Streptogramins targets###
###The RPLD-1 marker for MQD-1 resistance may contain mosiac regions and must be extracted from the genome assembly###
if ($miscR_Type{"RPLD-1"}) {
    my $RPLD1_ref = extractFastaByID("11__RPLD__RPLD-1__11",$misc_DB);
    `echo "$RPLD1_ref" > TEMP_RPLD1_Ref.fna`;
    system("LoTrac_target.pl -1 $fastq1 -2 $fastq2 -q TEMP_RPLD1_Ref.fna -S 2.2M -f -n $justName -L 0.8 -I 0.8");
    my @TEMP_RPLD1_extract = glob("EXTRACT*RPLD-1*.fasta");
    my $RPLD1_fileName = $TEMP_RPLD1_extract[0];
    my $RPLD1_fragName = `cat $RPLD1_fileName | grep ">" | tail -n1 | sed 's/>//g'`;
    chomp($RPLD1_fileName);
    chomp($RPLD1_fragName);
    my $RPLD1_seq = extractFastaByID($RPLD1_fragName,$RPLD1_fileName);
    my $RPLD1_aaSeq = sixFrame_Translate($RPLD1_seq,1);
    if ($RPLD1_aaSeq ne "KPWRQKGTGRAR") {
        print "RPLD1 INFO: RPLD1 name: $RPLD1_fragName | RPLD1 filename: $RPLD1_fileName | RPLD1 seq: $RPLD1_seq | RPLD aa seq: $RPLD1_aaSeq\n";
        my @miscR_value = split(':',$miscR_Type{"RPLD-1"});
        print $fh "RPLD-1\timperfect\tPossible_MQD_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"RPLD-1"} = $RPLD1_seq;
    }
}

if ($miscR_Type{"RPLD-2"}) {
    my $RPLD2_seq = freebayes_prior_fix($MISC_bam, $misc_DB, "12__RPLD__RPLD-2__12");
    my $RPLD2_aaSeq = sixFrame_Translate($RPLD2_seq,1);
    if ($RPLD2_aaSeq eq "DAVFGIEPNKSVVFDVI") {
        print "RPLD2 seq: $RPLD2_seq\n";
        my @miscR_value = split(':',$miscR_Type{"RPLD-2"});
        print $fh "RPLD-2\timperfect\tPossible_MQD-2_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"RPLD-2"} = $RPLD2_seq;
    }
}

if ($miscR_Type{"RPLV-1"}) {
    my $RPLV1_seq = freebayes_prior_fix($MISC_bam, $misc_DB, "13__RPLV__RPLV-1__13");
    my $RPLV1_aaSeq = sixFrame_Translate($RPLV1_seq,1);
    if ($RPLV1_aaSeq ne "KRTAHITVA") {
        my @miscR_value = split(':',$miscR_Type{"RPLV-1"});
        print $fh "RPLV-1\timperfect\tPossible_MQD-3_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"RPLV-1"} = $RPLV1_seq;
    }
}

if ($miscR_Type{"RPLV-2"}) {
    my $RPLV2_seq = freebayes_prior_fix($MISC_bam, $misc_DB, "14__RPLV__RPLV-2__14");
    my $RPLV2_aaSeq = sixFrame_Translate($RPLV2_seq,1);
    if ($RPLV2_aaSeq ne "PTMKRFRPRA") {
        print "RPLV2 seq: $RPLV2_seq\n";
        my @miscR_value = split(':',$miscR_Type{"RPLV-2"});
        print $fh "RPLV-2\timperfect\tPossible_MQD-4_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"RPLV-2"} = $RPLV2_seq;
    }
}
#################################################################


#################################################################
###Type the 23S ribosomal RNA MLS resistance target###
if ($miscR_Type{"23SWT1-1"}) {
    my @two3SWT1_output;
    my $two3SWT1_seq = freebayes_prior_fix($MISC_bam, $misc_DB, "1__23SWT1__23SWT1-1__1");
    my $two3SWT1_pos21 = substr($two3SWT1_seq,20,1);
    my $two3SWT1_pos20 = substr($two3SWT1_seq,19,1);
    if ($two3SWT1_pos21 eq "G") {
        my $pos_out = "pos21=G";
        push(@two3SWT1_output,$pos_out);
    } elsif ($two3SWT1_pos20 eq "C") {
        my $pos_out = "pos20=C";
        push(@two3SWT1_output,$pos_out);
    }
    my @miscR_value = split(':',$miscR_Type{"23SWT1-1"});
    if (@two3SWT1_output) {
        my $two3SWT1_final = join(':',@two3SWT1_output);
        print $fh "23SWT1-1\t$two3SWT1_final\tMLS\t$miscR_value[0]\n";
    } elsif ($miscR_value[1] eq "yes") {
        print $fh "23SWT1-1\timperfect\tPossible_MLS_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"23SWT1-1"} = $two3SWT1_seq;
    }
}

if ($miscR_Type{"23SWT2-1"}) {
    my @two3SWT2_output;
    my $two3SWT2_seq = freebayes_prior_fix($MISC_bam, $misc_DB, "2__23SWT2__23SWT2-1__2");
    my $two3SWT2_pos14 = substr($two3SWT2_seq,13,1);
    if ($two3SWT2_pos14 eq "G") {
        my $pos_out = "pos14=G";
        push(@two3SWT2_output,$pos_out);
    }
    my @miscR_value = split(':',$miscR_Type{"23SWT2-1"});
    if (@two3SWT2_output) {
        my $two3SWT2_final = join(':',@two3SWT2_output);
        print $fh "23SWT2-1\t$two3SWT2_final\tMLS\t$miscR_value[0]\n";
    } elsif ($miscR_value[1] eq "yes") {
        print $fh "23SWT2-1\timperfect\tPossible_MLS_(Extract Seq)\t$miscR_value[0]\n";
        $miscR_extract{"23SWT2-1"} = $two3SWT2_seq;
    }
}
#################################################################


#################################################################
###Type the FOLP-1 Sulfonamide resistance target###
my $target = "7__FOLP__FOLP-1__7";
my @TEMP_VANC_bam = glob("VANC_*\.sorted\.bam");
my $bamFile = $TEMP_VANC_bam[0];
(my $samFile = $bamFile) =~ s/\.bam/\.sam/g;
system("samtools view -h $bamFile > $samFile");
system("cat $samFile | grep -E \"^\@HD|^\@SQ.*$target|^\@PG\" > FOLP_target_seq.sam");
system("awk -F'\t' '\$3 == \"$target\" {print \$0}' $samFile >> FOLP_target_seq.sam");
system("samtools view -bS FOLP_target_seq.sam > FOLP_target_seq.bam");
system("samtools index FOLP_target_seq.bam FOLP_target_seq.bai");
$REF_seq = extractFastaByID("$target","$vanc_DB");
open(my $rf,'>',"FOLP_target_ref.fna");
print $rf "$REF_seq\n";
close $rf;
system("freebayes -q 20 -p 1 -f FOLP_target_ref.fna FOLP_target_seq.bam -v FOLP_target_seq.vcf");
open(MYINPUTFILE, "FOLP_target_seq.vcf");
my %srst2_seroT;
while(<MYINPUTFILE>) {
    my $line = $_;
    chomp($line);
    if ($line =~ /^7__FOLP__FOLP-1__7/) {
	my @FOLP_line = split('\t', $line);
	my $ref_allele = $FOLP_line[3];
	#my $ref_len = length($ref_allele);
	my $alt_allele = $FOLP_line[4];
	#my $alt_len = length($alt_allele);
	#my $lociDP = ($FOLP_line[7] =~ /DP=([0-9]+);/);
	#$FOLP_line[7] =~ /DPB=([0-9]+);/;
	$FOLP_line[7] =~ /DPB=(\d+\.?\d*);/;
	print "FOLP DP: $1 | ref allele: $ref_allele | alt allele: $alt_allele\n";
	if (length($ref_allele) != length($alt_allele)) {
	    push(@sxR_tmR_array,"yes");
	    print $fh "FOLP\tIndel=nt$FOLP_line[1]\tsxR\t$1\n";
	}
    }
}

###BASH!!!! FOLP-1 script###
#if [[ -e $(ls VANC*_Final.sorted.vcf) ]]
#then
#    cat VANC*_Final.sorted.vcf | grep "FOLP" > TEMP_vanc-FOLP.vcf
#    if [[ -s TEMP_vanc-FOLP.vcf ]]
#    then
#        while read snpLine
#        do
#            ref_allele=$(echo "$snpLine" | awk -F"\t" '{print $4}')
#            ref_length=$(echo ${#ref_allele})
#            alt_allele=$(echo "$snpLine" | awk -F"\t" '{print $5}')
#            alt_length=$(echo ${#alt_allele})
#            lociDP=$(echo "$snpLine" | awk -F"\t" '{print $8}' | sed 's/.*;DP=\([0-9]\+\);.*/\1/g')
#            lociPos=$(echo "$snpLine" | awk -F"\t" '{print $2}')
#
#            if [[ "$ref_length" -ne "$alt_length" ]]
#            then
#                sxR_tmRarray+=("yes")
#                printf "FOLP\tIndel=nt$lociPos\tsxR\t$lociDP\n" >> Misc_Resistance_results.txt
#            fi
#        done < TEMP_vanc-FOLP.vcf
#    fi
#fi
#################################################################


#################################################################
###Type the Co-trimoxazole resistance target (presence of tmR (FOLA-1) and sxR (FOLP-1))###
my $sxR_tmR_len = scalar(@sxR_tmR_array);
if ($sxR_tmR_len == 2) {
    print $fh "FOLA/FOLP\ttmpR+/sxR+\tcotR\tNA\n";
}
#################################################################





###Print out all of the extracted sequence to the 'Check_Target_Sequence.txt' file###
if (%miscR_extract) {
    open ( my $exFile_out, ">>", 'Check_Target_Sequence.txt' ) or die "Could not open file 'Check_Target_Sequence.txt': $!";
    print $exFile_out '#' x 56;
    print $exFile_out "--Misc. Resistance Extracted Sequence--";
    print $exFile_out '#' x 56;
    print $exFile_out "\n";
    foreach my $key (keys(%miscR_extract)) {
        print $exFile_out ">$key\n";
        print $exFile_out "$miscR_extract{$key}\n";
        print $exFile_out "\n";
}
print $exFile_out '#' x 151;
print $exFile_out "\n\n";
close $exFile_out;
}

###Delete TEMP files and close open files###
#unlink(TEMP*);
close $fh;
