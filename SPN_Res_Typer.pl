#!/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;
use File::Copy qw(copy);
use Env;
use lib $ENV{MODULESHOME}."/init";
use perl;

###MODULE LOAD###
#module load samtools/0.1.18
#module load bowtie2/2.1.0
#module load Python/2.7
#module load freebayes/0.9.21

sub checkOptions {
    my %opts;
    getopts('h1:2:d:r:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $ref_dir, $res_DB, $outDir, $outName);

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

    if($opts{d}) {
        $ref_dir = $opts{d};
        $ref_dir =~ s/\/$//g;
        if (-d $ref_dir) {
            print "Directory containing the  resistance reference sequences: $ref_dir\n";
        } else {
            print "The directory containing the resistance references is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/ref_dir).\n";
            help();
        }
    } else {
        print "The resistance reference sequence directory (including full path) has not been given.\n";
        help();
    }

    if($opts{r}) {
        $res_DB = "$ref_dir/$opts{r}";
        if ($res_DB) {
            print "The resistance reference sequence file: $opts{r}\n";
        } else {
            print "The resistance reference sequence file is not in the correct format or doesn't exist.\n";
            #print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The resistance reference sequence file has not been given.\n";
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
	chomp($outName);
        print "The output file name prefix: $outName\n";
    } else {
        $outName = `echo "$fastq1" | awk -F"/" '{print \$(NF)}' | sed 's/_S[0-9]\\+_L[0-9]\\+_R[0-9]\\+.*//g'`;
	chomp($outName);
        print "The default output file name prefix is: $outName";
    }

    return ($help, $fastq1, $fastq2, $ref_dir, $res_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
SPN_Res_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -d <ref directory: dir> -r <resistance seq: file> -o <output directory name: string> -n <output name prefix: string> [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -d   reference sequence directory (including full path)
    -r   resistance reference sequence file
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $ref_dir, $res_DB, $outDir, $outName) = checkOptions( @ARGV );





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
my $Res_output = "OUT_Res_Results.txt";
open(my $fh,'>',$Res_output) or die "Could not open file '$Res_output' $!";
my $BIN_res_out = "BIN_Res_Results.txt";
open(my $bh,'>',$BIN_res_out) or die "Could not open file '$BIN_res_out' $!";
my @Bin_Res_arr = (0) x 19;
#print $fh "Resistance_Group\tTarget\n";
#=pod
my $outNameRES = "RES_".$outName;
my $out_nameARG = "ARG_".$outName;
my $out_nameRESFI = "RESFI_".$outName;
my $out_nameFOLP = "FOLP_".$outName;
print "resistance db: $res_DB\n";
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $outNameRES --log --save_scores --min_coverage 99.9 --max_divergence 5 --gene_db $res_DB");
###Type ARG-ANNOT Resistance Genes###
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $out_nameARG --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db $ref_dir/ARGannot_r1.fasta");
###Type ResFinder Resistance Genes###
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $out_nameRESFI --log --save_scores --min_coverage 70 --max_divergence 30 --gene_db $ref_dir/ResFinder.fasta");
###Type FOLP Resistance Gene###
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $out_nameFOLP --log --save_scores --min_coverage 95.0 --max_divergence 15 --gene_db $ref_dir/SPN_FOLP_Gene-DB_Final.fasta");
#=cut

my @TEMP_RES_bam = glob("RES_*\.sorted\.bam");
my @TEMP_RES_fullgene = glob("RES_*__fullgenes__*__results\.txt");
my $RES_bam = $TEMP_RES_bam[0];
my $RES_full_name = $TEMP_RES_fullgene[0];
print "res bam is: $RES_bam || res full gene $RES_full_name\n";
(my $RES_vcf = $RES_bam) =~ s/\.bam/\.vcf/g;
(my $RES_bai = $RES_bam) =~ s/\.bam/\.bai/g;

my @TEMP_ARG_fullgene = glob("ARG_*__fullgenes__*__results\.txt");
my $ARG_full_name = $TEMP_ARG_fullgene[0];
my @TEMP_RESFI_fullgene = glob("RESFI_*__fullgenes__*__results\.txt");
my $RESFI_full_name = $TEMP_RESFI_fullgene[0];
my $merged_net = "ARG-RESFI_fullgenes_results.txt";
#copy $ARG_full_name, $merged_net;
if ($ARG_full_name) {
    system("tail -n+2 $ARG_full_name > $merged_net");
}
if ($RESFI_full_name) {
    system("tail -n+2 $RESFI_full_name >> $merged_net");
}

my %drugRes_Col = (
    "TET" => "neg",
    "EC" => "neg",
    "FQ" => "neg",
    "COT" => "neg",
    "OTHER" => "neg",
    );

my %Res_Targets = (
    "RPOB" => "neg",
    "FOLA" => "neg",
    "PARC" => "neg",
    "GYRA" => "neg",
    "RPLD1" => "neg",
    "RPLD2" => "neg",
    "RPLV1" => "neg",
    "RPLV2" => "neg",
    "R23S1" => "neg",
    "R23S2" => "neg",
    "FOLP" => "neg",
    "ERMB" => "neg",
    "ERMBUP" => "neg",
    "ERMBS" => "neg",
    "ERMTR" => "neg",
    "TET" => "neg",
    "CAT" => "neg",
    );

###Type the Presence/Absence Targets###
open(MYINPUTFILE, "$RES_full_name");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @miscR_fullgene;
    @miscR_fullgene = split('\t',$line);
    if ($miscR_fullgene[5] >= 10) {
        if ($miscR_fullgene[3] =~ m/(ERMTR|MEF)/) {
            if ($drugRes_Col{"EC"} eq "neg") {
                $drugRes_Col{"EC"} = $miscR_fullgene[2];
            } else {
                my $new_val = $drugRes_Col{"EC"}.":".$miscR_fullgene[2];
                $drugRes_Col{"EC"} = $new_val;
            }
        }
        if ($miscR_fullgene[3] =~ m/TET/) {
	    $Res_Targets{"TET"} = "pos";
            if ($drugRes_Col{"TET"} eq "neg") {
                $drugRes_Col{"TET"} = $miscR_fullgene[2];
            } else {
                my $new_val = $drugRes_Col{"TET"}.":".$miscR_fullgene[2];
                $drugRes_Col{"TET"} = $new_val;
            }
        }
        if ($miscR_fullgene[3] =~ m/CAT/) {
	    $Res_Targets{"CAT"} = "pos";
            if ($drugRes_Col{"OTHER"} eq "neg") {
                $drugRes_Col{"OTHER"} = $miscR_fullgene[2];
            } else {
                my $new_val = $drugRes_Col{"OTHER"}.":".$miscR_fullgene[2];
                $drugRes_Col{"OTHER"} = $new_val;
            }
        }

        if ($miscR_fullgene[3] =~ m/RPOB/) {
            $Res_Targets{"RPOB"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/FOLA/) {
            $Res_Targets{"FOLA"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/PARC/) {
            $Res_Targets{"PARC"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/GYRA/) {
            $Res_Targets{"GYRA"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/RPLD1/) {
            $Res_Targets{"RPLD1"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/RPLD2/) {
            $Res_Targets{"RPLD2"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/RPLV1/) {
            $Res_Targets{"RPLV1"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/RPLV2/) {
            $Res_Targets{"RPLV2"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/R23S1/) {
            $Res_Targets{"R23S1"} = "pos";
        } elsif ($miscR_fullgene[3] =~ m/R23S2/) {
            $Res_Targets{"R23S2"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/ERMB-1/) {
	    $Res_Targets{"ERMB"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/ERMBUP/ && $miscR_fullgene[6]) {
	    $Res_Targets{"ERMBUP"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/ERMBS/ && ! $miscR_fullgene[6]) {
	    $Res_Targets{"ERMBS"} = "pos";
	} elsif ($miscR_fullgene[3] =~ m/ERMTR/) {
	    $Res_Targets{"ERMTR"} = "pos";
	}
    }
}

###Print out ERMB targets###
if ($Res_Targets{"ERMB"} eq "pos" && $Res_Targets{"ERMTR"} eq "neg") {
    my $ERMB_out = "ERMB";
    if ($Res_Targets{"ERMBS"} eq "pos") {
        $ERMB_out = $ERMB_out.":ERMBS";
    }
    if ($Res_Targets{"ERMBUP"} eq "pos") {
        $ERMB_out = $ERMB_out.":ERMBUP";
    }

    if ($drugRes_Col{"EC"} eq "neg") {
        $drugRes_Col{"EC"} = "$ERMB_out";
    } else {
        my $new_val = $drugRes_Col{"EC"}.":$ERMB_out";
        $drugRes_Col{"EC"} = $new_val;
    }
}

while (my ($key, $val) = each %Res_Targets) {
    my @val_arr = split(':',$val);
    my @val_sort = sort(@val_arr);
    my $val_out = join(':',@val_sort);
    print "$key\t$val_out\n";
}
print "\n";
print Dumper(\%drugRes_Col);

###############################################################################################
###ARG-ANNOT and ResFinder Safety Net###
open(MYINPUTFILE, "$merged_net");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @miscR_fullgene;
    @miscR_fullgene = split('\t',$line);
    if ($miscR_fullgene[5] >= 10) {
        if ($miscR_fullgene[3] =~ m/ERM/i) {
            if ($Res_Targets{"ERMB"} eq "neg" && $Res_Targets{"ERMTR"} eq "neg") {
                if ($drugRes_Col{"EC"} eq "neg") {
                    $drugRes_Col{"EC"} = "ERM";
                } else {
                    my $new_val = $drugRes_Col{"EC"}.":ERM";
                    $drugRes_Col{"EC"} = $new_val;
                }
                $Res_Targets{"ERMB"} = "pos";
            }
        #} elsif ($miscR_fullgene[3] =~ m/LNU/i) {
        #    if ($Res_Targets{"LNUB"} eq "neg") {
        #        if ($drugRes_Col{"EC"} eq "neg") {
        #            $drugRes_Col{"EC"} = "LNU";
        #        } else {
        #            my $new_val = $drugRes_Col{"EC"}.":LNU";
        #            $drugRes_Col{"EC"} = $new_val;
        #        }
        #        $Res_Targets{"LNUB"} = "pos";
        #    }
        #} elsif ($miscR_fullgene[3] =~ m/LSA/i) {
        #    if ($Res_Targets{"LSA"} eq "neg") {
        #        if ($drugRes_Col{"EC"} eq "neg") {
        #            $drugRes_Col{"EC"} = "LSA";
        #        } else {
        #            my $new_val = $drugRes_Col{"EC"}.":LSA";
        #            $drugRes_Col{"EC"} = $new_val;
        #        }
        #        $Res_Targets{"LSA"} = "pos";
        #    }
        } elsif ($miscR_fullgene[3] =~ m/MEF/i) {#&& $Res_Targets{"MEF"} eq "neg") {
            if ($Res_Targets{"MEF"} eq "neg") {
                if ($drugRes_Col{"EC"} eq "neg") {
                    $drugRes_Col{"EC"} = "MEF";
                } else {
                    my $new_val = $drugRes_Col{"EC"}.":MEF";
                    $drugRes_Col{"EC"} = $new_val;
                }
                $Res_Targets{"MEF"} = "pos";
            }
        } elsif ($miscR_fullgene[3] =~ m/TET/i) { #&& $Res_Targets{"TET"} eq "neg") {
            if ($Res_Targets{"TET"} eq "neg") {
                if ($drugRes_Col{"TET"} eq "neg") {
                    $drugRes_Col{"TET"} = "TET";
                } else {
                    my $new_val = $drugRes_Col{"TET"}.":TET";
                    $drugRes_Col{"TET"} = $new_val;
                }
                $Res_Targets{"TET"} = "pos";
            }
        } elsif ($miscR_fullgene[3] =~ m/CAT/i) {#&& $Res_Targets{"CAT"} eq "neg") {
            if ($Res_Targets{"CAT"} eq "neg") {
                if ($drugRes_Col{"OTHER"} eq "neg") {
                    $drugRes_Col{"OTHER"} = "CAT";
                } else {
                    my $new_val = $drugRes_Col{"OTHER"}.":CAT";
                    $drugRes_Col{"OTHER"} = $new_val;
                }
                $Res_Targets{"CAT"} = "pos";
            }
        } else {
            if ($drugRes_Col{"OTHER"} eq "neg") {
                $drugRes_Col{"OTHER"} = $miscR_fullgene[2];
            } else {
                my $new_val = $drugRes_Col{"OTHER"}.":".$miscR_fullgene[2];
                $drugRes_Col{"OTHER"} = $new_val;
            }
        }
    }
}
###############################################################################################


###Now Type the Non-Presence/Absence Targets###
#my @sxR_tmR_array;
###############################################################################################
###Type the RPOB-1 Rifampicin targets###
if ($Res_Targets{"RPOB"} eq "pos") {
    my $RPOB_output;
    my $RPOB_seq = freebayes_prior_fix($RES_bam, $res_DB, "15__RPOB__RPOB-1__15");
    my $RPOB_aaSeq = sixFrame_Translate($RPOB_seq,1);
    my $RPOB_aaRef = "FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL";
    print "RPOB ref: $RPOB_aaRef || RPOB seq: $RPOB_aaSeq\n";
    if ($RPOB_aaSeq ne "FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL") {
        my $mask = $RPOB_aaSeq ^ $RPOB_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($RPOB_aaRef,$-[0],1), ' ', substr($RPOB_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
	    my $diff_element = substr($RPOB_aaRef,$-[0],1).($-[0]+1).substr($RPOB_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "RPOB seq: $RPOB_seq\n";
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[6] = $bin_out;
        my $RPOB_out = "RPOB1-".$diff_output;
        if ($drugRes_Col{"OTHER"} eq "neg") {
            $drugRes_Col{"OTHER"} = $RPOB_out;
        } else {
            my $new_val = $drugRes_Col{"OTHER"}.":".$RPOB_out;
            $drugRes_Col{"OTHER"} = $new_val;
        }
    }
}
###############################################################################################


###############################################################################################
###Type the FOLA-1 Trimethoprim targets###
###The FOLA-1 marker for tmR resistance may contain mosiac regions and must be extracted from the genome assembly###
my $REF_seq = extractFastaByID("7__FOLA__FOLA-1__7","$res_DB");
`echo "$REF_seq" > TEMP_FOLA_Ref.fna`;
module "unload perl/5.22.1";
module "load perl/5.16.1-MT";
system("LoTrac_target.pl -1 $fastq1 -2 $fastq2 -q TEMP_FOLA_Ref.fna -S 2.2M -f -n $outName -o $outDir");
module "unload perl/5.16.1-MT";
module "load perl/5.22.1";
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
    my $FOLA_aaRef = "QDVQSVLDWYQDQEKNLYII";
    print "FOLA Seq: $FOLA_aaRef || $FOLA_aaSeq || $FOLA_seq\n";
    if ($FOLA_aaSeq !~ /QDVQSVL[D|G]WYQ[D|V|A]QEKNLYII/) {
	#push(@sxR_tmR_array,"yes");
        my $mask = $FOLA_aaSeq ^ $FOLA_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($FOLA_aaRef,$-[0],1), ' ', substr($FOLA_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($FOLA_aaRef,$-[0],1).($-[0]+1).substr($FOLA_aaSeq,$-[0],1);
	    if ($diff_element !~ /D12[A|V]/ && $diff_element ne "D8G") {
		push(@seq_diffs,$diff_element);
	    }
        }
        print "FOLA seq: $FOLA_seq\n";
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[6] = $bin_out;
        my $FOLA_out = "FOLA-".$diff_output;
        if ($drugRes_Col{"COT"} eq "neg") {
            $drugRes_Col{"COT"} = $FOLA_out;
        } else {
            my $new_val = $drugRes_Col{"COT"}.":".$FOLA_out;
            $drugRes_Col{"COT"} = $new_val;
        }
    }
}
###############################################################################################


###############################################################################################
###Type the PARC and GYRA FLQ Targets###
if ($Res_Targets{"PARC"} eq "pos") {
    my @PARC_output;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    my $PARC_seq = freebayes_prior_fix($RES_bam, $res_DB,"10__PARC__PARC-12__10");
    my $PARC_aaSeq = sixFrame_Translate($PARC_seq,1);
    my $PARC_aaRef = "DSSIYDAMVRMSQNWKNREILVEMHGNNGSMDGDPPAAMRYTEARLSEIAGYLLQDIEKKTVPFAWNFDD";
    print "PARCgbs sequence: $PARC_seq || $PARC_aaSeq\n";
    if ($PARC_aaSeq !~ /DSSIYDAMVRMSQNWKN[R|C]EILVEMHGNNGSMDGDPPAAMRYTEARLSEIA[G|D]YLLQDIEK[K|N]TVPFAWNFDD/) {
        my $mask = $PARC_aaSeq ^ $PARC_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($PARC_aaRef,$-[0],1), ' ', substr($PARC_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($PARC_aaRef,$-[0],1).($-[0]+1).substr($PARC_aaSeq,$-[0],1);
            if ($diff_element ne "G51D" && $diff_element ne "K60N" && $diff_element ne "R18C") {
		push(@seq_diffs,$diff_element);
	    }
        }
        print "PARC seq: $PARC_seq\n";
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[14] = $bin_out;
        my $PARC_out = "PARC-".$diff_output;
        if ($drugRes_Col{"FQ"} eq "neg") {
            $drugRes_Col{"FQ"} = $PARC_out;
        } else {
            my $new_val = $drugRes_Col{"FQ"}.":".$PARC_out;
            $drugRes_Col{"FQ"} = $new_val;
        }
    }
}

if ($Res_Targets{"GYRA"} eq "pos") {
    my @GYRA_output;
    #my @miscR_value = split(':',$miscR_Type{"PARCGBS-1"});
    #my $PARC_seq = extractFastaByID("7__PARCGBS__PARCGBS-1__7","TEMP_miscR_consensus.fna");
    my $GYRA_seq = freebayes_prior_fix($RES_bam, $res_DB,"8__GYRA__GYRA-1__8");
    my $GYRA_aaSeq = sixFrame_Translate($GYRA_seq,1);
    my $GYRA_aaRef = "VMGKYHPHGDSSIYEAMVRMAQWWSY";
    print "GYRA sequence: $GYRA_seq || $GYRA_aaSeq\n";
    if ($GYRA_aaSeq ne "VMGKYHPHGDSSIYEAMVRMAQWWSY") {
        my $mask = $GYRA_aaSeq ^ $GYRA_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($GYRA_aaRef,$-[0],1), ' ', substr($GYRA_aaSeq,$-[0],1), ' ', $-[0], "\n";
            #my $diff_element = "pos".($-[0]+1).":".substr($RPOB_aaRef,$-[0],1)."->".substr($RPOB_aaSeq,$-[0],1);
            my $diff_element = substr($GYRA_aaRef,$-[0],1).($-[0]+1).substr($GYRA_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        print "GYRA seq: $GYRA_seq\n";
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[10] = $bin_out;
        my $GYRA_out = "GYRA-".$diff_output;
        if ($drugRes_Col{"FQ"} eq "neg") {
            $drugRes_Col{"FQ"} = $GYRA_out;
        } else {
            my $new_val = $drugRes_Col{"FQ"}.":".$GYRA_out;
            $drugRes_Col{"FQ"} = $new_val;
        }
    }
}
###############################################################################################


###############################################################################################
###Type the RPLD and RPLV Macrolide and Streptogramins targets###
###The RPLD-1 marker for MQD-1 resistance may contain mosiac regions and must be extracted from the genome assembly###
if ($Res_Targets{"RPLD1"} eq "pos") {
    my $RPLD1_ref = extractFastaByID("11__RPLD1__RPLD1-1__11",$res_DB);
    `echo "$RPLD1_ref" > TEMP_RPLD1_Ref.fna`;
    module "unload perl/5.22.1";
    module "load perl/5.16.1-MT";
    system("LoTrac_target.pl -1 $fastq1 -2 $fastq2 -q TEMP_RPLD1_Ref.fna -S 2.2M -f -n $outName -L 0.8 -I 0.8");
    module "unload perl/5.16.1-MT";
    module "load perl/5.22.1";
    my $RPLD1_file = glob("EXTRACT*RPLD1-1*.fasta");
    my $RPLD1_error = glob("ERROR*RPLD1-1*.fasta");
    my @RPLD1_output;
    #my $RPLD1_fileName = $TEMP_RPLD1_extract[0];
    if ($RPLD1_file && ! $RPLD1_error) {
	my $RPLD1_seq_hedr = `cat $RPLD1_file | grep ">" | tail -n1 | sed 's/>//g'`;
	chomp($RPLD1_seq_hedr);
	my $RPLD1_seq = extractFastaByID($RPLD1_seq_hedr,$RPLD1_file);
	chomp($RPLD1_seq);
	my $RPLD1_aaSeq = sixFrame_Translate($RPLD1_seq,1);
	my $RPLD1_aaRef = "KPWRQKGTGRAR";
	print "RPLD1 Seq: $RPLD1_aaRef || $RPLD1_aaSeq || $RPLD1_seq\n";
	if ($RPLD1_aaSeq ne "KPWRQKGTGRAR") {
	    my $mask = $RPLD1_aaSeq ^ $RPLD1_aaRef;
	    my @seq_diffs;
	    while ($mask =~ /[^\0]/g) {
		print substr($RPLD1_aaRef,$-[0],1), ' ', substr($RPLD1_aaSeq,$-[0],1), ' ', $-[0], "\n";
		my $diff_element = substr($RPLD1_aaRef,$-[0],1).($-[0]+1).substr($RPLD1_aaSeq,$-[0],1);
		push(@seq_diffs,$diff_element);
	    }
	    print "RPLD1 seq: $RPLD1_seq\n";
	    my $diff_output = join(',',@seq_diffs);
            #my $bin_out = join(':',@seq_diffs);
            #$Bin_Res_arr[6] = $bin_out;
	    my $RPLD1_out = "RPLD1-".$diff_output;
	    if ($drugRes_Col{"EC"} eq "neg") {
		$drugRes_Col{"EC"} = $RPLD1_out;
	    } else {
		my $new_val = $drugRes_Col{"EC"}.":".$RPLD1_out;
		$drugRes_Col{"EC"} = $new_val;
	    }
	}
    }
}

if ($Res_Targets{"RPLD2"} eq "pos") {
    my $RPLD2_seq = freebayes_prior_fix($RES_bam, $res_DB, "12__RPLD2__RPLD2-1__12");
    my $RPLD2_aaSeq = sixFrame_Translate($RPLD2_seq,1);
    my $RPLD2_aaRef = "DAVFGIEPNESVVFDVI";
    if ($RPLD2_aaSeq !~ /DAVFGIEPN[E|K]SV[V|E]FDVI/) {
        print "RPLD2 seq: $RPLD2_seq\n";
	my $mask = $RPLD2_aaSeq ^ $RPLD2_aaRef;
	my @seq_diffs;
	while ($mask =~ /[^\0]/g) {
	    print substr($RPLD2_aaRef,$-[0],1), ' ', substr($RPLD2_aaSeq,$-[0],1), ' ', $-[0], "\n";
	    my $diff_element = substr($RPLD2_aaRef,$-[0],1).($-[0]+1).substr($RPLD2_aaSeq,$-[0],1);
	    if ($diff_element ne "E10K" && $diff_element ne "V13E") {
		push(@seq_diffs,$diff_element);
	    }
        }
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[6] = $bin_out;
	my $RPLD2_out = "RPLD2-".$diff_output;
	if ($drugRes_Col{"EC"} eq "neg") {
	    $drugRes_Col{"EC"} = $RPLD2_out;
	} else {
	    my $new_val = $drugRes_Col{"EC"}.":".$RPLD2_out;
	    $drugRes_Col{"EC"} = $new_val;
	}
    }
}

if ($Res_Targets{"RPLV1"} eq "pos") {
    my $RPLV1_seq = freebayes_prior_fix($RES_bam, $res_DB, "13__RPLV1__RPLV1-1__13");
    my $RPLV1_aaSeq = sixFrame_Translate($RPLV1_seq,1);
    my $RPLV1_aaRef = "KRTAHITVA";
    if ($RPLV1_aaSeq ne "KRTAHITVA") {
	print "RPLV1 seq: $RPLV1_seq\n";
	my $mask = $RPLV1_aaSeq ^ $RPLV1_aaRef;
	my @seq_diffs;
	while ($mask =~ /[^\0]/g) {
	    print substr($RPLV1_aaRef,$-[0],1), ' ', substr($RPLV1_aaSeq,$-[0],1), ' ', $-[0], "\n";
	    my $diff_element = substr($RPLV1_aaRef,$-[0],1).($-[0]+1).substr($RPLV1_aaSeq,$-[0],1);
	    push(@seq_diffs,$diff_element);
	}
	my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[6] = $bin_out;
	my $RPLV1_out = "RPLV1-".$diff_output;
	if ($drugRes_Col{"EC"} eq "neg") {
	    $drugRes_Col{"EC"} = $RPLV1_out;
	} else {
	    my $new_val = $drugRes_Col{"EC"}.":".$RPLV1_out;
	    $drugRes_Col{"EC"} = $new_val;
	}
    }
}

if ($Res_Targets{"RPLV2"} eq "pos") {
    my $RPLV2_seq = freebayes_prior_fix($RES_bam, $res_DB, "14__RPLV2__RPLV2-1__14");
    my $RPLV2_aaSeq = sixFrame_Translate($RPLV2_seq,1);
    my $RPLV2_aaRef = "PTMKRFRPRA";
    if ($RPLV2_aaSeq ne "PTMKRFRPRA") {
        print "RPLV2 seq: $RPLV2_seq\n";
        my $mask = $RPLV2_aaSeq ^ $RPLV2_aaRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($RPLV2_aaRef,$-[0],1), ' ', substr($RPLV2_aaSeq,$-[0],1), ' ', $-[0], "\n";
            my $diff_element = substr($RPLV2_aaRef,$-[0],1).($-[0]+1).substr($RPLV2_aaSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[6] = $bin_out;
        my $RPLV2_out = "RPLV2-".$diff_output;
        if ($drugRes_Col{"EC"} eq "neg") {
            $drugRes_Col{"EC"} = $RPLV2_out;
        } else {
            my $new_val = $drugRes_Col{"EC"}.":".$RPLV2_out;
            $drugRes_Col{"EC"} = $new_val;
        }
    }
}
###############################################################################################


###############################################################################################
###Type the 23S ribosomal RNA MLS resistance target###
if ($Res_Targets{"R23S1"} eq "pos") {
    my $R23S1_ntFNA = freebayes_prior_fix($RES_bam, $res_DB, "1__R23S1__R23S1-1__1");
    my @R23S1_ntArr = split(/\n/,$R23S1_ntFNA);
    my $R23S1_ntSeq = $R23S1_ntArr[1];
    my $R23S1_ntRef = "GTTACCCGCGACAGGACGGAAAGACCCCATGGAG";
    print "R23S1 seq: $R23S1_ntSeq || R23S1 ref: $R23S1_ntRef\n";
    if ($R23S1_ntSeq ne "GTTACCCGCGACAGGACGGAAAGACCCCATGGAG") {
        print "R23S1 seq: $R23S1_ntSeq || R23S1 ref: $R23S1_ntRef\n";
        my $mask = $R23S1_ntSeq ^ $R23S1_ntRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($R23S1_ntRef,$-[0],1), ' ', substr($R23S1_ntSeq,$-[0],1), ' ', $-[0], "\n";
            my $diff_element = substr($R23S1_ntRef,$-[0],1).($-[0]+1).substr($R23S1_ntSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[6] = $bin_out;
        my $R23S1_out = "R23S1-".$diff_output;
        if ($drugRes_Col{"EC"} eq "neg") {
            $drugRes_Col{"EC"} = $R23S1_out;
        } else {
            my $new_val = $drugRes_Col{"EC"}.":".$R23S1_out;
            $drugRes_Col{"EC"} = $new_val;
        }
    }
}

if ($Res_Targets{"R23S2"} eq "pos") {
    my $R23S2_ntFNA = freebayes_prior_fix($RES_bam, $res_DB, "2__R23S2__R23S2-1__2");
    my @R23S2_ntArr = split(/\n/,$R23S2_ntFNA);
    my $R23S2_ntSeq = $R23S2_ntArr[1];
    my $R23S2_ntRef = "AGACAGTTCGGTCCCTATCCGTCGCGGGCG";
    if ($R23S2_ntSeq ne "AGACAGTTCGGTCCCTATCCGTCGCGGGCG") {
        print "R23S2 seq: $R23S2_ntSeq || R23S2 ref: $R23S2_ntRef\n";
        my $mask = $R23S2_ntSeq ^ $R23S2_ntRef;
        my @seq_diffs;
        while ($mask =~ /[^\0]/g) {
            print substr($R23S2_ntRef,$-[0],1), ' ', substr($R23S2_ntSeq,$-[0],1), ' ', $-[0], "\n";
            my $diff_element = substr($R23S2_ntRef,$-[0],1).($-[0]+1).substr($R23S2_ntSeq,$-[0],1);
            push(@seq_diffs,$diff_element);
        }
        my $diff_output = join(',',@seq_diffs);
        #my $bin_out = join(':',@seq_diffs);
        #$Bin_Res_arr[6] = $bin_out;
        my $R23S2_out = "R23S2-".$diff_output;
        if ($drugRes_Col{"EC"} eq "neg") {
            $drugRes_Col{"EC"} = $R23S2_out;
        } else {
            my $new_val = $drugRes_Col{"EC"}.":".$R23S2_out;
            $drugRes_Col{"EC"} = $new_val;
        }
    }
}
###############################################################################################


###############################################################################################
###Type the FOLP-1 Sulfonamide resistance target###
my $target = "1__FOLP__FOLP-1__1";
my @TEMP_FOLP_bam = glob("FOLP_*\.sorted\.bam");
my $bamFile = $TEMP_FOLP_bam[0];
(my $samFile = $bamFile) =~ s/\.bam/\.sam/g;
system("samtools view -h $bamFile > $samFile");
system("cat $samFile | grep -E \"^\@HD|^\@SQ.*$target|^\@PG\" > FOLP_target_seq.sam");
system("awk -F'\t' '\$3 == \"$target\" {print \$0}' $samFile >> FOLP_target_seq.sam");
system("samtools view -bS FOLP_target_seq.sam > FOLP_target_seq.bam");
system("samtools index FOLP_target_seq.bam FOLP_target_seq.bai");
$REF_seq = extractFastaByID("$target","$ref_dir/SPN_FOLP_Gene-DB_Final.fasta");
open(my $rf,'>',"FOLP_target_ref.fna");
print $rf "$REF_seq\n";
close $rf;
system("freebayes -q 20 -p 1 -f FOLP_target_ref.fna FOLP_target_seq.bam -v FOLP_target_seq.vcf");
open(MYINPUTFILE, "FOLP_target_seq.vcf");
my %srst2_seroT;
while(<MYINPUTFILE>) {
    my $line = $_;
    chomp($line);
    if ($line =~ /^1__FOLP__FOLP-1__1/) {
        my @FOLP_line = split('\t', $line);
        my $ref_allele = $FOLP_line[3];
        #my $ref_len = length($ref_allele);
        my $alt_allele = $FOLP_line[4];
        #my $alt_len = length($alt_allele);
        #my $lociDP = ($FOLP_line[7] =~ /DP=([0-9]+);/);
        #$FOLP_line[7] =~ /DPB=([0-9]+);/;
        $FOLP_line[7] =~ /DPB=(\d+\.?\d*);/;
        print "FOLP DP: $1 | ref allele: $ref_allele | alt allele: $alt_allele\n";
	my $FOLP_dp = $1;
	my $FOLP_loc = $FOLP_line[1];
        if (length($ref_allele) != length($alt_allele) && $FOLP_dp >= 2) {
	    my $FOLP_out;
            if (length($ref_allele) > length($alt_allele)) {
		$FOLP_out = "FOLP_".$FOLP_loc."-del";
	    } else {
		$FOLP_out = "FOLP_".$FOLP_loc."-ins";
            }
	    if ($drugRes_Col{"COT"} eq "neg") {
		$drugRes_Col{"COT"} = $FOLP_out;
	    } else {
		my $new_val = $drugRes_Col{"COT"}.":".$FOLP_out;
		$drugRes_Col{"COT"} = $new_val;
	    }
        }
    }
}
###############################################################################################


###############################################################################################
###Print Drug Resistance Output###
while (my ($key, $val) = each %drugRes_Col) {
    my @val_arr = split(':',$val);
    #print "@val_arr\n";
    my @val_sort = sort { "\L$a" cmp "\L$b" } @val_arr;
    #print "@val_sort\n";
    my $val_out = join(':',@val_sort);
    print "$key\t$val_out\n";
    $drugRes_Col{"$key"} = $val_out;
    #print $fh "$key\t$val_out\n";
}

print $fh "EC\t$drugRes_Col{'EC'}\n";
print $fh "COT\t$drugRes_Col{'COT'}\n";
print $fh "TET\t$drugRes_Col{'TET'}\n";
print $fh "FQ\t$drugRes_Col{'FQ'}\n";
print $fh "OTHER\t$drugRes_Col{'OTHER'}\n";
###############################################################################################
