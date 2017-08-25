#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

###MODULE LOAD###
#module load perl/5.12.3
#module load ncbi-blast+/2.2.29
#module load BEDTools/2.17.0
#module load Python/2.7

sub checkOptions {
    my %opts;
    getopts('h1:2:r:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $sero_DB, $outDir, $outName);

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
        $sero_DB = $opts{r};
        if (-e $sero_DB) {
            print "The serotype reference database sequence: $sero_DB\n";
        } else {
            print "The serotype reference sequence location is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The serotype reference sequence location (including full path) has not been given.\n";
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

    return ($help, $fastq1, $fastq2, $sero_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
GBS_serotyper.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <reference databases directory: file path> -o <output directory name: string> -n <output name prefix: string>  [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   serotype reference sequence directory (including full path)
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $sero_DB, $outDir, $outName) = checkOptions( @ARGV );




###Subroutines###
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


##Start Doing Stuff##
chdir "$outDir";
my $serotype_output = "OUT_SeroType_Results.txt";
open(my $fh,'>',$serotype_output) or die "Could not open file '$serotype_output' $!";
my $serotype_extract = "Serotype_Extraction_Sequence.fna";
open(my $se,'>',$serotype_extract) or die "Could not open file '$serotype_extract' $!";
print $fh "Target\tDiffs\tAvgDepth\tSerotype\n";

###Detect GAS serotype sequence###
my $sero_outName = "SERO_$outName";
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $sero_outName --log --save_scores --min_coverage 99.9 --max_divergence 5 --gene_db $sero_DB");

###mpileup the 'SERO_.*.sorted.bam and create the called variants file with freebayes.
opendir(DIR, ".") or die "Couldn't open directory for reading: $!";
my @seroTargets = grep (/SERO_.*__fullgenes__.*\.txt/,readdir(DIR));
closedir(DIR);
my @seroT_bam = glob("SERO_*\.sorted\.bam");
my $sero_bam = $seroT_bam[0];

my $seroT_output = $seroTargets[0];
#print "\nseroT_output is: $seroT_output\n";
open(MYINPUTFILE, "$seroT_output");
#my @srst2ARR;
my %srst2;
while(<MYINPUTFILE>) {
    my $line = $_;
    chomp($line);
    my @seroT_line = split('\t', $line);
    if ($seroT_line[2] eq "gene") {
        next;
    } else {
        if ($seroT_line[5] > 10) {
            #my $newLine = "$seroT_line[2]:$seroT_line[3]:$seroT_line[4]:$seroT_line[5]:$seroT_line[6]:$seroT_line[7]";
            #$srst2_seroT{$seroT_line[2]} = $newLine;
	    my @newLine = ("$seroT_line[2]", "$seroT_line[3]", "$seroT_line[4]", "$seroT_line[5]", "$seroT_line[6]", "$seroT_line[7]");
	    #push(@srst2ARR,\@newLine);
	    $srst2{$seroT_line[2]} = \@newLine;
        }
    }
}
#print Dumper(\%srst2_seroT);
#print Dumper(\%srst2HSH);

my %singlTarg95 = ( 
    'WZY1' => '1', 'WZY2' => '2', 'WCHE3' => '3', 'WZY4' => '4', 'WZY5' => '5', 'WZY8' => '8', 'WCRG10A' => '10A', 'GTF10F' => '10F', 'WZY11A' => '11A', 'WZY11B:11C' => '11B/C', 
    'WCII12A:12B:46' => '12A/B:46', 'WCII12F' => '12F', 'WZY13' => '13', 'WZY14' => '14', 'WZY16F' => '16F', 'GTF17F' => '17F', 'WZY20' => '20', 'WZY21' => '21', 
    'WZY23A' => '23A', 'WZY23B' => '23B', 'WZY23F' => '23F', 'WZY24F:24A:24B' => '24A/B/F', 'WCYE25A:25F' => '25A/F', 'WZY28A' => '28A', 'WZY31' => '31', 'TTS' => '37', 
    'WZY34' => '34', 'WZY35A' => '35A', 'WZY35C:42' => '35C/42', 'WCRO35F' => '35F', 'WHAI47F' => '47F', 'WCYV38' => '38', 'RRGA' => 'PI-1', 'PITB' => 'PI-2',
);

my %singlTarg100 = ( 
    'WZY18A' => '18A', 'WZY18F' => '18F', 'WZY19A' => '19A', 'WZY19AVAR' => '19A', 'WZY19F' => '19F', 'WZY19FVAR' => '19F',
);

foreach my $k (keys %srst2) {
    #print "TARGET $k: $srst2HSH{$k}|| @{$srst2HSH{$k}}[0]\n";
    my @targARR = @{$srst2{$k}};
    my $target = $targARR[0];
    my $diff = $targARR[4];
    my $depth = $targARR[3];
    print "Target: $targARR[0] | Depth: $targARR[3] | Diff: $targARR[4]\n";

    if (exists $singlTarg95{$target}) {
        print "Found Single Target 95% Match!! || $target\t$diff\t$depth\n";
        if (! $diff) {
            print "perfect match\n";
            $diff = "-";
        }
        print $fh "$target\t$diff\t$depth\t$singlTarg95{$target}\n";
    } elsif (exists $singlTarg100{$target} && ! $diff) {
        print "Found Single Target 100% Match || $target\t$diff\t$depth\n";
        print $fh "$target\t-\t$depth\t$singlTarg100{$target}\n";
    }
}

if (exists $srst2{WZY35B}) {
    if (exists $srst2{WCIG35B} && @{$srst2{WCIG35B}}[4]) {
	#print "Found 35B with imperfect WCIG35B: @{$srst2{WCIG35B}}[0]\n";
	print $fh "@{$srst2{WZY35B}}[0]:@{$srst2{WCIG35B}}[0]\t-:@{$srst2{WCIG35B}}[4]\t@{$srst2{WZY35B}}[3]:@{$srst2{WCIG35B}}[3]\t35B:35D\n";
    } else {
	#print "Found 35B with identical SCIG35B\n";
	print $fh "@{$srst2{WZY35B}}[0]\t-\t@{$srst2{WZY35B}}[3]\t35B\n";
    }
}

if (exists $srst2{WCIP6AC} && ! @{$srst2{WCIP6AC}}[4]) {
    if (exists $srst2{WCIN6AB}) {
	if (! @{$srst2{WCIN6AB}}[4]) {
	    print "Found 6A Match || @{$srst2{WCIN6AB}}[0]: None\n";
	    print $fh "@{$srst2{WCIP6AC}}[0]:@{$srst2{WCIN6AB}}[0]\t-:-\t@{$srst2{WCIP6AC}}[3]:@{$srst2{WCIN6AB}}[3]\t6A\n";
	} else {
            print "Found 6A Match || @{$srst2{WCIN6AB}}[0]: @{$srst2{WCIN6AB}}[4]\n";
            print $fh "@{$srst2{WCIP6AC}}[0]:@{$srst2{WCIN6AB}}[0]\t-:@{$srst2{WCIN6AB}}[4]\t@{$srst2{WCIP6AC}}[3]:@{$srst2{WCIN6AB}}[3]\t6A\n";
	}
    }
    elsif (exists $srst2{WCIN6CD}) {
	if (! @{$srst2{WCIN6CD}}[4]) {
	    print "Found 6C Match || @{$srst2{WCIN6CD}}[0]: None\n";
	    print $fh "@{$srst2{WCIP6AC}}[0]:@{$srst2{WCIN6CD}}[0]\t-:-\t@{$srst2{WCIP6AC}}[3]:@{$srst2{WCIN6CD}}[3]\t6C\n";
	} else {
	    print "Found 6C Match || @{$srst2{WCIN6CD}}[0]: @{$srst2{WCIN6CD}}[4]\n";
	    print $fh "@{$srst2{WCIP6AC}}[0]:@{$srst2{WCIN6CD}}[0]\t-:@{$srst2{WCIN6CD}}[4]\t@{$srst2{WCIP6AC}}[3]:@{$srst2{WCIN6CD}}[3]\t6C\n";
	}
    }
}

if (exists $srst2{WCIP6BD} && ! @{$srst2{WCIP6BD}}[4]) {
    if (exists $srst2{WCIN6AB}) {
        if (! @{$srst2{WCIN6AB}}[4]) {
            print "Found 6B Match || @{$srst2{WCIN6AB}}[0]: None\n";
            print $fh "@{$srst2{WCIP6BD}}[0]:@{$srst2{WCIN6AB}}[0]\t-:-\t@{$srst2{WCIP6BD}}[3]:@{$srst2{WCIN6AB}}[3]\t6B\n";
        } else {
            print "Found 6B Match || @{$srst2{WCIN6AB}}[0]: @{$srst2{WCIN6AB}}[4]\n";
            print $fh "@{$srst2{WCIP6BD}}[0]:@{$srst2{WCIN6AB}}[0]\t-:@{$srst2{WCIN6AB}}[4]\t@{$srst2{WCIP6BD}}[3]:@{$srst2{WCIN6AB}}[3]\t6B\n";
        }
    }
    elsif (exists $srst2{WCIN6CD}) {
        if (! @{$srst2{WCIN6CD}}[4]) {
            print "Found 6B Match || @{$srst2{WCIN6CD}}[0]: None\n";
            print $fh "@{$srst2{WCIP6AC}}[0]:@{$srst2{WCIN6CD}}[0]\t-:-\t@{$srst2{WCIP6AC}}[3]:@{$srst2{WCIN6CD}}[3]\t6D\n";
        } else {
            print "Found 6D Match || @{$srst2{WCIN6CD}}[0]: @{$srst2{WCIN6CD}}[4]\n";
            print $fh "@{$srst2{WCIP6AC}}[0]:@{$srst2{WCIN6CD}}[0]\t-:@{$srst2{WCIN6CD}}[4]\t@{$srst2{WCIP6AC}}[3]:@{$srst2{WCIN6CD}}[3]\t6D\n";
        }
    }
}

if (exists $srst2{WZY7C}) {
    if (exists $srst2{WCHF7C}) {
	print $fh "@{$srst2{WZY7C}}[0]:@{$srst2{WCHF7C}}[0]\tND:ND\t@{$srst2{WZY7C}}[3]:@{$srst2{WCHF7C}}[3]\t7C\n";
    } else {
	print $fh "@{$srst2{WZY7C}}[0]:@{$srst2{WCHF7C}}[0]\tND:ND\t@{$srst2{WZY7C}}[3]:@{$srst2{WCHF7C}}[3]\t7B\n";
    }
}

if (exists $srst2{WZY7F}) {
    if ($srst2{WCWD7F} && ! @{$srst2{WCWD7F}}[4]) {
	print $fh "@{$srst2{WZY7F}}[0]:@{$srst2{WCWD7F}}[0]\tND:-\t@{$srst2{WZY7F}}[3]:@{$srst2{WCWD7F}}[3]\t7F\n";
    } elsif ($srst2{WCWD7F} && @{$srst2{WCWD7F}}[4]) {
	print $fh "@{$srst2{WZY7F}}[0]:@{$srst2{WCWD7F}}[0]\tND:@{$srst2{WCWD7F}}[4]\t@{$srst2{WZY7F}}[3]:@{$srst2{WCWD7F}}[3]\t7A:7F\n";
	my $sero_seq = freebayes_prior_fix($sero_bam, $sero_DB,"35__WCWD7F__WCWD7F-1__45");
	print $se "$sero_seq\n";
    }
}

if (exists $srst2{WZY9N}) {
    if (exists $srst2{WCJA9L} && ! @{$srst2{WCJA9L}}[4]) {
        print $fh "@{$srst2{WZY9N}}[0]:@{$srst2{WCJA9L}}[0]\tND:-\t@{$srst2{WZY9N}}[3]:@{$srst2{WCJA9L}}[3]\t9L\n";
    } elsif (exists $srst2{WCJA9N} && ! @{$srst2{WCJA9N}}[4]) {
        print $fh "@{$srst2{WZY9N}}[0]:@{$srst2{WCJA9N}}[0]\tND:-\t@{$srst2{WZY9N}}[3]:@{$srst2{WCJA9N}}[3]\t9N\n";
    }
}

if (exists $srst2{WZY9V}) {
    if ($srst2{WCJE9V} && ! @{$srst2{WCJE9V}}[4]) {
        print $fh "@{$srst2{WZY9V}}[0]:@{$srst2{WCJE9V}}[0]\tND:-\t@{$srst2{WZY9V}}[3]:@{$srst2{WCJE9V}}[3]\t9V\n";
    } elsif ($srst2{WCJE9V} && @{$srst2{WCJE9V}}[4]) {
        print $fh "@{$srst2{WZY9V}}[0]:@{$srst2{WCJE9V}}[0]\tND:@{$srst2{WCJE9V}}[4]\t@{$srst2{WZY9V}}[3]:@{$srst2{WCJE9V}}[3]\t9A:9V\n";
        my $sero_seq = freebayes_prior_fix($sero_bam, $sero_DB,"40__WCJE9V__WCJE9V-1__51");
        print $se "$sero_seq\n";
    }
}

if (exists $srst2{WZY15A}) {
    if ($srst2{WCIZ15F} && ! @{$srst2{WCIZ15F}}[4]) {
	print $fh "@{$srst2{WZY15A}}[0]:@{$srst2{WCIZ15F}}[0]\tND:-\t@{$srst2{WZY15A}}[3]:@{$srst2{WCIZ15F}}[3]\t15F\n";
    } elsif (! $srst2{WCIZ15F}) {
	print $fh "@{$srst2{WZY15A}}[0]\tND\t@{$srst2{WZY15A}}[3]\t15A\n";
    }
}

if (exists $srst2{WZY15B}) {
    if ($srst2{WCIZ15B} && ! @{$srst2{WCIZ15B}}[4]) {
	print $fh "@{$srst2{WZY15B}}[0]:@{$srst2{WCIZ15B}}[0]\tND:-\t@{$srst2{WZY15B}}[3]:@{$srst2{WCIZ15B}}[3]\t15B\n";
    } elsif (! $srst2{WCIZ15B} || ($srst2{WCIZ15B} && @{$srst2{WCIZ15B}}[4])) {
	if ($srst2{WCIZ15B}) {
	    print $fh "@{$srst2{WZY15B}}[0]:@{$srst2{WCIZ15B}}[0]\tND:@{$srst2{WCIZ15B}}[4]\t@{$srst2{WZY15B}}[3]:@{$srst2{WCIZ15B}}[3]\t15C\n";
	} else {
	    print $fh "@{$srst2{WZY15B}}[0]:-\tND:-\t@{$srst2{WZY15B}}[3]:-\t15C\n";
        }
    }
}

if (exists $srst2{WZY18C} && ! @{$srst2{WZY18C}}[4]) {
    if ($srst2{WCIX18C} && ! @{$srst2{WCIX18C}}[4]) {
	print $fh "@{$srst2{WZY18C}}[0]:@{$srst2{WCIX18C}}[0]\t-:-\t@{$srst2{WZY18C}}[3]:@{$srst2{WCIX18C}}[3]\t18C\n";
    } elsif ($srst2{WCIX18C} && @{$srst2{WCIX18C}}[4]) {
	print $fh "@{$srst2{WZY18C}}[0]:@{$srst2{WCIX18C}}[0]\t-:@{$srst2{WCIX18C}}[4]\t@{$srst2{WZY18C}}[3]:@{$srst2{WCIX18C}}[3]\t18B:18C\n";
        my $sero_seq = freebayes_prior_fix($sero_bam, $sero_DB,"58__WCIX18C__WCIX18C-1__69");
        print $se "$sero_seq\n";
    }
}

if (exists $srst2{WZY22F} && exists $srst2{WCWA22A}) {
    print $fh "@{$srst2{WZY22F}}[0]:@{$srst2{WCWA22A}}[0]\tND:ND\t@{$srst2{WZY22F}}[3]:@{$srst2{WCWA22A}}[3]\t22A\n";
} elsif (exists $srst2{WZY22F} && exists $srst2{WCWA22F}) {
    print $fh "@{$srst2{WZY22F}}[0]:@{$srst2{WCWA22F}}[0]\tND:ND\t@{$srst2{WZY22F}}[3]:@{$srst2{WCWA22F}}[3]\t22F\n";
}

if (exists $srst2{WZY33F}) {
    if (exists $srst2{WCJE33A} && ! @{$srst2{WCJE33A}}[4]) {
	print $fh "@{$srst2{WZY33F}}[0]:@{$srst2{WCJE33A}}[0]\tND:ND\t@{$srst2{WZY33F}}[3]:@{$srst2{WCJE33A}}[3]\t33A\n";
    } else {
	print $fh "@{$srst2{WZY33F}}[0]\tND\t@{$srst2{WZY33F}}[3]\t33F\n";
    }
}

close $fh;
close $se;



###JUNK YARD###
=pop
if (exists $srst2{WZY22F} && ! @{$srst2{WZY22F}}[4]) {
    print "We found a perfect match: @{$srst2{WZY22F}}[0]\n";
}
elsif (exists $srst2{WZY22F} && @{$srst2{WZY22F}}[4]) {
    print "We have a imperfect match: @{$srst2{WZY22F}}[0]\n";
}

if (exists $srst2{WCJE33A} && ! @{$srst2{WCJE33A}}[4]) {
    print "We found a perfect match: @{$srst2{WCJE33A}}[0]\n";
}
elsif (exists $srst2{WCJE33A} && @{$srst2{WCJE33A}}[4]) {
    print "We have a imperfect match: @{$srst2{WCJE33A}}[0]\n";
}

foreach my $row (@srst2ARR) {
    my $target = @$row[0];
    my $diff = @$row[4];
    my $depth = @$row[3];
    if (exists $singlTarg95{$target}) {
        print "Found Single Target 95% Match!! || $target\t$diff\t$depth\n";
        if (! $diff) {
            print "perfect match\n";
            $diff = "NF";
        }
        print $fh "$target\t$diff\t$depth\t$singlTarg95{$target}\n";
    } elsif (exists $singlTarg100{$target} && ! $diff) {
        print "Found Single Target 100% Match || $target\t$diff\t$depth\n";
        print $fh "$target\tNF\t$depth\t$singlTarg100{$target}\n";
    }
}
=cut
