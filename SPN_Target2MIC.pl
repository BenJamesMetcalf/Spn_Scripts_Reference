#!/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;
use File::Copy qw(copy);
use Env;
#use lib $ENV{MODULESHOME}."/init";
use lib "/usr/share/Modules/init/";
use perl;

###Start Doing Stuff###
my $Res_output = "RES-MIC_".$ARGV[1];
open(my $fh,'>',$Res_output) or die "Could not open file '$Res_output' $!";
#print Dumper \@ARGV;
print "Output file name is: $Res_output\n";

my %Res_hash;
my $RES_full_name = $ARGV[0];
open(MYINPUTFILE, "$RES_full_name");
while(<MYINPUTFILE>) {
    #next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @res_arr;
    @res_arr = split('\t',$line);
    $Res_hash{$res_arr[0]} = $res_arr[1];
}
close MYINPUTFILE;

while (my ($key, $val) = each %Res_hash) {
    my @val_arr = split(':',$val);
    my @val_sort = sort(@val_arr);
    my $val_out = join(':',@val_sort);
    print "$key\t$val_out\n";
}
print "\n";
print Dumper \%Res_hash;

my %drug;
my %Out_hash;
###ER_CL Category###
$drug{ERY} = "=,0.06,S";
$drug{CLI} = "=,0.06,S";
$drug{SYN} = "<=,1.0,S";
$drug{LZO} = "<=,2.0,S";
$drug{ERY_CLI} = "neg";
my @Res_targs = split(':',$Res_hash{EC});
my $left_side = "no";
if ($Res_hash{"EC"} eq "neg") {
    #print "ER_CL,$Res_hash{EC},".$drug{ERY}.",".$drug{CLI}.",".$drug{SYN}.",".$drug{LZO}.",".$drug{ERY_CLI}."\n";
    $Out_hash{EC} = "$Res_hash{EC},$drug{ERY},$drug{CLI},$drug{SYN},$drug{LZO},$drug{ERY_CLI}";
} elsif (!grep(/R23S1.*(A(21|29)G|C(5|9)T)/i,@Res_targs) && !grep(/RPLD1.*(5K|G9R|Q5K|K6E)/i,@Res_targs)) {
    if (grep (/(R23S1|RPLD1|RPLV)/i,@Res_targs)) { 
	print "I'm in the right side of the EC decision tree\n";
	#flag everything
	#Run the ERMB+ERMBS code to overite flags
	$drug{ERY} = "Flag,Flag,Flag";
	$drug{CLI} = "Flag,Flag,Flag";
	$drug{SYN} = "Flag,Flag,Flag";
	$drug{LZO} = "Flag,Flag,Flag";
	$drug{ERY_CLI} = "Flag";
	if (grep (/ERM/i,@Res_targs) && !grep(/ERMBS/i,@Res_targs)) {
	    $drug{ERY} = ">,32,R";
	    $drug{CLI} = ">,2,R";
	    $drug{ERY_CLI} = "pos";
	    if (grep (/LSA/i,@Res_targs)) {
		$drug{SYN} = "Flag,Flag,Flag";
	    }
	}
	#print "ER_CL,$Res_hash{EC},".$drug{ERY}.",".$drug{CLI}.",".$drug{SYN}.",".$drug{LZO}.",".$drug{ERY_CLI}."\n";
	$Out_hash{EC} = "$Res_hash{EC},$drug{ERY},$drug{CLI},$drug{SYN},$drug{LZO},$drug{ERY_CLI}";
    } else {
	$left_side = "yes";
    }
} else {
    $left_side = "yes";
}

if ($left_side eq "yes") {
    print "I'm in the left side of the EC decision tree\n";
    #Go to left side of decision tree
    #Check for 3 known R23S Mutations (not R23SA21G)
    if (grep (/R23S1.*(A29G|C(5|9)T)/i,@Res_targs)) {
	print "Found R23S1\n";
	$drug{ERY} = "<=,0.25,S";
        $drug{CLI} = "<=,0.25,S";
    } 
    #Check for 4 known RPLD mutations
    if (grep (/RPLD1.*(5K|G9R|Q5K|K6E)/i,@Res_targs)) {
	print "Found RPLD1\n";
	$drug{ERY} = "<=,0.25,S";
        $drug{CLI} = "<=,0.25,S";
    }
    #Check for MEF
    if (grep(/MEF/i,@Res_targs)) {
        print "Found MEF\n";
        $drug{ERY} = "=,8,R";
    }
    #Check for LSA/LNU
    if (grep(/LSA/i,@Res_targs)) {
        print "Found LSA\n";
        $drug{CLI} = ">=,1,R";
    }
    if (grep(/LNU/i,@Res_targs)) {
        print "Found LNU\n";
        $drug{CLI} = "Flag,Flag,Flag";
    }
    #Check for R32SA21G and ERMB+ERMBS targets
    if (grep (/R23S1.*A21G/i,@Res_targs)) {
	print "Found R23S1-A21G\n";
	$drug{ERY} = ">,32.0,R";
        $drug{CLI} = "=,1.0,R";
        $drug{ERY_CLI} = "pos";
    }
    if (grep(/ERM/i,@Res_targs) && !grep(/ERMBS/i,@Res_targs)) {
	print "Found ERM and no ERMBS\n";
	$drug{ERY} = ">,32,R";
	$drug{CLI} = ">,2,R";
	$drug{ERY_CLI} = "pos";
	if (grep (/LSA/i,@Res_targs)) {
	    $drug{SYN} = "Flag,Flag,Flag";
	}
    }
    #print "ER_CL,$Res_hash{EC},".$drug{ERY}.",".$drug{CLI}.",".$drug{SYN}.",".$drug{LZO}.",".$drug{ERY_CLI}."\n";
    $Out_hash{EC} = "$Res_hash{EC},$drug{ERY},$drug{CLI},$drug{SYN},$drug{LZO},$drug{ERY_CLI}";
}

###COT Category###
$drug{SXT} = "=,0.5,S";
@Res_targs = split(':',$Res_hash{COT});
my $newTarg = "no";
if ($Res_hash{"COT"} eq "neg") {
    $Out_hash{"COT"} = "$Res_hash{COT},$drug{SXT}";
} else {
    #if (grep(/^FOLA-I20L$/i,@Res_targs)) {
    if (grep(/FOLA.*I20L/i,@Res_targs)) {
	if (grep(/FOLP_.*-ins/i,@Res_targs)) {
	    print "Found I20L AND FOLP insert\n";
            $drug{SXT} = ">=,4,R";
	} else {
	    print "Found I20L\n";
	    $drug{SXT} = "=,2,I";
	}
    #} elsif (grep(/FOLP_.*-ins/i,@Res_targs) && (!grep(/FOLA/i,@Res_targs) || grep(/^FOLA-(I20L|D12N)$/i,@Res_targs))) {
    } elsif (grep(/FOLP_.*-ins/i,@Res_targs) && (!grep(/FOLA/i,@Res_targs) || grep(/^FOLA-D12N$/i,@Res_targs))) {
	print "Found FOLP insert\n";
	$drug{SXT} = "=,2,I";
    } elsif (!grep(/^FOLA-D12N$/i,@Res_targs)) {
	print "Found new COT target\n";
        $drug{SXT} = "Flag,Flag,Flag";
    }
    $Out_hash{"COT"} = "$Res_hash{COT},$drug{SXT}";
}

###TET Category###
$drug{TET} = "<=,0.25,S";
$drug{DOX} = "<=,0.25,S";
if ($Res_hash{"TET"} eq "neg") {
    $Out_hash{"TET"} = "$Res_hash{TET},$drug{TET},$drug{DOX}";
} else {
    my @Res_targs = split(':',$Res_hash{TET});
    if ( grep( /TET/i, @Res_targs ) ) {
        $drug{TET} = ">,8,R";
        $drug{DOX} = ">=,1,R";
        $Out_hash{"TET"} = "$Res_hash{TET},$drug{TET},$drug{DOX}";
    }
}

###FQ Category###
$drug{CIP} = "NA,NA,NA";
$drug{LFX} = "<=,2,S";
@Res_targs = split(':',$Res_hash{FQ});
if ($Res_hash{"FQ"} eq "neg") {
    $Out_hash{"FQ"} = "$Res_hash{FQ},$drug{CIP},$drug{LFX}";
} elsif (!grep(/GYRA/i,@Res_targs) && grep(/PARC-S2[F|Y]/i,@Res_targs)) {
    $Out_hash{"FQ"} = "$Res_hash{FQ},$drug{CIP},$drug{LFX}";
} else {
    $drug{LFX} = "Flag,Flag,Flag";
    $Out_hash{"FQ"} = "$Res_hash{FQ},$drug{CIP},$drug{LFX}";
}

###OTHER Category###
$drug{CHL} = "<=,2,S";
$drug{RIF} = "<=,1,S";
$drug{VAN} = "=,0.5,S";
$drug{DAP} = "NA,NA,NA";
if ($Res_hash{"OTHER"} eq "neg") {
    print "OTHER,$Res_hash{OTHER},".$drug{CHL}.",".$drug{RIF}.",".$drug{VAN}.",".$drug{DAP}."\n";
    $Out_hash{"OTHER"} = "$Res_hash{OTHER},$drug{CHL},$drug{RIF},$drug{VAN},$drug{DAP}";
} else {
    my @Res_targs = split(':',$Res_hash{OTHER});
    if ($Res_hash{OTHER} ne "ant(6)-Ia:Ant6-Ia_AGly:aph(3')-III:Aph3-III_AGly:Sat4A_Agly" && $Res_hash{OTHER} ne "aph(3')-III:Aph3-III_AGly:Sat4A_AGly" && $Res_hash{OTHER} ne "msr(D):MsrD_MLS") {
        foreach my $target (@Res_targs) {
	    if ($target =~ m/CAT|RPOB|MSR/i) {
		if ($target =~ m/CAT/i) {
		    print "Found CAT\n";
		    $drug{CHL} = ">=,8,R";
		} elsif ($target =~ m/RPOB/i) {
		    print "Found RPOB\n";
		    $drug{RIF} = "Flag,Flag,Flag";
		}
            } elsif ($target !~ m/CAT|RPOB|MSR/i) {
                print "Found an ARGANNOT/RESFINDER target. Flag everything\n";
                $drug{CHL} = "Flag,Flag,Flag";
                $drug{RIF} = "Flag,Flag,Flag";
                $drug{VAN} = "Flag,Flag,Flag";
                last;
            }
        }
    }
    print "OTHER,$Res_hash{OTHER},".$drug{CHL}.",".$drug{RIF}.",".$drug{VAN}.",".$drug{DAP}."\n";
    $Out_hash{"OTHER"} = "$Res_hash{OTHER},$drug{CHL},$drug{RIF},$drug{VAN},$drug{DAP}";
}

print $fh $Out_hash{EC}.",". $Out_hash{COT}.",".$Out_hash{TET}.",".$Out_hash{FQ}.",".$Out_hash{OTHER}."\n";
print "$Out_hash{EC}||$Out_hash{COT}||$Out_hash{TET}||$Out_hash{FQ}||$Out_hash{OTHER}\n";
