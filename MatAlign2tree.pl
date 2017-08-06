#!/usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Carp;
use File::Spec::Functions;
#use File::Temp qw/ tempfile tempdir /;
use File::Copy;
#use Data::Dumper::Simple;
use Cwd qw(abs_path);
use File::Basename;
#use IO::All;

my $curr_path 		= dirname(abs_path($0));
my $pcm_path		= 'PCMs';
my $inDir			= '.';
my $extension		= 'pcm';
my $format			= 'pcm';
my $outDir			= 'output';
my $fixnames        = 1;
my $adjustNames     = 1;
my $verbose			= 0;
my $man				= 0;
my $help			= 0;
my $neighbor_path	= "$curr_path/app/neighbor.app/Contents/MacOS/neighbor";
my $matalign_path	= "$curr_path/app/matalign-v4a";
my $logolist;
my $treeformat      = "UPGMA";
=head1 NAME

    MatAlign2tree.pl

=head1 SYNOPSIS

    perl MatAlign2tree.pl [options]

    Required:
	
    Optional:
    --in          <dir>    Input folder containing pcm and pfm files path, default '.'
    --out         <dir>    Output dir name, Default (output)
    --pcmpath     <path>   relative pcm path to input folder.
                           Default('PCMs')
    --matalign    <path>   Path to matalign-v4a program. Default('this_script_path/app/matalign-v4a')
    --neighbor    <path>   Path to neighbor program. Default('this_script_path/app/neighbor.app/Contents/MacOS/neighbor')
    --ext         <str>    Extension of matrix files, Default (pcm)
    --format      <str>    format of matrix, pcm or pfm, Default (pcm)
    --fixnames             Assign each name a number and pass this number to the Neighbor
                           program.  Then replace numbers with names later.
                           Default (on) to turn off, --nofixnames
    --logolist    <file>   a logo list file.
    --tree        <str>    default:UPGMA. could be UPGMA and NJ
                           
    --help                 Print help message
    --man                  Print manual page  

=head1 DESCRIPTION

    Given a directory of matrices, run MatAlign, produce a distance matrix
    and run the Neighbor tree estimation program on this distance matrix to 
    generate a NJ tree. Write all of the generate files to subdir.
    
    Requires that phylip and matalign-v4a be installed!

=cut

my $argv = GetOptions	# get command line options
(
	"in=s"      	=>      \$inDir,
	"out=s"   		=>      \$outDir,
	"pcmpath=s"		=>		\$pcm_path,
	"matalign=s"	=>		\$matalign_path,
	"neighbor=s"	=>		\$neighbor_path,
	"fixnames!"		=>		\$fixnames,
	"ext=s"     	=>      \$extension,
	"format=s"		=>		\$format,
	"logolist=s"	=>		\$logolist,
    "tree=s"        =>      \$treeformat,
	"verbose"       =>      \$verbose,
	help            =>	sub { pod2usage(1) },	
	man				=>	sub { pod2usage(-verbose => 2) },
) or pod2usage(2);
pod2usage("--in is a required parameter\n") unless(defined $inDir);

die "A file named 'outfile' already exists, which will cause errors\n" unless(! -e 'outfile');
die "A file named 'outtree' already exists, which will cause errors\n" unless(! -e 'outtree');
die "format must be 'pcm' or 'pfm'\n" unless($format eq "pcm" || $format eq "pfm");

pod2usage("The path to matalign($matalign_path) is invalid.\n") unless(-e $matalign_path);

mkdir($outDir) unless(-d $outDir);

$pcm_path = catdir($inDir,$pcm_path);

my %mx;
my %namesH;

my @dirFiles;
opendir(DIR, $pcm_path) or die "Couldn't open directory \'$pcm_path\': $!";
@dirFiles = grep{-f "$pcm_path/$_" && /\.$extension$/} readdir(DIR);
close DIR;

if(defined $logolist){
	my @lists;
	open LIST, "<", "$logolist" or die "Couldn't open \"$logolist\": $!";
	while (my $list = <LIST>)
	{
		chomp $list;
		push @lists,$list.".$extension" if(grep $_ eq $list.".$extension", @dirFiles);
	}
	close LIST;
	@dirFiles = @lists;
}

die "The directory \"$pcm_path\" contains no matrices with file extension \"$extension\"\nPlease check the file extension and the directory path\n" unless(@dirFiles);

#calculate pfm to generate alignments
foreach my $fm (@dirFiles){
	open PFM,"<",catfile($pcm_path,$fm) or die "Couldn't open $fm pcm file: $!";
	my @motif = ();
	while(my $pfm = <PFM>){
		chomp $pfm;
		if($pfm=~/^[ACGT]\s+\|\s+(\d.+\d)/){
			my @v = split(/\s+/, $1);
			push @motif,[@v];
		}
	}
	close PFM;
	if($adjustNames==1){
		$fm =~ s/[\(\):]+/-/g;
	}
	print $fm if($verbose);
    check($fm, @motif);
	if($format eq "pfm"){#convert pfm to pcm for matalign
		my @motifPCM = PFM2PCM(map{[@$_]} @motif);
		my $outPCMDir = catdir($outDir, "pcm");
		mkdir($outPCMDir) unless(-d $outPCMDir);
		open OUTPFM,">",catfile($outPCMDir,$fm) or die "Couldn't write to $fm pcm file: $!";
		my @a=qw(A C G T);
		for(my $i=0; $i<scalar(@motifPCM); $i++){
			print OUTPFM $a[$i]," | ",join(" ",@{$motifPCM[$i]}),"\n";
		}
		close OUTPFM;
	}
	print Dumper(@motif) if($verbose);
}

# write f1 file containing names of all matrices to compare for megaMatalign
my $f1_file = catfile($outDir, "f1.matrix_paths.txt");
open(F1, ">$f1_file") or die "Couldn't open \"$f1_file\": $!";
if($format eq "pcm"){ print F1 File::Spec->rel2abs( $pcm_path ), "\n";}
else{ print F1 File::Spec->rel2abs( catdir($outDir, "pcm") ), "\n";}
foreach my $file(@dirFiles)
{
        print F1 "$file\n"
} 
close F1;
print "checkpoint matalign0" if($verbose);
# run matalign program
my $matalign_out_file = catfile($outDir, 'matalign.out');
system("$matalign_path -A A:T C:G -f1 $f1_file > $matalign_out_file");
print "checkpoint matalign1" if($verbose);

open(IN, $matalign_out_file) or die "Couldn't open \"$matalign_out_file\": $!";
my $startState = 0;
while(<IN>)
{
        if($startState)
        {
                #Matrix_1, Consensus, Matrix_2, Consensus, ALLR, Overlap, Distance, Aligned_Dist, Shared_Consensus, E_value,            #P_value
                my @r = split(/\s+/);
                my $m1 = $r[0];
                $m1 =~ s/\.$extension$//g;

                my $m2 = $r[2];
                $m2 =~ s/\.$extension$//g;
				
				if($adjustNames eq 1)
				{
					$m1 =~ s/[\(\):]+/-/g;
					$m2 =~ s/[\(\):]+/-/g;
				}
				
				my $tmp = $r[6];
				$tmp = 0.00 if($tmp < 0);
                $mx{$m1}{$m2}  = $tmp;
                $mx{$m2}{$m1}  = $tmp;
        }
        elsif($startState || /^Matrix_1\s+Consensus\s+Matrix_2\s+Consensus\s+ALLR\s+Overlap\s+Distance\s+Aligned_Dist\s+Shared_Consensus\s+E_value\s+P_value/)
        {
                $startState = 1;
        }
}
close IN;

#sort all names
my @names = sort keys %mx;

my %num2name;
my %name2num;
if($fixnames eq 1)
{	
	my $cntr = 1;
	# output a table indicating the number assigned to every name
	my $num2namesPath = catfile($outDir, 'num_2_names.tsv');
	open my $N2N, ">", $num2namesPath or croak "Couldn't open \"$num2namesPath\": $!";
	foreach my $name(@names)
	{
		my $num = sprintf("N%06iN", $cntr);
		print $N2N "$name\t$num\n";
		$num2name{$num} = $name;
		$name2num{$name} = $num;
		$cntr++;
	}	
	close $N2N;
} 

print "checkpoint 1" if($verbose);
my $distMX_file = catfile($outDir, 'matalign.distMX.txt');
open(OUT, ">$distMX_file") or die "Couldn't open \"$distMX_file\": $!";
print OUT scalar @names, "\n";
foreach my $n1(@names)
{
        if($fixnames eq 1)
        {
        	printf(OUT "%-13s", $name2num{$n1});
        }
        else
        {
	        my $tmpName = $n1;
	        $tmpName =~ s/^(.+)\..+?W(\d+)\.motif_(\d+)/$1.$2.$3/g;    
	        printf(OUT "%-13s", substr($tmpName, 0, 10));
	#       printf OUT("%-13s", substr($n1, 0, 10));
        }

        foreach my $n2(@names)
        {
                if($n1 eq $n2)
                {
                        printf(OUT "% 8.4f ", 0.00);
                }
                elsif(defined $mx{$n1}{$n2})
                {
                        printf(OUT "% 8.4f ", $mx{$n1}{$n2})
                }
                else
                {
                        printf(OUT "% 8.4f ", $mx{$n2}{$n1})
                }
        }
        print OUT "\n";
}
close OUT;

# run neighbor NJ tree building program
#my ($FH, $inDir_path) = tempfile();
#print $FH "$distMX_file\nY\n";
print "checkpoint 2" if($verbose);
my $cmdfile_path = catfile($outDir, 'infile');
$cmdfile_path = File::Spec->rel2abs( $cmdfile_path ); # make the output path absolute
open my $INFILE, ">", $cmdfile_path or croak "Couldn't open \"$cmdfile_path\": $!";
if($treeformat eq "UPGMA") {
    print $INFILE "$distMX_file\nN\nY\n";
}else{
    print $INFILE "$distMX_file\nY\n";
}
close $INFILE;

system("$neighbor_path < $cmdfile_path > /dev/null");
print "checkpoint 3" if($verbose);
if($fixnames eq 1)
{
	# replace all of the numbers in outtree with the proper names
	my $outTreeFile = catfile($outDir,"NJ.matalign.distMX.nwk");
	open my $TREE, ">", $outTreeFile or croak "Couldn't open \"$outTreeFile\": $!";
	
	open my $IN, "<", "outtree" or die "Couldn't open \"outtree\": $!";
	while(my $linein=<$IN>)
	{
		foreach my $num(sort keys %num2name)
		{
			$linein =~ s/$num/$num2name{$num}/g
		}
		print $TREE $linein;
	}
	close $IN;
	close $TREE;
	
	move("outtree", catfile($outDir,"outtree"));
	move("outfile", catfile($outDir,"outfile"));
}
else
{
	move_message("outtree", catfile($outDir,"NJ.matalign.distMX.nwk"));
	move_message("outfile", catfile($outDir,"NJ.matalign.distMX.neighbor.out"));
}

unlink(catfile($outDir,"matalign.distMX.txt"));
unlink(catfile($outDir,"f1.matrix_paths.txt"));
unlink(catfile($outDir,"infile"));

################################   subroutines   ################################################
sub move_message
{
        my($from, $to) = @_;
        print "file name \"$from\" moved to \"$to\"\n";
        move($from, $to);
}

sub PFM2PCM{
	my (@M) = @_;#motif rows is ACGT
	my $total = 1000;
	my $len=scalar(@{$M[0]});
	for(my $i=0; $i<$len; $i++){
		for(my $j=0; $j<4; $j++){
			$M[$j][$i]=int($M[$j][$i]*$total+0.49);
		}
	}
	return @M;
}

sub check{
	my ($name, @motif) = @_;#motif rows is ACGT
	if(scalar(@motif)!=4) {
		print Dumper(@motif);
		die("motif $name format is not good: $!");
	}
}