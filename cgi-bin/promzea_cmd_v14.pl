#!/usr/bin/perl -w
# cd /srv/www/public_html/cgi-bin
# perl promzea_cmd_v13.pl -i (anthocyanin_pathway.fas,genelist.txt or affymetrixprobeslist.txt) -p (200, 500 or 1000) -t (fasta, gene, affy or bed) -s (1 or 2) -o (ZM or AT)
my @timeData = localtime(time);
print ("$timeData[2]",':',"$timeData[1]",':',"$timeData[0]\n"); 
	
##
# The url prefix for linking to results in the email 
##
my $URL_FOR_EMAIL = 'http://131.104.92.180/results/';

##########################
######## imports #########
##########################
use strict;
use warnings;
use lib "/var/www/cgi-bin";  
use PromzeaModule;
# 	perl modules				  
use CGI qw(:standard );
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Fcntl qw(:flock :seek);
use Data::Dumper;
use Email::Valid;
use List::Util qw(max);
use MIME::Lite;
use Net::SMTP;						  
use Getopt::Long;


##########################
########  Main 	 #########
##########################
log_time("START");

#__________________________________________________________________________________________

#	Count file for display variable incrementation
open(IN,"+<","../html/count.txt")or die("Can't open counts for reading: $!");
flock(IN,LOCK_EX);	# lock the file (exclusive lock)
seek(IN,0,SEEK_SET);	# rewind it to the beginning
my $count_execution = <IN>;       # read only the first line.
$count_execution = $count_execution + 1;    # increment the counter
truncate(IN,0);         # this erases (truncates) the file to length=0
seek(IN,0,SEEK_SET);	# rewind it to the beginning
print IN "$count_execution\n";   # write out the new count
close(IN);             # close the file.
my $display = sprintf("%08d", $count_execution);
#__________________________________________________________________________________________

#	Help initialization
printhelp() unless @ARGV;

my ($in,$promlength,$type,$organism,$stringence);
GetOptions(
	"p=s"=>\$promlength,
	"t=s"=>\$type,
	"i=s"=>\$in,
	"o=s"=>\$organism, # AT:A.thaliana, ZM: Zea mays B73
	"s=s"=>\$stringence
);
foreach (@ARGV) {
    help() if (/^-h\s/i);
}
help("1") unless ($in);
help("2") unless ($promlength =~ /200|500|1000|intron|utr5|utr3|bed/);
help("3") unless ($type =~ /fasta/ or $type =~ /gene/ or $type =~ /affy/ or $type=~/bed/ or $type=~/seq/);
my ($seqfile_cDNA,@file,$option);
#__________________________________________________________________________________________

# 	Conversion of user data to promoter ID and file			
open(READIN,"$in");
my @seqfile_cDNA = <READIN>;
close READIN;

if($type =~ /fasta/i){
	$seqfile_cDNA = "../html/temp/$display.fas";
	open (OUT, "+>>", $seqfile_cDNA);
	foreach(@seqfile_cDNA){
		print OUT  $_;
		}
	close OUT;
	$option = 1;
}
elsif($type =~ /affy/i){
		$seqfile_cDNA = "..html/temp/$display.txt";
		open (OUT, ">>",$seqfile_cDNA);
		foreach(@seqfile_cDNA){
			chomp;
			print OUT "$_\n";
		}
		close OUT;
		$option = 2;
}
elsif($type =~ /gene/i){
		$seqfile_cDNA = "../html/temp/$display.txt";
		open (OUT, ">>",$seqfile_cDNA);
		foreach(@seqfile_cDNA){
			chomp;
			print OUT "$_\n";
		}
		close OUT;
		$option = 3;
}
elsif($type =~ /bed/i){
		$seqfile_cDNA = "../html/temp/$display.bed";
		open (OUT, ">>",$seqfile_cDNA);
		foreach(@seqfile_cDNA){
			chomp;
			print OUT "$_\n";
		}
		close OUT;
		$option = 4;
}
elsif($type =~ /seq/i){
    $seqfile_cDNA = "../html/temp/$display.fasta";
    open (OUT, ">>",$seqfile_cDNA);
    foreach(@seqfile_cDNA){
        chomp;
        print OUT "$_\n";
    }
    close OUT;
    $option = 5;
}
else{die("Missing or incorrect Input")};
#__________________________________________________________________________________________
#	Perl CGI core script 
#my $home = $ENV{PWD};
#my @home = split("/",$home);
#pop @home;
my $sitehome = "/var/www";
my $site_results = "$sitehome/html/results/"."$display";
my $siteimage = "$sitehome/html/results/"."$display";
my $image_dir = ".";
my($motifs, $metric_HoA,$consensus_seq,$go_desc,@motif,$stamp);
my(%metric_HoA,%metric,%consensus_seq,%go_desc,%stamp,$i);

##
# test01 | test02
# up to here is quick 0.018439 | 0.018173  seconds
# this next line takes 83.376374 | 74.219779  seconds
# after this next line is quick 0.002417 | 0.249264 seconds
##
log_time("maize_motif start");

#__________________________________________________________________________________________
#	Apply motif discovery tools, filtering and graphic generation: 
($motifs, $metric_HoA,$consensus_seq,$go_desc,$stamp) = maize_motif($organism,$promlength,$stringence, $seqfile_cDNA,$option,$sitehome,$siteimage);
@motif = @$motifs;
%metric_HoA = %$metric_HoA;
%consensus_seq = %$consensus_seq;
%go_desc = %$go_desc;
%stamp = %$stamp;
log_time("maize_motif end");

#	Rank motifs using metric results: 
my (@bgd,@fgd,$metric,$motif_name,$nmotif);
foreach $motif_name(keys (%metric_HoA)){
	@bgd=();
	@fgd=();
	for (my $i=0; $i<=$#{$metric_HoA{$motif_name}[0]};$i++){
		push @fgd, $metric_HoA{$motif_name}[0][$i];
	}
	for (my $i=0; $i<=$#{$metric_HoA{$motif_name}[1]};$i++){
		push @bgd, $metric_HoA{$motif_name}[1][$i];
	}
	$metric = MNCP(\@fgd,\@bgd);
	$metric{$motif_name}=$metric;
}

#	Define the height of html table in function of motif number
$nmotif=0;
foreach $motif_name(keys(%metric)){$nmotif++;}
my $n_motif_height = 100*$nmotif;

#	HTML ouput and plug-in files
promzea_html_output(\@motif,\%metric,\%stamp,\%consensus_seq,\%go_desc,$site_results,$image_dir,$display,$n_motif_height,$in);


if (@motif){
	my @answers = run("zip -r $site_results $site_results",".");
	#convert to pdf and attach to email
	@answers = run("htmldoc --webpage --size a4 --pagemode outline -f $site_results/$display.pdf $site_results/$display.html",".");
}else{ die("Promzea is not able to determine any motif.")};
#__________________________________________________________________________________________________________
sub help{

    print STDERR "Input file is required\n" if (defined $_[0] && $_[0] eq "1");
	print STDERR "Promoter length is inaccurate choose a promoter size 200,500 or 1000\n" if (defined $_[0] && $_[0] eq "2");
	print STDERR "type file is required \n\n" if (defined $_[0] && $_[0] eq "3");


    print STDERR <<HELP;
$0:
-i The input file containig cDNA fsya sequences or a list affymetrix or gramene ID.
-p The promoter length to be tested by Promzea 200,500 or 1000
-t The type of file given fasta, affy or gramene
-h prints the help describing parameter to enter in promzea command line.

HELP

    exit;
}



log_time("END");

exit 0;


