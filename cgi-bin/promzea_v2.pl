#!/usr/bin/perl

##########################
######## settings ########
##########################

##
# verbosity determons the loging level of this script
#	0=nothing
#	1=errors only
#	2=timeings
#	3=progress output (calls to 3rd party tools)
#	4=debugging information
##
my $VERBOSITY = 1;

##
# The url prefix for linking to results in the email 
##
#my $URL_FOR_EMAIL = 'http://192.168.1.110/promzea/html/results/';
#my $URL_FOR_EMAIL = 'http://131.104.92.180/results/';
my $URL_FOR_EMAIL = 'http://promzea.cs.uoguelph.ca/results/';

##
# The url blast version
##
#my $blastdr = "bin/blast-2.2.25/bin/blastall";
my $blastdr = 'bin/blast-2.2.21/bin/blastall';

##########################
######## imports #########
##########################

use strict;
use warnings; 
use PromzeaModule;
use Cwd;
use Cwd qw(abs_path);
use Data::Dumper;
use Email::Valid;
use Fcntl qw(:flock :seek);
use File::Path;
use IO::CaptureOutput qw(capture_exec capture);
use List::Util qw(max);
use MIME::Lite;
use Net::SMTP;
use Statistics::Basic qw(mean);
use Time::HiRes qw(gettimeofday);

##########################
########## main ##########
##########################

log_time("START");

my $email;
my $in;
my $promlength;
my $seq;
my $organism;
my $stringence; 
my $cdnafile; 
my $matrix; 
my $seq_microarray; 
my $email_verif;

#	get the request
our $workfile=$ARGV[0];
open (FILE, "<$workfile");
while (my $line = <FILE>) {
	my $aPrefix=$line;
	$aPrefix =~ s/^([^\=]+)\=.*$/$1/g;
	$aPrefix =~ s/\n+//g;
	$line =~ s/^[^\=]+\=(.*)$/$1/g;
	$line =~ s/\n+//g;
	$line =~ s/\\n/\n/g;
	if($aPrefix eq "email" ){
		$email = $line;
	}elsif($aPrefix eq "editForm" ){
		$in = $line;
	}elsif($aPrefix eq "promlength" ){
		$promlength = $line;
	}elsif($aPrefix eq "seq" ){
		$seq = $line;
		# $seq =~ s/[\n\r]+/~/g;
		# $seq =~ s/(>[^>~]+)~/\n$1\n/g;
		# $seq =~ s/~//g;
	}elsif($aPrefix eq "cdnafile" ){
		$cdnafile = $line;
	}elsif($aPrefix eq "organism" ){
		$organism = $line;
	}elsif($aPrefix eq "stringence" ){
		$stringence = $line;
	}elsif($aPrefix eq "matrix" ){
		$matrix = $line;
	}elsif($aPrefix eq "seq_microarray" ){
		$seq_microarray = $line;
	}elsif($aPrefix eq "email_verification" ){
		$email_verif = $line;
	}else{
		dieemail("Unknown aPrefix=($aPrefix)");
	}
}
close(FILE);


# Count file for display variable incrementation
open(IN,"+<","../html/count.txt") or dieemail("Can't open counts for reading: $!"); 
flock(IN,LOCK_EX);
seek(IN,0,SEEK_SET); 
my $count_execution = <IN>; 
$count_execution = $count_execution + 1;
truncate(IN,0);
seek(IN,0,SEEK_SET);
print IN "$count_execution\n";
close(IN);

my $display = sprintf("%08d", $count_execution);
# 	Variable initialization									  
my $seqfile_cDNA;
# 	CGI variable initialization	
if(!defined($in)){
	$in="Untitled";
}
my $seq_config =  "";
my $file_config = "a-zA-Z0-9_.-";
my $ariz_config = "";

# 	Conversion of user data to promoter ID and file			
my ($option,@seqfile_cDNA);
if(!defined($promlength)){
	dieemail("Missing or Input; promlength.");
}
if( defined($seq) && $seq=~ m/(\>[0-9A-Z\.\_]+[\n\b\r\s]+[ATGC\n\r\s\b]+)+/i ){
	@seqfile_cDNA = $seq;
	open (OUT, ">>","../html/temp/$display.fas");
	foreach (@seqfile_cDNA){ print OUT  $_;}
	close OUT;
	$seqfile_cDNA = "../html/temp/$display.fas";
	$option = 1;
}elsif( defined($seq) && length($seq)>3 ){
	dieemail("Incorrect Input. seq needs a \">tag\", got(".$seq.")");
}elsif( defined($cdnafile) ){
	$seqfile_cDNA = "../html/temp/$display.fas";
	open (OUT,"+>","$seqfile_cDNA") or dieemail("Can't open outfile for writing on : $!");
	flock( OUT, LOCK_EX );
	print OUT $cdnafile;
	close(OUT);
	$option = 1;
}elsif( defined($matrix) && defined($seq_microarray) && ($matrix =~ 'Affy46k') and ($seq_microarray=~ m/[(Zm|ZmAffx).(\.).(\d+).(\.).(\d).(\.).(\w+)]+/i) ){
	@seqfile_cDNA = $seq_microarray; 
	my $malength=scalar(@seqfile_cDNA);
	if($malength<2){
		dieemail("seq_microarray must be of length 2 to 100 but we have $malength");
	}
	open (OUT, "+>>../html/temp/$display.txt");
	foreach(@seqfile_cDNA){print OUT  $_;}
	close OUT;
	$seqfile_cDNA = "../html/temp/$display.txt";
	$option = 2;
}elsif( defined($matrix) && defined($seq_microarray) && ($matrix =~ 'Gramene') and ($seq_microarray=~ m/[(GRMZM2G|GRMZM5G).(\d+).(\w+)]+/i) ){
	@seqfile_cDNA = $seq_microarray ; 
		open (OUT, "+>>../html/temp/$display.txt");
		foreach (@seqfile_cDNA){
			print OUT  $_;
		}
		close OUT;
		$seqfile_cDNA = "../html/temp/$display.txt";
		$option = 3;
}else {
	dieemail("Missing or incorrect Input. [seq(".$seq.") | cdnafile(".$cdnafile.") | matrix(".$matrix.") & seq_microarray(".$seq_microarray.") ");
};

#	Perl CGI core script 
#------------------------------------------------------------------------
#	Perl CGI core script 
my $home = 	$ENV{'PWD'};
my @home = split("/",$home);
pop @home;
my $sitehome = "/var/www";
my $site_results = "$sitehome/html/results/"."$display";
my $siteimage = "$sitehome/html/results/"."$display";
my $image_dir = ".";
my($motifs, $metric_HoA,$consensus_seq,$go_desc,@motif,$stamp);
my(%metric_HoA,%metric,%consensus_seq,%go_desc,%stamp,$i);

#log_time("maize_motif start");
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


#	send the e-mail results, transform results in zip and pdf files										 
if (@motif){
	#system("htmldoc --webpage -f $site_results/$display.pdf $site_results/$display.html");
	my $URL = "$URL_FOR_EMAIL$display/$display.html";
	my @answers = run("zip -r $site_results $site_results",".");
	#convert to pdf and attach to email
	@answers = run("htmldoc --webpage --size a4 --pagemode outline -f $site_results/$display.pdf $site_results/$display.html",".");
	@answers = run("/usr/local/bin/sendEmail -a $site_results/$display.pdf -a $site_results.zip -f \"test\@promzea.cs.uoguelph.ca\" -t \"$email\" -bcc \"timothy.d.f.l\@gmail.com\" -s mail.uoguelph.ca:587 -xu cliseron -xp cliseron -u \"Promzea results $display\" -m \"Results from Promzea query are visible at $URL\" -o tls=auto -q",".");
}else{ dieemail("Promzea is not able to determine any motif.")};


unlink ("../html/temp/$display.fas","../html/temp/$display.txt");
unlink("$workfile");
# maybe rename/move workfile to the results?
#rmtree ("$site_results");
#unlink ("$site_results.zip");
exit 0;

# 	Subroutine

sub dieemail {
	my ($msg) = @_;
	$msg  =~ s/[\r\n]+/\\n/g;
	$msg  =~ s/\"/&quot;/g;
	my @answers = run("/usr/local/bin/sendEmail -a $workfile -f \"test\@promzea.cs.uoguelph.ca\" -t \"timothy.d.f.l\@gmail.com\" -s mail.uoguelph.ca:587 -xu cliseron -xp cliseron -u \"Promzea error\" -m \"$msg\" -o tls=auto -q",".");
	unlink("$workfile");
	#rmtree ("$site_results");
	#unlink ("$site_results.zip");
	#print "$msg\n";
	exit -1;
}

