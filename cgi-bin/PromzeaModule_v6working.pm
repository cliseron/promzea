package PromzeaModule_v6;
use strict;
use warnings;
use Exporter qw(import);
my @subs = ("run","promzea_html_output","v_log","log_time","maize_motif","fasta_to_fastaline","MNCP");
our $VERSION = 1.00;
our @ISA = qw(Exporter);
our @EXPORT=@subs;
our %EXPORT_TAGS = ( 	DEFAULT => \@subs);

use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Fcntl qw(:flock :seek);
use Data::Dumper;
use Email::Valid;
use List::Util qw(max);
use MIME::Lite;
use Net::SMTP;
use Statistics::Basic qw(mean);
use Cwd;
use Cwd qw(abs_path);
use IO::CaptureOutput qw(capture_exec);
use Time::HiRes qw(gettimeofday);
use IO::CaptureOutput qw(capture_exec capture);
use Bio::SeqIO;
use Bio::SeqIO::fasta;
use Bio::DB::Fasta;
use Bio::SearchIO;
use Bio::PrimarySeq;
use Bio::Perl;
use Bio::SeqFeature::Generic;
use Parallel::ForkManager;

##
# verbosity determons the loging level of this script
#	0=nothing
#	1=errors only
#	2=timeings
#	3=progress output (calls to 3rd party tools)
#	4=debugging information
##
my $VERBOSITY = 0;

# run a 3rd party command
# @return ($stdout, $stderr, $success, $exit_code)
sub run{
        my ($command,$run_from)= @_;
        v_log(3,$command);
        my ($s, $usec) = gettimeofday();
                my $dir = getcwd;
                if($run_from ne "."){
                                chdir $run_from;
                }
                my @answers = capture_exec($command);
                chdir $dir;
                $answers[3]=$answers[3]/256;
        if($answers[2]==0){
                v_log(1,"$command\n".$answers[1]);
        }
        my ($s2, $usec2) = gettimeofday();
        my $time = $s2-$s;
        v_log(3,"time=$time");
        v_log(4,"err=".$answers[1]);
        if ($answers[3]!=0 && $command !~ m/.*pscan.*/ ){
                if( $command =~ m/.*Email.*/ ){
                        die("3rd party app($command) failed:\nout=$answers[1]\nerr=$answers[2]\n");
                }
        }
        return @answers;
}

#______________________________________________________________________________________________________
##
# logs messages depending on the verbosity level
sub v_log{
	my ($level,$message)= @_;
	if($level<=$VERBOSITY){
		print STDERR "LOG_".$level.": ".$message."\n";
	}
}
#__________________________________________________________________________________________
sub log_time{
	my ($message) = @_;
	my ($s, $usec) = gettimeofday();
	v_log(2, " TIMER: $s.$usec\t$message");
}
#__________________________________________________________________________________________
sub promzea_html_output{
	my($motif,$metric,$stamp,$consensus_seq,$go_desc,$site_results,$image_dir,$display,$n_motif_height,$in) = @_;
	my($motif_name,$n,$pictogram,$value,$motif_progs,$seq_consensus),
	my($input_motif_freq,$genomefreq,$GOpiechart,$goview,$golist);
	my @motif= @$motif;
	my %metric = %$metric;
	my %stamp = %$stamp;
	my %consensus_seq = %$consensus_seq;
	my %go_desc = %$go_desc;
	my $stamp_file = "inputstamp.txt";
	if(@motif){
		# file for stamp plug in
		$n=1;
		open(OUT,"+>", "$site_results/inputstamp.txt");
		foreach $motif_name(sort {$metric{$b} <=> $metric{$a}}(keys(%metric))) {
			print OUT ">Motif$n\n".$stamp{$motif_name};
			$n++;
		}
		close OUT;	
		#	HTML result page
		#	Send the user to create HTML page created.
		open (HTML, "+>", $site_results."/".$display.".html");
		print HTML qq(<head> <title> results - $display </title>
		<meta http-equiv="Content-Type" content="text/html; charset=utf-8" /><style type="text/css"></style><link href="../../stylesheet2.css" rel="stylesheet" type="text/css" />
		</head><body class="background"><div class="maincontent" id="results">
			<div class="header"><img src="../../images/home_menu.png" alt="home direction" width="57" height="23" border="0" align="baseline" usemap="#Map" />
			<map name="Map" id="Map"><area shape="rect" coords="4,5,717,99" href="../../index.html" alt="home_direction" />
			<area shape="rect" coords="52,19,68,37" href="#" /></map> </div><div class="results" id="results_wrap"> );
		#----------------------------------------------------------------------
		print HTML qq(<div class="results" id="results_wrap"><p> <h2><center>Results Summary</center></h2></br ><h2><center>$in</center></h2></br ><center>Promzea - $display</center></p></div>
		<table width="500" height="$n_motif_height" align="center" cellpadding="0" cellspacing="0"></td></tr>);
		$n=1;
		foreach $motif_name(sort {$metric{$b} <=> $metric{$a}}(keys(%metric))) {
			$pictogram = "$image_dir/"."$motif_name".'logo.png';
			$value = $metric{$motif_name};
			$motif_progs = substr($motif_name, 0, -1);
			print HTML qq(<tr><td width="300" height="200" align="center" valign="middle">&nbsp;<img  src="$pictogram"alt= "No Pictogram"></img></td>
				<td width="100" align="center" valign="middle">&nbsp; <b> Motif$n </b><br />(from $motif_progs) <br /> <b>score: $value </b> </td> </tr>);
			$n++;
		}
		print HTML qq(</table>);
		print HTML qq(<p align= "center">Compare your motifs to known promoter motif databases using STAMP website</br >
		<a HREF=\"$stamp_file\"target="_blank"> motif file to copy in STAMP website</a></br ></br ></p>
		Open the above link, copy content of the newly open file and paste in STAMP program link below</br >
					In STAMP, under "Similarity Matching", we suggest selecting the plant motif databases:</br >
					Athamap, AGRIS, PLACE, TRANSFAC; then submit</br ></p>); 
		print HTML qq(<p align= "center"><a HREF=\"http://www.benoslab.pitt.edu/stamp/\" target="_blank"> STAMP website</a></br ></br ></p>); 
					
		$n =1;
		foreach $motif_name (sort {$metric{$b} <=> $metric{$a}}(keys(%metric))) {
			$pictogram = "$image_dir/"."$motif_name".'logo.png';
			$input_motif_freq = "$image_dir/sequence"."$motif_name".'.png';
			$genomefreq = "$image_dir/graph"."$motif_name".'.png';
			$GOpiechart = "$image_dir/GO"."$motif_name".'.png';
			$goview = "./go/"."$motif_name"."_gofile.txt";
			$golist = "./go/"."$motif_name"."_golist.txt";
			$seq_consensus = $consensus_seq{$motif_name};
			print HTML qq(<hr  />  <hr  /> <table width="720" height="940" border="0" align="center" cellpadding="0" cellspacing="0">);
			print HTML qq(<tr><td height="20" colspan="2" align="center" valign="middle"><h2>Motif$n</h2></td> </tr>);
			print HTML qq(
	<tr><td width="200" align="center" valign="middle">&nbsp;<img  src="$pictogram"alt= "No Pictogram" width="200" height="150"></img></td>
		<td rowspan="2" width="490" height="380" align="center" valign="middle">&nbsp; <img src="$input_motif_freq" alt= "No frequency graph"></img></td> </tr>);
			print HTML qq(<tr><td align="center"><b>$seq_consensus</b></td></tr> );
#			if($#{$go_desc{$motif_name}} != -1){
				print HTML qq(<tr><td height="20" colspan="2" align="center" valign="middle"> <hr /> Motif$n annotation in the genome </td> </tr>);
				print HTML qq(<tr><td colspan="2" align="center" height="400" valign="middle">&nbsp;<img src="$GOpiechart" alt="No specific annotation found"></img></td> </tr>);
				if (defined($go_desc{$motif_name})){
					print HTML qq(<tr><td colspan="2" align="center" height="60" valign="middle"><b>Annotation complete description</b> <div class="scroll"><p>);
					#print @{$go_desc{$motif_name}}."\n";
					if(@{$go_desc{$motif_name}} <= 6){
						for my $i(0 .. $#{$go_desc{$motif_name}}) {
							my $GO_description = "${${$go_desc{$motif_name}}[$i]}[0] => ${${$go_desc{$motif_name}}[$i]}[1] ";
							print HTML qq($GO_description</br >);
						}
					}else{
						for my $i(0..6) {
							my $GO_description = "${${$go_desc{$motif_name}}[$i]}[0] => ${${$go_desc{$motif_name}}[$i]}[1] ";
							print HTML qq($GO_description</br >);
						}
					}  
					print HTML qq(</p></div></td></tr>);
				}
			print HTML qq(<tr><td colspan="2" align="center" height="20" valign="middle"> <A HREF=\"$golist\">Genome-wide Motif$n search results</a></td></tr>);
			print HTML qq(<tr><td colspan="2" align="center" height="20" valign="middle"> <A HREF=\"$goview\">Motif$n gene list of over-represented annotation(s)</a></td></tr>
			</table>);
			$n++;
		}
		print HTML qq(<p>&nbsp;</p> <p><font size="-1"><center>Sequence logo generated by <a HREF="http://weblogo.berkeley.edu/">weblogo</a>
		<br><a HREF="http://search.cpan.org/dist/Chart-Clicker/">Graphic generated with Chart::Clicker,Perl module </a>
		<br><a HREF="http://www.plant.uoguelph.ca/research/homepages/raizada/index.html">Promzea program from the Raizada lab</a></center></font></p>
		</div></div></div></body></html>);
		close HTML;
		unlink ("./html/temp/$display.fas","./html/temp/$display.txt","./temp/error$display.txt");
	}else{
		diehtml("Promzea is not able to determine any motif.")
	};
}
#______________________________________________________________________________________________________
sub maize_motif{
	my ($organism, $promlength, $stringence, $seqfile_cDNA,$option,$sitehome,$siteimage)= @_;
	use POSIX qw(ceil);
	use Parallel::ForkManager;
	use Sys::Info; 
	use Sys::Statistics::Linux::MemStats;
	use 5.010;
	use List::Util qw(min); 
	#number of processors used 
    my $cpu_info = Sys::Info->new;
	my $cpu = $cpu_info->device('CPU');
	my $MAX_PROCESSES_FOR_CPU = $cpu->count;
	my $MAX_PROCESSES = $MAX_PROCESSES_FOR_CPU-1;
	#	Initialize the name for the directory of temp files
	my($display,$usage,$infile,@return_motif_results,%consensus_seq,%go_descrip);
	my @file = split (/\//, $seqfile_cDNA);
	my $intermediate = $file[$#file] ;
	@file = split (/\./, $intermediate);
	$display = $file[0];
	my $seqfile =  "$display".'_upst.fasta';
	#	Create a directory for the processing file 
	mkdir("$sitehome/html/results/$display", 0777); #|| print $!;
	mkdir("$sitehome/html/results/$display/go", 0777); #|| print $!;
	mkdir("$sitehome/html/results/$display/logo", 0777); #|| print $!;
	#	Initialization of the root paths
	my $cgi_path_temp = "$sitehome/html/temp";	
	my $dbpath = "$sitehome/html/data/"."$organism";
	my $results = "$sitehome"."/html/results/"."$display";
	my $results2 = "html/results/"."$display";
	my $affy = "affyMaize_probe_to_genome.txt";

	#	Data folder paths and file
	my $memebkg;	
	if($promlength =~ /intron|utr5|utr3/){
		$memebkg = "$dbpath/"."$organism"."\_"."200"."\_meme_bg2\.txt";
	}else{
		$memebkg = "$dbpath/"."$organism"."\_"."$promlength"."\_meme_bg2\.txt";
	}
	#my $cDNA_flat_db = "$organism\_"."premRNA.fasta";
	my $cDNA_flat_db = "$organism\_"."translations.fasta";
	my $genome_flat_db = "$dbpath"."/"."$organism\_"."genome.fas";
	my $prom_flat_db  = "$organism"."\_"."$promlength"."\_"."upstream.fa";
	my $prom_file_id = "$dbpath"."/"."$prom_flat_db";
	my $rand_prom_file = "$dbpath"."/"."random"."$prom_flat_db";
	my @answers = run("cp $rand_prom_file $results/",".");
	my $rand_genome = "$results"."/"."random"."$prom_flat_db";
	my $annotation_count_total = "$dbpath/"."$organism"."\_annotation_count\.txt";
	my $annotation_txt = "$dbpath/"."$organism"."\_GOannot.txt";
	my (@AoS,@mtfPWMs,$motifname,$skipline,$n,@in,$site,$score,$nb_seq,$seqmt,$mfconsensus,$mtstrand,$mtstart,$mtend);
	my ($width,$i,$j,@A,@C,@G,@T,$seq,$p_value,$motif,$A,$C,$G,$T,$threshold,$consensus);
	# hash containing annotation of all the genes in the genome
	my $gene_to_annot = annotation_hash($annotation_txt);
	my %gene_to_annot = %$gene_to_annot;
	# 	Loading of Promoter sequences from cDNA file	 	     
	if ($option == 1){
		cDNA_to_promoter($organism,$stringence,$results,$cgi_path_temp,$dbpath,$cDNA_flat_db,$prom_flat_db,\%gene_to_annot,$intermediate);
		}
	elsif ($option == 2){
		affy_to_promoter($results,$cgi_path_temp,$dbpath,$cDNA_flat_db,$prom_flat_db,\%gene_to_annot,$intermediate);
		}
	elsif ($option == 3){
		geneid_to_promoter($results,$cgi_path_temp,$dbpath,$cDNA_flat_db,$prom_flat_db,\%gene_to_annot,$intermediate);
		}
	elsif ($option == 4){
		chip_seq_peaks($results,$organism,$genome_flat_db,$cgi_path_temp,$intermediate);
		}
	else{ print "No option"}
	my $file = "$results/$seqfile";

	#	Creation of input file for BioProspector		       
	my $BioPrin = "$results/BioPr$display.txt";
	my $BioProut = "$results/BioPr$display.out";
	my $random_geno = "$results/promoter_random$display.fa";	
	

	#	Create one SeqIO object to read in,and another to write out
	my $list_id = fasta_to_BioPfasta($BioPrin,$file);
	my @list_id = @$list_id;

	#	Transform in pscan compatible files 
	my $scanin = "$results/inBioPr$display.txt";	
	my $count = fasta_to_fastaline($scanin,$file);
	$nb_seq = $count;

	# 	Gather all the promoters into an array - except those from the input set of genes
	my $nb_random = 4000;
	my $min_promoter_length = 100;
	my %promoters;
	if ($option == 4){
		my $promoters = seq_selection($file);
		%promoters = %$promoters;
	}else{
		my $promoters = promoter_selection($prom_file_id,$min_promoter_length);
		%promoters = %$promoters;
	}
	#while( my ($k, $v) = each %promoters ) {
        #print "key: $k,"."length:".length($v)."\n";
    	#}
	@answers = run("fasta-subsample $prom_file_id $nb_random > $random_geno",".");

	#	MEME motif search	
	my $AoS_ref;
	my $mtfPWMs_ref;
	($AoS_ref,$mtfPWMs_ref) = meme_run_parse(\@AoS,\@mtfPWMs,$results,$display,$file,$memebkg,$nb_seq,$scanin,$random_geno,$nb_random,$MAX_PROCESSES);
	@AoS = @$AoS_ref;
	@mtfPWMs = @$mtfPWMs_ref;

	#	BioProspector motif search(BioP) 
	($AoS_ref,$mtfPWMs_ref) = biop_run_parse(\@AoS,\@mtfPWMs,$results,$display,$BioPrin,$nb_seq,$scanin,$random_geno,$nb_random);
	@AoS = @$AoS_ref;
	@mtfPWMs = @$mtfPWMs_ref;
	
	#	Weeder motif search			    		
	($AoS_ref,$mtfPWMs_ref) = weeder_run_parse($organism,\@AoS,\@mtfPWMs,$results,$display,$file,$scanin,$nb_seq,$random_geno,$nb_random);
	@AoS = @$AoS_ref;
	@mtfPWMs = @$mtfPWMs_ref;
	unlink($BioPrin,$BioProut,$seqfile,$scanin,"$file.mix","$file.html");

	# 	STAMP plugging motif file writing
	my %stamp;
	for $i( 0 .. $#mtfPWMs){
		$stamp{"$mtfPWMs[$i][0]"} = "A $mtfPWMs[$i][1][0]\nC $mtfPWMs[$i][1][1]\nG $mtfPWMs[$i][1][2]\nT $mtfPWMs[$i][1][3]\n";
	}
	
	#	Motif filtering and transform data format		 	 
	my (@name,$numCols,$row,$sum);
	my ($nbseq,%motifs,@AoA,@sortedAoA,$length,@Array);
	my (@rank,$range,$range1,$startpos,$pos,$seq_length,$nb,$key);
	my ($line,$nline,$name,@goinput,@rowmtrx,@GOanalysis,$Goinput,$gofile);
	
	#	Initialization of hash to create graph ouput
	for $i ( 0 .. $#mtfPWMs){
		$motifname = $mtfPWMs[$i][0];
		$motifs{$motifname} = "motif$i";
	}
	
	my ($Aline,$Gline,$Cline,$Tline);
	for $i ( 0 .. $#mtfPWMs){ 
		#	JASPAR file generation: transform PWM from cluster output to input file for Clover program
		$motifname = $mtfPWMs[$i][0];
		$Aline = $mtfPWMs[$i][1][0];@A = split (/\s+/,$Aline);
		$Cline = $mtfPWMs[$i][1][1];@C = split (/\s+/,$Cline);
		$Gline = $mtfPWMs[$i][1][2];@G = split (/\s+/,$Gline);
		$Tline = $mtfPWMs[$i][1][3];@T = split (/\s+/,$Tline);
		my $jaspar = "$results/"."jaspar$motifname.txt";
		open (OUT, "+>", $jaspar);
		print OUT ">$motifname\n";
		foreach $i(0..$#A){
			print OUT "$A[$i]\t$C[$i]\t$G[$i]\t$T[$i]\n";
		}	
		close OUT;
		#	create hash of the motif consensus by transforming PWM from cluster output
		my $consensus_seq_ref = PWM_to_ConsensusHash($motifname,\%consensus_seq,\@A,\@C,\@G,\@T);
		%consensus_seq = %$consensus_seq_ref;
	}
	
	# 	Add the length of the promoter sequences in array of array in case of promoter having a gap at their starting point
	#	Determine motif occurence and representation in the genome			  			
		
	my @idAoS = sort{ $a->[1] cmp $b->[1] } @AoS;
	my @nbseq = ();
	#print Dumper (\@idAoS);
	for $nb(0..$#idAoS){
		if($option == 4){
			$key = $idAoS[$nb][1];
			$length = length($promoters{$key});
			my $per = "$idAoS[$nb][2]"*100/"$length";
			$per += 0.5; 
			$idAoS[$nb][2] =  int($per);
			push (@nbseq,$key);
			push (@{$idAoS[$nb]},"$length");
		}else{
			$key = $idAoS[$nb][1];
			#print "$key\n";
			$length = length($promoters{$key});
			push (@nbseq,$key);
			if (exists($promoters{$key})){ 
				push (@{$idAoS[$nb]},"$length");
			}
		}
	}
	@nbseq = uniq(@nbseq);
	$nbseq = @nbseq;
	#print Dumper(\@idAoS);	
	
	#	Graph preparation: create axis 
	my ($x_axis_ref, $y_axis_ref, $coordinates_ref, $pt);
###
    my ($min,$max);
    my @intron;
    if($promlength =~ /intron|utr5|utr3/){
        for $i (0 .. $#idAoS){
            my $temp_intron = $idAoS[$i][8];
            push(@intron,$temp_intron);
        }
        $max = max(@intron);
        $min = 0;
    }
    else{
        $min = -$promlength;
        $max = 0;
    }
    my $binRange;
    if ($promlength =~ /200|500|1000/){  
        $binRange = "$promlength"/10;
    }else{
        $binRange = max(@intron)/10;
    }
        
    my ($hash_ref, $rank_order) = axis_generation($min,$max,$binRange);
    my %x_axis = %$hash_ref;
    my @x_rank_order = @$rank_order;
    #my $pm = new Parallel::ForkManager($MAX_PROCESSES);
####	
	foreach $motifname(keys %motifs){
		#my $pid = $pm->start and next;
		my @ycoord = ();
        #@rank = (0) x $#coordinates;
		$name= $motifname;
		#	Weblogo files preparation
		my $logo = "$results/logo/"."$motifname"."logo.fasta";
		open (LOGO, "+>", $logo);
		#$line=0 represent the file headers
		for $line(1..$#idAoS){
			$nline= $idAoS[$line][0];
			if ("$nline" !~ "$name"){}
			else{
				#	Control if the line has the good motif name to include it in graph
				$length = $idAoS[$line][8];
				$startpos = $idAoS[$line][2];
                my $location;
                if($promlength =~ /intron|utr5|utr3/){
                    $location = $startpos;
                }else{
                    $location = $startpos - $length;
                }
                push(@ycoord,$location);
				##	Seqlogo fasta file:
				print LOGO '>'."$idAoS[$line][1]\n$idAoS[$line][6]\n";
			}
		}
		close LOGO;
        #print "y-coord\n";print $_."\t" foreach (@ycoord);print "\n";
        my ($x_axis,$y_axis) = FrequencyTable(\%x_axis,\@ycoord,\@x_rank_order,$binRange);
        my @x_axis = @$x_axis;
        my @y_axis = @$y_axis;
		
		#	WebLogo grahics execution
		my $lpng = "$siteimage/"."$motifname"."logo";
		my @answers = run("seqlogo -F PNG -E 1 -f $logo -k 1 -o $lpng -abceMnY",".");
		# 	Motif frequency Graphic 
		my $graph = "$siteimage/sequence$motifname.png";
		motif_position_freq_graph_gnu(\@x_axis,\@y_axis,$graph,$binRange);
###
        #$pm->finish; # Terminates the child process
###		
	}
	#$pm->wait_all_children;
	my (%arrayMNCP,@bgd_count, @sortedbdg_count, @MNCP_bgd);
	my (@fgd_count, @sortedfdg_count, @MNCP_fgd);
	
	#	load gene description file for annoationgene list file 
	
	my $nbseq_genome = count_seq($prom_file_id);
    
    
	foreach $motifname(keys %motifs){
		# Forks and returns the pid for the child:

		#	Foreground arrays for motif metric - MnCP score 
		#	Run and Parse clover outputs
		@MNCP_bgd = ();
		@MNCP_fgd =();
		my $jaspar = "$results/"."jaspar$motifname.txt";
		my ($AoU,$nmot) = clover_parsing ("$jaspar","$file");
		my @AoU= @$AoU;
		my @idAoU = sort{ $a->[1] cmp $b->[1] } @AoU;
		$nbseq = $nb_seq;
		for $line(0..$#idAoU){
			#	MNCP metric preparation
			push(@fgd_count,$idAoU[$line][1]);
		}
		$line = 1;
		@sortedfdg_count = sort(@fgd_count);
		$count=1; 
		for $i(0..$#sortedfdg_count){
			if (defined($sortedfdg_count[$i+1])){
				if("$sortedfdg_count[$i]" eq "$sortedfdg_count[$i+1]"){	$count++;}
				else{	push (@MNCP_fgd,$count);	$count=1;}
			}
			else{# last line
				if("$sortedfdg_count[$i]" eq "$sortedfdg_count[$i-1]"){ push (@MNCP_fgd,$count);}
				else{	$count = 1;	push (@MNCP_fgd,$count);}
			}
		}
		my $nomatch = $nbseq - $nmot;
		my @nomatch = (("0")x$nomatch);
		if(@nomatch){
			@MNCP_fgd = (@MNCP_fgd, @nomatch);
		}
		push (@{$arrayMNCP{$motifname}}, [@MNCP_fgd]);
		unlink($AoU,@AoU);
		#	Annotation Background arrays for motif metric - MnCP score 
		#	Run and Parse clover outputs
		$jaspar = "$results/"."jaspar$motifname.txt";
		$nmot = 0;
		my $nbfile_clover = ceil($nbseq_genome/5000); 
		for my $i(1..$nbfile_clover){ #number of files promoter sequence are split into
		# Clover can processs more than 5000 sequences at the time because of memory allocation problem			
			my $prom_genome = $dbpath."\/$organism\_$promlength\_up\_$i\.fasta";
			my($AoA_temp,$nmot_temp) = clover_parsing("$jaspar","$prom_genome");
			#($AoA,$nmot) = clover_parsing("$jaspar","$rand_genome");
			#my @AoA = @$AoA;
			my @AoA_temp = @$AoA_temp;
			$nmot = $nmot + $nmot_temp;
			push(@AoA,@AoA_temp);
		}
		#	Initialization for motif MNCP metric parameters
		@bgd_count=();
		@sortedbdg_count = ();
		@MNCP_bgd = ();
		#	Return array
		push (@return_motif_results,$motifname);
		#	GO annotation input reading
		my $golist = "$results/go/"."$motifname"."_golist.txt";
		open (GOLIST,"+>",$golist)|| die ( "could not open GO output");
		print GOLIST "If the sequence of the promoter is shorter than 200-500-1000 pb,the motif localization can be approximate\.\nIt is advised to check the motif real coordinates in the promoters of genes of interest\.The motif coordinates are accurate for first-intron\.\n
		\nGene ID\tCoord.\tScore\tAnnotation\t\n";
		for my $line(0..$#AoA){		
			#	GO annotation array for the motif 
			@rowmtrx = ($AoA[$line][0],$AoA[$line][1], $AoA[$line][2], $AoA[$line][3]);
			push ( @goinput,[@rowmtrx]);
			#	MNCP metric preparation
			push(@bgd_count, $AoA[$line][1]);	
		}
		#	create Go annotation
		#-----------------------------------------------
		@goinput = sort{ $b->[0] cmp $a->[0] }@goinput;
		my ($gene,$location,$locationend);
		for my $row(0..$#goinput){ 
			$gene = $goinput[$row][1];
			if(($promlength =~ /intron|utr5|utr3/) or ($option == 4)){
				$location = $goinput[$row][2];
				$locationend = $goinput[$row][3];
			}
			else{$location = $goinput[$row][2]-"$promlength";
				$locationend = $goinput[$row][3]-"$promlength";
			} 
			my @des = (); 
			if (defined($gene_to_annot{$gene})) {	 
				for my $j(0..$#{$gene_to_annot{$gene}}){
					push @des,$gene_to_annot{$gene}[$j][1];
				}
				my $temp = join( '|', @des );
				print GOLIST "$gene\t$location\t$locationend\t$temp\n"
			}
			else{print GOLIST "$gene\t$location\t$locationend\tNo Description\n"}
		}
		@goinput = ();
		close GOLIST;
		
		#	Generation of GO overepresentation files
		$gofile = "$results"."/go/"."$motifname"."_gofile.txt";
		# All genome id for GO annotations search separate by 5000 gene as Clover can not process everything at once:		   
		
		my($group, $group_total, $pval,$logp,@pvals,@groups,@go_descrip_mid,@GOfile);
		my @gene_to_motif = uniq(@bgd_count);
		my $nb_gene_w_mf = @gene_to_motif;#number of genes output of clover
		my $annotation_count_total = $annotation_count_total;
		my $annotation_count = annotation_counting(\%gene_to_annot,$annotation_count_total,\@gene_to_motif);
		my %annotation_count = %$annotation_count;
		my $sig = GO_hypergeometric(\%annotation_count,$nb_gene_w_mf,$nbseq_genome);
		my @sig = @$sig;
		# GO annotation pie chart: legend and data preparation
		open (GFILE, "+>",$gofile);
		print GFILE "# GO term\tp-value\tAnnotation description\n"; 		
		for my $i(0..$#sig){
			$group = $sig[$i][0];
			$pval = $sig[$i][1];
			$group_total= $sig[$i][2];
			if($pval == 0){$pval = 1/800;}					
			$logp= - log($pval);
			push (@go_descrip_mid, ["$group","$group_total"]);
			push (@groups, $group);
			push (@pvals, $logp);
			my $GOfile = $annotation_count{$group}[3];
			#print "$GOfile\n";
			@GOfile = @$GOfile;
			my $val = sprintf '%.2e',$sig[$i][1];
			print GFILE "##".$sig[$i][0]."\t".$val."\t".$sig[$i][2]."\n";
			print GFILE "## Gene ID:\n";
			for my $j(0..$#GOfile){
				print GFILE $GOfile[$j]."\n";
			}
			
		}
		close GFILE;
		$go_descrip{$motifname}= [ @go_descrip_mid ];
		@go_descrip_mid =();
		#	GO annotation pie chart: graphic execution
		if($#groups>=6){
			for $n(0..6){@groups = ( $groups[0],$groups[1],$groups[2],$groups[3],$groups[4],$groups[5],$groups[6],"others");}
			$sum = 0;
			for(@pvals){$sum += $_;}
			my $rem_logp = ($sum*6)/100;
			@pvals = ( $pvals[0],$pvals[1],$pvals[2],$pvals[3],$pvals[4],$pvals[5],$pvals[6],$rem_logp);
		} 
		else{};
		if(@groups){GO_annotation_graph($siteimage,$motifname,\@pvals,\@groups);}
		else{ 
			@groups = ("No description");
			@pvals = ("1");
			GO_annotation_graph($siteimage,$motifname,\@pvals,\@groups);
		}
		#	MNCP score background preparation
		#-----------------------------------------------
		#	Background arrays for motif metric
		@sortedbdg_count = sort(@bgd_count); 
		$count=1; 
		for $i(0..$#sortedbdg_count){
			if (defined($sortedbdg_count[$i+1])){
				if("$sortedbdg_count[$i]" eq "$sortedbdg_count[$i+1]")	{ $count++; }
				else{push (@MNCP_bgd,$count); $count=1;}
			}
			else{
				#	Last line procesing
				if("$sortedbdg_count[$i]" eq "$sortedbdg_count[$i-1]")	{push (@MNCP_bgd,$count);}
				else{$count = 1; push (@MNCP_bgd,$count);}
			}
		}	
		
		$nomatch = $nbseq_genome - $nmot;
		@nomatch = (("0")x$nomatch);
		if (@nomatch){@MNCP_bgd = (@MNCP_bgd, @nomatch);}
		push(@{$arrayMNCP{$motifname}}, [@MNCP_bgd]);
		unlink(@AoA,$nbseq,$jaspar,$seqfile_cDNA);	
	}
	unlink($random_geno);
	return ( \@return_motif_results, \%arrayMNCP,\%consensus_seq,\%go_descrip,\%stamp);	
}
#__________________________________________________________________________________________
sub cDNA_to_promoter{
	    my ($organism,$stringence,$results,$input_path,$path,$db,$db_prom,$descri_gene,@seqfiles)= @_;
		#------------------------------------------
		#       Variable Initialization
		#------------------------------------------
	    use List::Util qw(sum);
	    my ($organsim, $filename, @params, $seq_file,$cDNA_id,$format, $output_file, $prefix_out, $seq_out, @names, $seq);
	    my ($seqfiles,@idgene,$gene_id,@Aoblast,@stAoblast,@line,$line,@ID_list,@ID_list_gene,@unique_gene_ID,@des);
	    my ($eval,$perIdentity);
	    my $site_temp = "../html/temp/";
		#------------------------------------------
		#       Flat database loading
		#------------------------------------------
		#print qq(<p>Input was BLASTed against full-length cDNAs(maizegenome.org)</p>);

		#       cDNA sequences files
	    chomp(@seqfiles);
	    $filename = "@seqfiles";
	    @names = split /\./,"$filename";
	    
		#------------------------------------------
		#       Output file Initialization
		#------------------------------------------
	    $prefix_out = shift(@names);
	    $format = '_upst.fasta';
	    $output_file = $prefix_out.$format;
	    my $filepath = "$results/$output_file";
	    #$seq_out = Bio::SeqIO->new( -format => 'fasta', -file => ">$filepath");
	    open (ID, "+>>", "$results/idfile.txt");

		#-------------------------------------------------
		#       BLAST cDNA against fl-cDNA  maizesequence.org
		#-------------------------------------------------
	    my %selected_genes;
	    my %gene_to_annot = %$descri_gene;
	    $seq_file = new Bio::SeqIO(-format => 'fasta', -file => "$input_path/$filename");
	    while ($seq = $seq_file->next_seq){
		my $probe = $seq->display_id;
		my $nucl_seq = $seq ->seq;
		#       Determine query sequence lenght
		my $size = $seq->length();
		my $probe_seq = ">$probe\n$nucl_seq\n";
		open ( SEQ , "+>", "$results/$probe.txt");
        print SEQ $probe_seq;
		close SEQ;
		my $dirdb = "$path/$db";
		my $temp = "temp"."$probe".".out";
		#my @answers = run("blastall -p blastn -d $dirdb -i ../html/results/$prefix_out/$probe.txt -o ../html/results/$prefix_out/$temp -e 1e-35 -m8",".");
		my @answers = run("blastall -p blastx -d $dirdb -i ../html/results/$prefix_out/$probe.txt -o ../html/results/$prefix_out/$temp -e 10 -m8",".");
		#@answers = run("blastall -p blastx -d $dirdb -i ../html/results/$prefix_out/$probe.txt -o ../outtest.txt -e 10 -m8",".");
		open (IN,"<","$results/$temp") ; #|| dieemail("impossible to open my temporary blast output($results/$temp)");
		while(my $line =<IN>){
		    if (defined($line)){
			@line = split (/\t/,$line);
			push (@Aoblast, [@line]);
		    }
		}
		close IN;
				
		foreach my $i(0..$#Aoblast){
		# queryId(0), subjId(1), %Identity(2), alnLength(3), mismatchCount(4), gapOpenCount(5), 
		# queryStart(6), queryEnd(7), subjStart(8), subjEnd(9), eVal(10), bitScore(11)
			#print "$Aoblast[$i][0]\t$Aoblast[$i][1]\t$Aoblast[$i][2]\t$Aoblast[$i][10]\n";
			$perIdentity  = $Aoblast[$i][2];
			$eval = $Aoblast[$i][10];
			#       Compare alignment to 90% of query identity and length equal to 70%
		  	if($stringence = 2 && $perIdentity > 85 && $eval < 1e-50){
				#print "identity $perIdentity; evalue $eval\n";
				$cDNA_id =  $Aoblast[$i][1];
				if ($organism =~ /ZM/){
					if ($cDNA_id =~ /\_P/){
					    @idgene = split (/\_/, $cDNA_id);
					    $gene_id = $idgene[0];
					}
					elsif ($cDNA_id =~ /\_FGP/){
					    $gene_id =~ s/\_FGP/\_FG/;
					}
				}
				elsif ($organism =~ /AT/){
					@idgene = split (/\./, $cDNA_id);
					$gene_id = $idgene[0];
				}
				$selected_genes{$gene_id} = 1;		
				push(@ID_list,$gene_id);
				push(@ID_list_gene,$gene_id);
				@unique_gene_ID = uniq(@ID_list_gene);
			}
			elsif ($stringence = 1 && $perIdentity > 50 && $eval < 1e-10){
				$cDNA_id =  $Aoblast[$i][1];
				if ($organism =~ /ZM/){
					if ($cDNA_id =~ /\_P/){
					    @idgene = split (/\_/, $cDNA_id);
					    $gene_id = $idgene[0];
					}
					elsif ($cDNA_id =~ /\_FGP/){
					    $gene_id =~ s/\_FGP/\_FG/;
					}
				}
				elsif ($organism =~ /AT/){
					@idgene = split (/\./, $cDNA_id);
					$gene_id = $idgene[0];
				}
				$selected_genes{$gene_id} = 1;		
				push(@ID_list,$gene_id);
				push(@ID_list_gene,$gene_id);
				@unique_gene_ID = uniq(@ID_list_gene);
				
			}
		}
		my $temp1 = join(" ",@unique_gene_ID); 
		my $temp2 = "";
		for my $i(0..$#unique_gene_ID){
			my $gene = $unique_gene_ID[$i];
			if (defined($gene_to_annot{$gene})) {	 
				for my $j(0..$#{$gene_to_annot{$gene}}){
					push @des,$gene_to_annot{$gene}[$j][1];
				}
				@des = uniq(@des);					
				$temp2 = join( ';', @des );
			}
		}
		if (defined ($temp1) && defined ($temp2)){
			print ID $seq->display_id,"\t$temp1\t",$temp2,"\r\n";
		}elsif(defined ($temp1)){
			print ID $seq->display_id,"\t$temp1\t","No description","\r\n";
		}else{  print ID $seq->display_id,"\tNo match found\t"," ","\r\n";
		}
		@ID_list_gene=();
		@unique_gene_ID=();
		my @nothing = uniq();
		@des = ();		
		@Aoblast = ();
		unlink ("$results/$probe.txt","$results/$temp",$temp1, $temp);
	}
	my @uniqueID = uniq(@ID_list);
	my $path_dbprom = "$path/$db_prom";
	my $seq_in = &fasta_to_hash("$path_dbprom");
	open (OUT,"+>",$filepath);
	for my $i(0..$#uniqueID){
		my $k = $uniqueID[$i];
		if( exists $$seq_in{$k}){
			my $seq = $$seq_in{$k};
			if(length($seq)>50){
				print OUT "\>$k\n$seq\n";}
			else{}
		}
	}
	close OUT;	
	return \%selected_genes;
}
#__________________________________________________________________________________________
sub affy_to_promoter{
	my ($results,$input_path,$path,$db,$db_prom,$descri_gene,@seqfiles)= @_;
	use Bio::DB::Fasta;
	use Bio::SeqIO;
	use Bio::SeqIO::fasta;
	use Bio::SearchIO;
	#	Variable Initialization
	my ($filename, $dbt,@params, $seq_file,$factory,$blast_report,$result,$cDNA_id,$hit,$evalue);
	my ($format, $output_file, $prefix_out, $seq_out, @names, $seq,$seqfiles,@idgene,$gene_id);
	#	cDNA sequences files
	chomp (@seqfiles);
	$filename = "@seqfiles";
	@names = split /\./,"$filename";
	#	Output file Initialization
	$prefix_out = shift(@names);
	$format = '_upst.fasta';
	$output_file = $prefix_out.$format;
	my $filepath = "$results/$output_file";
	#$seq_out = Bio::SeqIO->new( -format => 'fasta', -file => ">$filepath");
	open(ARRAY,"<","$path/affyMaize_probe_to_genome.txt");
	my @array_to_genome = <ARRAY>;
	close (ARRAY);
	open (ID, "+>>", "$results/idfile.txt");
	open ( IN,"<","$input_path/$filename");
	my ($line,@array_genome_line,$array,$n,@ID_list);
	while ($line= <IN>){
		if (defined($line)){
			chomp $line;
			$line =~ s/\cM|\r|\n//;
			foreach my $n (0..$#array_to_genome){
				$array = "$array_to_genome[$n]";
				chomp $array;
				$array =~ s/\cM|\r|\n//;
				@array_genome_line = split (/\t/,$array);
				if ($line =~ /$array_genome_line[0]/i){
					#	Iterate over each hit on the query sequence
					$cDNA_id = $array_genome_line[1];
					if ($cDNA_id =~ /\_T/){
						@idgene = split (/\_/, $cDNA_id);
						$gene_id = $idgene[0];
					}
					else { 
						$gene_id = $cDNA_id ;
						$gene_id =~ s/\_FGT/_FG/;
					}
					print ID "$line\t$gene_id\r\n";
					push(@ID_list,$gene_id);
				}
			}
		}	
	}
	my @uniqueID = uniq(@ID_list);
	my $path_dbprom = "$path/$db_prom";
        my $seq_in = &fasta_to_hash("$path_dbprom");
	open (OUT,"+>",$filepath);
	for my $i(0..$#uniqueID){
		my $k = $uniqueID[$i];
		if( exists $$seq_in{$k}){
			my $seq = $$seq_in{$k};
			if(length($seq)>50){
				print OUT "\>$k\n$seq\n";}
			else{}
		}
	}
}
#__________________________________________________________________________________________
sub geneid_to_promoter{
	my ($results,$input_path,$path,$db,$db_prom,$descri_gene,@seqfiles)= @_;
	use Bio::DB::Fasta;
	use Bio::SeqIO;
	use Bio::SeqIO::fasta;
	use Bio::SearchIO;
	#	Variable Initialization
	my ($filename,$cDNA_id,$format, $output_file, $prefix_out, $seq_out, @names, $seq,$seqfiles,@idgene,$gene_id);
	#	cDNA sequences files
	chomp (@seqfiles);
	$filename = "@seqfiles";
	@names = split /\./,"$filename";
	#	Output file Initialization
	$prefix_out = shift(@names);
	$format = '_upst.fasta';
	$output_file = $prefix_out.$format;
	my $filepath = "$results/$output_file";
	#$seq_out = Bio::SeqIO->new( -format => 'fasta', -file => ">$filepath");
	open (ID, "+>>", "$results/idfile.txt");
	open (IN,"<", "$input_path/$filename");
	my ($line,@array_genome_line,$array,$n,@ID_list);
	while ($line= <IN>){
		if (defined($line)){
			chomp $line;
			$line =~ s/\cM|\r|\n//;
			$cDNA_id = $line;
			if ($cDNA_id =~ /\_T/){
				@idgene = split (/\_/, $cDNA_id);
				$gene_id = $idgene[0];
			}
			else { 
				$gene_id = $cDNA_id ;
				$gene_id =~ s/\_FGT/_FG/;
			}
			print ID "$line\t$gene_id\r\n";
			push(@ID_list,$gene_id);
		}
	}	
	my @uniqueID = uniq(@ID_list);
	my $path_dbprom = "$path/$db_prom";
        my $seq_in = &fasta_to_hash("$path_dbprom");
	open (OUT,"+>",$filepath);
	for my $i(0..$#uniqueID){
		my $k = $uniqueID[$i];
		if( exists $$seq_in{$k}){
			my $seq = $$seq_in{$k};
			if(length($seq)>50){
				print OUT "\>$k\n$seq\n";}
			else{}
		}
	}
}
#__________________________________________________________________________________________
sub chip_seq_peaks{
	my ($results,$organism,$genome_flat_db,$cgi_path_temp,@seqfiles) = @_; 
	use Bio::DB::Fasta;
	use Bio::SeqIO;
	use Bio::Seq;
	use Bio::PrimarySeq;
	chomp (@seqfiles);
	my $filename = "@seqfiles";
	my @names = split /\./,"$filename";
	#	Output file Initialization
	my $prefix_out = shift(@names);
	my $format = '_upst.fasta';
	my $output_file = $prefix_out.$format;
	my $filepath = "$results/$output_file";	
	my $BEDfile= "$cgi_path_temp/$filename";
	
	my $db  = Bio::SeqIO->new(-file => $genome_flat_db , '-format' => 'Fasta');
	my $out = Bio::SeqIO->new(-file => ">$filepath" , '-format' => 'Fasta');

	my @array;
	open(IN,"<",$BEDfile);
	while(<IN>) {
	    chomp;
	    push @array,$_;
	}
	close IN;

	my $fasta_id = "";
	
	while( my $chr = $db->next_seq() ) {
	    my $longid = $chr->display_id;
	    my @temp = ();
	    #print $longid."\n";
	    if($organism =~ "ZM"){
		@temp = split(/\:/, $longid);
		$fasta_id = $temp[2];
		#print "$fasta_id\n";
	    }
	    elsif($organism =~ "AT"){
		@temp = split(/Chr/, $longid);
		$fasta_id = $temp[1];
		#print "$fasta_id\n";
	    }
	    elsif($organism =~ "OS"){
		@temp = split(/\s+/, $longid);
		$fasta_id = $temp[0];
		#print "$fasta_id\n";
	    }
	    for my $i(0..$#array){
		my ($chrname,$start,$stop)= split ("\t",$array[$i]);
		if("$fasta_id" =~ m/^$chrname$/){
		    #print "peak on chr: $chrname\t$longid\n";
		    my $end =  $chr->length;
		    if ($stop > $end){$stop = $end;}
		    my $seq = $chr->subseq($start,$stop);
		    my $display = "peak\_$i\|$chrname\:$start";
		    my $seqobj = Bio::PrimarySeq->new( -seq => $seq, -id  => $display);
		    $out->write_seq($seqobj);
		}
	    }
	}
}
#__________________________________________________________________________________________
sub fasta_to_fastaline{
	#Transform in pscan compatible files 
	my ($scanin,$file)= @_;
	open (OUT,"+>",$scanin);
	open (IN, $file);
	my $count = 0;
	open (IN, $file);
	while (<IN>){
		s/\cM|\r|^$//;
		if(!m/^$/){
			if(m/^(\>)/){$count++;print OUT;}
			else{tr/a-z/A-Z/;print OUT;}
		}
	}
	close IN;
	close OUT;
	return $count;
}
#__________________________________________________________________________________________
sub promoter_selection {
	my($prom_file_id,$min_promoter_length) = @_;
	my %promoter = (name => '', sequence => ''); # name, sequence
	my(@promoters,$selected_genes,%list_promoters,$name,$promoter_key);
	open (IN, "<",$prom_file_id);
	while (<IN>) {
		if (/^>(\S+)/) {
			$name = $_;
			chomp $name;
			$promoter_key = substr($name,1);
			if ($promoter{name}) {
				if (not defined $selected_genes->{$promoter{name}} 
					and length $promoter{sequence} >= $min_promoter_length
					#and $promoter{sequence} !~ /N{50,}/
					 ) {
					push @promoters, \%promoter;	
				}
			}
			$promoter{name} = $1;
			$promoter{sequence} = '';
		}
		elsif(/(\w+)/) {
			$promoter{sequence} .= uc($1);
			$list_promoters{$promoter_key}= $promoter{sequence};
		}
	}
	close IN;
	return \%list_promoters;
}
#__________________________________________________________________________________________
sub seq_selection {
	my($prom_file_id) = @_;
	my %promoter = (name => '', sequence => ''); # name, sequence
	my(@promoters,$selected_genes,%list_promoters,$name,$promoter_key);
	open (IN, "<",$prom_file_id);
	while (<IN>) {
		if (/^>(\S+)/) {
			$name = $_;
			chomp $name;
			$promoter_key = substr($name,1);
			if ($promoter{name}) {
				if (not defined $selected_genes->{$promoter{name}}) {
					push @promoters, \%promoter;	
				}
			}
			$promoter{name} = $1;
			$promoter{sequence} = '';
		}
		elsif(/(\w+)/) {
			$promoter{sequence} .= uc($1);
			$list_promoters{$promoter_key}= $promoter{sequence};
		}
	}
	close IN;
	return \%list_promoters;
}

#__________________________________________________________________________________________
sub fasta_to_BioPfasta{
	my($BioPrin,$file) = @_;
	use Bio::SeqIO;
	use Bio::SeqIO::fasta;
	my $seq_in = Bio::SeqIO->new('-file' => "<$file",'-format' => 'fasta');
	my (@list_id,$list_id_name,$line);
	#	Write each entry in the input file to the output file
	open (OUT,"+>",$BioPrin);
	while (my $inseq = $seq_in->next_seq){
		#	Adapt the file for Bioprospector input type
		$line = ">".$inseq->id."\n".$inseq->seq();
		$list_id_name =$inseq->id;
		push(@list_id,$list_id_name);
		print OUT "$line\n";
	}
	close OUT;
	return \@list_id;
}
#__________________________________________________________________________________________
sub meme_run_parse {
	my($Astring,$motifPWMs,$results,$display,$file,$memebkg,$nbseq,$scanin,$random_geno,$nb_random,$MAX_PROCESSES) = @_;		
	my @Astring = @$Astring;
	my @motifPWMs = @$motifPWMs;
	my @output = run("meme $file -nmotifs 10 -dna -maxw 10 -maxsize 10000000000 -p $MAX_PROCESSES -bfile $memebkg -revcomp -nostatus -text", ".");
	my $infile = $output[0];
###	
	#my $outme = "memeresults.txt";
	#open(MEOUT,">$outme");
	#print MEOUT "$output[0]";
	#close MEOUT;
###

	my(@response,$A,$C,$G,$T,@background,$n,@Temp,$motifname,$width,$site,$score,$threshold,$consensus,@A,@C,@G,@T);
	my($seqmt,$mtstrand,$mtend,$mfconsensus,@Base,@temp_meme,$mtstart,$p_value,$j);
	my $nsites = 0;
	my $i = 0;
	@response = split("\n",$infile);
	for(my $i=0;$i<=$#response;$i++){
		if($response[$i]=~/^Background/){
			$i++;
			@background = split /A\s+|C\s+|G\s+|T\s+/,$response[$i];
			chomp $background[4];	
		}
		else{}
	}
	$A=$background[1]; $C=$background[2]; $G=$background[3]; $T=$background[4];
	$n = 1;
	for(my $i=0;$i<=$#response;$i++){	
		next unless $response[$i]=~ /^MOTIF/;
		@Temp = split /MOTIF\s+|\s+width\s+\=\s+|\s+sites\s+\=\s+|\s+llr\s+\=\s|\s+E\-value\s+\=\s/,$response[$i];
		$motifname = "MEME".$n;
		$width = $Temp[2];
		$site = $Temp[3];
		$score = $Temp[5];
		chomp $score;
		$threshold = $site/$nbseq;
		do{$i++;}while($response[$i]!~ /^Multilevel/);
		@Temp = split (/\s+/,$response[$i]);
		$consensus = $Temp[1];
		#	Obtain the motif size
		do{$i++;}while($response[$i] !~ m/^Sequence\sname\s+Strand\s+Start/);
		$i++;$i++;
		@A= (("0")x$width); @C= (("0")x$width); @G= (("0")x$width); @T= (("0")x$width);
		do {	#	Read line with motif information for
			@Temp= split (/\s+|\t/,$response[$i]); 
			$seqmt = $Temp[0]; 
			$mfconsensus = $Temp[5];
			chomp $response[$i];
			#---------------------------------
			#	Use in case of single strand search
			#$mfconsensus = $_[5]; chomp; $mtstart = $_[1]; $mtstrand = "plus";$mtend = $_[1] + $width;
			#---------------------------------
			#	Use in case of dble strand (revcom)	
			if($Temp[1] =~ /\-/){
				$mtstrand = '-';	
				$mtstart = $Temp[2] - $width;
				$mtend = $Temp[2];
				#$Temp[5] =~ tr/ATGCatgc/TACGtacg/;				
				$mfconsensus = $Temp[5];
			}
			elsif($Temp[1] =~ /\+/){
				$mtstrand = '+';	
				$mtend = $Temp[2] + $width;
				$mtstart = $Temp[2];
				
			}
			#---------------------------------
			@Base = split (//,$mfconsensus);
			for my $f(0..($width-1)){
				if($Base[$f] eq 'A') {$A[$f]="$A[$f]"+1} 
				elsif($Base[$f] eq 'C'){$C[$f]="$C[$f]"+1} 
				elsif($Base[$f] eq 'G'){$G[$f]="$G[$f]"+1} 	
				elsif($Base[$f] eq 'T') {$T[$f]="$T[$f]"+1} 
			}
			@Temp = ($motifname,$seqmt,$mtstart,$mtend,$score,$n,$mfconsensus,$mtstrand);
			#	Unfiltered meme results
			push (@temp_meme,[@Temp]);
			$nsites++;
			$i++;
		} while($response[$i] !~ /^\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-/);
		my $fimoinput = "$results/"."$n"."fimo$display".".txt";
		
		my $fimores = input_MEMEtoFIMO($fimoinput,$A,$C,$G,$T,$motifname,\@A,\@C,\@G,\@T,$score);
		#	Calculate the p-value of binomial distribtion
		$p_value = hypergeometic_motif_MEME($results,$display,$scanin,$nbseq,$random_geno,$nb_random,$fimoinput,$memebkg);		
		if (($threshold >= 0.5)and($p_value < 0.05)){
			my @Temp1 = ("@A","@C","@G","@T");
			my @Temp2 = ("$motifname",\@Temp1);
			#@_ = ("$motifname",["@A","@C","@G","@T"]);
			#push (@motifPWMs,[@_]);
			push @motifPWMs,\@Temp2;
			for $j(0..$#temp_meme){push(@Astring,$temp_meme[$j]);}
			#push (@Astring,"@temp_meme");
			}	
		(@A,@C,@G,@T,@temp_meme) = ();
		$n++;
		unlink($fimoinput);
	}
 	return (\@Astring,\@motifPWMs);
}	
#__________________________________________________________________________________________
sub input_MEMEtoFIMO{
	my($fimoinput,$A,$C,$G,$T,$motifname,$Aarray,$Carray,$Garray,$Tarray,$score) = @_;	
	my($numCols,$sum,$width,$row,$colA,$colC,$colG,$colT);
	my(@A) = @$Aarray;
	my(@C) = @$Carray;
	my(@G) = @$Garray;
	my(@T) = @$Tarray;
	# 	FIMO input file generation: transform pwm from cluster output to input file for FIMO program
	open(FOUT, "+>",$fimoinput);
	print FOUT "MEME version 4.6.1\r\n\r\nALPHABET= ACGT\r\n\r\nstrands: + -\r\n\r\n";
	print FOUT "Background letter frequencies (from) \r\nA $A C $C G $G T $T \r\n\r\n";
	#print FOUT "MEME version 4.3.0\n\nALPHABET= ACGT\n\nstrands: + -\n\n";
	#print FOUT "Background letter frequencies\nA $A C $C G $G T $T\n\n";
	$numCols = 4;
	$sum = $A[1] + $C[1] + $G[1] + $T[1];
	$width = scalar(@A);
	print FOUT "MOTIF $motifname\r\nletter-probability matrix: alength= 4 w= $width nsites= $sum E= $score \r\n";
	#print FOUT "MOTIF $motifname\nBL   MOTIF $motifname width=$width seqs=0\nletter-probability matrix: alength= 4 w=$width nsites=$sum E=$score\n";
	foreach my $row(0..$#A) {
		if(defined $row){		
			if ($A[$row]==0){$colA = sprintf ("%.6f",0)} 
			else{$colA = sprintf ("%.6f",($A[$row]/$sum))};
			if ($C[$row]==0){$colC = sprintf ("%.6f",0)} 
			else{$colC = sprintf ("%.6f",($C[$row]/$sum))};
			if ($G[$row]==0){$colG = sprintf ("%.6f",0)} 
			else{$colG = sprintf ("%.6f",($G[$row]/$sum))};
			if ($T[$row]==0){$colT = sprintf ("%.6f",0)} 
			else{$colT = sprintf ("%.6f",($T[$row]/$sum))};
			print FOUT " $colA  $colC  $colG  $colT \r\n";
		}
	}		
	close(FOUT);
	my $res = "fimo_sucess\n";
	return ($res);
}
#__________________________________________________________________________________________
sub biop_run_parse{
	my($Astring,$motifPWMs,$results,$display,$BioPrin,$nb_seq,$scanin,$random_geno,$nb_random)= @_;
	my @Astring = @$Astring;
	my @motifPWMs = @$motifPWMs;
	my($infile,$n,@temp_biop,@answers,@response,@Temp,$motifname,$consensus,@A,@C,@G,@T,@in,$width,$score,$site,$threshold);
	my($seqmt,$mfconsensus,$mtstrand,$mtstart,$mtend,@Base,$p_value,$j);
	@answers = run("BioProspector -i $BioPrin -n 10 -r 5",".");
	if ($answers[2]==0){
		diehtml("BioProspector @answers failed: $?");
	}
	unlink("errors$display");
	if ( $? == -1 ){
		print "BioProspector failed: $!\n";
	}
	else{
		$infile = $answers[1];
		@response = split("\n",$infile);
		$n = 1;
		for(my $i=0;$i<=$#response;$i++){
			next unless  $response[$i]=~ /^Motif #/;
			$motifname = "BioP".$n;
			@Temp = split /\:\s+\(|\//,$response[$i];
			$consensus = $Temp[1];
			$i++;$i++;
			@Temp = split /\;/,$response[$i];
			#	Obtain the motif size
			@in = split /\s+\(|\,/,$Temp[0];
			$width =$in[1];
			@A= (("0")x$width); @C= (("0")x$width); @G= (("0")x$width); @T= (("0")x$width);
			#	Obtain the motif score
			@in = split / /,$Temp[2];
			$score = $in[2];
			@in = split / /,$Temp[3];
			$site = $in[2];
			$threshold = $site/$nb_seq;
			do{$i++;}while($response[$i] !~ /^\>/);
			do{
				@Temp= split (/\s+|\t/,$response[$i]);
				$seqmt = substr $Temp[0],1; 
				$i++;
				$mfconsensus =$response[$i];
				chomp $mfconsensus;
				chomp $response[$i];
				if ( $Temp[5]=~ /f/){ 
					$mtstrand = '+';
					$mtstart = $Temp[6];
					$mtend = $Temp[6] + length($mfconsensus);
				}
				elsif ($Temp[5]=~ /r/){ 
					$mtstrand = '-';
					$mtstart = $Temp[6] - length($mfconsensus);
					$mtend  = $Temp[6];
					#$mfconsensus =~ tr/ATGCatgc/TACGtacg/;
				}
				#	Create PWM:
				@Base = split (//,$mfconsensus);
				for my $f(0..($width-1)){
					if($Base[$f] eq 'A') {$A[$f]="$A[$f]"+1;} 
					elsif($Base[$f] eq 'C'){$C[$f]="$C[$f]"+1} 
					elsif($Base[$f] eq 'G'){$G[$f]="$G[$f]"+1} 	
					elsif($Base[$f] eq 'T') {$T[$f]="$T[$f]"+1} 
				}				
				@Temp = ($motifname,$seqmt,$mtstart,$mtend,$score,$n,$mfconsensus,$mtstrand);
				#	Unfiltered Bioprospector
				push (@temp_biop,[@Temp]);
				$i++;
			}
			while($response[$i]!~ /^\*\*\*/);
			my @Temp1 = ("@A","@C","@G","@T");
			my @Temp2 = ("$motifname",\@Temp1);
			#@_ = ("$motifname",["@A","@C","@G","@T"]);		
			my $pscan_input = "$results/"."$n"."pscaninput"."$display".".wil";
			open(OUT, "+>",$pscan_input);
			print OUT ">$motifname\n@A\n@C\n@G\n@T\n";
			close OUT;
			#	Calculate the p-value of binomial distribution of a motif	
			$p_value = binomial_Bioprospector_motif($results,$scanin,$nb_seq,$random_geno,$nb_random,$pscan_input);		
			if($p_value < 0.7){
				push (@motifPWMs,\@Temp2);
				for $j(0..$#temp_biop){ push(@Astring, $temp_biop[$j]);}
			}
			(@A,@C,@G,@T,@temp_biop) = ();
			$n++;
			unlink($pscan_input);
		}
		#close IN;
	}
	return (\@Astring,\@motifPWMs);			
}
#__________________________________________________________________________________________
sub weeder_run_parse{
	my($organism_bckg,$Astring,$motifPWMs,$results,$display,$file,$scanin,$nb_seq,$random_geno,$nb_random)= @_;
	my($seqin,$line,$pattern,$seqline,$gcount,$n_seq,@seq,@line,@motif,@temp_weeder,$infile,$skipline,@background,$i,$j,$p_value);
	my(@Base,@A,@C,@G,@T,$consensus,$threshold,$width,$mfconsensus,$mtstrand,$mtend,$mtstart,$seqmt,$score,$motifname,@Temp);
	my @Astring = @$Astring;
	my @motifPWMs = @$motifPWMs;
	chdir "/var/www/cgi-bin/";
    my @weedercommand = ("/var/www/cgi-bin/weederlauncher.out","$file","$organism_bckg","small","S","M");
    #print STDERR "/var/www/cgi-bin/weederlauncher.out $file $organism_bckg small S M";
	capture { qx(@weedercommand) };#, \undef, \undef;
        chdir "/var/www/html/";
	$infile = "$file.wee";
	#	Weeder parser
	open (IN,"<",$infile) or diehtml("Unable($?) to open file($infile).");
	my $n = 1;
	while(<IN>){
		next unless m/^Your/;
		$skipline = <IN> ;
		#	Put sequence in array to rename them in final array
		do{ 	$line=<IN>;
			chomp $line;
		    	@Temp = split /Sequence\s+|\s\:\s\>/,$line;
		    	push (@seq,[@Temp]);} while( $line !~ m/^$/);
		pop(@seq);
		do{$skipline = <IN>;chomp $skipline;} while($skipline !~ m/\*\*\*\s+Interesting\smotifs/);
		MOTIF:
		$skipline = <IN>;
		if (defined($skipline)){
			while($skipline !~ m/^A|^G|^C|^T/){
				$skipline = <IN>;chomp $skipline;
			}
			chomp $skipline;
			$consensus = $skipline;
			$threshold = low_complexity($consensus);
			$width = length($skipline);
			@A= (("0")x$width); @C= (("0")x$width); @G= (("0")x$width); @T= (("0")x$width);
			do{$skipline = <IN>;
				chomp $skipline;}
				while($skipline !~ m/^Seq\s+St\s+oligo/);
			do{ $line=<IN>;
				chomp $line;
				$line =~ s/^\s+//;
				if ( $line !~ m/^$/){
					@Temp = split (/\s+|\]\s+/,$line);
						$mfconsensus = $Temp[2];
						if ($mfconsensus =~ m/^\[/){
							$mfconsensus = substr($mfconsensus,1,)}
						chomp $mfconsensus;
						@Base = split (//,$mfconsensus);
						for $i(0..($width-1)){
							if($Base[$i] eq 'A') {$A[$i]="$A[$i]"+1;} 
							elsif($Base[$i] eq 'C'){$C[$i]="$C[$i]"+1} 
							elsif($Base[$i] eq 'G'){$G[$i]="$G[$i]"+1} 	
							elsif($Base[$i] eq 'T') {$T[$i]="$T[$i]"+1}
						}
						if($Temp[1] =~ m/\-/){
							$mtstrand = '-';	
							$mtend = $Temp[3];
							$mtstart = $mtend - "$width";
						}
						elsif($Temp[1] =~ m/\+/){
							$mtstrand = '+';	
							$mtstart = $Temp[3];
							$mtend = $mtstart + "$width";
						}
						for $i(0..$#seq){
							
							if($seq[$i][1] =~ $Temp[0]){
								$seqmt = $seq[$i][2];
								chomp $seqmt;
								$seqmt =~ s/[\n\r\s]+//g;
							}
							else {};
						}
					$score = 1;
					$motifname = "Weeder$n";
					@line = ($motifname,$seqmt,$mtstart,$mtend,$score,$n,$mfconsensus,$mtstrand);
					#print "$motifname,$seqmt,$mtstart,$mtend,$score,$n,$mfconsensus,$mtstrand\n";
					push (@temp_weeder,[@line]);			
				}
			} while( $line !~ m/^$/);
			#	Calculate the p-value		
			my $pscan_input = "$results/"."$n"."pscaninput"."$display".".wil";
			open(OUT, "+>", $pscan_input);
			print OUT ">$motifname\n@A\n@C\n@G\n@T\n";
			close OUT;
			$p_value = binomial_motif($results,$scanin,$nb_seq,$random_geno,$nb_random,$pscan_input);		
			if ($threshold < 0.8 and $p_value < 0.3){
				my @Temp1 = ("@A","@C","@G","@T");
				my @Temp2 = ("$motifname",\@Temp1);
				push @motifPWMs,\@Temp2;
				for $j(0..$#temp_weeder){push(@Astring, $temp_weeder[$j]);
				}
			}	
			$n++;
			@temp_weeder = ();
			unlink($pscan_input);
			do{$skipline = <IN>;chomp $skipline;}while($skipline !~ m/^\=\=\=\=\=/);
			#	Redundant motif not the best
			$skipline = <IN>;$skipline = <IN>;
			if (defined	($skipline)) {goto MOTIF}
			else{}
		}
		else{}
	}
	close IN;
	unlink ("$file.mix","$file.html","$file.wee");
	return (\@Astring,\@motifPWMs);
}
#__________________________________________________________________________________________
sub count_seq {
	my ($file) = @_ ;	
	my $nseq = 0;
	open(CIN,"<", $file);
	while(<CIN>){
		if(/\>/){$nseq++;}
	}
	close CIN;
	return $nseq;
}
#__________________________________________________________________________________________
sub clover_parsing {
	my ($mymotifs,$myseqs) = @_; 
	my @answers = run("clover -u 9.21 -t 0.00 $mymotifs $myseqs",".");
	my @results=split(/[\n\r]/,$answers[0],-1);
	my (@AoA,$val,$x);
	my $nseq = 0;
	my $n = 0;
	for $x(0..$#results) {
		$val = $results[$x];
		chomp $val;
		my $seqmt = $val;
		$seqmt =~ s/\>//;
		if($val =~ /^\>/){
			$n++;
			$x++;
			until($val =~ /^$/){
				$val = $results[$x];
				chomp $val;
				$val=~ s/\r//;
				if($val !~ m/^$/){
					if ($val =~ /^\s+/){
                        			$val =~ s/\s+//;
                    			}
					my($motif,$start,$no,$end,$strand,$mfconsensus,$score) = split(/\s+/, $val);
					#print "$sequence\n";
					push( @AoA, [$motif,$seqmt,$start,$end,"1",$score,$mfconsensus,$strand]);
				}
				$x++;
			}
		}
	}
return (\@AoA,$n);
}
#__________________________________________________________________________________________
sub PWM_to_ConsensusHash {	
	my ($motifname,$consensus_seq,$A,$C,$G,$T) = @_;
	
	use List::Util qw(max);
	
	my $numCols = 4;
	my($colA,$colC,$colG,$colT,$consensus);
	my %consensus_seq = %$consensus_seq;
	my @A = @$A;
	my @C = @$C;
	my @G = @$G;
	my @T = @$T;
	my $sum = $A[1] + $C[1] + $G[1] + $T[1];
	for my $row(0..$#A) {
		if ($A[$row]==0){$colA = 0;} 
		else{$colA = $A[$row]/$sum};
		if ($C[$row]==0){$colC = 0;} 
		else{$colC = $C[$row]/$sum};
		if ($G[$row]==0){$colG = 0;} 
		else{$colG = $G[$row]/$sum};
		if ($T[$row]==0){$colT = 0;} 
		else{$colT = $T[$row]/$sum};
		my @rowMatrix = ($colA, $colC,$colG, $colT);
		if ( max(@rowMatrix) <= 0.25){$consensus = $consensus."N";}
		elsif($colA == max(@rowMatrix) && $colA == $colC){$consensus = $consensus."M";}
		elsif($colA == max(@rowMatrix) && $colA == $colT){$consensus = $consensus."W";}
		elsif($colA == max(@rowMatrix) && $colA == $colG){$consensus = $consensus."R";}
		elsif($colC == max(@rowMatrix) && $colC == $colG){$consensus = $consensus."S";}
		elsif($colC == max(@rowMatrix) && $colC == $colT){$consensus = $consensus."Y";}
		elsif($colT == max(@rowMatrix) && $colT == $colG){$consensus = $consensus."K";}
		elsif($colA == max(@rowMatrix)){$consensus = $consensus."A";}
		elsif($colC == max(@rowMatrix)){$consensus = $consensus."C";}
		elsif($colG == max(@rowMatrix)){$consensus = $consensus."G";}
		elsif($colT == max(@rowMatrix)){$consensus = $consensus."T";}
	}
	$consensus_seq{$motifname} = "$consensus";
	return(\%consensus_seq);		
}
#__________________________________________________________________________________________
sub axis_generation{
    my ($min,$max,$range) = @_;
    my (%hash,@orderAxis,$i); 
    #print "$max\t$min\n";
    for($i = $min; $i<= $max;$i = $i+$range){
        #print $i."\n";
        $hash{$i} = 0;
        push(@orderAxis, $i);
    }
    return(\%hash,\@orderAxis);
}
#__________________________________________________________________________________________
sub FrequencyTable{
    my ($xaxis_ref,$ycoord,$rank_order_ref,$binRange)=@_;
    my %hash = %$xaxis_ref;
    my @ycoord = @$ycoord;
    my @rank_order = @$rank_order_ref;
    
    #foreach (@ycoord){
    #    print $_."\n";
    #}
    
    my(@x_axis,@y_axis);
    foreach my $coord(@ycoord){
        #print "coordinates\t: $coord\n";
        foreach my $range(keys %hash){
            if ($coord >= $range and $coord < ($range+$binRange)){
                $hash{$range}+= 1;
            }
            #print "$range\t$hash{$range}\n";
        }
    }
    
    foreach my $rank(@rank_order){
        push (@x_axis,$rank);
        push (@y_axis, $hash{$rank});
    }
    return(\@x_axis,\@y_axis);
}
#__________________________________________________________________________________________
sub motif_position_freq_graph{
	my ($x_axis_ref,$y_values_ref,$graph,$range) = @_;
	use Chart::Clicker;
	use Chart::Clicker::Context;
	use Chart::Clicker::Data::DataSet;
	use Chart::Clicker::Data::Marker;
	use Chart::Clicker::Data::Series;
	use Chart::Clicker::Decoration::Legend::Tabular;
	use Chart::Clicker::Renderer::Area;
	use Chart::Clicker::Renderer::Pie;
	use Geometry::Primitive::Rectangle;
	use Geometry::Primitive::Circle;
	use Graphics::Color::RGB;
	my @x_axis = @$x_axis_ref;
    my @y_values = @$y_values_ref;
	#####
    my (@x_middle_range,$middle_range);   
    for my $i(0..($#x_axis)){
        $middle_range = $x_axis[$i] + ($range/2);
        push(@x_middle_range,$middle_range);
    }
    #push (@x_middle_range, "0");
    my @y_axis_range = (0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150);
    my $y = 0;
    my $ymax;
    foreach $y(@y_axis_range){
        $ymax = $y;
        last if(max(@y_values) < $y);
    }
    my $binRange = $ymax/5;
    my $min =0;
    my $max = $ymax;
    my ($hash_ref, $rank_order) = axis_generation($min,$max,$binRange);
    my @y_axis = @$rank_order;
    foreach my $i(0..$#x_middle_range){print "$x_middle_range[$i]\t"};print "\n"; #########################
    foreach my $i(0..$#y_values){print "$y_values[$i]\t"};print "\n";#########################
    #####
    #keys => \@x_axis,
    my $cc = Chart::Clicker->new(width => 480, height => 370, format => 'png');
	my $series1 = Chart::Clicker::Data::Series->new(
                        keys => \@x_middle_range,
                        values => \@y_values,
						name => ""
						);
	my $ds = Chart::Clicker::Data::DataSet->new(series => [$series1]);
	$cc->title->text("Motif frequency in user input promoters");
	$cc->title->padding->bottom(5);
	$cc->add_to_datasets($ds);
	$cc->border->width(0);
	my $defctx = $cc->get_context('default');
	my $area = Chart::Clicker::Renderer::Area->new(opacity => .4);
	$area->brush->width(1);
	$defctx->renderer($area);
	$defctx->domain_axis->format("%d");
    $defctx->domain_axis->tick_values(\@x_axis);
	$defctx->range_axis->format("%d ");
	$defctx->range_axis->tick_values(\@y_axis);
	$defctx->range_axis->label('Frequency');
	$defctx->domain_axis->label('Coordinates');
	$defctx->renderer->brush->width(1);
	$cc->write_output("$graph");
}	
#__________________________________________________________________________________________
sub motif_position_freq_graph_gnu{
    my ($x_axis_ref,$y_values_ref,$graph,$range) = @_;
    use Chart::Graph::Gnuplot qw(gnuplot);
    use List::Util qw(max);
    my @x_axis = @$x_axis_ref;
    my @y_values = @$y_values_ref;
    #####
    my (@x_middle_range,$middle_range);   
    for my $i(0..($#x_axis)){
        $middle_range = $x_axis[$i] + ($range/2);
        push(@x_middle_range,$middle_range);
    }
    #push (@x_middle_range, "0");
    my @y_axis_range = (0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,250,300,350,450,500,600,650,700,750,800,850,900,950,1000);
    my $y = 0;
    my $ymax;
    foreach $y(@y_axis_range){
        $ymax = $y;
        last if(max(@y_values) < $y);
    }
    my $binRange = $ymax/5;
    my $min =0;
    my $max = $ymax;
    my ($hash_ref, $rank_order) = axis_generation($min,$max,$binRange);
    my @y_axis = @$rank_order;
    my $ymin = 0;
    #Debugging aid - save the temporary files if desired
    #$Chart::Graph::save_tmpfiles = 1;
    #Debugging aid - turn on extra debugging messages
    #$Chart::Graph::debug = 1; 
    # Call and "usual" global parameters
    
    gnuplot({'title' => 'Motif frequency',
        'x-axis label' => 'Coordinates',
        'y-axis label' => 'Frequency',
        'yrange' => [$ymin,$ymax],
        'size' => [0.7,0.5],
        'output type' => 'png',
        'output file' => "$graph",
        # Setting date/time specific options.
        'xtics' => \@x_axis,
        'ytics' => \@y_axis,
        # Set output range - note quoting of date string
        'extra_opts' => 'set grid',
    },
    # Data for when stock opened
    #type +> filledcurves
    [{'title' => 'motif position frequencies',
        'type' => 'columns',
        'style' => 'lines'}, 
    	\@x_middle_range , 
    	\@y_values
    ]
    );
    my @answers = run ("chcon -v --type=httpd_sys_content_t $graph",".");
}
#__________________________________________________________________________________________
sub GO_annotation_graph{
	my ($siteimage,$motifname,$pvalue,$groups) = @_ ; 
	use Chart::Clicker;
	use Chart::Clicker::Context;
	use Chart::Clicker::Data::DataSet;
	use Chart::Clicker::Data::Marker;
	use Chart::Clicker::Data::Series;
	use Chart::Clicker::Decoration::Legend::Tabular;
	use Chart::Clicker::Renderer::Area;
	use Chart::Clicker::Renderer::Pie;
	use Geometry::Primitive::Rectangle;
	use Geometry::Primitive::Circle;
	use Graphics::Color::RGB;
	my $cc = Chart::Clicker->new(width => 300, height => 300);
	my ($series, $series1,$series2,$series3,$series4,$series5,$series6,$series7);
	my ($i, $value, $data, @real_pval);
	for $i(0..$#$pvalue){
		$value = sprintf '%.2e', exp(-$$pvalue[$i]);
		push(@real_pval, $value);
	}
	if(defined($$pvalue[0])){
		$series1 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[0],0],
						name => "$$groups[0]");
	}
	if(defined($$pvalue[1])){
		$series2 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[1],0],
						name => "$$groups[1]");
	}	
	if(defined($$pvalue[2])){
		$series3 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[2],0],
						name => "$$groups[2]");
	}	
	if(defined($$pvalue[3])){
		$series4 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[3],0],
						name => "$$groups[3]");
	}
	if(defined($$pvalue[4])){
		$series5 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[4],0],
						name => "$$groups[4]");
	}
	if(defined($$pvalue[5])){
		$series5 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[5],0],
						name => "$$groups[5]");
	}
	if(defined($$pvalue[6])){
		$series5 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[6],0],
						name => "$$groups[6]");
	}
	if(defined($$pvalue[7])){
		$series5 = Chart::Clicker::Data::Series->new(
						keys => [ 1,2],
						values => [$$pvalue[7],0],
						name => "$$groups[7]");
	}
	if (defined($series1 && $series2 && $series3 && $series4 && $series5 && $series6 && $series7)){
		$series = [ $series1, $series2, $series3, $series4, $series5, $series6, $series7]
	}
	elsif (defined($series1 && $series2 && $series3 && $series4 && $series5 && $series6)){
		$series = [ $series1, $series2, $series3, $series4, $series5, $series6]
	}
	elsif (defined($series1 && $series2 && $series3 && $series4 && $series5)){
		$series = [ $series1, $series2, $series3, $series4, $series5]
	}
	elsif (defined($series1 && $series2 && $series3 && $series4)){
		$series = [ $series1, $series2, $series3, $series4]
	}
	elsif (defined($series1 && $series2 && $series3)){
		$series = [ $series1, $series2, $series3]
	}
	elsif (defined($series1 && $series2)){
		$series = [ $series1, $series2]
	}
	elsif (defined($series1)){
		$series = [ $series1]
	}
	my $ds = Chart::Clicker::Data::DataSet->new(series => $series);
	if (defined($$pvalue[7])){
		$data = [[$real_pval[0]], [$real_pval[1]],[$real_pval[2]], [$real_pval[3]], [$real_pval[4]], [$real_pval[5]], [$real_pval[6]], [$real_pval[7]]];
	}
	elsif(defined($$pvalue[6])){
		$data = [[$real_pval[0]], [$real_pval[1]],[$real_pval[2]], [$real_pval[3]], [$real_pval[4]], [$real_pval[5]], [$real_pval[6]]];
	}
	elsif(defined($$pvalue[5])){
		$data = [[$real_pval[0]], [$real_pval[1]],[$real_pval[2]], [$real_pval[3]], [$real_pval[4]], [$real_pval[5]]];
	}
	elsif(defined($$pvalue[4])){
		$data = [[$real_pval[0]], [$real_pval[1]],[$real_pval[2]], [$real_pval[3]], [$real_pval[4]]];
	}
	elsif(defined($$pvalue[3])){
		$data = [[$real_pval[0]], [$real_pval[1]], [$real_pval[2]], [$real_pval[3]]];
	}
	elsif(defined($$pvalue[2])){
		$data = [[$real_pval[0]], [$real_pval[1]], [$real_pval[2]]];
	}
	elsif(defined($$pvalue[1])){
		$data = [[$real_pval[0]] ,[$real_pval[1]]];
	}
	elsif(defined($$pvalue[0])){
		$data = [[$real_pval[0]]];
	}
	else{ };
	$cc->legend(Chart::Clicker::Decoration::Legend::Tabular->new( 
				header => [ qw(Annotation p-value) ],
				data =>  $data ));
	$cc->add_to_datasets($ds);
	$cc->border->width(0);
	my $defctx = $cc->get_context('default');
	my $ren = Chart::Clicker::Renderer::Pie->new;
	$ren->border_color(Graphics::Color::RGB->new(red => 1, green => 1, blue => 1));
	$ren->brush->width(2);
	$ren->gradient_color(Graphics::Color::RGB->new(red => 1, green => 1, blue => 1, alpha => .3));
	$ren->gradient_reverse(1);
	$defctx->renderer($ren);
	$defctx->domain_axis->hidden(1);
	$defctx->range_axis->hidden(1);
	$cc->plot->grid->visible(0);
	$cc->write_output("$siteimage/GO$motifname.png");
	$series="a";$series1="a";$series2="a";$series3="a";$series4="a";$series5="a";
	unlink($series, $series1,$series2,$series3,$series4,$series5,$cc,$defctx);
}
#__________________________________________________________________________________________
 sub random_int_between {
	my($min, $max) = @_;
	#	Assumes that the two arguments are integers themselves!
	return $min if $min == $max;
	($min, $max) = ($max, $min) if $min > $max;
	return $min + int rand(1 + $max - $min);
}
#__________________________________________________________________________________________
sub maxim {my $max = pop(@_);
	foreach (@_)
		{ $max = $_ if $_ > $max;
		}
	$max;
}
#__________________________________________________________________________________________
sub low_complexity{
	my ($consensus) = @_;
	my @DNAconsensus = split( '', $consensus );
	#	Initialize the counts.
	my $count_of_A =0;
	my $count_of_C= 0;
	my $count_of_G= 0;
	my $count_of_T= 0;
	#	Increment the appropriate count.
	foreach my $base (@DNAconsensus) {
		if     ( $base eq 'A' ) {
			$count_of_A++;
		} elsif ( $base eq 'C' ) {
			$count_of_C++;
		} elsif ( $base eq 'G' ) {
			$count_of_G++;
		} elsif ( $base eq 'T' ) {
			$count_of_T++;
		} else {}
	}
	my $maximum = max $count_of_A,$count_of_C,$count_of_G,$count_of_T;
	my $sum = $count_of_A+$count_of_C+$count_of_G+$count_of_T;
	my $percent = $maximum/$sum;
	return $percent;
}
#____________________________________________________________________________________________________
sub binomial_motif{
	my ($results,$user_seq,$nb_useq, $rand_seq, $nb_rseq,$motif_matrix) = @_;
	my (@Array, $seqmt, $score,@Aseq,@Arand );
	#	Run pscan on sequence to analyse
	#$pscanexec -q $user_seq -M $motif_matrix -ss;
	my @answers = run("pscan -q $user_seq -M $motif_matrix",".");
	my $pscan_user_seq = "$user_seq".".ris";  
	open (SUBIN,"<",$pscan_user_seq) or die("can't open file");
	while ( <SUBIN> ) {
		chomp;
		if(m/^$/){ }
		else{
			if(!m/^>/){ }
			else{
				chomp;
				@Array = split (/\t/,$_);
				$seqmt = substr $Array[0],1;
				$score = $Array[1];
				if ($score >= 0.97){
					push (@Aseq,$seqmt);
				}
			}
		}
	}
	close SUBIN;
	@Aseq= uniq(@Aseq);
	my $nb_umotif = @Aseq;
	# 	Open .ris and count number of sequences by motifs (sequence that contains a motif)
	@answers = run("pscan -q $rand_seq -M  $motif_matrix","."); 
	#	Group same motif results and sorting them
	my $pscan_rand_seq = "$rand_seq".".ris"; 
	open (SUBIN,"<",$pscan_rand_seq);
	while ( <SUBIN> ) {
		chomp;
		if(m/^$/){}
		else{
			if(!m/^>/){ }
			else{
				@Array = split (/\t/,$_);
				$seqmt = substr $Array[0],1;
				$score = $Array[1];
				if ($score >= 0.965){
					push (@Arand,$seqmt);
				}
			}
		}
	}
	close SUBIN;
	@Arand= uniq(@Arand);
	$nb_rseq = $nb_rseq;
	my $g = @Arand;
	my $p = $g/$nb_rseq;
	unlink($pscan_rand_seq,$pscan_user_seq);
	my $res = binomialp_value($nb_useq,$nb_umotif,$p);	
	return $res;
}
#____________________________________________________________________________________________________
sub hypergeometic_motif_MEME{
	my ($results,$display,$user_seq,$nb_useq,$rand_seq,$nb_rseq,$motif_matrix,$bcgd) = @_;
	my (@Array, $seqmt, $score,@Aseq,@Arand,$line );
	my $fimotxt = "$results/fimo$display.txt";
	#run FIMO on sequence to analyse
	my @fimotxt =qx(fimo --bgfile $bcgd --output-pthresh 1e-4 --no-qvalue --text --verbosity 1 $motif_matrix $user_seq);
	for my $l(1..($#fimotxt)){
			@Array = split (/\t/,$fimotxt[$l]);
			$seqmt = $Array[1];
			push (@Aseq,$seqmt);
	}	
	@Aseq= uniq(@Aseq);
	my $nb_umotif = @Aseq;
	unlink(@fimotxt);
	@fimotxt =qx(fimo --bgfile $bcgd --output-pthresh 1e-4 --no-qvalue --text --verbosity 1 $motif_matrix $rand_seq);
	for my $l(1..($#fimotxt)){
			@Array = split (/\t/,$fimotxt[$l]);		
			$seqmt = $Array[1];
			#	FIMO version 4.6.0 single strand 
			#if( $Array[4] =~ m/^\+/){
			#	FIMO version 4.3.0 single strand 
			#if( $Array[0] =~ m/^\+/){
				push (@Arand,$seqmt);
			#}
	}	
	@Arand= uniq(@Arand);
	$nb_rseq = $nb_rseq;
	my $g = @Arand;
	unlink(@fimotxt);
	#print "$nb_rseq,$nb_useq,$nb_umotif,$g\n";
	my $res = hypergeom($g,$nb_rseq,$nb_useq,$nb_umotif);
	return $res;
}
#__________________________________________________________________________________________
sub binomial_Bioprospector_motif{
	my ($results,$user_seq,$nb_useq, $rand_seq, $nb_rseq,$motif_matrix) = @_;
	my (@Array, $seqmt, $score,@Aseq,@Arand );
	#	Run pscan on sequence to analyse
	my @answers = run("pscan -q $user_seq -M $motif_matrix",".");
	#if($answers[2]==0){
	#	diehtml("can't run pscan");
	#}
	my $pscan_user_seq = "$user_seq".".ris"; 
	$pscan_user_seq =~ s/\s+//;
	open (SUBIN,"$pscan_user_seq") or die("can't open Bioprospector PSCAN");
	while ( <SUBIN> ) {
		chomp;
		if(m/^$/){ }
		else{
			if(!m/^>/){ }
			else{
				chomp;
				@Array = split (/\t/,$_);
				$seqmt = substr $Array[0],1;
				$score = $Array[1];
				if ($score >= 0.9){
					push (@Aseq,$seqmt);
					#print "$seqmt\n";
				}
			}
		}
	}
	close SUBIN;
	@Aseq= uniq(@Aseq);
	my $nb_umotif = @Aseq;
	#	Open .ris and count number of sequences by motifs (sequence that contains a motif)
	@answers = run("pscan -q $rand_seq -M  $motif_matrix",".");	
	#	Group same motif results and sorting them
	my $pscan_rand_seq = "$rand_seq".".ris"; 
	open (SUBIN,"<",$pscan_rand_seq);
	while ( <SUBIN> ) {
		chomp;
		if(m/^$/){}
		else{
			if(!m/^>/){ }
			else{
				@Array = split (/\t/,$_);
				$seqmt = substr $Array[0],1;
				$score = $Array[1];
				if ($score >= 0.9){
					push (@Arand,$seqmt);
				}
			}
		}
	}
	close SUBIN;
	@Arand= uniq(@Arand);
	$nb_rseq = $nb_rseq;
	my $g = @Arand;
	my $p = $g/$nb_rseq;
	unlink($pscan_rand_seq,$pscan_user_seq);
	my $res = binomialp_value($nb_useq,$nb_umotif,$p);	
	return $res;
}
#____________________________________________________________________________________________________
sub uniq{
    my %seen = ();
    my @r = ();
    foreach my $a (@_) {
        unless ($seen{$a}) {
            push @r, $a;
            $seen{$a} = 1;
        }
    }
    return @r;
}
#______________________________________________________________________________________________________
=pod
sub hypergeo_p_val{
	my ($n,$k,$Genome,$g) = @_;
	my ($i,$c,$prob_i,$p_value);
	if($n < $g){
		$c = $n;
	}
	else { $c = $g }
	$p_value = 0;
	for $i($k..$c){
		$prob_i = (bc($n,$i)*bc(($Genome-$n),($g-$i)))/ bc($Genome,$g);
		$p_value += $prob_i;	
	}
	if ($p_value>1){
		$p_value = 1;
	}
	return $p_value;
}
=cut
#______________________________________________________________________________________________________
sub binomialp_value{
	my ($n,$k,$p)=@_;
	my $p_value;
	if($n == 0){
		$p_value = 1}
	elsif($k==0){
		$p_value = 1}
	elsif($p == 0){
		$p_value = 0}
	elsif($p == 1){
		$p_value = 1}
	else{$p_value = bc($n,$k)*($p**$k)*((1-$p)**($n-$k));}
	return $p_value;
}
#__________________________________________________________________________________________
sub GO_hypergeometric{
	my ($annotation_count,$nb_gene_w_mf,$nb_gene_genome) = @_;
	my %annotation_count = %$annotation_count;
	my (@Sigf,$pval,$n,$k,$gene_in_genome_annot);

	foreach my $key (keys %annotation_count){
		if ($annotation_count{$key}[0] > 0){
			#total number of gene in list $nb_gene_w_mf
			$k =  $annotation_count{$key}[0]; #gene with motif
			$gene_in_genome_annot = $annotation_count{$key}[1];
			#$pval = hypergeometric($nb_gene_w_mf,$k,$nb_gene_genome,$gene_in_genome_annot);
			my $not_genome_annotated = $nb_gene_genome-$gene_in_genome_annot;
			$pval = hypergeom($gene_in_genome_annot,$not_genome_annotated,$nb_gene_w_mf,$k); #100,40,1000,300 
			#print $pval."\n";
			if($pval<0.01){
				# put significant GO annotation, p-value, description in array of arrays			
				my $descr = $annotation_count{$key}[2];
				push(@Sigf,[$key,$pval,$descr]);
			}
		}
	}
	my @sorted_Sigf = sort { $a->[1] <=> $b->[1] } @Sigf;
	return (\@sorted_Sigf);
}
#__________________________________________________________________________________________
sub annotation_counting{
	my($gene_to_annot,$annotation_count_total,$gene_to_motif) = @_;
	my(%annotation_count);

	my @gene_to_motif = @$gene_to_motif;
	my %gene_to_annot = %$gene_to_annot;
	open(AN,"<","$annotation_count_total") or die("Can't open file");
	while (<AN>){	
		chomp;
		s/^\s+//;
		#print;   	
		my($val,$key) = split /\s/;
		#print "$val\;$key\n" ;
		my $ar = [];	
	  	$annotation_count{$key} = [0, $val,"na",$ar];
	}
	close(AN);
	# Create array gene and annotation
=pod
	open(AN,"<","$annotation_txt")or die("can't open file");
	while (<AN>){
		chomp;
		my($gene,$annot,$desc)= split /\t/;
		if (defined ($gene_to_annot{$gene})){
			push(@{$gene_to_annot{$gene}},[$annot, $desc]);
		}
		else { 
			$gene_to_annot{$gene} = [[$annot, $desc]];
		}
	}
	close (AN);
=cut
	#Create a hash of array with key annotation from Specie_annotation.txt file
	my ($geneid,$GO, $description);
	for my $i(0..$#gene_to_motif){
		#print $gene_to_motif[$i]."\n";
		$geneid = $gene_to_motif[$i];
		if (defined($gene_to_annot{$geneid})) {	 
			for my $j(0..$#{$gene_to_annot{$geneid}}){
				$GO = $gene_to_annot{$geneid}[$j][0];
				$description = $gene_to_annot{$geneid}[$j][1];
				$annotation_count{$GO}[0]++;
				if ("$annotation_count{$GO}[2]" =~ /na/){			
					$annotation_count{$GO}[2]= $description;
				}					
				if(!defined ($annotation_count{$GO}[3])){
					@{$annotation_count{$GO}[3]}= $geneid;
				}else{
					push(@{$annotation_count{$GO}[3]},$geneid);
				}
	 		}
		}	
	}
	#for my $i (0..$#{$annotation_count{$GO}[3]}){
	#	print "in subroutine".$annotation_count{$GO}[3][$i]."\#\#\#\n";
	#}
	return (\%annotation_count);
}
#__________________________________________________________________________________________
sub annotation_hash{
	my($annotation_txt) = @_;
	my %gene_to_annot;	
	open(AN,"<","$annotation_txt")or die("can't open file");
		while (<AN>){
			chomp;
			my($gene,$annot,$desc)= split /\t/;
			if (defined ($gene_to_annot{$gene})){
				push(@{$gene_to_annot{$gene}},[$annot, $desc]);
			}
			else { 
				$gene_to_annot{$gene} = [[$annot, $desc]];
			}
		}
	close (AN);
	return (\%gene_to_annot);
}

#__________________________________________________________________________________________
sub hypergeom {
   # There are m "bad" and n "good" balls in an urn.
   # Pick N of them. The probability of i or more successful selections:
   # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
   my ($n, $m, $N, $i) = @_;
   #print "$n\, $m\, $N\, $i\n";	
   my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
   my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
   return exp($loghyp1 - $loghyp2);
}
#__________________________________________________________________________________________
sub logfact {
   return gammln(shift(@_) + 1.0);
}
#__________________________________________________________________________________________
sub gammln {
  my $xx = shift;
  my @cof = (76.18009172947146, -86.50532032941677,
             24.01409824083091, -1.231739572450155,
             0.12086509738661e-2, -0.5395239384953e-5);
  my $y = my $x = $xx;
  my $tmp = $x + 5.5;
  $tmp -= ($x + .5) * log($tmp);
  my $ser = 1.000000000190015;
  for my $j (0..5) {
     $ser += $cof[$j]/++$y;
  }
  -$tmp + log(2.5066282746310005*$ser/$x);
}
#__________________________________________________________________________________________
sub bc { 
	my $r=1;
	$r*=$_/($_[0]-$_+1)for(1+pop..$_[0]);
	if($r== "inf"){
		$r=1;
	}
	return $r;
}
#__________________________________________________________________________________________
sub MNCP {
	my ($fgd,$bgd)= @_;
	my @fgd = @$fgd;
	my @bgd = @$bgd;
	my(@sortedfgd,@st_total,$n,$rank,$i, $j,$sumrank,$count,$average);
	my ($n_fgd,$n_bgd,$n_total,@total,$slp,@slps,$MNCP);
	@sortedfgd =  sort {$b <=> $a} @$fgd;
	$n = @sortedfgd;
	my (@rankfgd) = (0) x $n;
	$sumrank = 0;
	$count = 0;
	for ($i=0; $i<=$n+1; $i++){	
		$sumrank += $i;
		$count += 1;
		if (defined($sortedfgd[$i+1])){
			if( ($i == $n-1) || (($sortedfgd[$i]) != ($sortedfgd[$i+1])) ){
				$average = ($sumrank/$count)+1;
				for ($j=($i-$count+1); $j<=($i+1);$j++){
					$rankfgd[$j]= $average;
				}
				$sumrank = 0;
				$count = 0;
			}
		}
		else{
			$average = $sumrank/$count+1;
			for ($j=($i-$count+1); $j<$n;$j++){
				$rankfgd[$j]= $average;
			}
		}
	}
	push(@total, (@fgd, @bgd));
	$n = @total;
	@st_total =  sort {$b <=> $a} @total;
	my (@rank_st_total) = (0) x $n;
	$sumrank = 0;
	$count = 0;
	for ($i=0; $i<=$n+1; $i++){
		$sumrank += $i;
		$count += 1;
		if (defined($st_total[$i+1])){
			if( ($i == $n-1) || (($st_total[$i]) != ($st_total[$i+1])) ){
				$average = $sumrank/$count+1;
				for ($j=($i-$count+1); $j<=($i+1);$j++){
					$rank_st_total[$j]= $average;
				}
				$sumrank = 0;
				$count = 0;
			}
		}
		else{
			$average = ($sumrank/$count)+1;
			for ($j=($i-$count+1); $j<$n;$j++){
				$rank_st_total[$j]= $average;
			}
		}
	}
	$n_fgd = @fgd;
	$n_total = @st_total;
	my $stop = 0;
	for ($i=0;$i<=$#sortedfgd ;$i++){
		for ($j=0;$j<=$#st_total;$j++){
			my $i_fdg = int($sortedfgd[$i]);
			my $j_bgd = int($st_total[$j]);
			if($i_fdg == $j_bgd){
				while($stop==0){
					$slp = ($rankfgd[$i] / $n_fgd ) / ( $rank_st_total[$j]/ $n_total);
					push(@slps,$slp);
					$stop = 1;
				};
			}
			else{}
		}
		$stop=0;	
	}
	$MNCP = mean(@slps);
	return $MNCP;
}
#__________________________________________________________________________________________
sub diehtml {
    my ($msg) = @_;
    print qq(<h2>Error</h2>$msg</body></html>);
    exit;
}
#__________________________________________________________________________________________
sub fasta_to_hash{  
	open(FILE, "<@_") or die("Cannot open FASTA file.\n");
	my %seqs;
	my $line;
	my $hash_key;
	open(FILE,"<@_")or die"\n\ncould not open fata file";
	while($line=<FILE>){
		chomp $line;  
		$line=~ s/^\s+//;  
		if(ord$line==62){        
			#if line starts with a >    
			$line =~ m/^>(\w+)/; 
			#matches a word starting with > (Fasta)    
			$hash_key = $1;
			#print "$1\n";
			my $sequence='';    
			while($line=<FILE>){      
				chomp $line;      
				$line=~ s/^\s+//;      
				next if $line =~ m/^#/; 
				#discard comments      
				last if ord$line==62;      
				$line=~ s/\s+$//;      
				$sequence.= $line;         
				#if the line doesnt contain a FASTA header the line is added to the value of the sequence    
			}    
			$seqs{$hash_key}=$sequence unless exists$seqs{$hash_key};    
			redo unless eof(FILE);  
		}
	}
	close(FILE);
	return \%seqs;
} 

1;

