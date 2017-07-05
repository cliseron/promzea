# Promzea
Promzea is an ensemble de-novo cis-elements tools written in perl scripts. The aim is to find cis-elements from maize, rice or Arabidopsis co-regulated genes.
This is the command line version of the web interface promzea.org

# Installation
## Install perl-5.12.5 
wget  http://www.cpan.org/src/5.0/perl-5.12.5.tar.gz
tar –xvf  perl-5.12.5.tar.gz
cd perl-5.12.5
sh Configure -de
make
make test
make install
sudo mv /usr/bin/perl /usr/bin/perlbroken
sudo ln -s /usr/local/bin/perl /usr/bin/perl
cd ..

## Install perl and Bioperl modules
### Bioperl installation
In the home user directory
sudo cpan
install Bundle::CPAN
	> install Algorithm::Munkres
	> install Array::Compare
	> install Clone
	> install Convert::Binary::C 
	> install GD 
	> install Graph
> install GraphViz
> force install Bio::Perl

### Perl module
> sudo cpan
> install FCGI
> install CGI
> install Email::Valid
> install List::Util
> install MIME::Lite
> install Net::SMTP			   
> install Getopt::Long
> install Sys::Info
> install Parallel::ForkManager
> install Sys::Statistics::Linux::MemStats
> install Statistics::Basic’
use Cwd;
> install SMUELLER/PathTools-3.33.tar.gz 
> install IO::CaptureOutput
> install Geometry::Primitive
> install Exception::Class
> install Graphics::Color::RGB
> install Math::Trig 
> install Forest
plus the other Chart-Clicker dependencies 
> install Chart::Clicker 

## Dependent Program 
Please check the Licence of each program before any useage. We disclaiming any responsability in any improper use of the following tiers-program.
### Clover
> wget zlab.bu.edu/~mfrith/downloads/clover-2011-10-24.tar.gz
> tar -xvf clover-2011-10-24.tar.gz
> cd clover-2011-10-24
> sudo make
> sudo cp clover /usr/local/bin
> sudo ln -sfn /usr/local/bin/clover /bin/\
###	pscan 
>wget http://159.149.160.51/pscan/Source/pscan.tar.gz
>tar -xvf pscan.tar.gz
> sudo  vi ~/.bash_profile
	PATH=$PATH:$HOME/bin:/usr/local/include:/usr/local/share:/usr/local/lib:/bin/
	Save  :wq
>sudo g++ pscan.cpp -o pscan -O3 -lgsl –lgslcblas
>sudo ln -sfn /usr/local/bin/pscan /bin/

### weeder
>wget http://159.149.160.51/modtools/weeder1.4.2.tar.gz
>tar xvf weeder1.4.2.tar
>cd Weeder1.4.2
>./compileall
>sudo mv adviser.out FreqFiles weederlauncher.out weederTFBS.out /var/www/cgi-bin/
> perl count_word.pl  6 /var/www/html/data/OS/OS_1000_upstream.fa OS.6.freq
> perl count_word.pl  8 /var/www/html/data/OS/OS_1000_upstream.fa OS.8.freq

### MEME and FIMO
> tar zxf meme_4.6.1.tar.gz
> cd meme_4.6.1
> sudo ./configure --prefix=/
> sudo make
> sudo make test
> sudo make install

###	BioPropector
> wget http://motif.stanford.edu/distributions/bioprospector/BioProspector.2004.zip
> unzip BioProspector.2004.zip
> mv BioProspector.linux BioProspector
> sudo chmod 555 BioProspector
> sudo mv BioProspector /bin
